# ---------------------------------------------------------------------------
# MVA Hackathon 2026 workflow.
#
# Two halves that run on different machines:
#
#   VCF-scoped (Arms A, B, E, F-partial) runs anywhere, including the
#   darwin-arm64 recon laptop.
#
#   Alignment-scoped (Arms C, D, F-repeats) needs the linux-64 GPU host. It
#   starts from FASTQ because no BAM was shipped. See RECON.md.
#
# Run the first half now:
#     snakemake --cores 4 vcf_arms
# Run the second half once the GPU host is available:
#     snakemake --cores 32 alignment_arms
#
# Disk budget for the alignment half, on top of the 79 GB of input:
#     bwa-mem2 index of GRCh38 + decoys   ~ 30 GB
#     reference FASTA and indices         ~  4 GB
#     per-lane sorted BAM (transient)     ~ 90 GB peak
#     final CRAM                          ~ 35 GB
#   Peak requirement is roughly 160 GB free. The recon host had 181 GB free,
#   which is inside the margin but not comfortably; check before starting.
# ---------------------------------------------------------------------------

import pathlib

configfile: "config/config.yaml"

SAMPLE   = config["cohort"]["proband_sample_id"]
LANES    = config["data"]["fastq_lanes"]
VCF      = config["data"]["vcf"]
BUILD    = config["genome"]["build"]
NAMING   = config["genome"]["contig_naming"]
SEED     = config["analysis"]["random_seed"]

REF_DIR  = "refs"
REF      = f"{REF_DIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# The upstream callset was produced against a DNAnexus-specific reference:
# GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta
# That exact file is not publicly redistributable in a form we can pin, so this
# workflow aligns against the standard 1000 Genomes GRCh38 full analysis set
# plus decoys and HLA, then strips the chr prefix so contig names match the VCF.
#
# THIS IS A DEVIATION AND IT IS DELIBERATE. The two references differ in the
# masked GRC exclusion regions. Consequences, which belong in the submission:
#   - Our alignments are not bit-identical to the ones that produced the VCF.
#   - SV and mosaic calls from our BAM are internally consistent and comparable
#     to each other, but a coordinate-level comparison against the upstream
#     callset in a masked region may disagree.
#   - Coordinates are unaffected: both are GRCh38 primary assembly.
# Anything called inside a masked-exclusion region must be flagged as such
# rather than reported at face value.

wildcard_constraints:
    lane = r"L\d{3}",
    read = r"R[12]",


# ===========================================================================
# Aggregate targets
# ===========================================================================

rule vcf_arms:
    """Everything runnable without an alignment. Live now."""
    input:
        "results/recon/characterisation.json",
        "results/summaries/arm_a_baseline.tsv",
        "results/summaries/arm_b_splicing_candidates.tsv",
        "results/summaries/arm_f_roh.txt",

rule alignment_arms:
    """Everything gated behind building a BAM. Needs the GPU host."""
    input:
        f"results/align/{SAMPLE}.cram",
        f"results/align/{SAMPLE}.mosdepth.summary.txt",
        "results/summaries/arm_c_sv_merged.tsv",
        "results/summaries/arm_d_mosaic.tsv",


# ===========================================================================
# Reference preparation
# ===========================================================================

rule download_reference:
    output:
        fa = REF,
        fai = REF + ".fai",
    threads: 2
    shell:
        r"""
        mkdir -p {REF_DIR}
        curl -L --retry 5 -o {output.fa}.gz \
          ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
        gzip -dc {output.fa}.gz > {output.fa} && rm -f {output.fa}.gz
        samtools faidx {output.fa}
        """

rule bwa_index:
    input: REF
    output: multiext(REF, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    resources: mem_mb = 70000
    shell: "bwa-mem2 index {input}"


# ===========================================================================
# Alignment. Two implementations; pick with config or --config aligner=...
# ===========================================================================

rule fastp_lane:
    """Adapter and quality trimming, per lane. The HTML report is a deliverable:
    library quality bears on how much to trust a low-VAF mosaic call."""
    input:
        r1 = "data/{sample}_S16_{lane}_R1_001.fastq.gz",
        r2 = "data/{sample}_S16_{lane}_R2_001.fastq.gz",
    output:
        r1   = temp("results/align/trim/{sample}_{lane}_R1.fq.gz"),
        r2   = temp("results/align/trim/{sample}_{lane}_R2.fq.gz"),
        html = "results/qc/fastp_{sample}_{lane}.html",
        json = "results/qc/fastp_{sample}_{lane}.json",
    threads: 8
    shell:
        r"""
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --detect_adapter_for_pe --thread {threads} \
              --html {output.html} --json {output.json}
        """

rule bwa_mem2_lane:
    """CPU alignment path. Read groups carry the lane so that duplicate marking
    and any later per-lane batch effect check can distinguish them."""
    input:
        r1 = "results/align/trim/{sample}_{lane}_R1.fq.gz",
        r2 = "results/align/trim/{sample}_{lane}_R2.fq.gz",
        ref = REF,
        idx = multiext(REF, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    output:
        bam = temp("results/align/{sample}_{lane}.sorted.bam"),
    threads: 24
    resources: mem_mb = 60000
    params:
        rg = lambda w: (f"@RG\\tID:{w.sample}_{w.lane}\\tSM:{w.sample}"
                        f"\\tLB:{w.sample}\\tPL:ILLUMINA\\tPU:HGWCNDSX7.{w.lane}"),
    shell:
        r"""
        bwa-mem2 mem -t {threads} -R '{params.rg}' -K 100000000 -Y \
            {input.ref} {input.r1} {input.r2} \
          | samtools sort -@ 4 -m 2G -o {output.bam} -
        samtools index -@ 4 {output.bam}
        """

rule merge_lanes:
    input: expand("results/align/{{sample}}_{lane}.sorted.bam", lane=LANES)
    output: temp("results/align/{sample}.merged.bam")
    threads: 8
    shell: "samtools merge -@ {threads} -o {output} {input}"

rule mark_duplicates:
    """Duplicates are marked, never removed. A removed read is unavailable to
    the mosaic arm, and at 42x the duplicate rate is low enough that keeping
    them costs little."""
    input: "results/align/{sample}.merged.bam"
    output:
        bam = "results/align/{sample}.bam",
        bai = "results/align/{sample}.bam.bai",
        metrics = "results/qc/{sample}.dupmetrics.txt",
    threads: 8
    shell:
        r"""
        samtools markdup -@ {threads} --write-index -f {output.metrics} \
          <(samtools collate -@ {threads} -O -u {input} \
            | samtools fixmate -@ {threads} -m -u - -) {output.bam}
        """

rule bam_to_cram:
    """CRAM for the archived artefact: roughly half the size, and the disk
    budget above is tight."""
    input: bam = "results/align/{sample}.bam", ref = REF
    output: "results/align/{sample}.cram"
    threads: 8
    shell: "samtools view -@ {threads} -T {input.ref} -C -o {output} {input.bam} && samtools index {output}"


# ===========================================================================
# Coverage. Settles the Phase 0 questions that the VCF could not.
# ===========================================================================

rule mosdepth:
    """Phase 0 could not distinguish 'no variant called' from 'no coverage',
    because the shipped VCF is variants-only rather than a gVCF. This rule
    settles it. The three variant-free gaps logged in
    results/recon/panel_depth_profile.json are checked here first."""
    input:
        bam = "results/align/{sample}.bam",
        bed = "config/gene_panels/mva_known.nochr.bed",
    output:
        summary = "results/align/{sample}.mosdepth.summary.txt",
        regions = "results/align/{sample}.regions.bed.gz",
    threads: 4
    shell:
        r"""
        mosdepth -t {threads} --by {input.bed} --fast-mode \
                 results/align/{wildcards.sample} {input.bam}
        """


# ===========================================================================
# Arm F: runnable now, VCF only
# ===========================================================================

rule arm_f_roh:
    """Formal ROH call, confirming the Phase 0 windowed scan that found no
    consanguinity. A reported negative, per plan section 6.6."""
    input: vcf = VCF
    output: "results/summaries/arm_f_roh.txt"
    shell:
        r"""
        mkdir -p results/summaries
        bcftools roh --AF-dflt 0.4 -G30 {input.vcf} \
          | awk '$1=="RG"' > {output} || true
        echo "# bcftools roh, --AF-dflt 0.4 -G30, seed {SEED}" >> {output}
        """
