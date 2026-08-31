# PROVENANCE.md

Everything needed to establish what went in, what produced the outputs, and at
what version. Regenerate the input table with `scripts/07_write_provenance.py`.

## 1. Input files

Source: the challenge data distribution. Not redistributed by this repository.

<!-- BEGIN INPUT TABLE -->

*Regenerated 2026-08-31 12:19 UTC by `scripts/07_write_provenance.py`.*

| File | Bytes | SHA256 |
|---|---:|---|
| `Challenge_Clinical_Phenotype_1.docx` | 16,865 | `0b8129496f239d766c09506f86c23ca2392d7f5f1db13ee8ee513688fea584b7` |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` | 315,153,971 | `bec64c68d3ad10d05b8071d07cdf87d217013d66d3ba02e24b505035c9d42175` |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` | 2,343,376 | `6f8fed62f11c475fc63a8e2b50925ffc7be33b6930225e821cb977951761a2e0` |
| `WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz` | 10,467,073,561 | `5ec244c0648552f2b23ae5b1a0b1350e8e464eb5cccd40178a109c44c174939f` |
| `WGS_EX2312012_HGWCNDSX7_S16_L001_R2_001.fastq.gz` | 11,035,873,766 | `b80f2cc9f29c229a2bfe557988a4b93c74494fc0d9108b3752bc04cf57239cab` |
| `WGS_EX2312012_HGWCNDSX7_S16_L002_R1_001.fastq.gz` | 10,375,801,767 | `0f60e090d9873a18a6914d374fa070a20502aa546d05b51ec520bd9e06ae7e71` |
| `WGS_EX2312012_HGWCNDSX7_S16_L002_R2_001.fastq.gz` | 10,911,120,034 | `992ee50db8a89709b7adb327daabcacbf58261503f8c439bb18a3c682e43cd91` |
| `WGS_EX2312012_HGWCNDSX7_S16_L003_R1_001.fastq.gz` | 10,195,738,985 | `bbeb51349d4411a661c5922e44650b357d968551ef85735ebb39e6c6071e5d46` |
| `WGS_EX2312012_HGWCNDSX7_S16_L003_R2_001.fastq.gz` | 10,621,354,344 | `f3f08f04cd17b089ceb8a8598634f3eb28bc299dc1c00cbfa8e24ba2adb9017b` |
| `WGS_EX2312012_HGWCNDSX7_S16_L004_R1_001.fastq.gz` | 10,306,364,082 | `ebc3a6853b9afb5b839ba1e4ff7966eb4f716fba5fee5a39c526d1468079f438` |
| `WGS_EX2312012_HGWCNDSX7_S16_L004_R2_001.fastq.gz` | 10,755,107,565 | `b16eeae1f84c57fdb889349c25d625eb96a9cb265059004081630131e1037d8b` |
| **11 files** | **84,985,948,316** | |

<!-- END INPUT TABLE -->

TODO(source): record the HuggingFace dataset revision hash for the distribution
these files came from. It was not captured at download time and is needed for a
complete provenance chain.

## 2. Upstream provenance of the callset

Read from the VCF header, preserved at
`results/recon/WGS_EX2312012_HGWCNDSX7.vcf.gz.header.txt`.

| Stage | Tool | Version | Date |
|---|---|---|---|
| Per-sample calling, GVCF mode | Sentieon Haplotyper | sentieon-genomics-202308.02 | 2025-02-04 |
| Joint genotyping | Sentieon GVCFtyper | sentieon-genomics-202308.02 | 2025-02-05 |
| Hard filtering | GATK VariantFiltration | 4.2.4.0 | 2025-02-05 |

Reference: `GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta`
Known sites: `Homo_sapiens_assembly38_no_chr.dbsnp138.vcf.gz` (dbSNP build 138)
Executed on DNAnexus. Calling thresholds `--call_conf 10 --emit_conf 10`.

## 3. Phase 0 tooling

| Tool | Version | Role |
|---|---|---|
| bcftools | 1.24 (htslib 1.24) | VCF headers, streaming queries, statistics |
| samtools | 1.24 (htslib 1.24) | Installed; no alignments present to use it on |
| pandoc | 3.x | Clinical document to text |
| Python | 3.x (Homebrew, arm64) | Characterisation and verification scripts |

Installed via Homebrew on 31 August 2026.

TODO(source): pin exact pandoc and Python patch versions once the container is
built in Phase 1, and record them from inside the container rather than the host.

## 4. External resources queried

| Resource | Endpoint | Date | What was sent |
|---|---|---|---|
| Ensembl REST | `/lookup/symbol/homo_sapiens` (batch POST) | 2026-08-31 | Nine public gene symbols. No patient data. |

Raw response preserved at `results/recon/ensembl_mva_genes.json`.

TODO(source): Ensembl release number was not captured in the response body.
Record it before the coordinates are used in any published table.

## 4a. Analysis tooling on the compute node, 1 September 2026

| Tool | Version | Role |
|---|---|---|
| Ensembl VEP | 116.1, merged cache, GRCh38 | Consequence annotation over all 5,012,204 records |
| bwa-mem2 | 2.3 | FASTQ to BAM, 4 lanes, 14 threads, 4h10m |
| samtools / bcftools / htslib | 1.24 | Sorting, duplicate marking, pileup, region queries |
| mosdepth | 0.3.14 | Panel coverage |
| SpliceAI | 1.3.1 | Splicing prediction, with a documented NumPy 2 compatibility shim |
| Reference (alignment and VEP) | Ensembl 115 GRCh38 primary assembly, no-chr naming | Matches the callset convention |
| gnomAD | v4.1 exomes and genomes sites, remote tabix slices | Allele frequency and homozygote counts |
| gnomAD constraint | v4.1 (autosomes), v2.1.1 (chrX/Y, which v4.1 omits) | Gene constraint |

**Reference deviation, stated because it affects reproducibility.** The upstream
callset was produced against a DNAnexus-specific reference
(`GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr`).
That exact file is not publicly pinnable, so our alignment used the Ensembl 115
primary assembly. Both are GRCh38 primary assembly so coordinates are unaffected,
but the two differ in masked regions and decoy content, and our BAM is therefore
not bit-identical to the one that produced the supplied VCF.

## 5. Databases not yet snapshotted

Phase 1 work. Each needs a version and a snapshot date in
`config/db_versions.yaml` before it is used: VEP cache, gnomAD, ClinVar,
AlphaMissense, CADD, REVEL, SpliceAI model weights hash, Pangolin weights,
Exomiser data, HPO release, SIGNOR, OmniPath, Reactome, Open Targets, ChEMBL,
DrugBank, Broad Repurposing Hub, LINCS/CMap.

## 6. Phase 0 outputs and the script that produced each

| Output | Script |
|---|---|
| `results/recon/manifest.tsv`, `sha256.txt`, `filetypes.tsv`, headers, stats | `scripts/00_inventory.sh` |
| `results/recon/characterisation.json` | `scripts/01_characterise.py` |
| `results/recon/Challenge_Clinical_Phenotype_1.docx.md` | `scripts/02_extract_clinical.sh` |
| `results/recon/roh_proxy.json` | `scripts/03_roh_proxy.py` |
| index verification (stdout) | `scripts/04_verify_panel_tally.py` |
| `results/recon/panel_depth_profile.json` | `scripts/05_panel_depth_profile.py` |
| gap background rates (stdout) | `scripts/06_gap_background_rate.py` |
| `config/gene_panels/mva_known.tsv`, `.nochr.bed` | Ensembl REST, see `config/gene_panels/panel_provenance.md` |
| `config/config.yaml` | Phase 0, hand-assembled from the above |

Random seed for all downstream analysis: `20261024`, set in
`config/config.yaml` under `analysis.random_seed`.
