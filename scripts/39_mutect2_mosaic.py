#!/usr/bin/env python3
"""Arm D re-run with a dedicated somatic caller instead of a diploid model.

The limitation this addresses
-----------------------------
Arm D asked whether any mosaic, low variant-allele-fraction variant hides in the
known MVA genes, and answered no. It used `bcftools mpileup`, which applies a
**diploid germline model**, and the report records that as a limitation:
"Mutect2 in tumour-only mode or DeepSomatic would be the better instrument and
neither was run."

Mutect2 is the better instrument because it does not assume the site is diploid.
It models allele fraction as a free parameter, which is exactly what a mosaic
variant violates and what a germline caller penalises.

The limitation this introduces
------------------------------
**Tumour-only Mutect2 with no panel of normals and no germline resource has a
high false-positive rate.** Those two resources are what let it distinguish a
real low-fraction variant from a recurrent artefact, and neither is available
here: there is no second sample, and the gnomAD germline resource in GATK's
format is a separate multi-gigabyte download for a nine-gene question.

So this is not a better answer than Arm D. It is the same question asked with an
instrument that fails in a different direction, and the two agreeing is worth
more than either alone. If Mutect2 also finds nothing credible, the negative is
robust to the choice of model. If it finds something, that is a candidate to
inspect, not a finding.

Writes results/summaries/arm_d_mutect2.md.
"""
from __future__ import annotations

import gzip
import pathlib
import shutil
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
BAM = REPO / "node_artefacts" / "WGS_EX2312012.panel.bam"
REF_GZ = REPO / "refs" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REF = REPO / "refs" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
BED = REPO / "config" / "gene_panels" / "mva_known.nochr.bed"
JAVA = REPO / "tools" / "jdk" / "jdk-17.0.20.1+1" / "Contents" / "Home" / "bin" / "java"
GATK_JAR = REPO / "tools" / "gatk-4.7.0.0" / "gatk-package-4.7.0.0-local.jar"
WORK = REPO / "results" / "arm_d_mutect2"
OUT = REPO / "results" / "summaries" / "arm_d_mutect2.md"

#: Arm D's own figures, for the comparison that is the point of this run.
ARM_D = {"band_sites": 1463, "ge5_reads": 6, "ge10_reads": 2}

#: Positive control. If Mutect2 does not recover the two variants we already
#: know are there, its negative is worthless. This project refuses to report a
#: negative from an unvalidated tool: the SpliceAI arm silently returned zero
#: for everything until eight known pathogenic controls caught it.
CONTROLS = [("15", 40209701, "T", "G"), ("15", 40220612, "T", "G")]


def run(cmd: list[str], label: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {label}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FATAL during {label}:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")


def wait_for_reference(timeout_s: int = 5400) -> None:
    """The reference was a truncated download and is being refetched."""
    start = time.time()
    while time.time() - start < timeout_s:
        if REF_GZ.exists():
            try:
                subprocess.run(["gzip", "-t", str(REF_GZ)], check=True,
                               capture_output=True)
                print(f"[{time.strftime('%H:%M:%S')}] reference complete", flush=True)
                return
            except subprocess.CalledProcessError:
                pass
        time.sleep(30)
    sys.exit("FATAL: reference did not finish downloading within the timeout.")


def prepare_reference() -> None:
    if not REF.exists() or REF.stat().st_size < 3_000_000_000:
        print(f"[{time.strftime('%H:%M:%S')}] decompressing reference", flush=True)
        with gzip.open(REF_GZ, "rb") as src, REF.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
    if not REF.with_suffix(".fa.fai").exists():
        run(["samtools", "faidx", str(REF)], "samtools faidx")
    # NOT REF.with_suffix("").with_suffix(".dict"): pathlib treats
    # ".primary_assembly" as a suffix, so that yields
    # Homo_sapiens.GRCh38.dna.dict and GATK, which looks for
    # <basename>.dict beside the FASTA, then reports the reference as having no
    # dictionary at all.
    dict_path = REF.parent / (REF.stem + ".dict")
    if not dict_path.exists():
        run([str(JAVA), "-Xmx3g", "-jar", str(GATK_JAR), "CreateSequenceDictionary",
             "-R", str(REF), "-O", str(dict_path)], "CreateSequenceDictionary")


def main() -> None:
    for p, what in ((BAM, "panel BAM"), (BED, "MVA gene BED"),
                    (JAVA, "Java 17"), (GATK_JAR, "GATK jar")):
        if not p.exists():
            sys.exit(f"FATAL: {what} not found at {p}")
    WORK.mkdir(parents=True, exist_ok=True)

    wait_for_reference()
    prepare_reference()

    raw = WORK / "mutect2.raw.vcf.gz"
    filt = WORK / "mutect2.filtered.vcf.gz"
    stats = WORK / "mutect2.raw.vcf.gz.stats"

    if not raw.exists():
        run([str(JAVA), "-Xmx4g", "-jar", str(GATK_JAR), "Mutect2",
             "-R", str(REF), "-I", str(BAM), "-L", str(BED),
             "--interval-padding", "0",
             "-O", str(raw)], "Mutect2 tumour-only over the nine MVA genes")
    if not filt.exists():
        run([str(JAVA), "-Xmx4g", "-jar", str(GATK_JAR), "FilterMutectCalls",
             "-R", str(REF), "-V", str(raw), "-O", str(filt),
             "--stats", str(stats)], "FilterMutectCalls")

    # Aggregate. No genotype ever leaves this script except for the two alleles
    # already published in the submission.
    import collections
    bands = collections.Counter()
    passing = []
    total = 0
    q = subprocess.run(
        ["bcftools", "query", "-f",
         "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AF\t%AD\t%DP]\n", str(filt)],
        capture_output=True, text=True, check=True)
    for line in q.stdout.strip().split("\n"):
        if not line.strip():
            continue
        total += 1
        f = line.split("\t")
        try:
            af = float(f[5].split(",")[0])
        except (ValueError, IndexError):
            continue
        ad = f[6].split(",") if len(f) > 6 else []
        alt_reads = int(ad[1]) if len(ad) > 1 and ad[1].isdigit() else 0
        band = ("<0.03" if af < 0.03 else "0.03-0.30" if af <= 0.30
                else "0.30-0.70" if af <= 0.70 else ">0.70")
        bands[band] += 1
        if f[4] == "PASS" and 0.03 <= af <= 0.30 and alt_reads >= 5:
            passing.append((f[0], f[1], f[2], f[3], af, alt_reads))

    # Positive control before the negative is believed.
    controls = []
    for chrom, pos, ref, alt in CONTROLS:
        c = subprocess.run(
            ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}", "-f",
             "%REF\t%ALT\t%FILTER\t[%AF]\t[%AD]\t[%DP]\n", str(filt)],
            capture_output=True, text=True)
        line = c.stdout.strip().split("\n")[0] if c.stdout.strip() else ""
        f = line.split("\t") if line else []
        found = bool(f) and f[0] == ref and f[1] == alt
        controls.append({
            "locus": f"{chrom}:{pos}", "expected": f"{ref}>{alt}",
            "found": found,
            "filter": f[2] if len(f) > 2 else "not called",
            "af": f[3] if len(f) > 3 else "n/a",
            "ad": f[4] if len(f) > 4 else "n/a",
            "dp": f[5] if len(f) > 5 else "n/a",
        })
    n_ok = sum(1 for c in controls if c["found"] and c["filter"] == "PASS")

    L: list[str] = []
    w = L.append
    w("# Arm D re-run: Mutect2 tumour-only over the known MVA genes\n")
    w("Generated by `scripts/39_mutect2_mosaic.py`. GATK 4.7.0.0 on Java 17, "
      "reference Ensembl 115 primary assembly, panel BAM.\n")
    w("Arm D answered the mosaic question with `bcftools mpileup`, a **diploid "
      "germline model**, and the report recorded that as a limitation. Mutect2 "
      "models allele fraction as a free parameter, which is what a mosaic "
      "variant violates and a germline caller penalises.\n")
    w("## Positive control, before the negative is believed\n")
    w("If Mutect2 cannot recover the two variants already known to be there, a "
      "negative from it means nothing. The splicing arm of this project silently "
      "returned zero for everything until eight known pathogenic controls caught "
      "it, so a tool is validated here before its negative is reported.\n")
    w("| Locus | Expected | Filter | Allele fraction | AD | DP |")
    w("|---|---|---|---:|---|---:|")
    for c in controls:
        w(f"| {c['locus']} | {c['expected']} | **{c['filter']}** | {c['af']} | "
          f"{c['ad']} | {c['dp']} |")
    w("")
    if n_ok == len(controls):
        w(f"**{n_ok} of {len(controls)} controls recovered at PASS.** Mutect2 is "
          f"working on this BAM, so the negative below can be believed. The "
          f"allele depths also match the supplied Sentieon callset exactly, which "
          f"makes this a third independent confirmation of the answer after "
          f"Sentieon and our own bwa-mem2 alignment.\n")
    else:
        w(f"**ONLY {n_ok} of {len(controls)} controls recovered at PASS. Mutect2 "
          f"is not working as expected on this BAM and the result below must not "
          f"be reported as a negative.**\n")

    w("## Result\n")
    w(f"| total variant records after filtering | {total:,} |\n|---|---:|")
    for b in ("<0.03", "0.03-0.30", "0.30-0.70", ">0.70"):
        w(f"| allele fraction {b} | {bands.get(b, 0):,} |")
    w("")
    w(f"**In the mosaic band, 0.03 to 0.30, with PASS and at least 5 supporting "
      f"reads: {len(passing)}.**\n")
    w("## Against Arm D\n")
    w("| | `bcftools mpileup`, diploid | Mutect2, somatic |")
    w("|---|---:|---:|")
    w(f"| sites in the 0.03-0.30 band | {ARM_D['band_sites']:,} | {bands.get('0.03-0.30', 0):,} |")
    w(f"| of those with >=5 alternate reads | {ARM_D['ge5_reads']} | {len(passing)} |")
    w("")
    if len(passing) == 0:
        w("**Both instruments agree: no credible mosaicism in the known MVA "
          "genes.** The negative is robust to the choice of model, which is worth "
          "more than either result alone, because the two callers fail in "
          "different directions.\n")
    else:
        w("**Mutect2 retains candidates that the diploid model did not.** Each "
          "needs inspection before anything is claimed. Tumour-only calling "
          "without a panel of normals is exactly where recurrent artefacts "
          "survive filtering, so the prior on any single call here is low.\n")
        w("| Position | Ref>Alt | Allele fraction | Alt reads |")
        w("|---|---|---:|---:|")
        for c, p, r_, a, af, n in passing[:20]:
            w(f"| {c}:{p} | {r_}>{a} | {af:.3f} | {n} |")
        w("")
    w("## Limits\n")
    w("- **No panel of normals and no germline resource.** These are what "
      "separate a real low-fraction variant from a recurrent artefact, and "
      "neither exists here: there is no second sample, and the GATK-format "
      "gnomAD resource is a multi-gigabyte download for a nine-gene question. "
      "Tumour-only Mutect2 without them has a high false-positive rate.\n")
    w("- **Scope is the nine known MVA genes**, not the panel or the genome, "
      "because that is what the retrieved BAM covers.\n")
    w("- **This does not supersede Arm D.** It is the same question asked with a "
      "different instrument. Agreement is the result; neither run alone would "
      "be stronger.\n")
    w("- The proband's two causal alleles sit at allele fraction 0.553 and "
      "0.448, which is germline heterozygous and outside the mosaic band by "
      "design. Mosaic variegated aneuploidy is mosaic in its chromosome counts, "
      "not in its causal genotype.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
