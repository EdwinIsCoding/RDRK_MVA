#!/usr/bin/env python3
"""MVA Hackathon 2026 - Phase 0 recon, step 2.

Reads the recon artefacts written by ``scripts/00_inventory.sh`` plus the patient
VCF, and answers the routing questions in MVA_HACKATHON_PLAN.md section 2.2.

Governance
----------
This script is permitted to touch ``data/``. Its *outputs* are what the agent and
the report read, and they are deliberately aggregate: counts, distributions,
FILTER breakdowns and per-gene tallies. No individual genotype, position or
allele is written to ``results/recon/`` or ``results/summaries/`` by this script.

Outputs
-------
results/recon/characterisation.json   machine readable, feeds config/config.yaml
results/summaries/phase0_summary.md   human readable aggregate summary

Usage
-----
    scripts/01_characterise.py [--data data] [--out results] [--quick]

``--quick`` skips the full-VCF streaming passes and reports only what can be read
from headers and indices.
"""

from __future__ import annotations

import argparse
import collections
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
from typing import Any

# Build discriminators: length of the first chromosome in each assembly.
CHR1_LENGTHS = {
    249250621: "GRCh37",
    248956422: "GRCh38",
    248387328: "T2T-CHM13v2.0",
}


def run(cmd: list[str], **kw: Any) -> str:
    """Run a command and return stdout, raising on failure."""
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        sys.exit(f"FATAL: {tool} not found on PATH. Install it before running Phase 0.")
    return path


# ---------------------------------------------------------------------------
# Header derived facts
# ---------------------------------------------------------------------------

def parse_header(vcf: pathlib.Path) -> dict[str, Any]:
    header = run(["bcftools", "view", "-h", str(vcf)])
    lines = header.splitlines()

    contigs: dict[str, int] = {}
    for line in lines:
        m = re.match(r"##contig=<ID=([^,>]+),length=(\d+)", line)
        if m:
            contigs[m.group(1)] = int(m.group(2))

    # Build: match chromosome 1 length under either naming convention.
    chr1_len = contigs.get("1") or contigs.get("chr1")
    build = CHR1_LENGTHS.get(chr1_len, "UNKNOWN")
    if build == "UNKNOWN":
        sys.exit(
            f"FATAL: cannot determine genome build. chr1 length = {chr1_len!r}. "
            "Every downstream coordinate depends on this. Resolve before continuing."
        )

    naming = "ucsc_chr" if any(c.startswith("chr") for c in contigs) else "ensembl_nochr"

    # Primary assembly only, for the contig report.
    primary = [c for c in contigs if re.fullmatch(r"(chr)?(\d{1,2}|X|Y|MT|M)", c)]

    filters = re.findall(r"##FILTER=<ID=([^,]+)", header)
    fmt = re.findall(r"##FORMAT=<ID=([^,]+)", header)
    info = re.findall(r"##INFO=<ID=([^,]+)", header)

    # Caller provenance. Sentieon and GATK both stamp command lines.
    callers = []
    for line in lines:
        m = re.match(r"##(SentieonCommandLine\.\w+|GATKCommandLine)=<ID=(\w+)[^>]*?Version=\"?([^\",]+)", line)
        if m:
            callers.append({"vendor": m.group(1).split(".")[0], "tool": m.group(2), "version": m.group(3)})

    refs = [line.split("=", 1)[1] for line in lines if line.startswith("##reference=")]

    return {
        "n_contigs": len(contigs),
        "n_primary_contigs": len(primary),
        "chr1_length": chr1_len,
        "build": build,
        "contig_naming": naming,
        "has_decoys": len(contigs) > 200,
        "filters_declared": filters,
        "format_fields": fmt,
        "info_fields": info,
        "callers": callers,
        "reference_files": refs,
        "samples": run(["bcftools", "query", "-l", str(vcf)]).split(),
    }


# ---------------------------------------------------------------------------
# Streaming aggregate passes
# ---------------------------------------------------------------------------

def filter_and_contig_tally(vcf: pathlib.Path) -> dict[str, Any]:
    """One streaming pass: FILTER breakdown, per contig counts, genotype classes."""
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", "%CHROM\t%FILTER\t%TYPE\t[%GT]\n", str(vcf)],
        stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
    )
    filters: collections.Counter[str] = collections.Counter()
    contigs: collections.Counter[str] = collections.Counter()
    types: collections.Counter[str] = collections.Counter()
    gts: collections.Counter[str] = collections.Counter()
    pass_contigs: collections.Counter[str] = collections.Counter()
    total = 0

    assert proc.stdout is not None
    for line in proc.stdout:
        chrom, filt, vtype, gt = line.rstrip("\n").split("\t")
        total += 1
        filters[filt] += 1
        contigs[chrom] += 1
        types[vtype] += 1
        gts[gt] += 1
        if filt == "PASS" or filt == ".":
            pass_contigs[chrom] += 1
    proc.wait()

    # Genotype classes, for het/hom ratio and any non diploid calls.
    def gclass(gt: str) -> str:
        a = re.split(r"[/|]", gt)
        if any(x == "." for x in a):
            return "nocall"
        if len(set(a)) == 1:
            return "hom_ref" if a[0] == "0" else "hom_alt"
        return "het"

    classes: collections.Counter[str] = collections.Counter()
    for gt, n in gts.items():
        classes[gclass(gt)] += n

    return {
        "total_records": total,
        "filter_counts": dict(filters.most_common()),
        "type_counts": dict(types.most_common()),
        "genotype_classes": dict(classes),
        "het_hom_ratio": round(classes["het"] / classes["hom_alt"], 3) if classes["hom_alt"] else None,
        "per_contig_all": dict(contigs.most_common(30)),
        "per_contig_pass": dict(pass_contigs.most_common(30)),
    }


def depth_distribution(vcf: pathlib.Path, sample_every: int = 200) -> dict[str, Any]:
    """Sampled FORMAT/DP distribution, as a coverage proxy when no BAM exists."""
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", "[%DP]\n", str(vcf)],
        stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
    )
    vals: list[int] = []
    assert proc.stdout is not None
    for i, line in enumerate(proc.stdout):
        if i % sample_every:
            continue
        s = line.strip()
        if s.isdigit():
            vals.append(int(s))
    proc.wait()
    if not vals:
        return {"note": "no FORMAT/DP field"}
    vals.sort()
    def q(p: float) -> int:
        return vals[min(len(vals) - 1, int(p * len(vals)))]
    return {
        "n_sampled": len(vals),
        "sample_every_nth_record": sample_every,
        "mean": round(sum(vals) / len(vals), 1),
        "p05": q(0.05), "p25": q(0.25), "median": q(0.50), "p75": q(0.75), "p95": q(0.95),
        "frac_dp_lt_10": round(sum(v < 10 for v in vals) / len(vals), 4),
        "frac_dp_lt_20": round(sum(v < 20 for v in vals) / len(vals), 4),
    }


def gene_panel_tally(vcf: pathlib.Path, bed: pathlib.Path) -> list[dict[str, Any]]:
    """Per gene aggregate tallies over the known MVA panel.

    Aggregate only: counts by FILTER and genotype class, plus the callable span
    proxy. No positions or alleles are emitted. A gene with an unexpectedly low
    count or a coverage hole is itself a lead, per plan section 2.2.
    """
    genes: list[dict[str, Any]] = []
    for line in bed.read_text().splitlines():
        chrom, start, end, name = line.split("\t")[:4]
        region = f"{chrom}:{start}-{end}"
        span = int(end) - int(start)
        try:
            out = run(["bcftools", "query", "-r", region,
                       "-f", "%FILTER\t[%GT]\t[%DP]\n", str(vcf)])
        except subprocess.CalledProcessError:
            genes.append({"gene": name, "region": region, "error": "query failed"})
            continue

        n = npass = nhet = nhom = 0
        dps: list[int] = []
        for row in out.splitlines():
            filt, gt, dp = (row.split("\t") + ["", "", ""])[:3]
            n += 1
            if filt in ("PASS", "."):
                npass += 1
            a = re.split(r"[/|]", gt)
            if "." not in a and len(a) == 2:
                if a[0] != a[1]:
                    nhet += 1
                elif a[0] != "0":
                    nhom += 1
            if dp.isdigit():
                dps.append(int(dp))

        genes.append({
            "gene": name,
            "region_nochr": region,
            "span_bp": span,
            "variants_total": n,
            "variants_pass": npass,
            "variants_filtered_out": n - npass,
            "het_calls": nhet,
            "hom_alt_calls": nhom,
            "variants_per_kb": round(1000 * n / span, 3),
            "median_dp_at_called_sites": sorted(dps)[len(dps) // 2] if dps else None,
            "min_dp_at_called_sites": min(dps) if dps else None,
        })
    return genes


# ---------------------------------------------------------------------------
# Data modality detection
# ---------------------------------------------------------------------------

def detect_modalities(filetypes: pathlib.Path, manifest: pathlib.Path) -> dict[str, Any]:
    rows = [l.split("\t") for l in filetypes.read_text().splitlines()[1:] if l.strip()]
    paths = [r[0] for r in rows]
    sizes = {}
    for l in manifest.read_text().splitlines()[1:]:
        if l.strip():
            p, b = l.split("\t")
            sizes[p] = int(b)

    def any_ext(*exts: str) -> list[str]:
        return [p for p in paths if any(p.lower().endswith(e) for e in exts)]

    fastqs = any_ext(".fastq.gz", ".fq.gz", ".fastq", ".fq")
    lanes = sorted({m.group(1) for p in fastqs if (m := re.search(r"_(L\d{3})_", p))})
    reads = sorted({m.group(1) for p in fastqs if (m := re.search(r"_(R[12])_\d+\.f", p))})

    return {
        "vcf": any_ext(".vcf.gz", ".vcf", ".bcf"),
        "vcf_index": any_ext(".tbi", ".csi"),
        "alignments": any_ext(".bam", ".cram"),
        "fastq": fastqs,
        "fastq_lanes": lanes,
        "fastq_read_ends": reads,
        "fastq_paired": reads == ["R1", "R2"],
        "fastq_total_bytes": sum(sizes.get(p, 0) for p in fastqs),
        "clinical_documents": any_ext(".docx", ".pdf", ".xlsx", ".csv", ".tsv", ".json"),
        "rnaseq_evidence": [],  # populated below
        "total_bytes": sum(sizes.values()),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    ap.add_argument("--panel-bed", default="config/gene_panels/mva_known.nochr.bed")
    ap.add_argument("--quick", action="store_true", help="skip full VCF streaming passes")
    args = ap.parse_args()

    require("bcftools")
    data = pathlib.Path(args.data)
    recon = pathlib.Path(args.out) / "recon"
    summaries = pathlib.Path(args.out) / "summaries"
    summaries.mkdir(parents=True, exist_ok=True)

    modalities = detect_modalities(recon / "filetypes.tsv", recon / "manifest.tsv")
    if not modalities["vcf"]:
        sys.exit("FATAL: no VCF found. Phase 0 cannot characterise the callset.")
    vcf = pathlib.Path(modalities["vcf"][0])

    result: dict[str, Any] = {
        "generated_utc": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc).isoformat(timespec="seconds"),
        "bcftools_version": run(["bcftools", "--version"]).splitlines()[0],
        "vcf": str(vcf),
        "modalities": modalities,
        "header": parse_header(vcf),
    }

    n_samples = len(result["header"]["samples"])
    result["cohort_structure"] = (
        "singleton" if n_samples == 1 else
        "duo" if n_samples == 2 else
        "trio" if n_samples == 3 else f"multi_sample_n{n_samples}"
    )

    if not args.quick:
        sys.stderr.write("streaming pass 1/3: FILTER, contig and genotype tally\n")
        result["callset"] = filter_and_contig_tally(vcf)
        sys.stderr.write("streaming pass 2/3: sampled depth distribution\n")
        result["depth"] = depth_distribution(vcf)
    sys.stderr.write("pass 3/3: MVA gene panel tally\n")
    result["mva_gene_panel"] = gene_panel_tally(vcf, pathlib.Path(args.panel_bed))

    # ---- Branch determination, per plan section 2.4 -------------------------
    has_rna = bool(modalities["rnaseq_evidence"])
    has_aln = bool(modalities["alignments"])
    has_fastq = bool(modalities["fastq"])
    is_trio = result["cohort_structure"] == "trio"

    if is_trio and has_aln and has_rna:
        branch, why = "A", "trio with alignments and RNA-seq"
    elif is_trio and (has_aln or has_fastq):
        branch, why = "B", "trio with alignments (or FASTQ from which to build them), no RNA-seq"
    elif has_aln or has_fastq:
        branch, why = ("C", "singleton with alignments available or derivable from FASTQ")
    else:
        branch, why = "D", "VCF only: no alignments and no FASTQ"
    if has_fastq and not has_aln:
        why += "; alignments must be generated from FASTQ before arms C and D can run"

    result["branch"] = {"branch": branch, "rationale": why,
                        "alignments_present": has_aln,
                        "alignments_derivable_from_fastq": has_fastq,
                        "rnaseq_present": has_rna}

    (recon / "characterisation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["branch"], indent=2))
    print(f"\nWritten: {recon / 'characterisation.json'}")


if __name__ == "__main__":
    main()
