#!/usr/bin/env python3
"""Arm B: splicing prediction over the Arm A shortlist.

Plan section 6.2, the highest-prior arm. This runs the shortlist rather than the
panel, which is what makes it feasible without a GPU: SpliceAI at a plus or minus
500 bp window over 408 genes is days of compute, but over roughly a dozen
variants it is minutes even on CPU.

The distance parameter
----------------------
``-D 500`` rather than the SpliceAI default of 50. Plan section 6.2 is explicit
that the default window "is precisely why these variants get missed clinically":
a cryptic acceptor created 200 bp into an intron is invisible at plus or minus
50 bp. Every shortlisted variant here is intronic or untranslated, so the wide
window is the whole point.

What a score here does and does not mean
----------------------------------------
A high SpliceAI delta is a **prediction** that a variant alters splicing. It is
not an observation. There is no RNA-seq for this proband, so no aberrant
junction can be confirmed, and plan section 0.4 and CLAUDE.md both forbid
presenting the prediction as though it were a finding. Every output row carries
that qualification.

Thresholds, from SpliceAI's own publication rather than chosen here:
    0.2  permissive, high recall
    0.5  recommended
    0.8  high precision
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import shutil
import subprocess
import sys

sys.path.insert(0, "src")

SHORTLIST = "results/arm_a_shortlist.tsv"
FASTA = "refs/mva_chroms.fa.gz"
ANNOT = "grch38"


def build_vcf(shortlist: pathlib.Path, out: pathlib.Path, fai: pathlib.Path) -> int:
    """Write the shortlist as a minimal VCF for SpliceAI to consume."""
    contigs = []
    if fai.exists():
        for line in fai.read_text().splitlines():
            f = line.split("\t")
            contigs.append(f"##contig=<ID={f[0]},length={f[1]}>")
    n = 0
    with out.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in contigs:
            fh.write(c + "\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        rows = list(csv.DictReader(shortlist.open(newline=""), delimiter="\t"))
        rows.sort(key=lambda r: (r["chrom"], int(r["pos"])))
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\t.\t.\n")
            n += 1
    return n


def parse_spliceai(vcf: pathlib.Path) -> list[dict]:
    """Parse the SpliceAI INFO field.

    Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    The four DS_ fields are delta scores for acceptor gain, acceptor loss,
    donor gain and donor loss. The maximum of the four is the headline number.
    """
    out = []
    for line in vcf.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        info = f[7]
        if "SpliceAI=" not in info:
            continue
        for ann in info.split("SpliceAI=")[1].split(";")[0].split(","):
            p = ann.split("|")
            if len(p) < 10:
                continue
            try:
                ds = [float(x) if x not in ("", ".") else 0.0 for x in p[2:6]]
                dp = [int(x) if x not in ("", ".") else 0 for x in p[6:10]]
            except ValueError:
                continue
            labels = ["acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss"]
            best = max(range(4), key=lambda i: ds[i])
            out.append({
                "chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": p[0],
                "gene": p[1],
                "ds_acceptor_gain": ds[0], "ds_acceptor_loss": ds[1],
                "ds_donor_gain": ds[2], "ds_donor_loss": ds[3],
                "max_delta": ds[best], "max_type": labels[best],
                "max_offset_bp": dp[best],
            })
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shortlist", default=SHORTLIST)
    ap.add_argument("--fasta", default=FASTA)
    ap.add_argument("--distance", type=int, default=500)
    ap.add_argument("--outdir", default="results/arm_b")
    args = ap.parse_args()

    short = pathlib.Path(args.shortlist)
    if not short.exists():
        sys.exit(f"FATAL: {short} not found. Run scripts/18_arm_a_shortlist.py first.")
    fasta = pathlib.Path(args.fasta)
    if not fasta.exists():
        sys.exit(f"FATAL: reference {fasta} not found. Run 'make downloads'.")
    if shutil.which("spliceai") is None:
        sys.exit("FATAL: spliceai not on PATH. Install per pyproject.toml [splicing].")

    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    in_vcf, out_vcf = out / "shortlist.vcf", out / "shortlist.spliceai.vcf"

    n = build_vcf(short, in_vcf, pathlib.Path(str(fasta) + ".fai"))
    print(f"{n} shortlisted variants -> {in_vcf}")

    cmd = ["spliceai", "-I", str(in_vcf), "-O", str(out_vcf),
           "-R", str(fasta), "-A", ANNOT, "-D", str(args.distance)]
    print("running: " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"spliceai failed:\n{r.stderr[-3000:]}")

    rows = parse_spliceai(out_vcf)
    json.dump(rows, open(out / "spliceai_scores.json", "w"), indent=1)

    # Aggregate summary only; the variant-level file stays gitignored.
    bands = {">=0.8 (high precision)": 0, "0.5-0.8 (recommended)": 0,
             "0.2-0.5 (permissive)": 0, "<0.2 (no signal)": 0}
    by_gene: dict[str, float] = {}
    for r_ in rows:
        d = r_["max_delta"]
        k = (">=0.8 (high precision)" if d >= 0.8 else
             "0.5-0.8 (recommended)" if d >= 0.5 else
             "0.2-0.5 (permissive)" if d >= 0.2 else "<0.2 (no signal)")
        bands[k] += 1
        by_gene[r_["gene"]] = max(by_gene.get(r_["gene"], 0.0), d)

    L = ["# Arm B: splicing prediction over the Arm A shortlist\n"]
    L.append(f"Generated by `scripts/19_arm_b_splicing.py`. "
             f"SpliceAI at `-D {args.distance}`, against the SpliceAI default of 50.\n")
    L.append("## Result\n")
    L.append(f"{n} shortlisted variants scored, {len(rows)} gene-level annotations.\n")
    L.append("| SpliceAI delta band | n |")
    L.append("|---|---:|")
    for k, v in bands.items():
        L.append(f"| {k} | {v} |")
    L.append("")
    L.append("Maximum delta per gene:\n")
    L.append("| gene | max delta |")
    L.append("|---|---:|")
    for g, d in sorted(by_gene.items(), key=lambda kv: -kv[1]):
        L.append(f"| {g} | {d:.3f} |")
    L.append("")
    hits = [r_ for r_ in rows if r_["max_delta"] >= 0.2]
    L.append("## Reading this\n")
    if not hits:
        L.append("**No shortlisted variant reaches even the permissive 0.2 threshold.**\n")
        L.append("That is a real negative and worth stating plainly: within the nine known")
        L.append("MVA genes, no rare variant in this proband is predicted to alter splicing.")
        L.append("It does not exclude the cryptic-allele hypothesis, because the shortlist")
        L.append("covers only the known genes and only variants the caller emitted. It does")
        L.append("mean the hypothesis is not supported where it was most expected.\n")
    else:
        L.append(f"{len(hits)} variant(s) reach the permissive 0.2 threshold.\n")
        L.append("**These are predictions, not observations.** There is no RNA-seq for this")
        L.append("proband, so no aberrant junction has been seen. The confirmatory experiment")
        L.append("is RT-PCR across the affected exon junction on patient fibroblast RNA,")
        L.append("which the MVA Society may be able to provide.\n")
    L.append("## Limits\n")
    L.append("- Scope is the nine known MVA genes, not the 408-gene panel. A cryptic allele")
    L.append("  in a novel gene would not appear here.")
    L.append("- Only variants present in the callset are scored. A variant the caller never")
    L.append("  emitted, including anything inside a structural variant, is invisible.")
    L.append("- SpliceAI scores splicing, not branch points or polypyrimidine tracts")
    L.append("  directly. Plan section 6.2 also asks for BPP or LaBranchoR and UTRannotator,")
    L.append("  which are not run here.")

    pathlib.Path("results/summaries").mkdir(parents=True, exist_ok=True)
    pathlib.Path("results/summaries/arm_b_splicing.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
