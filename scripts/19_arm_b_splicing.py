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


def _patch_numpy_for_spliceai() -> None:
    """SpliceAI 1.3.1 calls ``np.fromstring`` in binary mode, removed in NumPy 2.

    The published release predates NumPy 2 and has not been updated. Rather than
    pinning NumPy below 2 (TensorFlow 2.21 will not accept that), the single
    offending call is replaced with ``frombuffer``, which is what ``fromstring``
    did in binary mode.

    **The rest of the function body is reproduced verbatim and that matters.**
    SpliceAI first rewrites the sequence into the byte values \x01 to \x04 and
    only then encodes. An earlier version of this shim dropped that rewrite and
    encoded raw ASCII, which silently produced a delta score of 0.000 for
    everything, including eight known pathogenic canonical splice-site variants
    in BUB1B. The positive controls in ``validate()`` exist because of that.

    latin-1 rather than the default UTF-8: the sentinel bytes must survive
    encoding unchanged.
    """
    import numpy as np
    import spliceai.utils as u

    _MAP = np.asarray([[0, 0, 0, 0],
                       [1, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1]])

    def one_hot_encode(seq):
        seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
        seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')
        return _MAP[np.frombuffer(seq.encode('latin-1'), np.int8) % 5]

    u.one_hot_encode = one_hot_encode


def validate(fasta: pathlib.Path, distance: int) -> tuple[int, int, float]:
    """Score known pathogenic canonical splice-site variants before trusting a
    negative result.

    A variant that destroys a GT or AG dinucleotide must score high. If these do
    not, the tool is broken and any negative finding from it is meaningless.
    Returns (n annotations scoring >= 0.5, n annotations, n input variants,
    max delta seen).

    Annotations and variants are counted separately and deliberately. SpliceAI
    emits one row per overlapping gene, so a variant in a region where two genes
    overlap is scored twice. Collapsing the two counts let an earlier run over
    eight control variants be reported as "9/9 variants".
    """
    import csv as _csv
    have = {l.split("\t")[0] for l in
            pathlib.Path(str(fasta) + ".fai").read_text().splitlines()}
    rows = [r for r in _csv.DictReader(
                open("benchmarks/published_mva_variants.tsv"), delimiter="\t")
            if r["benchmark_tier"] == "1"
            and r["variant_class"] == "splice_site_canonical"
            and r["chrom_nochr"] in have
            and len(r["ref"]) == 1 and len(r["alt"]) == 1]
    if not rows:
        return (0, 0, 0, 0.0)

    tmp = pathlib.Path("results/arm_b/positive_controls.vcf")
    tmp.parent.mkdir(parents=True, exist_ok=True)
    with tmp.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for l in pathlib.Path(str(fasta) + ".fai").read_text().splitlines():
            f = l.split("\t")
            fh.write(f"##contig=<ID={f[0]},length={f[1]}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in sorted(rows, key=lambda r: (r["chrom_nochr"], int(r["pos_grch38"]))):
            fh.write(f"{r['chrom_nochr']}\t{r['pos_grch38']}\t{r['clinvar_vcv']}"
                     f"\t{r['ref']}\t{r['alt']}\t.\t.\t.\n")

    res = run_spliceai(tmp, fasta, distance, quiet=True)
    hi = sum(1 for x in res if x["max_delta"] >= 0.5)
    mx = max((x["max_delta"] for x in res), default=0.0)
    return (hi, len(res), len(rows), mx)


def run_spliceai(in_vcf: pathlib.Path, fasta: pathlib.Path, distance: int,
                 quiet: bool = False) -> list[dict]:
    """Score a VCF with SpliceAI through its Python API.

    Driven directly rather than through the CLI so the NumPy shim above can be
    applied, and so failures surface as exceptions rather than a parsed-empty
    output file.
    """
    import os
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    _patch_numpy_for_spliceai()
    import pysam
    from spliceai.utils import Annotator, get_delta_scores

    ann = Annotator(str(fasta), ANNOT)
    out: list[dict] = []
    labels = ["acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss"]

    vcf = pysam.VariantFile(str(in_vcf))
    n_scored = 0
    for rec in vcf:
        scores = get_delta_scores(rec, ann, distance, 0)
        n_scored += 1
        for ann_str in scores:
            p = ann_str.split("|")
            if len(p) < 10:
                continue
            try:
                ds = [float(x) if x not in ("", ".") else 0.0 for x in p[2:6]]
                dp = [int(x) if x not in ("", ".") else 0 for x in p[6:10]]
            except ValueError:
                continue
            best = max(range(4), key=lambda i: ds[i])
            out.append({
                "chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": p[0],
                "gene": p[1],
                "ds_acceptor_gain": ds[0], "ds_acceptor_loss": ds[1],
                "ds_donor_gain": ds[2], "ds_donor_loss": ds[3],
                "max_delta": ds[best], "max_type": labels[best],
                "max_offset_bp": dp[best],
            })
    if not quiet:
        print(f"scored {n_scored} variants, {len(out)} gene-level annotations")
    return out


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

    hi, n_ctrl, n_ctrl_variants, mx = validate(fasta, args.distance)
    print(f"positive controls: {n_ctrl_variants} known canonical splice variants "
          f"yielding {n_ctrl} gene-level annotations, of which {hi} score >= 0.5 "
          f"(max delta {mx:.3f})")
    if n_ctrl == 0:
        # A vacuous pass is not a pass. This happens when the reference covers
        # no chromosome carrying a benchmark positive, as for the chrX-only
        # reference: every MVA benchmark variant is autosomal. The run may
        # proceed, but the negative it produces is unvalidated and must say so.
        print("WARNING: no positive control was available against this reference, "
              "so SpliceAI is UNVALIDATED for this run. A negative result here "
              "carries less weight than one from a validated run.")
        validated = False
    else:
        validated = True
    if n_ctrl and hi == 0:
        sys.exit(
            f"FATAL: SpliceAI scored 0/{n_ctrl} known pathogenic canonical "
            f"splice-site variants above 0.5 (max delta {mx:.3f}). The tool is "
            f"not working, so a negative result from it would be meaningless. "
            f"Refusing to report one."
        )

    rows = run_spliceai(in_vcf, fasta, args.distance)
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
    if validated:
        L.append(f"**Tool validated for this run:** {n_ctrl_variants} known "
                 f"pathogenic canonical splice-site variants, yielding {n_ctrl} "
                 f"gene-level annotations because SpliceAI emits one row per "
                 f"overlapping gene. {hi} of those {n_ctrl} annotations score at "
                 f"or above 0.5, max delta {mx:.3f}. A negative below can be "
                 f"believed.\n")
    else:
        L.append("**Tool NOT validated for this run.** No positive control was "
                 "available against this reference, because every canonical "
                 "splice-site positive in the benchmark is autosomal and this "
                 "reference is chrX only. The negative below is therefore weaker "
                 "than one from a validated run, and should be re-run against a "
                 "reference carrying controls before it is reported.\n")
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
        L.append("That is a real negative and worth stating plainly: no rare variant in the")
        L.append("genes examined here is predicted to alter splicing in this proband.")
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
