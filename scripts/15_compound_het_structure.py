#!/usr/bin/env python3
"""Compound-heterozygote structure across the mitotic panel.

Why now, before VEP and gnomAD
------------------------------
The leading hypothesis (RECON.md class 1) is a compound heterozygote with one
cryptic allele. Whether that is even structurally possible in a given gene is
answerable from the callset alone: it needs at least two heterozygous variants,
and ideally evidence that two of them sit on opposite haplotypes.

In a singleton there are no parents to phase against, so the only native phasing
signal is HaplotypeCaller's ``PGT``/``PID`` tags, recorded at Phase 0 as
load-bearing for this branch. Two heterozygous variants sharing a ``PID`` are in
the same phasing group; if their ``PGT`` values differ (``0|1`` against ``1|0``)
they are on **opposite haplotypes**, which is exactly the *trans* configuration a
recessive compound heterozygote requires.

This is a structural screen, not a diagnosis. Without gnomAD frequencies most of
these variants will be common polymorphisms. The point is to find which genes can
carry the hypothesis at all, and to rank where the splicing arm should look first
once SpliceAI is available.

Output is aggregate by design: per-gene counts, phasing-group structure and class
distributions. No individual variant position, allele or genotype is written, in
keeping with CLAUDE.md rule 1.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, "src")

from mva.track1.regions import GeneModel  # noqa: E402

VCF = "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"


def gene_records(vcf: str, contig: str, start: int, end: int) -> list[dict]:
    """Stream one gene's records with the phasing tags."""
    out = subprocess.run(
        ["bcftools", "query", "-r", f"{contig}:{start}-{end}",
         "-f", "%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t[%DP]\t[%PGT]\t[%PID]\n", vcf],
        capture_output=True, text=True, check=True).stdout
    recs = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        recs.append({"pos": int(f[0]), "ref": f[1], "alt": f[2], "filter": f[3],
                     "gt": f[4], "dp": f[5], "pgt": f[6], "pid": f[7]})
    return recs


def analyse_gene(recs: list[dict], model: GeneModel, gene_sym: str) -> dict:
    hets, homs = [], []
    for r in recs:
        a = r["gt"].replace("|", "/").split("/")
        if len(a) != 2 or "." in a:
            continue
        (hets if a[0] != a[1] else homs if a[0] != "0" else []).append(r)

    # Phasing groups: PID identifies a group, PGT the haplotype assignment.
    groups: dict[str, list[dict]] = collections.defaultdict(list)
    for r in hets:
        if r["pid"] and r["pid"] != ".":
            groups[r["pid"]].append(r)

    trans_pairs = 0
    cis_pairs = 0
    for members in groups.values():
        pgts = [m["pgt"] for m in members if m["pgt"] and m["pgt"] != "."]
        for i in range(len(pgts)):
            for j in range(i + 1, len(pgts)):
                if pgts[i] != pgts[j]:
                    trans_pairs += 1
                else:
                    cis_pairs += 1

    n_phased = sum(len(v) for v in groups.values())
    return {
        "gene": gene_sym,
        "n_variants": len(recs),
        "n_het": len(hets),
        "n_hom_alt": len(homs),
        "n_het_phased": n_phased,
        "n_het_unphased": len(hets) - n_phased,
        "n_phasing_groups": len(groups),
        "largest_phasing_group": max((len(v) for v in groups.values()), default=0),
        "het_pairs_in_trans": trans_pairs,
        "het_pairs_in_cis": cis_pairs,
        # Can this gene carry a compound heterozygote at all?
        "can_carry_comp_het": len(hets) >= 2,
        "has_trans_evidence": trans_pairs > 0,
        "n_filtered_out": sum(1 for r in recs if r["filter"] not in ("PASS", ".")),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", default="refs/Homo_sapiens.GRCh38.115.gtf.gz")
    ap.add_argument("--panel", default="config/gene_panels/mitotic_extended.tsv")
    ap.add_argument("--out", default="results/summaries/compound_het_structure.md")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.panel, newline=""), delimiter="\t"))
    core = {r["symbol"] for r in rows if r["in_core_panel"] == "yes"}
    known = {r["symbol"] for r in rows if r["known_mva_gene"] == "yes"}
    meta = {r["symbol"]: r for r in rows}

    sys.stderr.write("building gene model...\n")
    gm = GeneModel.from_gtf(args.gtf, symbols=core, flank=2000)

    results = []
    for i, g in enumerate(gm.genes, 1):
        if i % 50 == 0:
            sys.stderr.write(f"  {i}/{len(gm)}\r")
        try:
            recs = gene_records(VCF, g.contig, g.start, g.end)
        except subprocess.CalledProcessError:
            continue
        a = analyse_gene(recs, gm, g.symbol)
        a["span_kb"] = round((g.end - g.start) / 1000, 1)
        a["variants_per_kb"] = round(a["n_variants"] / max(1, (g.end - g.start) / 1000), 3)
        a["known_mva"] = g.symbol in known
        a["panel_tier"] = meta.get(g.symbol, {}).get("panel_tier", "")
        results.append(a)
    sys.stderr.write("\n")

    pathlib.Path("results/summaries").mkdir(parents=True, exist_ok=True)
    json.dump(results, open("results/recon/compound_het_structure.json", "w"), indent=1)

    dens = sorted(r["variants_per_kb"] for r in results)
    median_dens = dens[len(dens) // 2]

    L = ["# Compound-heterozygote structure across the mitotic panel\n"]
    L.append("Generated by `scripts/15_compound_het_structure.py`. Aggregate counts only.\n")
    L.append("## Method\n")
    L.append("For each panel gene, count heterozygous calls and inspect HaplotypeCaller's")
    L.append("`PGT`/`PID` physical phasing tags. Two heterozygous variants sharing a `PID`")
    L.append("with differing `PGT` lie on opposite haplotypes, which is the *trans*")
    L.append("configuration a recessive compound heterozygote requires. In a singleton")
    L.append("this is the only native phasing signal available.\n")
    L.append("This is a structural screen. Without gnomAD frequencies most of these")
    L.append("variants are common polymorphisms; the question here is only which genes")
    L.append("*can* carry the hypothesis, and where the splicing arm should look first.\n")

    L.append("## Known MVA genes\n")
    L.append("| gene | span kb | variants | var/kb | het | hom | phasing groups | het in trans | can carry comp het |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in sorted([x for x in results if x["known_mva"]], key=lambda x: -x["n_het"]):
        L.append(f"| {r['gene']} | {r['span_kb']} | {r['n_variants']} | {r['variants_per_kb']} | "
                 f"{r['n_het']} | {r['n_hom_alt']} | {r['n_phasing_groups']} | "
                 f"{r['het_pairs_in_trans']} | {'yes' if r['can_carry_comp_het'] else 'NO'} |")
    L.append("")

    no_ch = [r["gene"] for r in results if r["known_mva"] and not r["can_carry_comp_het"]]
    trans = [r for r in results if r["known_mva"] and r["has_trans_evidence"]]
    L.append("### Reading\n")
    if no_ch:
        L.append(f"**{', '.join(no_ch)} cannot carry a compound heterozygote** in this")
        L.append("callset: fewer than two heterozygous calls in the gene body. Unless a")
        L.append("second allele is a structural variant or sits outside the called region,")
        L.append("the compound-heterozygote hypothesis is not available for these genes.")
        L.append("Note that this is exactly what a homozygous deletion of one allele would")
        L.append("also look like, which Arm C must check once a BAM exists.\n")
    if trans:
        L.append("Genes with read-backed *trans* evidence between heterozygous calls: "
                 + ", ".join(f"**{r['gene']}** ({r['het_pairs_in_trans']} pairs)" for r in trans)
                 + ".\n")
    else:
        L.append("**No known MVA gene shows a read-backed *trans* pair.** Physical phasing")
        L.append("only resolves variants close enough to share a read or assembly region,")
        L.append("so this is expected for variants further apart than a read pair, and is")
        L.append("not evidence against a compound heterozygote. It does mean phase for any")
        L.append("candidate pair will have to come from population-based phasing or from")
        L.append("targeted long-read or Sanger work.\n")

    L.append("## Panel-wide context\n")
    L.append(f"- genes analysed: **{len(results)}**")
    L.append(f"- median variant density across the panel: **{median_dens} per kb**")
    L.append(f"- genes able to carry a compound heterozygote: "
             f"**{sum(1 for r in results if r['can_carry_comp_het'])}**")
    L.append(f"- genes with read-backed *trans* evidence: "
             f"**{sum(1 for r in results if r['has_trans_evidence'])}**")
    L.append(f"- total heterozygous *trans* pairs panel-wide: "
             f"**{sum(r['het_pairs_in_trans'] for r in results)}**\n")

    L.append("### Where the known MVA genes sit in the panel density distribution\n")
    L.append("Phase 0 noted that BUB1B, CEP57, TRIP13 and BUB1 carry far fewer variants per")
    L.append("kb than BUB3, CEP192 or SMC5, and left it unexplained. Ranking every panel")
    L.append("gene by density puts that in context:\n")
    L.append("| gene | var/kb | percentile within the panel |")
    L.append("|---|---:|---:|")
    for r in sorted([x for x in results if x["known_mva"]], key=lambda x: x["variants_per_kb"]):
        pct = 100 * sum(1 for d in dens if d < r["variants_per_kb"]) / len(dens)
        L.append(f"| {r['gene']} | {r['variants_per_kb']} | {pct:.0f} |")
    L.append("")

    top = sorted(results, key=lambda r: -r["het_pairs_in_trans"])[:10]
    if top and top[0]["het_pairs_in_trans"]:
        L.append("### Panel genes with the most read-backed *trans* pairs\n")
        L.append("Not candidates by themselves: these are gene bodies dense enough in")
        L.append("heterozygous variation for physical phasing to resolve pairs. Listed so")
        L.append("the splicing arm knows where phase will and will not be recoverable.\n")
        L.append("| gene | het | trans pairs | known MVA |")
        L.append("|---|---:|---:|---|")
        for r in top:
            if r["het_pairs_in_trans"]:
                L.append(f"| {r['gene']} | {r['n_het']} | {r['het_pairs_in_trans']} | "
                         f"{'yes' if r['known_mva'] else ''} |")
        L.append("")

    pathlib.Path(args.out).write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
