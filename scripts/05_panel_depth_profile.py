#!/usr/bin/env python3
"""10 kb resolution variant-density and depth profile across each panel gene
plus 300 kb of flanking sequence.

Purpose: separate a sharp gene-specific coverage collapse (which in a recessive
chromosomal-instability disorder is a candidate homozygous or hemizygous
deletion) from a broad regional depression (mappability or GC) and from a short
run of homozygosity. Aggregate counts only.
"""
import collections, json, subprocess

VCF = "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
BED = "config/gene_panels/mva_known.nochr.bed"
WIN, FLANK = 10_000, 300_000

genes = []
for line in open(BED):
    c, s, e, name = line.split()[:4]
    genes.append((name, c, int(s), int(e)))

report = {}
for name, c, s, e in genes:
    lo, hi = max(0, s - FLANK), e + FLANK
    out = subprocess.run(
        ["bcftools", "query", "-r", f"{c}:{lo}-{hi}", "-f", "%POS\t[%GT]\t[%DP]\n", VCF],
        capture_output=True, text=True).stdout

    n = collections.Counter(); dps = collections.defaultdict(list); het = collections.Counter()
    for row in out.splitlines():
        pos, gt, d = (row.split("\t") + ["", ""])[:3]
        w = int(pos) // WIN
        n[w] += 1
        a = gt.replace("|", "/").split("/")
        if len(a) == 2 and "." not in a and a[0] != a[1]:
            het[w] += 1
        if d.isdigit():
            dps[w].append(int(d))

    flank_w = [w for w in n if not (s // WIN <= w <= e // WIN)]
    gene_w  = [w for w in range(s // WIN, e // WIN + 1)]
    def med(v): 
        v = sorted(v); return v[len(v)//2] if v else None
    flank_density = med([n[w] for w in flank_w]) if flank_w else None
    flank_dp = med([d for w in flank_w for d in dps[w]])
    gene_dp = med([d for w in gene_w for d in dps[w]])

    empty = [w for w in gene_w if n[w] == 0]
    # Longest contiguous run of empty 10 kb windows inside the gene body.
    runs, cur = [], []
    for w in gene_w:
        if n[w] == 0: cur.append(w)
        elif cur: runs.append(cur); cur = []
    if cur: runs.append(cur)
    longest = max(runs, key=len) if runs else []

    report[name] = {
        "region": f"{c}:{s}-{e}",
        "gene_windows": len(gene_w),
        "empty_gene_windows": len(empty),
        "longest_empty_run_kb": len(longest) * WIN // 1000,
        "longest_empty_run_region": (f"{c}:{longest[0]*WIN}-{(longest[-1]+1)*WIN}" if longest else None),
        "median_variants_per_10kb_flank": flank_density,
        "median_variants_per_10kb_gene": med([n[w] for w in gene_w]),
        "median_DP_flank": flank_dp,
        "median_DP_gene": gene_dp,
        "dp_ratio_gene_vs_flank": round(gene_dp / flank_dp, 2) if gene_dp and flank_dp else None,
    }

json.dump(report, open("results/recon/panel_depth_profile.json", "w"), indent=2)

cols = ["gene", "gene_win", "empty", "run_kb", "var/10kb_flank", "var/10kb_gene", "DP_flank", "DP_gene", "DP_ratio"]
print("".join(f"{c:>16}" for c in cols))
for g, r in report.items():
    print("".join(f"{str(v):>16}" for v in [
        g, r["gene_windows"], r["empty_gene_windows"], r["longest_empty_run_kb"],
        r["median_variants_per_10kb_flank"], r["median_variants_per_10kb_gene"],
        r["median_DP_flank"], r["median_DP_gene"], r["dp_ratio_gene_vs_flank"]]))
print()
for g, r in report.items():
    if r["longest_empty_run_kb"] and r["longest_empty_run_kb"] >= 20:
        print(f"  GAP  {g}: {r['longest_empty_run_kb']} kb with no called variant at {r['longest_empty_run_region']}")
