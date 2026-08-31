#!/usr/bin/env python3
"""Windowed heterozygosity scan: distinguishes coverage holes from runs of
homozygosity (ROH). Aggregate output only, 1 Mb windows.

A coverage hole depletes het AND hom calls together. An ROH depletes het while
hom calls stay normal. The distinction reroutes the recessive model: extensive
ROH raises the prior on a homozygous-by-descent single allele over a compound
heterozygote with a cryptic second hit.
"""
import collections, subprocess, sys, json

WIN = 1_000_000
vcf = sys.argv[1] if len(sys.argv) > 1 else "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
autosomes = [str(i) for i in range(1, 23)]

het = collections.Counter(); hom = collections.Counter(); tot = collections.Counter()
p = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%POS\t%FILTER\t[%GT]\n", vcf],
                     stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in p.stdout:
    chrom, pos, filt, gt = line.rstrip("\n").split("\t")
    if chrom not in autosomes or filt not in ("PASS", "."):
        continue
    key = (chrom, int(pos) // WIN)
    tot[key] += 1
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        continue
    if a[0] != a[1]:
        het[key] += 1
    elif a[0] != "0":
        hom[key] += 1
p.wait()

rows = []
for key in sorted(tot, key=lambda k: (int(k[0]), k[1])):
    h, m, t = het[key], hom[key], tot[key]
    rows.append({"chrom": key[0], "win_mb": key[1], "het": h, "hom_alt": m,
                 "total": t, "het_frac": round(h / t, 4) if t else None})

# Genome-wide reference values from windows with adequate variant density.
dense = [r for r in rows if r["total"] >= 300]
med_het_frac = sorted(r["het_frac"] for r in dense)[len(dense) // 2]
med_total = sorted(r["total"] for r in dense)[len(dense) // 2]

json.dump({"window_bp": WIN, "n_windows": len(rows),
           "median_total_per_window": med_total,
           "median_het_fraction": med_het_frac,
           "windows": rows}, open("results/recon/roh_proxy.json", "w"))

print(f"windows: {len(rows)}  median variants/Mb: {med_total}  median het fraction: {med_het_frac}")

# ROH candidates: normal density, strongly depleted heterozygosity.
roh = [r for r in rows if r["total"] >= 0.5 * med_total and r["het_frac"] < 0.35 * med_het_frac]
low = [r for r in rows if r["total"] < 0.25 * med_total]
print(f"\nROH-like windows (normal density, het fraction < 35% of median): {len(roh)}")
runs = []
for r in sorted(roh, key=lambda r: (int(r['chrom']), r['win_mb'])):
    if runs and runs[-1][0] == r['chrom'] and r['win_mb'] == runs[-1][2] + 1:
        runs[-1][2] = r['win_mb']
    else:
        runs.append([r['chrom'], r['win_mb'], r['win_mb']])
for c, s, e in sorted(runs, key=lambda x: -(x[2] - x[1]))[:20]:
    print(f"   chr{c}:{s}-{e+1} Mb   ({e-s+1} Mb)")
print(f"\nLow-density windows (< 25% of median variants, i.e. coverage/mappability holes): {len(low)}")
for r in sorted(low, key=lambda r: r['total'])[:20]:
    print(f"   chr{r['chrom']}:{r['win_mb']}-{r['win_mb']+1} Mb  total={r['total']} het={r['het']} hom={r['hom_alt']}")
