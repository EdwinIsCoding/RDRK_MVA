#!/usr/bin/env python3
"""Background rate of empty 10 kb windows and of >=20 kb variant-free runs.

Calibrates the panel-gene gaps found by 05_panel_depth_profile.py. Without this
the gaps are uninterpretable: runs of homozygosity of this length are common in
every genome, so an uncalibrated 'gap in BUB1B' claim is not evidence.
Restricted to chromosomes 5, 11 and 15 (which carry TRIP13, CEP57 and BUB1B)
and to windows with mappable flanking sequence.
"""
import collections, subprocess

VCF = "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
WIN = 10_000
CHROMS = ["5", "11", "15"]

for chrom in CHROMS:
    n = collections.Counter()
    p = subprocess.Popen(["bcftools", "query", "-r", chrom, "-f", "%POS\t%FILTER\n", VCF],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        pos, filt = line.rstrip("\n").split("\t")
        if filt in ("PASS", "."):
            n[int(pos) // WIN] += 1
    p.wait()

    if not n:
        continue
    lo, hi = min(n), max(n)
    windows = list(range(lo, hi + 1))
    # Exclude sparse regions (centromere, telomere) by requiring the local
    # 1 Mb neighbourhood to be reasonably dense, so we count gaps inside
    # otherwise well-covered sequence only.
    dense = []
    for w in windows:
        nb = sum(n[x] for x in range(w - 50, w + 51))
        if nb >= 800:
            dense.append(w)
    empty = [w for w in dense if n[w] == 0]

    runs, cur = [], []
    for w in dense:
        if n[w] == 0:
            cur.append(w)
        else:
            if cur: runs.append(len(cur))
            cur = []
    if cur: runs.append(len(cur))
    ge2 = sum(1 for r in runs if r >= 2)

    print(f"chr{chrom}: {len(dense):,} evaluable 10 kb windows | "
          f"empty {len(empty):,} ({100*len(empty)/len(dense):.2f}%) | "
          f"variant-free runs >=20 kb: {ge2} | longest run: {max(runs) * 10 if runs else 0} kb")
