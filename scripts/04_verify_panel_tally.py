#!/usr/bin/env python3
"""Cross-check indexed region queries against a full streaming pass.

The shipped .tbi lacks the count metadata that newer bcftools writes, and
warns that it predates the data file. Before any 'coverage hole' claim is made
about a panel gene, the indexed tally must be shown to agree with an
index-independent stream. A disagreement means the index is untrustworthy and
every region-scoped result so far is void.
"""
import collections, subprocess, sys

VCF = "data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
BED = "config/gene_panels/mva_known.nochr.bed"

intervals = collections.defaultdict(list)
for line in open(BED):
    c, s, e, name = line.split()[:4]
    intervals[c].append((int(s), int(e), name))

stream = collections.Counter()
dp = collections.defaultdict(list)
p = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%POS\t[%DP]\n", VCF],
                     stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in p.stdout:
    c, pos, d = line.rstrip("\n").split("\t")
    iv = intervals.get(c)
    if not iv:
        continue
    pos = int(pos)
    for s, e, name in iv:
        if s <= pos <= e:
            stream[name] += 1
            if d.isdigit():
                dp[name].append(int(d))
p.wait()

print(f"{'gene':10} {'streamed':>9} {'indexed':>9} {'agree':>7} {'med_DP':>7} {'min_DP':>7}")
bad = 0
for line in open(BED):
    c, s, e, name = line.split()[:4]
    region = f"{c}:{s}-{e}"
    out = subprocess.run(["bcftools", "query", "-r", region, "-f", ".\n", VCF],
                         capture_output=True, text=True)
    idx = len(out.stdout.splitlines())
    ok = "YES" if idx == stream[name] else "NO"
    if ok == "NO":
        bad += 1
    d = sorted(dp[name])
    print(f"{name:10} {stream[name]:>9} {idx:>9} {ok:>7} "
          f"{(d[len(d)//2] if d else '-'):>7} {(d[0] if d else '-'):>7}")

print()
if bad:
    print(f"*** INDEX UNRELIABLE: {bad} genes disagree. All region-scoped results are void. ***")
else:
    print("Index agrees with stream on every panel interval.")
