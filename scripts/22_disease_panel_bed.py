#!/usr/bin/env python3
"""CDS-plus-splice-region BED for the widened disease panel.

Coding sequence with 20 bp of flanking intron, not whole gene bodies. Two
reasons: gnomAD exomes covers coding sequence, so a whole-gene interval would
mostly return nothing; and the widened search is looking for the ordinary
coding variant that an "achievable" answer implies, having already searched
intronic space in the known MVA genes and found nothing.

The 20 bp flank keeps canonical splice sites and the near-splice region in
scope, since those are coding-adjacent and commonly causal.
"""
from __future__ import annotations
import collections, csv, gzip, pathlib, re, sys

GTF = "refs/Homo_sapiens.GRCh38.115.gtf.gz"
PANEL = "config/gene_panels/disease_genes.tsv"
FLANK = 20
PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

want = {r["symbol"] for r in csv.DictReader(open(PANEL, newline=""), delimiter="\t")}
iv = collections.defaultdict(list)
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "CDS" or f[0] not in PRIMARY:
            continue
        m = re.search(r'gene_name "([^"]+)"', f[8])
        if not m or m.group(1) not in want:
            continue
        iv[f[0]].append((max(0, int(f[3]) - 1 - FLANK), int(f[4]) + FLANK, m.group(1)))

merged, n_genes = [], set()
for c, rows in iv.items():
    rows.sort()
    cs, ce, names = None, None, set()
    for s, e, name in rows:
        if cs is None:
            cs, ce, names = s, e, {name}
        elif s <= ce:
            ce = max(ce, e); names.add(name)
        else:
            merged.append((c, cs, ce, ",".join(sorted(names)))); cs, ce, names = s, e, {name}
    if cs is not None:
        merged.append((c, cs, ce, ",".join(sorted(names))))
    n_genes |= {n for _, _, n in rows}

merged.sort(key=lambda r: (len(r[0]), r[0], r[1]))
out = pathlib.Path("config/gene_panels/disease_cds.nochr.bed")
with out.open("w") as fh:
    for c, s, e, n in merged:
        fh.write(f"{c}\t{s}\t{e}\t{n}\n")
bp = sum(e - s for _, s, e, _ in merged)
print(f"{len(merged):,} merged CDS intervals over {len(n_genes):,} genes "
      f"({len(want):,} requested), {bp/1e6:.1f} Mb -> {out}")
missing = want - n_genes
if missing:
    print(f"not found in the GTF: {len(missing)} e.g. {sorted(missing)[:8]}")
