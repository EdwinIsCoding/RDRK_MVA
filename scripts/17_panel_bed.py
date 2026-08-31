#!/usr/bin/env python3
"""Write BED intervals for the core mitotic panel from the pinned Ensembl GTF.

Built from the GTF rather than the live REST API. The REST lookup overstated
BUB3 by 9.9-fold (158,779 bp against 16,072 bp) because it returns whatever the
current release spans, including readthrough annotation. The GTF is versioned,
matches the gene model the pipeline uses, and is reproducible.
"""
from __future__ import annotations
import csv, gzip, pathlib, re, sys

GTF = "refs/Homo_sapiens.GRCh38.115.gtf.gz"
PANEL = "config/gene_panels/mitotic_extended.tsv"
OUT = "config/gene_panels/mitotic_extended.nochr.bed"

core = {r["symbol"] for r in csv.DictReader(open(PANEL, newline=""), delimiter="\t")
        if r["in_core_panel"] == "yes"}
rows = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "gene":
            continue
        m = re.search(r'gene_name "([^"]+)"', f[8])
        if m and m.group(1) in core:
            rows[m.group(1)] = (f[0], int(f[3]), int(f[4]))

with open(OUT, "w") as fh:
    for sym, (c, s, e) in sorted(rows.items(), key=lambda kv: (kv[1][0], kv[1][1])):
        fh.write(f"{c}\t{max(0, s-5000)}\t{e+5000}\t{sym}\n")
print(f"wrote {len(rows)}/{len(core)} core panel genes to {OUT}")
missing = sorted(core - set(rows))
if missing:
    print(f"absent from the GTF ({len(missing)}): {', '.join(missing[:12])}")
