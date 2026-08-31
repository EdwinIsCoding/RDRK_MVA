#!/usr/bin/env python3
"""Regenerate the input-file table in PROVENANCE.md from recon artefacts.

Keeps checksums, sizes and the manifest in one place and reproducible, so the
provenance table is never hand-edited and never drifts from the data.
"""
import datetime, pathlib, subprocess

recon = pathlib.Path("results/recon")
sha = {}
for line in (recon / "sha256.txt").read_text().splitlines():
    if not line.strip() or line.startswith("SKIPPED"):
        continue
    h, p = line.split(None, 1)
    sha[p.strip()] = h

rows = []
for line in (recon / "manifest.tsv").read_text().splitlines()[1:]:
    if line.strip():
        p, b = line.split("\t")
        rows.append((p, int(b), sha.get(p, "PENDING")))

total = sum(b for _, b, _ in rows)
out = ["| File | Bytes | SHA256 |", "|---|---:|---|"]
for p, b, h in rows:
    out.append(f"| `{pathlib.Path(p).name}` | {b:,} | `{h}` |")
out.append(f"| **{len(rows)} files** | **{total:,}** | |")

table = "\n".join(out)
marker_a, marker_b = "<!-- BEGIN INPUT TABLE -->", "<!-- END INPUT TABLE -->"
prov = pathlib.Path("PROVENANCE.md")
text = prov.read_text()
pre, rest = text.split(marker_a, 1)
_, post = rest.split(marker_b, 1)
stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
prov.write_text(f"{pre}{marker_a}\n\n*Regenerated {stamp} by `scripts/07_write_provenance.py`.*\n\n{table}\n\n{marker_b}{post}")
missing = sum(1 for _, _, h in rows if h == "PENDING")
print(f"PROVENANCE.md updated: {len(rows)} files, {missing} checksums pending")
