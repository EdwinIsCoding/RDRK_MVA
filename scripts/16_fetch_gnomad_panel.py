#!/usr/bin/env python3
"""Fetch gnomAD v4.1 allele frequencies for the panel regions, by remote slicing.

Why this is possible without a 60 GB download
---------------------------------------------
gnomAD publishes per-chromosome sites VCFs with tabix indexes over HTTPS.
``bcftools`` can range-request a region directly, so pulling every gnomAD record
across the mitotic panel costs a few hundred megabytes of transfer instead of
the full release. Slicing BUB1B returned 9,149 records in under three seconds.

Why it does not breach the data governance
------------------------------------------
The only thing transmitted is a gene interval: "send me chr15:40161023-40221136".
That reveals an interest in BUB1B, which is public knowledge and is written in
this repository's gene panel. **No proband coordinate, allele or genotype leaves
the machine.** The join against the proband happens locally, afterwards.

This is a deliberate distinction from the gnomAD GraphQL API, which would
require sending the proband's own variant coordinates to a third party. That is
not done, and must not be.

What is fetched
---------------
``AF``            overall allele frequency
``AF_grpmax``     maximum frequency across genetic ancestry groups, the field
                  that replaced ``AF_popmax`` in v4. This is the one to filter
                  on: a variant rare overall but common in one ancestry group is
                  not a plausible cause of a severe recessive condition.
``nhomalt``       count of homozygous individuals observed. For a severe
                  recessive paediatric phenotype this is the sharpest single
                  filter available: a variant with homozygotes in a population
                  reference cannot be causal in the homozygous state.
``AN``            allele number, so coverage-poor sites can be distinguished
                  from genuinely absent ones.

Both exomes and genomes are fetched. A panel gene's introns are covered by the
genomes callset only, and the introns are where Arm B is looking.

Output: refs/gnomad_panel/{exomes,genomes}.chr{N}.tsv.gz plus a merged
        refs/gnomad_panel/panel_af.tsv.gz keyed by chrom:pos:ref:alt
"""

from __future__ import annotations

import argparse
import csv
import gzip
import pathlib
import subprocess
import sys
import collections

BASE = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
FMT = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_grpmax\t%INFO/nhomalt\t%INFO/AN\n"


def panel_regions(panel_tsv: pathlib.Path, gtf_bed: pathlib.Path,
                  flank: int = 2000) -> dict[str, list[tuple[int, int, str]]]:
    """Gene intervals grouped by chromosome, in chr-prefixed form.

    gnomAD uses UCSC naming; the proband callset does not. The rename happens
    here, once, rather than being rediscovered at every join.
    """
    by_chrom: collections.defaultdict[str, list] = collections.defaultdict(list)
    for line in gtf_bed.read_text().splitlines():
        if not line.strip():
            continue
        c, s, e, name = line.split("\t")[:4]
        by_chrom[f"chr{c}" if not c.startswith("chr") else c].append(
            (max(1, int(s) - flank), int(e) + flank, name))
    return dict(by_chrom)


def slice_gnomad(kind: str, chrom: str, regions: list[tuple[int, int, str]],
                 out: pathlib.Path, timeout: int = 1800) -> int:
    """Range-request one chromosome's gnomAD file for the given intervals."""
    url = f"{BASE}/{kind}/gnomad.{kind}.v4.1.sites.{chrom}.vcf.bgz"
    region_args = ",".join(f"{chrom}:{s}-{e}" for s, e, _ in regions)
    cmd = ["bcftools", "query", "-r", region_args, "-f", FMT, url]
    n = 0
    with gzip.open(out, "wt") as fh:
        fh.write("chrom\tpos\tref\talt\taf\taf_grpmax\tnhomalt\tan\n")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, bufsize=1 << 20)
        assert proc.stdout is not None
        for line in proc.stdout:
            fh.write(line)
            n += 1
        proc.wait()
        if proc.returncode != 0:
            err = (proc.stderr.read() if proc.stderr else "").strip()
            raise RuntimeError(f"{kind} {chrom} failed: {err[:400]}")
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default="config/gene_panels/mitotic_extended.nochr.bed")
    ap.add_argument("--panel", default="config/gene_panels/mitotic_extended.tsv")
    ap.add_argument("--outdir", default="refs/gnomad_panel")
    ap.add_argument("--kinds", nargs="+", default=["exomes", "genomes"])
    ap.add_argument("--chroms", nargs="*", default=None,
                    help="restrict to these chromosomes, e.g. chr15 chr11")
    args = ap.parse_args()

    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    bed = pathlib.Path(args.bed)
    if not bed.exists():
        sys.exit(f"FATAL: {bed} not found. Build it with scripts/17_panel_bed.py")

    regions = panel_regions(pathlib.Path(args.panel), bed)
    chroms = args.chroms or sorted(regions, key=lambda c: (len(c), c))

    total = 0
    for kind in args.kinds:
        for chrom in chroms:
            if chrom not in regions:
                continue
            dest = out / f"{kind}.{chrom}.tsv.gz"
            if dest.exists() and dest.stat().st_size > 100:
                sys.stderr.write(f"  {kind} {chrom}: cached\n")
                continue
            try:
                n = slice_gnomad(kind, chrom, regions[chrom], dest)
            except Exception as exc:  # noqa: BLE001
                sys.stderr.write(f"  {kind} {chrom}: FAILED {exc}\n")
                dest.unlink(missing_ok=True)
                continue
            total += n
            sys.stderr.write(f"  {kind} {chrom}: {n:,} records "
                             f"({len(regions[chrom])} genes)\n")

    # Merge into one lookup keyed on the proband's own naming convention
    # (no chr prefix), so the join needs no rename at query time.
    merged = out / "panel_af.tsv.gz"
    seen: dict[tuple, tuple] = {}
    for f in sorted(out.glob("*.chr*.tsv.gz")):
        with gzip.open(f, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = (row["chrom"].removeprefix("chr"), row["pos"], row["ref"], row["alt"])
                af = row["af_grpmax"] if row["af_grpmax"] not in (".", "") else row["af"]
                prev = seen.get(key)
                # Exomes and genomes both cover coding sequence. Keep the higher
                # frequency: under-reporting rarity would falsely promote a
                # variant, which is the more dangerous error here.
                try:
                    afv = float(af) if af not in (".", "") else 0.0
                except ValueError:
                    afv = 0.0
                nh = row["nhomalt"] if row["nhomalt"] not in (".", "") else "0"
                if prev is None or afv > prev[0]:
                    seen[key] = (afv, int(nh) if nh.isdigit() else 0)
                elif prev is not None:
                    seen[key] = (prev[0], max(prev[1], int(nh) if nh.isdigit() else 0))

    with gzip.open(merged, "wt") as fh:
        fh.write("key\taf_grpmax\tnhomalt\n")
        for (c, p, r, a), (af, nh) in seen.items():
            fh.write(f"{c}:{p}:{r}:{a}\t{af:.8g}\t{nh}\n")

    print(f"\nfetched {total:,} gnomAD records")
    print(f"merged to {merged}: {len(seen):,} unique variants")


if __name__ == "__main__":
    main()
