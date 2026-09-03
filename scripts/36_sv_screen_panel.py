#!/usr/bin/env python3
"""Targeted structural variant screen over the known MVA genes.

The gap this closes, and the one it does not
--------------------------------------------
Arm C's genome-wide SV calling was never completed: Delly was launched twice,
failed once on a renamed subcommand, and the compute booking ended before the
second run returned. The report states that plainly and this does not change it.

What this does instead is ask the narrower question that actually bears on the
case: **is there structural variant evidence over the nine known MVA genes?**
That is answerable from the panel BAM we retrieved, without a new booking.

What it can and cannot see
--------------------------
It can see a breakpoint whose reads lie inside the panel regions: discordant
read pairs, and split reads carrying an `SA` supplementary alignment.

It cannot see:

- a breakpoint outside the panel intervals, because those reads were never
  retrieved;
- a balanced rearrangement whose partner is outside the panel;
- copy-number change without a breakpoint in the panel, though panel depth is
  reported alongside and covers that partly.

**The mate-outside-panel artefact is the trap here.** A read whose mate fell
outside the retrieved regions looks discordant when it is merely orphaned by the
extraction. Pairs are therefore only counted when **both** mates are present in
this BAM, and the orphan count is reported separately rather than folded in.
Without that, a panel-scoped BAM manufactures thousands of false discordant
pairs and the screen reports a structural variant over every gene.

Writes results/summaries/arm_c_sv_screen.md.
"""
from __future__ import annotations

import bisect
import collections
import csv
import pathlib
import statistics
import sys

import pysam

BAM = "node_artefacts/WGS_EX2312012.panel.bam"
REGIONS_BED = pathlib.Path("node_artefacts/regions.bed")
PANEL = pathlib.Path("config/gene_panels/mva_known.tsv")
OUT = pathlib.Path("results/summaries/arm_c_sv_screen.md")

MIN_MAPQ = 20
#: A pair is discordant if its insert size is far from the library mode, or its
#: orientation is not the expected forward/reverse. The insert threshold is
#: derived from this library rather than assumed.
INSERT_SD_MULTIPLE = 6
#: Minimum reads supporting a putative breakpoint before it is worth reporting.
MIN_SUPPORT = 3
#: From config/config.yaml analysis.random_seed, so the background sample is
#: reproducible (CLAUDE.md rule 5).
SEED = 20261024


def genes() -> list[tuple[str, str, int, int]]:
    rows = list(csv.DictReader(PANEL.open(newline=""), delimiter="\t"))
    return [(r["symbol"], r["chrom_nochr"], int(r["start"]), int(r["end"]))
            for r in rows]


def load_regions(path: pathlib.Path) -> dict[str, list[tuple[int, int]]]:
    """The intervals the panel BAM was extracted over, for the mate test.

    Asking the BAM whether each mate is present costs a query per read and does
    not finish in any useful time. Whether a mate was retrieved is decided
    instead by whether its position falls inside an extracted interval, which is
    a bisect over a sorted list.
    """
    by_contig: dict[str, list[tuple[int, int]]] = {}
    if not path.exists():
        return by_contig
    with path.open() as fh:
        for line in fh:
            f = line.split()
            if len(f) < 3:
                continue
            by_contig.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    for c in by_contig:
        by_contig[c].sort()
    return by_contig


def in_regions(by_contig: dict[str, list[tuple[int, int]]],
               contig: str | None, pos: int | None) -> bool:
    if contig is None or pos is None:
        return False
    iv = by_contig.get(contig)
    if not iv:
        return False
    i = bisect.bisect_right(iv, (pos, float("inf"))) - 1
    return i >= 0 and iv[i][0] <= pos < iv[i][1]


def library_insert_stats(bam: pysam.AlignmentFile, regions) -> tuple[float, float]:
    """Insert-size centre and spread, measured from this library."""
    sizes = []
    for _, chrom, start, end in regions:
        for rd in bam.fetch(chrom, start, min(end, start + 200_000)):
            if (rd.is_proper_pair and not rd.is_secondary and not rd.is_supplementary
                    and rd.mapping_quality >= MIN_MAPQ and rd.template_length > 0):
                sizes.append(rd.template_length)
            if len(sizes) > 200_000:
                break
    if len(sizes) < 1000:
        return 0.0, 0.0
    med = statistics.median(sizes)
    # Median absolute deviation, robust to the long tail a few chimaeras give.
    mad = statistics.median([abs(s - med) for s in sizes]) * 1.4826
    return med, mad


def scan_region(bam, retrieved, chrom: str, start: int, end: int,
                lo: float, hi: float) -> dict:
    """Breakpoint evidence in one interval. Used for the MVA genes and, with
    identical settings, for the background regions they are judged against."""
    n_reads = n_disc = n_split = n_orphan = 0
    breakpoints: collections.Counter = collections.Counter()
    for rd in bam.fetch(chrom, start, end):
        if rd.is_secondary or rd.is_supplementary or rd.mapping_quality < MIN_MAPQ:
            continue
        n_reads += 1
        if rd.has_tag("SA"):
            n_split += 1
            breakpoints[("split", rd.reference_name, rd.reference_start // 1000)] += 1
        if not rd.is_paired or rd.mate_is_unmapped or rd.is_unmapped:
            continue
        if not in_regions(retrieved, rd.next_reference_name, rd.next_reference_start):
            n_orphan += 1
            continue
        tlen = abs(rd.template_length)
        wrong_orientation = rd.is_reverse == rd.mate_is_reverse
        diff_contig = rd.reference_name != rd.next_reference_name
        if diff_contig or wrong_orientation or (tlen and (tlen > hi or tlen < lo)):
            n_disc += 1
            breakpoints[("pair", rd.reference_name, rd.reference_start // 1000)] += 1
    clusters = [(k, v) for k, v in breakpoints.items() if v >= MIN_SUPPORT]
    kb = max(1.0, (end - start) / 1000)
    return {"reads": n_reads, "disc": n_disc, "split": n_split,
            "orphan": n_orphan, "clusters": clusters,
            "clusters_per_kb": len(clusters) / kb,
            "max_cluster": max((v for _, v in clusters), default=0)}


def background(bam, retrieved, mva_contigs, lo, hi, n_regions=400):
    """The same scan over non-MVA panel regions, to calibrate what "a cluster"
    means here.

    The first run of this screen reported clustered breakpoint evidence over
    seven of the nine MVA genes, which would be an extraordinary finding and is
    in fact the ordinary rate of split-read artefact in short-read data. This
    project has twice closed a lead by calibrating it against a genome-wide
    background (`scripts/06_gap_background_rate.py`), and the same discipline
    applies here: a cluster count means nothing until you know how often one
    occurs by chance.
    """
    import random
    rng = random.Random(SEED)
    pool = [(c, s, e) for c, ivs in retrieved.items() for s, e in ivs
            if (e - s) >= 2000 and c not in mva_contigs]
    if not pool:
        return []
    sample = rng.sample(pool, min(n_regions, len(pool)))
    out = []
    for c, s, e in sample:
        try:
            out.append(scan_region(bam, retrieved, c, s, e, lo, hi))
        except (ValueError, KeyError):
            continue
    return out


def main() -> None:
    if not pathlib.Path(BAM).exists():
        sys.exit(f"FATAL: {BAM} not present. Nothing to screen.")
    bam = pysam.AlignmentFile(BAM, "rb")
    regions = genes()
    retrieved = load_regions(REGIONS_BED)
    if not retrieved:
        sys.exit(f"FATAL: {REGIONS_BED} is absent. Without the extraction "
                 f"intervals the mate test cannot run, and every mate that was "
                 f"simply never retrieved would be counted as a discordant "
                 f"pair. That would manufacture a structural variant over every "
                 f"gene, so this refuses to run rather than report one.")

    med, mad = library_insert_stats(bam, regions)
    if med == 0:
        sys.exit("FATAL: could not measure insert size; the BAM looks wrong.")
    hi = med + INSERT_SD_MULTIPLE * mad
    lo = max(0, med - INSERT_SD_MULTIPLE * mad)

    L: list[str] = []
    w = L.append
    w("# Arm C: structural variant screen over the known MVA genes\n")
    w("Generated by `scripts/36_sv_screen_panel.py` from the panel BAM.\n")
    w("**This is not the genome-wide SV calling that Arm C never completed.** "
      "It is the narrower question that bears on the case: is there breakpoint "
      "evidence over the nine known MVA genes? Read the limits at the end before "
      "the result.\n")
    w("## Library, measured rather than assumed\n")
    w(f"| median insert size | {med:.0f} bp |\n|---|---|\n"
      f"| robust spread (MAD-derived) | {mad:.0f} bp |\n"
      f"| discordant if insert outside | {lo:.0f} to {hi:.0f} bp |\n"
      f"| minimum reads to report a breakpoint | {MIN_SUPPORT} |\n")

    mva_contigs = {c for _, c, _, _ in regions}
    bg = background(bam, retrieved, mva_contigs, lo, hi)
    bg_rates = sorted(r["clusters_per_kb"] for r in bg)
    bg_max = sorted(r["max_cluster"] for r in bg)

    def pct(sorted_vals, x):
        import bisect as _b
        if not sorted_vals:
            return float("nan")
        return 100.0 * _b.bisect_right(sorted_vals, x) / len(sorted_vals)

    w("## Calibration first\n")
    w("A cluster count means nothing until you know how often one arises by "
      f"chance. The identical scan was run over **{len(bg)}** randomly sampled non-MVA "
      "panel regions, seeded from `config/config.yaml`.\n")
    if bg_rates:
        w("| background, clusters per kb | value |")
        w("|---|---:|")
        w(f"| median | {statistics.median(bg_rates):.3f} |")
        w(f"| 90th percentile | {bg_rates[int(0.90 * (len(bg_rates) - 1))]:.3f} |")
        w(f"| 99th percentile | {bg_rates[int(0.99 * (len(bg_rates) - 1))]:.3f} |")
        w(f"| regions with at least one cluster | "
          f"{100.0 * sum(1 for r in bg_rates if r > 0) / len(bg_rates):.0f}% |")
        w("")
        w(f"**{100.0 * sum(1 for r in bg_rates if r > 0) / len(bg_rates):.0f}% of "
          f"ordinary panel regions contain at least one cluster by this "
          f"definition.** Clustered split reads are the background, not the "
          f"signal, which is what the first run of this screen got wrong.\n")

    w("## Result by gene\n")
    thresh = bg_rates[int(0.99 * (len(bg_rates) - 1))] if bg_rates else 0.0
    max_thresh = bg_max[int(0.99 * (len(bg_max) - 1))] if bg_max else 0
    w(f"Judged against the background 99th percentile: **{thresh:.3f} clusters "
      f"per kb**, and a largest single cluster of **{max_thresh} reads**. "
      f"Percentile ranks are not quoted per gene because the background is "
      f"zero-inflated, {100.0 * sum(1 for r in bg_rates if r == 0) / len(bg_rates):.0f}% "
      f"of regions having no cluster at all, so most genes would tie at the same "
      f"rank and the number would imply a precision it does not have.\n")
    w("| Gene | Region | Reads | Discordant pairs | Split reads | Orphaned mates "
      "(excluded) | Clusters/kb | Largest cluster | Verdict |")
    w("|---|---|---:|---:|---:|---:|---:|---:|---|")

    outliers = []
    for sym, chrom, start_, end_ in regions:
        r = scan_region(bam, retrieved, chrom, start_, end_, lo, hi)
        unusual = (r["clusters_per_kb"] > thresh) or (r["max_cluster"] > max_thresh)
        if unusual:
            outliers.append((sym, r, pct(bg_rates, r["clusters_per_kb"]),
                             pct(bg_max, r["max_cluster"])))
        verdict = "**exceeds background**" if unusual else "within background"
        w(f"| {sym} | {chrom}:{start_}-{end_} | {r['reads']:,} | {r['disc']} | "
          f"{r['split']} | {r['orphan']:,} | {r['clusters_per_kb']:.3f} | "
          f"{r['max_cluster']} | {verdict} |")
    w("")

    w("## Reading this\n")
    if not outliers:
        w("**No known MVA gene carries breakpoint evidence beyond the background "
          "rate.** Every gene sits inside the distribution measured over ordinary "
          "panel regions on the same reads with the same settings.\n")
        w("Combined with the panel depth already reported, 42 to 51x with no gene "
          "under-covered, this is a **reported negative for structural variation "
          "over the known MVA genes**, at the sensitivity a panel-scoped BAM "
          "allows.\n")
        w("**The uncalibrated version of this screen said the opposite.** It "
          "flagged clustered evidence over seven of nine genes, which would have "
          "been a remarkable finding and was the ordinary rate of split-read "
          "artefact. The calibration is the analysis; the raw counts are not.\n")
    else:
        w("**Genes exceeding the 99th percentile of the background:**\n")
        for sym, r, p_rate, p_max in outliers:
            w(f"- `{sym}`: {r['clusters_per_kb']:.3f} clusters/kb "
              f"({p_rate:.0f}th percentile), largest cluster {r['max_cluster']} "
              f"reads ({p_max:.0f}th)")
            for (kind, contig, kb), n in sorted(r["clusters"], key=lambda x: -x[1])[:5]:
                w(f"  - {n} {kind} reads near {contig}:{kb * 1000:,}")
        w("\nAn outlier is a hypothesis, not a call. Short-read breakpoint "
          "evidence in repetitive or paralogous sequence is unreliable and there "
          "is no matched control here, so this warrants inspection rather than "
          "belief.\n")

    w("## Limits, which are the point\n")
    w("- **Not genome-wide.** Only the nine known MVA genes, because only they "
      "are in the retrieved BAM. A causal SV elsewhere is invisible here and "
      "Arm C's genome-wide gap remains open.\n")
    w("- **Mate-outside-panel pairs are excluded, not counted.** The orphan "
      "column above shows how many. Counting them as discordant would "
      "manufacture a structural variant over every gene, which is the failure "
      "mode a panel-scoped BAM invites.\n")
    w("- **No matched control.** Discordant-pair rates vary by region and there "
      "is no second sample to normalise against, so only clustering is used as "
      "evidence, never a raw rate.\n")
    w("- **A balanced rearrangement with a partner outside the panel is "
      "undetectable**, and the parental recurrent-miscarriage history is exactly "
      "the signal such a rearrangement would give. That hypothesis remains open "
      "and untested.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
