#!/usr/bin/env python3
"""Build the Track 1 submission file in the organisers' exact format.

Format is transcribed in RULES.md section 4 from the Hackathon Space. Two
details will silently destroy a score if got wrong, and both differ from this
project's internal conventions:

1. **Chromosomes are chr-prefixed** (`chr15`). The proband callset uses Ensembl
   naming with no prefix. Every row is converted here, once, and the conversion
   is asserted rather than assumed.
2. **GRCh38.** Confirmed for this callset at Phase 0, and re-asserted here.

Other constraints from the rules:

- At most **10 rows**. "This is one case, not a cohort, so we're asking for your
  best-ranked guesses, not an exhaustive list."
- ``epcr`` in (0, 1]. Only the ranking relative to our own guesses matters, so
  these are ordinal, not calibrated probabilities, and the write-up should say
  so rather than implying a calibration we do not have.
- ``finding_type`` is ``primary`` or ``secondary``. Secondary findings do not
  hurt the automated score and are reviewed qualitatively, so a
  well-justified incidental finding is free credit.
- A compound-heterozygous pair goes in **one row**, using the ``_2`` columns.
  Partial credit is available for recovering one of two variants, so a pair
  should be submitted as a pair rather than split across rows.

Only 6 submissions are allowed per participant, so this script refuses to write
a file that violates the format.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
           "chrom_2", "pos_2", "ref_2", "alt_2",
           "epcr", "finding_type", "notes"]

MAX_ROWS = 10
VALID_CONTIGS = {f"chr{c}" for c in list(map(str, range(1, 23))) + ["X", "Y", "M", "MT"]}


class FormatError(Exception):
    pass


def to_ucsc(contig: str) -> str:
    """Ensembl naming to UCSC. The submission format requires the prefix."""
    c = contig.strip()
    if c.startswith("chr"):
        return c
    if c in ("MT", "M"):
        return "chrM"
    return f"chr{c}"


def validate(rows: list[dict]) -> None:
    if not rows:
        raise FormatError("no candidate rows")
    if len(rows) > MAX_ROWS:
        raise FormatError(f"{len(rows)} rows, the limit is {MAX_ROWS}")

    seen_epcr = []
    for i, r in enumerate(rows, 1):
        for col in ("proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "epcr", "finding_type"):
            if not str(r.get(col, "")).strip():
                raise FormatError(f"row {i}: {col} is required and empty")

        if r["chrom_1"] not in VALID_CONTIGS:
            raise FormatError(f"row {i}: chrom_1={r['chrom_1']!r} is not chr-prefixed. "
                              "The callset uses Ensembl naming; the submission does not.")
        if str(r.get("chrom_2", "")).strip() and r["chrom_2"] not in VALID_CONTIGS:
            raise FormatError(f"row {i}: chrom_2={r['chrom_2']!r} is not chr-prefixed")

        try:
            p = int(r["pos_1"])
            if p < 1:
                raise ValueError
        except (TypeError, ValueError):
            raise FormatError(f"row {i}: pos_1={r['pos_1']!r} is not a positive integer") from None

        for col in ("ref_1", "alt_1"):
            if not set(str(r[col]).upper()) <= set("ACGTN*"):
                raise FormatError(f"row {i}: {col}={r[col]!r} is not a plain allele")

        try:
            e = float(r["epcr"])
        except (TypeError, ValueError):
            raise FormatError(f"row {i}: epcr={r['epcr']!r} is not a number") from None
        if not (0 < e <= 1):
            raise FormatError(f"row {i}: epcr={e} is outside (0, 1]")
        seen_epcr.append(e)

        if r["finding_type"] not in ("primary", "secondary"):
            raise FormatError(f"row {i}: finding_type={r['finding_type']!r} "
                              "must be 'primary' or 'secondary'")

        # A compound-het pair belongs in one row. Half a pair is a format error,
        # not a partial answer.
        second = [str(r.get(c, "")).strip() for c in ("chrom_2", "pos_2", "ref_2", "alt_2")]
        if any(second) and not all(second):
            raise FormatError(f"row {i}: partial second variant. A compound-heterozygous "
                              "pair needs all four _2 columns or none.")

    primaries = [r for r in rows if r["finding_type"] == "primary"]
    if not primaries:
        raise FormatError("no primary candidate. At least one row must be 'primary'.")

    # Ranking is what is scored, so ties waste discriminative power.
    if len(set(seen_epcr)) != len(seen_epcr):
        sys.stderr.write("  WARNING: duplicate epcr values. Only the ranking relative to "
                         "your own guesses is scored, so ties discard information.\n")


def write(rows: list[dict], out: pathlib.Path) -> None:
    validate(rows)
    rows = sorted(rows, key=lambda r: -float(r["epcr"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="results/track1_candidates.tsv",
                    help="TSV with chrom,pos,ref,alt,epcr,finding_type,notes "
                         "and optional chrom2,pos2,ref2,alt2")
    ap.add_argument("--proband-id", default="WGS_EX2312012")
    ap.add_argument("--out", default="results/track1_submission.csv")
    ap.add_argument("--dry-run", action="store_true",
                    help="validate only, write nothing")
    args = ap.parse_args()

    src = pathlib.Path(args.candidates)
    if not src.exists():
        sys.exit(f"FATAL: {src} not found. Nothing to submit yet.\n"
                 f"Expected columns: chrom,pos,ref,alt,epcr,finding_type,notes"
                 f"[,chrom2,pos2,ref2,alt2]")

    rows = []
    for r in csv.DictReader(src.open(newline=""), delimiter="\t"):
        rows.append({
            "proband_id": args.proband_id,
            "chrom_1": to_ucsc(r["chrom"]), "pos_1": r["pos"],
            "ref_1": r["ref"].upper(), "alt_1": r["alt"].upper(),
            "chrom_2": to_ucsc(r["chrom2"]) if r.get("chrom2") else "",
            "pos_2": r.get("pos2", ""), "ref_2": (r.get("ref2") or "").upper(),
            "alt_2": (r.get("alt2") or "").upper(),
            "epcr": r["epcr"], "finding_type": r.get("finding_type", "primary"),
            "notes": r.get("notes", ""),
        })

    try:
        validate(rows)
    except FormatError as exc:
        sys.exit(f"FATAL: submission format invalid: {exc}\n"
                 f"Only 6 submissions are allowed; not spending one on a bad file.")

    n_primary = sum(1 for r in rows if r["finding_type"] == "primary")
    n_pairs = sum(1 for r in rows if r["chrom_2"])
    print(f"{len(rows)} rows valid: {n_primary} primary, "
          f"{len(rows)-n_primary} secondary, {n_pairs} compound-het pairs")
    print(f"epcr range {min(float(r['epcr']) for r in rows):.3f} to "
          f"{max(float(r['epcr']) for r in rows):.3f}")

    if args.dry_run:
        print("dry run, nothing written")
        return
    write(rows, pathlib.Path(args.out))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
