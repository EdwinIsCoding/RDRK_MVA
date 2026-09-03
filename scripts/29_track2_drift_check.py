#!/usr/bin/env python3
"""Has the evidence moved since the Track 2 report was written?

Why this exists
---------------
The Track 2 report quotes counts from two live sources. ClinicalTrials.gov and
ChEMBL both change without notice, and the report's central premise is a set of
four zeros: no registered trial has this proband's disease as a condition and
none has their tumour as a prevention target. A judge running the pipeline in
November could see different numbers from the ones the report states, and would
have no way to tell whether we were wrong or the world moved.

So the counts are pinned here with the date they were taken, and this script
re-queries and reports the difference. It changes nothing on its own.

If a zero has become non-zero, that is not a number to patch. It means a
chemoprevention trial now exists for this disease, section 4.1 of the report is
no longer true, and the section has to be rewritten rather than adjusted.

Usage
-----
    python scripts/29_track2_drift_check.py            # compare against the pin
    python scripts/29_track2_drift_check.py --update   # re-pin to today
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import pathlib
import sys

sys.path.insert(0, "src")

from mva.track2.chemoprevention import (
    MESH_HEREDITARY_NEOPLASTIC,
    count_trials,
)

PIN = pathlib.Path("config/track2_evidence_pin.json")
CACHE = pathlib.Path("results/track2/cache_drift")

QUERIES = {
    "mva_any_trial":
        'AREA[ConditionSearch]"Mosaic Variegated Aneuploidy"',
    "bub1b_any_trial":
        'AREA[ConditionSearch]"BUB1B" OR AREA[InterventionName]"BUB1B"',
    "rhabdomyosarcoma_prevention":
        'AREA[ConditionSearch]"Rhabdomyosarcoma" AND AREA[DesignPrimaryPurpose]PREVENTION',
    "rhabdomyosarcoma_prevention_drug":
        'AREA[ConditionSearch]"Rhabdomyosarcoma" AND AREA[DesignPrimaryPurpose]PREVENTION '
        'AND AREA[InterventionType]DRUG',
    "hereditary_prevention_drug":
        f'AREA[ConditionSearch]"{MESH_HEREDITARY_NEOPLASTIC}" '
        f'AND AREA[DesignPrimaryPurpose]PREVENTION AND AREA[InterventionType]DRUG',
}

#: Counts whose value is load-bearing for a claim in the report, and what breaks
#: if they move. A drift in one of these is a rewrite, not an edit.
LOAD_BEARING = {
    "mva_any_trial": "Track 2 report section 4.1: the evidence base is empty",
    "rhabdomyosarcoma_prevention": "Track 2 report section 4.1: the evidence base is empty",
}


def observe() -> dict[str, int]:
    CACHE.mkdir(parents=True, exist_ok=True)
    return {k: count_trials(q, CACHE) for k, q in QUERIES.items()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true",
                    help="re-pin the counts to today rather than comparing")
    args = ap.parse_args()

    now = observe()

    if args.update or not PIN.exists():
        PIN.parent.mkdir(parents=True, exist_ok=True)
        PIN.write_text(json.dumps({
            "taken": dt.date.today().isoformat(),
            "source": "ClinicalTrials.gov API v2",
            "queries": QUERIES,
            "counts": now,
            "note": ("Pinned so that drift is detectable. See "
                     "scripts/29_track2_drift_check.py. A load-bearing count "
                     "that moves is a rewrite of the claim it supports, not an "
                     "edit to the number."),
        }, indent=2) + "\n")
        print(f"pinned {len(now)} counts to {PIN}")
        for k, v in sorted(now.items()):
            print(f"  {k}: {v}")
        return 0

    pin = json.loads(PIN.read_text())
    was, taken = pin["counts"], pin["taken"]
    print(f"Comparing against the pin taken {taken}.\n")

    drifted = []
    for k in sorted(QUERIES):
        a, b = was.get(k), now.get(k)
        if a == b:
            print(f"  unchanged  {k}: {b}")
        else:
            drifted.append((k, a, b))
            print(f"  CHANGED    {k}: {a} -> {b}")

    if not drifted:
        print("\nNo drift. Every count in the Track 2 report still holds.")
        return 0

    print(f"\n{len(drifted)} count(s) have moved since {taken}.")
    breaking = [(k, a, b) for k, a, b in drifted if k in LOAD_BEARING]
    for k, a, b in breaking:
        print(f"\n  LOAD-BEARING: {k} went from {a} to {b}.")
        print(f"  This supports: {LOAD_BEARING[k]}")
        if a == 0 and b:
            print("  A claim that the evidence base is empty is no longer true. "
                  "Rewrite the section; do not adjust the number.")
    if not breaking:
        print("\nNone of the moved counts is load-bearing. Update the report's "
              "figures and re-pin with --update.")
    return 1 if breaking else 0


if __name__ == "__main__":
    raise SystemExit(main())
