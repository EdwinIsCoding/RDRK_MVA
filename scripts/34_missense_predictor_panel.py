#!/usr/bin/env python3
"""Consult the full in silico predictor panel for the two BUB1B alleles.

Why this exists
---------------
The Track 1 report applied ACMG criterion PP3 on the strength of "two concordant
in silico predictors", SIFT and PolyPhen-2, and listed AlphaMissense, CADD and
REVEL among its limitations as "not available". They are available, for a single
variant, through a public aggregator of dbNSFP. Leaving a stated limitation
unchecked when one API call closes it is not a limitation, it is an omission.

The answer changes the claim in both directions at once, which is why it was
worth asking: far more predictors exist than we cited, and they do **not** all
agree.

What is transmitted, and the decision behind it
-----------------------------------------------
A genomic coordinate and its alleles. No genotype, no sample identifier, nothing
linking the query to the proband. Both coordinates are already published in
`submission/track1_submission.csv` and in the Track 1 report, under CC BY 4.0 and
in a public repository, so the query reveals nothing a reader of the submission
could not already look up themselves.

This is nonetheless a step beyond the position taken at Phase 0, where gnomAD was
queried by gene interval **specifically** to avoid sending proband variants to a
third party. It is recorded as a decision in `ETHICS.md` section 2a rather than
taken quietly. The offline alternative, a local dbNSFP installation, is tens of
gigabytes for two lookups.

On thresholds
-------------
Where a predictor supplies its own categorical call, that call is reported. Where
it supplies only a score, the score is reported and **no threshold is applied**,
because the thresholds in common use would be recited from memory and CLAUDE.md
rule 2 exists to stop exactly that. The disagreement between predictors is
visible without needing a single threshold, which is the point.

Writes results/summaries/missense_predictor_panel.md.
"""
from __future__ import annotations

import datetime as dt
import json
import pathlib
import urllib.parse
import urllib.request

OUT = pathlib.Path("results/summaries/missense_predictor_panel.md")
CACHE = pathlib.Path("results/track2/cache_predictors")
API = "https://myvariant.info/v1/variant"

#: Both alleles, as published in the submission.
ALLELES = [
    ("Allele 1, p.Leu737Ter, nonsense", "chr15:g.40209701T>G"),
    ("Allele 2, p.Asn1002Lys, missense", "chr15:g.40220612T>G"),
]

#: Predictors that supply a categorical call, and how to read it.
CATEGORICAL = [
    ("AlphaMissense", "alphamissense", {"P": "pathogenic", "B": "benign", "A": "ambiguous"}),
    ("ClinPred", "clinpred", {"D": "damaging", "T": "tolerated"}),
    ("DEOGEN2", "deogen2", {"D": "damaging", "T": "tolerated"}),
    ("FATHMM", "fathmm", {"D": "damaging", "T": "tolerated"}),
    ("LIST-S2", "list-s2", {"D": "damaging", "T": "tolerated"}),
    ("M-CAP", "m-cap", {"D": "damaging", "T": "tolerated"}),
    ("MetaLR", "metalr", {"D": "damaging", "T": "tolerated"}),
    ("MetaRNN", "metarnn", {"D": "damaging", "T": "tolerated"}),
    ("MetaSVM", "metasvm", {"D": "damaging", "T": "tolerated"}),
    ("MutationAssessor", "mutationassessor",
     {"H": "high", "M": "medium", "L": "low", "N": "neutral"}),
    ("MutationTaster", "mutationtaster",
     {"A": "disease causing automatic", "D": "disease causing",
      "N": "polymorphism", "P": "polymorphism automatic"}),
    ("PROVEAN", "provean", {"D": "damaging", "N": "neutral"}),
    # PolyPhen-2 nests its call under hdiv/hvar rather than exposing "pred" at
    # the top level, so a generic reader misses it entirely. It was missing from
    # the first version of this panel, which reported 15 predictors instead of
    # 16 and, worse, omitted one of the two the report had originally cited.
    ("PolyPhen-2 HDIV", "polyphen2.hdiv", {"D": "probably damaging",
                                           "P": "possibly damaging",
                                           "B": "benign"}),
    ("PrimateAI", "primateai", {"D": "damaging", "T": "tolerated"}),
    ("SIFT", "sift", {"D": "deleterious", "T": "tolerated"}),
    ("SIFT4G", "sift4g", {"D": "deleterious", "T": "tolerated"}),
]

#: Scores reported without applying any threshold of ours.
SCORES = [
    ("REVEL", ("revel", "score")),
    ("CADD phred", ("cadd", "phred")),
    ("VEST4", ("vest4", "score")),
    ("MVP", ("mvp", "score")),
    ("MPC", ("mpc", "score")),
    ("GERP++ RS", ("gerp++", "rs")),
]

#: Calls that count as evidence for damage, for the tally.
DAMAGING = {"pathogenic", "damaging", "deleterious", "disease causing",
            "disease causing automatic", "high", "probably damaging"}
BENIGN = {"benign", "tolerated", "neutral", "polymorphism",
          "polymorphism automatic"}
#: "possibly damaging" is deliberately in neither set: it is PolyPhen-2 hedging,
#: and forcing it into one column would overstate whichever column it joined.


def fetch(hgvs: str) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    key = CACHE / (hgvs.replace(":", "_").replace(">", "_") + ".json")
    if key.exists():
        return json.loads(key.read_text())
    url = f"{API}/{urllib.parse.quote(hgvs)}?assembly=hg38"
    req = urllib.request.Request(url, headers={"User-Agent": "mva-hackathon-2026"})
    with urllib.request.urlopen(req, timeout=90) as fh:
        d = json.load(fh)
    key.write_text(json.dumps(d))
    return d


def first(v):
    """dbNSFP returns per-transcript lists. Collapse, keeping disagreement visible."""
    if isinstance(v, list):
        uniq = sorted({str(x) for x in v if x is not None})
        return "/".join(uniq) if uniq else None
    return v


def main() -> None:
    L: list[str] = []
    w = L.append
    w("# In silico predictor panel for the two BUB1B alleles\n")
    w(f"Generated by `scripts/34_missense_predictor_panel.py` on "
      f"{dt.date.today().isoformat()}. Source: MyVariant.info, serving dbNSFP.\n")
    w("The Track 1 report applied ACMG PP3 on two predictors and listed "
      "AlphaMissense, CADD and REVEL as unavailable. They are available. This is "
      "what consulting them shows.\n")
    w("Where a predictor supplies its own categorical call, that call is shown. "
      "Where it supplies only a score, the score is shown and **no threshold of "
      "ours is applied**, because the thresholds in common use would be recited "
      "from memory (`CLAUDE.md` rule 2). The disagreement is visible without "
      "one.\n")

    for label, hgvs in ALLELES:
        try:
            d = fetch(hgvs).get("dbnsfp", {})
        except Exception as exc:
            w(f"## {label}\n\nLookup failed: {exc}\n")
            continue
        w(f"## {label}\n")
        w(f"`{hgvs}` (GRCh38)\n")
        if not d:
            w("No dbNSFP record. Expected for a nonsense variant, which the "
              "missense predictors do not score.\n")

        rows, dmg, ben, other = [], 0, 0, 0
        for name, key, legend in CATEGORICAL:
            blk = d
            for part in key.split("."):
                blk = blk.get(part) if isinstance(blk, dict) else None
            if not isinstance(blk, dict):
                continue
            pred = first(blk.get("pred"))
            if pred is None:
                continue
            words = sorted({legend.get(p, p) for p in str(pred).split("/")})
            meaning = "/".join(words)
            score = first(blk.get("score"))
            rows.append((name, meaning, score))
            if any(x in DAMAGING for x in words):
                dmg += 1
            elif any(x in BENIGN for x in words):
                ben += 1
            else:
                other += 1

        if rows:
            w("### Predictors that make a call\n")
            w("| Predictor | Call | Score |")
            w("|---|---|---|")
            for name, meaning, score in rows:
                w(f"| {name} | **{meaning}** | {score if score is not None else 'n/a'} |")
            w("")
            def v(n):
                return "calls" if n == 1 else "call"
            w(f"**{dmg} {v(dmg)} it damaging, {ben} {v(ben)} it tolerated or "
              f"benign, and {other} {'sits' if other == 1 else 'sit'} in "
              f"between.**\n")

        srows = []
        for name, path in SCORES:
            blk = d
            for p in path:
                blk = blk.get(p) if isinstance(blk, dict) else None
            if blk is not None:
                srows.append((name, first(blk)))
        if srows:
            w("### Scores, with no threshold applied\n")
            w("| Predictor | Score |")
            w("|---|---|")
            for name, val in srows:
                w(f"| {name} | {val} |")
            w("")

    w("## What this does to the PP3 claim\n")
    w("The Track 1 report said PP3 rested on **two concordant** in silico "
      "predictors. Consulting the full panel corrects that in two directions at "
      "once, and both belong in the report.\n")
    w("**More evidence exists than we cited.** AlphaMissense, the strongest "
      "single modern missense predictor, calls the variant pathogenic, and a "
      "clear majority of the panel calls it damaging.\n")
    w("**The predictors are not concordant.** Several tolerate it, including two "
      "meta-predictors that combine many of the others, and REVEL sits in the "
      "middle of its range rather than near either end. Describing the "
      "computational evidence as concordant was true only of the two predictors "
      "we happened to have run.\n")
    w("PP3 is retained at supporting weight, which is where it already was, and "
      "the disagreement is now stated wherever it is applied. The honest summary "
      "is that computational evidence favours a damaging effect and does not "
      "establish one, which is the same conclusion the report reached by a "
      "narrower route.\n")
    w("None of this changes the call. Allele 2 is not the allele carrying the "
      "diagnosis; allele 1 is a ClinVar-pathogenic nonsense variant in a gene "
      "where loss of function is the established mechanism. Allele 2's role is "
      "to be the second hit, and the argument for that rests on the recessive "
      "mechanism and on its rarity, not on any predictor score.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
