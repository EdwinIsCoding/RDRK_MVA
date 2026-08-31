#!/usr/bin/env python3
"""Join gnomAD v4.1 constraint metrics onto the mitotic panel.

LOEUF (``lof.oe_ci.upper``, the upper bound of the 90% confidence interval on
the observed/expected ratio for loss-of-function variants) is the field to use,
not pLI. pLI saturates: it cannot distinguish a moderately constrained gene from
an extremely constrained one, and it is unreliable for short genes where the
expected count is small. LOEUF is continuous and carries its own uncertainty.
Both are emitted, because pLI is what most readers recognise.

Lower LOEUF means more constrained. gnomAD's own guidance is that LOEUF < 0.6
marks a gene intolerant of heterozygous loss of function.

Important caveat, recorded here because it bears directly on this project:
**constraint is a weak prior for an autosomal recessive disorder.** LOEUF
measures selection against heterozygous loss of function. A gene causing a
recessive condition can be entirely unconstrained in heterozygotes and still be
the answer, which is exactly the situation for many recessive disease genes.
BUB1B itself is a case in point. Constraint is therefore recorded and reported,
but it must not be used as a hard filter, and it should carry little weight in
the additive score of plan section 6.7.

Transcript selection: MANE Select where available, otherwise the canonical
transcript, otherwise the transcript with the largest expected LoF count.

Usage:
    python3 scripts/12_join_constraint.py \
        [--constraint refs/gnomad.v4.1.constraint_metrics.tsv]
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys

CONSTRAINT_DEFAULT = "refs/gnomad.v4.1.constraint_metrics.tsv"
GNOMAD_VERSION = "gnomAD v4.1"


def load_constraint(path: pathlib.Path) -> dict[str, dict]:
    """One record per gene symbol, choosing the best transcript."""
    best: dict[str, dict] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue

            def num(key: str) -> float | None:
                v = (row.get(key) or "").strip()
                if v in ("", "NA", "NaN", "nan"):
                    return None
                try:
                    return float(v)
                except ValueError:
                    return None

            rank = (
                2 if (row.get("mane_select") or "").lower() == "true" else
                1 if (row.get("canonical") or "").lower() == "true" else 0
            )
            rec = {
                "transcript": row.get("transcript", ""),
                "gene_id": row.get("gene_id", ""),
                "loeuf": num("lof.oe_ci.upper"),
                "oe_lof": num("lof.oe"),
                "pli": num("lof.pLI"),
                "lof_exp": num("lof.exp") or 0.0,
                "mis_z": num("mis.z_score"),
                "syn_z": num("syn.z_score"),
                "_rank": rank,
            }
            prev = best.get(gene)
            if prev is None or (rec["_rank"], rec["lof_exp"]) > (prev["_rank"], prev["lof_exp"]):
                best[gene] = rec
    return best


def fmt(v: float | None, places: int = 3) -> str:
    return "NA" if v is None else f"{v:.{places}f}"


def annotate(panel: pathlib.Path, constraint: dict[str, dict]) -> tuple[int, int]:
    lines = panel.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:] if l.strip()]

    matched = 0
    for r in rows:
        c = constraint.get(r["symbol"])
        if c is None:
            r["ensembl_gene_id"] = "NA"
            r["gnomad_loeuf"] = "NA"
            r["gnomad_pli"] = "NA"
            r["gnomad_oe_lof"] = "NA"
            r["gnomad_mis_z"] = "NA"
            r["gnomad_transcript"] = "NA"
            r["constraint_bin"] = "no_constraint_data"
            continue
        matched += 1
        r["ensembl_gene_id"] = c["gene_id"] or "NA"
        r["gnomad_loeuf"] = fmt(c["loeuf"])
        r["gnomad_pli"] = fmt(c["pli"])
        r["gnomad_oe_lof"] = fmt(c["oe_lof"])
        r["gnomad_mis_z"] = fmt(c["mis_z"], 2)
        r["gnomad_transcript"] = c["transcript"] or "NA"
        loeuf = c["loeuf"]
        r["constraint_bin"] = (
            "no_constraint_data" if loeuf is None else
            "highly_constrained" if loeuf < 0.35 else
            "constrained" if loeuf < 0.60 else
            "moderate" if loeuf < 1.00 else
            "unconstrained"
        )

    out_cols = [c for c in header if c not in
                ("ensembl_gene_id", "gnomad_loeuf", "gnomad_pli")]
    out_cols += ["ensembl_gene_id", "gnomad_transcript", "gnomad_loeuf",
                 "gnomad_oe_lof", "gnomad_pli", "gnomad_mis_z", "constraint_bin"]

    with panel.open("w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in out_cols) + "\n")
    return matched, len(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--constraint", default=CONSTRAINT_DEFAULT)
    ap.add_argument("--panel", default="config/gene_panels/mitotic_extended.tsv")
    args = ap.parse_args()

    cpath = pathlib.Path(args.constraint)
    if not cpath.exists():
        sys.exit(
            f"FATAL: {cpath} not found.\n"
            "Download it with:\n"
            "  curl -o refs/gnomad.v4.1.constraint_metrics.tsv \\\n"
            "    https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
            "constraint/gnomad.v4.1.constraint_metrics.tsv"
        )

    constraint = load_constraint(cpath)
    sys.stderr.write(f"loaded constraint for {len(constraint):,} gene symbols ({GNOMAD_VERSION})\n")

    matched, total = annotate(pathlib.Path(args.panel), constraint)
    print(f"panel: {matched}/{total} genes matched to {GNOMAD_VERSION} constraint")

    # Report the constraint profile of the known MVA genes, since that is the
    # calibration that matters: if the known answers are unconstrained, the
    # score must not lean on constraint.
    lines = pathlib.Path(args.panel).read_text().splitlines()
    header = lines[0].split("\t")
    print(f"\n{'gene':10} {'LOEUF':>8} {'pLI':>8}  bin")
    for line in lines[1:]:
        r = dict(zip(header, line.split("\t")))
        if r.get("known_mva_gene") == "yes":
            print(f"{r['symbol']:10} {r['gnomad_loeuf']:>8} {r['gnomad_pli']:>8}  {r['constraint_bin']}")


if __name__ == "__main__":
    main()
