#!/usr/bin/env python3
"""Phenotype-driven gene prior from HPO annotations. Arm E, plan section 6.5.

Why not Exomiser or LIRICAL
---------------------------
Both are the right tools and both are named in the plan. Neither runs here:
they need Java 17 or later and this host has Java 11, and their data bundles are
several gigabytes competing with downloads already saturating the connection.
They should be run on the GPU host, and this module's output compared against
theirs as a cross-check rather than replaced by it.

What this computes instead
--------------------------
The same underlying signal, from the HPO annotation file directly: how specific
the overlap is between the proband's eight terms and the terms annotated to each
gene's known diseases.

Specificity matters more than count. "Autosomal recessive inheritance"
(HP:0000007) is annotated to thousands of genes and carries almost no
information; "Rhabdomyosarcoma" (HP:0002859) is annotated to very few and is
close to diagnostic on its own. So each matched term is weighted by its
information content, -log(p), where p is the fraction of annotated genes
carrying it. This is the Resnik-style weighting that Exomiser's phenotype
scoring is also built on.

Ancestor closure is applied where hp.obo is available, so a gene annotated with
a child term of one of the proband's terms still matches. Without the ontology
the match is exact-term only, which under-counts, and that is reported rather
than hidden.

**This is a gene prior, not evidence of causation.** A high score means the
gene's known phenotypes resemble this child's, which is a reason to look, not a
finding. It carries a modest weight in the additive score for that reason.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import pathlib
import sys

PROBAND_HPO = {
    "HP:0002859": "Rhabdomyosarcoma",
    "HP:0000121": "Nephrocalcinosis",
    "HP:0004322": "Short stature",
    "HP:0001508": "Failure to thrive",
    "HP:0003202": "Skeletal muscle atrophy",
    "HP:0001622": "Premature birth",
    "HP:0001518": "Small for gestational age",
    "HP:0200067": "Recurrent spontaneous abortion",
}

#: Terms describing inheritance or clinical modifiers rather than phenotype.
#: They are annotated to enormous numbers of genes and only add noise.
UNINFORMATIVE_PREFIXES = ("HP:0000005",)   # Mode of inheritance subtree root


def load_obo(path: pathlib.Path) -> dict[str, set[str]]:
    """term -> set of ancestors (including itself). Empty if the file is absent."""
    if not path.exists():
        return {}
    parents: dict[str, set[str]] = collections.defaultdict(set)
    cur = None
    for line in path.read_text().splitlines():
        if line == "[Term]":
            cur = None
        elif line.startswith("id: HP:"):
            cur = line[4:].strip()
        elif line.startswith("is_a:") and cur:
            parents[cur].add(line.split()[1])
    # Transitive closure, memoised.
    anc: dict[str, set[str]] = {}

    def climb(t: str, seen: frozenset = frozenset()) -> set[str]:
        if t in anc:
            return anc[t]
        if t in seen:
            return set()
        out = {t}
        for p in parents.get(t, ()):
            out |= climb(p, seen | {t})
        anc[t] = out
        return out

    for t in list(parents):
        climb(t)
    return anc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--g2p", default="refs/hpo/genes_to_phenotype.txt")
    ap.add_argument("--obo", default="refs/hpo/hp.obo")
    ap.add_argument("--panel", default="config/gene_panels/disease_genes.tsv")
    ap.add_argument("--out", default="config/gene_panels/phenotype_prior.tsv")
    args = ap.parse_args()

    g2p = pathlib.Path(args.g2p)
    if not g2p.exists():
        sys.exit(f"FATAL: {g2p} not found.")

    anc = load_obo(pathlib.Path(args.obo))
    if anc:
        sys.stderr.write(f"ontology: {len(anc):,} terms with ancestor closure\n")
    else:
        sys.stderr.write("ontology: hp.obo absent, falling back to exact-term matching. "
                         "This under-counts genes annotated with child terms.\n")

    gene_terms: dict[str, set[str]] = collections.defaultdict(set)
    with g2p.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym, term = row.get("gene_symbol"), row.get("hpo_id")
            if sym and term:
                gene_terms[sym].add(term)
    sys.stderr.write(f"annotations: {len(gene_terms):,} genes\n")

    # Expand each gene's terms to include ancestors, so a specific annotation
    # satisfies a more general proband term.
    if anc:
        gene_terms = {g: set().union(*(anc.get(t, {t}) for t in ts)) if ts else set()
                      for g, ts in gene_terms.items()}

    # Information content from the annotation corpus itself.
    n_genes = len(gene_terms)
    freq: collections.Counter = collections.Counter()
    for ts in gene_terms.values():
        freq.update(ts)
    ic = {t: -math.log(c / n_genes) for t, c in freq.items() if c}

    # The proband's terms, with the ancestor sets used for matching.
    query = {}
    for t, label in PROBAND_HPO.items():
        query[t] = {"label": label, "ic": ic.get(t)}
    sys.stderr.write("\nproband term specificity (information content, higher is rarer):\n")
    for t, d in sorted(query.items(), key=lambda kv: -(kv[1]["ic"] or 0)):
        n = freq.get(t, 0)
        sys.stderr.write(f"  {t}  {d['label']:32} IC {d['ic'] or 0:5.2f}  "
                         f"annotated to {n:,} genes\n")

    panel = {r["symbol"]: r for r in
             csv.DictReader(open(args.panel, newline=""), delimiter="\t")}

    rows = []
    max_possible = sum(v["ic"] or 0 for v in query.values())
    for gene in set(gene_terms) | set(panel):
        ts = gene_terms.get(gene, set())
        matched = [t for t in PROBAND_HPO if t in ts]
        score = sum(ic.get(t, 0.0) for t in matched)
        if not matched:
            continue
        rows.append({
            "symbol": gene,
            "phenotype_score": round(score, 3),
            "fraction_of_max": round(score / max_possible, 4) if max_possible else 0,
            "n_terms_matched": len(matched),
            "matched_terms": ";".join(f"{t}({PROBAND_HPO[t]})" for t in matched),
            "in_disease_panel": "yes" if gene in panel else "no",
        })
    rows.sort(key=lambda r: (-r["phenotype_score"], -r["n_terms_matched"], r["symbol"]))

    out = pathlib.Path(args.out)
    cols = list(rows[0].keys())
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\n{len(rows):,} genes share at least one term with the proband -> {out}")
    print(f"  maximum attainable score (all 8 terms): {max_possible:.2f}\n")
    print(f"  {'gene':12}{'score':>7}{'terms':>7}  matched")
    for r in rows[:20]:
        short = ";".join(t.split("(")[1].rstrip(")") for t in r["matched_terms"].split(";"))
        print(f"  {r['symbol']:12}{r['phenotype_score']:>7.2f}{r['n_terms_matched']:>7}  {short[:70]}")

    known = {"BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192", "SMC5", "CENATAC"}
    print("\n  where the known MVA genes rank:")
    for i, r in enumerate(rows, 1):
        if r["symbol"] in known:
            print(f"    {i:5}. {r['symbol']:10} score {r['phenotype_score']:.2f}  "
                  f"{r['n_terms_matched']} terms")


if __name__ == "__main__":
    main()
