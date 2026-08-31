#!/usr/bin/env python3
"""Assemble config/gene_panels/mitotic_extended.tsv, per plan section 6.5.

Every gene carries the specific source that nominated it, so the panel is
reproducible and auditable rather than a list someone once wrote down. A gene
nominated by several independent sources is more defensible than one nominated
by a single loose GO annotation, so the nomination count is retained and used
for tiering.

Sources
-------
GO        QuickGO annotation search, descendants included, human only, with
          evidence codes recorded. Weak annotations propagated from electronic
          inference are kept but marked.
Reactome  Pathway participants via the ContentService.
STRING    First-shell physical and functional partners of the known MVA genes
          above a high confidence threshold.

CORUM is named in the plan but is not included: it has no stable unauthenticated
API and its bulk download terms need checking first. Recorded as a gap rather
than silently dropped, and noted in panel_provenance.md.

Constraint (gnomAD LOEUF and pLI) is joined by scripts/12_join_constraint.py,
which needs a bulk file download and is kept separate so this script stays
runnable offline-ish and fast.

Usage:
    python3 scripts/11_build_mitotic_panel.py
"""

from __future__ import annotations

import collections
import json
import pathlib
import sys
import time
import urllib.parse
import urllib.request

OUT_TSV = pathlib.Path("config/gene_panels/mitotic_extended.tsv")
RAW = pathlib.Path("results/recon/panel_sources")

# Plan section 6.5 names these five GO terms explicitly. Descendants are
# included, because the specific child terms are where the informative
# annotations sit.
GO_TERMS = {
    "GO:0007094": "mitotic spindle assembly checkpoint signaling",
    "GO:0000776": "kinetochore",
    "GO:0007098": "centrosome cycle",
    "GO:0007062": "sister chromatid cohesion",
    "GO:0007059": "chromosome segregation",
    # Added beyond the plan: the SMC5/6 and DNA-repair axis that Atelis syndrome
    # implicates, and the APC/C which the plan names in prose but not as a term.
    "GO:0000819": "sister chromatid segregation",
    "GO:0051983": "regulation of chromosome segregation",
}

REACTOME_PATHWAYS = {
    "R-HSA-69618": "Mitotic Spindle Checkpoint",
    "R-HSA-2500257": "Resolution of Sister Chromatid Cohesion",
    "R-HSA-380259": "Loss of Nlp from mitotic centrosomes",
    "R-HSA-380270": "Recruitment of mitotic centrosome proteins and complexes",
    "R-HSA-141424": "Amplification of signal from the kinetochores",
    "R-HSA-174184": "Cdh1:APC/C mediated degradation of Cdc20 and other targets",
    "R-HSA-2467813": "Separation of Sister Chromatids",
    "R-HSA-68877": "Mitotic Prometaphase",
}

# Seeds for the STRING first shell: the established and candidate MVA genes.
STRING_SEEDS = ["BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192", "SMC5", "CENATAC"]
STRING_MIN_SCORE = 700          # 0.7, STRING's "high confidence" cutoff
STRING_MAX_PARTNERS = 60        # per seed

# Evidence codes that indicate an annotation was inferred electronically rather
# than asserted by a curator reading an experiment.
WEAK_EVIDENCE = {"IEA", "ISS", "ISO", "ISA", "ISM", "IBA", "IRD", "IKR", "RCA"}

# GO terms differ enormously in how much they narrow the search. "chromosome
# segregation" annotates 549 human symbols and says little; "mitotic spindle
# assembly checkpoint signaling" annotates 75 and is close to the disease
# mechanism. Treating them as equal nominations produced a 1,296 gene panel
# whose bulk was single weak annotations to the broad terms. These are the
# terms specific enough that one curated annotation is worth acting on.
SPECIFIC_GO = {"GO:0007094", "GO:0000776", "GO:0007098", "GO:0007062", "GO:0051983"}


def load_hgnc(path: pathlib.Path = pathlib.Path("refs/hgnc_complete_set.txt")):
    """Return (approved symbols, alias-or-previous -> approved symbol).

    GO and Reactome both carry withdrawn symbols: CASC5 for KNL1, KNTC2 for
    NDC80, SGOL1 for SGO1, WAPAL for WAPL, CPAP for CENPJ. Left unmapped these
    silently fail every downstream join, which is how a gene quietly drops out
    of a panel it belongs in.
    """
    if not path.exists():
        sys.stderr.write(f"  WARNING: {path} absent, symbols will not be normalised\n")
        return set(), {}
    import csv as _csv
    approved, alias = set(), {}
    with path.open(newline="") as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            sym = row["symbol"]
            approved.add(sym)
            for field in ("alias_symbol", "prev_symbol"):
                for a in (row.get(field) or "").split("|"):
                    a = a.strip()
                    # Never let an alias override a symbol that is itself
                    # approved for a different gene.
                    if a and a not in alias:
                        alias[a] = sym
    return approved, alias


def fetch(url: str, accept: str = "application/json", retries: int = 4) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": accept,
                                               "User-Agent": "mva-hackathon-2026"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90) as fh:
                return fh.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed after {retries} attempts: {url}\n{last}")


# ---------------------------------------------------------------------------

def from_go() -> dict[str, list[dict]]:
    """QuickGO annotation search, paged, human, descendants included."""
    nominations: dict[str, list[dict]] = collections.defaultdict(list)
    for go_id, label in GO_TERMS.items():
        seen_pages, page, total = 0, 1, None
        rows = []
        while True:
            q = urllib.parse.urlencode({
                "goId": go_id, "goUsage": "descendants", "taxonId": "9606",
                "limit": "200", "page": str(page),
            })
            data = json.loads(fetch(f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?{q}"))
            total = data.get("numberOfHits", 0)
            results = data.get("results", [])
            rows.extend(results)
            if not results or len(rows) >= total or page >= 25:
                break
            page += 1
            seen_pages += 1
            time.sleep(0.2)

        for r in rows:
            sym = r.get("symbol")
            if not sym:
                continue
            nominations[sym].append({
                "source": "GO",
                "detail": f"{go_id} ({label})",
                "evidence": r.get("goEvidence", ""),
                "weak": r.get("goEvidence", "") in WEAK_EVIDENCE,
            })
        sys.stderr.write(f"  GO {go_id} {label}: {total} annotations, "
                         f"{len({r.get('symbol') for r in rows})} symbols\n")
        (RAW / f"go_{go_id.replace(':', '_')}.json").write_text(json.dumps(rows[:5000]))
        time.sleep(0.3)
    return nominations


def from_reactome() -> dict[str, list[dict]]:
    nominations: dict[str, list[dict]] = collections.defaultdict(list)
    for stid, label in REACTOME_PATHWAYS.items():
        try:
            data = json.loads(fetch(
                f"https://reactome.org/ContentService/data/participants/{stid}/referenceEntities"))
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"  Reactome {stid} FAILED: {exc}\n")
            continue
        syms = set()
        for ent in data:
            # geneName is an ordered alias list whose first element is the
            # primary symbol. Adding every element treats aliases as separate
            # genes: MAD1L1 also arrived as MAD1 and TXBP181, TAOK1 as
            # KIAA1361, MAP3K16 and MARKK. Take the primary symbol only.
            if ent.get("schemaClass") not in ("ReferenceGeneProduct", "ReferenceIsoform"):
                continue
            if ent.get("databaseName") != "UniProt":
                continue
            names = ent.get("geneName") or []
            if names:
                syms.add(names[0])
        for s in syms:
            nominations[s].append({"source": "Reactome", "detail": f"{stid} ({label})",
                                   "evidence": "pathway_participant", "weak": False})
        sys.stderr.write(f"  Reactome {stid} {label}: {len(syms)} gene symbols\n")
        (RAW / f"reactome_{stid}.json").write_text(json.dumps(sorted(syms)))
        time.sleep(0.3)
    return nominations


def from_string() -> dict[str, list[dict]]:
    nominations: dict[str, list[dict]] = collections.defaultdict(list)
    for seed in STRING_SEEDS:
        q = urllib.parse.urlencode({
            "identifiers": seed, "species": "9606",
            "required_score": str(STRING_MIN_SCORE), "limit": str(STRING_MAX_PARTNERS),
        })
        try:
            text = fetch(f"https://string-db.org/api/tsv/interaction_partners?{q}",
                         accept="text/plain").decode()
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"  STRING {seed} FAILED: {exc}\n")
            continue
        lines = text.strip().splitlines()
        if len(lines) < 2:
            continue
        head = lines[0].split("\t")
        i_b, i_score = head.index("preferredName_B"), head.index("score")
        partners = []
        for line in lines[1:]:
            f = line.split("\t")
            partners.append((f[i_b], float(f[i_score])))
        for name, score in partners:
            nominations[name].append({
                "source": "STRING", "detail": f"first shell of {seed} (score {score:.3f})",
                "evidence": f"combined_score={score:.3f}", "weak": False})
        sys.stderr.write(f"  STRING {seed}: {len(partners)} partners at score >= "
                         f"{STRING_MIN_SCORE / 1000:.2f}\n")
        (RAW / f"string_{seed}.tsv").write_text(text)
        time.sleep(0.4)
    return nominations


# ---------------------------------------------------------------------------

def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    approved, alias = load_hgnc()
    sys.stderr.write(f"HGNC: {len(approved):,} approved symbols, {len(alias):,} aliases\n")

    sys.stderr.write("GO:\n")
    raw_noms = collections.defaultdict(list)
    for src in (from_go(), from_reactome(), from_string()):
        for sym, entries in src.items():
            raw_noms[sym].extend(entries)

    # Normalise to approved HGNC symbols and drop anything that is not a gene.
    noms = collections.defaultdict(list)
    dropped, renamed = [], []
    for sym, entries in raw_noms.items():
        if approved and sym not in approved:
            mapped = alias.get(sym)
            if mapped:
                renamed.append((sym, mapped))
                sym = mapped
            else:
                dropped.append(sym)
                continue
        noms[sym].extend(entries)
    sys.stderr.write(f"normalisation: {len(renamed)} symbols renamed, "
                     f"{len(dropped)} unrecognised symbols dropped\n")
    (RAW / "symbol_normalisation.json").write_text(json.dumps(
        {"renamed": sorted(renamed), "dropped": sorted(dropped)}, indent=1))

    # Known MVA genes are always in the panel, whatever the sources returned.
    known = {}
    for line in pathlib.Path("config/gene_panels/mva_known.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        known[f[0]] = f[9]  # tier
        noms[f[0]].append({"source": "known_MVA_gene", "detail": f[10],
                           "evidence": "curated", "weak": False})

    rows = []
    for sym, entries in noms.items():
        sources = sorted({e["source"] for e in entries})
        strong = [e for e in entries if not e["weak"]]
        go_terms = sorted({e["detail"] for e in entries if e["source"] == "GO"})
        specific_go = sorted({e["detail"] for e in entries
                              if e["source"] == "GO" and not e["weak"]
                              and e["detail"].split()[0] in SPECIFIC_GO})
        reactome = sorted({e["detail"] for e in entries if e["source"] == "Reactome"})
        string_seeds = sorted({e["detail"].split()[3] for e in entries
                               if e["source"] == "STRING"})

        # Tiering. Known genes, then converging evidence, then a curated
        # annotation to a mechanism-specific GO term, then the long tail.
        # Only tiers 1 to 3 form the core panel that the analysis arms use;
        # tier 4 is retained in the file so the cut is visible and revisable,
        # not silently discarded.
        if sym in known:
            tier = 1
        elif len(sources) >= 3 and strong:
            tier = 2
        elif (len(sources) >= 2 and strong) or len(specific_go) >= 2:
            tier = 3
        elif specific_go:
            tier = 4
        else:
            tier = 5

        rows.append({
            "symbol": sym,
            "panel_tier": tier,
            "n_sources": len(sources),
            "sources": ",".join(sources),
            "n_nominations": len(entries),
            "has_experimental_evidence": "yes" if strong else "no",
            "in_core_panel": "",   # set below, once tiers are known
            "n_specific_go_terms": len(specific_go),
            "specific_go_terms": "; ".join(specific_go),
            "go_terms": "; ".join(go_terms),
            "reactome_pathways": "; ".join(reactome),
            "string_seeds": ",".join(string_seeds),
            "known_mva_gene": "yes" if sym in known else "no",
            # Filled by scripts/12_join_constraint.py.
            "ensembl_gene_id": "TODO(source)",
            "gnomad_loeuf": "TODO(source)",
            "gnomad_pli": "TODO(source)",
        })

    # The core panel is tiers 1 to 4, which lands in the 300 to 500 range the
    # plan asks for. Tier 5 is the long tail of single weak annotations to broad
    # terms and is excluded from analysis but kept in the file.
    for r in rows:
        r["in_core_panel"] = "yes" if r["panel_tier"] <= 4 else "no"

    rows.sort(key=lambda r: (r["panel_tier"], -r["n_sources"], r["symbol"]))
    cols = list(rows[0].keys())
    with OUT_TSV.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    by_tier = collections.Counter(r["panel_tier"] for r in rows)
    by_src = collections.Counter(s for r in rows for s in r["sources"].split(","))
    print(f"\nWrote {len(rows)} genes to {OUT_TSV}")
    print("\nby tier:")
    for t in sorted(by_tier):
        print(f"  tier {t}: {by_tier[t]}")
    print("\ngenes nominated by each source:")
    for s, n in by_src.most_common():
        print(f"  {s:16} {n}")
    multi = sum(1 for r in rows if r["n_sources"] >= 2)
    core = sum(1 for r in rows if r["in_core_panel"] == "yes")
    print(f"\nnominated by 2 or more independent sources: {multi}")
    print(f"CORE PANEL (tiers 1-4, used by the analysis arms): {core} genes")
    print(f"long tail (tier 5, retained but not used): {len(rows) - core} genes")


if __name__ == "__main__":
    main()
