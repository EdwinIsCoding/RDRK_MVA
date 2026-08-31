#!/usr/bin/env python3
"""Build the widened disease-gene panel for Arm E.

Why widen, and why not to the whole genome
------------------------------------------
The mitotic panel (408 genes) encodes an assumption: that the answer lies in
chromosome-segregation machinery. Two things argue against relying on that.

First, the proband's phenotype does not match canonical MVA closely. There is no
microcephaly, no seizures, no developmental delay. Rhabdomyosarcoma with IUGR,
nephrocalcinosis and parental recurrent miscarriage overlaps several cancer
predisposition and growth syndromes that have nothing to do with the spindle
assembly checkpoint.

Second, and decisively, the organisers describe Track 1 as "a foundational track
- it's designed to be achievable", scored against a confirmed clinical answer.
An achievable answer is far more likely to sit in an established disease gene
than in a novel one.

So the widening is to **genes with curated disease associations**, not to the
whole genome. A genome-wide search over 5 million variants without a gene prior
is not more thorough, it is less discriminating: it would bury a real answer
under thousands of rare variants in genes nobody has ever linked to disease.

Sources, all public and versioned
---------------------------------
ClinGen      gene-disease validity classifications. Definitive, Strong and
             Moderate are kept; Limited, Disputed and Refuted are recorded but
             tiered down, since a disputed association is a weak prior.
PanelApp     Genomics England curated panels, green (diagnostic-grade) entries.
             Panels selected for relevance to this phenotype plus the broad
             intellectual disability and cancer predisposition panels.
gene2phenotype  DDG2P developmental disorder gene panel from EBI.

Output: config/gene_panels/disease_genes.tsv
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import pathlib
import sys
import time
import urllib.request

OUT = pathlib.Path("config/gene_panels/disease_genes.tsv")
RAW = pathlib.Path("results/recon/disease_panel_sources")

#: PanelApp panels chosen for this phenotype, each matched to the proband's
#: HPO terms. **Every id here was resolved by searching the PanelApp panel index
#: by name, not guessed.** An earlier revision hard-coded ids from memory and
#: six of the nine were wrong: 391 is "Adult solid tumours" not rhabdomyosarcoma,
#: 478 is "Fetal anomalies" not cancer predisposition, 115 is head and neck
#: cancer not Fanconi anaemia. The gene content fetched was real but the
#: recorded provenance was false, which is worse than a gap.
#:
#: The name is verified against the API response at fetch time (see
#: from_panelapp), so a future id drift fails loudly instead of silently
#: relabelling.
PANELAPP_PANELS = {
    # Rhabdomyosarcoma, HP:0002859, the presenting event
    "290": "Familial rhabdomyosarcoma",
    "243": "Childhood solid tumours",
    "259": "Childhood solid tumours cancer susceptibility",
    "1320": "DICER1-related cancer predisposition",
    # Chromosome instability, the MVA differential
    "508": "Fanconi anaemia or Bloom syndrome",
    # Nephrocalcinosis, HP:0000121, present since birth
    "149": "Nephrocalcinosis or nephrolithiasis",
    "487": "Cystic renal disease",
    # Growth restriction and failure to thrive, HP:0001518 and HP:0001508
    "38": "Beckwith-Wiedemann syndrome (BWS) and other congenital overgrowth disorders",
    "162": "Severe microcephaly",
    # Broad developmental panel, kept as a wide net
    "285": "Intellectual disability",
    # Secondary findings, paediatric. The hackathon rules explicitly welcome
    # incidental findings and they do not hurt the automated score.
    "930": "Additional findings health related - children",
}

STRONG = {"definitive", "strong", "moderate"}


def fetch(url: str, retries: int = 4, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={
        "User-Agent": "mva-hackathon-2026", "Accept": "application/json, text/csv, */*"})
    last = None
    for a in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (a + 1))
    raise RuntimeError(f"failed: {url}\n{last}")


def from_clingen() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    try:
        raw = fetch("https://search.clinicalgenome.org/kb/gene-validity/download")
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"  ClinGen FAILED: {exc}\n")
        return out
    (RAW / "clingen.csv").write_bytes(raw)
    text = raw.decode("utf-8", errors="replace")
    # The download carries several preamble lines before the real header.
    lines = text.splitlines()
    # The header row is CSV-quoted, so it begins with a double quote and a
    # bare startswith("GENE SYMBOL") never matches. Strip quotes first.
    start = next((i for i, l in enumerate(lines)
                  if l.lstrip('"').upper().startswith("GENE SYMBOL")), 0)
    rdr = csv.reader(io.StringIO("\n".join(lines[start:])))
    header = next(rdr, [])
    idx = {h.strip().upper(): i for i, h in enumerate(header)}
    g_i = idx.get("GENE SYMBOL")
    d_i = idx.get("DISEASE LABEL")
    c_i = idx.get("CLASSIFICATION")
    if g_i is None:
        sys.stderr.write("  ClinGen: unexpected header, skipped\n")
        return out
    n = 0
    for row in rdr:
        if len(row) <= g_i or not row[g_i].strip():
            continue
        sym = row[g_i].strip()
        # The banner rows use runs of '+' as separators; skip them.
        if set(sym) <= {"+"} or sym.upper() == "GENE SYMBOL":
            continue
        cls = (row[c_i].strip() if c_i is not None and len(row) > c_i else "")
        dis = (row[d_i].strip() if d_i is not None and len(row) > d_i else "")
        out.setdefault(sym, []).append({
            "source": "ClinGen", "detail": f"{dis} ({cls})".strip(),
            "strong": cls.lower() in STRONG})
        n += 1
    sys.stderr.write(f"  ClinGen: {n} assertions, {len(out)} genes\n")
    return out


def from_panelapp() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for pid, name in PANELAPP_PANELS.items():
        try:
            d = json.loads(fetch(
                f"https://panelapp.genomicsengland.co.uk/api/v1/panels/{pid}/?format=json"))
        except Exception as exc:  # noqa: BLE001
            sys.stderr.write(f"  PanelApp {pid} {name}: FAILED {exc}\n")
            continue
        # Verify the panel is the one intended. Ids are not stable across
        # PanelApp releases and a silent substitution would poison provenance.
        actual = d.get("name", "")
        if actual.strip().lower() != name.strip().lower():
            sys.stderr.write(
                f"  PanelApp {pid}: FATAL name mismatch. Expected {name!r}, "
                f"API returned {actual!r}. Skipping rather than mislabelling.\n")
            continue
        green = 0
        for g in d.get("genes", []):
            sym = (g.get("gene_data") or {}).get("gene_symbol")
            conf = str(g.get("confidence_level", ""))
            if not sym:
                continue
            # Confidence 3 is green: diagnostic-grade. 2 amber, 1 red.
            out.setdefault(sym, []).append({
                "source": "PanelApp", "detail": f"{name} (confidence {conf})",
                "strong": conf == "3"})
            if conf == "3":
                green += 1
        sys.stderr.write(f"  PanelApp {name}: {green} green of {len(d.get('genes', []))}\n")
        (RAW / f"panelapp_{pid}.json").write_text(json.dumps(
            [(g.get('gene_data') or {}).get('gene_symbol') for g in d.get("genes", [])]))
        time.sleep(0.4)
    return out


def from_g2p() -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    try:
        # The old /downloads/DDG2P.csv.gz path now returns the React app.
        # The live CSV endpoint is under the API, verified 31 August 2026.
        raw = fetch("https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/")
    except Exception as exc:  # noqa: BLE001
        sys.stderr.write(f"  gene2phenotype FAILED: {exc}\n")
        return out
    (RAW / "ddg2p.csv").write_bytes(raw)
    try:
        text = gzip.decompress(raw).decode("utf-8", errors="replace")
    except Exception:
        text = raw.decode("utf-8", errors="replace")
    rdr = csv.DictReader(io.StringIO(text))
    n = 0
    for row in rdr:
        sym = (row.get("gene symbol") or row.get("gene_symbol")
               or row.get("hgnc symbol") or "").strip()
        if not sym:
            continue
        conf = (row.get("confidence category") or row.get("DDD category")
                or row.get("confidence_category") or "").strip().lower()
        dis = (row.get("disease name") or row.get("disease_name") or "").strip()
        out.setdefault(sym, []).append({
            "source": "gene2phenotype", "detail": f"{dis} ({conf})".strip(),
            "strong": conf in ("definitive", "strong", "moderate", "confirmed", "probable")})
        n += 1
    sys.stderr.write(f"  gene2phenotype: {n} records, {len(out)} genes\n")
    return out


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)

    merged: dict[str, list[dict]] = {}
    sys.stderr.write("sources:\n")
    for src in (from_clingen(), from_panelapp(), from_g2p()):
        for sym, entries in src.items():
            merged.setdefault(sym, []).extend(entries)

    # Carry the mitotic panel through so the widened set is a superset.
    mito: dict[str, str] = {}
    p = pathlib.Path("config/gene_panels/mitotic_extended.tsv")
    if p.exists():
        for r in csv.DictReader(p.open(newline=""), delimiter="\t"):
            if r["in_core_panel"] == "yes":
                mito[r["symbol"]] = r["panel_tier"]
                merged.setdefault(r["symbol"], []).append({
                    "source": "mitotic_panel", "detail": f"tier {r['panel_tier']}",
                    "strong": r["panel_tier"] in ("1", "2")})

    rows = []
    for sym, entries in merged.items():
        sources = sorted({e["source"] for e in entries})
        strong = [e for e in entries if e["strong"]]
        rows.append({
            "symbol": sym,
            "tier": 1 if strong and len(sources) >= 2 else 2 if strong else 3,
            "n_sources": len(sources),
            "sources": ",".join(sources),
            "has_strong_assertion": "yes" if strong else "no",
            "in_mitotic_panel": "yes" if sym in mito else "no",
            "assertions": " | ".join(sorted({e["detail"] for e in entries})[:6])[:400],
        })
    rows.sort(key=lambda r: (r["tier"], -r["n_sources"], r["symbol"]))

    cols = list(rows[0].keys())
    with OUT.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    import collections
    by_tier = collections.Counter(r["tier"] for r in rows)
    print(f"\nwrote {len(rows)} genes to {OUT}")
    for t in sorted(by_tier):
        print(f"  tier {t}: {by_tier[t]}")
    print(f"  with a strong assertion: {sum(1 for r in rows if r['has_strong_assertion']=='yes')}")
    print(f"  also on the mitotic panel: {sum(1 for r in rows if r['in_mitotic_panel']=='yes')}")


if __name__ == "__main__":
    main()
