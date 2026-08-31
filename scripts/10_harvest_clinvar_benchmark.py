#!/usr/bin/env python3
"""Harvest published MVA-causal variants from ClinVar into the positive-control
benchmark required by MVA_HACKATHON_PLAN.md section 5.1.

Why ClinVar rather than reading papers directly: every field emitted here is
traceable to a VCV accession and, where ClinVar records one, a PMID. Nothing is
recalled from memory. CLAUDE.md rule 2 forbids inventing an identifier, and a
benchmark built from LLM recall of the MVA literature would be exactly that
failure mode, dressed up as ground truth.

The benchmark is used two ways (plan section 5.2 and 5.3):
  1. Spike each variant into a background genome and check the pipeline ranks
     its gene in the top N. The `canonical_spdi` field carries the exact
     reference and alternate alleles needed to construct a spike-in record.
  2. Re-run with these same ClinVar records masked, since the pipeline would
     otherwise be scored on its ability to look up the answer.

Output: benchmarks/published_mva_variants.tsv
        benchmarks/clinvar_raw/<gene>.json   (unmodified API responses)

Usage:
    python3 scripts/10_harvest_clinvar_benchmark.py [--email you@example.com]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Genes from plan section 0.1, plus CENATAC which causes MVA and is absent from
# the plan's table, and CEP57L1 as a paralogue candidate.
GENES = ["BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192", "SMC5", "CENATAC", "CEP57L1"]

# A ClinVar trait counts as MVA-relevant if its name matches any of these.
# Deliberately broader than "mosaic variegated aneuploidy" so that the Atelis
# and near-tetraploidy phenotypes, which overlap MVA clinically, are captured.
TRAIT_PATTERNS = [
    r"mosaic variegated aneuploidy",
    r"\bMVA\b",
    r"atelis",
    r"premature chromatid separation",
    r"aneuploidy",
    r"chromosome (?:instability|breakage)",
    r"tetraploid",
]

# OMIM numbers of the established MVA loci, used as a second route to relevance
# in case a trait name is phrased unusually.
MVA_OMIM = {"257300", "614114", "617598", "620185", "620184"}

PATHOGENIC = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}


def api(endpoint: str, params: dict, retries: int = 4) -> dict:
    """Call an E-utilities endpoint, returning parsed JSON. Retries on transient
    failure, because NCBI throttles and a half-populated benchmark is worse than
    a slow one.

    ``doseq`` matters: elink merges comma-separated ids into a single linkset,
    which silently attributes every PMID to the first uid. Repeated ``id=``
    parameters return one linkset per uid, which is what the caller needs.
    """
    params = {**params, "retmode": "json"}
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params, doseq=True)
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return json.load(fh)
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"E-utilities call failed after {retries} attempts: {url}\n{last}")


def trait_is_mva(trait_set: list[dict]) -> tuple[bool, str, str]:
    """Return (relevant, joined trait names, joined OMIM/MONDO identifiers)."""
    names, xrefs = [], []
    relevant = False
    for trait in trait_set or []:
        name = trait.get("trait_name", "")
        names.append(name)
        for xref in trait.get("trait_xrefs", []) or []:
            ident = f"{xref.get('db_source')}:{xref.get('db_id')}"
            xrefs.append(ident)
            if xref.get("db_source") == "OMIM" and xref.get("db_id") in MVA_OMIM:
                relevant = True
        if any(re.search(p, name, re.I) for p in TRAIT_PATTERNS):
            relevant = True
    return relevant, "; ".join(dict.fromkeys(names)), "; ".join(dict.fromkeys(xrefs))


def classify_variant(hgvs: str, obj_type: str, consequences: list) -> str:
    """Assign a coarse variant class. Recall is reported broken down by this
    field, because a pipeline that recovers coding loss of function but misses
    every deep-intronic positive has a specific, actionable weakness that an
    aggregate recall figure would hide (plan section 5.2)."""
    text = " ".join(str(c) for c in (consequences or [])).lower()
    c_part = hgvs.split(":c.")[-1] if ":c." in hgvs else ""

    if "frameshift" in text or re.search(r"(del|dup|ins)", c_part) and obj_type in (
        "Deletion", "Duplication", "Insertion"):
        base = "frameshift_or_indel"
    else:
        base = None

    # Intronic positions are written as c.<exonpos><+|-><offset>.
    m = re.search(r"[+-](\d+)", c_part.split("_")[0])
    if m and not c_part.startswith("*"):
        offset = int(m.group(1))
        if offset <= 2:
            return "splice_site_canonical"
        if offset <= 10:
            return "splice_region_near"
        return "deep_intronic"

    if c_part.startswith("*") or c_part.startswith("-"):
        return "utr_or_promoter"
    if "nonsense" in text or re.search(r"p\.\w+\d+(Ter|\*)", hgvs):
        return "nonsense"
    if "synonymous" in text:
        return "synonymous"
    if "missense" in text or re.search(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", hgvs):
        return "missense"
    if base:
        return base
    if obj_type and obj_type != "single nucleotide variant":
        return obj_type.replace(" ", "_").lower()
    return "unclassified"


def harvest_gene(gene: str, raw_dir: pathlib.Path, email: str | None) -> list[dict]:
    common = {"db": "clinvar"}
    if email:
        common["email"] = email

    # Every significance, not only pathogenic. The variant classes this project
    # most needs to be tested on, deep intronic and cryptic splice alleles, are
    # under-represented among confidently pathogenic ClinVar records precisely
    # because they are hard to find. Restricting the search to P/LP produces a
    # benchmark made of easy nonsense and frameshift variants, which a pipeline
    # can score well on while remaining blind to the case in hand. Rows are
    # tiered by confidence instead, and recall is reported per class.
    term = f'{gene}[gene]'
    search = api("esearch.fcgi", {**common, "term": term, "retmax": "3000"})
    uids = search["esearchresult"].get("idlist", [])
    sys.stderr.write(f"  {gene}: {len(uids)} ClinVar records total\n")
    if not uids:
        return []

    # esummary in batches, then one elink call for all UIDs to get PMIDs.
    summaries: dict[str, dict] = {}
    for i in range(0, len(uids), 100):
        batch = uids[i:i + 100]
        res = api("esummary.fcgi", {**common, "id": ",".join(batch)})
        summaries.update({k: v for k, v in res["result"].items() if k != "uids"})
        time.sleep(0.4)

    pmids: dict[str, list[str]] = {}
    for i in range(0, len(uids), 50):
        batch = uids[i:i + 50]
        # cmd=neighbor_history would lose the per-UID mapping, so request
        # per-id linking explicitly.
        res = api("elink.fcgi", {"dbfrom": "clinvar", "db": "pubmed",
                                 "id": batch, "linkname": "clinvar_pubmed"})
        for ls in res.get("linksets", []):
            uid = ls["ids"][0] if ls.get("ids") else None
            links = []
            for db in ls.get("linksetdbs", []) or []:
                links.extend(db.get("links", []))
            if uid:
                pmids[uid] = links
        time.sleep(0.4)

    (raw_dir / f"{gene}.json").write_text(
        json.dumps({"query": term, "uids": uids, "summaries": summaries, "pmids": pmids}, indent=1))

    rows = []
    for uid, rec in summaries.items():
        germ = rec.get("germline_classification") or {}
        sig = (germ.get("description") or "").strip()

        relevant, trait_names, trait_xrefs = trait_is_mva(germ.get("trait_set", []))
        if not relevant:
            continue

        # Tier 1 is the scoring set: confidently pathogenic, so a miss is
        # unambiguously a pipeline failure. Tier 2 is reported but uncertain,
        # useful for measuring sensitivity to the hard classes without letting
        # an uncertain call count as ground truth.
        tier = 1 if sig.lower() in PATHOGENIC else 2

        vset = (rec.get("variation_set") or [{}])[0]
        loc38 = next((l for l in vset.get("variation_loc", [])
                      if l.get("assembly_name") == "GRCh38"), {})
        spdi = vset.get("canonical_spdi", "")
        ref = alt = ""
        if spdi.count(":") == 3:
            _, _, ref, alt = spdi.split(":")

        hgvs = vset.get("variation_name", rec.get("title", ""))
        rows.append({
            "gene": gene,
            "clinvar_vcv": rec.get("accession_version", rec.get("accession", "")),
            "clinvar_uid": uid,
            "hgvs_c": hgvs,
            "hgvs_p": rec.get("protein_change", "") or "",
            "variant_class": classify_variant(hgvs, vset.get("variant_type", ""),
                                              rec.get("molecular_consequence_list")),
            "variant_type": vset.get("variant_type", ""),
            "build": "GRCh38",
            "chrom_nochr": loc38.get("chr", ""),
            "pos_grch38": loc38.get("start", ""),
            "stop_grch38": loc38.get("stop", ""),
            "canonical_spdi": spdi,
            "ref": ref,
            "alt": alt,
            "clinical_significance": sig,
            "benchmark_tier": tier,
            "review_status": germ.get("review_status", ""),
            "last_evaluated": germ.get("last_evaluated", ""),
            "trait_names": trait_names,
            "trait_xrefs": trait_xrefs,
            "pmids": ";".join(pmids.get(uid, [])),
            # Plan section 5.1 asks whether the source paper reported that a
            # standard pipeline missed the variant. ClinVar does not record
            # this, so it must be filled by reading the cited paper.
            "missed_by_standard_pipeline": "TODO(source)",
            "notes": "",
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", default=None, help="contact address for the NCBI API")
    ap.add_argument("--out", default="benchmarks/published_mva_variants.tsv")
    args = ap.parse_args()

    raw_dir = pathlib.Path("benchmarks/clinvar_raw")
    raw_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    for gene in GENES:
        all_rows.extend(harvest_gene(gene, raw_dir, args.email))
        time.sleep(0.5)

    # Order so that the classes the pipeline is most likely to miss sort first.
    priority = {"deep_intronic": 0, "splice_region_near": 1, "utr_or_promoter": 2,
                "synonymous": 3, "splice_site_canonical": 4}
    all_rows.sort(key=lambda r: (priority.get(r["variant_class"], 9), r["gene"], r["hgvs_c"]))

    cols = list(all_rows[0].keys()) if all_rows else []
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

    print(f"\nWrote {len(all_rows)} rows to {out}")
    by_class: dict[str, int] = {}
    by_gene: dict[str, int] = {}
    no_pmid = 0
    for r in all_rows:
        by_class[r["variant_class"]] = by_class.get(r["variant_class"], 0) + 1
        by_gene[r["gene"]] = by_gene.get(r["gene"], 0) + 1
        if not r["pmids"]:
            no_pmid += 1
    print("\nby variant class:")
    for k, v in sorted(by_class.items(), key=lambda kv: -kv[1]):
        print(f"  {k:24} {v}")
    print("\nby gene:")
    for k, v in sorted(by_gene.items(), key=lambda kv: -kv[1]):
        print(f"  {k:24} {v}")
    print(f"\nrows with no linked PMID: {no_pmid}/{len(all_rows)}")


if __name__ == "__main__":
    main()
