#!/usr/bin/env python3
"""Harvest a deep-intronic and cryptic-splice mechanism control set.

Why this exists
---------------
``scripts/10_harvest_clinvar_benchmark.py`` produced 108 confidently pathogenic
variants across the MVA genes. Their class breakdown is:

    nonsense 49, frameshift 35, canonical splice 15, deletion 5, missense 3

**Zero deep intronic. Zero near-splice-region. Zero synonymous. Zero UTR.**

That is a problem, because the leading hypothesis for this proband
(MVA_HACKATHON_PLAN.md section 0.2 class 1, and RECON.md) is a cryptic second
allele that is deep intronic, at a branch point, or in a UTR. Scoring the
pipeline only on nonsense and frameshift positives would produce a flattering
recall figure that says nothing about whether the pipeline can find the kind of
variant we are actually looking for.

The reason ClinVar has none is not an oversight. These alleles are hard to find,
so they are under-ascertained, and when found they are often deposited as
uncertain rather than pathogenic. The very scarcity is the point.

What this set is, and what it is not
------------------------------------
This harvests confidently pathogenic **deep intronic and near-splice variants
from any gene and any disease**, then asks whether the splicing arm recovers
them. It tests the machinery: the plus or minus 500 bp window, the SpliceAI and
Pangolin invocation, the branch point scoring, the ranking.

It is **not** MVA-specific and must never be reported as though it were. A good
score here means "the splicing arm can find a deep intronic pathogenic variant";
it does not mean "the splicing arm can find the MVA allele". That distinction
belongs in the submission in those words.

Genes matching the MVA panel are excluded, so this set stays independent of the
benchmark it is meant to complement.

Output: benchmarks/splice_mechanism_controls.tsv
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERY = (
    '("intron variant"[molecular consequence]) '
    'AND ("pathogenic"[clinsig] OR "likely pathogenic"[clinsig]) '
    'AND ("criteria provided, multiple submitters, no conflicts"[review status] '
    'OR "reviewed by expert panel"[review status] '
    'OR "criteria provided, single submitter"[review status])'
)

#: Minimum distance from the nearest exon boundary, in bases, for a variant to
#: be informative here. Beyond 10 bp the canonical splice-site motif no longer
#: explains the effect. The bands below matter more than the cutoff: the whole
#: point of the plan's plus or minus 500 bp window is to catch what a default
#: plus or minus 50 bp window misses, so recall must be reported per band.
DEEP_INTRONIC_MIN_OFFSET = 11

#: Offset bands for reporting. A tool run at the default plus or minus 50 bp
#: cannot see anything past the first band, by construction.
OFFSET_BANDS = [
    (11, 50, "intronic_11_50"),          # a default window still reaches these
    (51, 100, "deep_intronic_51_100"),   # needs a widened window
    (101, 500, "deep_intronic_101_500"), # the regime the plan targets
    (501, 10**9, "beyond_500bp_window"), # expected misses, kept to bound recall
]

ACCEPTED_SIGNIFICANCE = {
    "pathogenic",
    "likely pathogenic",
    "pathogenic/likely pathogenic",
}


def api(endpoint: str, params: dict, retries: int = 4) -> dict:
    params = {**params, "retmode": "json"}
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params, doseq=True)
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return json.load(fh)
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed after {retries} attempts: {url}\n{last}")


def intron_offset(hgvs: str) -> int | None:
    """Distance from the nearest exon boundary, from the HGVS c. description.

    ``c.1735-12T>C`` is 12 bases into the intron from the acceptor side;
    ``c.1234+56A>G`` is 56 bases from the donor side. Returns None where the
    variant is not intronic.
    """
    if ":c." not in hgvs:
        return None
    c = hgvs.split(":c.", 1)[1]
    # Ignore UTR coordinates, which also use +/- but mean something else.
    if c.startswith("*") or c.startswith("-"):
        return None
    m = re.match(r"^\d+([+-])(\d+)", c)
    return int(m.group(2)) if m else None


def band(offset: int) -> str:
    for lo, hi, name in OFFSET_BANDS:
        if lo <= offset <= hi:
            return name
    return "unbanded"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=9000,
                    help="how many ClinVar records to pull before filtering")
    ap.add_argument("--email", default=None)
    ap.add_argument("--out", default="benchmarks/splice_mechanism_controls.tsv")
    args = ap.parse_args()

    common: dict = {"db": "clinvar"}
    if args.email:
        common["email"] = args.email

    # Exclude the MVA panel so this set is independent of the main benchmark.
    panel_genes = set()
    for line in pathlib.Path("config/gene_panels/mitotic_extended.tsv").read_text().splitlines()[1:]:
        if line.strip():
            panel_genes.add(line.split("\t")[0])

    head = api("esearch.fcgi", {**common, "term": QUERY, "retmax": "1"})
    total = int(head["esearchresult"]["count"])

    # esearch returns UIDs in a fixed order, so taking the first N samples one
    # end of the database rather than the database. Walk the whole result set
    # in strides instead.
    uids: list[str] = []
    stride = 5000
    for start in range(0, min(total, args.sample * 4), stride):
        page = api("esearch.fcgi", {**common, "term": QUERY,
                                    "retstart": str(start), "retmax": str(stride)})
        uids.extend(page["esearchresult"].get("idlist", []))
        if len(uids) >= args.sample * 4:
            break
        time.sleep(0.35)
    # Deterministic thinning to the requested sample size, seeded so the
    # control set is reproducible.
    import random
    random.Random(20261024).shuffle(uids)
    uids = uids[:args.sample]
    sys.stderr.write(f"{total} matching records in ClinVar; sampling {len(uids)}\n")

    summaries: dict[str, dict] = {}
    for i in range(0, len(uids), 100):
        res = api("esummary.fcgi", {**common, "id": ",".join(uids[i:i + 100])})
        summaries.update({k: v for k, v in res["result"].items() if k != "uids"})
        time.sleep(0.4)
        sys.stderr.write(f"  fetched {len(summaries)}/{len(uids)}\r")
    sys.stderr.write("\n")

    pmids: dict[str, list[str]] = {}
    for i in range(0, len(uids), 50):
        res = api("elink.fcgi", {"dbfrom": "clinvar", "db": "pubmed",
                                 "id": uids[i:i + 50], "linkname": "clinvar_pubmed"})
        for ls in res.get("linksets", []):
            if ls.get("ids"):
                pmids[ls["ids"][0]] = [l for db in ls.get("linksetdbs", []) or []
                                       for l in db.get("links", [])]
        time.sleep(0.4)

    rows = []
    skipped_panel = 0
    for uid, rec in summaries.items():
        germ = rec.get("germline_classification") or {}
        if (germ.get("description") or "").strip().lower() not in ACCEPTED_SIGNIFICANCE:
            continue
        vset = (rec.get("variation_set") or [{}])[0]
        hgvs = vset.get("variation_name", rec.get("title", ""))
        offset = intron_offset(hgvs)
        if offset is None or offset < DEEP_INTRONIC_MIN_OFFSET:
            continue

        genes = [g.get("symbol") for g in (rec.get("genes") or [])]
        if any(g in panel_genes for g in genes):
            skipped_panel += 1
            continue

        loc38 = next((l for l in vset.get("variation_loc", [])
                      if l.get("assembly_name") == "GRCh38"), {})
        spdi = vset.get("canonical_spdi", "")
        ref = alt = ""
        if spdi.count(":") == 3:
            _, _, ref, alt = spdi.split(":")
        if not (loc38.get("chr") and loc38.get("start") and ref and alt):
            continue

        traits = "; ".join(t.get("trait_name", "")
                           for t in (germ.get("trait_set") or []))
        rows.append({
            "gene": genes[0] if genes else "",
            "clinvar_vcv": rec.get("accession_version", ""),
            "clinvar_uid": uid,
            "hgvs_c": hgvs,
            "intron_offset_bp": offset,
            "offset_band": band(offset),
            "variant_class": "deep_intronic" if offset > 50 else "splice_region_near",
            "build": "GRCh38",
            "chrom_nochr": loc38["chr"],
            "pos_grch38": loc38["start"],
            "canonical_spdi": spdi,
            "ref": ref,
            "alt": alt,
            "clinical_significance": germ.get("description", ""),
            "review_status": germ.get("review_status", ""),
            "trait_names": traits,
            "pmids": ";".join(pmids.get(uid, [])),
            "control_set": "mechanism_only_not_MVA_specific",
        })

    rows.sort(key=lambda r: (-int(r["intron_offset_bp"]), r["gene"]))
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\nWrote {len(rows)} deep-intronic pathogenic controls to {out}")
    print(f"  excluded because the gene is on the MVA panel: {skipped_panel}")
    offs = [int(r["intron_offset_bp"]) for r in rows]
    print(f"  intron offset: min {min(offs)} median {sorted(offs)[len(offs)//2]} max {max(offs)} bp")
    print("\n  by offset band:")
    for _, _, name in OFFSET_BANDS:
        n = sum(1 for r in rows if r["offset_band"] == name)
        print(f"    {name:26} {n}")
    print(f"  distinct genes: {len({r['gene'] for r in rows})}")
    print("\n  top diseases represented:")
    for t, n in collections.Counter(
            r["trait_names"].split(";")[0].strip() for r in rows).most_common(8):
        print(f"    {t[:60]:60} {n}")


if __name__ == "__main__":
    main()
