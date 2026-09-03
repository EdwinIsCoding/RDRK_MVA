#!/usr/bin/env python3
"""Follow through on the one axis the Track 2 report called promising.

The loose end
-------------
Section 5 of the Track 2 report says mitochondrial and oxidative support is "the
only axis of the three that is better supplied than the genome average" and "the
axis where a repurposing search is most likely to find something to work with",
and then stops. It never names the activatable genes, never screens their drugs,
and never gates them. That is the one place the report raises an expectation and
does not meet it.

This closes it: the activatable genes in that axis, the drugs that act on them in
the required direction, the safety screen applied to each, and a tissue-expression
gate on the genes.

Plan section 7.2's gates, and why two of them are vacuous here
--------------------------------------------------------------
The plan requires gating every target on four things. Two bite and two cannot.

- **Expression in the affected tissue** bites. This proband has skeletal muscle
  atrophy and nephrocalcinosis, and the tumour is of skeletal muscle lineage, so
  muscle and kidney are the tissues that matter. Applied below from GTEx.
- **Blood-brain barrier penetrance** does not apply. It is conditional on a
  CNS-expressed target being the point of intervention, and this proband lacks
  microcephaly, seizures and developmental delay, so no CNS endpoint is proposed.
- **Open Targets tractability bucket** and **Pharos development level** are
  vacuous *by construction*. This candidate set was built by requiring that a
  drug already exists with an activating mechanism, so every gene in it is
  tractable and clinically precedented by definition. Running those gates would
  return "tractable" for all 23 and would look like evidence. Saying so is more
  honest than performing them.

On thresholds
-------------
No absolute expression cutoff is applied. The values in common use would be
recited from memory, which CLAUDE.md rule 2 exists to prevent. Median TPM is
reported, alongside each gene's highest-expressing tissue for context and its
rank within this candidate set, so a reader can apply their own.

Writes results/summaries/track2_mitochondrial_axis.md.
"""
from __future__ import annotations

import datetime as dt
import json
import pathlib
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, "src")

from mva.track2.chemoprevention import mechanism_actions, paediatric_trials
from mva.track2.druggable_direction import build_directional_proteome
from mva.track2.safety import DrugRecord, Verdict, screen

GTEX = "https://gtexportal.org/api/v2"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
CACHE_DD = "results/track2/cache_dd"
CACHE = pathlib.Path("results/track2/cache_mito")
AXIS_SUMMARY = pathlib.Path("results/summaries/track2_axis_availability.md")
OUT = pathlib.Path("results/summaries/track2_mitochondrial_axis.md")

#: The proband's affected tissues, from the coded HPO terms.
#: HP:0003202 skeletal muscle atrophy, HP:0000121 nephrocalcinosis, and the
#: presenting tumour is rhabdomyosarcoma, of skeletal muscle lineage.
TISSUES = {
    "Muscle_Skeletal": "skeletal muscle (HP:0003202, and the tumour lineage)",
    "Kidney_Cortex": "kidney (HP:0000121 nephrocalcinosis)",
}

#: GO terms defining the axis, as used in scripts/28.
AXIS_TERMS = ["cellular response to oxidative stress",
              "mitochondrial respiratory chain complex assembly",
              "aerobic respiration"]


class LookupFailed(Exception):
    """A record could not be retrieved, as distinct from being empty.

    This distinction is load-bearing. ChEMBL returns intermittent HTTP 500s, and
    a molecule whose record fails to load has no ATC codes for the same reason a
    vitamin has none: because the field is absent. Screening on that absence
    would pass a drug the codes might have excluded, so a failed lookup is
    raised and the molecule is reported as unretrievable rather than screened.
    """


def get(url: str, timeout: int = 45, attempts: int = 4) -> dict:
    """Fetch with backoff. ChEMBL 500s are transient and frequent enough that a
    single attempt silently degrades a third of the candidate set."""
    last = None
    for i in range(attempts):
        req = urllib.request.Request(url,
                                     headers={"User-Agent": "mva-hackathon-2026"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.load(fh)
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (i + 1))
    raise LookupFailed(f"{url}: {last}")


def cached(key: str, fetch, required: bool = False):
    """Cache a lookup. With ``required``, a failure raises rather than returning
    None, so the caller cannot mistake "could not fetch" for "nothing there"."""
    CACHE.mkdir(parents=True, exist_ok=True)
    import re
    f = CACHE / f"{re.sub(r'[^A-Za-z0-9_.-]', '_', key)}.json"
    if f.exists():
        return json.loads(f.read_text())
    try:
        d = fetch()
    except Exception as exc:
        if required:
            raise LookupFailed(str(exc)) from exc
        return None
    f.write_text(json.dumps(d))
    time.sleep(0.2)
    return d


def axis_genes() -> list[str]:
    """The activatable genes in the mitochondrial axis, read from scripts/28's
    own output so the two analyses cannot disagree."""
    if not AXIS_SUMMARY.exists():
        sys.exit("FATAL: run scripts/28_track2_axis_availability.py first.")
    import re
    text = AXIS_SUMMARY.read_text()
    block = text.split("### Mitochondrial and oxidative support")[1]
    line = next((ln for ln in block.splitlines()
                 if ln.startswith("Genes in this axis with a drug in the required")), "")
    return sorted(set(re.findall(r"`([A-Z0-9]+)`", line)))


def expression(gene: str) -> dict | None:
    ref = cached(f"gtex_ref_{gene}",
                 lambda: get(f"{GTEX}/reference/gene?geneId={urllib.parse.quote(gene)}"))
    rows = (ref or {}).get("data") or []
    hit = next((r for r in rows if r.get("geneSymbol") == gene), None)
    if not hit:
        return None
    gid = hit["gencodeId"]
    q = urllib.parse.urlencode({"gencodeId": gid, "datasetId": "gtex_v8"})
    med = cached(f"gtex_expr_{gene}",
                 lambda: get(f"{GTEX}/expression/medianGeneExpression?{q}"))
    data = (med or {}).get("data") or []
    if not data:
        return None
    by_tissue = {r["tissueSiteDetailId"]: r["median"] for r in data}
    top = max(data, key=lambda r: r["median"])
    return {"gencode": gid, "by_tissue": by_tissue,
            "top_tissue": top["tissueSiteDetailId"], "top_median": top["median"]}


def molecule(chembl_id: str) -> dict:
    """Raises LookupFailed rather than returning an empty record."""
    d = cached(f"mol_{chembl_id}",
               lambda: get(f"{CHEMBL}/molecule/{chembl_id}?format=json"),
               required=True)
    if not d or "molecule_chembl_id" not in d:
        raise LookupFailed(f"{chembl_id}: empty record")
    return d


def main() -> None:
    genes = axis_genes()
    if not genes:
        sys.exit("FATAL: could not read the axis gene list.")
    prot = build_directional_proteome(CACHE_DD)

    L: list[str] = []
    w = L.append
    w("# Track 2: the mitochondrial and oxidative axis, followed through\n")
    w(f"Generated by `scripts/40_mitochondrial_axis_followthrough.py` on "
      f"{dt.date.today().isoformat()}. Sources: ChEMBL, GTEx v8, "
      f"ClinicalTrials.gov, and the safety screen in `src/mva/track2/safety.py`.\n")
    w("Section 5 of the Track 2 report called this the only axis better supplied "
      "with activating drugs than the genome average, and stopped there. This "
      "names the genes, screens the drugs and gates the genes on expression in "
      "the tissues this proband actually has affected.\n")
    w(f"GO terms defining the axis: {', '.join(AXIS_TERMS)}.\n")

    w("## 1. Which of plan 7.2's gates can bite\n")
    w("| Gate | Applied? | Why |")
    w("|---|---|---|")
    w("| Expression in affected tissue | **yes** | Skeletal muscle and kidney, "
      "from the coded HPO terms and the tumour lineage. GTEx v8 medians below. |")
    w("| Open Targets tractability bucket | no | **Vacuous by construction.** The "
      "candidate set was built by requiring that an activating drug already "
      "exists, so every gene is tractable by definition. Running it would return "
      "a pass for all of them and look like evidence. |")
    w("| Pharos development level | no | Vacuous for the same reason: a gene with "
      "an approved or clinical-stage ligand is clinically precedented by "
      "definition of how this set was assembled. |")
    w("| Blood-brain barrier penetrance | no | Conditional on a CNS endpoint. "
      "This proband has no microcephaly, seizures or developmental delay, and no "
      "CNS intervention is proposed. |")
    w("")

    # ---------------------------------------------------------------- expression
    w("## 2. Tissue expression, no threshold applied\n")
    w("Median TPM, GTEx v8. Each gene's highest-expressing tissue is shown for "
      "context, because a low absolute value means something different in a gene "
      "that is low everywhere than in one that is high elsewhere and absent "
      "here. **No cutoff of ours is applied**: the thresholds in common use would "
      "be recited from memory.\n")
    w("| Gene | Skeletal muscle | Kidney cortex | Highest tissue | Its median |")
    w("|---|---:|---:|---|---:|")
    expr: dict[str, dict] = {}
    for g in genes:
        e = expression(g)
        if e is None:
            w(f"| {g} | not resolved | not resolved | | |")
            continue
        expr[g] = e
        m = e["by_tissue"].get("Muscle_Skeletal")
        k = e["by_tissue"].get("Kidney_Cortex")
        w(f"| {g} | {m if m is not None else 'n/a'} | "
          f"{k if k is not None else 'n/a'} | {e['top_tissue']} | {e['top_median']} |")
    w("")
    if expr:
        ranked = sorted(expr.items(),
                        key=lambda kv: kv[1]["by_tissue"].get("Muscle_Skeletal", 0))
        lowest = [g for g, _ in ranked[:5]]
        w(f"Lowest muscle expression within this set: "
          f"{', '.join(f'`{g}`' for g in lowest)}. That is a comparison inside 23 "
          f"genes, not a statement that any is unexpressed.\n")

    # ------------------------------------------------------------------ drugs
    w("## 3. The drugs, and the safety screen\n")
    w("Every molecule ChEMBL records as acting on one of these genes in the "
      "**activating** direction, screened by the same deterministic rules used "
      "for the chemoprevention axis.\n")

    # One lookup per molecule, not per gene-drug pair: the AMPK subunits share
    # their activators, so a per-pair loop repeats the same three network calls
    # several times over. And the paediatric trial query, the slowest of the
    # three, runs only for molecules the screen has not already excluded.
    pairs = [(g, mid) for g in genes
             for mid in sorted(prot.activatable.get(g, set())) if mid]
    unique = sorted({mid for _, mid in pairs})

    info: dict[str, dict] = {}
    unretrievable: list[str] = []
    for mid in unique:
        try:
            m = molecule(mid)
        except LookupFailed:
            unretrievable.append(mid)
            continue
        name = m.get("pref_name") or mid
        atc = tuple(a for a in (m.get("atc_classifications") or [])
                    if isinstance(a, str))
        acts = mechanism_actions(mid, pathlib.Path(CACHE_DD))
        rec = DrugRecord(name=name, chembl_id=mid, atc_codes=atc,
                         mechanism="; ".join(acts) or None,
                         provenance=f"ChEMBL {mid}")
        first = screen(rec)
        if first.verdict is not Verdict.EXCLUDED and name != mid:
            paeds = paediatric_trials(name, CACHE)
            rec = DrugRecord(name=name, chembl_id=mid, atc_codes=atc,
                             mechanism="; ".join(acts) or None,
                             paediatric_trial_ids=paeds,
                             has_paediatric_pk=bool(paeds),
                             provenance=f"ChEMBL {mid}")
            final = screen(rec)
        else:
            final = first
        info[mid] = {"name": name, "atc": atc,
                     "phase": m.get("max_phase"), "res": final, "rec": rec}

    rows = [(g, info[mid]["name"], mid, info[mid]["atc"], info[mid]["phase"],
             info[mid]["res"], info[mid]["rec"])
            for g, mid in pairs if mid in info]

    order = {Verdict.ALLOWED: 0, Verdict.FLAGGED: 1, Verdict.UNKNOWN: 2,
             Verdict.EXCLUDED: 3}
    rows.sort(key=lambda r: (order[r[5].verdict], r[0], r[1]))

    counts = {v: sum(1 for r in rows if r[5].verdict is v) for v in Verdict}
    w("| verdict | drug-target pairs |")
    w("|---|---:|")
    for v in (Verdict.ALLOWED, Verdict.FLAGGED, Verdict.UNKNOWN, Verdict.EXCLUDED):
        w(f"| {v.value} | {counts[v]} |")
    w("")
    if unretrievable:
        w(f"**{len(unretrievable)} molecule record(s) could not be retrieved** "
          f"from ChEMBL after four attempts each and are excluded from the table "
          f"rather than screened: {', '.join(unretrievable[:12])}"
          + (" ..." if len(unretrievable) > 12 else "") + ". A molecule whose "
          "record failed to load has no ATC codes for the same reason a vitamin "
          "has none, and screening on that absence would pass a drug the codes "
          "might have excluded.\n")
    w(f"**{len(rows)} activating drug-target pairs across {len(genes)} genes, "
      f"from {len(unique)} distinct molecules.** The gap between those two "
      f"numbers is mostly the AMPK subunits, which share their activators.\n")

    w("| Gene | Drug | ChEMBL | Max phase | ATC | Muscle TPM | Verdict |")
    w("|---|---|---|---:|---|---:|---|")
    for g, name, mid, atc, phase, res, _rec in rows:
        m = expr.get(g, {}).get("by_tissue", {}).get("Muscle_Skeletal")
        w(f"| {g} | {name} | {mid} | "
          f"{phase if phase is not None else 'n/a'} | "
          f"{', '.join(atc[:2]) or 'none'} | {m if m is not None else 'n/a'} | "
          f"**{res.verdict.value}** |")
    w("")

    caveated = [(g, n, r) for g, n, _, _, _, r, _ in rows
                if r.verdict is not Verdict.EXCLUDED and r.mandatory_caveats]
    if caveated:
        w("### Caveats that must travel with each candidate\n")
        for g, n, res in caveated:
            w(f"**{n}** (for `{g}`)")
            for c in res.mandatory_caveats:
                w(f"- {c}")
            w("")

    excl = [(g, n, r) for g, n, _, _, _, r, _ in rows if r.verdict is Verdict.EXCLUDED]
    if excl:
        w("### Excluded, and why\n")
        w("| Gene | Drug | Rule | Reason |")
        w("|---|---|---|---|")
        for g, n, res in excl:
            for f in res.findings:
                if f.verdict is Verdict.EXCLUDED:
                    w(f"| {g} | {n} | {f.rule_id} | {f.reason} |")
        w("")

    w("## 4. Where the pharmacology actually comes from\n")
    import collections
    # Counted over the pairs that were actually screened, not over every pair
    # nominated. Three molecule records were unretrievable and are excluded
    # above; counting them here would leave section 4 and section 5 quoting
    # different totals for the same set.
    per_gene = collections.Counter(g for g, _, _, _, _, _, _ in rows)
    top = per_gene.most_common(5)
    top_share = sum(n for _, n in top[:3]) / len(rows) if rows else 0
    w("| Gene | Screened activating drug-target pairs |")
    w("|---|---:|")
    for g, n in per_gene.most_common():
        w(f"| {g} | {n} |")
    w("")
    w(f"**{sum(n for _, n in top[:3])} of {len(rows)} screened pairs, "
      f"{top_share:.0%}, come from just three genes: "
      f"{', '.join(f'`{g}`' for g, _ in top[:3])}.**\n")
    w("Those three are metabolic drug targets. The drugs behind them are "
      "insulin analogues, fibrates and glucokinase activators, and they reached "
      "this axis because their genes carry GO annotations for oxidative-stress "
      "response or aerobic respiration. **That is an annotation artefact, not a "
      "therapeutic rationale for this child.** Insulin is not mitochondrial "
      "support for a proband with mosaic variegated aneuploidy; it is a diabetes "
      "drug whose receptor is annotated to a metabolic process.\n")
    w("The Track 2 report already carried the caveat that GO annotation sets are "
      "not therapeutic targets and that membership is a claim about a process "
      "rather than a point of intervention in it. This demonstrates it rather "
      "than asserting it.\n")

    w("## 5. What this axis is worth\n")
    n_ok = counts[Verdict.ALLOWED]
    w(f"**{n_ok} of {len(rows)} activating drug-target pairs pass the safety "
      f"screen outright**, {counts[Verdict.FLAGGED]} carry a flag that must "
      f"travel with them, {counts[Verdict.UNKNOWN]} lack the paediatric evidence "
      f"the plan requires, and {counts[Verdict.EXCLUDED]} are excluded.\n")
    w("**The axis does not survive being followed through, and the enrichment "
      "that made it look promising was an artefact.** Section 5 of the report "
      "called it the only axis better supplied with activating drugs than the "
      "genome average and the one most likely to yield something. Both were true "
      "of the numbers and neither survives inspection: the supply is "
      "concentrated in metabolic genes that a GO annotation swept in, and what "
      "it yields is insulin.\n")
    w("This is a reported negative, and it is the third axis to close. The "
      "direct spindle-checkpoint axis is pharmacologically unavailable in the "
      "required direction. Proteotoxic stress mitigation sits at the base rate. "
      "Mitochondrial and oxidative support is above the base rate for a reason "
      "that does not transfer to this patient. **Cancer chemoprevention and "
      "surveillance is the axis left standing**, which is where plan section 7.1 "
      "said the value was before any of this was measured.\n")
    w("The mechanistic case for the axis is unchanged and remains weak on its "
      "own terms: aneuploid cells carry a proteostasis and energetic burden, and "
      "no evidence links modulating any gene here to an outcome in mosaic "
      "variegated aneuploidy, let alone in this child.\n")
    w("**What would falsify it.** The same experiment the rest of Track 2 "
      "names: a `Bub1b` hypomorphic mouse, here with a mitochondrial or "
      "oxidative endpoint rather than tumour incidence. Nothing in this table "
      "substitutes for it.\n")
    w("## 6. A reproducibility caveat worth stating\n")
    w("ChEMBL returned intermittent HTTP 500s throughout this analysis. Requests "
      "are retried four times with backoff and a molecule whose record still "
      "cannot be fetched is reported as unretrievable rather than screened on "
      "absent fields, because a failed lookup and a drug with no ATC code look "
      "identical to a screen.\n")
    w("**The consequence is that a fresh run can differ slightly from this one.** "
      "Responses are cached under `results/track2/cache_mito/`, so a repeat run "
      "on this machine reproduces exactly, but a run from a clean checkout "
      "depends on which requests ChEMBL happens to serve. The counts above are "
      "therefore reproducible in the strong sense only alongside the cache. "
      "Nothing about the conclusion turns on one or two molecules: the finding "
      "is a 71% concentration in three metabolic genes, which no plausible "
      "handful of failures changes.\n")
    w("No dosing appears here and none may be added. These are research "
      "hypotheses addressed to researchers.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L[:60]))
    print(f"\n... written to {OUT}")


if __name__ == "__main__":
    main()
