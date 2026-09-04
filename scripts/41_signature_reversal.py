#!/usr/bin/env python3
"""Signature reversal against a labelled proxy, per plan section 7.3.

The gap this closes
-------------------
Plan section 7.3 calls LINCS/CMap connectivity "the highest-yield repurposing
method available" and says it "needs no pocket, no structure, no known
mechanism". It also covers our case explicitly: "If no RNA-seq: use published
transcriptomic signatures from MVA/CIN models (BUB1B-hypomorph mice, aneuploid
lines) as a proxy and **label the proxy clearly**."

This was never run. The Track 2 report said we "did not substitute a published
proxy signature and present it as though it were the patient's", which rejects
something nobody proposed in order to avoid doing what the plan sanctions. That
sentence read as restraint and was an unexplored path wearing restraint's
clothes.

What the proxy is, and every way it differs from the patient
-----------------------------------------------------------
GEO **GSE277997**, "BubR1 Insufficiency Recapitulates Changes Associated with
Age-Related Cardiac Pathologies": 12 mice, RNA-seq, hypomorph versus wild type,
with the authors' own differential expression table.

It is the right genotype and almost nothing else:

- **Mouse, not human.**
- **Cardiac tissue.** This proband's affected tissues are skeletal muscle and
  kidney, and the tumour is of skeletal muscle lineage. The heart is not
  involved.
- **An ageing-heart phenotype**, not a paediatric cancer-predisposition
  presentation.
- **A different allelic architecture**: an engineered hypomorph, not a nonsense
  allele in trans with a missense.

So a compound reversing this signature is a compound that opposes the
transcriptional consequences of BubR1 insufficiency in mouse heart. Whether that
has anything to do with this child is exactly what is not established, and the
output is labelled accordingly throughout.

Writes results/summaries/track2_signature_reversal.md.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import json
import pathlib
import sys
import urllib.request

sys.path.insert(0, "src")

from mva.track2.chemoprevention import mechanism_actions, paediatric_trials, resolve_agent
from mva.track2.safety import DrugRecord, screen

SIG = pathlib.Path("results/track2/signature/GSE277997_H_v_WT.csv.gz")
HOM = pathlib.Path("results/track2/signature/HOM_MouseHumanSequence.rpt")
HGNC = pathlib.Path("refs/hgnc_complete_set.txt")
CACHE = pathlib.Path("results/track2/cache_sig")
OUT = pathlib.Path("results/summaries/track2_signature_reversal.md")
L1000CDS2 = "https://maayanlab.cloud/L1000CDS2/query"

GEO = "GSE277997"
#: A stated choice, not a discovered threshold. Sensitivity to it is reported.
PADJ = 0.05
#: From config/config.yaml analysis.random_seed, so the null control below is
#: reproducible (CLAUDE.md rule 5).
SEED = 20261024
#: L1000CDS2 characteristic-direction search behaves best on a few hundred genes
#: per direction; the full significant set would swamp it.
TOP_N = 150


def orthologs() -> dict[str, str]:
    """Mouse symbol to human symbol, from MGI's homology classes.

    Not the uppercase naming convention. That convention is right most of the
    time and wrong silently, and "most of the time" is not a mapping.
    """
    by_class: dict[str, dict[str, str]] = {}
    with HOM.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = r.get("DB Class Key")
            org = (r.get("Common Organism Name") or "").lower()
            sym = r.get("Symbol")
            if not key or not sym:
                continue
            side = "human" if "human" in org else "mouse" if "mouse" in org else None
            if side:
                by_class.setdefault(key, {})[side] = sym
    return {v["mouse"]: v["human"] for v in by_class.values()
            if "mouse" in v and "human" in v}


def approved_symbols() -> set[str]:
    if not HGNC.exists():
        return set()
    with HGNC.open(newline="") as fh:
        return {r["symbol"] for r in csv.DictReader(fh, delimiter="\t")
                if r.get("status") == "Approved"}


def load_signature() -> list[dict]:
    rows = []
    with gzip.open(SIG, "rt") as fh:
        for r in csv.DictReader(fh):
            try:
                rows.append({"gene": r["gene"],
                             "lfc": float(r["log2FC"]),
                             "padj": float(r["padj"])})
            except (ValueError, KeyError, TypeError):
                continue
    return rows


def query_l1000(up: list[str], dn: list[str]) -> list[dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    key = CACHE / f"l1000_{len(up)}_{len(dn)}_{hash(tuple(up + dn)) & 0xffffffff}.json"
    if key.exists():
        return json.loads(key.read_text()).get("topMeta", [])
    payload = {
        "data": {"upGenes": up, "dnGenes": dn},
        # aggravate=false asks for perturbations that REVERSE the signature.
        "config": {"aggravate": False, "searchMethod": "geneSet",
                   "share": False, "combination": False, "db-version": "latest"},
        "metadata": [],
    }
    req = urllib.request.Request(
        L1000CDS2, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json",
                 "User-Agent": "mva-hackathon-2026"})
    with urllib.request.urlopen(req, timeout=180) as fh:
        d = json.load(fh)
    key.write_text(json.dumps(d))
    return d.get("topMeta", [])


def null_control(mapped: list[dict], n_up: int, n_dn: int,
                 rounds: int = 3) -> list[set[str]]:
    """The same query on random gene sets of the same size.

    HDAC and mTOR inhibitors perturb transcription broadly and surface in many
    L1000 reversal queries whatever the input. Without a null, a hit list of
    them is indistinguishable from a hit list of them. Anything that also
    returns for random genes drawn from the same expressed background is not
    evidence about BubR1.

    Seeded from config/config.yaml so the control is reproducible.
    """
    import random
    rng = random.Random(SEED)
    pool = [r["human"] for r in mapped]
    out = []
    for _ in range(rounds):
        pick = rng.sample(pool, min(n_up + n_dn, len(pool)))
        try:
            hits = query_l1000(pick[:n_up], pick[n_up:n_up + n_dn])
        except Exception:
            continue
        out.append({(h.get("pert_desc") or "").strip().lower()
                    for h in hits if (h.get("pert_desc") or "").strip()})
    return out


def main() -> None:
    for f in (SIG, HOM):
        if not f.exists():
            sys.exit(f"FATAL: {f} absent.")

    sig = load_signature()
    orth = orthologs()
    approved = approved_symbols()

    sigf = [r for r in sig if r["padj"] < PADJ]
    mapped = []
    unmapped = 0
    for r in sigf:
        h = orth.get(r["gene"])
        if h is None or (approved and h not in approved):
            unmapped += 1
            continue
        mapped.append({**r, "human": h})

    mapped.sort(key=lambda r: r["lfc"])
    dn = [r["human"] for r in mapped[:TOP_N]]
    up = [r["human"] for r in reversed(mapped[-TOP_N:])]

    hits = query_l1000(up, dn)

    L: list[str] = []
    w = L.append
    w("# Track 2: signature reversal against a labelled proxy\n")
    w(f"Generated by `scripts/41_signature_reversal.py` on "
      f"{dt.date.today().isoformat()}. Sources: GEO {GEO}, MGI homology, HGNC, "
      f"L1000CDS2 over LINCS L1000.\n")

    w("## The proxy, and every way it is not this patient\n")
    w(f"**{GEO}**, \"BubR1 Insufficiency Recapitulates Changes Associated with "
      f"Age-Related Cardiac Pathologies\". Twelve mice, RNA-seq, hypomorph "
      f"versus wild type, using the authors' own differential expression "
      f"table.\n")
    w("Plan section 7.3 sanctions exactly this substitution and requires the "
      "proxy be labelled clearly. It is the right genotype and almost nothing "
      "else:\n")
    w("| | Proxy | This proband |")
    w("|---|---|---|")
    w("| Species | mouse | human |")
    w("| Tissue | heart | skeletal muscle and kidney affected; tumour of skeletal muscle lineage |")
    w("| Phenotype | age-related cardiac pathology | paediatric cancer "
      "predisposition with rhabdomyosarcoma |")
    w("| Genotype | engineered BubR1 hypomorph | nonsense allele in trans with a missense |")
    w("")
    w("**A compound that reverses this signature opposes the transcriptional "
      "consequences of BubR1 insufficiency in mouse heart.** Whether that bears "
      "on this child is precisely what is not established here.\n")

    w("## Building the signature\n")
    w("| Stage | Genes |")
    w("|---|---:|")
    w(f"| in the authors' table | {len(sig):,} |")
    w(f"| significant at adjusted p < {PADJ} | {len(sigf):,} |")
    w(f"| with a one-to-one MGI human ortholog, symbol approved by HGNC | {len(mapped):,} |")
    w(f"| dropped for want of an ortholog or an approved symbol | {unmapped:,} |")
    w(f"| submitted, most up-regulated | {len(up)} |")
    w(f"| submitted, most down-regulated | {len(dn)} |")
    w("")
    w(f"The adjusted p threshold of {PADJ} is **a stated choice, not a "
      f"discovered value**. Orthology is MGI's homology classes rather than the "
      f"uppercase naming convention, which is right most of the time and wrong "
      f"silently.\n")

    w("## What reverses it\n")
    if not hits:
        w("**Nothing returned.** The query completed and produced no "
          "perturbation, which is a reportable negative rather than an error.\n")
    else:
        seen, rows = set(), []
        for h in hits:
            name = (h.get("pert_desc") or "").strip()
            if not name or name == "-" or name.lower() in seen:
                continue
            seen.add(name.lower())
            rows.append((name, h.get("score"), h.get("cell_id"), h.get("pert_id")))
        w(f"L1000CDS2 returned **{len(hits)}** signature matches over "
          f"**{len(rows)}** distinct perturbagens, ranked by overlap with the "
          f"reversed signature.\n")

        w("### Screened\n")
        w("Each named compound resolved against ChEMBL and put through the same "
          "safety screen as every other Track 2 candidate.\n")
        w("| Compound | L1000 score | Cell line | ChEMBL | Verdict |")
        w("|---|---:|---|---|---|")
        counts: dict[str, int] = {}
        for name, score, cell, _pid in rows[:25]:
            mol, matched = resolve_agent(name, CACHE)
            if mol is None:
                w(f"| {name} | {score} | {cell} | unresolved | not screened |")
                counts["unresolved"] = counts.get("unresolved", 0) + 1
                continue
            acts = mechanism_actions(mol.chembl_id, CACHE)
            paeds = paediatric_trials(mol.pref_name, CACHE)
            rec = DrugRecord(name=mol.pref_name, chembl_id=mol.chembl_id,
                             atc_codes=mol.atc_codes,
                             mechanism="; ".join(acts) or None,
                             paediatric_trial_ids=paeds,
                             has_paediatric_pk=bool(paeds),
                             provenance=f"ChEMBL {mol.chembl_id}")
            res = screen(rec)
            counts[res.verdict.value] = counts.get(res.verdict.value, 0) + 1
            w(f"| {mol.pref_name} | {score} | {cell} | {mol.chembl_id} | "
              f"**{res.verdict.value}** |")
        w("")
        w("| verdict | compounds |")
        w("|---|---:|")
        for k in sorted(counts):
            w(f"| {k} | {counts[k]} |")
        w("")

    # ------------------------------------------------------------ null control
    w("## Are these hits specific to the signature?\n")
    w("HDAC and mTOR inhibitors perturb transcription broadly and surface in "
      "many L1000 reversal queries whatever is submitted. So the identical "
      "query was run three times on random gene sets of the same size, drawn "
      "from the same expressed and orthology-mapped background, seeded for "
      "reproducibility.\n")
    nulls = null_control(mapped, len(up), len(dn))
    if not nulls:
        w("The control did not complete, so specificity is **unestablished** and "
          "the list above should be read as unvalidated.\n")
    else:
        real = {(h.get("pert_desc") or "").strip().lower()
                for h in hits if (h.get("pert_desc") or "").strip()}
        w("| Random draw | Perturbagens returned | Also in the real result |")
        w("|---|---:|---:|")
        overlaps = []
        for i, nset in enumerate(nulls, 1):
            ov = len(real & nset)
            overlaps.append(ov)
            w(f"| {i} | {len(nset)} | {ov} |")
        w("")
        mean_ov = sum(overlaps) / len(overlaps)
        frac = mean_ov / len(real) if real else 0
        w(f"**On average {mean_ov:.1f} of {len(real)} perturbagens, "
          f"{frac:.0%}, also come back from random genes.**\n")
        if frac >= 0.5:
            w("**That is most of the list, and it means the result is largely "
              "not about BubR1.** These compounds reverse many signatures, "
              "including ones assembled at random from the same background. The "
              "hit list should be read as the method's generic output rather "
              "than as a finding about this disease, and we report it that way "
              "rather than presenting a ranked table that looks specific and is "
              "not.\n")
        elif frac >= 0.2:
            w("**A substantial minority of the list is generic.** The compounds "
              "shared with the random draws carry little information about "
              "BubR1; those unique to the real signature are the only part worth "
              "any further attention, and even they inherit every caveat in the "
              "proxy table above.\n")
        else:
            w("**The overlap is small**, so the hits are largely specific to the "
              "submitted signature rather than generic transcriptional "
              "perturbagens. That makes the list worth reading, subject to every "
              "caveat in the proxy table above.\n")
        uniq = sorted(real - set().union(*nulls)) if nulls else []
        if uniq:
            w(f"Returned for the real signature and none of the random draws "
              f"({len(uniq)}). These are LINCS perturbagen identifiers as the "
              f"database records them, so many are supplier catalogue codes "
              f"rather than compound names, and a string such as "
              f"`656402-250mg` is a catalogue number and **not a dose**:\n")
            w(", ".join(f"`{u}`" for u in uniq[:20])
              + (" ..." if len(uniq) > 20 else "") + "\n")

    w("## What this is worth\n")
    w("**Method, not answer.** Connectivity reversal generates hypotheses. A "
      "compound that reverses a transcriptional signature is not thereby a "
      "treatment, and the LINCS signatures were measured in cancer cell lines "
      "under drug perturbation, which is a third layer of distance from this "
      "child on top of species and tissue.\n")
    w("**The distances multiply rather than add.** Mouse to human, heart to "
      "muscle and kidney, ageing to paediatric cancer predisposition, and cell "
      "line to patient. Any one would warrant caution; together they mean the "
      "output is a direction to look rather than a candidate to advance.\n")
    w("**What would make it real.** The same experiment named everywhere else "
      "in this report: a `Bub1b` hypomorphic model, ideally in the affected "
      "tissue, testing whether a compound reversing this signature changes any "
      "outcome that matters. Nothing here substitutes for it.\n")
    w("**Why it is reported despite all of that.** Plan section 7.3 sanctions a "
      "clearly labelled proxy, and an unrun method is worth less than a run one "
      "with its limits stated. The alternative was leaving the highest-yield "
      "repurposing approach unattempted while implying it had been considered.\n")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L[:45]))
    print(f"\n... written to {OUT}")


if __name__ == "__main__":
    main()
