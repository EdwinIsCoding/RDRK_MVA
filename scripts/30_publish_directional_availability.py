#!/usr/bin/env python3
"""Publish the directional availability table as a reusable resource.

Why publish it
--------------
Our Track 2 result is a negative: the direct spindle-checkpoint axis requires
activating ten targets and no activating drug exists for any of them. Measuring
the base rate to keep that claim honest produced something more generally
useful than the claim itself, which is the answer to "can this protein be pushed
in the direction I need?" for the whole druggable proteome.

Any team proposing to activate a target can check it here in one lookup, instead
of discovering after the fact that every ligand for their target is an
inhibitor. That is the reusable part of a result that is otherwise specific to
one child.

What is published, and why only this
------------------------------------
ChEMBL is released under CC BY-SA 3.0 and this project's outputs are CC BY 4.0.
Redistributing ChEMBL content under a licence that drops share-alike would be a
licence question we are not in a position to answer. So this publishes **gene
symbols and counts**, which are facts derived from the database rather than a
copy of it, with attribution. Anyone needing the underlying molecule records
should take them from ChEMBL directly, where they are properly licensed.

Writes resources/directional_availability/.
"""
from __future__ import annotations

import datetime as dt
import json
import pathlib
import sys

sys.path.insert(0, "src")

from mva.track2.druggable_direction import build_directional_proteome

CACHE = "results/track2/cache_dd"
OUT = pathlib.Path("resources/directional_availability")
DB_VERSIONS = pathlib.Path("config/db_versions.yaml")
HGNC = pathlib.Path("refs/hgnc_complete_set.txt")


def n_protein_coding() -> int:
    import csv
    if not HGNC.exists():
        raise SystemExit(
            f"FATAL: {HGNC} is absent. It supplies the denominator for every "
            f"rate in this analysis, and returning zero here would silently "
            f"turn the base rate into 0.00% and make the argument look absurd "
            f"rather than fail. Run `make downloads-track2` first."
        )
    with HGNC.open(newline="") as fh:
        return sum(1 for r in csv.DictReader(fh, delimiter="\t")
                   if r.get("locus_group") == "protein-coding gene"
                   and r.get("status") == "Approved")


def chembl_version() -> str:
    """Read the ChEMBL release from config/db_versions.yaml, never from memory."""
    try:
        import yaml
        d = yaml.safe_load(DB_VERSIONS.read_text())
        c = d["track2"]["chembl"]
        return f"{c['db_version']} (released {c['release_date']})"
    except Exception:
        return "TODO(source)"


def main() -> None:
    prot = build_directional_proteome(CACHE)
    OUT.mkdir(parents=True, exist_ok=True)
    n_genes = n_protein_coding()
    version = chembl_version()
    today = dt.date.today().isoformat()

    for label, table in (("activatable", prot.activatable),
                         ("inhibitable", prot.inhibitable)):
        f = OUT / f"{label}_genes.tsv"
        with f.open("w") as fh:
            fh.write("gene_symbol\tn_drug_mechanisms\n")
            for gene in sorted(table):
                fh.write(f"{gene}\t{len(table[gene])}\n")
        print(f"wrote {f} ({len(table):,} genes)")

    both = sorted(set(prot.activatable) & set(prot.inhibitable))
    (OUT / "both_directions_genes.tsv").write_text(
        "gene_symbol\n" + "\n".join(both) + "\n")

    summary = {
        "generated": today,
        "chembl_version": version,
        "activatable_genes": len(prot.activatable),
        "inhibitable_genes": len(prot.inhibitable),
        "genes_with_both_directions": len(both),
        "approved_protein_coding_genes_hgnc": n_genes,
        "activating_mechanism_records": prot.n_activating_mechanisms,
        "inhibiting_mechanism_records": prot.n_inhibiting_mechanisms,
        "activation_rate": round(len(prot.activatable) / n_genes, 5) if n_genes else None,
        "inhibition_rate": round(len(prot.inhibitable) / n_genes, 5) if n_genes else None,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    act_rate = summary["activation_rate"] or 0
    inh_rate = summary["inhibition_rate"] or 0
    (OUT / "README.md").write_text(f"""# Directional drug availability across the human proteome

Generated {today} by `scripts/30_publish_directional_availability.py`.

## What this answers

**Can this protein be pushed in the direction my hypothesis needs?**

A repurposing hypothesis that needs a target *activated* is worth very little if
every ligand for that target is an inhibitor. That is easy to discover late and
cheap to check early, and this table is the check.

| | Genes | Share of protein-coding genome |
|---|---:|---:|
| An **activating** drug mechanism exists | {len(prot.activatable):,} | {act_rate:.2%} |
| An **inhibiting** drug mechanism exists | {len(prot.inhibitable):,} | {inh_rate:.2%} |
| Both directions available | {len(both):,} | |
| Approved protein-coding genes (HGNC) | {n_genes:,} | |

**Inhibition is roughly {inh_rate / act_rate:.1f} times more available than
activation.** Any claim that a particular target set is "undruggable in the
required direction" has to be read against that, which is the reason we built
this: it weakened our own headline finding, and we would rather state the
weakened version than the flattering one.

## Files

| File | Contents |
|---|---|
| `activatable_genes.tsv` | Gene symbol, and how many drug mechanisms act to increase its activity |
| `inhibitable_genes.tsv` | The same for mechanisms that decrease it |
| `both_directions_genes.tsv` | Genes reachable either way |
| `summary.json` | The counts above, machine-readable |

## Method

Every ChEMBL drug-mechanism record with an activating action type (agonist,
partial agonist, activator, opener, positive allosteric or positive modulator,
stabiliser) or an inhibiting one (inhibitor, antagonist, blocker, negative
allosteric or negative modulator, disrupting agent, inverse agonist,
downregulator, sequestering agent) was retrieved, and its target resolved to
gene symbols. Built from {prot.n_activating_mechanisms:,} activating and
{prot.n_inhibiting_mechanisms:,} inhibiting mechanism records.

Source: {version}.

## Limits, which matter before you use this

- **Availability is not suitability.** A gene here has a drug that acts in that
  direction. Nothing is said about whether the drug is safe, reaches the tissue,
  or would help any particular patient.
- **Action-type annotation is ChEMBL's**, and it is incomplete for some targets
  and absent for others. A gene missing from both tables may have no drug, or
  may have drugs whose action type was never curated.
- **Gene symbol resolution is via ChEMBL target components.** Multi-protein
  targets and complexes contribute every component gene, so a symbol here may be
  reachable only as part of a complex.
- **This is a snapshot.** ChEMBL is versioned and moves; regenerate rather than
  citing these counts years later.

## Licence and attribution

Derived from **ChEMBL** ({version}), which is released under CC BY-SA 3.0. We
publish gene symbols and counts, which are facts derived from the database
rather than a copy of it, with attribution. If you need the underlying molecule
and mechanism records, take them from ChEMBL directly where they are properly
licensed: <https://www.ebi.ac.uk/chembl/>.

Gene counts are against HGNC approved protein-coding genes.

This file and the tables beside it are released under CC BY 4.0 as part of the
MVA Hackathon 2026 submission.
""")
    print(f"wrote {OUT}/README.md and summary.json")


if __name__ == "__main__":
    main()
