# Directional drug availability across the human proteome

Generated 2026-09-03 by `scripts/30_publish_directional_availability.py`.

## What this answers

**Can this protein be pushed in the direction my hypothesis needs?**

A repurposing hypothesis that needs a target *activated* is worth very little if
every ligand for that target is an inhibitor. That is easy to discover late and
cheap to check early, and this table is the check.

| | Genes | Share of protein-coding genome |
|---|---:|---:|
| An **activating** drug mechanism exists | 359 | 1.86% |
| An **inhibiting** drug mechanism exists | 1,541 | 7.99% |
| Both directions available | 214 | |
| Approved protein-coding genes (HGNC) | 19,297 | |

**Inhibition is roughly 4.3 times more available than
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
gene symbols. Built from 1,269 activating and
4,897 inhibiting mechanism records.

Source: ChEMBL_37 (released 2026-05-01).

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

Derived from **ChEMBL** (ChEMBL_37 (released 2026-05-01)), which is released under CC BY-SA 3.0. We
publish gene symbols and counts, which are facts derived from the database
rather than a copy of it, with attribution. If you need the underlying molecule
and mechanism records, take them from ChEMBL directly where they are properly
licensed: <https://www.ebi.ac.uk/chembl/>.

Gene counts are against HGNC approved protein-coding genes.

This file and the tables beside it are released under CC BY 4.0 as part of the
MVA Hackathon 2026 submission.
