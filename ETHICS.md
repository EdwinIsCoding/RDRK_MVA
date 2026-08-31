# ETHICS.md

This repository analyses whole genome sequencing and clinical phenotype data
from a living child with a suspected rare disease, shared by a family who
consented in order to advance rare disease research. This document records what
we did, what we deliberately did not do, and where a judgement call was made
that a reader might reasonably want to reverse.

Written at Phase 0, 31 August 2026. To be updated at each STOP checkpoint.

---

## 1. Consent scope as we understand it

The phenotype document supplied by the organisers states that it contains real
clinical information from a living patient and family who have consented to
share their story to advance rare disease research, and that any publication or
communication arising from the challenge must not include information that could
re-identify the patient or family beyond what is already publicly available
through the family's own blog posts.

We read that as three obligations:

1. Use the data for rare disease research within this challenge.
2. Publish nothing that could re-identify the patient or family.
3. Treat the boundary of what is already public as a limit we do not push
   against, not a licence to aggregate.

Where the hackathon rules are more specific than this, the rules govern. They
are transcribed in `RULES.md`.

## 2. What we did not do

- **We did not commit patient data to version control.** `.gitignore` excluding
  `data/` and every genomic file extension was written before the repository's
  first commit, when it had no commits at all. A pre-commit hook
  (`.githooks/pre-commit`) hard-fails on anything under `data/`, on genomic
  extensions anywhere in the tree, and on files above 5 MB. Git history is
  effectively permanent and this repository is intended to be public, so the
  guard is placed before the mistake rather than after it.

- **We did not put any patient genomic data into a hosted LLM context.** No
  variant, coordinate, allele, genotype or read from the proband has been read
  into an LLM. Every number in `DATA_CARD.md` is an aggregate produced by a
  script in `scripts/`, and the scripts are what touched the data. This
  constraint is written into `CLAUDE.md` as rule 1 so it binds future sessions.

- **We did not look for the family's blog.** The phenotype document mentions
  that the family have published their own posts. Locating them would enrich the
  phenotype and would also re-identify the family, and cross-referencing a
  de-identified genome against self-published family material is precisely the
  aggregation the consent language guards against. We did not search for it and
  the repository contains no reference to it. Anyone continuing this work should
  not either.

- **We did not send patient data to any third-party API.** Two external calls
  were made during Phase 0, both to the Ensembl REST API, both transmitting only
  gene symbols (`BUB1B`, `CEP57`, `TRIP13`, `BUB1`, `BUB3`, `CEP192`, `SMC5`,
  `CENATAC`, `CEP57L1`) to retrieve public reference coordinates. No patient
  information was included.

- **We did not report a lead we could not calibrate.** Two candidate findings
  were raised and closed during Phase 0 (`DATA_CARD.md` section 5). Reporting
  either would have been easy and would have been wrong.

## 3. A judgement call: we read the clinical phenotype document

**What we did.** `data/Challenge_Clinical_Phenotype_1.docx` was converted to text
by `scripts/02_extract_clinical.sh` and read in full, including by the LLM
assisting this work.

**Why.** The plan's own Phase 0 requires determining whether HPO terms are
available, what karyotype and aneuploidy percentages were reported, how many
tissues exist, and whether a candidate heterozygous variant has already been
named. Plan section 12 requires reporting the last two at STOP #1. None of that
is answerable without reading the document. The document is patient data but it
is not genomic data: it contains eight HPO terms, a gestational age, a birth
weight and a family history, and no sequence, coordinate or genotype. It is the
material the challenge distributes specifically to be reasoned over.

**The tension we are not going to paper over.** The scaffold's own default is
that patient data does not enter a hosted LLM API. This is an exception to that
default, made deliberately and recorded here rather than taken silently. A
reasonable person could set the line differently.

**How to reverse it.** If the hackathon rules or the team's own judgement place
the phenotype document off-limits to a hosted model, then: delete
`results/recon/Challenge_Clinical_Phenotype_1.docx.md`, remove the `hpo_terms`
and `phenotype` blocks from `config/config.yaml` and section 4 of
`DATA_CARD.md`, and drive Exomiser and LIRICAL from an HPO term list prepared by
a human. Nothing downstream depends on the free text, only on the eight coded
terms, so the cost of reversing is low. It gets higher the longer it is left.

## 4. What the outputs are, and are not

Everything this repository produces is a **research hypothesis**, generated by
computational prediction, for evaluation by qualified researchers and
clinicians.

- **Nothing here is a diagnosis.** A ranked candidate variant is a hypothesis
  requiring orthogonal confirmation. Every Track 1 candidate ships with the
  specific experiment that would falsify it, because a candidate without a
  falsification route is not a scientific claim.
- **Nothing here is clinical advice, and nothing here is a treatment
  recommendation.** Track 2 nominates repurposing hypotheses with mechanistic
  rationale and safety contraindications. No dosing, dose range or safety margin
  appears anywhere in this repository, and none may be added. That prohibition is
  rule 3 in `CLAUDE.md`.
- **Confidence is stated, including when it is low.** Where evidence is weak we
  say so rather than rounding up. The known limitations are enumerated in
  `config/config.yaml` under `limitations` and carried verbatim into the
  submission.
- **We do not have RNA-seq**, so every splicing result is a prediction and never
  an observation. Because that distinction is easy to lose between analysis and
  write-up, it is recorded as an anti-pattern in `CLAUDE.md`.

## 5. A note on the audience

The family may read the output. That has two consequences we intend to honour.

First, an overconfident result is not a neutral error here. A polished dossier
confers authority the evidence does not support, and a family reading a
confident wrong answer pays a cost that a leaderboard does not capture. Ranked
hypotheses with honest uncertainty are the right output even where a single
confident answer would score better.

Second, a negative result is worth stating plainly. "We excluded this, here is
how" is useful to a family and to the field. Several such statements are already
in `DATA_CARD.md` and `RECON.md`.

## 6. Licence

All outputs of this work are released under CC BY 4.0, per the challenge terms.
This does not extend to the input data, which is not ours to license and is not
redistributed.
