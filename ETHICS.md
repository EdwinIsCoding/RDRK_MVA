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

**Status: confirmed by the project owner on 31 August 2026.** The decision was
put to them explicitly, with the reversal procedure below, and they elected to
keep it. The exception does not generalise: it covers the clinical phenotype
document only, and CLAUDE.md rule 1 continues to bar patient *genomic* data from
an LLM context. Analyses of the callset are performed by scripts, which write
aggregate summaries; the summaries are what is read.

**How to reverse it.** If the hackathon rules or the team's own judgement place
the phenotype document off-limits to a hosted model, then: delete
`results/recon/Challenge_Clinical_Phenotype_1.docx.md`, remove the `hpo_terms`
and `phenotype` blocks from `config/config.yaml` and section 4 of
`DATA_CARD.md`, and drive Exomiser and LIRICAL from an HPO term list prepared by
a human. Nothing downstream depends on the free text, only on the eight coded
terms, so the cost of reversing is low. It gets higher the longer it is left.

## 3a. LLM use, measured against the organisers' own guidance

The organisers published guidance on LLM use in the Community discussion linked
from the Space, and updated the submission requirements on 28 August 2026. LLM
use is **permitted**, with Claude named explicitly, subject to two mandatory
conditions and a disclosure requirement.

### The two conditions

1. **No training on inputs or outputs**, and no rights taken by the provider in
   either.
2. **Retention limited in time and purpose.** Zero retention is not required;
   short-lived logs for abuse monitoring, debugging or service quality are
   acceptable.

The guidance frames this as whether the provider acts as a *processor* (a tool)
or a *recipient* (gains rights in the data). Only the former is acceptable.

Two cautions are attached: disable training settings and **do not rate outputs**,
because some providers use content you give feedback on or that is flagged for
safety review; and check the terms attached to any free credits, which may
override account defaults.

### Resolved, and it is a deviation that must be disclosed

The account is a **Claude Max subscription**, which is governed by Anthropic's
consumer terms rather than the commercial terms in the organisers' example.
The "help improve Claude" setting, which permits training on conversation
content, was **ON for the whole of this work until 31 August 2026**, when the
owner disabled it on discovering the requirement.

Measured against the organisers' first mandatory condition, "No Training on your
inputs or outputs, and no rights taken by the provider in either", that period
does not comply. Their framing is whether the provider acts as a *processor*, a
tool, or a *recipient* that gains rights in the data. With training enabled, a
consumer plan is closer to the latter, which is the case they said must not
arise.

This is recorded here rather than quietly corrected, for three reasons: the
Hackathon Rules require suspected unauthorised disclosure to be reported to Sage
Bionetworks' Privacy and Compliance Office; the methods description must carry an
accurate LLM disclosure in both tracks; and an inconsistency discovered by a
judge would be far more damaging than the deviation itself.

**Disclosure line for both tracks, to be used verbatim:**

> Anthropic Claude (Claude Code), Max subscription, consumer terms. The
> "help improve Claude" training setting was enabled until 31 August 2026 and
> disabled thereafter. No patient genomic data (VCF, BAM, variant tables or
> pasted variant blocks) entered the model context at any point; see the audit
> below. Clinical phenotype terms and aggregate callset statistics did.

**TODO(owner):** report to Sage Bionetworks' Privacy and Compliance Office via
Sage's Help Center, and ask Anthropic whether already-collected conversation
content can be deleted from training pipelines under the consumer terms. Neither
can be actioned from inside this repository.

### What actually reached the model, audited

The guidance distinguishes what must be deleted from what may be kept.

| Their category | What we sent |
|---|---|
| **Delete:** VCF/BAM files | Never in context. Read only by scripts. |
| **Delete:** variant tables with genotypes | Never in context. Written to gitignored files under `results/`, never read back. |
| **Delete:** "prompts containing pasted variant blocks" | **None.** This is precisely what CLAUDE.md rule 1 prohibits, and it held. |
| **Delete:** model weights trained on raw genomic data | None trained. |
| **Keep:** HPO terms | Yes, the eight coded terms. Explicitly on the keep list. |
| **Keep:** gene rankings, candidate rankings | Yes, gene-level and aggregate only. |
| **Keep:** code, reports without raw genomic data | Yes. |

Two things did reach the model that are patient-derived and worth naming rather
than glossing:

- **The clinical phenotype document**, in full: eight HPO terms, gestational age,
  approximate birth weight, and the family history of recurrent miscarriage. This
  was a deliberate decision, recorded in section 3 above and confirmed by the
  project owner. It is clinical rather than genomic, and HPO terms sit on the
  organisers' keep list.
- **A small number of runs-of-homozygosity interval coordinates**, five of them,
  printed while summarising Arm F. These are patient-derived genomic intervals.
  They are not variants, carry no genotypes and no alleles, and are not a
  "variant block" in the sense the guidance prohibits. Recorded here for
  completeness rather than because it is thought to be a breach.

The pseudonymous sample identifier `WGS_EX2312012`, assigned by the organisers,
also appears throughout.

### Assessment

The architecture adopted at Phase 0, in which scripts touch `data/` and the
model reads only aggregate summaries, was chosen before this guidance was read
and turns out to match its central requirement almost exactly. The one category
the guidance singles out, prompts containing pasted variant blocks, is the one
CLAUDE.md rule 1 was written to prevent.

What remains outstanding is not a property of this repository: it is whether the
account's plan disables training on customer content. That is the owner's to
verify.

## 3b. Deletion obligation

The Hackathon Rules require **all data to be deleted at the conclusion of the
Hackathon**. The guidance scopes this to "what you control: your machines, cloud
instances, notebooks, repositories, storage"; provider logs are out of scope.

This needs scheduling rather than remembering. Concretely, on conclusion:

- delete `data/` in full: the VCF, its index, all eight FASTQ files, and the
  clinical phenotype document
- delete `results/` in full, which holds the extracted phenotype text, the
  variant-level shortlists, the ROH and mtDNA raw tables, and the recon outputs
- delete the derived FASTA and gnomAD slices under `refs/` that were cached for
  this analysis, though these are public reference data and not patient data
- retain, per the guidance: candidate variant rankings, HPO terms, gene
  rankings, code and reports, none of which contain raw genomic data

The repository as tracked in git already contains none of the delete-list
material, which the pre-commit hook enforces.

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
