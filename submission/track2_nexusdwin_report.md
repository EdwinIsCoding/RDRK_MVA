# Track 2 report: drug repurposing hypotheses for BUB1B-related mosaic variegated aneuploidy

**Submitter:** NexusDwin
**Date:** 3 September 2026
**Licence:** CC BY 4.0
**Repository:** see the GitHub URL supplied with this submission

---

## 1. Summary

The proband is a compound heterozygote for two loss-of-function alleles in
`BUB1B`, established in Track 1 and read-level verified. Loss of function is the
mechanism, so the therapeutic direction is restoration or compensation, never
inhibition.

**Our central result is a negative, and it is the most useful thing we found.**
The obvious therapeutic idea, restoring spindle assembly checkpoint function, is
not merely unproven. Nominating targets by signed causal edges from the mosaic
variegated aneuploidy seed genes yields ten targets, **every one of which
requires activation**. ChEMBL holds 118 drug-mechanism records across them and
**none acts in that direction**. Six of the same ten have inhibitors.

We then did something we expect few entries to do: we measured whether that
finding means what it appears to mean. Across all of ChEMBL only 359 human genes
have a drug with an activating mechanism, against 1,541 with an inhibiting one,
from 19,297 approved protein-coding genes. Activation is available for 1.86% of
the genome. At that rate, the expected number of activatable targets among ten is
0.19 and the probability of seeing zero is 0.83.

**So the finding is real but weaker than it first looks, and we report the
weaker version.** What survives is stated in section 4.

The tractable axis is cancer chemoprevention and surveillance, because the
presenting event is rhabdomyosarcoma. There we found that the chemoprevention
evidence base for this child's disease and this child's tumour is not thin but
**empty**: zero registered trials in mosaic variegated aneuploidy, zero naming
`BUB1B`, zero prevention trials of any kind in rhabdomyosarcoma. Every candidate
we can offer is transferred from a different syndrome with a different driver
gene, and that assumption does more work than any filter we apply afterwards.

We propose no treatment. Every output is a research hypothesis addressed to
researchers, with the experiment that would falsify it.

---

## 2. Mechanism characterisation

Required by the submission rules as the basis for the repurposing rationale.

### 2.1 The variants and their direction of effect

| | Allele 1 | Allele 2 |
|---|---|---|
| Position (GRCh38) | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
| HGVS | `NM_001211.6:c.2210T>G` | `NM_001211.6:c.3006T>G` |
| Protein | p.Leu737Ter | p.Asn1002Lys |
| Consequence | nonsense | missense |
| Effect | **loss of function** | **presumed hypomorphic, inferred not measured** |

**This is loss of function, not gain of function.** Allele 1 introduces a
premature termination codon at residue 737 of 1,050, truncating the protein
before and within the kinase domain. It is classified Pathogenic/Likely
pathogenic in ClinVar (VCV000533901.9) against mosaic variegated aneuploidy
syndrome 1, OMIM 257300.

Allele 2 lies within the kinase domain. **That it is hypomorphic rather than null
is an inference and we label it as one.** The reasoning is that complete BubR1
nullity is not compatible with life, so a living compound heterozygote carrying
one clear null must retain residual function from the other allele. That is an
argument from viability, not a measurement. No functional assay of p.Asn1002Lys
exists, and the same protein change reached through a different nucleotide
substitution is classified **uncertain significance** in ClinVar
(VCV004600147.1). A reader should treat the hypomorph model as the most probable
of several rather than as established.

### 2.2 The pathway disrupted

BubR1, the product of `BUB1B`, is a core component of the mitotic checkpoint
complex. The checkpoint delays anaphase until every kinetochore is correctly
attached to the spindle. Reduced BubR1 function shortens or weakens that delay.

### 2.3 The downstream biological consequence

Chromosomes missegregate. Daughter cells inherit whole-chromosome gains and
losses, cell to cell, which is the *variegated* aneuploidy the disease is named
for. Two consequences follow and they are the two things a therapy could
plausibly address:

1. **Cancer predisposition.** Aneuploidy is a tumour-initiating state, and
   rhabdomyosarcoma is the tumour most characteristically reported in
   `BUB1B`-related mosaic variegated aneuploidy. It is this proband's presenting
   event, coded HP:0002859.
2. **Chronic cellular stress.** An aneuploid cell expresses genes at
   stoichiometries its protein complexes were not built for, which imposes a
   proteostasis and energetic burden.

### 2.4 What follows for repurposing, and what does not

**The causal lesion is not druggable and we do not propose to drug it.** Every
cell in this child's body carries both alleles from conception. A small molecule
cannot restore a truncated protein, and it cannot undo aneuploidy that has
already occurred in a differentiated tissue. Plan section 7.1 says so and
sections 3 and 4 below turn that assertion into a measurement.

What is left is downstream: the tumour risk, and the stress state.

---

## 3. The direct axis, and why we closed it

### 3.1 Method

Targets were nominated by **signed causal edges** from OmniPath, which
aggregates SIGNOR, Reactome and others, rather than by network proximity.
Proximity in an unsigned graph tells you which protein sits near the disrupted
biology and never what to do to it, and a candidate proposed without a direction
is not a hypothesis but a gene name.

Where the sign was contradictory or absent, the target was **rejected rather
than assigned a plausible direction**. Four were rejected on that basis and are
named rather than dropped: `GSK3B`, `KAT2B`, `MAPK14`, `MAPK8`.

Each surviving target was then put to ChEMBL with its required direction
attached.

### 3.2 Result

| | |
|---|---:|
| candidate modulators with a signed edge into a seed gene | 14 |
| direction resolved | 10 |
| of those, requiring **activation** | **10** |
| requiring inhibition | 0 |
| ChEMBL drug-mechanism records across the ten | 118 |
| acting in the required direction | **0** |

The ten are `ATM`, `AURKA`, `AURKB`, `CDK1`, `CENPE`, `EGFR`, `KNL1`,
`MAD2L1BP`, `PLK1` and `TTK`. Every drug that exists for them is an inhibitor,
antagonist or blocker.

### 3.3 We then tested our own finding, and it came out weaker

A result like the one above invites a conclusion it does not support. If
activating drugs are rare everywhere, then "no activator across ten targets" is
a fact about drug discovery rather than about the spindle assembly checkpoint.

So we measured the denominator.

| | genes | share of protein-coding genome |
|---|---:|---:|
| have an activating drug mechanism in ChEMBL | 359 | 1.86% |
| have an inhibiting drug mechanism in ChEMBL | 1,541 | 7.99% |
| approved protein-coding genes, HGNC | 19,297 | |

Built from 1,269 activating and 4,897 inhibiting mechanism records.

**Inhibition is 4.3 times more available than activation.** At the genome-wide
activation rate, the expected number of activatable targets among ten is
**0.19**, and the probability of observing zero is **0.83**.

Observing zero is therefore unremarkable. We say so.

### 3.4 What survives, stated at the strength the evidence supports

1. **The axis is unavailable in fact.** Whatever the reason, no activating drug
   exists for any of these targets today. A repurposing proposal aimed at them
   has nothing to repurpose. Repurposing is constrained by what exists, so this
   is decisive for this track even though it is uninformative about biology.
2. **Every nominated target requires the scarce direction and none requires the
   plentiful one.** An axis containing a mix would have the better-supplied
   direction open to it. This one does not, and six of the ten have inhibitors
   sitting ready to be reached for by mistake.
3. **The safety argument is independent of availability.** Activating mitotic
   kinases such as `PLK1`, `AURKB` or `CDK1` in a child with a cancer
   predisposition syndrome would be contraindicated if an activator existed. The
   axis stays closed even if one is developed.

We expect several entries to propose `PLK1`, `AURKA` or `TTK` modulation for
this proband. Anyone doing so should state which direction they intend and check
whether a drug achieves it, because for these ten targets none does, and the
inhibitors that do exist point the wrong way for a loss-of-function disease.

---

## 4. Cancer chemoprevention and surveillance

The plan calls this the highest-value output of the track. The proband's
presenting event is a malignancy, which makes it the most directly relevant of
the remaining axes.

### 4.1 The evidence base for this child is empty

Asked first, because the answer governs what everything after it can claim.

| Question, put to ClinicalTrials.gov | Trials |
|---|---:|
| any interventional trial in mosaic variegated aneuploidy | **0** |
| any trial naming `BUB1B` | **0** |
| any prevention trial in rhabdomyosarcoma | **0** |
| any prevention trial with a drug, in rhabdomyosarcoma | **0** |

Not a thin evidence base. An empty one.

**This is the dominant limitation of the axis and no downstream filtering
repairs it.** Every candidate below is transferred from a different hereditary
cancer predisposition syndrome, with a different driver gene, a different tumour
spectrum and a different mechanism of tumour initiation. Familial adenomatous
polyposis is driven by biallelic `APC` loss and adenoma formation; this proband's
tumour arises from whole-chromosome missegregation. An agent that suppresses
adenoma formation has no established reason to suppress that.

### 4.2 Derivation

Candidates are derived from the registry's own MeSH heading **"Neoplastic
Syndromes, Hereditary"**, not from a list of syndromes we chose. A list we chose
would be a statement about what an author recalled. The heading makes it a
statement about what the registry holds, and anyone can rerun the query.

```
AREA[ConditionSearch]"Neoplastic Syndromes, Hereditary"
  AND AREA[DesignPrimaryPurpose]PREVENTION
  AND AREA[InterventionType]DRUG
```

| Stage | Count |
|---|---:|
| prevention trials with a drug intervention | 39 |
| distinct agent names after dose and formulation text is stripped | 49 |
| resolved to a ChEMBL molecule | 25 |
| unresolved, reported as gaps rather than guessed at | 24 |

A name resolves only on an **exact** match against a ChEMBL preferred name or
synonym. Fuzzy search returns a nearest molecule for almost any string, and a
plausible wrong ChEMBL identifier in a report is worse than a gap because a
reader cannot tell it is wrong.

### 4.3 The safety screen

Deterministic rules over structured drug fields. No model judgement enters a
verdict, and every verdict names the rule that produced it and the field that
rule read. This is the part of Track 2 that carries the most scientific content
for this patient, and it is built as code rather than prose.

| Verdict | Agents | Meaning |
|---|---:|---|
| allowed | 16 | passes the screen; still only a hypothesis |
| flagged | 5 | a real tension that must be stated wherever the agent is proposed |
| unknown | 4 | required paediatric evidence absent; **not a pass** |
| excluded | 0 | categorical; no efficacy argument overrides |

`UNKNOWN` being distinct from `ALLOWED` is deliberate. A compound with no
paediatric exposure data is not thereby safe for a child.

**A screen that excludes nothing may be permissive or may be broken, and the
candidate set cannot tell you which**, because prevention trials rarely test
cytotoxics. So known cytotoxic agents are pushed through the identical code
path, named by the WHO ATC L01A alkylating-agent subgroup rather than by us:
bendamustine, busulfan, carboquone, carmustine and chlorambucil. **Five of five
are excluded.** The empty exclusion list above is therefore a property of the
candidate set and not of a broken screen.

### 4.4 Candidates

Ordered by screen verdict. This is a prioritisation for research evaluation and
not a ranking of clinical preference. **No dose appears anywhere in this
repository and none may be added.**

| Agent | ChEMBL | Prevention trials | Endpoint measured | Verdict |
|---|---|---:|---|---|
| SULINDAC | CHEMBL15770 | 1 | duodenal polyp burden | allowed |
| URSODIOL | CHEMBL1551 | 1 | number and size of duodenal adenomas | allowed |
| MESALAMINE | CHEMBL704 | 2 | occurrence of colorectal neoplasia | allowed |
| LETROZOLE | CHEMBL1444 | 1 | invasive breast cancer at 5 years | allowed |
| NOGAPENDEKIN ALFA | CHEMBL4297690 | 1 | cumulative incidence of adenomas | allowed |
| ASPIRIN | CHEMBL25 | 2 | Ki-67 and apoptosis, a surrogate | allowed |
| NAPROXEN | CHEMBL154 | 2 | PGE2 concentration, a surrogate | allowed |
| ATORVASTATIN | CHEMBL1487 | 1 | Ki-67 and apoptosis, a surrogate | allowed |
| **CELECOXIB** | CHEMBL118 | **5** | number and size of duodenal adenomas | **flagged** |
| ERLOTINIB | CHEMBL553 | 2 | duodenal polyp burden | flagged |
| PATIDEGIB | CHEMBL538867 | 1 | number of new basal cell carcinomas | unknown |

Celecoxib carries the most registry support of any agent in the set, at five
prevention trials, with an endpoint that counts lesions rather than a biomarker.

**Fourteen of the twenty-five resolved agents have at least one primary outcome
that counts tumours or lesions.** The rest reached the list through a surrogate
endpoint such as Ki-67 staining, or through an endpoint unrelated to cancer such
as seizure frequency in tuberous sclerosis. That makes them weaker candidates
rather than stronger ones, and the full table in
`results/summaries/track2_chemoprevention.md` prints each trial's primary
outcome text beside our classification of it so a reader can overrule us.

### 4.5 Caveats that travel with the flagged agents

These are not footnotes. They are attached to the candidate wherever it appears.

- **CELECOXIB** carries ATC L01XX33 and is classified antineoplastic. It is not
  categorically cytotoxic, which is why it is not excluded, but the
  classification travels with it. Paediatric exposure evidence exists:
  NCT00934739, NCT02876094, NCT00474773, NCT02934191, NCT00006299.
- **ERLOTINIB** carries ATC L01EB02, likewise antineoplastic and likewise not
  cytotoxic. Paediatric exposure: NCT02233049, NCT00360854, NCT00570232,
  NCT00077454, NCT01032070.
- **SIROLIMUS** is both antineoplastic-classified (L01EG04) and immunosuppressive
  (L04AH01). The immunosuppression is a **direct tension with tumour immune
  surveillance in a child already predisposed to cancer**, and it is the reason
  we do not propose a rapalog for the proteostasis axis despite the mechanistic
  appeal described in section 5. Naming that tension is the point; a proposal
  that omits it is the error a clinically trained reader will find first.
- **PIRFENIDONE** is immunosuppressive, with the same tension.

### 4.6 Surveillance

We deliberately do not propose a surveillance protocol. Producing an imaging or
monitoring schedule for a named child would be a clinical recommendation, which
`CLAUDE.md` rule 3 and `ETHICS.md` section 4 forbid this project from generating,
and which is properly the work of the clinical genetics and oncology services
already involved.

What we can say is what the analysis supports: the tumour risk is the dominant
modifiable component of this phenotype, no chemoprevention agent has been tested
against it, and surveillance is therefore the only intervention in this axis with
an established evidence base in cancer predisposition generally. A research
programme here would be better spent on surveillance protocol evaluation than on
any candidate in section 4.4.

---

## 5. The remaining axes, measured the same way

Both require increasing the activity of something, so both face the same
availability question as the direct axis. Gene sets are defined by GO terms
resolved through QuickGO by exact name match, and counted twice: over all
evidence codes and over experimental codes only, because electronically inferred
annotations dominate large GO terms without adding evidence.

| Axis | GO terms | Genes | With an activating drug | Expected at base rate | Ratio |
|---|---|---:|---:|---:|---:|
| Proteotoxic stress mitigation | GO:0034620, GO:0016236, GO:0044183 | 702 | 11 | 13.1 | **0.8** |
| Mitochondrial and oxidative support | GO:0034599, GO:0033108, GO:0009060 | 824 | 23 | 15.3 | **1.5** |

**Proteotoxic stress mitigation is no better supplied than the genome average.**
Eleven available targets from 702 genes is what picking genes at random would
give. It differs from the direct axis in being non-empty rather than in being
enriched, and we do not blur that distinction. The mechanistic story is
attractive, the pharmacology is ordinary, and the most obvious agent for it,
sirolimus, is flagged for immunosuppression in exactly the patient who can least
afford it.

**Mitochondrial and oxidative support is the only axis of the three that is
better supplied than the genome average**, at 1.5 times expectation. That makes
it the axis where a repurposing search is most likely to find something to work
with. It does not make it the axis most likely to help this child, and we have no
evidence bearing on that question.

A gene appearing in these counts means ChEMBL records a drug acting on it in the
required direction. It does not mean the drug is safe in a child, reaches the
relevant tissue, or that the GO annotation reflects the biology that matters
here. GO annotation sets are not therapeutic targets: they are broad, they
overlap, and membership is a claim about a process rather than about a point of
intervention in it.

---

## 6. What would falsify each hypothesis

A candidate without a falsification route is not a scientific claim.

| Hypothesis | The experiment that would falsify it |
|---|---|
| Chemoprevention transfers from other hereditary cancer syndromes to `BUB1B`-related MVA | A tumour-incidence study in a `Bub1b` hypomorphic mouse, which is the standard model for this disease. If an agent that suppresses adenoma formation in `Apc` mutants does not reduce tumour incidence there, the transfer fails. |
| The second allele is hypomorphic rather than null | A functional assay of p.Asn1002Lys: kinase activity, mitotic checkpoint competence, or protein stability in a cell line. The whole compensation framing rests on this and it has never been measured. |
| The direct axis is closed | Development of a selective activator of any of the ten targets, which would reopen availability while leaving the safety objection standing. |
| Proteostasis burden is therapeutically relevant here | Proteotoxic stress markers in patient-derived cells against controls. We have no patient cells and no patient transcriptome, so this is unaddressed rather than answered. |

---

## 7. Limitations

Stated plainly, and none of them is repaired elsewhere in this report.

1. **No chemoprevention evidence exists for this disease or this tumour.** Every
   candidate is a transfer across syndromes, and the transfer is the weakest
   link in the argument.
2. **No patient transcriptome.** There is no RNA-seq for this proband, so
   signature-reversal methods such as LINCS or CMap connectivity, which plan
   section 7.3 identifies as the highest-yield repurposing approach available,
   cannot be run against patient data at all. We did not substitute a published
   proxy signature and present it as though it were the patient's.
3. **No patient cells and no functional assay.** Every mechanistic statement
   about allele 2 is an inference from viability.
4. **The hypomorph model is unmeasured**, and the same protein change is a
   variant of uncertain significance in ClinVar.
5. **GO-derived gene sets are coarse** and their overlap with druggable space is
   a statement about annotation as much as about biology.
6. **The safety screen reads structured fields, not full drug labels.** It can
   only act on what the fields record. Absence of a flag is not evidence of
   safety, which is why `UNKNOWN` exists as a separate verdict.
7. **No dosing, dose range, safety margin or clinical recommendation appears
   anywhere in this work**, by rule and by design. Nothing here is addressed to
   a clinician treating this child.

---

## 8. Reproducibility

Every number in this report regenerates from the repository.

```bash
make track2        # direction audit, chemoprevention axis, axis availability
make test          # the full automated suite, which prints its own count
```

Outputs land in `results/summaries/`: `track2_direction_audit.md`,
`track2_chemoprevention.md` and `track2_axis_availability.md`. Sources are
OmniPath signed causal edges, ChEMBL, ClinicalTrials.gov API v2, QuickGO and
HGNC. Every identifier in this report was returned by one of those, and where a
lookup failed the name is carried as a gap rather than filled with a plausible
guess.

The safety screen and the derivation logic are covered by tests, including tests
that assert no dose reaches an output, that a near-miss ChEMBL search returns
nothing rather than the nearest molecule, and that the base-rate arithmetic which
weakens our own headline claim is correct.

---

## 9. Data handling and AI assistant disclosure

Required by the organisers' update of 28 August 2026.

> **Anthropic Claude (Claude Code), Max subscription, consumer terms. The "help
> improve Claude" training setting was enabled until 31 August 2026 and disabled
> thereafter. No patient genomic data (VCF, BAM, variant tables or pasted variant
> blocks) entered the model context at any point. Clinical phenotype terms and
> aggregate callset statistics did.**

We disclose a compliance gap rather than presenting a clean summary. For part of
this work the account permitted training on conversation content, which does not
meet the organisers' condition of no training on inputs or outputs. The setting
has been disabled. Notification to Sage Bionetworks' Privacy and Compliance
Office was outstanding at the time of writing. The full account is in
`ETHICS.md` section 3a and section 8 of the Track 1 report, and it is not
softened here.

No patient data of any kind was used in the Track 2 analysis above. It runs on
public databases and on the two variant coordinates established in Track 1.

---

## 10. What this report is

Research hypotheses addressed to researchers, with their evidence, their gaps
and the experiments that would falsify them. Nothing here is a diagnosis, a
treatment recommendation or clinical advice. A polished dossier confers authority
that this evidence does not support, and the family may read this.

The most useful sentence we can offer is not a drug name. It is that the direct
therapeutic axis is unavailable, that the chemoprevention evidence base for this
disease is empty, and that a research programme aimed at this child's tumour risk
would currently be better spent on surveillance than on any repurposing candidate
we can name.
