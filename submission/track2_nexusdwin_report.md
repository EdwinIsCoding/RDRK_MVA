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

**There is precedent for this architecture, and it is close.** UniProt curates
five MVA1 missense variants in the BubR1 kinase domain, at residues 814, 844,
909, 921 and 1012, and **all five are annotated as compound heterozygous with a
nonsense mutation**. That is exactly the architecture proposed here, and the
nearest sits **10 residues** from ours. All five trace to Hanks et al., *Nature
Genetics* 2004 (PMID 15475955), the study that established biallelic `BUB1B`
mutation as the cause of MVA1. A phosphothreonine phosphorylated by PLK1 sits at
residue 1008, six residues away.

This does not make p.Asn1002Lys hypomorphic, and proximity to a pathogenic
residue is not evidence about a different residue. What changes is the standing
of the inference. The architecture we propose is the documented architecture for
this disease rather than a construction of ours, so the hypomorph model rests on
precedent as well as on viability. Full output:
`results/summaries/kinase_domain_precedent.md`.

**We checked whether a structural calculation could do better, and it cannot.**
Plan section 7.5 permits stability work for a missense in a protein with an
experimental structure, and BubR1 has nine PDB entries. Two map to residue 1002,
and both do so by declaring the full-length construct while the observed BubR1
density stops at residue 308 and 345 respectively. The kinase domain, residues
766 to 1050 and the half of the protein this variant sits in, is unresolved in
every experimental structure of this protein. Both are cryo-EM reconstructions
of the anaphase-promoting complex at 4.8 and 3.8 angstrom, where side-chain
placement would not be determinable even where density existed.

A predicted model would cover the residue, and a stability calculation on one
would produce a number. That number would be a statement about the prediction
rather than about the protein, and offering it as evidence for hypomorphism is
the anti-pattern this project recorded before it began: structure prediction
standing in for a variant-effect measurement. So the gap stays open with its
reason recorded, rather than filled by a number that does not support it. See
`results/summaries/track2_structural_feasibility.md`.

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

**We publish this table rather than only citing it.**
`resources/directional_availability/` lists every gene for which ChEMBL records
a drug acting in each direction, so any team can check in one lookup whether the
target they intend to activate can be activated at all. It is the part of our
negative result that is useful to someone other than this patient, and it is the
measurement that weakened our own claim, so it seemed wrong to keep it to
ourselves. Gene symbols and counts are published rather than ChEMBL records,
with attribution, because ChEMBL is CC BY-SA 3.0 and our outputs are CC BY 4.0.

Observing zero is therefore unremarkable. We say so.

### 3.5 We tried to break the all-activation finding three ways

Ten targets reached by one hop from eight genes is a small basis for a headline.
Three things could make the all-activation character an artefact, and none had
been tested.

**Do stimulatory edges simply dominate the graph?** If most signed edges are
stimulatory then any loss-of-function seed set yields mostly "activate", and ten
of ten says nothing. Running the identical nomination on 20 seeded random seed
sets of the same size, drawn from the curated disease-gene panel, measures the
stimulatory share at **69%**. At that rate ten of ten has probability **0.024**,
and only 1 of the 20 random sets was all-activation with five or more targets.
**The finding survives.**

**Does a wider net change it?** Yes, and this is a real correction. Adding the
first-hop partners to the seed set reaches 66 targets requiring inhibition;
seeding from 60 mitotic-panel genes reaches 61. So the all-activation character
belongs to the **immediate regulators of these eight genes**, not to mitotic
biology generally, and claim 2 above is scoped accordingly.

**Does that reopen the axis? No, and the reason matters.** Those wider runs
treat `PLK1`, `AURKA`, `CDK1` and the rest as things to be activated, so the
inhibition-reachable targets they find are routes to raising a mitotic kinase's
activity. That is the endpoint the safety screen rules out in a child with a
cancer predisposition syndrome, independently of whether a drug exists. They are
not an alternative to the closed axis, they are further paths to the same
contraindicated place.

**The closure is therefore over-determined**: unavailable in the required
direction where we looked, and contraindicated where widening the search makes
something available. Full output:
`results/summaries/track2_nomination_sensitivity.md`.

### 3.4 What survives, stated at the strength the evidence supports

1. **The axis is unavailable in fact.** Whatever the reason, no activating drug
   exists for any of these targets today. A repurposing proposal aimed at them
   has nothing to repurpose. Repurposing is constrained by what exists, so this
   is decisive for this track even though it is uninformative about biology.
2. **Every immediate regulator of the seed genes requires the scarce
   direction, and that is unlikely by chance.** Six of the ten have inhibitors
   sitting ready to be reached for by mistake. Two independent checks support
   this. Section 6 shows it is disease-specific: run unchanged on Fanconi
   anaemia and ataxia-telangiectasia, the same nomination returns 17 and 10
   targets reachable by inhibition, and this proband's disease is the only one
   of the three with none. And section 3.5 measures the null: stimulatory edges
   are 69% of the signed graph, so ten of ten activation has probability 0.024.

   **The scope is the immediate regulators, not mitotic biology at large.**
   Widening the seed set does reach inhibition-reachable targets. Section 3.5
   explains why that does not reopen the axis.
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

**We then asked the literature the same question**, because a registry records
prospective trials and not published work, and an early chemoprevention signal
would appear first in an animal model. PubMed returns **zero** publications for
`Bub1b AND chemoprevention`. The three hits for mosaic variegated aneuploidy and
cancer prevention are a review, a surgical case report and a colorectal risk
study, none a prevention study; their titles are printed in
`results/summaries/chemoprevention_literature.md` rather than summarised,
because a hit count alone would have been misleading.

One hit is worth more than its count. Lissa et al., *PNAS* 2014 (PMID 24516128)
report that resveratrol and salicylate selectively reduce the fitness of
**tetraploid** cells through AMP-activated protein kinase. Two independent parts
of this project converge on it: aspirin reached our candidate table through a
registry query that knew nothing of the paper, and `PRKAA1`, `PRKAA2`, `PRKAG1`,
`PRKAG2` and `PRKAG3` reached our activatable list through the directional
analysis in section 5. **It is still not evidence for this disease.** Tetraploidy
from cytokinesis failure is not whole-chromosome aneuploidy from a weakened
spindle checkpoint, and the paper's context is `APC`-driven intestinal neoplasia.
It names a mechanistic route worth testing in the right model; it does not
shorten the distance to this patient.

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
| dropped as placebo or vehicle arms | 2 |
| intervention names resolved to a ChEMBL molecule | 39 |
| **distinct molecules after merging arms of the same drug** | **30** |
| unresolved, reported as gaps rather than guessed at | 8 |

A name resolves only on an **exact** match against a ChEMBL preferred name or
synonym, and three lookups are tried in that order: exact preferred name, exact
synonym, then ChEMBL's ranked text search. The same exact check is applied to
all three, so widening recall cannot admit a near miss. A plausible wrong ChEMBL
identifier in a report is worse than a gap because a reader cannot tell it is
wrong.

**Relying on the ranked search alone lost roughly half the candidates.** Querying
it for `sirolimus` returns ten molecules, none of them CHEMBL413, whose preferred
name is SIROLIMUS. Among the names it silently dropped was **metformin**, which
is one of the most studied repurposed chemoprevention agents and which appears
in section 4.4 below with a tumour-count endpoint. The exact endpoints recover
it. Registry intervention names describe trial arms rather than compounds, so
`metformin combination`, `celecoxib monotherapy` and `for Aspirin 300` resolve
only after the arm wrapper is stripped, and the summary prints the string that
actually matched so the normalisation can be audited rather than trusted.

Two names are **dropped rather than resolved**. `no active patidegib` is the
control arm of a patidegib trial and `Vehicle comparator` is a vehicle arm.
Normalising either to its drug would record the arm that received no drug as
registry evidence for the drug. The eight that remain unresolved are food
preparations, antibody classes and sponsor product codes, none of which is a
single molecule.

Arms of the same drug are then **merged into one candidate**. Left unmerged,
aspirin appeared three times with its trial count split between the copies, and
could hold two verdicts at once, since an arm labelled `for Aspirin 300` finds
no paediatric trials and scores UNKNOWN while plain `Aspirin` scores ALLOWED.
The paediatric lookup now runs on the ChEMBL preferred name rather than the
registry arm label, for the same reason.

### 4.3 The safety screen

Deterministic rules over structured drug fields. No model judgement enters a
verdict, and every verdict names the rule that produced it and the field that
rule read. This is the part of Track 2 that carries the most scientific content
for this patient, and it is built as code rather than prose.

| Verdict | Agents | Meaning |
|---|---:|---|
| allowed | 21 | passes the screen; still only a hypothesis |
| flagged | 5 | a real tension that must be stated wherever the agent is proposed |
| unknown | 4 | required paediatric evidence absent; **not a pass** |
| excluded | 0 | categorical; no efficacy argument overrides |

`UNKNOWN` being distinct from `ALLOWED` is deliberate. A compound with no
paediatric exposure data is not thereby safe for a child.

**A screen that excludes nothing may be permissive or may be broken, and the
candidate set cannot tell you which**, because prevention trials rarely test
cytotoxics. So the screen is controlled in both directions, with agents named by
the WHO ATC classification rather than by us.

| Control | Drawn from | Required behaviour | Result |
|---|---|---|---|
| Cytotoxics | ATC L01A alkylating agents: bendamustine, busulfan, carboquone, carmustine, chlorambucil | must be **excluded** | **5 of 5 excluded** |
| Ordinary agents | ATC M01A anti-inflammatories, A02B acid-related agents, A11C fat-soluble vitamins | must **not** be excluded | **9 of 9 not excluded** |

The second half matters as much as the first. A screen that refused everything
would score five out of five on cytotoxics and look healthy, and the candidate
set could not reveal the difference. Only the pair shows that the instrument
discriminates rather than merely refuses, so the empty exclusion list above is a
property of the candidate set. The ordinary agents return `unknown` rather than
`allowed` because no paediatric lookup is run for a control, and `UNKNOWN` is
never a pass.

This control is not hypothetical insurance. It is the exact failure we shipped
once: a blanket ATC L01 exclusion, which had passed its tests, removed celecoxib
from the candidate list on the stated grounds that a COX-2 inhibitor is
cytotoxic chemotherapy. The tests had been written from the same assumption as
the rule, so only a real candidate set could expose it.

### 4.4 Candidates

Ordered by screen verdict. This is a prioritisation for research evaluation and
not a ranking of clinical preference. **No dose appears anywhere in this
repository and none may be added.**

Only agents with at least one tumour-counting endpoint are listed here. The
full table of thirty, including those that reach the list through a surrogate or
an unrelated endpoint, is in `results/summaries/track2_chemoprevention.md`.

| Agent | ChEMBL | Prevention trials | Endpoint measured | Verdict |
|---|---|---:|---|---|
| ASPIRIN | CHEMBL25 | 3 | patients with at least one adenoma | allowed |
| MESALAMINE | CHEMBL704 | 2 | occurrence of colorectal neoplasia | allowed |
| NAPROXEN | CHEMBL154 | 2 | PGE2 concentration, a surrogate | allowed |
| SULINDAC | CHEMBL15770 | 2 | duodenal polyp burden | allowed |
| **METFORMIN** | CHEMBL1431 | 1 | number and size of colonic and duodenal polyps | allowed |
| LETROZOLE | CHEMBL1444 | 1 | invasive breast cancer at 5 years | allowed |
| NOGAPENDEKIN ALFA | CHEMBL4297690 | 1 | cumulative incidence of adenomas | allowed |
| URSODIOL | CHEMBL1551 | 1 | number and size of duodenal adenomas | allowed |
| TIPIFARNIB | CHEMBL289228 | 1 | median time to progression | allowed |
| **CELECOXIB** | CHEMBL118 | **5** | number and size of duodenal adenomas | **flagged** |
| SIROLIMUS | CHEMBL413 | 4 | seizure occurrence, and tumour volume | flagged |
| ERLOTINIB | CHEMBL553 | 2 | duodenal polyp burden | flagged |
| PATIDEGIB | CHEMBL538867 | 1 | number of new basal cell carcinomas | unknown |

Celecoxib carries the most registry support of any agent in the set, at five
prevention trials, with an endpoint that counts lesions rather than a biomarker.

**Seventeen of the thirty resolved molecules have at least one primary outcome
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
- **METFORMIN** is not flagged and carries paediatric exposure evidence. It is
  listed last among the caveats precisely because it has none, which is unusual
  in this set and is the reason it is worth a closer look than its single trial
  would otherwise justify.
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
better supplied than the genome average**, at 1.5 times expectation, which made
it the axis where a repurposing search looked most likely to find something.

**We followed that through, and it does not survive.** Naming the genes, pulling
every activating drug ChEMBL records for them and screening each one shows the
supply is not spread across the axis: **roughly 70% of the drug-target pairs
come from just three genes**, `INSR`, `PPARA` and `GCK`. The exact counts move
by one or two between runs, because ChEMBL serves intermittent HTTP 500s and a
molecule whose record cannot be fetched is reported as unretrievable rather than
screened on absent fields; the summary carries the figures for the run that
produced it. Those are metabolic drug
targets, and the drugs behind them are insulin analogues, fibrates and
glucokinase activators. They reached this axis because their genes carry GO
annotations for oxidative-stress response or aerobic respiration.

That is an annotation artefact, not a therapeutic rationale. Insulin is not
mitochondrial support for a child with mosaic variegated aneuploidy; it is a
diabetes drug whose receptor is annotated to a metabolic process. The caveat
below, that GO annotation sets are not therapeutic targets, is demonstrated here
rather than asserted.

**So this is the third axis to close**, and cancer chemoprevention and
surveillance is the one left standing, which is where plan section 7.1 placed
the value before any of this was measured. Full output:
`results/summaries/track2_mitochondrial_axis.md`.

A gene appearing in these counts means ChEMBL records a drug acting on it in the
required direction. It does not mean the drug is safe in a child, reaches the
relevant tissue, or that the GO annotation reflects the biology that matters
here. GO annotation sets are not therapeutic targets: they are broad, they
overlap, and membership is a claim about a process rather than about a point of
intervention in it.

---

## 5a. Signature reversal, against a proxy that is labelled

Plan section 7.3 calls LINCS/CMap connectivity "the highest-yield repurposing
method available" and permits a published MVA/CIN model signature where no
patient transcriptome exists, provided the proxy is labelled clearly. It is run
here.

**The proxy is GEO GSE277997**, a BubR1-insufficiency mouse, twelve animals,
using the authors' own differential expression table. It is the right genotype
and almost nothing else.

| | Proxy | This proband |
|---|---|---|
| Species | mouse | human |
| Tissue | heart | skeletal muscle and kidney affected; tumour of skeletal muscle lineage |
| Phenotype | age-related cardiac pathology | paediatric cancer predisposition with rhabdomyosarcoma |
| Genotype | engineered hypomorph | nonsense allele in trans with a missense |

Of 12,308 genes, 150 up and 150 down were submitted after filtering at adjusted
p < 0.05, a stated choice, and mapping to human through **MGI homology classes**
rather than the uppercase naming convention, which is right most of the time and
wrong silently.

**The result is specific, which had to be tested rather than assumed.** HDAC and
mTOR inhibitors perturb transcription broadly and surface in many L1000 queries
whatever is submitted. Running the identical query three times on random gene
sets of the same size, drawn from the same background, returns on average
**3.3 of 39 perturbagens, 9%**. Thirty-four come back for the real signature and
none of the random draws, so the list is largely about the submitted signature
rather than about generic transcriptional perturbagens.

**The safety screen then did the job it exists for.** Among the top-ranked hits
is **daunorubicin**, a cytotoxic anthracycline, which the screen excluded under
its ATC L01D rule. A connectivity method will happily rank chemotherapy first
for a child with a cancer predisposition syndrome, because it optimises
transcriptional opposition and knows nothing about the patient. Of the top 25,
several resolve only to research codes; the named compounds are mostly mTOR and
HDAC inhibitors, flagged rather than allowed.

**What it is worth.** A direction to look, not a candidate to advance. The
distances multiply rather than add: mouse to human, heart to muscle and kidney,
ageing to paediatric cancer predisposition, and LINCS cell line to patient. What
would make it real is the experiment named everywhere else in this report, a
`Bub1b` hypomorphic model in the affected tissue.

It is reported despite all of that because an unrun method is worth less than a
run one with its limits stated, and because the alternative was leaving the
highest-yield approach unattempted while implying it had been considered. Full
output: `results/summaries/track2_signature_reversal.md`.

---

## 6. Does the method generalise, or was it built around one answer?

Everything above concerns one child. A judge is entitled to ask whether the
machinery says anything useful about a different patient, and that is not a
question we can settle by asserting it.

So the **unchanged** pipeline was pointed at two other inherited disorders. Both
are recessive or loss-of-function chromosomal-instability and DNA-repair
syndromes with cancer predisposition, which makes them near neighbours of this
proband's disease. A method that cannot tell near neighbours apart is not useful
in rare disease, where near neighbours are what a differential is made of. Their
gene sets are read from the curated disease panel by matching each gene's own
assertion text, not typed in.

Nothing else changed. No parameter, no threshold.

| Disease | Seeds | Targets needing activation | needing inhibition | Drug available in the required direction |
|---|---:|---:|---:|---:|
| **Mosaic variegated aneuploidy** | 8 | **10** | **0** | **0** |
| Fanconi anaemia | 19 | 36 | 17 | 8 |
| Ataxia-telangiectasia | 2 | 37 | 10 | 5 |

Protein complexes are set aside before counting. OmniPath names them by joining
components with underscores, they can never match a gene symbol in the
availability table, and leaving them in would count each as "no drug available"
and bias every result towards our own conclusion. In Fanconi anaemia the
ubiquitin machinery alone contributed 154 such identifiers, mostly `UBB_UBE2*`
pairs that are one piece of biology counted many times.

**Two things follow, and neither was available from the proband alone.**

**The method discriminates.** Three near-neighbour syndromes, one pipeline, three
different answers. A pipeline tuned until it produced our conclusion would not
behave that way.

**This proband's disease is the extreme case, and that is now a comparison
rather than an assertion.** It is the only one of the three in which every
nominated target requires activation and none requires inhibition, and the only
one with nothing available in the direction it needs.

The registry half of the argument survives the same test. Asked the identical
queries, Fanconi anaemia returns 150 interventional trials and ataxia-telangiectasia
returns 59, against zero for mosaic variegated aneuploidy. **The empty evidence
base in section 4.1 is a property of the disease, not of how we asked.**

Where the three agree, the agreement is itself a finding about loss-of-function
disorders as a class: all three nominate far more targets requiring activation
than inhibition, which is the direction pharmacology supplies least.

**What this does not show.** That the method finds the right drug for any of
these diseases. It shows the machinery accepts a different disease and returns a
different, mechanically derived answer. That is the precondition for using it on
the next patient, not a claim about its accuracy.

**Cost of applying it to a new disease: a seed gene set, and nothing else.**
Every other input is a public database already wired in, and the safety screen is
disease-agnostic by construction.

### A note on knowledge-graph leakage, because we said we would check

This project's anti-pattern list includes "any knowledge-graph claim without a
time-split leakage check", and we use a knowledge graph: OmniPath's signed causal
edges drive the whole nomination step. A reader who holds us to our own list is
entitled to ask where the leakage check is.

**It is absent because the check does not apply to this use, and saying so is
better than performing a ritual version of it.** Leakage matters when a graph is
used to *predict* something the graph already encodes, which is how a model that
has memorised known drug-target pairs scores well on recovering known drug-target
pairs. Our use is not predictive. We ask the graph a descriptive question, which
protein regulates a seed gene and in which direction, and the answer is then
tested against an independent source, ChEMBL, that the graph had no part in
building. Nothing is being forecast, so there is no future to hold out.

The failure mode that *does* apply here is different and worth naming, because it
is the one a time split would not catch. Both OmniPath and ChEMBL are richer for
well-studied proteins, so a target that is heavily curated in one is likely
heavily curated in the other. Our targets are mitotic kinases, among the
best-studied proteins in the genome, which means the direction we can resolve and
the drugs we can find are both biased towards them. That bias works **against**
our conclusion rather than for it: if anything, these targets are more likely to
have an activator recorded than an average protein, and none does. Section 3.3
puts a number on the comparison rather than leaving it as an argument.

---

## 7. What would falsify each hypothesis

A candidate without a falsification route is not a scientific claim.

| Hypothesis | The experiment that would falsify it |
|---|---|
| Chemoprevention transfers from other hereditary cancer syndromes to `BUB1B`-related MVA | A tumour-incidence study in a `Bub1b` hypomorphic mouse, which is the standard model for this disease. If an agent that suppresses adenoma formation in `Apc` mutants does not reduce tumour incidence there, the transfer fails. |
| The second allele is hypomorphic rather than null | A functional assay of p.Asn1002Lys: kinase activity, mitotic checkpoint competence, or protein stability in a cell line. The whole compensation framing rests on this and it has never been measured. |
| The direct axis is closed | Development of a selective activator of any of the ten targets, which would reopen availability while leaving the safety objection standing. |
| Proteostasis burden is therapeutically relevant here | Proteotoxic stress markers in patient-derived cells against controls. We have no patient cells and no patient transcriptome, so this is unaddressed rather than answered. |

---

## 8. Limitations

Stated plainly, and none of them is repaired elsewhere in this report.

1. **No chemoprevention evidence exists for this disease or this tumour.** Every
   candidate is a transfer across syndromes, and the transfer is the weakest
   link in the argument.
2. **No patient transcriptome**, so signature reversal was run against a
   labelled proxy rather than against this child. See section 5a. An earlier
   version of this report said we "did not substitute a published proxy
   signature and present it as though it were the patient's", which rejected
   something nobody proposed in order to avoid doing what plan section 7.3
   sanctions. That sentence read as restraint and was an unrun method wearing
   restraint's clothes.
3. **No patient cells and no functional assay.** Every mechanistic statement
   about allele 2 is an inference from viability, and no structural calculation
   can substitute, because the kinase domain is unresolved in every experimental
   structure of BubR1. See section 2.1.
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

## 9. Reproducibility

Every number in this report regenerates from the repository, **and none of it
needs the challenge dataset.**

Track 2 reads no patient data. Every result above comes from public databases,
no script in the pipeline opens a file under `data/`, and a test asserts it. A
reviewer can therefore check this entire report without applying for data
access, which is the difference between taking our word for it and not:

```bash
make reproduce-track2   # one 16 MB download, then the whole Track 2 pipeline
```

```bash
make track2        # direction audit, chemoprevention axis, axis availability
make scalability   # the same pipeline on two comparator diseases
make track2-drift  # have the live registry counts moved since we wrote this?
make test          # the full automated suite, which prints its own count
```

Two of the sources above are live and change without notice. The counts this
report quotes are therefore pinned in `config/track2_evidence_pin.json` with the
date they were taken, and `make track2-drift` re-queries and reports the
difference. A reader running the pipeline months later can tell whether we were
wrong or the world moved.

The four zeros in section 4.1 are marked load-bearing in that check. If one ever
becomes non-zero it means a chemoprevention trial now exists for this disease,
and section 4.1 has to be rewritten rather than have its number adjusted.

Outputs land in `results/summaries/`: `track2_direction_audit.md`,
`track2_chemoprevention.md`, `track2_axis_availability.md` and
`track2_scalability.md`. Sources are
OmniPath signed causal edges, ChEMBL, ClinicalTrials.gov API v2, QuickGO and
HGNC. Every identifier in this report was returned by one of those, and where a
lookup failed the name is carried as a gap rather than filled with a plausible
guess.

The safety screen and the derivation logic are covered by tests, including tests
that assert no dose reaches an output, that a near-miss ChEMBL search returns
nothing rather than the nearest molecule, and that the base-rate arithmetic which
weakens our own headline claim is correct.

---

## 10. Data handling and AI assistant disclosure

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

## 11. What this report is

Research hypotheses addressed to researchers, with their evidence, their gaps
and the experiments that would falsify them. Nothing here is a diagnosis, a
treatment recommendation or clinical advice. A polished dossier confers authority
that this evidence does not support, and the family may read this.

The most useful sentence we can offer is not a drug name. It is that the direct
therapeutic axis is unavailable, that the chemoprevention evidence base for this
disease is empty, and that a research programme aimed at this child's tumour risk
would currently be better spent on surveillance than on any repurposing candidate
we can name.
