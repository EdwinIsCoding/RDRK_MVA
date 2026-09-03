# Verification of the Track 1 work

An adversarial re-check of the Track 1 claims against primary sources rather
than against the previous agent's prose. Performed 3 September 2026 by the
agent that took over the project.

Method: every claim below was tested against the thing it is a claim about.
Coordinates and genotypes against the original callset in `data/`, not against
the VEP output. Reference bases against the GRCh38 FASTA. ClinVar accessions
against NCBI E-utilities. Population frequencies against gnomAD v4.1. Read-level
figures against the panel BAM. Derived numbers by re-running the scripts that
produced them.

Two scripts were added and their outputs are the evidence for what follows:

- `scripts/25_verify_track1_claims.py` writes `results/summaries/verification_callset.md`
- `scripts/26_verify_readlevel.py` writes `results/summaries/verification_readlevel.md`

Both obey `CLAUDE.md` rule 1. They read `data/` and the BAM; their outputs are
aggregate counts plus the five loci already published in
`submission/track1_submission.csv` under CC BY 4.0.

---

## Summary

**The Track 1 call is correct and survives every check that bears on it.** Both
alleles are in the callset exactly as claimed, the reference bases are right,
the ClinVar accession is right and is listed against mosaic variegated
aneuploidy syndrome 1, and the read-level evidence supports two germline
heterozygous variants in uniquely mappable sequence. Nothing found here changes
the answer or the ranking.

**Nine claims around it are wrong.** One is material to the science, three are
material to the write-up, and five are inaccuracies of counting or wording. The
most serious is that the second allele is described as absent from gnomAD when
the project's own lookup table records a frequency for it. The most consequential
for the report is that the two variants set aside as unresolved are in fact
common polymorphisms and are now closed.

---

## 1. Claims confirmed

| Claim | Source checked against | Result |
|---|---|---|
| `chr15:40209701 T>G`, GT 0/1, DP 46, AD 21/25, GQ 99, PASS | original callset | **Exact** |
| `chr15:40220612 T>G`, GT 0/1, DP 28, AD 15/13, GQ 99, PASS | original callset | **Exact** |
| `chr22:20996720 C>G`, GT 0/1, DP 48, AD 24/24, GQ 99, PASS | original callset | **Exact** |
| REF base T at both BUB1B positions | GRCh38 FASTA, Ensembl 115 | **Confirmed** |
| REF base C at the LZTR1 position | ClinVar SPDI `NC_000022.11:20996719:C:G` | **Confirmed** |
| ClinVar VCV000533901.9 Pathogenic/Likely pathogenic | NCBI E-utilities | **Confirmed**, criteria provided, multiple submitters, no conflicts, last evaluated 9 October 2024 |
| That record is listed against MVA1 | NCBI E-utilities | **Confirmed**, trait *Mosaic variegated aneuploidy syndrome 1*, OMIM 257300, MONDO:0009759 |
| ClinVar VCV001409252.7 Pathogenic/Likely pathogenic | NCBI E-utilities | **Confirmed**, multiple submitters, no conflicts, last evaluated 24 September 2025 |
| Its three conditions | NCBI E-utilities | **Confirmed**: Noonan syndrome 10 (OMIM 616564), Noonan syndrome 2 (OMIM 605275), LZTR1-related schwannomatosis (OMIM 615670) |
| Alleles 10,911 bp apart | arithmetic | **Confirmed**, 40220612 − 40209701 = 10,911 |
| 5,012,204 callset records | streamed count | **Confirmed** |
| 4,950,283 VEP-annotated records | streamed count | **Confirmed** |
| The 61,921 difference is entirely decoy and unplaced contigs | per-contig counts of both files | **Confirmed**, and stronger than claimed: losses on primary contigs (1-22, X, Y, MT) are exactly **zero** |
| 415 scored variants reduce to 12 | re-ran `scripts/18_arm_a_shortlist.py` | **Confirmed**, and the variant-level shortlist regenerated bit-identical |
| Both causal alleles were in that shortlist | membership test on the regenerated shortlist | **Confirmed**, 2 of 2, BUB1B carries 3 of the 12 |
| BUB1B ranked 18th of 2,503 genes on the HPO prior | re-ran `scripts/24_phenotype_prior.py` | **Confirmed**, 18th, score 9.64, and the corpus is 2,503 genes |
| Panel coverage 42-51x against a genome mean of 43.8x | mosdepth output | **Confirmed**, BUB1B 48.13x, genome mean 43.79x |
| Splice-distance concordance 268/268 SNVs | `tests/test_evidence.py`, which asserts concordance of exactly 1.0 | **Confirmed**, the test passes in a suite of 140 |
| The Track 2 direction audit | re-ran `make track2` | **Reproduces byte-identical**: 10 targets, all requiring activation, 118 ChEMBL mechanisms, none activating |

### A claim strengthened rather than merely confirmed

The report says phase is *inferred, not demonstrated* because the alleles lie
beyond a read pair. That was an assertion. It is now a measurement.

Extracting every read name at each allele from the panel BAM and intersecting
the two sets gives **zero reads or read pairs in common**. The observed
template length at allele 1 has a median of 502 bp and a maximum of 1,103 bp,
against a required span of 10,911 bp. No short-read data of this library can
resolve the phase, and none of ours does.

---

## 2. Discrepancies found

### 2.1 The second allele is not absent from gnomAD

**Claimed**, in the report, the submission notes file, `config/config.yaml` and
`README.md`: `chr15:40220612 T>G` is "novel", "absent from gnomAD", gnomAD
popmax "absent".

**Found**: it is present in gnomAD v4.1 exomes.

| | value |
|---|---|
| allele frequency | 6.84 × 10⁻⁷ |
| group max frequency | 8.99 × 10⁻⁷ |
| allele number | 1,461,878 |
| homozygotes | 0 |

That is a single observed allele in roughly 1.46 million. It is absent from
gnomAD genomes, which is presumably where the wording came from.

This is not an outside source contradicting the pipeline. **The project's own
lookup table records it**: `refs/gnomad_mva_known/panel_af.tsv.gz` holds the key
`15:40220612:T:G` with a group max frequency of 8.99 × 10⁻⁷, and the Arm A
summary counts it among the 4 "rare" variants rather than the 7 "absent" ones.
The pipeline was right and the write-up did not follow it.

**Consequence.** The variant remains ultra-rare and the argument is untouched.
The ACMG criterion should be PM2_supporting for an ultra-rare allele, not PM2
for an absent one. The word "novel" should go.

### 2.2 The p.Asn1002Lys protein change is in ClinVar

**Claimed**: allele 2 has "no ClinVar record".

**Found**: true of the nucleotide change, misleading about the protein change.
Two ClinVar records exist at `chr15:40220612`:

| Accession | Change | Protein | Classification |
|---|---|---|---|
| VCV004600147.1 | `c.3006T>A` | **p.Asn1002Lys**, the same substitution | Uncertain significance, single submitter |
| VCV004842697.1 | `c.3006T>C` | p.Asn1002= , synonymous | Likely benign, single submitter |

Both AAG and AAA encode lysine, so `c.3006T>G` and `c.3006T>A` produce an
identical protein change. Someone has seen p.Asn1002Lys before and classified it
**uncertain**.

This belongs in the report. It does not support PS1, which requires the same
amino acid change with an established pathogenic classification, and a variant of
uncertain significance is not that. It is honest context that cuts slightly
against the missense, and omitting it while claiming "no ClinVar record" reads as
stronger than the evidence.

### 2.3 The reported "popmax" figures are global allele frequencies

**Claimed**: allele 1 gnomAD v4.1 popmax 7.87 × 10⁻⁵; LZTR1 popmax 1.4 × 10⁻⁶.

**Found**: both numbers are the global allele frequency. The group max
frequencies are higher.

| Variant | Reported as popmax | Actual global AF | Actual group max |
|---|---|---|---|
| `chr15:40209701 T>G` | 7.87 × 10⁻⁵ | 7.87 × 10⁻⁵ | **9.98 × 10⁻⁵** |
| `chr22:20996720 C>G` | 1.4 × 10⁻⁶ | 1.37 × 10⁻⁶ | **2.99 × 10⁻⁵** |

Both remain rare and the filtering decision is unchanged, but it is a systematic
mislabel across the deliverables. It is worth noticing that allele 1's true group
max of 9.98 × 10⁻⁵ sits just under the project's own filter threshold of
1.0 × 10⁻⁴ set in `config/config.yaml`. The correct answer survives that filter
with almost no margin, which is a genuine robustness observation and a better
thing to report than the mislabelled number.

### 2.4 PEX5 and CTU2 are common polymorphisms, and are now closed

This is the finding that most changes a section of the report.

**Claimed**: two "homozygous loss-of-function calls absent from gnomAD",
excluded from the submission, with the exclusion reason "open, not resolved" and
an earlier reason of "mis-calls in repetitive sequence" recorded as disproved by
mapping quality.

**Found**: both are common variants in gnomAD v4.1 exomes, and neither is a
loss-of-function change.

| Call | Change | gnomAD AF | group max | homozygotes in gnomAD |
|---|---|---|---|---|
| `PEX5` 12:7190512 | 45 bp deletion | **0.252** | 0.720 | **173,260** |
| `CTU2` 16:88714226 | 6 bp deletion | **0.767** | 0.925 | **425,713** |

The CTU2 deletion is the *major* allele in gnomAD. A variant carried
homozygously by 425,713 people is not a candidate for a severe recessive
paediatric phenotype, and neither is one carried homozygously by 173,260.

Three further points, each independently sufficient to exclude them:

1. **Neither is a loss-of-function change.** 45 bp and 6 bp are both multiples of
   three, so both are in-frame deletions. Describing them as loss-of-function
   was wrong.
2. **The PEX5 call is not supported by our own alignment at all.** At
   `12:7190512` the panel BAM shows 25 reads spanning the position, mean MAPQ
   58.6, and **zero** reads carrying a deletion or any other indel. The callset
   genotypes it 1/1. The two disagree completely.
3. **The CTU2 deletion is real but not homozygous.** 20 of 30 reads carry the
   6 bp deletion, an allele fraction of 0.67, against a callset genotype of 1/1
   with AD 0,27.

**Why the pipeline said "absent".** It did not. No local gnomAD slice covers
chromosome 16 at all, and none covers the PEX5 locus on chromosome 12, so the
annotator had no frequency for either. Given an empty contig its `_is_covered`
test returns false and it emits `gnomad_not_assayed`, which is the correct and
deliberately separate verdict introduced by commit `fe460cd`. The "absent from
gnomAD" description was added in the write-up, not produced by the analysis. This
is exactly the failure that the `not_assayed` bucket exists to prevent, appearing
downstream of the fix.

**And the reasoning that was overturned was closer to right than what replaced
it.** The original exclusion reason was "mis-calls in repetitive sequence". That
was overturned on the grounds that both loci map uniquely at MAPQ near 60. Unique
mappability does not address repeat-*length* polymorphism, which is what these
are: gnomAD carries a 90 bp allele at the PEX5 locus, twice the 45 bp unit, and a
ladder of alleles from 2 bp to 55 bp at the CTU2 locus. Both sit in tandem
repeats, and repeat-length genotyping from short reads is unreliable in exactly
this way. The correction was right that the mapping-quality argument fails and
wrong to conclude that the repeat argument fails with it.

**Both are now resolved and should be described as closed**, not open.

### 2.5 There is a phasing group in BUB1B

**Claimed**, in the report, the submission notes and `config/config.yaml`:
"HaplotypeCaller's `PGT`/`PID` physical phasing reported no phasing group
anywhere in `BUB1B`".

**Found**: one phasing group exists in the BUB1B region, `40216568_TA_T`,
carried by two records at `15:40216568` and `15:40216570`. Both are genotyped
1/1.

The conclusion is unaffected, for two independent reasons. Neither causal allele
carries a `PGT` or `PID` tag, so no phasing group links them; and homozygous
records carry no *trans* information in any case. The companion analysis in
`results/summaries/compound_het_structure.md` reports zero phasing groups for
BUB1B because it screens heterozygous calls only, which is a defensible method
choice.

The statement as written is nevertheless false, and it is the kind of absolute
claim a reviewer can check in one command. It should say that no phasing group
links the two causal alleles.

### 2.6 The SpliceAI positive controls are eight variants, not nine

**Claimed**: "9/9 known pathogenic canonical splice-site variants score at or
above 0.5". The report states in the same sentence that there were **eight**
controls, so the deliverable contradicts itself.

**Found**: eight control variants, nine gene-level annotations. `15:40218455` is
annotated to both BUB1B and PAK6 and is therefore scored twice. `validate()` in
`scripts/19_arm_b_splicing.py` returns `len(res)`, the number of annotation rows,
and the summary labels it as a number of variants.

Re-running the controls reproduces this exactly and confirms the substance:

| Variants | Annotations | All at or above 0.5 | Minimum delta | Maximum delta |
|---|---|---|---|---|
| 8 | 9 | yes | 0.60 | 1.000 |

The tool is validated and the negative can be believed. Only the count is
mislabelled. The correct phrasing is eight controls, nine annotations, all nine
at or above 0.5.

### 2.7 The strand-balance figures do not sum to their own alternate counts

**Claimed**, in `results/summaries/arm_c_readlevel_verification.md` and carried
into the report:

| | allele 1 | allele 2 |
|---|---|---|
| alt reads | 26 | 12 |
| strand balance | 15 fwd / 12 rev | 6 fwd / 10 rev |

15 + 12 is 27, not 26. 6 + 10 is 16, not 12. Both rows are internally
inconsistent.

**Found**, recomputing from the panel BAM at the stated thresholds of MAPQ ≥ 20
and BQ ≥ 20:

| | depth | ref | alt | VAF | alt fwd/rev | mean MAPQ of alt |
|---|---|---|---|---|---|---|
| allele 1 | 47 | 21 | 26 | 0.553 | **14 / 12** | 60.0 |
| allele 2 | **29** | **16** | **13** | 0.448 | **5 / 8** | 60.0 |
| LZTR1 | 47 | 24 | 23 | 0.489 | 10 / 13 | 60.0 |

Allele 1 agrees on every figure except the forward-strand count. Allele 2 is two
reads deeper in this recomputation, most likely a base-alignment-quality setting
difference. Every conclusion holds: both are heterozygous at a VAF near 0.5,
strand balanced, with mean MAPQ 60.0 on every alternate read. The published
numbers should be replaced with internally consistent ones.

### 2.8 No shortlist artefact is timestamped in the repository

**Claimed**: "Both causal alleles were present in the shortlist ... That artefact
is timestamped in the repository (commit `849bf98`)."

**Found**: nothing under `results/` is tracked by git, by design and by
`.gitignore`. Commit `849bf98` contains `pyproject.toml`,
`scripts/18_arm_a_shortlist.py`, `src/mva/track1/annotators.py` and
`tests/test_annotators.py`. It does not contain the shortlist.

What the repository does hold is the code, and a commit message written at the
time stating the outcome: 12 survivors, of which BUB1B carries three.

The underlying honesty claim holds on the evidence available:

| Time, 31 August 2026 | Event |
|---|---|
| 16:00:04 | commit `849bf98`, the Arm A code and its result described in the message |
| 16:21 | shortlist files written, per their file modification times |
| 16:35:51 | commit `d7c3549`, Arm B runs |
| 16:38 | Hackathon Space captured to `results/rules/` |
| 16:43:02 | commit `355b64c`, RULES.md transcribed from that capture |
| 19:47:43 | commit `0a64650`, the submission built |

The shortlist therefore existed roughly seventeen minutes before the Space was
first captured, and I confirmed by direct membership test that the regenerated
shortlist contains both causal alleles. The ordering supports the claim. The
sentence should say that the *code and its stated result* are timestamped, not
that the artefact is, because a reviewer who looks will not find the artefact.

One caveat I created and should declare: re-running Arm A overwrote
`results/arm_a_shortlist.tsv` and `results/summaries/arm_a_shortlist.md`, so
their modification times now reflect my reproduction run rather than the
original. The 16:21 times above were recorded before I overwrote them, and the
regenerated shortlist is bit-identical to the original.

### 2.9 `config/db_versions.yaml` does not exist

`CLAUDE.md` rule 5 requires that all database versions come from
`config/db_versions.yaml`. No such file existed. Database versions were recorded
in `PROVENANCE.md` instead, which was adequate in substance, but a hard rule
naming a file that was never created is a rule nobody can follow.

**Resolved.** The file now exists. Every value in it was read from a file
header, an API status endpoint or a tool's version output, with the date
recorded, and three entries stay `TODO(source)` because the resource publishes
no version through the interface used. Resources that were never used are listed
under `not_used`, so their absence is a decision rather than an oversight.

---

## 3. Items assessed rather than measured

### 3.1 `proband_id`: leave it as `PROBAND01`

The submission uses `PROBAND01`, taken from the organisers' own
`track1_submission_template.csv`, in which both example rows carry that value.
The field specification says the identifier is "provided in the dataset", and the
only identifier in the dataset is the VCF sample name `WGS_EX2312012`.

I cannot test this, because Track 2 has no leaderboard and I will not probe the
Track 1 one. My assessment is to leave it. The organisers state that this is one
case and not a cohort, which makes the field close to vestigial; the value they
shipped in the template is the value a scorer built against that template is most
likely to accept; and `WGS_EX2312012` is the proband's pseudonymous identifier,
which there is no reason to introduce into a file the organisers will publish.

**This remains the human's call.** If the submission has been made and scored
zero, change this first, exactly as the handover says.

### 3.2 The LZTR1 secondary finding: report it, but drop the modifier framing

Reporting it is right. It is an established pathogenic nonsense variant in a
tumour-suppressor gene, confirmed above against ClinVar, and the organisers
explicitly welcome secondary findings.

The weak part is the sentence proposing that it might modify the presentation,
resting on rhabdomyosarcoma, short stature and failure to thrive overlapping the
RASopathy phenotype. Nothing in this dataset supports a modifying contribution,
rhabdomyosarcoma is not among the three conditions ClinVar lists for this
variant, and a dual diagnosis is being raised only because two findings happen to
be present in one genome. "Cannot be excluded on this data" is true of a great
many things.

Recommendation: keep it as a secondary finding warranting clinical review for
tumour surveillance, and delete the speculation about a modifying contribution.
That is a smaller claim and a defensible one.

### 3.3 Arm C structural variant calling: do not rebuild it

The handover costs a rebuild at roughly seven hours of a fresh booking and
predicts the output is a reported negative over the panel.

My assessment is not to run it, and the reason is stronger than the cost. Panel
coverage is 42-51x with no gene under-covered, both causal alleles are read-level
verified in uniquely mappable sequence, and a compound heterozygote with a
ClinVar-pathogenic nonsense already explains the phenotype. A structural variant
call set would extend the negative-results table and cannot change the call.
Track 2 has no submission yet, three of its required artefacts do not exist, and
scientific rigour is 35% of a panel score there against nothing at all for
completeness here.

The gap is already disclosed in the report, in `RECON.md` and in
`results/summaries/arm_d_mosaic.md`. Disclosed and unfinished is an acceptable
position. Undisclosed would not be.

---

## 4. What was changed as a result

Recorded here rather than folded silently into the documents.

| Document | Change |
|---|---|
| `submission/track1_nexusdwin_report.md` | Allele 2 described as ultra-rare in gnomAD with its frequency, not "novel, absent". ACMG criterion changed to PM2_supporting. The ClinVar record for the same protein change added. Group max frequencies corrected and labelled. Phasing sentence narrowed to the two causal alleles. Positive-control count corrected to eight variants and nine annotations. Strand-balance figures replaced with internally consistent ones. Section 7 rewritten: PEX5 and CTU2 closed as common polymorphisms. LZTR1 modifier speculation removed. The shortlist provenance sentence corrected. |
| `submission/track1_submission.csv` | Notes field corrected on both rows for the same reasons. Coordinates, alleles, `epcr` and `finding_type` are unchanged. |
| `submission/arm_c_readlevel_verification.md` and its copy under `results/summaries/` | Strand-balance and depth figures corrected; PEX5 and CTU2 section replaced with the resolved finding. |
| `config/config.yaml` | `allele_2.gnomad_popmax` corrected from "absent". Phasing note narrowed. |
| `README.md` | The same three corrections, in the summary table and findings list. |
| `scripts/19_arm_b_splicing.py` | Reports controls as variants and annotations separately, so the count cannot be mislabelled again. |
| `tests/` | A regression test that fails if the strand counts in the read-level summary stop summing to the alternate-read counts, and one that fails if a candidate is described as absent from gnomAD while the lookup holds a frequency for it. |

**Nothing was changed about the answer.** The variants, their order, the `epcr`
values and the finding types are exactly as submitted.

---

## 5. Findings from the Track 2 work, recorded here for the same reason

The verification remit was Track 1, but building Track 2 turned up defects of
the same kind in the same places, and burying them in commit messages would
defeat the point of this document.

### 5.1 The safety screen categorically excluded the best chemoprevention agent

`rule_genotoxic` excluded any drug carrying an ATC code beginning `L01`, on the
stated grounds that cytotoxic chemotherapy is out of scope. Running the screen
over a real candidate set for the first time showed what that costs. **Celecoxib
carries `L01XX33`** alongside `M01AH01`, because of its familial adenomatous
polyposis indication. The rule excluded the single best-evidenced chemoprevention
agent in hereditary cancer predisposition, and excluded it on the reasoning that
a COX-2 inhibitor is cytotoxic chemotherapy.

The exclusion is now scoped to `L01A` to `L01D`, the cytotoxic subgroups, with
the boundaries read from the ChEMBL ATC endpoint and cached at
`refs/atc/atc_l01.json` rather than recalled. `L01E`, `L01F` and `L01X` are
flagged instead, so an antineoplastic classification still travels with a
candidate and never passes silently.

**The general lesson is that a rule written against an imagined candidate set
was wrong in a way that only a real one could reveal.** It had passed its tests,
because the tests were written from the same imagination.

### 5.2 Our own headline Track 2 claim was overstated

The direction audit found ten targets, all requiring activation, and no
activating drug among 118 ChEMBL mechanism records. That was reported as though
it were a property of the spindle assembly checkpoint.

Measuring the base rate shows it is mostly a property of pharmacology. Only 359
of 19,297 approved protein-coding genes have any activating drug, so at the
genome-wide rate the expected count among ten targets is 0.19 and **the
probability of observing zero is 0.83**. The observation is unremarkable.

The claim is now made at the strength the evidence supports, in section 3.4 of
the Track 2 report. Three narrower statements survive and they are enough.

### 5.3 Three pipeline defects, each of which produced a confident wrong number

- **A cache-key collision.** A count query asking for one study poisoned the
  cache for the full query behind it, so the pipeline derived candidates from a
  single trial while printing the true total of 39 beside them. The output
  looked internally consistent and was wrong.
- **Salt forms escaped the safety screen.** ChEMBL attaches no ATC code to a
  salt, so `ERLOTINIB` was flagged and `ERLOTINIB HYDROCHLORIDE` was allowed.
  Same active molecule, same concern, opposite verdict. ATC is now inherited
  from the ChEMBL parent.
- **An empty lookup reported as no lookup.** Agents whose paediatric trial search
  returned nothing were described as "paediatric exposure not looked up" when it
  had been looked up and found nothing. Those are different statements and the
  screen distinguishes them.

### 5.4 A heuristic of ours was confidently wrong before it shipped

A first attempt to flag non-tumour prevention endpoints tested keywords against
each trial's **condition** list. It reported that atorvastatin in Lynch syndrome
was not a chemoprevention candidate, because the string "Lynch Syndrome" contains
no tumour word. The test was replaced rather than patched: classification now
reads the trial's **primary outcome** text, and the output prints that text
beside our classification of it so a reader can overrule us.

### 5.5 What this run added to prevent recurrence

| Guard | Prevents |
|---|---|
| `tests/test_track2_report_matches_data.py` | Every candidate row, ChEMBL identifier, trial count and verdict in the Track 2 report is checked against the generated summary, and every cited NCT identifier must appear in the pipeline's own output. |
| Dose assertions on the real output files | A dose reaching any artefact, which `CLAUDE.md` rule 3 forbids. |
| Identifier-shape assertions | A malformed accession, which is usually an invented one. |
| Negative control in the chemoprevention run | A safety screen that excludes nothing because it is broken rather than because the candidates are clean. Five of five known cytotoxics are excluded. |
| `binom_zero` in the library with a test | The base-rate arithmetic that weakens our own claim being wrong in the direction that flatters us. |
| Pitch-length test | A three-minute video script that is not three minutes. The first draft ran to 593 spoken words, nearly four minutes. |
| `tests/test_track1_report_matches_data.py` | The Track 1 report drifting from the callset, from the verification artefacts, or from itself. One test per corrected claim in section 2. |
| `tests/test_track2_needs_no_patient_data.py` | A convenience read of a patient file quietly destroying the property that lets a reviewer reproduce Track 2 without data access. |

### 5.6 The correction in section 2.7 was applied to only one of two documents

Recorded because it is the same class of error the section documents.

Section 2.7 found that the published strand counts did not sum to their own
alternate-read counts, and recorded that the figures had been replaced with
internally consistent ones. They were replaced in
`arm_c_readlevel_verification.md` and only half replaced in the Track 1 report,
where the strand row was corrected while the depth, read and VAF rows kept the
superseded values. The report therefore carried 5 forward and 8 reverse
alternate reads beside an alternate count of 12, which is the exact
inconsistency section 2.7 exists to describe.

It was found by `tests/test_track1_report_matches_data.py` within a minute of
that test first running, which is the argument for the test. Hand verification
does not repeat itself; a test does. The report now reads depth 29, 16 reference
and 13 alternate reads, and VAF 0.448, matching the recomputation.

---

## 6. History rewrite, 3 September 2026

### 6.1 What prompted it, and what it found first

The open item said commit `4ce0d4c` contained cluster hostnames that had been
sanitised but remained in history. Preparing the rewrite showed that was only
half the problem. **Both real hostnames were still in the working tree**, in
`scripts/gpu/.local.sh.example`, which is tracked and public.

The August sanitisation commit replaced the account name and the SSH key path
with placeholders and left `JUMP_HOST` and `GPU_NODE` as real values. Its message
stated that the repository "was carrying the GPU node and jump-host names, the
account name and the SSH key path" and that topology now came from a gitignored
file. Two thirds of that was accurate, and nobody had checked in the three days
since. **A commit message asserting that something was removed is not evidence
that it was removed.**

### 6.2 What was redacted

Seven distinct strings: two university hostnames, one account name and four
paths derived from it. The list was not written by hand. It was derived from the
lines the sanitisation commit itself removed, then filtered to strings that
appear **only** under `scripts/gpu/`, so that no generic token could be replaced
globally and corrupt unrelated content. Nothing was rejected by that filter,
which is itself a check: had a common word appeared in the candidate list, it
would have been excluded rather than replaced.

| | |
|---|---|
| Tool | `git filter-repo --replace-text` |
| Commits rewritten | 54, all of them |
| Occurrences replaced | 154 across 7 files, all under `scripts/gpu/` |
| Real hostnames remaining in any reachable commit | **0**, verified by exhaustive scan |

### 6.3 What survived, which matters for the honesty claim

Only two commits changed identity, both after the topology was introduced:

| Cited as | Now |
|---|---|
| `5214f06` | `0a64650` |
| `4ce0d4c` | `24e6ae7` |

**The four commits cited as evidence elsewhere are unchanged**, because they
predate the affected files: `849bf98`, `d7c3549`, `355b64c` and `fe460cd`. That
includes `849bf98`, which the Track 1 report cites as the timestamp for the
claim that both causal alleles were shortlisted before the leaderboard was seen.
Had that SHA moved, the claim would have needed restating rather than
renumbering.

The two changed references have been updated in this document, in
`docs/HANDOVER_PROMPT.md` and in `submission/README.md`.

### 6.4 What this does not accomplish

Stated because a rewrite invites the belief that the exposure is undone, and it
is not.

- ~~The remote still has the old history until it is force-pushed.~~ **Done,
  3 September 2026.** `main` was force-pushed with an explicit lease on the
  known remote value, and the result was verified by cloning the public
  repository fresh and scanning every reachable commit: zero real hostnames.
  The `pre-sanitise-backup` tag was deliberately not pushed.
- **Anyone who cloned or forked the repository keeps the old objects.** Rewriting
  our copy does not reach theirs.
- **GitHub retains unreachable objects** and serves them by SHA for a period.
  Removing them entirely requires asking GitHub Support after the force-push.
- **The backup retains everything.** A bundle of the pre-rewrite history was
  taken before the rewrite and a `pre-sanitise-backup` tag was created. Both
  contain the unredacted values and both should be destroyed once the rewrite is
  confirmed good.

None of the redacted strings is a credential. They are third-party university
infrastructure names, which is why this is a tidy-up rather than an incident.

### 6.5 The guard

`tests/test_no_topology_leaks.py` asserts the property rather than trusting a
commit message: no tracked file may carry an `.ac.uk` hostname other than the EBI
API host and documentation placeholders, all four variables in the example file
must look like placeholders, and no reachable commit may carry one either. That
last test failed before the rewrite and passes after it.
