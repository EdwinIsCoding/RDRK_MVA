# MVA Hackathon 2026

Submission for **Rare Disease, Real Kid: The MVA Hackathon 2026**
(Sage Bionetworks, MVA Society, HuggingFace, BEACON).

- **Track 1** predicts the causal variant or variants in a proband with
  suspected mosaic variegated aneuploidy.
- **Track 2** proposes drug repurposing hypotheses against the disrupted
  biology.

Deadline 24 October 2026, code freeze 17 October 2026. All outputs CC BY 4.0.

## Status

**Track 1 solved.** Submission built, verified against primary sources, and
independently reviewed. Awaiting upload.
**Track 2 report and pitch script written**, pitch video not yet recorded.

Everything outstanding is listed in one place, `submission/README.md`, together
with the seven declared scientific limitations and which of them are closed.

**Answer: a compound heterozygote in `BUB1B`**, causing mosaic variegated
aneuploidy syndrome 1 (OMIM 257300).

| Allele | Change | Evidence |
|---|---|---|
| `chr15:40209701 T>G` | `c.2210T>G` **p.Leu737Ter**, nonsense | ClinVar **VCV000533901.9** Pathogenic/Likely pathogenic, listed against MVA1. gnomAD v4.1 exomes AF 7.9e-05, group max 1.0e-04 |
| `chr15:40220612 T>G` | `c.3006T>G` **p.Asn1002Lys**, missense | **Ultra-rare**: gnomAD v4.1 exomes AF 6.8e-07, one allele in 1,461,878, no homozygotes. No ClinVar record for this nucleotide change; the same protein change via `c.3006T>A` is **VCV004600147.1**, Uncertain significance. **16 predictors: 10 damaging, 5 tolerated, 1 intermediate.** AlphaMissense pathogenic 0.88-0.92, CADD 24.5, REVEL 0.472. Favours a damaging effect; does not establish one |

Both recovered by three analysis routes over one library: the supplied Sentieon
callset, our own `bwa-mem2` alignment, and GATK Mutect2 over that alignment.
VAF 0.553 and 0.448, strand balanced, mean MAPQ 60.0 on every alternate read.
That is robustness to aligner and caller, not independent confirmation: all
three share one FASTQ library, so a sample swap or contamination would be
reproduced by all of them.

**Phase is inferred, not demonstrated.** The alleles lie 10,911 bp apart,
neither carries a `PGT`/`PID` phasing tag, and no read or read pair in our
alignment touches both positions. Confirmation needs parental testing or long
reads.

| Document | Contents |
|---|---|
| `submission/` | The Track 1 predictions file and report, ready to upload |
| `MVA_HACKATHON_PLAN.md` | The governing plan. Six analysis arms, four STOP checkpoints |
| `DATA_CARD.md` | What the data is and what it supports |
| `RECON.md` | Hypothesis triage and the confirmed branch |
| `STOP2_STATUS.md` | Why the benchmark could not test the leading hypothesis |
| `ETHICS.md` | Consent scope, the LLM data-handling audit, and a disclosed compliance gap |
| `RULES.md` | Hackathon rules, transcribed from the official Space |
| `PROVENANCE.md` | Input checksums, tool versions, database snapshot dates |
| `CLAUDE.md` | Agent contract. Hard rules and anti-patterns |
| `docs/VERIFICATION.md` | Adversarial re-check of every Track 1 claim against primary sources, and the nine corrections it produced |
| `submission/track2_nexusdwin_report.md` | The Track 2 report. Mechanism, the closed direct axis, chemoprevention and the safety screen |
| `submission/track2_nexusdwin_pitch.md` | The 3-minute pitch script |
| `resources/directional_availability/` | Reusable: every gene ChEMBL can push in each direction, published for other teams |
| `submission/README.md` | The submission checklist, the seven declared limitations with status, and everything outstanding |
| `docs/sage_disclosure_draft.md` | Draft disclosure to Sage Bionetworks' Privacy and Compliance Office, for the owner to send |
| `config/db_versions.yaml` | Every database version, with the date and source each was read from |

### Findings

Each changed a decision, and each was measured rather than assumed.

- **The plan's central hypothesis was wrong.** It reasoned that a cryptic splice
  allele was most likely, since a standard coding pipeline would already have
  solved the case. SpliceAI at plus or minus 500 bp found nothing above even the
  permissive threshold. The second allele is an ordinary ultra-rare missense.
- **The benchmark cannot test the hypothesis it was built for.** Of 108
  confidently pathogenic MVA-gene variants in ClinVar, none is deep intronic,
  near-splice, synonymous or UTR. See `STOP2_STATUS.md`.
- **Every known MVA gene is unconstrained** (BUB1B LOEUF 0.707, pLI 0.000).
  Weighted as a dominant-disease pipeline would, constraint would have actively
  deprioritised the correct answer.
- **The proband is male**, determined from chrX heterozygosity of 0.062 against
  0.620 autosomal. Not stated in the clinical document, and it makes X-linked a
  single-hit hypothesis.
- **The direct therapeutic axis is pharmacologically unavailable.** Signed-edge
  nomination yields ten targets all requiring activation; ChEMBL holds 118 drug
  mechanisms across them and none is activating.
- **We then weakened that finding ourselves.** Only 359 of 19,297 human genes
  have an activating drug at all, so seeing zero across ten targets has
  probability 0.83 under the base rate. The claim we make is the narrower one.
- **The chemoprevention evidence base for this child is empty**: zero registered
  trials in mosaic variegated aneuploidy, zero naming `BUB1B`, zero prevention
  trials in rhabdomyosarcoma. The same queries return 150 trials for Fanconi
  anaemia, so the emptiness is the disease and not the question.
- **We tried to break our own headline and it survived, with its scope
  narrowed.** All-activation across ten targets has probability 0.024 against a
  measured stimulatory share of 69%, so it is not an artefact of the graph. But
  it belongs to the immediate regulators of the eight seed genes rather than to
  mitotic biology, and where a wider net does reach inhibition-reachable
  targets, those routes activate a mitotic kinase and the safety screen rules
  them out anyway.
- **Signature reversal, the method the plan calls highest-yield, is now run**
  against a labelled BubR1-hypomorph proxy. A null control puts the hits at 91%
  specific, and the safety screen excluded daunorubicin from the top of the
  list.
- **The method discriminates between near-neighbour diseases.** Run unchanged on
  Fanconi anaemia and ataxia-telangiectasia, it returns 17 and 10 targets
  reachable by inhibition against 0 for this proband, which is what makes the
  closed axis a finding rather than a habit.
- **The mitochondrial axis looked promising and did not survive follow-through.**
  It was the only axis better supplied with activating drugs than the genome
  average, but roughly 70% of that supply is `INSR`, `PPARA` and `GCK`, whose
  drugs are insulins, fibrates and glucokinase activators swept in by a GO
  annotation. Three axes are now closed and chemoprevention is the one standing.
- **No structural variant over any known MVA gene**, from a breakpoint screen
  calibrated against 400 sampled panel regions. Uncalibrated it flagged seven of
  nine genes, which was the background rate of split-read artefact.
- **No credible mosaicism, on two callers that fail in opposite directions.**
  GATK Mutect2 tumour-only agrees with `bcftools mpileup`, and recovers both
  causal alleles at PASS as its own positive control.
- **No consanguinity**, confirmed by two independent methods.
- **No coverage gap** over any MVA gene: 42-51x against a genome mean of 43.8x.
- **Two variants we had set aside as unresolved are common polymorphisms.** The
  `PEX5` and `CTU2` homozygous calls are carried homozygously by 173,260 and
  425,713 people in gnomAD. We had described them as absent from gnomAD, from a
  lookup that never assayed their chromosomes. See `docs/VERIFICATION.md`.

### Method validation

- SpliceAI silently returned 0.000 for everything, including eight known
  pathogenic canonical splice variants, because of a NumPy 2 compatibility shim
  of ours. Corrected, all eight score above 0.5. They produce nine gene-level
  annotations, because one control is annotated to both `BUB1B` and `PAK6`, and
  all nine are at or above 0.5. The runner now refuses to report a negative if
  its positive controls fail.
- Splice distances agree with ClinVar HGVS intron offsets for **268/268 SNVs**.
- **366 automated tests** across the evidence schema, scoring, annotators,
  submission format, and every claim corrected during verification. The suite
  passes clean.

## Reproducing Track 2 without any data access

**Track 2 reads no patient data at all.** Every result in
`submission/track2_nexusdwin_report.md` comes from public databases: ChEMBL,
ClinicalTrials.gov, OmniPath, QuickGO, UniProt, PDBe and HGNC. No script in the
Track 2 pipeline opens a file under `data/`, and a test asserts it
(`tests/test_track2_needs_no_patient_data.py`).

So a reviewer can check every Track 2 claim without applying for the challenge
dataset:

```bash
git clone <this repository> && cd RDRK_MVA
pip install -e .            # or: micromamba create -f environment.yml
make reproduce-track2       # one 16 MB download, then the whole pipeline
```

That fetches HGNC, which supplies the denominator for the base-rate argument,
runs the direction audit, the chemoprevention axis, the availability
measurement, the three-disease generalisation check and the structural
feasibility check, publishes the reusable availability table, and runs the test
suite. Outputs land in `results/summaries/` and should match the tables in the
Track 2 report.

Two of the sources are live registries. `make track2-drift` compares today's
counts against the ones pinned in `config/track2_evidence_pin.json` on the day
the report was written, so a reviewer can tell whether a difference is our error
or the world moving.

**Track 1 is different** and does need the challenge distribution, because it is
an analysis of this proband's genome. `make reproduce` is that path.

## Data

Not in this repository and never will be. `data/` is gitignored and a pre-commit
hook hard-fails on anything under it, on genomic file extensions anywhere, and on
files above 5 MB. See `ETHICS.md`.

## Reproducing Phase 0

Requires `bcftools` (1.24 used here) and `pandoc`, with the dataset in `data/`.

```bash
git config core.hooksPath .githooks     # install the guard first

scripts/00_inventory.sh                 # manifest, checksums, headers, file types
scripts/02_extract_clinical.sh          # phenotype document to text
python3 scripts/01_characterise.py      # routing facts -> results/recon/characterisation.json

# Verification and calibration. Each closes a question that would otherwise
# have been guessed at.
python3 scripts/03_roh_proxy.py           # consanguinity: 1 Mb heterozygosity scan
python3 scripts/04_verify_panel_tally.py  # is the shipped .tbi index trustworthy?
python3 scripts/05_panel_depth_profile.py # 10 kb depth and density across panel genes
python3 scripts/06_gap_background_rate.py # background rate of variant-free runs
```

`SKIP_SHA=1` skips checksumming, which otherwise reads all 79 GB.

## Layout

```
config/          config.yaml (build, cohort, branch, limitations), gene panels
scripts/         Phase 0 recon and verification. These may touch data/.
src/mva/         io, track1, track2, llm, report
tests/           positive control harness (Phase 2)
benchmarks/      published MVA causal variants (Phase 2)
results/recon/   Phase 0 outputs. Gitignored where they contain anything specific.
results/summaries/  aggregate summaries, the only directory an agent reads freely
```
