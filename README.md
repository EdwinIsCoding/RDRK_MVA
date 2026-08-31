# MVA Hackathon 2026

Submission for **Rare Disease, Real Kid: The MVA Hackathon 2026**
(Sage Bionetworks, MVA Society, HuggingFace, BEACON).

- **Track 1** predicts the causal variant or variants in a proband with
  suspected mosaic variegated aneuploidy.
- **Track 2** proposes drug repurposing hypotheses against the disrupted
  biology.

Deadline 24 October 2026, code freeze 17 October 2026. All outputs CC BY 4.0.

## Status: Phase 0-2 built. STOP #1 cleared, STOP #2 partial.

| Document | Contents |
|---|---|
| `MVA_HACKATHON_PLAN.md` | The governing plan. Six analysis arms, four STOP checkpoints. |
| `DATA_CARD.md` | What the data is, what state it is in, and what it will and will not support. |
| `RECON.md` | Hypothesis class triage and the confirmed branch. |
| `STOP2_STATUS.md` | **Why STOP #2 cannot yet be cleared**, and what was measured instead. |
| `ETHICS.md` | Consent scope, what we did and did not do, and one documented judgement call. |
| `RULES.md` | **Incomplete.** Hackathon rules transcription, needed before the first analysis commit. |
| `PROVENANCE.md` | Input checksums, tool versions, database snapshot dates. |
| `CLAUDE.md` | Agent contract. Hard rules and anti-patterns. |

**Confirmed branch: C, deferred.** Singleton, GRCh38, no RNA-seq, no alignments
shipped but FASTQ present. Arms A, B, E and part of F are live now; Arms C, D and
repeat expansion detection are blocked behind building a BAM from FASTQ.

### Findings so far

Each of these changed a design decision, and each was measured rather than assumed.

- **The benchmark cannot test the leading hypothesis.** Of 108 confidently
  pathogenic MVA-gene variants in ClinVar, none is deep intronic, near-splice,
  synonymous or UTR. `STOP2_STATUS.md` sets out the consequence and the
  substitute control set.
- **Every known MVA gene is unconstrained** (BUB1B LOEUF 0.707 pLI 0.000,
  CENATAC 1.227, and so on). Weighted conventionally, constraint would have
  deprioritised every correct answer.
- **The direct spindle-checkpoint-restoration axis is pharmacologically
  unavailable.** Signed-edge nomination yields ten targets all requiring
  activation; ChEMBL holds 118 drug mechanisms across them, none activating.
  See `results/summaries/track2_direction_audit.md`.
- **No consanguinity** in the proband (longest ROH-like run 2 Mb), which keeps
  the compound heterozygote as the leading model.
- **Splice distance is exact for SNVs** (268/268 against ClinVar HGVS) after
  restricting exon boundaries to the MANE Select transcript.

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
