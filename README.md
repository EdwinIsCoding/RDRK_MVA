# MVA Hackathon 2026

Submission for **Rare Disease, Real Kid: The MVA Hackathon 2026**
(Sage Bionetworks, MVA Society, HuggingFace, BEACON).

- **Track 1** predicts the causal variant or variants in a proband with
  suspected mosaic variegated aneuploidy.
- **Track 2** proposes drug repurposing hypotheses against the disrupted
  biology.

Deadline 24 October 2026, code freeze 17 October 2026. All outputs CC BY 4.0.

## Status: Track 1 solved. Submission built, awaiting upload.

**Answer: a compound heterozygote in `BUB1B`**, causing mosaic variegated
aneuploidy syndrome 1 (OMIM 257300).

| Allele | Change | Evidence |
|---|---|---|
| `chr15:40209701 T>G` | `c.2210T>G` **p.Leu737Ter**, nonsense | ClinVar **VCV000533901.9** Pathogenic/Likely pathogenic, listed against MVA1. gnomAD popmax 7.9e-05 |
| `chr15:40220612 T>G` | `c.3006T>G` **p.Asn1002Lys**, missense | **Novel**, absent from gnomAD, SIFT 0.01, PolyPhen 0.997 |

Both verified at read level against **our own alignment**, built from raw FASTQ
independently of the supplied callset: VAF 0.553 and 0.444, strand balanced,
mean MAPQ 60.0 on every alternate read.

**Phase is inferred, not demonstrated.** The alleles lie 10,911 bp apart, beyond
a read pair, and no `PGT`/`PID` phasing group exists in `BUB1B`. Confirmation
needs parental testing or long reads.

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

### Findings

Each changed a decision, and each was measured rather than assumed.

- **The plan's central hypothesis was wrong.** It reasoned that a cryptic splice
  allele was most likely, since a standard coding pipeline would already have
  solved the case. SpliceAI at plus or minus 500 bp found nothing above even the
  permissive threshold. The second allele is an ordinary novel missense.
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
- **No consanguinity**, confirmed by two independent methods.
- **No coverage gap** over any MVA gene: 42-51x against a genome mean of 43.8x.

### Method validation

- SpliceAI silently returned 0.000 for everything, including eight known
  pathogenic canonical splice variants, because of a NumPy 2 compatibility shim
  of ours. Corrected, those controls score **9/9** above 0.5. The runner now
  refuses to report a negative if its positive controls fail.
- Splice distances agree with ClinVar HGVS intron offsets for **268/268 SNVs**.
- **140 automated tests** across the evidence schema, scoring, annotators and
  submission format.

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
