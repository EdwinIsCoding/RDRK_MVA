# Submission artefacts

These are deliverables, not intermediate results, so they are tracked whereas
everything under `results/` is not.

Committing variant coordinates here is deliberate and consistent with the rules.
The organisers' data-handling guidance lists "candidate variant rankings" among
the things participants may keep, as distinct from VCF and BAM files, variant
tables with genotypes, and prompts containing pasted variant blocks, which must
be deleted. The submission itself is published under CC BY 4.0 and the
organisers state that all submitted results are made accessible to the research
community.

| File | Purpose |
|---|---|
| `track1_submission.csv` | The Track 1 file to upload. Format validated by `scripts/23_track1_submission.py`. |
| `track1_candidates.tsv` | Human-readable source with the reasoning per row. |
| `track1_nexusdwin_report.md` | The Track 1 methods write-up, judged separately. |
| `arm_c_readlevel_verification.md` | Read-level confirmation from our own alignment. |
| `arm_d_mosaic.md` | The mosaic arm, reported negative. |
| `track2_nexusdwin_report.md` | **The Track 2 report.** Mechanism characterisation, the closed direct axis, the chemoprevention axis and its safety screen. |
| `track2_nexusdwin_pitch.md` | **The 3-minute pitch script**, with production constraints. The video itself is not yet recorded. |

## Track 2 submission checklist

The rules require three artefacts. Two exist.

| Required | Status |
|---|---|
| Written report, PDF or Markdown, filename carrying the team name | **Done**, `track2_nexusdwin_report.md` |
| Public GitHub repository with reproducible code | Repository exists and is public. History was rewritten and force-pushed on 3 September 2026 to purge third-party cluster hostnames and an account name. Verified by cloning the public repository and scanning all history: zero remain. See `docs/VERIFICATION.md` section 6. |
| 3-minute pitch video on YouTube or Vimeo | **Not done.** The script is written and timed at 179 seconds; the recording and upload remain. |

Only 3 Track 2 submissions are permitted per team and only the latest is
reviewed, so there is no value in submitting early to probe a score.

## Verification performed before submitting

- Coordinates read from the **original callset**, not from the VEP output.
- REF alleles checked against the GRCh38 reference: `15:40209701` T, `15:40220612` T,
  `22:20996720` C all match.
- Chromosomes chr-prefixed, since the submission requires it and the callset
  does not use it.
- The compound-heterozygous pair is in a **single row** using the `_2` columns,
  which is what earns the partial credit the metric offers for recovering one of
  two variants.
- Format gate passed. It had already rejected an earlier hand-written file in
  which one missing tab put `epcr` in the `alt_2` column.
