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
| `track1_submission.csv` | The file to upload. Format validated by `scripts/23_track1_submission.py`. |
| `track1_candidates.tsv` | Human-readable source with the reasoning per row. |

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
