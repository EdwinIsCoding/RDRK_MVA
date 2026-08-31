# Project: MVA Hackathon 2026

Rare Disease, Real Kid: The MVA Hackathon 2026. Sage Bionetworks, MVA Society,
HuggingFace, BEACON. The governing document is `MVA_HACKATHON_PLAN.md`. Phase 0
is complete; its findings are in `DATA_CARD.md` and `RECON.md`.

## Hard rules

1. **NEVER read, `cat`, `head`, `grep` or otherwise load patient genomic data
   into context.** `data/` is off-limits. Work from `results/summaries/` and
   `results/recon/` only. Scripts under `scripts/` and `src/` may touch `data/`;
   you may write and run them, you may not read their inputs.
   The one documented exception is the clinical phenotype document, which was
   read at Phase 0 under the reasoning recorded in `ETHICS.md` section 3. It is
   not genomic data. That exception does not generalise.
2. **NEVER invent an identifier.** rsIDs, PMIDs, ClinVar VCV accessions, ChEMBL
   IDs, NCT numbers, HPO terms, coordinates and gene symbols come from a tool
   call or a file, never from memory. If you do not have a source, write
   `TODO(source)` and stop. A plausible-looking wrong accession is worse than a
   gap, because a reader cannot tell it is wrong.
3. **NEVER generate dosing, dose ranges, safety margins or clinical
   recommendations.** Outputs are research hypotheses addressed to researchers.
4. **Genome build is GRCh38 with `ensembl_nochr` contig naming** (`1`, not
   `chr1`), set in `config/config.yaml`. Every coordinate is build-tagged.
   Reject any function signature taking a bare position. Most annotation
   resources ship with the `chr` prefix, so a rename step in front of VEP,
   SpliceAI, AnnotSV and similar is mandatory, not optional.
5. **All randomness seeded** from `config/config.yaml` (`analysis.random_seed`).
   All database versions from `config/db_versions.yaml`.
6. **Halt at every STOP checkpoint** in `MVA_HACKATHON_PLAN.md` and report. Do
   not proceed past one on your own judgement.
7. **British English, no em dashes**, in every file and every message.
8. **No `Co-Authored-By` trailer and no Claude attribution** in commit messages
   or pull request bodies.

## What Phase 0 established (do not re-derive)

- Singleton, `WGS_EX2312012`. No parents. De novo and segregation models are off.
- GRCh38, no-chr naming, 2,580 contigs including hs38d1 decoys.
- 5,012,204 records, 94.58% PASS, median depth 42×. Filtering was
  non-destructive, so 271,414 filtered records are recoverable. `MQ40` alone is
  177,522 of them and is a search target, not noise.
- `PGT`/`PID` physical phasing tags are present and are the only native phasing
  signal available.
- No RNA-seq. No BAM or CRAM, but 84.7 GB of paired FASTQ across four lanes.
- No karyotype, no aneuploidy percentage, one tissue.
- Eight HPO terms, already coded by the organisers. No extraction needed.
- No consanguinity: longest ROH-like run is 2 Mb across 2,777 windows.
- No coverage collapse over any known MVA gene: depth ratio to flank is
  0.86-1.07. Variant-free gaps in `BUB1B`, `CEP57` and `TRIP13` were calibrated
  against a genome-wide background of 37-85 such runs per chromosome and are
  **not** leads.

## Anti-patterns, do not implement

- Structure-prediction RMSD between wild type and mutant as a variant-effect
  score. Use FoldX, Rosetta `cartesian_ddg`, ThermoMPNN or RaSP, and state the
  r≈0.5-0.7 accuracy honestly.
- Docking scores or DiffDock confidence presented as binding affinity. Docking
  is an enrichment tool with roughly 2 kcal/mol error.
- Network proximity used to infer whether to inhibit or activate a target.
  PrimeKG and Hetionet edges are largely unsigned. Use SIGNOR, OmniPath or
  Reactome for direction.
- Any knowledge-graph claim without a time-split leakage check.
- An LLM critic agent on the same base model used as a validity check.
- Filtering the callset to `PASS` before Arm A has looked at what was filtered.
- Presenting a splicing prediction as though it were an observed junction. There
  is no RNA-seq. The distinction must survive into the final report.

## Style

- Scoring functions return `(score, evidence: list[Evidence])` where each
  `Evidence` carries a resolvable source identifier. Never a bare float.
- Snakemake rules over ad-hoc scripts. Declare inputs and outputs.
- Positive-control tests pass before any ranking logic is trusted.
- Every recon or analysis script writes aggregates to `results/`, never
  individual patient genotypes.
- Report negative results. Closing a lead cheaply is worth as much as opening
  one, and the plan's judging criteria reward it explicitly.
