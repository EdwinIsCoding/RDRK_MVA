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

## The answer (established 1 September 2026)

`BUB1B` compound heterozygote causing MVA1 (OMIM 257300):
`chr15:40209701 T>G` p.Leu737Ter (ClinVar VCV000533901.9, Pathogenic/Likely
pathogenic) with `chr15:40220612 T>G` p.Asn1002Lys (ultra-rare, not absent:
gnomAD v4.1 exomes group max 8.99e-07, a single allele in 1,461,878, and absent
from gnomAD genomes only. The same protein change via c.3006T>A is ClinVar
VCV004600147.1, uncertain significance). Both read-level verified: VAF 0.553 and
0.448, strand balanced, MAPQ 60.
Phase is **inferred, not demonstrated**: 10,911 bp apart, no PGT/PID phasing
group. Secondary finding: `LZTR1` chr22:20996720 p.Tyr748Ter, ClinVar
VCV001409252.7.

Track 1 is effectively complete. **The remaining work is the write-ups and
Track 2**, where the competition actually is: 61 teams were already tied at a
perfect Track 1 score before we submitted.

## What Phase 0 established (do not re-derive)

- Singleton, `WGS_EX2312012`. No parents. De novo and segregation models are off.
- GRCh38, no-chr naming, 2,580 contigs including hs38d1 decoys.
- 5,012,204 records, 94.58% PASS, median depth 42×. Filtering was
  non-destructive, so 271,414 filtered records are recoverable. `MQ40` alone is
  177,522 of them and is a search target, not noise.
- `PGT`/`PID` physical phasing tags are present and are the only native phasing
  signal available.
- No RNA-seq. No BAM or CRAM shipped, but 84.7 GB of paired FASTQ; an alignment
  was built from it on 1 September (61 GB BAM, bwa-mem2, 4h10m).
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
