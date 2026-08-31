# DATA_CARD.md

Phase 0 characterisation of the challenge dataset. Written by
`scripts/00_inventory.sh`, `scripts/01_characterise.py` and the verification
scripts `03`-`06`. Machine-readable source: `results/recon/characterisation.json`.

Generated 31 August 2026. Tooling: bcftools 1.24 / htslib 1.24, pandoc.

---

## 1. Dataset at a glance

| Property | Value |
|---|---|
| Proband sample ID | `WGS_EX2312012` |
| Cohort structure | **Singleton.** One sample in the VCF. No parents, no siblings. |
| Assay | Whole genome sequencing, Illumina, paired-end |
| Flowcell / lanes | `HGWCNDSX7`, sample index S16, lanes L001-L004 |
| Reference build | **GRCh38** (chromosome 1 = 248,956,422 bp) |
| Reference assembly | `GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta` |
| Contig naming | **Ensembl style, no `chr` prefix** (`1`, `2`, … `X`, `Y`, `MT`) |
| Contigs in header | 2,580 (primary assembly plus hs38d1 decoys) |
| Total input volume | 79 GB across 11 files |

### File manifest

| File | Bytes | Role |
|---|---|---|
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` | 315,153,971 | Small variant callset |
| `WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi` | 2,343,376 | Tabix index |
| `WGS_EX2312012_HGWCNDSX7_S16_L00{1..4}_R{1,2}_001.fastq.gz` | 8 files, 84.7 GB total | Raw reads, 4 lanes, paired |
| `Challenge_Clinical_Phenotype_1.docx` | 16,865 | Proband phenotype, HPO-coded |

SHA256 checksums for all files are in `results/recon/sha256.txt` and mirrored
into `PROVENANCE.md`.

**No BAM or CRAM was shipped.** The VCF header shows one existed upstream
(`WGS_EX2312012_HGWCNDSX7.bam`) but it was not distributed. Alignments must be
regenerated from the FASTQ. See section 6.

**No RNA-seq of any kind.** No FASTQ read group, counts matrix or expression
file indicates a transcriptome assay.

---

## 2. Variant calling provenance

The callset was produced on DNAnexus in February 2025:

| Stage | Tool | Version | Date |
|---|---|---|---|
| Per-sample calling | Sentieon Haplotyper (GATK HaplotypeCaller equivalent), GVCF mode | sentieon-genomics-202308.02 | 2025-02-04 |
| Joint genotyping | Sentieon GVCFtyper | sentieon-genomics-202308.02 | 2025-02-05 |
| Filtering | GATK VariantFiltration | 4.2.4.0 | 2025-02-05 |

dbSNP build 138 (`Homo_sapiens_assembly38_no_chr.dbsnp138.vcf.gz`) supplied the
`DB` flag. Calling confidence thresholds were `--call_conf 10 --emit_conf 10`,
which is permissive and retains more low-quality and low-allele-fraction
evidence than the GATK default of 30. That is favourable for the mosaic arm.

### Filtering is hard-threshold, not VQSR, and is non-destructive

`VariantFiltration` applied five hard filters and **tagged** rather than removed
the failing records. Every filtered variant is still present in the file and is
recoverable.

| FILTER | Records | Share |
|---|---|---|
| `PASS` | 4,740,790 | 94.58% |
| `MQ40` (RMS mapping quality < 40) | 177,522 | 3.54% |
| `QD2` (quality by depth < 2) | 35,402 | 0.71% |
| `LowQual;QD2` | 17,000 | 0.34% |
| `MQ40;QD2` | 12,541 | 0.25% |
| all remaining combinations | 28,949 | 0.58% |
| **Total non-PASS** | **271,414** | **5.42%** |

`MQ40` alone accounts for two thirds of the discarded records. Low mapping
quality concentrates in segmental duplications, repeats and paralogous
sequence, which is precisely where a diagnostic pipeline under-calls and where
hypothesis classes 2 (structural) and 6 (regulatory) live. **Arm A must not
filter to `PASS` blindly.**

No VQSR or CNN recalibration was applied, so there is no `VQSLOD` or
`CNN_1D/2D` score to inherit or reason about.

---

## 3. Callset quality

| Metric | Value | Interpretation |
|---|---|---|
| Total records | 5,012,204 | Normal for single-sample WGS |
| SNPs | 4,082,023 | |
| Indels | 925,563 | |
| Mixed SNP/indel sites | 4,618 | |
| Multiallelic sites | 66,860 | Must be split before annotation |
| Ts/Tv | 1.94 | Normal (includes filtered records and indels) |
| Heterozygous calls | 3,113,937 | |
| Homozygous alternate calls | 1,898,267 | |
| Het / hom-alt ratio | 1.64 | Within the outbred range |
| Median depth (sampled every 200th record) | **42×** | Ample |
| Mean depth | 44.2× | |
| Depth 5th / 95th percentile | 16× / 62× | |
| Fraction of called sites with DP < 10 | 2.4% | |
| Fraction of called sites with DP < 20 | 7.2% | |

42× median WGS is comfortably deep. It supports mosaic detection to roughly
5-10% variant allele fraction with adequate read support, which matters because
MVA is by definition a mosaicism disorder.

### FORMAT fields available

`GT`, `AD`, `DP`, `GQ`, `PL`, plus **`PGT` and `PID`**.

`PGT`/`PID` are HaplotypeCaller physical phasing tags. In a singleton with no
parents these are the only native phasing signal available, and they will
partially resolve *cis* versus *trans* for candidate compound heterozygotes that
fall within the same assembly region. This is load-bearing for branch C and is
recorded here so the phasing arm is designed around it rather than rediscovering
it later.

---

## 4. Clinical phenotype

Source: `data/Challenge_Clinical_Phenotype_1.docx`, extracted to text by
`scripts/02_extract_clinical.sh`. Eight features, each already HPO-coded by the
challenge organisers. No free-text extraction or LLM normalisation is required.

| HPO ID | Term | Note |
|---|---|---|
| `HP:0002859` | Rhabdomyosarcoma | The primary oncological event that triggered urgent investigation |
| `HP:0000121` | Nephrocalcinosis | Present since birth |
| `HP:0004322` | Short stature | Below expectation for age and family background |
| `HP:0001508` | Failure to thrive | Documented repeatedly |
| `HP:0003202` | Skeletal muscle atrophy | Reduced muscle bulk |
| `HP:0001622` | Premature birth | 32 weeks gestation |
| `HP:0001518` | Small for gestational age | Approximately 1 kg at birth, consistent with IUGR |
| `HP:0200067` | Recurrent spontaneous abortion | Parental history, multiple losses before the proband |

The organisers flag the parental reproductive loss as phenotypic input rather
than background, and ask that the constellation be read as a whole.

### Sex: male, determined from the callset

The phenotype document does not state the proband's sex. It is determinable
directly: chrX heterozygosity is 0.062 against 0.620 on the autosomes, and chrY
carries 9,732 PASS variants. The proband is male.

This raises the X-linked model from the low prior the plan assigns it, because
for a male an X-linked recessive cause needs a single hit rather than a compound
heterozygote with a cryptic second allele. See
`results/summaries/proband_sex_and_x_linked.md`.

### What the phenotype document does **not** contain

These absences change the plan and are recorded deliberately.

- **No karyotype and no aneuploidy percentage.** There is no cytogenetic result
  of any kind. This removes the ground truth that plan section 6.4 assumed for
  correlating candidate variant allele fraction against per-tissue aneuploidy.
- **No candidate variant.** No gene, no heterozygote, no prior genetic finding
  is named. The project is therefore **not** branch E.
- **No tissue other than the sequenced sample.** Single specimen, so no
  cross-tissue VAF comparison is possible.
- **No microcephaly, seizures, developmental delay or intellectual
  disability**, which are prominent in the canonical MVA description. Their
  absence should be carried into phenotype-driven ranking honestly rather than
  assumed away.
- No Wilms tumour, no leukaemia.

### Phenotype-driven prior

Rhabdomyosarcoma plus IUGR plus growth failure plus nephrocalcinosis plus
parental recurrent miscarriage is a coherent chromosomal-instability and
cancer-predisposition cluster. Rhabdomyosarcoma in particular is the tumour most
characteristically reported in `BUB1B`-related MVA. However, the differential
must stay open: this constellation overlaps other cancer predisposition and
DNA-repair syndromes, and Arm E exists to test exactly that.

---

## 5. Known MVA gene panel: aggregate tally

Coordinates from Ensembl REST (`/lookup/symbol`, GRCh38), stored with
provenance in `config/gene_panels/mva_known.tsv`. Intervals carry 5 kb flanks.

| Gene | Region (no-chr) | Span | Variants | PASS | Filtered | Het | Hom alt |
|---|---|---|---|---|---|---|---|
| `BUB1B` | 15:40155984-40226137 | 70.2 kb | 17 | 16 | 1 | 10 | 7 |
| `CEP57` | 11:95784965-95842070 | 57.1 kb | 7 | 7 | 0 | 4 | 3 |
| `TRIP13` | 5:887849-924357 | 36.5 kb | 3 | 3 | 0 | 3 | 0 |
| `BUB1` | 2:110630468-110683098 | 52.6 kb | 15 | 15 | 0 | 9 | 6 |
| `BUB3` | 10:123149395-123175467 | 26.1 kb | 32 | 32 | 0 | 32 | 0 |
| `CENATAC` | 11:118993051-119020811 | 27.8 kb | 20 | 20 | 0 | 3 | 17 |
| `CEP192` | 18:12986283-13130053 | 143.8 kb | 189 | 189 | 0 | 74 | 115 |
| `SMC5` | 9:70253270-70359874 | 106.6 kb | 98 | 97 | 1 | 43 | 55 |
| `CEP57L1` | 6:109090105-109179418 | 89.3 kb | 21 | 21 | 0 | 15 | 6 |

Counts were cross-validated against an index-independent streaming pass
(`scripts/04_verify_panel_tally.py`); the shipped `.tbi` agrees exactly on every
interval, so region-scoped queries are trustworthy.

**Correction, 31 August 2026.** The `BUB3` row originally reported 199 variants
over a 168.8 kb span. Those coordinates came from the live Ensembl REST API,
which returned a 158,779 bp span for `BUB3` against 16,072 bp in the pinned
Ensembl 115 GTF, a 9.9-fold overstatement covering neighbouring sequence. The
row above is corrected and the panel is now built from the pinned GTF. The other
eight genes agreed within 10% and are unaffected. The conclusion drawn from the
original figure, that `BUB3` shows no coverage anomaly, survives the correction:
at the correct span its density is the highest of the MVA genes.

### Two candidate leads raised and both closed at Phase 0

Recording these because closing a lead cheaply is worth as much as opening one,
and both would otherwise have been re-raised in week 2.

**Lead 1: apparent coverage collapse over `TRIP13`. Closed, artefact of small
sample.** `TRIP13` initially showed median depth 7× against a genome median of
42×. That figure came from only three called sites. The 10 kb profile in
`scripts/05_panel_depth_profile.py` shows the depth ratio of gene body to
300 kb flank is **1.03** for `TRIP13`, and between 0.86 and 1.07 for every other
panel gene. There is no coverage collapse over any known MVA gene.

**Lead 2: variant-free gaps inside `BUB1B` (10 kb), `CEP57` (20 kb) and
`TRIP13` (20 kb). Closed, within background rate.** Calibration on chromosomes
5, 11 and 15 (`scripts/06_gap_background_rate.py`) gives 2.5-3.3% of evaluable
10 kb windows empty, 37 to 85 variant-free runs of 20 kb or longer per
chromosome, and longest runs of 70 to 260 kb. Gaps of this size are ordinary
genomic background, most plausibly short runs of homozygosity. They are not
evidence of deletion.

Both closures are provisional in one specific respect: a homozygous deletion
would also produce zero calls, and the VCF cannot distinguish "no variant" from
"no coverage" because it is a variants-only file, not a gVCF. **Definitive
exclusion requires read depth from an alignment**, which is Arm C work. The
gaps are logged in `results/recon/panel_depth_profile.json` as regions to check
first once the BAM exists, at a low prior.

### Runs of homozygosity: no evidence of consanguinity

A 1 Mb windowed heterozygosity scan (`scripts/03_roh_proxy.py`, 2,777 autosomal
windows) gives a median heterozygous fraction of 0.633 and a median 1,626 PASS
variants per Mb. Only 33 windows are ROH-like, and the longest contiguous run is
2 Mb. Genuine consanguinity produces multiple runs above 5-10 Mb. There is no
such signal.

**Consequence for the recessive model:** homozygosity by descent is not the
expected mechanism. The compound heterozygote remains the leading hypothesis,
which keeps plan section 0.2 class 1, the cryptic second allele, at the top of
the prior list. All low-density windows found by the scan are centromeric or
acrocentric (chr7:59-60 Mb, chr21:6-7 Mb, chr17:24-25 Mb and similar) and are
expected mappability holes.

---

## 6. Consequences for the plan

| Question from plan section 2.2 | Answer | Consequence |
|---|---|---|
| Trio, quad or singleton? | **Singleton** | De novo calling and segregation filtering are impossible. Phasing must be read-backed (`PGT`/`PID`) and population-based. |
| Reference build? | **GRCh38, no-chr contig naming** | Written to `config/config.yaml`. Every annotation resource must be GRCh38 and every tool needs the naming convention handled explicitly. |
| BAM or CRAM present? | **No, but FASTQ is** | Arms C and D are not dead, they are deferred behind an alignment step that this laptop cannot run. See section 7. |
| RNA-seq present? | **No** | FRASER2 and OUTRIDER are out. Splicing predictions stay predictions. Track 2 signature reversal must use a published proxy signature, clearly labelled. |
| How many tissues? | **One** | The VAF against aneuploidy-percentage correlation in plan section 6.4 is not possible. |
| Caller and pre-filtering? | **Sentieon 202308.02, GATK hard filters, non-destructive** | 271,414 filtered records are recoverable. `MQ40` regions are a deliberate search target. |
| Coverage over MVA genes? | **Normal, 0.86-1.07× of local flank** | No homozygous-deletion lead from depth. Provisional pending alignment. |
| Clinical data format? | **`.docx`, already HPO-coded, 8 terms** | Exomiser and LIRICAL can be driven directly. No LLM extraction needed. |
| Karyotype / aneuploidy %? | **Absent** | Plan section 6.4 loses its ground truth. Report as a data limitation. |
| Clinical data names a candidate het? | **No** | Not branch E. The search is open, not a second-allele hunt. |

---

## 7. Compute environment: a material constraint

The plan assumes an RTX 6000 with 48 GB VRAM. The machine this recon ran on is
**not that machine**:

| Property | Value |
|---|---|
| Host | Apple M2 MacBook Pro, `arm64` |
| RAM | **8 GB** |
| GPU | None usable (no CUDA) |
| Free space, data volume | 181 GB of 934 GB |
| Free space, internal disk | 5.3 GB of 228 GB |
| Docker | Installed, **daemon not running** |
| Bioinformatics tooling | Only `bcftools` and `samtools`, installed during this recon |

Aligning 84.7 GB of gzipped FASTQ needs roughly 250-400 GB of scratch plus the
reference and index. With 181 GB free and 8 GB of RAM this is not feasible here.
`bwa-mem2` alone wants about 30 GB of RAM for the human index.

**VCF-scoped work runs comfortably on this laptop. Alignment-dependent work must
move to the RTX 6000 host.** The 79 GB dataset is already on a portable external
SSD, so the transfer path exists.

---

## 8. Data governance as executed

- `.gitignore` excluding `data/` and every genomic extension was written
  **before** the first commit. The repository had no commits at that point.
- A pre-commit hook (`.githooks/pre-commit`, `core.hooksPath` set) hard-fails on
  anything under `data/`, on genomic extensions anywhere, and on files above
  5 MB.
- No genomic variant, coordinate, allele or genotype from the proband has been
  read into an LLM context. All figures in this document are aggregates
  produced by scripts.
- **The clinical phenotype document was read.** It is patient data though not
  genomic data. See `ETHICS.md` section 3 for the reasoning and for how to
  reverse the decision.
- The phenotype document states that the family have published their own blog
  posts. **No attempt has been made or should be made to locate them.**
  Enriching the phenotype from family-authored public material would
  re-identify the family and exceeds the consent as described.
