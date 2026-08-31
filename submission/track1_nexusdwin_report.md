# Track 1 report: BUB1B compound heterozygote in PROBAND01

**Submitter:** NexusDwin
**Repository:** see the GitHub URL supplied with this submission
**Date:** 31 August 2026
**Licence:** CC BY 4.0

---

## 1. The call

| | Allele 1 | Allele 2 |
|---|---|---|
| Position (GRCh38) | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
| HGVS | `NM_001211.6:c.2210T>G` | `NM_001211.6:c.3006T>G` |
| Protein | `p.Leu737Ter` | `p.Asn1002Lys` |
| Consequence | nonsense (stop gained) | missense |
| dbSNP | rs759242053 | none, novel |
| ClinVar | **VCV000533901.9**, Pathogenic/Likely pathogenic, multiple submitters, no conflicts, listed against *Mosaic variegated aneuploidy syndrome 1* | no record |
| gnomAD v4.1 popmax | 7.87 × 10⁻⁵ | **absent** |
| In silico | n/a for a nonsense | SIFT deleterious (0.01), PolyPhen-2 probably damaging (0.997) |
| Genotype | 0/1, DP 46, AD 21/25, GQ 99, PASS | 0/1, DP 28, AD 15/13, GQ 99, PASS |

Both fall on the MANE Select transcript `ENST00000287598.11`.

**Interpretation.** Biallelic loss of function in `BUB1B` causes mosaic variegated
aneuploidy syndrome 1 (OMIM 257300), an autosomal recessive chromosomal
instability disorder. A premature termination codon at p.Leu737, already
classified pathogenic for this exact condition, in trans with a novel missense
predicted damaging by two orthogonal methods, is the canonical MVA1 allelic
architecture: one clearly disruptive allele plus one hypomorphic allele.
Complete BUB1B nullity is not compatible with life, so a residual-function
second allele is expected rather than surprising.

The phenotype fits. Rhabdomyosarcoma (HP:0002859) is the tumour most
characteristically reported in BUB1B-related MVA; intrauterine growth
restriction, prematurity, failure to thrive and short stature are core features.

**ACMG/AMP criteria applied.** Allele 1: PVS1 (nonsense in a gene where loss of
function is the established mechanism), PM2_supporting (popmax 7.9 × 10⁻⁵),
PP5 (ClinVar pathogenic, multiple submitters). Allele 2: PM2 (absent from
gnomAD), PP3 (two concordant in silico predictors), PM3_supporting (in trans
with a pathogenic allele, *inferred rather than demonstrated*, see section 5).

---

## 2. Approach

The pipeline is a transparent additive evidence framework, not a learned ranker.
With n = 1 there is nothing to train on, and a black-box model over a single
proband would fit noise while sounding authoritative. Every weight is fixed
a priori, documented in `src/mva/track1/scoring.py`, and every score decomposes
into evidence items each carrying a resolvable identifier.

**Stage 1, characterisation.** The callset is GRCh38 with Ensembl contig naming
(`15`, not `chr15`), single sample, 5,012,204 records, median depth 42×, called
by Sentieon 202308.02 with GATK hard filters applied non-destructively so 271,414
filtered records remain recoverable.

**Stage 2, gene priors.** Three panels, each with per-gene provenance:
- 9 known MVA and candidate genes, coordinates from the pinned Ensembl 115 GTF.
- 408 mitotic genes from GO, Reactome and STRING, weighted by GO term
  specificity and normalised to approved HGNC symbols.
- 4,774 curated disease genes from ClinGen, Genomics England PanelApp (eleven
  panels matched to the proband's HPO terms) and EBI gene2phenotype.

**Stage 3, population filtering.** gnomAD v4.1 allele frequencies were obtained
by remote tabix range requests against the public per-chromosome sites files,
transmitting only gene intervals. Across the known MVA genes this reduced 415
scored variants to **12**, of which 403 were common polymorphisms.

**Stage 4, consequence.** Ensembl VEP 116.1 with the merged cache over all
5,012,204 records.

**Stage 5, phenotype.** An information-content weighted match between the
proband's eight HPO terms and the HPO gene-to-phenotype corpus.

---

## 3. What actually found the variants

Both causal alleles were present in the shortlist of **3 rare heterozygous BUB1B
variants** produced at stage 3, using only population frequency and gene-panel
membership. Consequence annotation was needed to *rank* them, not to *find*
them. That artefact is timestamped in the repository (commit `849bf98`).

The phenotype prior independently ranked `BUB1B` **18th of 2,503 genes** on the
HPO corpus alone, top 0.8%, without using any variant data.

---

## 4. Negative results

Reported because excluding an explanation is a result.

| Question | Method | Outcome |
|---|---|---|
| Consanguinity or homozygosity by descent? | `bcftools roh`, plus an independent 1 Mb windowed heterozygosity scan | **No.** Longest homozygous segment 1.10 Mb; nothing above 1.5 Mb against the multiple >10 Mb tracts a first-cousin union produces. Two methods, one answer. |
| A cryptic splice allele, the a priori leading hypothesis? | SpliceAI at `-D 500` over the rare variants in the known MVA genes | **No.** Nothing reached even the permissive 0.2 threshold; BUB1B maximum delta 0.030. |
| Coverage gap or deletion over an MVA gene? | 10 kb depth and density profiling | **No.** Gene-body to flank depth ratio 0.86–1.07 across all nine genes. |
| Mitochondrial cause? | Contig `M` analysis | **No**, but weakly: 13 homoplasmic variants at median 4,177×, and the diploid germline caller used cannot detect low-level heteroplasmy. |

The splicing negative is worth dwelling on. The project plan reasoned that a
cryptic second allele was the highest-prior hypothesis, since a standard coding
pipeline would already have solved the case. **That reasoning was wrong.** The
second allele is an ordinary novel missense that any competent pipeline reaches
once it has consequence annotation. The searching was not wasted, because the
negative bounds the answer, but the prior was mistaken and is reported as such.

---

## 5. Limitations, stated plainly

**Trans configuration is inferred, not demonstrated.** The two alleles lie
**10,911 bp apart**, well beyond a read pair. HaplotypeCaller's `PGT`/`PID`
physical phasing reported no phasing group anywhere in `BUB1B`, which is exactly
what that separation predicts. No amount of short-read depth resolves this;
confirmation requires parental testing or long reads. We infer *trans* from the
recessive mechanism and from neither allele appearing on a shared haplotype in
population data, and we do not claim to have shown it.

**No RNA-seq exists for this proband**, so every splicing result above is a
prediction and never an observation of an aberrant junction.

**Singleton.** No de novo calling, no segregation filtering, no parental phasing.

**Missense pathogenicity rests on SIFT and PolyPhen-2 only.** AlphaMissense, CADD
and REVEL were not available in this run.

**We were not blind.** The Track 1 leaderboard was saturated with 61 perfect
scores before we submitted, and filenames on it identify the gene. Our shortlist
containing both alleles predates that (section 3, commit `849bf98`), but the
consequence annotation that ranked them was examined afterwards, and we knew the
gene when we looked. Reporting this is the honest description of how the result
was reached.

---

## 6. Validation

Method validation was treated as a first-class deliverable, not a footnote.

- **A benchmark of 1,551 ClinVar variants** across nine MVA genes, 1,541 with
  linked PMIDs, tiered by classification confidence. Its most useful property is
  a limitation it exposes: of the 108 confidently pathogenic variants, **none is
  deep intronic, near-splice, synonymous or UTR**. The benchmark cannot test the
  cryptic-allele hypothesis, so a separate mechanism control set of 130
  confidently pathogenic deep-intronic variants from other genes was built.
- **Tools are validated before their negatives are believed.** SpliceAI initially
  returned 0.000 for everything, including eight known pathogenic canonical
  splice-site variants in `BUB1B`. The cause was a NumPy 2 compatibility shim of
  ours that dropped a base-encoding step. Corrected, those controls score 9/9
  above 0.5 with a maximum delta of 1.000. The runner now **refuses to report a
  negative** if its positive controls fail.
- **Calibration that changed the design.** Every known MVA gene is unconstrained
  (BUB1B LOEUF 0.707, pLI 0.000; CENATAC 1.227). Weighted as a dominant-disease
  pipeline would, constraint would have actively deprioritised the correct answer.
  It contributes a small bonus and never a penalty. On chrX, where the male
  proband is hemizygous and selection is direct, it is weighted higher.
- **Region annotation checked against an independent source**: splice distances
  agree with ClinVar HGVS intron offsets for **268/268 SNVs**.
- **140 automated tests** covering the evidence schema, scoring, annotators and
  submission format.

---

## 7. Secondary finding

**`LZTR1` chr22:20996720 C>G, `NM_006767.4:c.2244C>G`, p.Tyr748Ter**, heterozygous
nonsense, rs1682503990, gnomAD popmax 1.4 × 10⁻⁶, DP 48, GQ 99, PASS.
**ClinVar VCV001409252.7, Pathogenic/Likely pathogenic**, multiple submitters, no
conflicts, against Noonan syndrome 10 (OMIM 616564), Noonan syndrome 2
(OMIM 605275) and LZTR1-related schwannomatosis (OMIM 615670).

Reported because it is an established pathogenic variant with tumour
surveillance implications warranting clinical review, and because
rhabdomyosarcoma, short stature and failure to thrive overlap the RASopathy
phenotype, so a dual diagnosis or modifying contribution cannot be excluded on
this data. It is not proposed as the primary cause; MVA1 explains the
chromosomal instability presentation and the BUB1B pair is the better fit.

Two homozygous loss-of-function calls absent from gnomAD (`PEX5`, `CTU2`) were
**deliberately excluded**. A true homozygous PEX5 knockout causes Zellweger
spectrum disease, which this child plainly does not have, so these are almost
certainly mis-calls in repetitive sequence. Read-level verification is pending.

---

## 8. Data handling and AI assistant disclosure

Required by the organisers' update of 28 August 2026.

> **Anthropic Claude (Claude Code), Max subscription, consumer terms. The "help
> improve Claude" training setting was enabled until 31 August 2026 and disabled
> thereafter.**

We disclose a compliance gap rather than presenting a clean summary. For part of
this work the account permitted training on conversation content, which does not
meet the organisers' condition of no training on inputs or outputs. The setting
has been disabled, account data has been purged, and Sage Bionetworks has been
notified.

**What the model was and was not exposed to.** The project's first rule, written
before any analysis and enforced by a pre-commit hook, was that patient genomic
data never enters an LLM context: scripts read `data/`, the model reads only
aggregate summaries. Measured against the organisers' own delete list, **no VCF
or BAM content, no variant table with genotypes, and no prompt containing a
pasted variant block** was ever sent. What was sent falls on their keep list:
HPO terms, gene-level and aggregate statistics, and code. Two patient-derived
items are named rather than glossed: the clinical phenotype document, and five
runs-of-homozygosity interval coordinates.

That architecture is what bounds the consequences of the gap. Had the work been
done the obvious way, by pasting variant tables into prompts, the exposed
material would have been a child's variant-level genome instead of eight
phenotype terms.

**Deletion.** All challenge data will be deleted at the conclusion of the
hackathon and `MVAHackathon2026@synapse.org` notified, per the data access terms.

---

## 9. Reproducibility

Every result is regenerable. The repository contains pinned environments
(`environment.yml`, `pyproject.toml`, `Dockerfile`), a Snakemake workflow, and
`make reproduce`. Reference resources are versioned: Ensembl 115 GTF, VEP 116.1
merged cache, gnomAD v4.1 (v2.1.1 for chrX constraint, which v4.1 omits), HPO
and ClinVar snapshots dated in `PROVENANCE.md`. SHA256 checksums of all eleven
input files are recorded. Randomness is seeded at 20261024.

`data/` and all derived results are excluded from version control by
`.gitignore` and by a pre-commit hook that rejects genomic file extensions,
anything under `data/` or `results/`, the proband sample identifier in a
filename, and files above 5 MB.
