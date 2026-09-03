# STOP #2 status: partial

Plan section 5.4 requires a recall table before the pipeline touches the
proband, and specifies the gate: *"If masked recall@20 for splice-class
positives is below ~50%, the splicing arm needs rework before it touches the
real data."*

**That gate cannot yet be evaluated, and the reason is not a scheduling
delay.** This document says what was built, what was measured, and what
specifically blocks the number the gate asks for.

Written 31 August 2026. **Updated 1 September 2026** with the outcome; see the
closing section.

---

## 1. The gate cannot be evaluated as written

Plan section 5.4 asks for masked recall@20 **for splice-class positives**
against the MVA benchmark. Those positives do not exist.

`benchmarks/published_mva_variants.tsv` holds every ClinVar variant across the
nine MVA and candidate genes: 1,551 rows, 1,541 with linked PMIDs, each
traceable to a VCV accession. Restricted to tier 1, meaning confidently
pathogenic and therefore usable as ground truth, 108 remain:

| Variant class | Tier 1 count |
|---|---:|
| nonsense | 49 |
| frameshift or indel | 35 |
| canonical splice site | 15 |
| deletion | 5 |
| missense | 3 |
| **deep intronic** | **0** |
| **near splice region** | **0** |
| **synonymous** | **0** |
| **UTR or promoter** | **0** |

The four classes at zero are precisely the cryptic second allele that
`RECON.md` ranks as hypothesis class 1, the leading hypothesis for this
proband.

This is not an artefact of how the harvest was run. It is a property of the
evidence base: cryptic alleles are hard to find, so they are under-ascertained,
and when found they are usually deposited as uncertain rather than pathogenic.
The scarcity is the same fact that makes the proband's case unsolved.

**The consequence matters.** A recall figure computed on the 108 available
positives would be dominated by nonsense and frameshift variants, which any
competent pipeline recovers. It would look reassuring and would say nothing
about whether the pipeline can find the kind of variant we are actually looking
for. Reporting it without this qualification would be the most misleading
number in the submission.

A regression test (`tests/test_positive_controls.py::test_documents_the_missing_variant_classes`)
fails if tier-1 rows in these classes ever appear, so the limitation cannot
silently go stale.

## 2. What was built to cover the gap

`benchmarks/splice_mechanism_controls.tsv`: 130 confidently pathogenic deep
intronic and near-splice variants across 100 genes and many diseases, harvested
from ClinVar with a seeded, reproducible sample and with MVA-panel genes
excluded so the set stays independent.

| Offset band from the nearest exon boundary | n | Why the band matters |
|---|---:|---|
| 11 to 50 bp | 44 | A default ±50 bp window still reaches these |
| 51 to 100 bp | 9 | Needs a widened window |
| 101 to 500 bp | 22 | The regime the plan's ±500 bp setting exists for |
| beyond 500 bp | 55 | Expected misses; they bound achievable recall honestly |

**This set tests the machinery, not the disease.** A good score means the
splicing arm can find a deep intronic pathogenic variant. It does not mean the
splicing arm can find the MVA allele. That distinction is written into the
harvest script's docstring and belongs in the submission in those words.

The 31 variants in the 51 to 500 bp band are the ones that specifically justify
running SpliceAI at ±500 bp rather than at its default.

## 3. What has been measured

### Region annotation is exact for SNVs

Splice distance drives the entire variant-class stratification, so it was
validated against an independent source of truth: ClinVar's own HGVS intron
offsets.

| | Concordance |
|---|---|
| SNVs | **268/268, 100.0%** |
| indels | 8/33, 24.2% |

The indel figure is a comparison artefact, not an error. HGVS 3'-shifts an
indel within a repeat while VCF left-aligns it, so the conventions legitimately
disagree by the repeat length. Recorded, with the consequence that indel splice
distances near a class boundary should not be trusted.

Reaching 100% required a fix. Pooling exon boundaries across all annotated
transcripts gave 89.4%, because a minor isoform's boundary could sit closer than
the canonical one's, promoting genuine deep intronic variants into the
canonical-splice class. Splice distance is now measured against MANE Select,
falling back to Ensembl canonical; all 736 panel genes resolve one.

### Constraint is the wrong prior for this disease

Measured across the known MVA genes with gnomAD v4.1:

| Gene | LOEUF | pLI | | Gene | LOEUF | pLI |
|---|---|---|---|---|---|---|
| BUB1B | 0.707 | 0.000 | | BUB3 | 0.770 | 0.001 |
| CEP57 | 0.740 | 0.000 | | CEP192 | 0.610 | 0.000 |
| TRIP13 | 0.592 | 0.192 | | SMC5 | 0.716 | 0.000 |
| BUB1 | 0.727 | 0.000 | | CENATAC | 1.227 | 0.000 |

Not one is constrained. This is what recessive biology predicts: LOEUF measures
selection against heterozygous loss of function, and MVA carriers are healthy.
Weighted as it would be in a dominant-disease pipeline, constraint would have
actively deprioritised every known answer. It now contributes a small bonus and
never a penalty (`src/mva/track1/scoring.py`).

Separately: **gnomAD v4.1 constraint covers autosomes only.** X-linked panel
genes (SMC1A, STAG2, ATRX, HDAC8 and others) can never receive a value and must
not be penalised for its absence.

## 4. What blocks the real number

The pipeline runs end to end, but three annotators raise rather than returning
neutral values, because a stub that quietly returned "no evidence" would let an
incomplete run masquerade as a complete one:

| Annotator | Needs | Consequence while absent |
|---|---|---|
| `VepAnnotator` | VEP cache, ~25 GB | No protein consequence. Missense, nonsense and frameshift classes unavailable. |
| `GnomadAnnotator` | gnomAD sites VCF, ~60 GB | No allele frequency. Absence of a frequency is not evidence of rarity. |
| `SpliceAiAnnotator` | SpliceAI and Pangolin weights, GPU | **Arm B cannot run at all.** This is the highest-prior arm. |

None can be installed on the recon host: 8 GB RAM, no CUDA, and the internal
disk filled once already during this work. All three are pinned in
`environment.yml` and install on the GPU host.

**Until SpliceAI runs, any recall figure measures gene-level triage, not variant
prioritisation.** The tests state this in their assertions rather than in a
footnote.

## 5. Recommendation

Do not treat this as STOP #2 cleared. The honest position is:

1. The harness is built, tested and validated where validation was possible.
2. The gate's specific question cannot be answered by the MVA benchmark, and
   that is a finding about the evidence base rather than a gap in the work.
3. The substitute measurement is ready and will run the moment SpliceAI does.

When the GPU host is available, the order is: install the toolchain, run recall
against both `published_mva_variants.tsv` and `splice_mechanism_controls.tsv`,
masked and unmasked, then report both tables with the qualification in section 2
attached to the second.

If the GPU host slips past roughly 15 September, plan section 10's buffer
discipline applies: cut Arms C and D, and go deep on Arm B.


---

## 10. Outcome, added 1 September 2026

The gate could not be evaluated as specified, and the case was solved anyway.
Both facts belong here.

**The answer is a `BUB1B` compound heterozygote**: `chr15:40209701 T>G`
(p.Leu737Ter, ClinVar VCV000533901.9, Pathogenic/Likely pathogenic for MVA1)
with `chr15:40220612 T>G` (p.Asn1002Lys, ultra-rare rather than novel: a
single allele in 1,461,878 in gnomAD v4.1 exomes, absent from genomes only).

### What this says about the gate

Section 1 argued that a recall figure computed on the 108 available positives
would be dominated by nonsense and frameshift variants and would say nothing
about finding a cryptic allele. That reasoning was sound but the conclusion drawn
from it was too pessimistic in one direction and too optimistic in another.

**Too pessimistic:** the benchmark's bias towards nonsense and frameshift
positives turned out to match the actual answer. One of the two causal alleles
is a nonsense variant, exactly the class the benchmark over-represents. A recall
figure on that set would have been more informative about this case than
predicted.

**Too optimistic:** the missing classes were assumed to matter because the
cryptic-allele hypothesis was ranked highest. It was wrong. The second allele is
an ordinary novel missense, and the benchmark contains only three missense
positives. So the benchmark was unrepresentative, but not in the direction the
analysis expected.

### What the blocked annotators cost

Section 4 listed VEP, gnomAD and SpliceAI as blocking a real recall figure. Two
of the three turned out not to need a GPU at all, only a route around a large
download: gnomAD frequencies came from remote tabix range requests against the
public per-chromosome sites files, and consequence annotation from VEP 116.1 on
a cluster node. **gnomAD alone reduced 415 scored variants in the known MVA genes
to 12, both causal alleles among them.** Consequence annotation was needed to
rank those 12, not to find them.

The honest summary is that the population-frequency filter did the work, the
consequence annotation did the ranking, and the splicing arm, which the plan
ranked highest, contributed a reported negative.
