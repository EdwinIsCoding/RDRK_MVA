# RECON.md — hypothesis class triage

Companion to `DATA_CARD.md`. States, for each of the six hypothesis classes in
`MVA_HACKATHON_PLAN.md` section 0.2, whether it is testable with the data
actually shipped, and drops the rest with a reason.

Written 31 August 2026 at STOP #1.

---

> **Resolved 1 September 2026.** The alignment was built from FASTQ in 4h10m,
> so the deferral below no longer applies and every arm is unblocked. The case
> is solved: a `BUB1B` compound heterozygote, `chr15:40209701 T>G`
> (p.Leu737Ter) with `chr15:40220612 T>G` (p.Asn1002Lys). The triage below is
> left as written, because how the reasoning stood before the answer was known
> is the honest record, and because one of its central calls was wrong.
>
> **Class 1, the cryptic second allele, was ranked highest and was wrong.** The
> second allele is an ordinary novel missense, not a cryptic splice variant.
> SpliceAI at plus or minus 500 bp found nothing above even the permissive 0.2
> threshold. The prior was reasonable given what was known, and it was mistaken.

## Confirmed branch: **C**, with one important amendment

Plan section 2.4 defines branch C as *singleton plus BAMs*. What we have is
**singleton plus FASTQ**, no alignments.

This is better than branch D (VCF only) and worse than branch C as written:

- Better than D, because 84.7 GB of raw paired reads across four lanes means
  every alignment-dependent method is technically available. Nothing is lost
  permanently.
- Worse than C, because those methods are gated behind a full WGS alignment that
  costs real time and cannot run on the recon machine (`DATA_CARD.md` section 7).

So: **branch C (deferred)**. Arms C and D are live but blocked on an alignment
step that must be scheduled explicitly and early, because everything downstream
of it queues behind it.

The three structural facts that set this branch:

1. **Singleton.** One sample, `WGS_EX2312012`. No parents. De novo detection and
   segregation filtering are both impossible.
2. **No RNA-seq.** No functional confirmation of any splicing prediction, and no
   patient-derived signature for Track 2.
3. **No prior genetic finding in the clinical document.** Not branch E. The
   search is open rather than a hunt for a known second allele.

---

## Hypothesis class triage

### Class 1 — Cryptic second allele in a known gene · **TESTABLE, highest prior**

Deep intronic, branch point, polypyrimidine tract, UTR, promoter, or synonymous
variants creating a cryptic splice site.

**Why it survives Phase 0 as the leading hypothesis.** The heterozygosity scan
found no consanguinity (longest ROH-like run 2 Mb), so homozygosity by descent
is not the expected mechanism and a compound heterozygote is. Published `BUB1B`
MVA families repeatedly show one coding allele plus one cryptic allele, and a
cryptic allele is exactly what a diagnostic pipeline reports as negative. The
callset was produced by a standard clinical workflow, which is the situation in
which this class is missed.

**What we can do.** SpliceAI and Pangolin at `-D 500` over the introns, UTRs and
2 kb promoters of the extended mitotic panel. Branch point scoring. uORF
annotation. All of it runs from the VCF alone.

**The honest limitation, stated up front.** With no RNA-seq we cannot observe an
aberrant junction. Every result in this arm is a *prediction*, and the report
must say so in those words. The strongest evidence this arm can produce is a
converging prediction plus a proposed RT-PCR experiment, not a demonstration.

**Phasing.** For a compound heterozygote in a singleton we cannot phase by
inheritance. Two partial substitutes exist and both should be used: the
`PGT`/`PID` physical phasing tags already present in the VCF, which resolve
variants within the same assembly region, and population-based phasing against a
reference panel. Where neither resolves the phase, say so rather than assuming
*trans*.

---

### Class 2 — Structural variant · **TESTABLE ONLY AFTER ALIGNMENT**

Deletions, duplications, inversions, mobile element insertions, complex
rearrangements.

**Status: blocked, not dropped.** Manta, Delly, CNVnator, ExpansionHunter and
MELT all consume alignments. We have none. We do have the FASTQ they would be
built from.

**Why it matters more here than usual.** 177,522 records carry the `MQ40` filter,
concentrated in repeats and segmental duplications, which is where short-read SV
calling fails and where diagnostic labs under-call. The parental recurrent
miscarriage history is also, independently, the kind of signal that a balanced
rearrangement produces, although we have no parental samples to test that.

**What Phase 0 already settled.** Depth across all nine panel genes is 0.86 to
1.07 times the local flanking depth, so there is no large homozygous deletion
over a known MVA gene detectable from the callset. Three variant-free gaps were
found (`BUB1B` 10 kb, `CEP57` 20 kb, `TRIP13` 20 kb) and calibration showed 37 to
85 such runs per chromosome genome-wide, so they are background. They are logged
for a depth check once a BAM exists, at low prior.

---

### Class 3 — Novel gene in mitotic machinery · **TESTABLE**

Kinetochore, SAC, centrosome, cohesin, condensin, SMC5/6, APC/C.

**Runs from the VCF.** The eight HPO terms are already coded by the organisers,
so Exomiser, LIRICAL and AMELIE can be driven directly with no extraction step.

**One caution carried forward from the phenotype.** The proband does *not* have
microcephaly, seizures or developmental delay, which are prominent in canonical
MVA. Phenotype-driven tools weight those heavily for MVA genes. The absence
should be handled explicitly rather than silently, and it widens the
differential beyond MVA proper towards other cancer predisposition and
DNA-repair syndromes. Arm E should not be scoped to mitotic genes alone.

---

### Class 4 — Mosaic / low variant allele fraction · **PARTIALLY TESTABLE, THEN BLOCKED**

**What is possible now.** The callset was generated with `--call_conf 10
--emit_conf 10`, which is permissive and retains more low-allele-fraction
evidence than the GATK default of 30. `AD` and `DP` are present on every record.
Median depth is 42×. So candidate low-VAF calls that *survived into the VCF* can
be examined immediately, and 42× supports interrogating allele fractions down to
roughly 5-10% with adequate read support.

**What is blocked.** Re-genotyping with Mutect2 tumour-only or DeepSomatic, and
force-calling panel positions with `bcftools mpileup`, all need alignments.
Variants that the caller never emitted cannot be recovered from the VCF at all.

**What is permanently lost.** Plan section 6.4 proposes correlating candidate VAF
against reported aneuploidy percentage per tissue. **This is not possible.**
There is one tissue and the clinical document contains no karyotype and no
aneuploidy percentage. The plan identified this as "strong, distinctive evidence
few teams will have"; it is unavailable, and the submission should say so
plainly rather than substitute a weaker analysis dressed up as the same thing.

---

### Class 5 — Uniparental disomy, imprinting, repeat expansion · **SPLIT**

- **UPD and runs of homozygosity: already partly done.** The 1 Mb scan found no
  extensive ROH, longest run 2 Mb. A formal `bcftools roh` pass will confirm and
  produce a citable negative. Cheap. **Note that UPD detection is substantially
  weaker in a singleton**, since the usual signal is Mendelian inconsistency
  against parents; what remains is homozygosity-tract evidence, which here is
  absent.
- **Repeat expansions: blocked.** ExpansionHunter and ExpansionHunter Denovo
  need alignments.

---

### Class 6 — Non-coding regulatory · **TESTABLE, LOWEST YIELD**

Enhancer disruption affecting a mitotic gene.

Runs from the VCF, but interpretation is weak: no RNA-seq to show an expression
consequence, no patient-tissue chromatin data, and enhancer-to-gene assignment is
itself uncertain. Keep this arm scoped to promoters and annotated regulatory
elements of the panel genes, treat any hit as a hypothesis of low confidence,
and do not let it consume time that class 1 should have.

---

## Arm status summary

| Arm | Plan section | Status as of 1 September 2026 | Outcome |
|---|---|---|---|
| A — Baseline annotation and filtering | 6.1 | **Complete** | Found the answer. 5.0M records to 12 rare variants in the known MVA genes, both causal alleles among them |
| B — Splicing and cryptic second allele | 6.2 | **Complete** | Reported negative. Nothing above 0.2; tool validated 9/9 on positive controls first |
| C — Structural variants | 6.3 | **Alignment done**, SV calling not run | 61 GB BAM built; read-level verification and coverage complete |
| D — Mosaic / low VAF | 6.4 | **Alignment done**, re-genotyping not run | Both causal alleles have VAF near 0.5, so not mosaic |
| E — Novel gene, phenotype-driven | 6.5 | **Complete** | HPO prior ranked BUB1B 18th of 2,503 genes without variant data |
| F — Completeness checks | 6.6 | **ROH and mtDNA complete**, repeat expansions not run | Both reported negatives |

SV calling, mosaic re-genotyping and repeat expansion detection are no longer
blocked; they were simply not reached before the answer was established and the
booking window closed. They would strengthen the negative results section rather
than change the call.

---

## What changes in the plan, and why

1. **Add an alignment stage that the plan does not contain.** The plan assumed
   BAMs. Producing them from 84.7 GB of FASTQ is a scheduled task with a real
   cost, and Arms C, D and part of F queue behind it. It should start as early
   as possible on the RTX 6000 host because it gates a third of the analysis.

2. **The recon machine is not the plan's machine.** 8 GB of RAM, no CUDA, 181 GB
   free, Docker daemon down. VCF-scoped work is comfortable here; alignment is
   not possible here. This needs resolving before week 1 (`DATA_CARD.md`
   section 7).

3. **Delete the VAF against aneuploidy correlation from the deliverables.** No
   karyotype, no aneuploidy percentage, one tissue. Replace it with an explicit
   statement of the limitation.

4. **Track 2 section 7.3 falls back to proxy signatures.** No patient
   transcriptome, so LINCS/CMap connectivity must be driven from published
   MVA/CIN model signatures, labelled as a proxy in every table where it appears.

5. **Arm E should not be scoped to mitotic genes only.** The phenotype lacks
   microcephaly, seizures and developmental delay. Constraining the search to
   the MVA gene neighbourhood risks assuming the answer.

6. **Weight the `MQ40` fraction as a search target, not as noise.** 177,522
   records, 3.5% of the callset, sit in low-mapping-quality sequence. Arm A must
   not filter to `PASS` blindly.

7. **Track 2 gains a firmer anchor than expected.** Rhabdomyosarcoma is the
   presenting event and cancer predisposition is explicit in the phenotype. Plan
   section 7.1 already identifies chemoprevention and surveillance as the
   highest-value and least-contested output; the phenotype supports that
   directly.

---

## Immediately actionable, in order

1. Resolve the compute question. Confirm access to the RTX 6000 host and move
   the external SSD to it, or decide that alignment-dependent arms are out of
   scope and say so in the submission.
2. Start the FASTQ to BAM alignment as soon as that host is available. It gates
   Arms C, D and F.
3. Build the extended mitotic panel and the positive-control benchmark. Both are
   pure desk work, need no patient data, and are on the critical path to
   STOP #2.
4. Run Arm A and Arm B on the VCF. Both are unblocked today.
