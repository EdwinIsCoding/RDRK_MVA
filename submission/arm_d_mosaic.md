# Arm D: mosaic and low variant allele fraction. Reported negative.

> **Arm C (structural variants) was NOT completed.** Delly was launched twice.
> The first run failed instantly because Delly 2.6 renamed its `call`
> subcommand to `sr`. The second run started correctly but the compute booking
> ended and `/scratch0` was wiped before the results could be retrieved. No
> genome-wide SV call set exists. What we do have for Arm C is read-level
> verification and panel coverage (`arm_c_readlevel_verification.md`), which is
> not the same thing and should not be described as if it were.

Plan section 6.4. MVA is by definition a mosaicism disorder, so allele fractions
that a germline caller discards are exactly the signal worth re-examining.

Method: `bcftools mpileup` re-genotyping across the nine known MVA genes against
our own alignment, with **no variant allele fraction filtering at all**, then
inspection of raw allele fractions.

## Result

| Threshold | Sites |
|---|---:|
| VAF between 0.03 and 0.30 at depth >= 20 | 1,463 |
| of those, supported by **>= 5** alternate reads | **6** |
| of those, supported by **>= 10** alternate reads | **2** |

**The headline number is misleading and should not be quoted alone.** Of 1,463
sites in the low-VAF band, 1,457 rest on between one and four alternate reads.
At 45x depth that is the expected rate of sequencing and alignment error, not
mosaicism. Six sites have five or more supporting reads and two have ten or
more, and those are the only candidates worth any further attention.

## Interpretation

**No credible mosaicism in the known MVA genes.** This is consistent with the
answer: both causal alleles sit at VAF 0.553 and 0.444, which is germline
heterozygous, not mosaic. The proband's chromosomal instability is a downstream
consequence of biallelic BUB1B loss of function, not of a mosaic causal variant.

That distinction is worth stating explicitly because the disease name invites the
opposite assumption. Mosaic variegated aneuploidy is mosaic in its *chromosome
counts*, cell to cell. The causal genotype is constitutional.

## Limits

- Scope is the nine known MVA genes, not the whole panel or genome.
- `bcftools mpileup` with a diploid model is not a dedicated somatic caller.
  **Addressed 3 September 2026.** Mutect2 4.7.0.0 tumour-only was run over the
  same nine genes on the same BAM: it reduces the 0.03-0.30 band from 1,463
  sites to 64 and leaves zero passing with 5 or more supporting reads, against 6
  here. Both instruments agree. See `results/summaries/arm_d_mutect2.md`.
  DeepSomatic was not run.
- The plan's stronger design, correlating candidate VAF against reported
  aneuploidy percentage per tissue, remains impossible: one tissue, and no
  karyotype or aneuploidy percentage exists in the clinical data.
