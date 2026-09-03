# Read-level verification from our own alignment

FASTQ to BAM built independently of the supplied callset: bwa-mem2 2.3 against
GRCh38 (Ensembl 115 primary assembly), four lanes merged, duplicates marked,
61 GB final BAM. Completed in 4h10m on 14 threads.

This matters because it is **orthogonal evidence**. The candidate variants were
called by Sentieon on someone else's pipeline; these reads were aligned by us
from raw FASTQ. Agreement is independent confirmation, not a restatement.

## The two BUB1B alleles

| | chr15:40209701 T>G | chr15:40220612 T>G |
|---|---|---|
| | p.Leu737Ter | p.Asn1002Lys |
| Depth (MQ>=20, BQ>=20) | 47 | 29 |
| Ref / alt reads | 21 / 26 | 16 / 13 |
| **VAF** | **0.553** | **0.448** |
| Strand balance of alt reads | 14 fwd / 12 rev | 5 fwd / 8 rev |
| Mean MAPQ of alt reads | **60.0** | **60.0** |

Recomputed 3 September 2026 by `scripts/26_verify_readlevel.py`. An earlier
version of this table gave allele 2 as depth 27 with 15/12 reads, and gave strand
splits of 15/12 and 6/10 that did not sum to their own alternate-read counts.
The figures above are internally consistent and were produced by a script that is
in the repository. Every conclusion below is unchanged.

Every indicator is what a true germline heterozygote looks like and none is what
an artefact looks like:

- **VAF near 0.5** on both. Not a mosaic fraction, not the low VAF of a
  sequencing artefact.
- **Strand balanced.** Strand bias is the commonest signature of a false
  positive and neither shows it.
- **MAPQ 60 on every alt read.** Both sit in uniquely mappable sequence, so
  neither is a paralogue mis-mapping, which is the failure mode that matters
  most in a gene family context.

## Panel coverage

| Gene | Mean depth | | Gene | Mean depth |
|---|---|---|---|---|
| CENATAC | 51.4x | | BUB1 | 46.8x |
| BUB1B | 48.1x | | CEP57L1 | 45.9x |
| BUB3 | 47.2x | | SMC5 | 45.0x |
| CEP192 | 45.7x | | CEP57 | 43.1x |
| | | | TRIP13 | 42.1x |

Genome mean 43.8x. **No known MVA gene is under-covered**, so no causal allele
is hiding in a coverage hole.

This closes a lead raised at Phase 0 for good. `TRIP13` then appeared to have a
6-fold depth deficit, which turned out to be an artefact of averaging over only
three called sites. The VCF could not settle it, because a variants-only file
cannot distinguish "no variant" from "no coverage". Real depth settles it:
42.1x against a genome mean of 43.8x. The Phase 0 closure was right, and it is
now demonstrated rather than inferred.

## PEX5 and CTU2: closed, 3 September 2026

Two homozygous calls, `PEX5` chr12:7190512 and `CTU2` chr16:88714226, were
excluded from the Track 1 submission. This section previously recorded them as
"open, not resolved". They are now resolved, and both descriptions we had given
them were wrong.

**They are common polymorphisms.** Checked against gnomAD v4.1 exomes by remote
range request:

| Call | Change | AF | group max AF | homozygotes in gnomAD |
|---|---|---|---|---|
| `PEX5` 12:7190512 | 45 bp deletion | 0.252 | 0.720 | 173,260 |
| `CTU2` 16:88714226 | 6 bp deletion | 0.767 | 0.925 | 425,713 |

The `CTU2` deletion is the major allele in the population. Neither is a
candidate for a severe recessive paediatric phenotype and neither ever was.

**Neither is loss of function.** 45 bp and 6 bp are both multiples of three, so
both are in-frame deletions.

**Our own alignment does not support the PEX5 genotype.** At `12:7190512` the
panel BAM shows 25 reads spanning the position, mean MAPQ 58.6, and **zero**
reads carrying a deletion or any other indel, against a callset genotype of 1/1.
The `CTU2` deletion is real but not homozygous: 20 of 30 reads carry it, an
allele fraction of 0.67, against a callset genotype of 1/1 with AD 0,27.

**Why they were ever described as absent from gnomAD.** They were not, by the
pipeline. No gnomAD slice fetched for this project covers chromosome 16, and none
covers the `PEX5` locus on chromosome 12, so `GnomadFrequencyAnnotator` returned
`gnomad_not_assayed` for both. That verdict is correct and is deliberately
distinct from "absent". The phrase entered in the write-up.

**And the reasoning we overturned was closer to right than what replaced it.**
Both loci are tandem repeats. gnomAD carries a 90 bp allele at the `PEX5` locus,
exactly twice the 45 bp unit, and a ladder of alleles from 2 bp to 55 bp at the
`CTU2` locus. Mapping quality near 60 shows the reads map uniquely; it says
nothing about repeat-*length* polymorphism, which is what these are and which
short reads genotype unreliably. The mapping-quality rebuttal was sound about
mapping and never touched the mechanism.
