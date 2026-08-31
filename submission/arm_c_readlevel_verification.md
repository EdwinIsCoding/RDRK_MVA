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
| Depth (MQ>=20, BQ>=20) | 47 | 27 |
| Ref / alt reads | 21 / 26 | 15 / 12 |
| **VAF** | **0.553** | **0.444** |
| Strand balance of alt reads | 15 fwd / 12 rev | 6 fwd / 10 rev |
| Mean MAPQ of alt reads | **60.0** | **60.0** |

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

## Where a stated hypothesis of ours was wrong

Two homozygous loss-of-function calls absent from gnomAD, `PEX5` chr12:7190512
and `CTU2` chr16:88714226, were excluded from the Track 1 submission on the
reasoning that they were "almost certainly mis-calls in repetitive sequence".

**The alignment does not support that reasoning.**

| | PEX5 12:7190512 | CTU2 16:88714226 |
|---|---|---|
| Depth | 25 | 29 |
| Mean MAPQ | 58.7 | 60.0 |
| MAPQ-0 (multi-mapping) reads | 0 | 0 |
| Soft-clipped reads | 8 of 28 | 0 of 36 |

Both lie in uniquely mappable sequence. Neither shows the multi-mapping
signature that "repetitive sequence" would predict. The exclusion may still be
correct, but the stated reason for it was not.

What remains true is the clinical argument, which is independent of mapping: a
genuine homozygous `PEX5` knockout causes Zellweger spectrum disease, which this
child does not have. So the call is more likely wrong than the gene is, and the
question is *how* it is wrong. The 8 soft-clipped reads at `PEX5` hint at a
structural feature that the caller may have represented as a large homozygous
deletion; `CTU2` has no soft-clipping at all and is harder to dismiss.

**Open, not resolved.** Both need direct inspection of the read alignments and a
check of whether the indel representation is stable under left-alignment before
anything is claimed about them either way.
