# Gene panel provenance

## `mva_known.tsv` and `mva_known.nochr.bed`

Nine genes: the seven named in `MVA_HACKATHON_PLAN.md` section 0.1, plus
`CENATAC` and `CEP57L1`.

### Assembly method

> **Superseded, 3 September 2026.** The coordinates below came from the REST
> API and are **not** the ones in `mva_known.tsv` today. The live API returned a
> 158,779 bp span for `BUB3` against 16,072 bp in the pinned Ensembl 115 GTF, so
> the panel was rebuilt from the GTF. `BUB3` now spans 123,154,395 to
> 123,170,467. The section below is the record of the original step, kept
> because the correction is part of the provenance rather than a replacement for
> it. See `DATA_CARD.md` section 5 and `PROVENANCE.md` section 4.

Coordinates were retrieved from the Ensembl REST API in a single batch call on
31 August 2026:

```bash
curl -s https://rest.ensembl.org/lookup/symbol/homo_sapiens \
  -H 'Content-Type: application/json' -H 'Accept: application/json' \
  -d '{"symbols":["BUB1B","CEP57","TRIP13","BUB1","BUB3","CEP192","SMC5","CENATAC","CEP57L1"]}'
```

The raw response is preserved at `results/recon/ensembl_mva_genes.json`. No
coordinate in this panel was typed from memory or from a paper; every one is
traceable to that response.

Build reported by the API for all nine genes: **GRCh38**, matching
`config/config.yaml`. `seq_region_name` is used verbatim, which gives the
`ensembl_nochr` naming the patient VCF also uses, so no rename is needed for
this file. Files derived from resources that use the `chr` prefix will need one.

`mva_known.nochr.bed` applies a **5 kb flank** on each side, to include promoter
sequence, untranslated regions and the proximal intronic space where cryptic
splice alleles are found. A wider window is appropriate for the splicing arm and
should be constructed there rather than by widening this file, which is also
used for depth and density tallies where the flank size matters.

### Tiers

Tiers are an editorial judgement about strength of prior, not a property of the
data, and they are used only to order reporting.

| Tier | Genes | Basis |
|---|---|---|
| 1 | `BUB1B`, `CEP57`, `TRIP13` | Established MVA loci MVA1, MVA2 and MVA3 |
| 2 | `BUB1`, `BUB3`, `CENATAC`, `CEP192`, `SMC5` | Reported in MVA or in the overlapping near-tetraploidy and Atelis phenotypes, fewer cases |
| 3 | `CEP57L1` | Paralogue of a known gene, candidate only, no established disease association |

TODO(source): attach an OMIM number or PMID to each tier assignment. The OMIM
numbers in `mva_known.tsv` for `BUB1B`, `CEP57` and `TRIP13` came from the plan
document and were never verified against OMIM itself.

**This gap did not reach a deliverable.** The one OMIM number that appears in the
Track 1 and Track 2 reports, 257300 for mosaic variegated aneuploidy syndrome 1,
was independently confirmed against ClinVar via E-utilities during the
verification pass, along with the three OMIM numbers cited for the LZTR1
conditions. See `docs/VERIFICATION.md` section 1. The tier assignments themselves
are editorial and order reporting only; they do not enter scoring.

### Built, 31 August 2026

`mitotic_extended.tsv` exists and this section previously said it did not. It was
built by `scripts/11_build_mitotic_panel.py` and enriched with constraint by
`scripts/12_join_constraint.py`, and ten scripts and modules read it, including
the Arm A shortlist that produced the Track 1 answer. The note below was left
stale for three days; corrected 3 September 2026.

| | |
|---|---:|
| Genes | 738 |
| In the core panel (`in_core_panel = yes`) | 408 |
| Carrying a gnomAD LOEUF value | 699 |

Nominating sources, per gene and recorded in the `sources` column:

| Source | Genes |
|---|---:|
| GO | 586 |
| Reactome | 263 |
| STRING | 188 |
| Known MVA gene | 9 |

**CORUM and MGI were not used**, though plan section 6.5 lists them. The panel
was built from GO, Reactome and STRING only, and the plan's list should not be
read as a description of what was done.

The 408 figure is the one quoted as "408 mitotic genes" in the Track 1 report,
and it is the core subset rather than the whole file.
