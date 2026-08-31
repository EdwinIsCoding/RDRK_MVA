# Gene panel provenance

## `mva_known.tsv` and `mva_known.nochr.bed`

Nine genes: the seven named in `MVA_HACKATHON_PLAN.md` section 0.1, plus
`CENATAC` and `CEP57L1`.

### Assembly method

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

TODO(source): attach an OMIM number or PMID to each tier assignment before the
panel is used in scoring. The OMIM numbers in `mva_known.tsv` for `BUB1B`,
`CEP57` and `TRIP13` came from the plan document and have not yet been verified
against OMIM itself.

### Not yet done

`mitotic_extended.tsv`, the 300 to 500 gene panel described in plan section 6.5,
has not been built. It is Phase 3 work (kickoff prompt P3) and requires GO,
Reactome, CORUM, STRING and MGI sources with per-gene attribution of which
source nominated it, plus gnomAD LOEUF and pLI.
