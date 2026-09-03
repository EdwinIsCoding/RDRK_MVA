"""Unit tests for the evidence and candidate schemas.

These guard the two CLAUDE.md rules that are enforced in types rather than by
discipline: no bare genomic positions, and no unresolvable source identifiers.
Both failures are silent when they happen, which is why they are tested rather
than trusted.
"""

from __future__ import annotations

import pytest

from mva.evidence import (
    Candidate,
    Evidence,
    GenomicPosition,
    SourceType,
    VariantClass,
    Zygosity,
)


def pos(**kw) -> GenomicPosition:
    base = dict(build="GRCh38", naming="ensembl_nochr", contig="15",
                pos=40206172, ref="T", alt="C")
    return GenomicPosition(**{**base, **kw})


class TestGenomicPosition:
    def test_round_trips_between_naming_conventions(self):
        p = pos()
        assert p.contig == "15"
        q = p.to_naming("ucsc_chr")
        assert q.contig == "chr15"
        assert q.pos == p.pos, "renaming must not move the coordinate"
        assert q.to_naming("ensembl_nochr") == p

    def test_naming_must_match_the_contig_string(self):
        # The failure this prevents: a GRCh38 position labelled with the wrong
        # convention silently fails every downstream join.
        with pytest.raises(ValueError, match="ucsc_chr"):
            pos(naming="ucsc_chr", contig="15")
        with pytest.raises(ValueError, match="ensembl_nochr"):
            pos(naming="ensembl_nochr", contig="chr15")

    def test_rejects_non_nucleotide_alleles(self):
        with pytest.raises(ValueError):
            pos(ref="X")
        with pytest.raises(ValueError):
            pos(alt="hello")

    def test_accepts_symbolic_alleles_for_structural_variants(self):
        assert pos(alt="<DEL>").alt == "<DEL>"

    def test_position_is_one_based(self):
        with pytest.raises(ValueError):
            pos(pos=0)

    def test_is_hashable_and_frozen(self):
        p = pos()
        assert {p, pos()} == {p}
        with pytest.raises(Exception):
            p.pos = 5  # type: ignore[misc]


class TestEvidenceSourceValidation:
    @pytest.mark.parametrize("source_type,good", [
        (SourceType.CLINVAR, "VCV004799479.1"),
        (SourceType.PUBMED, "28492532"),
        (SourceType.ENSEMBL, "ENSG00000156970"),
        (SourceType.GO, "GO:0007094"),
        (SourceType.REACTOME, "R-HSA-69618"),
        (SourceType.OMIM, "257300"),
        (SourceType.HPO, "HP:0002859"),
        (SourceType.MONDO, "MONDO:0009759"),
        (SourceType.CHEMBL, "CHEMBL25"),
        (SourceType.DRUGBANK, "DB00945"),
        (SourceType.CLINICALTRIALS, "NCT12345678"),
    ])
    def test_accepts_well_formed_identifiers(self, source_type, good):
        e = Evidence(criterion="c", statement="s", source_type=source_type,
                     source_id=good, weight=1.0)
        assert e.source_id == good

    @pytest.mark.parametrize("source_type,bad", [
        (SourceType.CLINVAR, "VCV123"),              # too short
        (SourceType.CLINVAR, "4799479"),             # missing VCV prefix
        (SourceType.PUBMED, "not-a-pmid"),
        (SourceType.PUBMED, "123"),                  # too short to be a PMID
        (SourceType.ENSEMBL, "ENSG123"),
        (SourceType.GO, "GO:12"),
        (SourceType.OMIM, "25730"),                  # five digits
        (SourceType.HPO, "HP:123"),
        (SourceType.CLINICALTRIALS, "NCT123"),
    ])
    def test_rejects_malformed_identifiers(self, source_type, bad):
        # CLAUDE.md rule 2. A plausible-looking wrong accession is worse than a
        # gap, because a reader cannot tell it is wrong without looking it up.
        with pytest.raises(ValueError, match="does not look like"):
            Evidence(criterion="c", statement="s", source_type=source_type,
                     source_id=bad, weight=1.0)

    def test_rejects_todo_placeholder_as_an_identifier(self):
        with pytest.raises(ValueError):
            Evidence(criterion="c", statement="s", source_type=SourceType.PUBMED,
                     source_id="TODO(source)", weight=1.0)

    @pytest.mark.parametrize("code", ["PVS1", "PS3", "PM2", "PP3", "BA1", "BS1", "BP7", "PM3_strong"])
    def test_accepts_real_acmg_codes(self, code):
        e = Evidence(criterion="c", statement="s", source_type=SourceType.GO,
                     source_id="GO:0007094", weight=1.0, acmg_code=code)
        assert e.acmg_code == code

    @pytest.mark.parametrize("code", ["PXX9", "PVS2", "PM7", "pathogenic", "PP6"])
    def test_rejects_invented_acmg_codes(self, code):
        with pytest.raises(ValueError, match="ACMG"):
            Evidence(criterion="c", statement="s", source_type=SourceType.GO,
                     source_id="GO:0007094", weight=1.0, acmg_code=code)

    def test_clinvar_url_strips_zero_padding(self):
        e = Evidence(criterion="c", statement="s", source_type=SourceType.CLINVAR,
                     source_id="VCV004799479.1", weight=1.0)
        assert e.url == "https://www.ncbi.nlm.nih.gov/clinvar/variation/4799479/"

    def test_local_computation_sources_have_no_url(self):
        e = Evidence(criterion="c", statement="s", source_type=SourceType.CALLSET,
                     source_id="proband VCF FORMAT/DP", weight=-1.0)
        assert e.url is None


class TestCandidate:
    def _ev(self, w: float, st: SourceType = SourceType.GO, sid: str = "GO:0007094") -> Evidence:
        return Evidence(criterion="c", statement="s", source_type=st, source_id=sid, weight=w)

    def test_score_is_derived_not_stored(self):
        c = Candidate(gene="BUB1B", arm="A_baseline")
        assert c.score == 0.0
        c.add(self._ev(2.0)).add(self._ev(-0.5, SourceType.CALLSET, "FILTER=MQ40"))
        assert c.score == 1.5
        # There is no settable score field, so a score can never drift out of
        # step with the evidence justifying it.
        assert not hasattr(c, "_score")

    def test_negative_weights_subtract(self):
        c = Candidate(gene="X", arm="A_baseline")
        c.add(self._ev(3.0)).add(self._ev(-4.0, SourceType.GNOMAD, "gnomAD v4.1"))
        assert c.score == -1.0

    def test_collects_acmg_codes_deduplicated_and_sorted(self):
        c = Candidate(gene="X", arm="A_baseline")
        for code in ("PP3", "PM2", "PP3"):
            c.add(Evidence(criterion="c", statement="s", source_type=SourceType.GO,
                           source_id="GO:0007094", weight=1.0, acmg_code=code))
        assert c.acmg_codes == ["PM2", "PP3"]

    def test_not_reportable_without_a_falsifying_experiment(self):
        # Plan section 9: every candidate ships with the experiment that would
        # falsify it. A hypothesis with no falsification route is not a claim.
        c = Candidate(gene="BUB1B", arm="B_splicing")
        c.add(self._ev(5.0))
        ok, why = c.is_reportable()
        assert not ok and "falsifying" in why

        c.falsifying_experiment = "RT-PCR on patient fibroblast RNA"
        assert c.is_reportable() == (True, "ok")

    def test_not_reportable_with_no_evidence(self):
        c = Candidate(gene="X", arm="A_baseline", falsifying_experiment="Sanger")
        ok, why = c.is_reportable()
        assert not ok and "no evidence" in why

    def test_confidence_requires_independent_sources_not_just_score(self):
        # Four readings of one tool must not reach the same confidence as four
        # independent lines of evidence at the same total score.
        correlated = Candidate(gene="X", arm="A")
        for _ in range(4):
            correlated.add(self._ev(2.0))
        independent = Candidate(gene="Y", arm="A")
        for st, sid in [(SourceType.GO, "GO:0007094"), (SourceType.CLINVAR, "VCV004799479.1"),
                        (SourceType.PUBMED, "28492532"), (SourceType.GNOMAD, "gnomAD v4.1")]:
            independent.add(self._ev(2.0, st, sid))

        assert correlated.score == independent.score == 8.0
        assert correlated.confidence() == "low"
        assert independent.confidence() == "high"

    def test_variant_class_and_zygosity_default_to_unknown(self):
        c = Candidate(gene="X", arm="A")
        assert c.variant_class is VariantClass.UNCLASSIFIED
        assert c.zygosity is Zygosity.UNKNOWN


# ---------------------------------------------------------------------------
# Region annotation, validated against an independent source of truth.
# ---------------------------------------------------------------------------

class TestRegionAnnotationAgainstClinVarHGVS:
    """The splice distance drives the variant-class stratification that the
    whole benchmark is reported by, so it is validated against ClinVar's own
    HGVS intron offsets rather than trusted."""

    @staticmethod
    def _hgvs_intron_offset(hgvs: str) -> int | None:
        import re
        if ":c." not in hgvs:
            return None
        c = hgvs.split(":c.", 1)[1]
        if c[:1] in "*-":
            return None
        m = re.match(r"^\d+([+-])(\d+)", c)
        return int(m.group(2)) if m else None

    def _compare(self, gene_model):
        import csv

        from mva.evidence import GenomicPosition
        agree = disagree = 0
        for row in csv.DictReader(
                open("benchmarks/published_mva_variants.tsv"), delimiter="\t"):
            off = self._hgvs_intron_offset(row["hgvs_c"])
            if off is None or not (row["chrom_nochr"] and row["ref"] and row["alt"]):
                continue
            if len(row["ref"]) != len(row["alt"]):
                continue  # indels: HGVS 3'-shifts, VCF left-aligns
            try:
                p = GenomicPosition(build="GRCh38", naming="ensembl_nochr",
                                    contig=row["chrom_nochr"], pos=int(row["pos_grch38"]),
                                    ref=row["ref"], alt=row["alt"])
            except Exception:
                continue
            _, _, dist = gene_model.classify(p)
            if dist is None:
                continue
            if abs(dist - off) <= 1:
                agree += 1
            else:
                disagree += 1
        return agree, disagree

    @pytest.mark.slow
    def test_snv_splice_distances_match_clinvar_exactly(self, gene_model):
        agree, disagree = self._compare(gene_model)
        assert agree + disagree >= 100, "too few comparable rows to validate"
        concordance = agree / (agree + disagree)
        # Was 89.4% when exon boundaries were pooled across all transcripts;
        # restricting to MANE Select / Ensembl canonical made it exact.
        assert concordance == 1.0, (
            f"SNV splice-distance concordance dropped to {concordance:.3f} "
            f"({disagree} disagreements). Something changed in transcript "
            f"selection or in the GTF release."
        )

    @pytest.mark.slow
    def test_every_gene_resolves_a_representative_transcript(self, gene_model):
        missing = [g.symbol for g in gene_model.genes if not g.has_representative_transcript]
        assert not missing, (
            f"{len(missing)} genes fall back to pooled transcripts, which "
            f"systematically understates splice distance: {missing[:10]}"
        )
