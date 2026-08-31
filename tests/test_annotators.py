"""Tests for the GPU-free annotators.

Two behaviours matter most and are easy to get wrong silently:

* A missing resource must raise, not return "no evidence". Absence of a gnomAD
  frequency is not evidence of rarity, and treating it as such would promote
  every unannotated variant.
* A missense variant must be classified but not scored, because no
  pathogenicity predictor is available on this host. Scoring it on class alone
  would treat every missense as equally suspicious.
"""

from __future__ import annotations

import gzip

import pytest

from mva.evidence import Candidate, GenomicPosition, VariantClass, Zygosity
from mva.track1.annotators import ConsequenceAnnotator, GnomadFrequencyAnnotator
from mva.track1.pipeline import NotAvailableHere


def cand(pos=40206172, contig="15", ref="T", alt="C", gene="BUB1B",
         zyg=Zygosity.HET) -> Candidate:
    return Candidate(
        gene=gene, arm="A_baseline", zygosity=zyg,
        position=GenomicPosition(build="GRCh38", naming="ensembl_nochr",
                                 contig=contig, pos=pos, ref=ref, alt=alt))


@pytest.fixture
def gnomad_table(tmp_path):
    p = tmp_path / "panel_af.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("key\taf_grpmax\tnhomalt\n")
        fh.write("15:40206172:T:C\t0\t0\n")            # absent
        fh.write("15:40206173:A:G\t5e-06\t0\n")        # ultra-rare
        fh.write("15:40206174:C:T\t0.004\t120\n")      # common, many homozygotes
        fh.write("15:40206175:G:A\t2e-05\t7\n")        # rare but homozygotes exist
    return p


class TestGnomadFrequencyAnnotator:
    def test_missing_table_raises_rather_than_returning_nothing(self, tmp_path):
        a = GnomadFrequencyAnnotator(table=tmp_path / "absent.tsv.gz")
        with pytest.raises(NotAvailableHere, match="NOT"):
            a.annotate(cand(), {})

    def test_absent_from_gnomad_scores_as_rare(self, gnomad_table):
        a = GnomadFrequencyAnnotator(table=gnomad_table)
        # A variant not in the table at all is absent from gnomAD.
        ev = a.annotate(cand(pos=40299999), {})
        assert any(e.criterion == "absent_from_gnomad" for e in ev)
        assert sum(e.weight for e in ev) > 0

    def test_common_variant_is_penalised(self, gnomad_table):
        a = GnomadFrequencyAnnotator(table=gnomad_table)
        ev = a.annotate(cand(pos=40206174, ref="C", alt="T"), {})
        assert sum(e.weight for e in ev) < 0
        assert any(e.acmg_code == "BS1" for e in ev)

    def test_homozygous_proband_at_a_site_with_gnomad_homozygotes_is_heavily_penalised(
            self, gnomad_table):
        """The sharpest filter available for a severe recessive phenotype: a
        variant seen homozygous in a population reference cannot be causal in
        the homozygous state."""
        a = GnomadFrequencyAnnotator(table=gnomad_table)
        ev = a.annotate(cand(pos=40206174, ref="C", alt="T", zyg=Zygosity.HOM_ALT), {})
        hom = [e for e in ev if e.criterion == "homozygotes_in_gnomad"]
        assert hom and hom[0].weight <= -5.0
        assert hom[0].acmg_code == "BS2"

    def test_heterozygous_proband_gets_a_milder_homozygote_penalty(self, gnomad_table):
        a = GnomadFrequencyAnnotator(table=gnomad_table)
        ev = a.annotate(cand(pos=40206175, ref="G", alt="A", zyg=Zygosity.HET), {})
        hom = [e for e in ev if e.criterion == "homozygotes_in_gnomad"]
        assert hom and -5.0 < hom[0].weight < 0

    def test_rare_variant_with_no_homozygotes_is_not_penalised(self, gnomad_table):
        a = GnomadFrequencyAnnotator(table=gnomad_table)
        ev = a.annotate(cand(pos=40206173, ref="A", alt="G"), {})
        assert not any(e.criterion == "homozygotes_in_gnomad" for e in ev)
        assert sum(e.weight for e in ev) > 0


class TestConsequenceAnnotator:
    def test_missing_csq_run_raises(self):
        with pytest.raises(NotAvailableHere, match="not be reported as complete"):
            ConsequenceAnnotator().annotate(cand(), {})

    def test_loss_of_function_scores_and_cites_pvs1_supporting_only(self):
        a = ConsequenceAnnotator({"15:40206172:T:C": (VariantClass.NONSENSE, "stop_gained|BUB1B")})
        c = cand()
        ev = a.annotate(c, {})
        assert c.variant_class is VariantClass.NONSENSE
        assert ev and ev[0].weight > 0
        # PVS1 proper needs last-exon and NMD-escape checks we have not made.
        assert ev[0].acmg_code == "PVS1_supporting"

    def test_missense_is_classified_but_not_scored(self):
        a = ConsequenceAnnotator({"15:40206172:T:C": (VariantClass.MISSENSE, "missense|BUB1B")})
        c = cand()
        ev = a.annotate(c, {})
        assert c.variant_class is VariantClass.MISSENSE
        assert ev and ev[0].weight == 0.0, (
            "missense must not be scored: no pathogenicity predictor is "
            "available, so scoring on class alone would treat every missense "
            "as equally suspicious"
        )

    def test_csq_overrides_the_coordinate_derived_class(self):
        a = ConsequenceAnnotator({"15:40206172:T:C": (VariantClass.SYNONYMOUS, "synonymous|BUB1B")})
        c = cand()
        c.variant_class = VariantClass.DEEP_INTRONIC   # from regions.py
        a.annotate(c, {})
        assert c.variant_class is VariantClass.SYNONYMOUS

    def test_unannotated_position_returns_nothing_without_raising(self):
        a = ConsequenceAnnotator({"1:1:A:C": (VariantClass.MISSENSE, "x")})
        assert a.annotate(cand(), {}) == []

    def test_occupies_the_vep_slot_so_completeness_reporting_works(self):
        # PipelineResult.REQUIRED_FOR_COMPLETE keys on annotator name.
        assert ConsequenceAnnotator().name == "vep"
