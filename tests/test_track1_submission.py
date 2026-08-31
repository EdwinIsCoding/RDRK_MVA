"""Tests for the Track 1 submission format.

Only 6 submissions are allowed per participant and the leaderboard scores
automatically, so a format error is not a warning, it is a wasted attempt
against a hard cap. The two conventions most likely to go wrong are chromosome
prefixing (the callset has none, the submission requires it) and compound-het
pairs (which belong in one row, since partial credit is available for
recovering one of two).
"""
from __future__ import annotations

import csv
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "t1", pathlib.Path(__file__).resolve().parents[1] / "scripts/23_track1_submission.py")
t1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t1)

import pytest


def row(**kw):
    base = dict(proband_id="WGS_EX2312012", chrom_1="chr15", pos_1="40206172",
                ref_1="T", alt_1="C", chrom_2="", pos_2="", ref_2="", alt_2="",
                epcr="0.9", finding_type="primary", notes="")
    return {**base, **kw}


class TestChromosomeConversion:
    def test_ensembl_naming_is_converted(self):
        assert t1.to_ucsc("15") == "chr15"
        assert t1.to_ucsc("X") == "chrX"

    def test_already_prefixed_is_left_alone(self):
        assert t1.to_ucsc("chr15") == "chr15"

    def test_mitochondrion_normalises(self):
        assert t1.to_ucsc("MT") == "chrM"
        assert t1.to_ucsc("M") == "chrM"

    def test_unprefixed_chromosome_is_rejected_by_validation(self):
        # The exact failure that would silently score zero.
        with pytest.raises(t1.FormatError, match="not chr-prefixed"):
            t1.validate([row(chrom_1="15")])


class TestFormatValidation:
    def test_a_well_formed_row_passes(self):
        t1.validate([row()])

    def test_row_limit_enforced(self):
        rows = [row(epcr=str(1 - i / 100)) for i in range(11)]
        with pytest.raises(t1.FormatError, match="limit is 10"):
            t1.validate(rows)

    @pytest.mark.parametrize("epcr", ["0", "0.0", "1.5", "-0.2", "abc"])
    def test_epcr_must_be_in_the_half_open_unit_interval(self, epcr):
        with pytest.raises(t1.FormatError):
            t1.validate([row(epcr=epcr)])

    def test_epcr_of_exactly_one_is_allowed(self):
        t1.validate([row(epcr="1.0")])

    def test_finding_type_is_constrained(self):
        with pytest.raises(t1.FormatError, match="primary"):
            t1.validate([row(finding_type="maybe")])

    def test_at_least_one_primary_required(self):
        with pytest.raises(t1.FormatError, match="no primary"):
            t1.validate([row(finding_type="secondary")])

    def test_secondary_findings_are_allowed_alongside_a_primary(self):
        t1.validate([row(epcr="0.9"),
                     row(chrom_1="chr17", pos_1="43044295", ref_1="A", alt_1="G",
                         epcr="0.2", finding_type="secondary",
                         notes="well-established pathogenic variant, unrelated to phenotype")])

    def test_partial_compound_het_pair_is_rejected(self):
        # Half a pair loses the partial credit the metric offers and is
        # ambiguous to the scorer.
        with pytest.raises(t1.FormatError, match="partial second variant"):
            t1.validate([row(chrom_2="chr15", pos_2="40206200")])

    def test_complete_compound_het_pair_passes(self):
        t1.validate([row(chrom_2="chr15", pos_2="40206200", ref_2="G", alt_2="A")])

    def test_alleles_must_be_plain(self):
        with pytest.raises(t1.FormatError, match="not a plain allele"):
            t1.validate([row(ref_1="<DEL>")])

    def test_position_must_be_a_positive_integer(self):
        with pytest.raises(t1.FormatError, match="positive integer"):
            t1.validate([row(pos_1="0")])
        with pytest.raises(t1.FormatError, match="positive integer"):
            t1.validate([row(pos_1="chr15:40206172")])

    def test_empty_submission_rejected(self):
        with pytest.raises(t1.FormatError, match="no candidate rows"):
            t1.validate([])


class TestOutput:
    def test_rows_are_written_in_descending_epcr(self, tmp_path):
        out = tmp_path / "sub.csv"
        t1.write([row(epcr="0.3"), row(epcr="0.9", pos_1="40206180"),
                  row(epcr="0.6", pos_1="40206190")], out)
        rows = list(csv.DictReader(out.open(newline="")))
        assert [r["epcr"] for r in rows] == ["0.9", "0.6", "0.3"]

    def test_header_matches_the_specification_exactly(self, tmp_path):
        out = tmp_path / "sub.csv"
        t1.write([row()], out)
        assert next(csv.reader(out.open(newline=""))) == t1.COLUMNS
