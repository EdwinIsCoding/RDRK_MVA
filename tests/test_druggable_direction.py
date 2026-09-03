"""Tests for the axis-availability instrument.

The instrument exists to keep a Track 2 claim at the strength its evidence
supports. Its most important output is a base rate that *weakens* the project's
headline finding, so the arithmetic behind that has to be right.
"""
from __future__ import annotations

import pytest

from mva.track2 import druggable_direction as dd


class TestBaseRateArithmetic:
    def test_zero_successes_probability(self):
        # Ten draws at a 1.86% rate leaves observing zero entirely ordinary.
        assert dd.binom_zero(0.0186, 10) == pytest.approx(0.829, abs=0.01)

    def test_certainty_at_the_extremes(self):
        assert dd.binom_zero(0.0, 10) == 1.0
        assert dd.binom_zero(1.0, 10) == 0.0

    def test_more_draws_lower_the_probability_of_seeing_nothing(self):
        p = 0.05
        assert dd.binom_zero(p, 50) < dd.binom_zero(p, 10) < dd.binom_zero(p, 1)


class TestGoTermsAreNeverGuessed:
    RESPONSE = {"results": [
        {"id": "GO:0006986", "name": "response to unfolded protein"},
        {"id": "GO:0034620", "name": "cellular response to unfolded protein"},
    ]}

    @pytest.fixture(autouse=True)
    def _stub(self, monkeypatch):
        monkeypatch.setattr(dd, "_cached", lambda cache, key, fetch: self.RESPONSE)

    def test_exact_name_resolves_even_when_it_is_not_first(self, tmp_path):
        """The first result is a different, plausible term. Taking it would put a
        real GO identifier for the wrong process into a report."""
        assert dd.resolve_go_term("cellular response to unfolded protein",
                                  tmp_path) == "GO:0034620"

    def test_case_and_whitespace_are_tolerated(self, tmp_path):
        assert dd.resolve_go_term("  Cellular Response To Unfolded Protein ",
                                  tmp_path) == "GO:0034620"

    def test_a_near_miss_returns_none(self, tmp_path):
        assert dd.resolve_go_term("response to unfolded proteins", tmp_path) is None
        assert dd.resolve_go_term("unfolded protein response", tmp_path) is None


class TestEvidenceCodeFiltering:
    ANN = {"results": [
        {"symbol": "HSPA5", "goEvidence": "IDA"},
        {"symbol": "DNAJB1", "goEvidence": "IEA"},
        {"symbol": "ATF4", "goEvidence": "IMP"},
        {"symbol": "XBP1", "goEvidence": "IEA"},
    ], "pageInfo": {"total": 1, "current": 1}}

    @pytest.fixture(autouse=True)
    def _stub(self, monkeypatch):
        monkeypatch.setattr(dd, "_cached", lambda cache, key, fetch: self.ANN)

    def test_all_evidence_returns_everything(self, tmp_path):
        assert dd.go_gene_symbols("GO:0034620", tmp_path) == {
            "HSPA5", "DNAJB1", "ATF4", "XBP1"}

    def test_experimental_only_drops_electronic_annotations(self, tmp_path):
        """IEA annotations are unreviewed and dominate large GO terms. Pooling
        them silently inflates a gene set without adding evidence."""
        assert dd.go_gene_symbols("GO:0034620", tmp_path,
                                  experimental_only=True) == {"HSPA5", "ATF4"}

    def test_iea_is_not_in_the_experimental_set(self):
        assert "IEA" not in dd.EXPERIMENTAL_EVIDENCE
        assert {"IDA", "IMP", "EXP"} <= dd.EXPERIMENTAL_EVIDENCE
