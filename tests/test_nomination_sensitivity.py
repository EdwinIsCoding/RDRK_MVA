"""Guards for the sensitivity analysis on the all-activation finding.

A sensitivity analysis that is only reported when it confirms is not a
sensitivity analysis. These tests defend both halves of what it found: that the
statistical claim survives, and that the scope claim had to narrow.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "summaries" / "track2_nomination_sensitivity.md"
SCRIPT = REPO / "scripts" / "42_nomination_sensitivity.py"
T2 = REPO / "submission" / "track2_nexusdwin_report.md"


@pytest.fixture(scope="module")
def out():
    if not OUT.exists():
        pytest.skip("run scripts/42_nomination_sensitivity.py first")
    return OUT.read_text()


class TestTheNullIsMeasuredNotAssumed:
    def test_the_stimulatory_share_is_reported(self, out):
        m = re.search(r"stimulatory share of the graph, measured\*?\*? \| \*\*(\d+)%", out)
        assert m, "the measured stimulatory share is no longer reported"
        share = int(m.group(1))
        assert 0 < share < 100

    def test_the_probability_is_reported_and_consistent(self, out):
        m = re.search(r"stimulatory share of (\d+)%, seeing (\d+) of (\d+) require "
                      r"activation has probability (\d+\.\d+)", out)
        assert m, "the null probability is no longer reported"
        share, n, _, p = int(m.group(1)) / 100, int(m.group(2)), None, float(m.group(4))
        assert abs(share ** n - p) < 0.01, (
            f"the reported probability {p} does not follow from {share} ** {n}")

    def test_the_control_is_seeded(self):
        assert "SEED" in SCRIPT.read_text()

    def test_both_outcomes_are_written(self):
        """The script must be able to report that the claim failed."""
        src = " ".join(SCRIPT.read_text().split())
        src = re.sub(r'"\s*"', "", src)
        assert "the claim must be weakened" in src, (
            "the branch that weakens the claim on a high null was removed")
        assert "The claim stands as written" in src


class TestTheScopeCorrectionSurvives:
    """Widening the seed set reaches inhibition targets. The report must say the
    finding belongs to the immediate regulators, not to mitotic biology."""

    def test_the_report_scopes_the_claim(self):
        text = T2.read_text()
        low = text.lower()
        assert "immediate regulator" in low, (
            "the report no longer scopes the all-activation claim to the "
            "immediate regulators of the seed genes")
        assert "not to mitotic biology" in low or "not mitotic biology" in low

    def test_the_report_does_not_claim_no_mix_exists(self):
        text = T2.read_text()
        assert "This one does not, and six of the ten" not in text, (
            "the unscoped version of claim 2 has returned")

    def test_the_wider_targets_are_not_offered_as_an_alternative(self, out):
        low = out.lower()
        assert "contraindicated place" in low or "same contraindicated" in low, (
            "the summary must say the wider inhibition targets lead to "
            "activating a mitotic kinase, which safety rules out")

    def test_the_report_explains_why_the_axis_stays_closed(self):
        low = T2.read_text().lower()
        assert "over-determined" in low, (
            "the report must say the closure holds on both availability and "
            "safety, not only availability")


class TestMalformedReferencesDoNotBecomeFindings:
    """OmniPath serves a placeholder where a PubMed identifier should be, and
    "123" appears in the wild. The Evidence validator refuses it. A malformed
    reference in one edge must not abandon the analysis, and must not be
    silently replaced with a plausible number either."""

    def test_the_nomination_survives_a_bad_reference(self):
        src = (REPO / "src" / "mva" / "track2" / "targets.py").read_text()
        assert "ValidationError" in src
        assert "malformed" in src

    def test_an_edge_with_no_usable_reference_is_refused(self):
        src = (REPO / "src" / "mva" / "track2" / "targets.py").read_text()
        assert "direction=None" in src, (
            "an edge whose every reference is malformed has no citable basis "
            "and its direction must not be asserted")
