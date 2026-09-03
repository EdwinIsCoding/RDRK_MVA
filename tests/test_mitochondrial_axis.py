"""The mitochondrial axis follow-through must not overclaim.

Section 5 of the Track 2 report called this axis the most promising of the three
and then stopped. Following through is only an improvement if the result is
reported at the strength it deserves, which for this axis is low: being better
supplied with activating drugs than the genome average is a weak thing to be.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "summaries" / "track2_mitochondrial_axis.md"
AXIS = REPO / "results" / "summaries" / "track2_axis_availability.md"
T2 = REPO / "submission" / "track2_nexusdwin_report.md"
SCRIPT = REPO / "scripts" / "40_mitochondrial_axis_followthrough.py"


@pytest.fixture(scope="module")
def out():
    if not OUT.exists():
        pytest.skip("run scripts/40_mitochondrial_axis_followthrough.py first")
    return OUT.read_text()


class TestTheGeneListComesFromTheOtherAnalysis:
    def test_genes_are_read_not_retyped(self):
        """If the list were retyped, the two analyses could disagree silently."""
        src = SCRIPT.read_text()
        assert "AXIS_SUMMARY" in src
        assert "track2_axis_availability" in src

    def test_the_gene_lists_agree(self, out):
        if not AXIS.exists():
            pytest.skip("axis availability not run")
        block = AXIS.read_text().split("### Mitochondrial and oxidative support")[1]
        line = next((ln for ln in block.splitlines()
                     if ln.startswith("Genes in this axis with a drug")), "")
        upstream = set(re.findall(r"`([A-Z0-9]+)`", line))
        assert upstream, "could not parse the upstream gene list"
        missing = sorted(g for g in upstream if g not in out)
        assert not missing, f"genes dropped between the two analyses: {missing}"


class TestVacuousGatesAreDeclaredNotPerformed:
    """Open Targets tractability and Pharos development level cannot fail for a
    set assembled by requiring an existing drug. Running them would return a
    pass for everything and look like evidence."""

    def test_the_vacuous_gates_are_named_and_explained(self, out):
        low = out.lower()
        assert "tractability" in low and "pharos" in low
        assert "vacuous" in low, (
            "the two gates that cannot bite must be declared vacuous, with the "
            "reason, rather than quietly skipped or performed for show")

    def test_the_gate_that_does_bite_is_applied(self, out):
        assert "GTEx" in out
        assert "Muscle_Skeletal" in out or "Skeletal muscle" in out
        assert "Kidney" in out

    def test_bbb_is_declined_with_a_reason(self, out):
        low = out.lower()
        assert "blood-brain" in low
        assert "no cns" in low or "seizures" in low


class TestNoThresholdIsInvented:
    def test_expression_reports_values_not_a_verdict(self, out):
        assert "No cutoff of ours is applied" in out or "no threshold" in out.lower(), (
            "an expression cutoff recited from memory is exactly what "
            "CLAUDE.md rule 2 forbids")

    def test_the_highest_tissue_is_shown_for_context(self, out):
        assert "Highest tissue" in out or "highest-expressing" in out


class TestTheResultIsNotOversold:
    def test_the_low_bar_is_stated(self, out):
        low = out.lower()
        assert "low bar" in low or "weak" in low, (
            "following through on this axis is only honest if the weakness of "
            "the mechanistic case survives into the conclusion")

    def test_a_falsification_route_is_given(self, out):
        low = out.lower()
        assert "falsif" in low
        assert "bub1b" in low or "hypomorphic mouse" in low

    def test_no_dosing(self, out):
        dose = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|IU)\b", re.I)
        hits = [line.strip() for line in out.splitlines() if dose.search(line)]
        assert not hits, f"dosing reached the output: {hits[:3]}"

    def test_no_em_dash(self, out):
        assert "\u2014" not in out

    def test_excluded_drugs_carry_their_rule(self, out):
        if "### Excluded, and why" not in out:
            pytest.skip("nothing excluded")
        assert "SAFE-" in out, "an exclusion must name the rule that produced it"


class TestTheReportFollowsThrough:
    def test_section_5_no_longer_stops_at_promising(self):
        text = T2.read_text()
        if "mitochondrial_axis" not in text and "section 5a" not in text.lower():
            pytest.skip("report not yet updated")
        assert "track2_mitochondrial_axis" in text, (
            "the report should point at the follow-through it promised")


class TestTheReportCitesTheStableFinding:
    """Exact pair counts move between runs because ChEMBL serves intermittent
    500s. The concentration in three genes does not, so that is what the report
    quotes and what this test checks."""

    def test_the_three_genes_are_named_in_both(self, out):
        report = T2.read_text()
        if "mitochondrial_axis" not in report:
            pytest.skip("report not yet updated")
        for gene in ("INSR", "PPARA", "GCK"):
            assert gene in out, f"{gene} missing from the summary"
            assert gene in report, f"{gene} missing from the report"

    def test_the_concentration_claim_matches_the_data(self, out):
        report = T2.read_text()
        if "roughly 70%" not in report:
            pytest.skip("report does not make the approximate claim")
        m = re.search(r"\*\*(\d+) of (\d+) screened pairs, (\d+)%", out)
        assert m, "could not parse the concentration from the summary"
        share = int(m.group(3))
        assert 60 <= share <= 80, (
            f"the summary now reports {share}% concentration; the report says "
            f"roughly 70% and one of them needs changing")

    def test_the_instability_is_disclosed(self, out):
        low = out.lower()
        assert "500" in out or "intermittent" in low
        assert "cache" in low, (
            "reproducibility depends on the cache and the summary must say so")
