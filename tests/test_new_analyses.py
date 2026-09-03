"""Guards for the three analyses added on 3 September 2026.

Each closed, or partly closed, a limitation the reports had declared. The first
of them produced a wrong answer before it was calibrated, which is the reason
these exist.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
SV = REPO / "results" / "summaries" / "arm_c_sv_screen.md"
PRECEDENT = REPO / "results" / "summaries" / "kinase_domain_precedent.md"
LIT = REPO / "results" / "summaries" / "chemoprevention_literature.md"
T1 = REPO / "submission" / "track1_nexusdwin_report.md"
T2 = REPO / "submission" / "track2_nexusdwin_report.md"


class TestSvScreenIsCalibrated:
    """Uncalibrated, this screen reported clustered breakpoint evidence over
    seven of nine MVA genes. That was the ordinary rate of split-read artefact.
    The calibration is the analysis."""

    def test_the_script_refuses_to_run_without_the_extraction_intervals(self):
        src = (REPO / "scripts" / "36_sv_screen_panel.py").read_text()
        assert "REGIONS_BED" in src
        assert "sys.exit" in src.split("retrieved = load_regions")[1][:600], (
            "the screen must refuse to run without the extraction intervals; "
            "without them every unretrieved mate counts as a discordant pair")

    def test_a_background_is_computed(self):
        src = (REPO / "scripts" / "36_sv_screen_panel.py").read_text()
        assert "def background(" in src
        assert "SEED" in src, "the background sample must be seeded (CLAUDE.md rule 5)"

    def test_the_summary_reports_a_background_and_a_verdict(self):
        if not SV.exists():
            pytest.skip("run scripts/36_sv_screen_panel.py first")
        text = SV.read_text()
        assert "Calibration first" in text
        assert "99th percentile" in text
        assert "within background" in text or "exceeds background" in text

    def test_no_gene_is_silently_called_an_outlier(self):
        if not SV.exists():
            pytest.skip("screen not run")
        text = SV.read_text()
        if "exceeds background" in text:
            assert "hypothesis, not a call" in text, (
                "an outlier must be presented as a hypothesis needing inspection")

    def test_the_genome_wide_gap_is_still_declared(self):
        if not SV.exists():
            pytest.skip("screen not run")
        assert "Not genome-wide" in SV.read_text()
        assert "genome-wide" in T1.read_text().lower(), (
            "the Track 1 report must still say genome-wide SV calling was never "
            "completed; the panel screen does not replace it")


class TestKinaseDomainPrecedent:
    def test_the_precedent_is_sourced_not_recalled(self):
        if not PRECEDENT.exists():
            pytest.skip("run scripts/37_kinase_domain_precedent.py first")
        text = PRECEDENT.read_text()
        assert "PMID" in text, "the precedent must carry its citation"
        for rs in set(re.findall(r"rs\d+", text)):
            assert rs[2:].isdigit()

    def test_the_report_does_not_overclaim_from_proximity(self):
        text = T2.read_text()
        if "precedent" not in text.lower():
            pytest.skip("precedent not yet in the report")
        assert "does not make p.Asn1002Lys hypomorphic" in text, (
            "proximity to a pathogenic residue is not evidence about a "
            "different residue, and the report must say so")

    def test_the_hypomorph_claim_is_still_labelled_an_inference(self):
        text = T2.read_text()
        assert "inference" in text.lower()


class TestChemopreventionLiterature:
    def test_the_direct_question_is_asked_and_its_answer_recorded(self):
        if not LIT.exists():
            pytest.skip("run scripts/38_chemoprevention_literature.py first")
        text = LIT.read_text()
        assert "Bub1b AND chemoprevention" in text
        assert "titles" in text.lower() or "PMID" in text, (
            "hit counts alone are misleading; the titles must be shown")

    def test_the_transfer_limitation_survives(self):
        text = T2.read_text()
        low = text.lower()
        assert "transfer" in low
        assert ("dominant limitation" in low or "weakest link" in low), (
            "the literature search confirms the transfer assumption rather than "
            "removing it, and the report must still say so")

    def test_the_tetraploidy_paper_is_not_overclaimed(self):
        text = T2.read_text()
        if "24516128" not in text:
            pytest.skip("paper not cited")
        assert "not evidence for this disease" in text.lower(), (
            "tetraploidy from cytokinesis failure is not whole-chromosome "
            "aneuploidy from a weakened checkpoint")
