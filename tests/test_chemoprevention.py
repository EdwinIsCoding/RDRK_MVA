"""Tests for the Track 2 chemoprevention derivation.

Three properties are load-bearing and each has bitten this project or its
predecessor:

1. **No dose reaches an output.** CLAUDE.md rule 3. Registry intervention names
   embed doses as a matter of course.
2. **No identifier is guessed.** CLAUDE.md rule 2. A fuzzy ChEMBL search returns
   a nearest molecule for almost any string, and a wrong ChEMBL ID in a report
   is worse than a gap because a reader cannot tell it is wrong.
3. **A heuristic is never presented as a registry fact.** The endpoint
   classification is ours, so it must be reproducible and printable.

No test here touches the network.
"""
from __future__ import annotations

import pathlib
import re

import pytest

from mva.track2 import chemoprevention as cp

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Anything that looks like a dose. Deliberately broader than the stripper, so
#: the test can fail on a pattern the stripper does not yet handle.
DOSE_PATTERN = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|µg|ug|g|ml|mL|IU|iu|units?)\b"
    r"|\b\d+(?:\.\d+)?\s*(?:mg|mcg|g)\s*/\s*(?:kg|m2|day)\b",
    re.IGNORECASE)


class TestNoDoseEscapes:
    @pytest.mark.parametrize("raw,expected", [
        ("Atorvastatin 20mg", "Atorvastatin"),
        ("Aspirin 325 mg", "Aspirin"),
        ("0.1% Uracil Cream", "Uracil"),
        ("Celecoxib 400 mg twice daily", "Celecoxib"),
        ("Sulindac 150 mg oral tablet", "Sulindac"),
        ("Metformin 500 mg/day", "Metformin"),
    ])
    def test_strip_dose_removes_the_dose(self, raw, expected):
        out = cp.strip_dose(raw)
        assert out.lower() == expected.lower(), f"{raw!r} -> {out!r}"
        assert not DOSE_PATTERN.search(out), f"a dose survived in {out!r}"

    @pytest.mark.parametrize("raw", [
        "Atorvastatin 20mg AND Aspirin 325 mg",
        "Aspirin 600 mg plus Metformin 850 mg",
    ])
    def test_combination_arms_split_into_agents(self, raw):
        parts = cp.split_combination(cp.strip_dose(raw))
        assert len(parts) == 2, parts
        for p in parts:
            assert not DOSE_PATTERN.search(p), p

    def test_generated_summary_carries_no_dose(self):
        """The whole point of the stripper, asserted on the real output."""
        out = REPO / "results" / "summaries" / "track2_chemoprevention.md"
        if not out.exists():
            pytest.skip("run scripts/27_track2_chemoprevention.py first")
        offenders = [line.strip() for line in out.read_text().splitlines()
                     if DOSE_PATTERN.search(line)]
        assert not offenders, (
            "a dose reached the chemoprevention output, which CLAUDE.md rule 3 "
            "forbids outright:\n  " + "\n  ".join(offenders[:5]))


class TestIdentifiersAreNeverGuessed:
    """resolve_molecule must reject a near miss rather than return its ChEMBL ID."""

    RESPONSE = {"molecules": [
        {"molecule_chembl_id": "CHEMBL999999", "pref_name": "DIMETHYL CELECOXIB",
         "molecule_synonyms": [], "atc_classifications": [], "max_phase": None,
         "molecule_hierarchy": {"parent_chembl_id": "CHEMBL999999"}},
        {"molecule_chembl_id": "CHEMBL118", "pref_name": "CELECOXIB",
         "molecule_synonyms": [{"molecule_synonym": "Celebrex"}],
         "atc_classifications": ["M01AH01", "L01XX33"], "max_phase": 4.0,
         "molecule_hierarchy": {"parent_chembl_id": "CHEMBL118"}},
    ]}

    @pytest.fixture(autouse=True)
    def _stub(self, monkeypatch):
        monkeypatch.setattr(cp, "_cached", lambda cache, key, fetch: self.RESPONSE)

    def test_exact_name_resolves(self, tmp_path):
        m = cp.resolve_molecule("celecoxib", tmp_path)
        assert m is not None and m.chembl_id == "CHEMBL118"
        assert m.atc_codes == ("M01AH01", "L01XX33")

    def test_synonym_resolves(self, tmp_path):
        m = cp.resolve_molecule("Celebrex", tmp_path)
        assert m is not None and m.chembl_id == "CHEMBL118"

    def test_near_miss_returns_none_rather_than_the_nearest_molecule(self, tmp_path):
        """The search response contains plausible molecules. None matches the
        query exactly, so the answer is a gap, not the first row."""
        assert cp.resolve_molecule("celecoxibb", tmp_path) is None
        assert cp.resolve_molecule("cele", tmp_path) is None


class TestSaltFormsInheritAtc:
    """A salt carrying no ATC of its own must not slip past the safety screen
    that flags its parent. ERLOTINIB was flagged and ERLOTINIB HYDROCHLORIDE
    was allowed, and they are the same active molecule."""

    SALT = {"molecules": [
        {"molecule_chembl_id": "CHEMBL1079742",
         "pref_name": "ERLOTINIB HYDROCHLORIDE", "molecule_synonyms": [],
         "atc_classifications": [], "max_phase": 4.0,
         "molecule_hierarchy": {"parent_chembl_id": "CHEMBL553"}}]}
    PARENT = {"molecule_chembl_id": "CHEMBL553", "pref_name": "ERLOTINIB",
              "atc_classifications": ["L01EB02"]}

    def test_salt_inherits_parent_atc(self, tmp_path, monkeypatch):
        def fake_cached(cache, key, fetch):
            return self.PARENT if key.startswith("molrec_") else self.SALT
        monkeypatch.setattr(cp, "_cached", fake_cached)

        m = cp.resolve_molecule("erlotinib hydrochloride", tmp_path)
        assert m is not None
        assert m.atc_codes == ("L01EB02",)
        assert m.atc_inherited_from == "CHEMBL553"

    def test_inherited_atc_reaches_the_safety_screen(self, tmp_path, monkeypatch):
        from mva.track2.safety import DrugRecord, Verdict, screen

        def fake_cached(cache, key, fetch):
            return self.PARENT if key.startswith("molrec_") else self.SALT
        monkeypatch.setattr(cp, "_cached", fake_cached)

        m = cp.resolve_molecule("erlotinib hydrochloride", tmp_path)
        r = screen(DrugRecord(name=m.pref_name, chembl_id=m.chembl_id,
                              atc_codes=m.atc_codes, provenance="ChEMBL"))
        assert r.verdict is Verdict.FLAGGED
        assert any("antineoplastic" in c for c in r.mandatory_caveats)


class TestEndpointClassification:
    @pytest.mark.parametrize("outcome,expected", [
        ("Number of new BCCs on the face at Month 12", "tumour"),
        ("Change in Duodenal Polyp Burden From Baseline to 6 Months", "tumour"),
        ("Reduction in the occurrence of any colorectal neoplasia", "tumour"),
        ("Proliferation (Ki-67) and apoptosis by immunohistochemical staining", "surrogate"),
        ("Ki67 in breast tissue of enrolled patients", "surrogate"),
        ("Whole Body Insulin Sensitivity", "other"),
        ("Bone mineral density", "other"),
    ])
    def test_endpoints_classify_as_expected(self, outcome, expected):
        assert cp.classify_endpoint(outcome) == expected

    def test_an_agent_with_no_recorded_outcome_is_not_silently_a_tumour_endpoint(self):
        ev = cp.TrialEvidence("X", ("NCT00000000",), ("Cancer",), (), ())
        assert cp.endpoint_classes(ev) == {"none recorded"}
