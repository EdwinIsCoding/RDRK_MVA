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


class TestArmLabelsAreNotCompounds:
    """Registry intervention names describe trial arms, not molecules.

    Leaving them unnormalised cost this pipeline roughly half its candidates,
    metformin among them. Normalising them carelessly is worse: it can turn a
    placebo arm into evidence for the drug it was controlling.
    """

    @pytest.mark.parametrize("raw,expected", [
        ("metformin combination", "metformin"),
        ("celecoxib monotherapy", "celecoxib"),
        ("Sulindac drug", "Sulindac"),
        ("Early Vigabatrin", "Vigabatrin"),
        ("Delayed Vigabatrin", "Vigabatrin"),
        ("TAVT-18 sirolimus", "sirolimus"),
        ("for Aspirin 300", "Aspirin"),
        ("100 for Aspirin 100", "Aspirin"),
        ("mesalamine 5-ASA", "mesalamine"),
    ])
    def test_the_drug_survives_the_arm_wrapper(self, raw, expected):
        cands = [c.lower() for c in cp.normalise_candidates(raw)]
        assert expected.lower() in cands, f"{raw!r} did not yield {expected!r}: {cands}"

    @pytest.mark.parametrize("raw", [
        "no active patidegib", "Vehicle comparator", "Placebo",
        "matching placebo", "placebo oral tablet",
    ])
    def test_control_arms_are_identified(self, raw):
        assert cp.is_control_arm(raw), f"{raw!r} was not recognised as a control arm"

    @pytest.mark.parametrize("raw", [
        "metformin combination", "celecoxib monotherapy", "Aspirin", "sirolimus",
    ])
    def test_real_agents_are_not_mistaken_for_control_arms(self, raw):
        assert not cp.is_control_arm(raw)

    def test_a_control_arm_is_refused_before_any_lookup(self, tmp_path, monkeypatch):
        """The placebo arm of a patidegib trial must not resolve to patidegib.

        Doing so would record the arm that received no drug as registry evidence
        for the drug.
        """
        called = []
        monkeypatch.setattr(cp, "_cached",
                            lambda cache, key, fetch: called.append(key) or None)
        mol, matched = cp.resolve_agent("no active patidegib", tmp_path)
        assert mol is None and matched is None
        assert not called, "a control arm triggered a ChEMBL lookup"


class TestExactLookupBeforeFuzzySearch:
    """ChEMBL's ranked text search has poor recall for common drug names.

    Searching it for "sirolimus" returns ten molecules, none of them CHEMBL413,
    whose preferred name is SIROLIMUS. Relying on it dropped metformin from the
    candidate set entirely.
    """

    def test_exact_endpoints_are_queried_before_the_search(self, tmp_path, monkeypatch):
        keys: list[str] = []

        def fake_cached(cache, key, fetch):
            keys.append(key)
            if key.startswith("molx_pref_name_"):
                return {"molecules": [{
                    "molecule_chembl_id": "CHEMBL413", "pref_name": "SIROLIMUS",
                    "molecule_synonyms": [], "atc_classifications": ["L04AH01"],
                    "max_phase": 4.0,
                    "molecule_hierarchy": {"parent_chembl_id": "CHEMBL413"}}]}
            return {"molecules": []}

        monkeypatch.setattr(cp, "_cached", fake_cached)
        m = cp.resolve_molecule("sirolimus", tmp_path)
        assert m is not None and m.chembl_id == "CHEMBL413"
        assert keys and keys[0].startswith("molx_pref_name_"), (
            f"the exact preferred-name lookup must come first, got {keys[:1]}")
        assert not any(k.startswith("mol_") and not k.startswith("molx_")
                       for k in keys), "the fuzzy search ran despite an exact hit"

    def test_exact_endpoint_results_still_face_the_exact_name_check(self, tmp_path, monkeypatch):
        """Defence in depth: even the exact endpoint's output is name-checked,
        so a filter that ever loosened could not admit a near miss."""
        monkeypatch.setattr(cp, "_cached", lambda cache, key, fetch: {"molecules": [{
            "molecule_chembl_id": "CHEMBL999", "pref_name": "SOMETHING ELSE",
            "molecule_synonyms": [], "atc_classifications": [], "max_phase": None,
            "molecule_hierarchy": {"parent_chembl_id": "CHEMBL999"}}]})
        assert cp.resolve_molecule("sirolimus", tmp_path) is None


class TestArmsOfOneDrugMergeIntoOneCandidate:
    def test_evidence_is_unioned(self):
        a = cp.TrialEvidence("Aspirin", ("NCT00000001",), ("Lynch Syndrome",),
                             ("PHASE2",), ("polyp count",))
        b = cp.TrialEvidence("for Aspirin 300", ("NCT00000002",), ("FAP",),
                             ("PHASE3",), ("adenoma burden",))
        m = cp.merge_evidence(a, b)
        assert m.nct_ids == ("NCT00000001", "NCT00000002")
        assert set(m.conditions) == {"Lynch Syndrome", "FAP"}
        assert set(m.outcomes) == {"polyp count", "adenoma burden"}
        assert m.n_trials == 2, "merging must not lose a trial"

    def test_the_shorter_agent_label_is_kept(self):
        a = cp.TrialEvidence("100 for Aspirin 100", ("NCT1",), (), (), ())
        b = cp.TrialEvidence("Aspirin", ("NCT2",), (), (), ())
        assert cp.merge_evidence(a, b).agent == "Aspirin"

    def test_generated_summary_lists_each_molecule_once(self):
        """Unmerged arms gave one molecule two verdicts at once: ALLOWED from the
        arm whose name found paediatric trials, UNKNOWN from the one that did not."""
        out = REPO / "results" / "summaries" / "track2_chemoprevention.md"
        if not out.exists():
            pytest.skip("run scripts/27_track2_chemoprevention.py first")
        ids = re.findall(r"^\| [A-Z][^|]*\| (CHEMBL\d+) \|", out.read_text(), re.M)
        dupes = {i for i in ids if ids.count(i) > 1}
        assert not dupes, f"these molecules appear more than once: {sorted(dupes)}"


class TestTheControlCannotPassVacuously:
    """"0 of 0 controls excluded" satisfies every equality test and reads clean.

    The ATC tables live under config/ specifically so a fresh clone has them.
    They used to live under refs/, which is gitignored, so anyone cloning this
    repository would have run the safety screen with no controls at all and been
    told it discriminates. Arm B already refuses to report an unvalidated
    negative; this asserts the same discipline here.
    """

    @pytest.fixture(scope="class")
    def runner(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "runner", REPO / "scripts" / "27_track2_chemoprevention.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_the_atc_tables_are_committed(self):
        """Not a style point. The safety screen's only validation depends on
        them, and it fails open rather than closed when they are absent."""
        import subprocess
        tracked = subprocess.run(["git", "ls-files", "config/atc"], cwd=REPO,
                                 capture_output=True, text=True).stdout.split()
        assert any(f.endswith("atc_l01.json") for f in tracked), (
            "the cytotoxic control table is not committed, so a fresh clone "
            "validates the safety screen against nothing")

    def test_control_sets_are_non_empty(self, runner):
        assert runner.negative_control_agents(), "no cytotoxic controls"
        assert runner.positive_control_agents(), "no ordinary controls"

    def test_a_missing_table_yields_no_controls_rather_than_a_silent_default(
            self, runner, tmp_path, monkeypatch):
        monkeypatch.setattr(runner, "ATC_DIR", tmp_path / "absent")
        assert runner.negative_control_agents() == []
        assert runner.positive_control_agents() == []

    def test_the_summary_refuses_to_claim_validation_without_controls(self):
        """The generated summary must say NOT VALIDATED when nothing was tested,
        never report a reassuring 0 of 0."""
        out = REPO / "results" / "summaries" / "track2_chemoprevention.md"
        if not out.exists():
            pytest.skip("run scripts/27_track2_chemoprevention.py first")
        text = out.read_text()
        assert "0 of 0" not in text, (
            "the control reported a vacuous pass: zero controls, zero failures, "
            "and a clean-looking result")
        assert "controls excluded" in text
