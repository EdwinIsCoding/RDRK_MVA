"""Tests for the Track 2 safety screen.

The screen is the part of Track 2 most likely to be got wrong quietly, because
a wrong ALLOWED looks exactly like a right one. These tests pin the behaviours
that plan section 7.4 calls load-bearing, using the specific compound classes
the plan names.
"""

from __future__ import annotations

import pytest

from mva.track2.safety import (
    DrugRecord,
    Verdict,
    screen,
)


class TestGenotoxicExclusion:
    """Categorical. No efficacy argument overrides these."""

    def test_explicit_genotoxic_flag_excludes(self):
        r = screen(DrugRecord(name="X", is_genotoxic=True, provenance="DrugBank"))
        assert r.verdict is Verdict.EXCLUDED
        assert not r.may_be_proposed

    @pytest.mark.parametrize("mechanism", [
        "alkylating agent",
        "topoisomerase II inhibitor",
        "DNA cross-linking agent",
        "platinum-based agent",
        "nucleoside analogue",
        "radiomimetic",
    ])
    def test_genotoxic_mechanisms_excluded_from_free_text(self, mechanism):
        r = screen(DrugRecord(name="X", mechanism=mechanism, provenance="ChEMBL"))
        assert r.verdict is Verdict.EXCLUDED, f"{mechanism!r} was not excluded"

    @pytest.mark.parametrize("atc,group", [
        ("L01AA01", "L01A alkylating"),
        ("L01BA01", "L01B antimetabolite"),
        ("L01CA01", "L01C plant alkaloid"),
        ("L01DB01", "L01D cytotoxic antibiotic"),
    ])
    def test_cytotoxic_antineoplastic_atc_excluded(self, atc, group):
        r = screen(DrugRecord(name="X", atc_codes=(atc,), provenance="WHO ATC"))
        assert r.verdict is Verdict.EXCLUDED, f"{group} was not excluded"

    @pytest.mark.parametrize("atc", ["L01EX01", "L01FF01", "L01XX33"])
    def test_non_cytotoxic_antineoplastic_atc_flagged_not_excluded(self, atc):
        """The blanket L01 rule excluded celecoxib, which carries L01XX33.

        Celecoxib is the best-evidenced chemoprevention agent in hereditary
        cancer predisposition, and the stated reason for excluding it was that
        cytotoxic chemotherapy is out of scope. A COX-2 inhibitor is not
        cytotoxic chemotherapy. The exclusion is now scoped to L01A-L01D.
        """
        r = screen(DrugRecord(name="X", atc_codes=(atc,), provenance="WHO ATC"))
        assert r.verdict is not Verdict.EXCLUDED
        assert r.may_be_proposed
        assert any("antineoplastic" in c for c in r.mandatory_caveats), (
            "a non-cytotoxic antineoplastic must still carry its classification "
            "as a caveat, never pass silently")

    def test_celecoxib_atc_pair_is_not_excluded(self):
        """The real ATC pair for CHEMBL118, retrieved from ChEMBL."""
        r = screen(DrugRecord(name="celecoxib", chembl_id="CHEMBL118",
                              atc_codes=("M01AH01", "L01XX33"),
                              provenance="ChEMBL"))
        assert r.verdict is not Verdict.EXCLUDED

    def test_in_vitro_chromosomal_instability_excluded(self):
        # This is the disease mechanism itself.
        r = screen(DrugRecord(name="X", causes_chromosomal_instability_in_vitro=True,
                              provenance="literature"))
        assert r.verdict is Verdict.EXCLUDED

    def test_exclusion_survives_strong_paediatric_evidence(self):
        # A compound with excellent paediatric data is still excluded if it is
        # genotoxic. Severity ordering must not let ALLOWED outrank EXCLUDED.
        r = screen(DrugRecord(
            name="X", is_genotoxic=True, has_paediatric_label=True,
            has_paediatric_pk=True, paediatric_trial_ids=("NCT12345678",),
            provenance="DrugBank"))
        assert r.verdict is Verdict.EXCLUDED


class TestImmunosuppressionIsFlaggedNotSilent:
    """Plan section 7.4: naming the tension is a strength, ignoring it is the
    error a clinically trained judge spots first."""

    @pytest.mark.parametrize("record", [
        DrugRecord(name="everolimus", mechanism="mTOR inhibitor", provenance="DrugBank"),
        DrugRecord(name="sirolimus", mechanism="rapalog", provenance="DrugBank"),
        DrugRecord(name="ciclosporin", mechanism="calcineurin inhibitor", provenance="DrugBank"),
        DrugRecord(name="X", is_immunosuppressant=True, provenance="DrugBank"),
        DrugRecord(name="Y", atc_codes=("L04AA18",), provenance="WHO ATC"),
    ])
    def test_immunosuppressants_are_flagged(self, record):
        r = screen(record)
        assert r.verdict is Verdict.FLAGGED

    def test_everolimus_may_be_proposed_but_never_silently(self):
        # The plan calls out everolimus by name: it is a plausible axis via
        # proteostasis, and it is immunosuppressive in a cancer-predisposed
        # child. It may be proposed, but the tension must travel with it.
        r = screen(DrugRecord(name="everolimus", mechanism="mTOR inhibitor",
                              has_paediatric_label=True, provenance="DrugBank"))
        assert r.may_be_proposed
        assert r.verdict is Verdict.FLAGGED
        caveats = " ".join(r.mandatory_caveats).lower()
        assert "immune surveillance" in caveats
        assert r.mandatory_caveats, "a flagged compound must carry its caveat"


class TestUnknownIsNotAPass:
    """The failure this guards: collapsing 'we did not look' into 'safe'."""

    def test_no_paediatric_data_looked_up_is_unknown(self):
        r = screen(DrugRecord(name="X", mechanism="HSF1 activator", provenance="ChEMBL"))
        assert r.verdict is Verdict.UNKNOWN
        assert r.may_be_proposed, "unknown is a caveat, not an exclusion"
        assert r.mandatory_caveats

    def test_looked_up_and_absent_is_still_unknown(self):
        r = screen(DrugRecord(name="X", mechanism="HSF1 activator",
                              has_paediatric_label=False, has_paediatric_pk=False,
                              provenance="DrugBank"))
        assert r.verdict is Verdict.UNKNOWN
        assert "unestablished" in " ".join(r.mandatory_caveats)

    def test_paediatric_evidence_reaches_allowed(self):
        r = screen(DrugRecord(name="X", mechanism="HSF1 activator",
                              has_paediatric_label=True, has_paediatric_pk=True,
                              paediatric_trial_ids=("NCT12345678",),
                              provenance="DrugBank+ClinicalTrials"))
        assert r.verdict is Verdict.ALLOWED
        assert not r.mandatory_caveats


class TestCnsPenetrance:
    def test_only_applies_when_the_target_is_cns(self):
        d = DrugRecord(name="X", crosses_bbb=False, has_paediatric_label=True,
                       provenance="DrugBank")
        assert screen(d, cns_target=False).verdict is Verdict.ALLOWED
        assert screen(d, cns_target=True).verdict is Verdict.EXCLUDED

    def test_unknown_bbb_with_cns_target_is_unknown(self):
        d = DrugRecord(name="X", has_paediatric_label=True, provenance="DrugBank")
        assert screen(d, cns_target=True).verdict is Verdict.UNKNOWN


class TestDeterminism:
    def test_same_input_gives_same_verdict_and_findings(self):
        d = DrugRecord(name="everolimus", mechanism="mTOR inhibitor",
                       has_paediatric_label=True, provenance="DrugBank")
        a, b = screen(d), screen(d)
        assert a == b, "the screen must be a pure function of its inputs"

    def test_every_finding_names_its_rule_and_field(self):
        # Auditability: a verdict nobody can trace to a rule and a field is not
        # usable in a submission that claims reproducibility.
        d = DrugRecord(name="everolimus", mechanism="mTOR inhibitor", provenance="DrugBank")
        for f in screen(d).findings:
            assert f.rule_id.startswith("SAFE-")
            assert f.field_read
            assert f.source
