"""Tests for directional target nomination and tractability.

The failure these guard is plan section 0.4's named anti-goal: using unsigned
network proximity to decide whether to inhibit or activate. That error is
invisible in the output, because a target with a confidently stated wrong
direction looks exactly like one with a right direction.
"""

from __future__ import annotations

from mva.track2.targets import Direction, SignedEdge, nominate, summarise
from mva.track2.tractability import (
    ACTIVATING_ACTIONS,
    INHIBITING_ACTIONS,
    Tractability,
    TractabilityResult,
)


def edge(src="PLK1", tgt="BUB1B", stim=False, inh=False,
         cons_stim=False, cons_inh=False, effort=5, pmids=("12345678",)):
    return SignedEdge(source=src, target=tgt, is_stimulation=stim, is_inhibition=inh,
                      consensus_stimulation=cons_stim, consensus_inhibition=cons_inh,
                      sources=("SIGNOR",), pmids=pmids, curation_effort=effort)


class TestEdgeSignResolution:
    def test_consensus_wins_over_raw_flags(self):
        # The raw flags are a union across databases and are routinely both
        # true; the consensus fields are the resolution.
        e = edge(stim=True, inh=True, cons_stim=False, cons_inh=True)
        assert e.sign is Direction.INHIBIT

    def test_contradictory_consensus_is_ambiguous(self):
        assert edge(cons_stim=True, cons_inh=True).sign is Direction.AMBIGUOUS

    def test_falls_back_to_raw_flags_when_consensus_silent(self):
        assert edge(stim=True).sign is Direction.ACTIVATE
        assert edge(inh=True).sign is Direction.INHIBIT

    def test_contradictory_raw_flags_are_ambiguous_not_a_coin_flip(self):
        assert edge(stim=True, inh=True).sign is Direction.AMBIGUOUS

    def test_no_sign_at_all_is_unsigned(self):
        assert edge().sign is Direction.UNSIGNED


class TestNomination:
    SEEDS = {"BUB1B", "BUB1"}

    def test_activator_of_seed_is_nominated_for_activation(self):
        n = nominate("PLK1", [edge(cons_stim=True)], self.SEEDS,
                     desired_effect_on_seed=Direction.ACTIVATE)
        assert n.direction is Direction.ACTIVATE
        assert n.is_nominable

    def test_inhibitor_of_seed_is_nominated_for_inhibition(self):
        # To raise BUB1B activity, inhibit what inhibits it.
        n = nominate("X", [edge(src="X", cons_inh=True)], self.SEEDS,
                     desired_effect_on_seed=Direction.ACTIVATE)
        assert n.direction is Direction.INHIBIT

    def test_unsigned_edge_is_not_nominable(self):
        # Plan section 0.4. Proximity is not a direction.
        n = nominate("X", [edge(src="X")], self.SEEDS)
        assert not n.is_nominable
        assert n.direction is Direction.UNSIGNED
        assert "unsigned" in n.blocking_reason

    def test_ambiguous_edge_is_not_nominable(self):
        n = nominate("X", [edge(src="X", cons_stim=True, cons_inh=True)], self.SEEDS)
        assert not n.is_nominable
        assert "contradict" in n.blocking_reason

    def test_gene_that_both_activates_and_inhibits_seeds_is_rejected(self):
        n = nominate("X", [edge(src="X", tgt="BUB1B", cons_stim=True),
                           edge(src="X", tgt="BUB1", cons_inh=True)], self.SEEDS)
        assert n.direction is Direction.AMBIGUOUS
        assert not n.is_nominable

    def test_no_edge_to_a_seed_is_rejected(self):
        n = nominate("X", [edge(src="X", tgt="UNRELATED", cons_stim=True)], self.SEEDS)
        assert n.direction is Direction.NO_EDGE

    def test_nomination_carries_resolvable_pmids(self):
        n = nominate("PLK1", [edge(cons_stim=True, pmids=("17376779", "17785528"))],
                     self.SEEDS)
        assert n.evidence
        for ev in n.evidence:
            assert ev.source_id.isdigit()
            assert ev.url.startswith("https://pubmed.ncbi.nlm.nih.gov/")

    def test_summary_reports_rejections_rather_than_hiding_them(self):
        noms = [
            nominate("A", [edge(src="A", cons_stim=True)], self.SEEDS),
            nominate("B", [edge(src="B")], self.SEEDS),
            nominate("C", [edge(src="C", cons_stim=True, cons_inh=True)], self.SEEDS),
        ]
        s = summarise(noms)
        assert s["n_total"] == 3
        assert s["n_nominable"] == 1
        assert s["rejected_unsigned"] == ["B"]
        assert s["rejected_ambiguous"] == ["C"]


class TestTractabilityClassification:
    def _r(self, direction, n_act, n_inh, n_mech=None, actions=None):
        n_mech = n_mech if n_mech is not None else n_act + n_inh
        return TractabilityResult(
            gene="X", chembl_target_id="CHEMBL1", required_direction=direction,
            status=(Tractability.NO_DRUGS if n_mech == 0 else
                    Tractability.AVAILABLE
                    if (n_act if direction is Direction.ACTIVATE else n_inh) > 0
                    else Tractability.WRONG_DIRECTION_ONLY),
            n_mechanisms=n_mech, n_activating=n_act, n_inhibiting=n_inh,
            max_phase=3.0, action_types=actions or {}, example_drugs=())

    def test_inhibitors_do_not_satisfy_a_requirement_to_activate(self):
        """The finding this encodes: 118 mechanisms across the SAC targets, all
        inhibitors, against a requirement to activate."""
        r = self._r(Direction.ACTIVATE, n_act=0, n_inh=118,
                    actions={"INHIBITOR": 118})
        assert r.status is Tractability.WRONG_DIRECTION_ONLY
        assert not r.is_actionable
        assert "wrong way" in r.explanation

    def test_activators_satisfy_a_requirement_to_activate(self):
        r = self._r(Direction.ACTIVATE, n_act=3, n_inh=0)
        assert r.is_actionable

    def test_inhibitors_satisfy_a_requirement_to_inhibit(self):
        r = self._r(Direction.INHIBIT, n_act=0, n_inh=9)
        assert r.is_actionable

    def test_no_drugs_is_distinct_from_wrong_direction(self):
        assert self._r(Direction.ACTIVATE, 0, 0).status is Tractability.NO_DRUGS

    def test_action_vocabularies_are_disjoint(self):
        # An action type counted as both would make every target look tractable.
        assert not (ACTIVATING_ACTIONS & INHIBITING_ACTIONS)


class TestPipelineCompletenessReporting:
    """An incomplete run must announce itself.

    The failure guarded here actually happened: a hardcoded spike-in INFO tag
    made bcftools error on the proband VCF, and the pipeline reported
    "0 records seen, 0 candidates, annotators ran fine" -- indistinguishable
    from a legitimate negative result.
    """

    def _result(self, ran):
        from mva.track1.pipeline import PipelineResult
        return PipelineResult(candidates=[], n_records_seen=0, n_records_in_panel=0,
                              annotators_run=ran, annotators_unavailable={})

    def test_default_annotator_set_is_reported_as_incomplete(self):
        r = self._result(["panel", "region", "quality"])
        assert not r.is_complete
        assert set(r.missing_evidence_classes) == {"vep", "gnomad", "spliceai"}
        note = r.completeness_note()
        assert note.startswith("INCOMPLETE RUN")
        assert "must not be presented as Track 1 findings" in note

    def test_full_annotator_set_is_complete(self):
        r = self._result(["panel", "region", "quality", "vep", "gnomad", "spliceai"])
        assert r.is_complete
        assert r.completeness_note().startswith("Complete run")

    def test_missing_spliceai_is_called_out_as_the_highest_prior_arm(self):
        r = self._result(["panel", "region", "quality", "vep", "gnomad"])
        assert "highest prior" in r.missing_evidence_classes["spliceai"]
