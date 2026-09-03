"""The deletion plan must be safe by default and must match ETHICS.md.

A destructive script is the one place where a bug is unrecoverable, and the
inputs here cannot be re-obtained once the challenge closes.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "33_delete_challenge_data.py"


@pytest.fixture(scope="module")
def mod():
    if not SCRIPT.exists():
        pytest.skip("deletion script absent")
    spec = importlib.util.spec_from_file_location("deletion", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class TestItIsSafeByDefault:
    def test_deletion_requires_an_explicit_flag(self):
        src = SCRIPT.read_text()
        assert '"--execute"' in src
        assert "if not args.execute:" in src, "there is no dry-run guard"

    def test_deletion_requires_a_typed_confirmation(self):
        src = SCRIPT.read_text()
        assert 'reply != "DELETE"' in src, (
            "a flag alone should not be enough to destroy 80 GB that cannot be "
            "re-obtained")

    def test_reference_data_is_not_deleted_by_default(self, mod):
        assert any(rel == "refs" for rel, _ in mod.OPTIONAL)
        assert not any(rel == "refs" for rel, _ in mod.DELETE)


class TestThePlanMatchesTheEthicsDocument:
    def test_every_delete_target_is_challenge_data_or_derived(self, mod):
        allowed_roots = {"data", "results", "node_artefacts"}
        for rel, _ in mod.DELETE:
            root = rel.split("/")[0]
            assert root in allowed_roots, (
                f"{rel} is outside the delete list in ETHICS.md section 3b")

    def test_nothing_on_the_keep_list_is_also_on_the_delete_list(self, mod):
        delete = {rel for rel, _ in mod.DELETE}
        for rel, _ in mod.KEEP:
            assert rel not in delete, f"{rel} is on both lists"

    def test_the_submission_is_kept(self, mod):
        """Candidate variant rankings are explicitly on the organisers' keep
        list. Deleting the submission would destroy the deliverable."""
        assert any(rel == "submission" for rel, _ in mod.KEEP)
        assert not any(rel.startswith("submission") for rel, _ in mod.DELETE)

    def test_aggregate_summaries_are_kept_but_raw_tables_are_not(self, mod):
        delete = {rel for rel, _ in mod.DELETE}
        assert any(rel == "results/summaries" for rel, _ in mod.KEEP)
        assert "results/summaries/arm_f_roh.raw.txt" in delete, (
            "the runs-of-homozygosity interval table is patient-derived and is "
            "named in ETHICS.md as one of the two items that reached a model")

    def test_the_manual_items_are_not_silently_omitted(self, mod):
        joined = " ".join(mod.MANUAL).lower()
        for phrase in ("synapse.org", "transcript", "compute node"):
            assert phrase in joined, (
                f"the manual obligation mentioning {phrase!r} is missing, which "
                f"would imply the obligation is discharged when it is not")

    def test_the_ethics_document_still_describes_this_obligation(self):
        text = (REPO / "ETHICS.md").read_text()
        assert "3b" in text and "delet" in text.lower()
        assert "MVAHackathon2026@synapse.org" in text
