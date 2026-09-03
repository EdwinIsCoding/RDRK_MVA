"""The published availability resource must agree with itself and with the report.

It is the one artefact of this project intended for reuse by other people, so a
row count that disagrees with its own summary, or a licence note that goes
missing, is worse here than in an internal file.
"""
from __future__ import annotations

import json
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
RES = REPO / "resources" / "directional_availability"


@pytest.fixture(scope="module")
def summary():
    f = RES / "summary.json"
    if not f.exists():
        pytest.skip("run scripts/30_publish_directional_availability.py first")
    return json.loads(f.read_text())


def _genes(name: str) -> list[str]:
    lines = (RES / name).read_text().splitlines()
    return [line.split("\t")[0] for line in lines[1:] if line.strip()]


class TestTheTablesMatchTheirSummary:
    def test_activatable_row_count(self, summary):
        assert len(_genes("activatable_genes.tsv")) == summary["activatable_genes"]

    def test_inhibitable_row_count(self, summary):
        assert len(_genes("inhibitable_genes.tsv")) == summary["inhibitable_genes"]

    def test_both_directions_is_the_actual_intersection(self, summary):
        a, i = set(_genes("activatable_genes.tsv")), set(_genes("inhibitable_genes.tsv"))
        both = {line for line in
                (RES / "both_directions_genes.tsv").read_text().splitlines()[1:]
                if line.strip()}
        assert both == a & i
        assert len(both) == summary["genes_with_both_directions"]

    def test_rates_are_consistent_with_the_counts(self, summary):
        n = summary["approved_protein_coding_genes_hgnc"]
        assert summary["activation_rate"] == pytest.approx(
            summary["activatable_genes"] / n, abs=1e-5)
        assert summary["inhibition_rate"] == pytest.approx(
            summary["inhibitable_genes"] / n, abs=1e-5)

    def test_inhibition_really_is_the_more_available_direction(self, summary):
        """The asymmetry this resource exists to make checkable. If it ever
        inverts, the Track 2 report's section 3.3 argument needs rewriting."""
        assert summary["inhibitable_genes"] > summary["activatable_genes"]


class TestGenesAreWellFormed:
    def test_no_empty_or_duplicated_symbols(self):
        for name in ("activatable_genes.tsv", "inhibitable_genes.tsv"):
            g = _genes(name)
            assert all(s.strip() for s in g), f"empty gene symbol in {name}"
            assert len(g) == len(set(g)), f"duplicate gene symbol in {name}"

    def test_the_closed_axis_targets_are_absent_from_the_activatable_table(self):
        """The Track 2 claim, checkable in one line by anyone."""
        a = set(_genes("activatable_genes.tsv"))
        closed = {"PLK1", "AURKA", "AURKB", "CDK1", "TTK", "CENPE",
                  "KNL1", "MAD2L1BP", "ATM"}
        assert not (a & closed), (
            f"{sorted(a & closed)} now has an activating drug. The Track 2 "
            f"report's central negative no longer holds as stated.")

    def test_the_wrong_direction_is_well_supplied(self):
        i = set(_genes("inhibitable_genes.tsv"))
        assert {"PLK1", "AURKA", "AURKB", "CDK1", "TTK"} <= i, (
            "the report states that six of the ten closed-axis targets have "
            "inhibitors; that is the contrast the negative depends on")


class TestAttributionSurvives:
    def test_readme_carries_chembl_attribution_and_licence(self):
        f = RES / "README.md"
        if not f.exists():
            pytest.skip("resource not generated")
        text = f.read_text()
        assert "ChEMBL" in text
        assert "CC BY-SA 3.0" in text, "the ChEMBL licence must be named"
        assert "CC BY 4.0" in text, "our own licence must be named"
        assert "ebi.ac.uk/chembl" in text, "attribution needs a resolvable link"

    def test_the_limits_section_is_present(self):
        """A reusable table without its caveats is the thing we were trying not
        to produce."""
        text = (RES / "README.md").read_text().lower()
        for phrase in ("availability is not suitability", "snapshot"):
            assert phrase in text, f"missing caveat: {phrase!r}"
