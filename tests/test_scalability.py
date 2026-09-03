"""The generalisation check must be capable of failing.

A comparison that returns the same answer whatever it is given demonstrates
nothing, and one whose gene sets were typed in by hand demonstrates only that
the author remembered some genes. Both properties are asserted here.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "summaries" / "track2_scalability.md"
PANEL = REPO / "config" / "gene_panels" / "disease_genes.tsv"


@pytest.fixture(scope="module")
def report():
    if not OUT.exists():
        pytest.skip("run scripts/31_track2_scalability.py first")
    return OUT.read_text()


def _rows(text: str, header_fragment: str) -> list[list[str]]:
    out, seen = [], False
    for line in text.splitlines():
        if header_fragment in line:
            seen = True
            continue
        if seen:
            if not line.strip().startswith("|"):
                break
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if "---" not in line:
                out.append(cells)
    return out


class TestComplexesAreSeparatedFromGenes:
    """OmniPath names complexes by joining components with underscores. A
    complex can never match a gene symbol in the availability table, so leaving
    it in counts it as 'no drug available' and biases every result towards our
    own conclusion."""

    def test_no_complex_identifier_reaches_a_gene_list(self, report):
        body = report.split("## 2.")[1].split("## 3.")[0]
        offenders = sorted({g for g in re.findall(r"`([A-Z0-9]+(?:_[A-Z0-9]+)+)`", body)})
        assert not offenders, (
            f"complex identifiers reached the per-gene lists: {offenders[:5]}")

    def test_complexes_are_counted_and_reported(self, report):
        assert "Complexes, set aside" in report, (
            "the count of discarded complexes must be shown, not silently dropped")


class TestTheComparisonDiscriminates:
    """If every disease returns the same verdict, the method measures nothing."""

    def test_three_diseases_are_compared(self, report):
        rows = _rows(report, "| Disease | Seeds | Partners |")
        assert len(rows) >= 3, f"expected three diseases, parsed {len(rows)}"

    def test_the_diseases_do_not_all_return_the_same_answer(self, report):
        rows = _rows(report, "| Disease | Seeds | Partners |")
        available = {r[0]: r[-1].replace("*", "") for r in rows}
        assert len(set(available.values())) > 1, (
            f"every disease returned the same availability count: {available}. "
            f"The method is not discriminating and the scalability claim fails.")

    def test_mva_is_the_only_all_activation_disease(self, report):
        """The Track 2 report's surviving claim 2, which rested on one disease
        until this comparison existed."""
        rows = _rows(report, "| Disease | Seeds | Partners |")
        inhibit = {r[0]: int(r[6]) for r in rows if r[6].isdigit()}
        assert inhibit, "could not parse the inhibit column"
        mva = [k for k in inhibit if k.lower().startswith("mosaic")]
        assert mva, "the proband's disease is missing from the comparison"
        assert inhibit[mva[0]] == 0
        assert all(v > 0 for k, v in inhibit.items() if k != mva[0]), (
            "the comparators no longer have inhibition-reachable targets, so "
            "MVA is no longer distinctive and report section 3.4 needs revising")

    def test_the_registry_emptiness_is_specific_to_this_disease(self, report):
        """Section 4.1 of the report claims the empty evidence base is a
        property of the disease and not of the query. That only holds while the
        comparators return trials to the identical queries."""
        rows = _rows(report, "| Disease | Any interventional trial |")
        counts = {r[0]: int(r[1]) for r in rows if r[1].isdigit()}
        assert counts, "could not parse the registry table"
        mva = [k for k in counts if k.lower().startswith("mosaic")]
        assert counts[mva[0]] == 0
        assert all(v > 0 for k, v in counts.items() if k != mva[0]), (
            "the comparators now return no trials either, so the four zeros may "
            "be an artefact of the query rather than a fact about the disease")


class TestSeedSetsAreDerivedNotTyped:
    def test_the_comparator_genes_come_from_the_committed_panel(self, report):
        """This project has shipped hand-written panel identifiers before, of
        which six of eleven were wrong. The comparator gene sets must be
        recoverable from the panel file rather than from an author."""
        if not PANEL.exists():
            pytest.skip("panel not present")
        panel_symbols = {line.split("\t")[0]
                         for line in PANEL.read_text().splitlines()[1:] if line.strip()}
        section = report.split("## 1. Seed sets")[1].split("## 2.")[0]
        listed = set(re.findall(r"`([A-Z][A-Z0-9]+)`", section))
        mva = {"BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CEP192", "SMC5", "CENATAC"}
        unknown = sorted((listed - mva) - panel_symbols)
        assert not unknown, (
            f"these comparator genes are not in the curated panel: {unknown}")
