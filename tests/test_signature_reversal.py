"""Guards for the signature-reversal analysis.

Plan section 7.3 calls this the highest-yield repurposing method and permits a
published proxy where no patient transcriptome exists, on one condition: that
the proxy is labelled clearly. Every test here defends that condition or the
null control that makes the result mean anything.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "results" / "summaries" / "track2_signature_reversal.md"
SCRIPT = REPO / "scripts" / "41_signature_reversal.py"
T2 = REPO / "submission" / "track2_nexusdwin_report.md"


@pytest.fixture(scope="module")
def out():
    if not OUT.exists():
        pytest.skip("run scripts/41_signature_reversal.py first")
    return OUT.read_text()


class TestTheProxyIsLabelled:
    """The plan's one condition. A proxy presented as the patient would be the
    single worst thing this analysis could do."""

    def test_every_difference_from_the_patient_is_stated(self, out):
        low = out.lower()
        for axis in ("mouse", "heart", "hypomorph"):
            assert axis in low, f"the proxy table no longer names {axis}"
        assert "not this patient" in low or "is not this patient" in low or \
               "every way it is not" in low, "the proxy is not labelled as a proxy"

    def test_the_report_carries_the_same_caveats(self):
        text = T2.read_text()
        if "5a" not in text:
            pytest.skip("section not in the report")
        low = text.lower()
        assert "mouse" in low and "heart" in low
        assert "direction to look, not a candidate to advance" in low, (
            "the report must not present reversal hits as candidates")

    def test_the_geo_accession_is_cited(self, out):
        assert re.search(r"GSE\d+", out), "the proxy dataset must be citable"


class TestOrthologyIsSourced:
    def test_mgi_homology_not_the_naming_convention(self):
        src = SCRIPT.read_text()
        assert "HOM_MouseHumanSequence" in src or "MGI" in src
        assert "uppercase naming convention" in src, (
            "the script must say why it does not simply uppercase mouse symbols")

    def test_symbols_are_validated_against_hgnc(self):
        assert "approved_symbols" in SCRIPT.read_text()


class TestTheNullControlExists:
    """Without it, a list of HDAC and mTOR inhibitors is indistinguishable from
    a list of HDAC and mTOR inhibitors."""

    def test_the_control_is_run_and_seeded(self):
        src = SCRIPT.read_text()
        assert "def null_control(" in src
        assert "SEED" in src, "the random draws must be reproducible"

    def test_the_summary_reports_the_overlap(self, out):
        assert "random genes" in out.lower()
        m = re.search(r"\*\*On average ([\d.]+) of (\d+) perturbagens, (\d+)%", out)
        assert m, "the null-control overlap is no longer reported"

    def test_a_high_overlap_would_be_called_out(self):
        """If the hits ever become mostly generic, the summary must say so
        rather than printing a table that looks specific."""
        src = SCRIPT.read_text()
        # Matched on the branch structure rather than on a prose string, which
        # a line break in the source would otherwise defeat.
        assert "frac >= 0.5" in src, (
            "the branch that disowns a mostly generic result was removed")
        assert "frac >= 0.2" in src, "the partial-overlap branch was removed"
        # Join adjacent Python string literals before matching, so a phrase
        # split across a line continuation still reads as one sentence.
        joined = re.sub(r'"\s*"', "", " ".join(src.split()))
        assert "largely not about BubR1" in joined, (
            "the high-overlap branch no longer disowns the result")


class TestTheSafetyScreenStillApplies:
    def test_reversal_hits_go_through_the_screen(self):
        src = SCRIPT.read_text()
        assert "from mva.track2.safety import" in src
        assert "screen(rec)" in src

    def test_a_cytotoxic_hit_is_excluded_not_listed(self, out):
        """A connectivity method will rank chemotherapy first for a child with a
        cancer predisposition syndrome, because it optimises transcriptional
        opposition and knows nothing about the patient."""
        if "DAUNORUBICIN" not in out.upper():
            pytest.skip("daunorubicin not in this run's hits")
        for line in out.splitlines():
            if "DAUNORUBICIN" in line.upper() and line.strip().startswith("|"):
                assert "excluded" in line.lower(), (
                    "daunorubicin is a cytotoxic anthracycline and must be "
                    f"excluded by the screen, not listed:\n  {line.strip()}")

    def test_no_dose_reaches_the_output(self, out):
        """CLAUDE.md rule 3. LINCS records some perturbagens by supplier
        catalogue number, and a code such as 656402-250mg contains a string that
        looks like a dose without being one, so backtick-quoted identifiers are
        excluded from the check and labelled as codes in the output itself."""
        dose = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|IU)\b", re.IGNORECASE)
        bad = []
        for line in out.splitlines():
            stripped = re.sub(r"`[^`]*`", "", line)
            if dose.search(stripped):
                bad.append(line.strip())
        assert not bad, f"a dose reached the signature output: {bad[:3]}"
        assert "not a dose" in out, (
            "the catalogue codes must be labelled so a reader does not read one "
            "as a dose")
