"""The Track 2 report must not drift from the data that produced it.

Every figure in submission/track2_nexusdwin_report.md was transcribed by hand
from a generated summary under results/summaries/. Hand transcription is exactly
how a wrong number reaches a deliverable, and this project has already shipped
several that way. These tests parse both and compare them.

They skip rather than fail when the generated summaries are absent, because
results/ is gitignored and a fresh clone will not have them until `make track2`
has run.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
REPORT = REPO / "submission" / "track2_nexusdwin_report.md"
PITCH = REPO / "submission" / "track2_nexusdwin_pitch.md"
CHEMOPREV = REPO / "results" / "summaries" / "track2_chemoprevention.md"
AVAIL = REPO / "results" / "summaries" / "track2_axis_availability.md"
DIRECTION = REPO / "results" / "summaries" / "track2_direction_audit.md"


def _rows(path: pathlib.Path) -> list[list[str]]:
    out = []
    for line in path.read_text().splitlines():
        if line.strip().startswith("|") and "---" not in line:
            out.append([c.strip() for c in line.strip().strip("|").split("|")])
    return out


@pytest.fixture(scope="module")
def report() -> str:
    if not REPORT.exists():
        pytest.skip("Track 2 report not written")
    return REPORT.read_text()


class TestCandidateTableMatchesTheGeneratedSummary:
    @pytest.fixture(scope="class")
    def generated(self):
        if not CHEMOPREV.exists():
            pytest.skip("run scripts/27_track2_chemoprevention.py first")
        by_agent = {}
        for r in _rows(CHEMOPREV):
            # candidate table: agent | ChEMBL | max phase | ATC | trials | paeds | verdict
            if len(r) == 7 and r[1].startswith("CHEMBL") and r[6].startswith("**"):
                by_agent[r[0].upper()] = {
                    "chembl": r[1],
                    "trials": int(r[4]),
                    "verdict": r[6].strip("*"),
                }
        assert by_agent, "no candidate rows parsed from the generated summary"
        return by_agent

    @pytest.fixture(scope="class")
    def claimed(self, generated):
        if not REPORT.exists():
            pytest.skip("Track 2 report not written")
        rows = {}
        for r in _rows(REPORT):
            # report table: Agent | ChEMBL | trials | endpoint | verdict
            if len(r) == 5 and r[1].startswith("CHEMBL"):
                name = r[0].replace("*", "").strip().upper()
                rows[name] = {"chembl": r[1],
                              "trials": int(r[2].replace("*", "").strip()),
                              "verdict": r[4].replace("*", "").strip()}
        assert rows, "no candidate rows parsed from the report"
        return rows

    def test_every_reported_agent_exists_in_the_generated_data(self, claimed, generated):
        missing = sorted(set(claimed) - set(generated))
        assert not missing, (
            f"the report names agents the pipeline did not produce: {missing}. "
            f"An agent that is not in the generated data is an invented one.")

    @pytest.mark.parametrize("field", ["chembl", "trials", "verdict"])
    def test_reported_values_match(self, claimed, generated, field):
        bad = {a: (v[field], generated[a][field]) for a, v in claimed.items()
               if a in generated and v[field] != generated[a][field]}
        assert not bad, (
            f"the report disagrees with the generated summary on {field}. "
            f"agent -> (report, data): {bad}")


class TestHeadlineNumbersMatch:
    def _num(self, text: str, pattern: str) -> str | None:
        m = re.search(pattern, text)
        return m.group(1) if m else None

    def test_base_rate_figures_match(self, report):
        if not AVAIL.exists():
            pytest.skip("run scripts/28_track2_axis_availability.py first")
        data = AVAIL.read_text()
        for label in ("359", "1,541", "19,297", "1.86%", "7.99%"):
            assert label in data, f"{label} is no longer in the generated data"
            assert label in report, (
                f"the report quotes a base-rate figure the data no longer has, "
                f"or has dropped {label}")

    def test_direct_axis_figures_match(self, report):
        if not DIRECTION.exists():
            pytest.skip("run scripts/14_track2_direction_audit.py first")
        data = DIRECTION.read_text()
        n_mech = self._num(data, r"\*\*(\d+) drug-mechanism records")
        assert n_mech, "could not parse the mechanism count from the audit"
        assert n_mech in report, (
            f"the audit reports {n_mech} drug-mechanism records and the report "
            f"does not quote that number")
        assert "0 act in the required direction" in data
        assert "**0**" in report or "none acts in that direction" in report

    def test_the_empty_evidence_base_is_still_empty(self, report):
        if not CHEMOPREV.exists():
            pytest.skip("run scripts/27_track2_chemoprevention.py first")
        data = CHEMOPREV.read_text()
        counts = re.findall(r"\| any trial in mosaic variegated aneuploidy \| \*\*(\d+)\*\*", data)
        assert counts == ["0"], (
            "a trial in mosaic variegated aneuploidy now exists in the registry. "
            "The report's central premise about the chemoprevention axis has "
            "changed and section 4.1 must be rewritten, not patched.")


class TestPitchScript:
    """The pitch is a submission artefact and is bound by the same rules.

    It also has a hard constraint the report does not: the video must be three
    minutes, and a script that overruns is a script that gets cut mid-argument.
    """

    @pytest.fixture(scope="class")
    def pitch(self):
        if not PITCH.exists():
            pytest.skip("pitch script not written")
        return PITCH.read_text()

    def spoken_words(self, pitch: str) -> int:
        return sum(len(l[1:].split()) for l in pitch.splitlines()
                   if l.startswith("> "))

    def test_fits_in_three_minutes(self, pitch):
        w = self.spoken_words(pitch)
        assert w > 200, f"only {w} spoken words; the script is a stub"
        assert w <= 420, (
            f"{w} spoken words is {round(w / 140 * 60)}s at 140 wpm, over the "
            f"three-minute limit. The video will be cut mid-argument.")

    def test_the_stated_word_count_is_true(self, pitch):
        """A script that misstates its own length is the kind of unchecked
        number this project keeps finding in its own documents."""
        m = re.search(r"\*\*Spoken length:\*\* (\d+) words", pitch)
        assert m, "the pitch no longer states its spoken length"
        assert int(m.group(1)) == self.spoken_words(pitch), (
            f"the pitch claims {m.group(1)} spoken words and actually has "
            f"{self.spoken_words(pitch)}")

    def test_no_dose_and_no_em_dash(self, pitch):
        assert not TestReportObeysTheHardRules.DOSE.search(pitch), \
            "a dose reached the pitch script"
        assert "\u2014" not in pitch, "an em dash reached the pitch script"

    def test_it_forbids_showing_patient_data_on_screen(self, pitch):
        """The production notes are the only place this constraint can live, and
        a video is the easiest way to leak a pileup screenshot."""
        low = pitch.lower()
        assert "do not show patient data" in low
        assert "vcf" in low and "bam" in low


class TestReportObeysTheHardRules:
    DOSE = re.compile(
        r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|µg|ug|ml|mL|IU|iu|units?)\b"
        r"|\b\d+(?:\.\d+)?\s*(?:mg|mcg|g)\s*/\s*(?:kg|m2|day)\b", re.IGNORECASE)

    def test_no_dose_appears(self, report):
        """CLAUDE.md rule 3."""
        bad = [l.strip() for l in report.splitlines() if self.DOSE.search(l)]
        assert not bad, "a dose reached the Track 2 report:\n  " + "\n  ".join(bad[:5])

    def test_no_em_dashes(self, report):
        """Project style, CLAUDE.md rule 7."""
        assert "—" not in report, "an em dash reached the Track 2 report"

    def test_identifiers_are_well_formed(self, report):
        """A malformed identifier is usually an invented one."""
        for nct in set(re.findall(r"\bNCT\d+", report)):
            assert len(nct) == 11, f"{nct} is not a valid NCT identifier shape"
        for go in set(re.findall(r"\bGO:\d+", report)):
            assert len(go) == 10, f"{go} is not a valid GO identifier shape"
        for vcv in set(re.findall(r"\bVCV\d+", report)):
            assert len(vcv) == 12, f"{vcv} is not a valid ClinVar accession shape"

    def test_nct_identifiers_are_traceable_to_the_generated_data(self, report):
        if not CHEMOPREV.exists():
            pytest.skip("run scripts/27_track2_chemoprevention.py first")
        data = CHEMOPREV.read_text()
        claimed = set(re.findall(r"\bNCT\d{8}\b", report))
        produced = set(re.findall(r"\bNCT\d{8}\b", data))
        invented = sorted(claimed - produced)
        assert not invented, (
            f"the report cites trial identifiers the pipeline never returned: "
            f"{invented}. CLAUDE.md rule 2.")
