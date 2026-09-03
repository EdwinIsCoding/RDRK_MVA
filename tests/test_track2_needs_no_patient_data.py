"""Track 2 must remain reproducible by anyone, with no data access.

This is a claim the README and the Track 2 report both make, and it is the
reason a judge can check our Track 2 results in an afternoon instead of applying
for the dataset. A single convenience read of a patient file would quietly
destroy it, so it is asserted rather than asserted-in-prose.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent

TRACK2_SCRIPTS = [
    "scripts/14_track2_direction_audit.py",
    "scripts/27_track2_chemoprevention.py",
    "scripts/28_track2_axis_availability.py",
    "scripts/29_track2_drift_check.py",
    "scripts/30_publish_directional_availability.py",
    "scripts/31_track2_scalability.py",
    "scripts/32_structural_feasibility.py",
]

#: Anything that would read the challenge distribution.
PATIENT_PATH = re.compile(r"""["'][^"']*\b(data/|WGS_EX|Challenge_Clinical|"""
                          r"""node_artefacts/)[^"']*["']""")


def _sources() -> list[pathlib.Path]:
    out = [REPO / s for s in TRACK2_SCRIPTS]
    out += sorted((REPO / "src" / "mva" / "track2").glob("*.py"))
    return [p for p in out if p.exists()]


class TestNoPatientDataInTheTrack2Pipeline:
    def test_every_track2_script_exists(self):
        missing = [s for s in TRACK2_SCRIPTS if not (REPO / s).exists()]
        assert not missing, (
            f"the pipeline moved and this test no longer covers it: {missing}")

    @pytest.mark.parametrize("name", TRACK2_SCRIPTS)
    def test_no_script_references_patient_data(self, name):
        p = REPO / name
        if not p.exists():
            pytest.skip(f"{name} absent")
        hits = [m.group(0) for m in PATIENT_PATH.finditer(p.read_text())]
        assert not hits, (
            f"{name} references patient data: {hits}. Track 2's independence "
            f"from the challenge distribution is what lets a judge reproduce it "
            f"without data access, and it is claimed in the README.")

    def test_no_track2_module_references_patient_data(self):
        bad = {}
        for p in sorted((REPO / "src" / "mva" / "track2").glob("*.py")):
            hits = [m.group(0) for m in PATIENT_PATH.finditer(p.read_text())]
            if hits:
                bad[p.name] = hits
        assert not bad, f"Track 2 modules reference patient data: {bad}"

    def test_the_only_required_download_is_declared(self):
        """If a Track 2 script starts needing a second reference file, the
        no-data path in the Makefile has to learn about it."""
        refs = set()
        for p in _sources():
            refs |= set(re.findall(r'["\'](refs/[^"\']+)["\']', p.read_text()))
        declared = {"refs/hgnc_complete_set.txt"}
        undeclared = sorted(refs - declared)
        assert not undeclared, (
            f"Track 2 now reads reference files the no-data path does not fetch: "
            f"{undeclared}. Add them to the downloads-track2 target.")


class TestTheClaimIsMadeWhereJudgesWillLookForIt:
    def test_readme_states_it(self):
        text = (REPO / "README.md").read_text().lower()
        assert "reproduce-track2" in text or "no patient data" in text, (
            "the README does not tell a reader Track 2 needs no data access")

    def test_the_makefile_offers_the_path(self):
        text = (REPO / "Makefile").read_text()
        assert "reproduce-track2:" in text
        assert "downloads-track2:" in text
