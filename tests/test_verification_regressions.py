"""Regressions for the errors found by the Track 1 verification pass.

Each test here corresponds to a numbered finding in docs/VERIFICATION.md. They
guard the deliverables rather than the library, because every one of these
errors reached a submitted document while the library underneath was correct.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent

DELIVERABLES = [
    REPO / "submission" / "track1_nexusdwin_report.md",
    REPO / "submission" / "track1_candidates.tsv",
    REPO / "submission" / "arm_c_readlevel_verification.md",
    REPO / "README.md",
    REPO / "config" / "config.yaml",
]


def _texts():
    return [(p, p.read_text()) for p in DELIVERABLES if p.exists()]


class TestAbsentFromGnomad:
    """docs/VERIFICATION.md section 2.1 and 2.4.

    Two separate deliverables described a variant as absent from gnomAD when the
    project's own lookup held a frequency for it, and two more described
    variants as absent when their chromosome was never assayed. "Absent from
    gnomAD" is the single most promoting piece of evidence a rare-disease filter
    has, so the phrase has to be earned.
    """

    #: Variants that are known NOT to be absent from gnomAD v4.1, with the
    #: group max allele frequency verified on 3 September 2026.
    NOT_ABSENT = {
        "40220612": 8.99e-07,   # BUB1B allele 2
        "7190512": 0.720,       # PEX5, in-frame 45 bp deletion
        "88714226": 0.925,      # CTU2, in-frame 6 bp deletion
    }

    @pytest.mark.parametrize("pos", sorted(NOT_ABSENT))
    def test_no_deliverable_calls_these_absent_from_gnomad(self, pos):
        for path, text in _texts():
            for line in text.splitlines():
                if pos not in line:
                    continue
                lowered = line.lower()
                if "absent" not in lowered:
                    continue
                # Permitted: an explicit correction, or the true statement that
                # BUB1B allele 2 is absent from gnomAD *genomes* specifically.
                if "genomes" in lowered or "not absent" in lowered:
                    continue
                if any(w in lowered for w in
                       ("earlier", "corrected", "wrong", "previously", "withdrawn")):
                    continue
                pytest.fail(
                    f"{path.relative_to(REPO)} describes position {pos} as absent "
                    f"from gnomAD. Its group max allele frequency is "
                    f"{self.NOT_ABSENT[pos]}. See docs/VERIFICATION.md.\n  {line.strip()}"
                )


class TestStrandCountsAddUp:
    """docs/VERIFICATION.md section 2.7.

    The published read-level table gave 15 forward and 12 reverse alternate
    reads beside an alternate count of 26, and 6 and 10 beside a count of 12.
    Neither sums. A reader checking one line of arithmetic found it.
    """

    ROW = re.compile(r"\|\s*(\d+)\s*fwd\s*/\s*(\d+)\s*rev\s*\|")
    ALT_ROW = re.compile(r"\|\s*Ref\s*/\s*alt reads\s*\|(.*)\|?\s*$", re.IGNORECASE)

    def test_read_level_strand_counts_sum_to_alt_counts(self):
        path = REPO / "submission" / "arm_c_readlevel_verification.md"
        text = path.read_text()

        alt_counts: list[int] = []
        strand_pairs: list[tuple[int, int]] = []
        for line in text.splitlines():
            if line.strip().startswith("| Ref / alt reads"):
                for cell in line.split("|")[2:]:
                    m = re.match(r"\s*(\d+)\s*/\s*(\d+)\s*$", cell)
                    if m:
                        alt_counts.append(int(m.group(2)))
            if "fwd" in line and "rev" in line and line.strip().startswith("|"):
                for cell in line.split("|")[2:]:
                    m = re.match(r"\s*(\d+)\s*fwd\s*/\s*(\d+)\s*rev\s*$", cell)
                    if m:
                        strand_pairs.append((int(m.group(1)), int(m.group(2))))

        assert alt_counts, "no alternate-read row found; the table shape changed"
        assert len(alt_counts) == len(strand_pairs), (
            f"{len(alt_counts)} alternate-read counts against "
            f"{len(strand_pairs)} strand splits in {path.name}"
        )
        for alt, (fwd, rev) in zip(alt_counts, strand_pairs):
            assert fwd + rev == alt, (
                f"strand split {fwd} fwd + {rev} rev = {fwd + rev}, but the same "
                f"column reports {alt} alternate reads. One of the two is wrong."
            )


class TestPositiveControlCounts:
    """docs/VERIFICATION.md section 2.6.

    ``validate()`` returned the number of SpliceAI annotation rows and the
    write-up called it a number of variants. Eight controls became "9/9".
    """

    def test_validate_reports_variants_and_annotations_separately(self):
        src = (REPO / "scripts" / "19_arm_b_splicing.py").read_text()
        body = src.split("def validate(")[1].split("\ndef ")[0]
        returns = [l.strip() for l in body.splitlines()
                   if l.strip().startswith("return (")]
        assert returns, "validate() no longer returns a tuple"
        for r in returns:
            assert r.count(",") >= 3, (
                "validate() must return annotation count and variant count "
                f"separately, so a run over N variants cannot be reported as "
                f"M/M variants. Got: {r}"
            )

    def test_no_deliverable_claims_nine_control_variants(self):
        for path, text in _texts():
            for line in text.splitlines():
                # Deliberately not requiring "splice" on the same line: the
                # original error lived on a line that only said "those controls
                # score 9/9", with the word splice two lines above it.
                if "9/9" in line:
                    pytest.fail(
                        f"{path.relative_to(REPO)} reports 9/9 control variants. "
                        f"There are eight controls yielding nine annotations, "
                        f"because 15:40218455 is annotated to both BUB1B and "
                        f"PAK6.\n  {line.strip()}"
                    )


class TestPhasingClaim:
    """docs/VERIFICATION.md section 2.5.

    A phasing group does exist in BUB1B (40216568_TA_T, two homozygous records).
    The claim that must survive is the narrower one: neither causal allele is in
    a phasing group.
    """

    def test_no_deliverable_claims_bub1b_has_no_phasing_group_at_all(self):
        bad = re.compile(
            r"no\s+(?:`?PGT`?/`?PID`?\s+)?phasing group\s+(?:exists\s+)?"
            r"(?:anywhere\s+)?in\s+`?BUB1B`?", re.IGNORECASE)
        for path, text in _texts():
            m = bad.search(text)
            if m:
                pytest.fail(
                    f"{path.relative_to(REPO)} claims no phasing group exists in "
                    f"BUB1B. One does: 40216568_TA_T. The defensible claim is "
                    f"that neither causal allele carries a phasing tag.\n"
                    f"  ...{text[max(0, m.start() - 60):m.end() + 60]}..."
                )
