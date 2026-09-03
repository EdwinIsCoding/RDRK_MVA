"""The Track 1 report must not drift from the callset or from its own tables.

Its numbers were verified by hand in September 2026 (docs/VERIFICATION.md) and
nine of them were wrong at the time. Hand verification does not repeat itself,
so the checks that can be automated are automated here: internal consistency,
agreement with the verification artefacts, and the identifier and style rules.

Checks needing the patient callset skip when the verification summaries are
absent, so this suite still runs on a clone with no data access.
"""
from __future__ import annotations

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
REPORT = REPO / "submission" / "track1_nexusdwin_report.md"
SUBMISSION = REPO / "submission" / "track1_submission.csv"
CALLSET = REPO / "results" / "summaries" / "verification_callset.md"
READLEVEL = REPO / "results" / "summaries" / "verification_readlevel.md"

#: The answer, as submitted. Any change here is a change to the call itself.
ALLELE_1 = ("chr15", "40209701", "T", "G")
ALLELE_2 = ("chr15", "40220612", "T", "G")
SECONDARY = ("chr22", "20996720", "C", "G")


@pytest.fixture(scope="module")
def report():
    if not REPORT.exists():
        pytest.skip("Track 1 report absent")
    return REPORT.read_text()


class TestTheCallIsInternallyConsistent:
    def test_the_submitted_csv_matches_the_report(self, report):
        if not SUBMISSION.exists():
            pytest.skip("submission absent")
        import csv
        rows = list(csv.DictReader(SUBMISSION.open()))
        primary = [r for r in rows if r["finding_type"] == "primary"]
        assert len(primary) == 1
        p = primary[0]
        assert (p["chrom_1"], p["pos_1"], p["ref_1"], p["alt_1"]) == ALLELE_1
        assert (p["chrom_2"], p["pos_2"], p["ref_2"], p["alt_2"]) == ALLELE_2
        for chrom, pos, ref, alt in (ALLELE_1, ALLELE_2):
            assert f"{chrom}:{pos} {ref}>{alt}" in report or pos in report, (
                f"{chrom}:{pos} is in the submission but not in the report")

    def test_the_allele_separation_arithmetic(self, report):
        stated = re.findall(r"10,?911", report)
        assert stated, "the report no longer states the allele separation"
        assert int(ALLELE_2[1]) - int(ALLELE_1[1]) == 10911, (
            "the coordinates no longer differ by the separation the report claims")

    def test_read_level_strand_counts_sum(self, report):
        """The published table once gave 15 fwd + 12 rev beside 26 alt reads."""
        alt, strand = [], []
        for line in report.splitlines():
            if line.strip().startswith("| Ref / alt reads"):
                for cell in line.split("|")[2:]:
                    m = re.match(r"\s*(\d+)\s*/\s*(\d+)\s*$", cell)
                    if m:
                        alt.append(int(m.group(2)))
            if "fwd" in line and "rev" in line and line.strip().startswith("|"):
                for cell in line.split("|")[2:]:
                    m = re.match(r"\s*(\d+)\s*fwd\s*/\s*(\d+)\s*rev\s*$", cell)
                    if m:
                        strand.append((int(m.group(1)), int(m.group(2))))
        if not alt or not strand:
            pytest.skip("read-level table not in its expected shape")
        assert len(alt) == len(strand)
        for a, (f, r) in zip(alt, strand, strict=True):
            assert f + r == a, f"{f} fwd + {r} rev is not {a} alternate reads"


class TestTheReportAgreesWithTheVerificationArtefacts:
    def test_genotypes_match_the_callset(self, report):
        if not CALLSET.exists():
            pytest.skip("run scripts/25_verify_track1_claims.py first")
        data = CALLSET.read_text()
        for pos in (ALLELE_1[1], ALLELE_2[1], SECONDARY[1]):
            assert pos in data, f"{pos} is no longer in the verified callset table"
        assert "| agrees | YES |" in data, (
            "the callset verification no longer agrees with the submission")

    def test_record_counts_match(self, report):
        """Only the counts the report actually quotes. It need not repeat every
        verified figure, but any it does repeat must still be true."""
        if not CALLSET.exists():
            pytest.skip("verification summary absent")
        data = CALLSET.read_text()
        for n in ("5,012,204", "4,950,283", "61,921"):
            assert n in data, (
                f"{n} is no longer produced by the verification script, so any "
                f"document still quoting it is stale")
        quoted = [n for n in ("5,012,204", "4,950,283", "61,921") if n in report]
        assert quoted, "the report quotes none of the verified record counts"

    def test_read_level_figures_match(self, report):
        if not READLEVEL.exists():
            pytest.skip("run scripts/26_verify_readlevel.py first")
        data = READLEVEL.read_text()
        for vaf in ("0.553", "0.448"):
            assert vaf in data, f"VAF {vaf} is no longer in the recomputation"
            assert vaf in report, f"the report does not carry VAF {vaf}"


class TestTheCorrectedClaimsStayCorrected:
    """One test per finding in docs/VERIFICATION.md that reached a deliverable."""

    def test_allele_2_is_not_called_absent_from_gnomad(self, report):
        for line in report.splitlines():
            if "40220612" in line and "absent" in line.lower():
                assert "genomes" in line.lower(), (
                    "allele 2 is described as absent from gnomAD. It is present "
                    f"in v4.1 exomes at 6.84e-07.\n  {line.strip()}")

    def test_the_protein_change_clinvar_record_is_still_cited(self, report):
        assert "VCV004600147" in report, (
            "the ClinVar record of uncertain significance for the same protein "
            "change was removed. It cuts against the missense and belongs here.")

    def test_popmax_is_not_mislabelled(self, report):
        assert "popmax" not in report.lower(), (
            "the report reintroduced 'popmax'. The figures quoted are global "
            "allele frequencies and group max separately; see VERIFICATION 2.3.")

    def test_no_claim_of_zero_phasing_groups_in_bub1b(self, report):
        bad = re.compile(r"no\s+(?:`?PGT`?/`?PID`?\s+)?phasing group\s+"
                         r"(?:exists\s+)?(?:anywhere\s+)?in\s+`?BUB1B`?", re.I)
        assert not bad.search(report), (
            "one phasing group does exist in BUB1B (40216568_TA_T). The "
            "defensible claim is that neither causal allele carries a tag.")

    def test_pex5_and_ctu2_are_described_as_closed(self, report):
        assert "7190512" in report and "88714226" in report
        assert "0.252" in report and "0.767" in report, (
            "the gnomAD frequencies that closed PEX5 and CTU2 were removed")
        assert "open rather than settled" not in report

    def test_the_shortlist_provenance_claim_is_accurate(self, report):
        assert "849bf98" in report, "the shortlist provenance commit is no longer cited"
        assert "not in git" in report or "nothing under `results/` is tracked" in report, (
            "the report claims an artefact is timestamped in the repository. "
            "Nothing under results/ is tracked; only the code and its commit "
            "message are.")


class TestTheReportObeysTheHardRules:
    def test_no_em_dash(self, report):
        assert "—" not in report

    def test_identifier_shapes(self, report):
        for acc in set(re.findall(r"\bVCV\d+", report)):
            assert len(acc) == 12, f"{acc} is not a valid ClinVar accession shape"
        for rs in set(re.findall(r"\brs\d+", report)):
            assert rs[2:].isdigit()

    def test_no_dosing_language(self, report):
        """CLAUDE.md rule 3, which binds both tracks."""
        dose = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mg|mcg|IU)\b", re.I)
        hits = [line.strip() for line in report.splitlines() if dose.search(line)]
        assert not hits, f"dosing language in the Track 1 report: {hits[:3]}"


class TestTheInSilicoPanelIsReportedHonestly:
    """The report once said PP3 rested on 'two concordant' predictors while
    listing AlphaMissense, CADD and REVEL as unavailable. They were available,
    and consulting them showed the panel is neither as narrow nor as concordant
    as either half of that sentence implied."""

    PANEL = REPO / "results" / "summaries" / "missense_predictor_panel.md"

    def test_pp3_no_longer_claims_two_concordant_predictors(self, report):
        low = report.lower()
        assert "two concordant in silico predictors" not in low, (
            "PP3 is again claiming two concordant predictors. Fifteen were "
            "consulted and five call the variant tolerated.")

    def test_the_report_states_the_disagreement(self, report):
        assert "tolerate" in report.lower(), (
            "the report no longer states that some predictors tolerate the "
            "variant, which is the honest half of this evidence")

    def test_alphamissense_is_no_longer_listed_as_unavailable(self, report):
        for line in report.splitlines():
            low = line.lower()
            if "alphamissense" in low and "not available" in low:
                pytest.fail(f"AlphaMissense is available and was consulted.\n  {line.strip()}")

    def test_the_counts_match_the_generated_panel(self, report):
        if not self.PANEL.exists():
            pytest.skip("run scripts/34_missense_predictor_panel.py first")
        text = self.PANEL.read_text()
        m = re.search(r"\*\*(\d+) calls? it damaging, (\d+) calls? it tolerated "
                      r"or benign, and (\d+) sits? in between\.\*\*", text)
        assert m, "could not parse the tally from the generated panel"
        # The missense allele is the second block; take the last match.
        tallies = re.findall(r"\*\*(\d+) calls? it damaging, (\d+) calls? it "
                             r"tolerated or benign, and (\d+) sits? in between\.\*\*", text)
        dmg, ben, mid = tallies[-1]
        total = int(dmg) + int(ben) + int(mid)
        assert f"{dmg} of {total}" in report or f"{dmg} damaging" in report, (
            f"the report does not carry the damaging count {dmg} of {total} that "
            f"the panel produced")
        assert f"{ben} of {total}" in report or f"{ben} tolerate" in report, (
            f"the report does not carry the tolerated count {ben}")


class TestProvenanceGapsStayClosed:
    def test_the_dataset_revision_is_recorded(self):
        text = (REPO / "PROVENANCE.md").read_text()
        assert "59e322d27f399006b398d366d33e703e48a29914" in text, (
            "the recovered HuggingFace dataset revision was removed")
        assert "SageBio/mva-hackathon-2026-data" in text

    def test_the_revision_is_labelled_as_inferred(self):
        """It was reconstructed from a last-modified date, not logged. Presenting
        it as a captured fact would be a stronger claim than the evidence."""
        text = (REPO / "PROVENANCE.md").read_text().lower()
        assert "inference, not a log entry" in text or "inferred" in text

    def test_the_ensembl_rest_coordinates_are_marked_superseded(self):
        text = (REPO / "config" / "gene_panels" / "panel_provenance.md").read_text()
        assert "Superseded" in text, (
            "panel_provenance.md again presents the REST coordinates as the "
            "source. The panel is built from the pinned Ensembl 115 GTF.")

    def test_bub3_carries_the_gtf_span_not_the_rest_span(self):
        """The concrete difference: 16,072 bp from the GTF against 158,779 bp
        from the live API."""
        panel = REPO / "config" / "gene_panels" / "mva_known.tsv"
        if not panel.exists():
            pytest.skip("panel absent")
        for line in panel.read_text().splitlines()[1:]:
            f = line.split("\t")
            if f and f[0] == "BUB3":
                span = int(f[7]) - int(f[6])
                assert span == 16072, (
                    f"BUB3 spans {span} bp. The pinned GTF gives 16,072; the "
                    f"live REST API gave 158,779 and was wrong.")
                return
        pytest.fail("BUB3 is not in the known-gene panel")
