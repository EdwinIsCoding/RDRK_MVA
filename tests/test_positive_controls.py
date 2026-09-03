"""Positive-control spike-in recall, per plan sections 5.2 and 5.3.

Read this before reading the numbers it produces
------------------------------------------------
Recall here is measured against a pipeline that is **not yet complete**. VEP,
gnomAD and SpliceAI are unavailable on the recon host (see
``mva.track1.pipeline``), so the only evidence in play is gene plausibility,
region class and quality. A recall figure from this configuration therefore
measures **gene-level triage**, not variant prioritisation, and the tests below
say so in their assertions rather than leaving it to a footnote.

A second and more serious limit, established in
``scripts/13_harvest_splice_controls.py``: of the 108 confidently pathogenic
MVA-gene variants in ClinVar, **none is deep intronic, near-splice, synonymous
or UTR**. The benchmark cannot, on its own, test the hypothesis class this
project is built around. ``benchmarks/splice_mechanism_controls.tsv`` exists to
partially cover that gap, and is explicitly not MVA-specific.

Neither limitation is a reason to skip the harness. Both are reasons to report
its output with the qualification attached.
"""

from __future__ import annotations

import collections

import pytest

from mva.evidence import VariantClass
from mva.track1.pipeline import Track1Pipeline
from mva.track1.spikein import load_benchmark, spike

# ---------------------------------------------------------------------------
# Benchmark integrity. These run anywhere and guard the ground truth itself.
# ---------------------------------------------------------------------------

class TestBenchmarkIntegrity:
    def test_benchmark_exists_and_is_non_trivial(self):
        rows = load_benchmark(tiers=())
        if not rows:
            pytest.skip("benchmark not built; run scripts/10_harvest_clinvar_benchmark.py")
        assert len(rows) >= 100

    def test_every_row_has_a_resolvable_clinvar_accession(self):
        rows = load_benchmark(tiers=())
        if not rows:
            pytest.skip("benchmark not built")
        bad = [r for r in rows if not r.clinvar_vcv.startswith("VCV")]
        assert not bad, f"{len(bad)} rows lack a VCV accession, e.g. {bad[:3]}"

    def test_positions_are_grch38_and_no_chr(self):
        # The proband callset is GRCh38 with Ensembl naming. A benchmark in a
        # different convention would spike into coordinates that do not exist.
        for r in load_benchmark(tiers=())[:200]:
            assert r.position.build == "GRCh38"
            assert not r.position.contig.startswith("chr")

    def test_tier1_is_confidently_pathogenic_only(self):
        for r in load_benchmark(tiers=(1,)):
            assert r.clinical_significance.lower() in {
                "pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"
            }, f"{r.label} is tier 1 but classified {r.clinical_significance!r}"

    def test_documents_the_missing_variant_classes(self):
        """The absence of cryptic-splice positives is a property of the
        benchmark that must stay visible. If this test starts failing because
        deep-intronic tier-1 rows appeared, that is good news and the
        limitation text in RECON.md and the submission should be updated."""
        tier1 = load_benchmark(tiers=(1,))
        if not tier1:
            pytest.skip("benchmark not built")
        classes = {r.variant_class for r in tier1}
        hard = {VariantClass.DEEP_INTRONIC, VariantClass.SPLICE_REGION,
                VariantClass.SYNONYMOUS, VariantClass.UTR_PROMOTER}
        present = classes & hard
        assert not present, (
            f"benchmark now contains tier-1 {[c.value for c in present]} positives. "
            "Update the stated limitation in RECON.md and in the submission: the "
            "benchmark can now test the cryptic-allele hypothesis directly."
        )


# ---------------------------------------------------------------------------
# Spike-in mechanics.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def spiked(tmp_path_factory, background_vcf):
    variants = load_benchmark(tiers=(1,))
    if not variants:
        pytest.skip("benchmark not built")
    out = tmp_path_factory.mktemp("spike") / "spiked.vcf.gz"
    spike(background_vcf, variants, out)
    return out, variants


class TestSpikeIn:
    def test_spiked_variants_are_present_and_findable(self, spiked):
        import subprocess
        out, variants = spiked
        found = subprocess.run(
            ["bcftools", "query", "-i", "INFO/SPIKED=1", "-f", "%INFO/SPIKEGENE\n", str(out)],
            capture_output=True, text=True, check=True).stdout.split()
        assert len(found) >= len(variants) * 0.9, (
            f"only {len(found)}/{len(variants)} spikes survived the merge"
        )

    def test_background_is_not_the_proband(self, background_vcf):
        # Spiking into the proband would contaminate the thing under
        # investigation and make every recall figure meaningless.
        import subprocess
        samples = subprocess.run(["bcftools", "query", "-l", str(background_vcf)],
                                 capture_output=True, text=True, check=True).stdout.split()
        assert "WGS_EX2312012" not in samples, "background VCF contains the proband"
        assert samples, "background VCF has no samples"


# ---------------------------------------------------------------------------
# Recall. Plan section 5.2.
# ---------------------------------------------------------------------------

def _recall_at(ranked_genes: list[str], gene: str, n: int) -> bool:
    return gene in ranked_genes[:n]


@pytest.mark.positive_control
@pytest.mark.slow
class TestRecall:
    def test_gene_level_triage_recovers_spiked_genes(self, spiked, gene_model, panel_tsv):
        """With annotation stubs unavailable this measures gene-level triage.

        The assertion is deliberately weak, because a strong assertion here
        would be measuring something the pipeline is not yet doing. It becomes
        meaningful once VEP, gnomAD and SpliceAI are wired in on the GPU host.
        """
        out, variants = spiked
        pipeline = Track1Pipeline(gene_model, panel_tsv=panel_tsv)

        # Restrict to the spiked genes' regions: a whole-genome pass over the
        # background is not what this test is measuring and costs minutes.
        regions = sorted({f"{v.position.contig}" for v in variants})
        result = pipeline.run(out, regions=regions)

        assert result.n_records_in_panel > 0, "no panel variants reached scoring"
        ranked = result.genes_ranked()
        spiked_genes = {v.gene for v in variants}
        recovered = spiked_genes & set(ranked)
        assert recovered, (
            f"no spiked gene appeared in the ranked list at all. "
            f"spiked={sorted(spiked_genes)} ranked_top10={ranked[:10]}"
        )

    def test_reports_recall_broken_down_by_variant_class(self, spiked, gene_model,
                                                         panel_tsv, capsys):
        """Plan section 5.2 requires recall per class, not in aggregate.

        A pipeline that recovers coding loss of function and misses every
        deep-intronic positive has a specific, actionable weakness that an
        aggregate figure hides.
        """
        out, variants = spiked
        pipeline = Track1Pipeline(gene_model, panel_tsv=panel_tsv)
        result = pipeline.run(out, regions=sorted({v.position.contig for v in variants}))
        ranked = result.genes_ranked()

        by_class: dict[str, list[bool]] = collections.defaultdict(list)
        for v in variants:
            by_class[v.variant_class.value].append(_recall_at(ranked, v.gene, 20))

        with capsys.disabled():
            print("\n\n  recall@20 by variant class")
            print(f"  {'class':26} {'n':>4} {'recall':>8}")
            for cls, hits in sorted(by_class.items()):
                print(f"  {cls:26} {len(hits):>4} {sum(hits) / len(hits):>8.2f}")
            print(f"\n  {result.completeness_note()}")
            print("\n  Read these numbers narrowly. Candidates are restricted to the")
            print("  mitotic panel and every spiked variant lies in a panel gene, so a")
            print("  recall of 1.00 shows the gene reaches the ranked list, not that the")
            print("  causal variant is prioritised within it. The figure becomes")
            print("  meaningful only once VEP, gnomAD and SpliceAI are wired in.\n")

        assert by_class, "no spiked variant was scored at all"
        # The caveat must travel with the number. A recall table published
        # without it would be the most misleading artefact in the submission.
        assert not result.is_complete, (
            "pipeline now reports a complete run; update the qualifying text "
            "above and in STOP2_STATUS.md, and re-read the recall figures"
        )
        assert "spliceai" in result.missing_evidence_classes
