"""Track 1 pipeline: turn a VCF into a ranked candidate list.

Scope, stated plainly
---------------------
What runs today, on any machine:
  * gene assignment and region classification from the Ensembl gene model
  * splice-distance stratification (near-splice, deep intronic, beyond window)
  * gene plausibility from the panels, with constraint as a bonus only
  * quality penalties from FILTER, depth and strand bias
  * the additive score and ranking of plan section 6.7

What is stubbed and must be filled on the GPU host before any result is
reported as an Arm A or Arm B finding:
  * ``VepAnnotator``       protein consequence: missense, nonsense, frameshift
  * ``GnomadAnnotator``    population allele frequency
  * ``SpliceAiAnnotator``  SpliceAI and Pangolin delta scores

Each stub raises rather than returning a neutral value. A stub that quietly
returned "no evidence" would let the pipeline produce a plausible ranked list
that silently omits the evidence class the whole project depends on, and nobody
reading the output would be able to tell. Failing loudly is the point.
"""

from __future__ import annotations

import abc
import csv
import dataclasses
import pathlib
from typing import Iterable, Iterator

from mva.evidence import Candidate, Evidence, GenomicPosition, SourceType, VariantClass, Zygosity
from mva.track1 import scoring
from mva.track1.regions import GeneModel, iter_vcf


class NotAvailableHere(RuntimeError):
    """Raised by an annotator that needs a resource this machine does not have.

    Deliberately not a silent no-op. See the module docstring.
    """


# ---------------------------------------------------------------------------
# Annotator interface
# ---------------------------------------------------------------------------

class Annotator(abc.ABC):
    """Adds evidence to a candidate. Every implementation must be able to state
    its own version, because that version goes into the evidence trail."""

    name: str = "abstract"

    @property
    @abc.abstractmethod
    def version(self) -> str: ...

    @abc.abstractmethod
    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]: ...


class PanelAnnotator(Annotator):
    """Gene plausibility from the locally built panels. Runs anywhere."""

    name = "panel"

    def __init__(self, panel_tsv: str | pathlib.Path = "config/gene_panels/mitotic_extended.tsv"):
        self.path = pathlib.Path(panel_tsv)
        self.panel: dict[str, dict[str, str]] = {}
        if self.path.exists():
            with self.path.open(newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    self.panel[row["symbol"]] = row

    @property
    def version(self) -> str:
        return f"panel/{self.path.name}/{len(self.panel)}genes"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        return scoring.gene_plausibility(candidate.gene, self.panel)


class QualityAnnotator(Annotator):
    """Artefact indicators from the callset itself. Runs anywhere."""

    name = "quality"

    @property
    def version(self) -> str:
        return "quality/1.0.0"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        return scoring.quality_penalties(
            filter_field=record.get("filter", ""),
            depth=record.get("dp"),
            strand_bias_fs=record.get("fs"),
            in_segdup=bool(record.get("in_segdup")),
        )


class RegionAnnotator(Annotator):
    """Splice distance and region class. Runs anywhere given a GTF.

    This is the arm that matters most (RECON.md: hypothesis class 1), so its
    evidence records the measured distance rather than only the class, letting
    a reader see how deep 'deep intronic' actually was.
    """

    name = "region"

    def __init__(self, model: GeneModel, gtf_version: str = "Ensembl/115"):
        self.model = model
        self.gtf_version = gtf_version

    @property
    def version(self) -> str:
        return f"regions/{self.gtf_version}/{len(self.model)}genes"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        dist = record.get("splice_distance")
        if dist is None or candidate.variant_class is VariantClass.UNCLASSIFIED:
            return []
        # No weight of its own: the region class conditions the splicing
        # prediction rather than being evidence for pathogenicity by itself.
        # An intron is not suspicious; a predicted cryptic site in one is.
        return [Evidence(
            criterion="region_class",
            statement=(f"{candidate.variant_class.value}, {dist} bp from the nearest "
                       f"exon boundary of {candidate.gene}"),
            source_type=SourceType.TOOL,
            source_id=f"ensembl-gtf/{self.gtf_version}",
            weight=0.0,
            detail={"splice_distance_bp": dist},
        )]


class VepAnnotator(Annotator):
    """Protein consequence. Needs a VEP cache of roughly 25 GB."""

    name = "vep"

    @property
    def version(self) -> str:
        return "ensembl-vep/116.1"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        raise NotAvailableHere(
            "VEP is not installed and no cache is present. Install per "
            "environment.yml on the GPU host. Until then, protein consequence "
            "classes (missense, nonsense, frameshift) are unavailable and Arm A "
            "results must not be reported as complete."
        )


class GnomadAnnotator(Annotator):
    """Population frequency. Needs the gnomAD v4.1 sites VCF."""

    name = "gnomad"

    @property
    def version(self) -> str:
        return "gnomAD v4.1 joint"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        af = record.get("gnomad_popmax_af", "__missing__")
        if af == "__missing__":
            raise NotAvailableHere(
                "gnomAD sites VCF is not present, so allele frequency is unknown. "
                "Absence of a frequency is NOT evidence of rarity: treating it as "
                "such would score every unannotated variant as if it were absent "
                "from gnomAD."
            )
        return scoring.frequency(af, source=self.version)


class SpliceAiAnnotator(Annotator):
    """SpliceAI and Pangolin delta scores. The highest-prior arm."""

    name = "spliceai"

    @property
    def version(self) -> str:
        return "spliceai/1.3.1+pangolin/1.0.2"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        sai, pan = record.get("spliceai_delta"), record.get("pangolin_delta")
        if sai is None and pan is None:
            raise NotAvailableHere(
                "Neither SpliceAI nor Pangolin has been run. This is the arm with "
                "the highest prior (RECON.md class 1); a ranked list produced "
                "without it is not an Arm B result."
            )
        return scoring.splice_effect(sai, pan, candidate.variant_class, self.version)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class PipelineResult:
    candidates: list[Candidate]
    n_records_seen: int
    n_records_in_panel: int
    annotators_run: list[str]
    annotators_unavailable: dict[str, str]
    contig_naming: str = "ensembl_nochr"
    n_records_unparseable: int = 0

    #: Evidence classes a complete Track 1 run would include, and the
    #: annotator that supplies each. Reported against what actually ran, so a
    #: result can never imply completeness it does not have.
    REQUIRED_FOR_COMPLETE = {
        "panel": "gene plausibility",
        "region": "splice distance and region class",
        "quality": "artefact indicators",
        "vep": "protein consequence (missense, nonsense, frameshift)",
        "gnomad": "population allele frequency",
        "spliceai": "splicing prediction (Arm B, highest prior)",
    }

    @property
    def missing_evidence_classes(self) -> dict[str, str]:
        """Evidence a complete run would have and this run does not.

        Distinct from ``annotators_unavailable``, which lists only annotators
        that were configured and then failed. An annotator absent from the
        configuration entirely produces no error at all, so reporting only
        failures would let a three-annotator run look complete.
        """
        ran = set(self.annotators_run)
        return {name: what for name, what in self.REQUIRED_FOR_COMPLETE.items()
                if name not in ran}

    @property
    def is_complete(self) -> bool:
        return not self.missing_evidence_classes

    def completeness_note(self) -> str:
        """One line for the top of any report built from this result."""
        if self.is_complete:
            return "Complete run: all evidence classes present."
        missing = self.missing_evidence_classes
        return ("INCOMPLETE RUN. Missing evidence: "
                + "; ".join(f"{k} ({v})" for k, v in missing.items())
                + ". Rankings from this run must not be presented as Track 1 findings.")

    def top(self, n: int = 20) -> list[Candidate]:
        return self.candidates[:n]

    def genes_ranked(self) -> list[str]:
        """Ranked gene symbols, deduplicated, keeping the best rank per gene."""
        seen: dict[str, None] = {}
        for c in self.candidates:
            if c.gene not in seen:
                seen[c.gene] = None
        return list(seen)

    def rank_of_gene(self, gene: str) -> int | None:
        genes = self.genes_ranked()
        return genes.index(gene) + 1 if gene in genes else None


class Track1Pipeline:
    """Assemble candidates from a VCF and rank them."""

    def __init__(
        self,
        gene_model: GeneModel,
        annotators: Iterable[Annotator] | None = None,
        panel_tsv: str | pathlib.Path = "config/gene_panels/mitotic_extended.tsv",
        restrict_to_panel: bool = True,
    ) -> None:
        self.gene_model = gene_model
        self.restrict_to_panel = restrict_to_panel
        self.panel_annotator = PanelAnnotator(panel_tsv)
        self.annotators: list[Annotator] = list(annotators) if annotators is not None else [
            self.panel_annotator,
            RegionAnnotator(gene_model),
            QualityAnnotator(),
        ]

    def _records(self, vcf: str | pathlib.Path, regions: list[str] | None) -> Iterator[dict]:
        if regions:
            for r in regions:
                yield from iter_vcf(vcf, region=r)
        else:
            yield from iter_vcf(vcf)

    @staticmethod
    def _detect_naming(vcf: str | pathlib.Path) -> str:
        """Read the contig convention from the VCF header rather than assuming it.

        The proband callset is no-chr; the GIAB background used for spike-in
        recall is chr-prefixed. Hardcoding either one makes GenomicPosition
        reject every record from the other, and the exception handler in run()
        turns that into a silent zero-candidate result. Detect, do not assume.
        """
        from mva.track1.spikein import detect_naming
        return detect_naming(vcf)

    @staticmethod
    def _regions_for_naming(regions: list[str] | None, naming: str) -> list[str] | None:
        """Region strings come from a no-chr BED; rewrite them if the VCF is
        chr-prefixed, otherwise every region query silently returns nothing."""
        if not regions:
            return regions
        if naming == "ucsc_chr":
            return [r if r.startswith("chr") else f"chr{r}" for r in regions]
        return [r[3:] if r.startswith("chr") else r for r in regions]

    def run(
        self,
        vcf: str | pathlib.Path,
        regions: list[str] | None = None,
        hpo: list[str] | None = None,
    ) -> PipelineResult:
        candidates: list[Candidate] = []
        unavailable: dict[str, str] = {}
        n_seen = n_panel = n_unparseable = 0

        naming = self._detect_naming(vcf)
        regions = self._regions_for_naming(regions, naming)

        for rec in self._records(vcf, regions):
            n_seen += 1
            # Multi-allelic ALT fields arrive comma separated; take each in turn.
            for alt in rec["alt"].split(","):
                if not alt or alt == ".":
                    continue
                try:
                    pos = GenomicPosition(
                        build="GRCh38", naming=naming,   # type: ignore[arg-type]
                        contig=rec["contig"], pos=rec["pos"],
                        ref=rec["ref"], alt=alt,
                    )
                except Exception:
                    # Symbolic alleles, decoy contigs and malformed records.
                    # Counted, not silently swallowed: if this equals the record
                    # count, something systematic is wrong.
                    n_unparseable += 1
                    continue

                gene, vclass, dist = self.gene_model.classify(pos)
                if gene is None:
                    continue
                if self.restrict_to_panel and gene not in self.panel_annotator.panel:
                    continue
                n_panel += 1

                gt = rec.get("gt", "")
                alleles = gt.replace("|", "/").split("/")
                zyg = (Zygosity.HOM_ALT if len(alleles) == 2 and alleles[0] == alleles[1]
                       and alleles[0] not in ("0", ".")
                       else Zygosity.HET if len(alleles) == 2 and alleles[0] != alleles[1]
                       else Zygosity.UNKNOWN)

                cand = Candidate(
                    gene=gene, position=pos, variant_class=vclass,
                    zygosity=zyg, arm="A_baseline",
                )
                rec_ext = {**rec, "splice_distance": dist}

                for ann in self.annotators:
                    try:
                        for ev in ann.annotate(cand, rec_ext):
                            cand.add(ev)
                    except NotAvailableHere as exc:
                        unavailable.setdefault(ann.name, str(exc))

                if cand.evidence:
                    candidates.append(cand)

        if n_seen and n_unparseable == n_seen:
            raise RuntimeError(
                f"every one of {n_seen} records failed to parse as a "
                f"{naming} GRCh38 position. This is a contig-convention or "
                f"build mismatch, not an empty result."
            )

        return PipelineResult(
            candidates=scoring.rank(candidates),
            n_records_seen=n_seen,
            n_records_in_panel=n_panel,
            annotators_run=[a.name for a in self.annotators if a.name not in unavailable],
            annotators_unavailable=unavailable,
            contig_naming=naming,
            n_records_unparseable=n_unparseable,
        )
