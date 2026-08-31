"""Concrete annotators that run without a GPU.

These replace the raising stubs in ``mva.track1.pipeline`` for two of the three
evidence classes that were previously treated as GPU-blocked. Neither actually
needed a GPU; both needed a route around a large download, and both have one.

``GnomadFrequencyAnnotator``
    Reads a locally cached table sliced from gnomAD's public per-chromosome
    sites files by remote range request (``scripts/16_fetch_gnomad_panel.py``).
    No proband coordinate was transmitted to obtain it.

``ConsequenceAnnotator``
    Wraps ``bcftools csq``, which computes protein consequence from a GFF3 and a
    reference FASTA. That is roughly 3 GB of reference rather than VEP's 25 GB
    cache, and bcftools is already a dependency.

What ``bcftools csq`` does not give, and VEP does: AlphaMissense, CADD and REVEL
scores. Those are missense pathogenicity predictors and they remain unavailable.
So a missense variant here gets a *class* but not a *pathogenicity prediction*,
and the scoring function's missense term stays unreachable. That is stated
rather than papered over.
"""

from __future__ import annotations

import gzip
import pathlib
import subprocess

from mva.evidence import Candidate, Evidence, SourceType, VariantClass
from mva.track1 import scoring
from mva.track1.pipeline import Annotator, NotAvailableHere

GNOMAD_TABLE = "refs/gnomad_panel/panel_af.tsv.gz"


class GnomadFrequencyAnnotator(Annotator):
    """Population allele frequency and, more usefully here, homozygote count."""

    name = "gnomad"

    #: Half-width, in bases, of the window used to decide whether a position is
    #: covered by the lookup at all. See ``_is_covered``.
    COVERAGE_WINDOW = 200

    def __init__(self, table: str | pathlib.Path = GNOMAD_TABLE,
                 version: str = "gnomAD v4.1 exomes+genomes"):
        self.path = pathlib.Path(table)
        self._version = version
        self.af: dict[str, float] = {}
        self.nhomalt: dict[str, int] = {}
        #: Sorted positions per contig, for the coverage test.
        self._positions: dict[str, list[int]] = {}
        if self.path.exists():
            import collections
            pos: dict[str, list[int]] = collections.defaultdict(list)
            with gzip.open(self.path, "rt") as fh:
                next(fh, None)
                for line in fh:
                    key, af, nh = line.rstrip("\n").split("\t")
                    try:
                        self.af[key] = float(af)
                        self.nhomalt[key] = int(nh)
                    except ValueError:
                        continue
                    c, p, *_ = key.split(":")
                    pos[c].append(int(p))
            self._positions = {c: sorted(v) for c, v in pos.items()}

    def _is_covered(self, contig: str, pos: int) -> bool:
        """Is this position inside a region the lookup actually assays?

        This distinction is load-bearing and was missing at first. An
        exomes-only lookup contains no intronic variants, so an intronic
        position is not in the table because it was never assayed, not because
        the allele is absent from the population. Treating the two the same
        turned 1,007 ordinary intronic variants across the rhabdomyosarcoma
        genes into "absent from gnomAD", which is the single most promoting
        piece of evidence the scoring function has.

        A position counts as covered if the lookup holds any variant within
        ``COVERAGE_WINDOW`` bases. gnomAD variant density in assayed sequence is
        far higher than one per 200 bp, so this separates assayed from
        unassayed reliably without needing the coverage tracks.
        """
        ps = self._positions.get(contig)
        if not ps:
            return False
        import bisect
        i = bisect.bisect_left(ps, pos)
        for j in (i - 1, i):
            if 0 <= j < len(ps) and abs(ps[j] - pos) <= self.COVERAGE_WINDOW:
                return True
        return False

    @property
    def version(self) -> str:
        return f"{self._version} ({len(self.af):,} panel variants)"

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        if not self.af:
            raise NotAvailableHere(
                f"{self.path} is absent or empty. Run "
                "scripts/16_fetch_gnomad_panel.py. Absence of a frequency is NOT "
                "evidence of rarity: scoring every unannotated variant as absent "
                "from gnomAD would promote exactly the wrong candidates."
            )
        p = candidate.position
        if p is None:
            return []
        key = f"{p.contig}:{p.pos}:{p.ref}:{p.alt}"

        popmax = self.af.get(key)
        if popmax is None and not self._is_covered(p.contig, p.pos):
            # Not assayed, so nothing can be said about its frequency. Zero
            # weight, and the candidate carries the gap visibly rather than
            # being rewarded for it.
            return [Evidence(
                criterion="gnomad_not_assayed",
                statement=("no gnomAD data at this position in the lookup used. "
                           "This is a coverage gap, not evidence of rarity, and no "
                           "frequency conclusion can be drawn."),
                source_type=SourceType.GNOMAD, source_id=self._version,
                weight=0.0,
                detail={"lookup": self.path.name},
            )]
        out = scoring.frequency(popmax, source=self._version)

        # Homozygote count. For a severe recessive paediatric phenotype this is
        # the sharpest single filter available: a variant observed homozygous in
        # a population reference, in people who are by construction not severely
        # affected as children, cannot be causal in the homozygous state.
        nh = self.nhomalt.get(key)
        if nh is not None and nh > 0:
            if candidate.zygosity.value == "hom_alt":
                out.append(Evidence(
                    criterion="homozygotes_in_gnomad",
                    statement=(f"{nh} homozygous individuals for this allele in gnomAD. "
                               f"The proband is homozygous here, so this allele cannot "
                               f"be causal in the homozygous state"),
                    source_type=SourceType.GNOMAD, source_id=self._version,
                    weight=-5.0, acmg_code="BS2",
                    detail={"nhomalt": nh},
                ))
            elif nh >= 5:
                out.append(Evidence(
                    criterion="homozygotes_in_gnomad",
                    statement=(f"{nh} homozygous individuals for this allele in gnomAD, "
                               f"which argues against a fully penetrant severe recessive "
                               f"allele even in trans"),
                    source_type=SourceType.GNOMAD, source_id=self._version,
                    weight=-1.5, acmg_code="BS2",
                    detail={"nhomalt": nh},
                ))
        return out


class ConsequenceAnnotator(Annotator):
    """Protein consequence via ``bcftools csq``.

    Runs once over the whole VCF rather than per variant: csq is haplotype-aware
    and needs to see a transcript's variants together to get frameshift and
    stop-gain right. The result is cached in a dict keyed on position.
    """

    name = "vep"   # occupies the same slot as VepAnnotator in the completeness check

    #: bcftools csq consequence strings mapped onto our classes.
    CSQ_MAP = {
        "stop_gained": VariantClass.NONSENSE,
        "stop_lost": VariantClass.NONSENSE,
        "start_lost": VariantClass.NONSENSE,
        "frameshift": VariantClass.FRAMESHIFT,
        "inframe_deletion": VariantClass.FRAMESHIFT,
        "inframe_insertion": VariantClass.FRAMESHIFT,
        "inframe_altering": VariantClass.FRAMESHIFT,
        "missense": VariantClass.MISSENSE,
        "synonymous": VariantClass.SYNONYMOUS,
        "splice_acceptor": VariantClass.SPLICE_CANONICAL,
        "splice_donor": VariantClass.SPLICE_CANONICAL,
        "splice_region": VariantClass.SPLICE_REGION,
        "5_prime_utr": VariantClass.UTR_PROMOTER,
        "3_prime_utr": VariantClass.UTR_PROMOTER,
        "intron": VariantClass.DEEP_INTRONIC,
        "non_coding": VariantClass.UNCLASSIFIED,
        "coding_sequence": VariantClass.UNCLASSIFIED,
    }

    #: Consequences that are loss of function in a gene where LoF is the
    #: established mechanism. PVS1 territory, though PVS1 proper also requires
    #: checking last-exon and NMD-escape rules, which is noted where used.
    LOF = {VariantClass.NONSENSE, VariantClass.FRAMESHIFT, VariantClass.SPLICE_CANONICAL}

    def __init__(self, csq_table: dict[str, tuple[VariantClass, str]] | None = None,
                 bcftools_version: str = "bcftools/1.24-csq"):
        self.table = csq_table or {}
        self._version = bcftools_version

    @property
    def version(self) -> str:
        return f"{self._version} ({len(self.table):,} annotated positions)"

    @classmethod
    def from_vcf(
        cls,
        vcf: str | pathlib.Path,
        gff3: str | pathlib.Path,
        fasta: str | pathlib.Path,
        regions_file: str | pathlib.Path | None = None,
    ) -> ConsequenceAnnotator:
        """Run ``bcftools csq`` once and cache the result."""
        for p, what in ((gff3, "GFF3"), (fasta, "reference FASTA")):
            if not pathlib.Path(p).exists():
                raise NotAvailableHere(
                    f"{what} not found at {p}. Run 'make downloads'. Without it "
                    "protein consequence is unavailable and Arm A is incomplete."
                )
        cmd = ["bcftools", "csq", "-f", str(fasta), "-g", str(gff3),
               "--local-csq", "-Ou", str(vcf)]
        if regions_file:
            cmd[2:2] = ["-R", str(regions_file)]
        query = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n"]

        csq = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = subprocess.Popen(query + ["-"], stdin=csq.stdout,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if csq.stdout:
            csq.stdout.close()

        table: dict[str, tuple[VariantClass, str]] = {}
        assert out.stdout is not None
        for line in out.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[4] in (".", ""):
                continue
            # BCSQ is a comma-separated list of pipe-delimited annotations;
            # the consequence is the first field of each.
            best_class, best_str = VariantClass.UNCLASSIFIED, f[4]
            for ann in f[4].split(","):
                cons = ann.split("|")[0].lstrip("@").split("&")[0]
                vc = cls.CSQ_MAP.get(cons)
                if vc is None:
                    continue
                # Keep the most severe consequence across transcripts.
                if best_class is VariantClass.UNCLASSIFIED or vc in cls.LOF:
                    best_class = vc
                    if vc in cls.LOF:
                        break
            table[f"{f[0]}:{f[1]}:{f[2]}:{f[3]}"] = (best_class, best_str[:200])
        out.wait()
        csq.wait()
        if out.returncode != 0:
            err = (out.stderr.read() if out.stderr else "")
            raise RuntimeError(f"bcftools csq failed: {err[:400]}")
        return cls(table)

    def annotate(self, candidate: Candidate, record: dict) -> list[Evidence]:
        if not self.table:
            raise NotAvailableHere(
                "bcftools csq has not been run, so protein consequence is "
                "unavailable. Missense, nonsense and frameshift classes cannot "
                "be assigned and Arm A must not be reported as complete."
            )
        p = candidate.position
        if p is None:
            return []
        hit = self.table.get(f"{p.contig}:{p.pos}:{p.ref}:{p.alt}")
        if hit is None:
            return []
        vclass, raw = hit
        # csq is authoritative for coding consequence; the coordinate-derived
        # class from regions.py stands only where csq says nothing.
        if vclass is not VariantClass.UNCLASSIFIED:
            candidate.variant_class = vclass

        out: list[Evidence] = []
        if vclass in self.LOF:
            out.append(Evidence(
                criterion="predicted_loss_of_function",
                statement=(f"{vclass.value} in {candidate.gene}. Loss of function is the "
                           f"established mechanism for the known MVA genes. PVS1 is cited "
                           f"as supporting rather than very strong: the last-exon and "
                           f"NMD-escape checks PVS1 proper requires have not been made."),
                source_type=SourceType.TOOL, source_id=self._version,
                weight=scoring.WEIGHTS.effect_lof_high_confidence,
                acmg_code="PVS1_supporting",
                detail={"bcsq": raw},
            ))
        elif vclass is VariantClass.MISSENSE:
            # Deliberately no weight. bcftools csq gives the class but no
            # pathogenicity prediction, and AlphaMissense, CADD and REVEL all
            # need VEP plugins. Scoring a missense variant on class alone would
            # treat every missense as equally suspicious.
            out.append(Evidence(
                criterion="missense_unscored",
                statement=("missense, but no pathogenicity predictor is available on "
                           "this host. AlphaMissense, CADD and REVEL need VEP plugins. "
                           "This variant is classified, not scored."),
                source_type=SourceType.TOOL, source_id=self._version,
                weight=0.0, detail={"bcsq": raw},
            ))
        return out
