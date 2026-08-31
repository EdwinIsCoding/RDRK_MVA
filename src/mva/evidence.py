"""Evidence and candidate schemas.

Two rules from CLAUDE.md are enforced here in types rather than left to
discipline, because both are easy to violate silently and expensive to detect
later:

1. **No bare float scores.** Every contribution to a candidate's score arrives
   attached to an ``Evidence`` carrying a resolvable source identifier. A number
   whose provenance nobody can reconstruct is not evidence, and a ranked list
   built from such numbers cannot be audited by a judge or acted on by a
   clinician.

2. **No bare genomic positions.** ``GenomicPosition`` requires a build. The
   proband's callset is GRCh38 with Ensembl-style contig naming; most annotation
   resources ship GRCh37 or GRCh38 with a ``chr`` prefix. A position passed
   around as a bare integer is how a project silently mixes coordinate systems
   and reports a variant in the wrong gene.
"""

from __future__ import annotations

import enum
import re
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

Build = Literal["GRCh37", "GRCh38", "T2T-CHM13v2.0"]
ContigNaming = Literal["ensembl_nochr", "ucsc_chr"]


class SourceType(str, enum.Enum):
    """Where a piece of evidence came from. Each member has a resolution rule in
    ``SOURCE_PATTERNS`` so that an identifier can be machine-checked rather than
    trusted."""

    CLINVAR = "clinvar"
    PUBMED = "pubmed"
    ENSEMBL = "ensembl"
    GNOMAD = "gnomad"
    HGNC = "hgnc"
    GO = "go"
    REACTOME = "reactome"
    STRING = "string"
    OMIM = "omim"
    MONDO = "mondo"
    HPO = "hpo"
    CHEMBL = "chembl"
    DRUGBANK = "drugbank"
    CLINICALTRIALS = "clinicaltrials"
    TOOL = "tool"          # a computation: tool name plus version plus params
    PANEL = "panel"        # membership of a locally built gene panel
    CALLSET = "callset"    # a property of the proband VCF itself


#: Patterns every source identifier must match. Deliberately strict: a
#: malformed accession that looks plausible is worse than a rejected one,
#: because a reader cannot tell it is wrong without going to look it up.
SOURCE_PATTERNS: dict[SourceType, re.Pattern[str]] = {
    SourceType.CLINVAR: re.compile(r"^VCV\d{9}(\.\d+)?$"),
    SourceType.PUBMED: re.compile(r"^\d{7,8}$"),
    SourceType.ENSEMBL: re.compile(r"^ENS[GTP]\d{11}(\.\d+)?$"),
    SourceType.HGNC: re.compile(r"^HGNC:\d+$"),
    SourceType.GO: re.compile(r"^GO:\d{7}$"),
    SourceType.REACTOME: re.compile(r"^R-HSA-\d+$"),
    SourceType.OMIM: re.compile(r"^\d{6}$"),
    SourceType.MONDO: re.compile(r"^MONDO:\d{7}$"),
    SourceType.HPO: re.compile(r"^HP:\d{7}$"),
    SourceType.CHEMBL: re.compile(r"^CHEMBL\d+$"),
    SourceType.DRUGBANK: re.compile(r"^DB\d{5}$"),
    SourceType.CLINICALTRIALS: re.compile(r"^NCT\d{8}$"),
    # Free-form sources still have to say what produced them and at what
    # version, so that a result can be regenerated.
    SourceType.TOOL: re.compile(r"^[\w.+-]+/[\w.+-]+"),          # name/version[/params]
    SourceType.GNOMAD: re.compile(r"^gnomAD (v[\d.]+)"),
    SourceType.STRING: re.compile(r"^9606\.ENSP\d{11}$|^STRING:"),
    # Free-form descriptors of a local computation. Spaces are allowed
    # ("proband VCF FORMAT/PGT,PID"), leading or trailing whitespace is not,
    # and the identifier must say enough to locate what produced it.
    SourceType.PANEL: re.compile(r"^\S(.*\S)?$"),
    SourceType.CALLSET: re.compile(r"^\S(.*\S)?$"),
}


class GenomicPosition(BaseModel):
    """A build-tagged position. There is deliberately no constructor that
    accepts a bare integer position."""

    model_config = ConfigDict(frozen=True)

    build: Build
    naming: ContigNaming
    contig: str
    pos: int = Field(ge=1, description="1-based, VCF convention")
    ref: str = Field(min_length=1)
    alt: str = Field(min_length=1)

    @field_validator("ref", "alt")
    @classmethod
    def _acgt(cls, v: str) -> str:
        if not re.fullmatch(r"[ACGTN]+|\*|<[^>]+>", v.upper()):
            raise ValueError(f"not a valid allele: {v!r}")
        return v.upper()

    @model_validator(mode="after")
    def _naming_matches_contig(self) -> GenomicPosition:
        has_chr = self.contig.startswith("chr")
        if self.naming == "ucsc_chr" and not has_chr:
            raise ValueError(f"naming=ucsc_chr but contig {self.contig!r} has no chr prefix")
        if self.naming == "ensembl_nochr" and has_chr:
            raise ValueError(f"naming=ensembl_nochr but contig {self.contig!r} has a chr prefix")
        return self

    def to_naming(self, naming: ContigNaming) -> GenomicPosition:
        """Rename the contig. Coordinates are identical between the two
        conventions for the primary assembly, so only the label changes."""
        if naming == self.naming:
            return self
        contig = self.contig[3:] if naming == "ensembl_nochr" else f"chr{self.contig}"
        return self.model_copy(update={"contig": contig, "naming": naming})

    @property
    def spdi_like(self) -> str:
        return f"{self.contig}:{self.pos}:{self.ref}:{self.alt}"

    def __str__(self) -> str:
        return f"{self.build}:{self.contig}:{self.pos}:{self.ref}>{self.alt}"


class Evidence(BaseModel):
    """One traceable reason to move a candidate's score.

    ``weight`` is the contribution to the additive score of plan section 6.7.
    It may be negative: artefact indicators and low mappability subtract.
    """

    model_config = ConfigDict(frozen=True)

    criterion: str = Field(min_length=1, description="short slug, e.g. 'panel_membership'")
    statement: str = Field(min_length=1, description="what this evidence asserts, in words")
    source_type: SourceType
    source_id: str = Field(min_length=1, description="resolvable identifier")
    weight: float
    #: ACMG/AMP criterion code where one genuinely applies (PVS1, PM2, PP3...).
    #: Left None rather than stretched: inventing a code that does not fit is
    #: worse than reporting the evidence without one.
    acmg_code: str | None = None
    detail: dict[str, str | float | int | bool | None] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _source_id_resolves(self) -> Evidence:
        pattern = SOURCE_PATTERNS.get(self.source_type)
        if pattern and not pattern.match(self.source_id):
            raise ValueError(
                f"source_id {self.source_id!r} does not look like a {self.source_type.value} "
                f"identifier (expected to match {pattern.pattern}). "
                "CLAUDE.md rule 2: write TODO(source) rather than an invented identifier."
            )
        return self

    @field_validator("acmg_code")
    @classmethod
    def _known_acmg(cls, v: str | None) -> str | None:
        if v is None:
            return v
        if not re.fullmatch(r"(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(_\w+)?", v):
            raise ValueError(f"{v!r} is not an ACMG/AMP criterion code")
        return v

    @property
    def url(self) -> str | None:
        """A link a reader can follow. None where the source is a local
        computation with nothing to point at."""
        sid = self.source_id
        # VCV accessions are zero-padded to nine digits; the web path is not.
        clinvar_num = sid.split(".")[0][3:].lstrip("0") or "0"
        base = {
            SourceType.CLINVAR: f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_num}/",
            SourceType.PUBMED: f"https://pubmed.ncbi.nlm.nih.gov/{sid}/",
            SourceType.ENSEMBL: f"https://ensembl.org/Homo_sapiens/Gene/Summary?g={sid}",
            SourceType.HGNC: f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{sid}",
            SourceType.GO: f"https://amigo.geneontology.org/amigo/term/{sid}",
            SourceType.REACTOME: f"https://reactome.org/content/detail/{sid}",
            SourceType.OMIM: f"https://omim.org/entry/{sid}",
            SourceType.MONDO: f"https://monarchinitiative.org/{sid}",
            SourceType.HPO: f"https://hpo.jax.org/browse/term/{sid}",
            SourceType.CHEMBL: f"https://www.ebi.ac.uk/chembl/explore/compound/{sid}",
            SourceType.DRUGBANK: f"https://go.drugbank.com/drugs/{sid}",
            SourceType.CLINICALTRIALS: f"https://clinicaltrials.gov/study/{sid}",
        }
        return base.get(self.source_type)


class VariantClass(str, enum.Enum):
    """Recall is reported broken down by this, per plan section 5.2."""

    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift_or_indel"
    MISSENSE = "missense"
    SYNONYMOUS = "synonymous"
    SPLICE_CANONICAL = "splice_site_canonical"
    SPLICE_REGION = "splice_region_near"
    DEEP_INTRONIC = "deep_intronic"
    UTR_PROMOTER = "utr_or_promoter"
    STRUCTURAL = "structural"
    REPEAT = "repeat_expansion"
    UNCLASSIFIED = "unclassified"


class Zygosity(str, enum.Enum):
    HET = "het"
    HOM_ALT = "hom_alt"
    HEMI = "hemizygous"
    MOSAIC = "mosaic"
    UNKNOWN = "unknown"


class Candidate(BaseModel):
    """One ranked hypothesis. Emitted by every arm into a common table."""

    model_config = ConfigDict(frozen=False)

    gene: str
    position: GenomicPosition | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    variant_class: VariantClass = VariantClass.UNCLASSIFIED
    zygosity: Zygosity = Zygosity.UNKNOWN
    arm: str = Field(description="which analysis arm produced this, e.g. 'B_splicing'")
    evidence: list[Evidence] = Field(default_factory=list)

    #: The experiment that would falsify this candidate. Required before a
    #: candidate may be reported: a hypothesis with no falsification route is
    #: not a scientific claim, and plan section 9 requires one per candidate.
    falsifying_experiment: str | None = None

    @property
    def score(self) -> float:
        """Sum of evidence weights. There is no separate score field, so a score
        can never drift out of step with the evidence that justifies it."""
        return sum(e.weight for e in self.evidence)

    @property
    def acmg_codes(self) -> list[str]:
        return sorted({e.acmg_code for e in self.evidence if e.acmg_code})

    def add(self, evidence: Evidence) -> Candidate:
        self.evidence.append(evidence)
        return self

    def confidence(self) -> Annotated[str, "qualitative, deliberately coarse"]:
        """Coarse confidence band. Deliberately not a probability: with n=1 and
        no functional validation, a calibrated probability would be a fiction,
        and a fictional probability is more misleading than a word."""
        s = self.score
        n_independent = len({e.source_type for e in self.evidence})
        if s >= 6.0 and n_independent >= 4:
            return "high"
        if s >= 3.5 and n_independent >= 3:
            return "moderate"
        if s >= 2.0:
            return "low"
        return "very_low"

    def is_reportable(self) -> tuple[bool, str]:
        """Gate applied before a candidate reaches the report."""
        if not self.evidence:
            return False, "no evidence attached"
        if self.falsifying_experiment is None:
            return False, "no falsifying experiment specified"
        unresolvable = [e.source_id for e in self.evidence if e.source_id.startswith("TODO")]
        if unresolvable:
            return False, f"unresolved sources: {unresolvable}"
        return True, "ok"
