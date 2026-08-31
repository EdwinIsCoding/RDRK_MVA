"""Transparent additive scoring, per plan section 6.7.

Deliberately not a learned model. With n=1 there is nothing to train on, and a
black-box ranker over a single proband would be fitting noise while sounding
authoritative. The framework here is additive, ACMG/AMP-shaped, and every weight
is justified in this file where a reader can argue with it.

**Weights are frozen before the proband is scored.** They were fixed from the
reasoning below and may be tuned only against the positive-control benchmark in
``benchmarks/``, never against the proband's own ranked list. Changing a weight
after seeing the proband's results converts this from a prior into a
post-hoc rationalisation.

The calibration that shaped these weights
-----------------------------------------
``scripts/12_join_constraint.py`` measured gnomAD v4.1 constraint for the known
MVA genes:

    BUB1B   LOEUF 0.707  pLI 0.000     CEP57   LOEUF 0.740  pLI 0.000
    BUB1    LOEUF 0.727  pLI 0.000     BUB3    LOEUF 0.770  pLI 0.001
    TRIP13  LOEUF 0.592  pLI 0.192     CEP192  LOEUF 0.610  pLI 0.000
    SMC5    LOEUF 0.716  pLI 0.000     CENATAC LOEUF 1.227  pLI 0.000

**Not one is constrained.** Every pLI is approximately zero and every LOEUF sits
in the moderate-to-unconstrained band. This is exactly what autosomal recessive
biology predicts: LOEUF measures selection against *heterozygous* loss of
function, and heterozygous carriers of MVA alleles are healthy.

The consequence is concrete. Had constraint been given the weight it usually
carries in a dominant-disease pipeline, this scoring function would actively
deprioritise every known answer. Constraint therefore contributes at most a
small positive nudge and **never a penalty**. This is the single most important
design decision in this file and it came from measurement, not from assumption.

Constraint on the X chromosome is a different matter
----------------------------------------------------
The proband is male (chrX heterozygosity 0.062 against 0.620 autosomal), so his
X is hemizygous. A loss-of-function variant on the X is therefore fully exposed
to selection in males, exactly as a dominant autosomal variant would be. LOEUF
measures precisely that, so on the X it is informative rather than misleading.

Measured across the panel:

    X-linked panel genes (gnomAD v2.1.1)  median LOEUF 0.297,  13/20 with pLI > 0.9
    known autosomal MVA genes (v4.1)      median LOEUF 0.727,   0/9  with pLI > 0.9

Same metric, opposite meaning, and the difference is mechanistic rather than
incidental: autosomal recessive carriers are healthy so selection is blind to
their heterozygous LoF, while a hemizygous male has nowhere to hide.

So constraint carries a larger bonus for X-linked genes than autosomal ones. It
is still never a penalty, because an unconstrained X gene is not thereby
excluded.

Note the release difference. gnomAD v4.1's constraint file is autosomes only, so
X values come from v2.1.1, which is a different build and an exome-only callset.
Values from the two releases are not directly comparable and the panel records
which release each came from.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva.evidence import Candidate, Evidence, SourceType, VariantClass, Zygosity

#: Bumped whenever a weight changes. Recorded on every scored candidate so a
#: ranked list can always be traced to the weights that produced it.
WEIGHTS_VERSION = "1.0.0-frozen-2026-08-31"


@dataclass(frozen=True)
class Weights:
    """Every field carries its justification in the comment above it."""

    # --- gene plausibility -------------------------------------------------
    # A known MVA gene is the strongest prior available. Three established loci
    # plus five reported ones; a variant in one of these starts well ahead.
    gene_known_mva_tier1: float = 3.0
    gene_known_mva_tier2: float = 2.0
    # Core panel membership from converging independent sources (GO, Reactome,
    # STRING). Tier 2 requires three sources, tier 3 two, tier 4 a single
    # curated annotation to a mechanism-specific GO term.
    gene_panel_tier2: float = 1.5
    gene_panel_tier3: float = 1.0
    gene_panel_tier4: float = 0.5
    # Constraint: small positive nudge only, never a penalty. See the module
    # docstring. An unconstrained gene gets zero, not a deduction.
    gene_constrained_bonus: float = 0.25
    # Constraint means something different on the X for a male proband, so it
    # is weighted differently. See CONSTRAINT_ON_X below.
    gene_constrained_bonus_hemizygous: float = 1.25

    # --- variant effect ----------------------------------------------------
    # Predicted loss of function in a gene where LoF is the known mechanism.
    effect_lof_high_confidence: float = 3.0        # PVS1 territory
    effect_splice_canonical: float = 2.5
    # SpliceAI and Pangolin delta scores. The 0.5 threshold is SpliceAI's own
    # recommended high-precision cutoff; 0.8 is high confidence.
    effect_spliceai_ge_0_8: float = 2.5
    effect_spliceai_ge_0_5: float = 1.5
    effect_spliceai_ge_0_2: float = 0.5
    # Agreement between two independently trained splice models is worth more
    # than either alone, because their failure modes are not identical.
    effect_splice_tools_agree: float = 1.0
    effect_missense_alphamissense_pathogenic: float = 1.5
    effect_missense_alphamissense_ambiguous: float = 0.25

    # --- frequency ---------------------------------------------------------
    # Absent from gnomAD entirely. PM2 supporting.
    freq_absent_gnomad: float = 1.5
    freq_ultra_rare: float = 1.0                   # popmax AF < 1e-5
    freq_rare: float = 0.5                         # popmax AF < 1e-4
    # A recessive allele can be present in gnomAD at low frequency, and MVA
    # carriers are healthy, so this penalty is mild rather than disqualifying.
    freq_too_common_penalty: float = -4.0          # popmax AF > 1e-3

    # --- inheritance and phasing ------------------------------------------
    # Singleton. De novo and segregation contribute nothing here; the fields
    # exist so the framework is complete, and are unreachable on this dataset.
    inh_de_novo_confirmed: float = 0.0             # impossible: no parents
    inh_segregates: float = 0.0                    # impossible: no parents
    # What IS available in a singleton: HaplotypeCaller physical phasing.
    inh_phased_in_trans: float = 2.0               # PM3-like, read-backed
    inh_phased_in_cis_penalty: float = -3.0        # two hits on one haplotype
    inh_homozygous_in_ar_gene: float = 1.5

    # --- functional --------------------------------------------------------
    # No RNA-seq exists for this proband, so these are structurally
    # unreachable. Kept, weighted high, and reported as unreachable, so that
    # the report states what the strongest available evidence would have been.
    func_rnaseq_aberrant_splicing: float = 4.0     # unreachable: no RNA-seq
    func_rnaseq_expression_outlier: float = 2.5    # unreachable: no RNA-seq

    # --- penalties ---------------------------------------------------------
    # 3.5% of the callset carries MQ40. Those regions are a deliberate search
    # target, so the penalty is a caution, not an exclusion.
    pen_low_mapping_quality: float = -1.0
    pen_segdup_or_repeat: float = -0.75
    pen_strand_bias: float = -1.0
    pen_low_depth: float = -1.0                    # DP < 10 against a median of 42
    pen_failed_hard_filter: float = -0.5


WEIGHTS = Weights()


#: Weights that cannot fire on this dataset, with the reason. Reported
#: explicitly rather than left as silent zeros, because a reader deserves to
#: know which evidence types were structurally unavailable rather than simply
#: absent.
UNREACHABLE: dict[str, str] = {
    "inh_de_novo_confirmed": "singleton: no parental samples",
    "inh_segregates": "singleton: no parental samples",
    "func_rnaseq_aberrant_splicing": "no RNA-seq in the distributed dataset",
    "func_rnaseq_expression_outlier": "no RNA-seq in the distributed dataset",
}


def gene_plausibility(
    gene: str,
    panel: dict[str, dict[str, str]],
) -> list[Evidence]:
    """Evidence from what the gene is, before looking at the variant."""
    out: list[Evidence] = []
    row = panel.get(gene)
    if row is None:
        return out

    if row.get("known_mva_gene") == "yes":
        tier = row.get("panel_tier", "1")
        w = WEIGHTS.gene_known_mva_tier1 if tier == "1" else WEIGHTS.gene_known_mva_tier2
        out.append(Evidence(
            criterion="known_mva_gene",
            statement=f"{gene} is an established or reported MVA gene",
            source_type=SourceType.PANEL,
            source_id=f"config/gene_panels/mva_known.tsv#{gene}",
            weight=w,
        ))
    else:
        tier_weight = {
            "2": WEIGHTS.gene_panel_tier2,
            "3": WEIGHTS.gene_panel_tier3,
            "4": WEIGHTS.gene_panel_tier4,
        }.get(row.get("panel_tier", ""), 0.0)
        if tier_weight:
            out.append(Evidence(
                criterion="mitotic_panel_membership",
                statement=(f"{gene} is on the extended mitotic panel at tier "
                           f"{row['panel_tier']}, nominated by {row.get('sources', '')}"),
                source_type=SourceType.PANEL,
                source_id=f"config/gene_panels/mitotic_extended.tsv#{gene}",
                weight=tier_weight,
                detail={"n_sources": row.get("n_sources", ""),
                        "go_terms": row.get("specific_go_terms", "")[:200]},
            ))

    # Constraint: bonus only, never a penalty, and weighted by whether the
    # gene is hemizygous in this proband. See the module docstring.
    source = row.get("constraint_source", "")
    hemizygous = source.startswith("gnomAD v2")   # the X/Y fill
    loeuf_str = row.get("gnomad_v2_loeuf" if hemizygous else "gnomad_loeuf", "NA")
    pli_str = row.get("gnomad_v2_pli" if hemizygous else "gnomad_pli", "NA")

    try:
        loeuf = float(loeuf_str)
    except (TypeError, ValueError):
        loeuf = None

    if loeuf is not None and loeuf < 0.60:
        band = "highly constrained" if loeuf < 0.35 else "constrained"
        if hemizygous:
            statement = (f"{gene} is {band} (LOEUF {loeuf:.3f}, pLI {pli_str}). "
                         f"The proband is male, so this X-linked gene is hemizygous "
                         f"and loss of function is fully exposed to selection; "
                         f"constraint is informative here in a way it is not for "
                         f"an autosomal recessive gene.")
            weight = WEIGHTS.gene_constrained_bonus_hemizygous
            source_id = "gnomAD v2.1.1 constraint (X/Y)"
        else:
            statement = f"{gene} is {band} (LOEUF {loeuf:.3f})"
            weight = WEIGHTS.gene_constrained_bonus
            source_id = "gnomAD v4.1 constraint metrics"
        out.append(Evidence(
            criterion="gnomad_constraint",
            statement=statement,
            source_type=SourceType.GNOMAD, source_id=source_id,
            weight=weight,
            detail={"loeuf": loeuf_str, "pli": pli_str, "hemizygous": hemizygous},
        ))
    return out


def splice_effect(
    spliceai_delta: float | None,
    pangolin_delta: float | None,
    variant_class: VariantClass,
    tool_versions: str = "spliceai/1.3.1+pangolin/1.0.2",
) -> list[Evidence]:
    """Evidence from splicing prediction. Arm B, the highest-prior arm.

    Note what this function does NOT claim. Without RNA-seq these are
    predictions of an effect on splicing, not observations of one. The wording
    of every statement here keeps that distinction, because it is the thing most
    likely to be lost between analysis and write-up.
    """
    out: list[Evidence] = []
    best = max([d for d in (spliceai_delta, pangolin_delta) if d is not None], default=None)
    if best is None:
        return out

    if best >= 0.8:
        w, band = WEIGHTS.effect_spliceai_ge_0_8, "high confidence"
    elif best >= 0.5:
        w, band = WEIGHTS.effect_spliceai_ge_0_5, "recommended threshold"
    elif best >= 0.2:
        w, band = WEIGHTS.effect_spliceai_ge_0_2, "permissive threshold"
    else:
        return out

    out.append(Evidence(
        criterion="splice_prediction",
        statement=(f"predicted to alter splicing, max delta {best:.2f} ({band}). "
                   f"Prediction only: no RNA-seq exists for this proband, so no "
                   f"aberrant junction has been observed."),
        source_type=SourceType.TOOL,
        source_id=tool_versions,
        weight=w,
        acmg_code="PP3" if best >= 0.5 else None,
        detail={"spliceai_delta": spliceai_delta, "pangolin_delta": pangolin_delta,
                "variant_class": variant_class.value},
    ))

    if (spliceai_delta is not None and pangolin_delta is not None
            and spliceai_delta >= 0.5 and pangolin_delta >= 0.5):
        out.append(Evidence(
            criterion="splice_tools_agree",
            statement=("SpliceAI and Pangolin agree above threshold; the two models "
                       "are independently trained and do not share failure modes"),
            source_type=SourceType.TOOL,
            source_id=tool_versions,
            weight=WEIGHTS.effect_splice_tools_agree,
        ))
    return out


def frequency(popmax_af: float | None, source: str = "gnomAD v4.1 joint") -> list[Evidence]:
    """Rarity. MVA is recessive and carriers are healthy, so a candidate allele
    may legitimately appear in gnomAD at low frequency."""
    if popmax_af is None:
        return [Evidence(
            criterion="absent_from_gnomad",
            statement="not observed in gnomAD",
            source_type=SourceType.GNOMAD, source_id=source,
            weight=WEIGHTS.freq_absent_gnomad, acmg_code="PM2",
        )]
    if popmax_af > 1e-3:
        return [Evidence(
            criterion="too_common",
            statement=(f"gnomAD popmax AF {popmax_af:.2e} is too common for a "
                       f"fully penetrant recessive allele in a severe phenotype"),
            source_type=SourceType.GNOMAD, source_id=source,
            weight=WEIGHTS.freq_too_common_penalty, acmg_code="BS1",
        )]
    if popmax_af < 1e-5:
        w, label = WEIGHTS.freq_ultra_rare, "ultra-rare"
    elif popmax_af < 1e-4:
        w, label = WEIGHTS.freq_rare, "rare"
    else:
        return []
    return [Evidence(
        criterion="rare_in_gnomad",
        statement=f"{label} in gnomAD, popmax AF {popmax_af:.2e}",
        source_type=SourceType.GNOMAD, source_id=source,
        weight=w, acmg_code="PM2" if popmax_af < 1e-5 else None,
    )]


def phasing(
    phase: str | None,
    zygosity: Zygosity,
    pid: str | None = None,
) -> list[Evidence]:
    """Phasing evidence available to a singleton.

    ``phase`` is one of "trans", "cis" or None. In a singleton it can only come
    from HaplotypeCaller's PGT/PID physical phasing tags, which resolve variants
    within one assembly region, or from population-based phasing. Where neither
    resolves it, this returns nothing rather than assuming *trans*: assuming the
    convenient phase is how a compound heterozygote gets reported that is
    actually two variants on the same haplotype.
    """
    if zygosity is Zygosity.HOM_ALT:
        return [Evidence(
            criterion="homozygous",
            statement="homozygous alternate in a gene with a recessive mechanism",
            source_type=SourceType.CALLSET, source_id="proband VCF FORMAT/GT",
            weight=WEIGHTS.inh_homozygous_in_ar_gene,
        )]
    if phase == "trans":
        return [Evidence(
            criterion="phased_in_trans",
            statement=("read-backed phasing places this allele in trans with the "
                       "second candidate allele, consistent with a compound heterozygote"),
            source_type=SourceType.CALLSET,
            source_id=f"proband VCF FORMAT/PGT,PID{f' PID={pid}' if pid else ''}",
            weight=WEIGHTS.inh_phased_in_trans, acmg_code="PM3",
        )]
    if phase == "cis":
        return [Evidence(
            criterion="phased_in_cis",
            statement=("read-backed phasing places both candidate alleles on the same "
                       "haplotype, which is not consistent with a recessive mechanism"),
            source_type=SourceType.CALLSET,
            source_id=f"proband VCF FORMAT/PGT,PID{f' PID={pid}' if pid else ''}",
            weight=WEIGHTS.inh_phased_in_cis_penalty, acmg_code="BP2",
        )]
    return []


def quality_penalties(
    filter_field: str,
    depth: int | None,
    strand_bias_fs: float | None = None,
    in_segdup: bool = False,
) -> list[Evidence]:
    """Artefact indicators. Cautions rather than exclusions: the low
    mapping-quality regions these flag are exactly where a diagnostic pipeline
    under-calls, and therefore where this project is deliberately looking."""
    out: list[Evidence] = []
    if filter_field and filter_field not in ("PASS", "."):
        tags = filter_field.split(";")
        if "MQ40" in tags:
            out.append(Evidence(
                criterion="low_mapping_quality",
                statement=("RMS mapping quality below 40; the region is repetitive or "
                           "paralogous. Flagged rather than excluded, since these regions "
                           "are a deliberate search target"),
                source_type=SourceType.CALLSET, source_id="proband VCF FILTER=MQ40",
                weight=WEIGHTS.pen_low_mapping_quality,
            ))
        other = [t for t in tags if t != "MQ40"]
        if other:
            out.append(Evidence(
                criterion="failed_hard_filter",
                statement=f"failed GATK hard filter(s): {','.join(other)}",
                source_type=SourceType.CALLSET,
                source_id=f"proband VCF FILTER={filter_field}",
                weight=WEIGHTS.pen_failed_hard_filter,
            ))
    if depth is not None and depth < 10:
        out.append(Evidence(
            criterion="low_depth",
            statement=f"read depth {depth} against a callset median of 42",
            source_type=SourceType.CALLSET, source_id="proband VCF FORMAT/DP",
            weight=WEIGHTS.pen_low_depth,
        ))
    if strand_bias_fs is not None and strand_bias_fs > 60:
        out.append(Evidence(
            criterion="strand_bias",
            statement=f"Fisher strand bias {strand_bias_fs:.1f}",
            source_type=SourceType.CALLSET, source_id="proband VCF INFO/FS",
            weight=WEIGHTS.pen_strand_bias,
        ))
    if in_segdup:
        out.append(Evidence(
            criterion="segmental_duplication",
            statement="within a segmental duplication or low-complexity region",
            source_type=SourceType.TOOL, source_id="ucsc/segdups-GRCh38",
            weight=WEIGHTS.pen_segdup_or_repeat,
        ))
    return out


def rank(candidates: list[Candidate]) -> list[Candidate]:
    """Rank by score, then by evidence-source diversity as a tie-break.

    Diversity as the tie-break is deliberate: two candidates on the same score,
    one supported by four independent kinds of evidence and one by four
    correlated readings of the same tool, are not equally believable.
    """
    return sorted(
        candidates,
        key=lambda c: (c.score, len({e.source_type for e in c.evidence})),
        reverse=True,
    )
