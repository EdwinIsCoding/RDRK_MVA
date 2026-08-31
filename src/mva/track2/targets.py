"""Directional target nomination, per plan section 7.2.

The anti-goal this module exists to avoid
-----------------------------------------
Plan section 0.4 forbids using unsigned network proximity to decide whether to
inhibit or activate a target. PrimeKG and Hetionet edges are largely undirected
and unsigned: proximity tells you *which* protein sits near the disrupted
biology, never what to do to it. A repurposing candidate proposed without a
direction is not a hypothesis, it is a gene name.

So direction here comes from signed causal edges (OmniPath, which aggregates
SIGNOR, Reactome and others), and a target with no signed edge is reported as
**undirected and not nominable**, rather than assigned a plausible-sounding
direction.

The ambiguous case is handled explicitly
----------------------------------------
Real signed databases contain contradictions. BUB1B to CDC20 is annotated both
stimulatory and inhibitory across sources, because different papers measured
different things. OmniPath exposes ``consensus_stimulation`` and
``consensus_inhibition`` for exactly this. Where the consensus itself is split
or absent, this module returns ``AMBIGUOUS`` and the target does not get a
direction. Picking the more convenient sign would be the section 0.4 error with
extra steps.
"""

from __future__ import annotations

import dataclasses
import enum
import json
import pathlib
import time
import urllib.parse
import urllib.request

from mva.evidence import Evidence, SourceType

OMNIPATH = "https://omnipathdb.org/interactions"


class Direction(str, enum.Enum):
    INHIBIT = "inhibit"
    ACTIVATE = "activate"
    AMBIGUOUS = "ambiguous"        # signed edges exist but contradict
    UNSIGNED = "unsigned"          # edge exists, no sign recorded
    NO_EDGE = "no_edge"

    @property
    def is_actionable(self) -> bool:
        """Only a resolved sign supports a therapeutic proposal."""
        return self in (Direction.INHIBIT, Direction.ACTIVATE)


@dataclasses.dataclass(frozen=True)
class SignedEdge:
    source: str
    target: str
    is_stimulation: bool
    is_inhibition: bool
    consensus_stimulation: bool
    consensus_inhibition: bool
    sources: tuple[str, ...]
    pmids: tuple[str, ...]
    curation_effort: int

    @property
    def sign(self) -> Direction:
        """Resolve the edge's sign, preferring the consensus fields.

        The raw ``is_*`` flags are the union across contributing databases and
        are routinely both true; the consensus fields are OmniPath's own
        resolution. Where consensus is split or silent, the edge is ambiguous.
        """
        if self.consensus_stimulation and not self.consensus_inhibition:
            return Direction.ACTIVATE
        if self.consensus_inhibition and not self.consensus_stimulation:
            return Direction.INHIBIT
        if self.consensus_stimulation and self.consensus_inhibition:
            return Direction.AMBIGUOUS
        # No consensus recorded. Fall back to the raw flags only when they
        # themselves are unambiguous.
        if self.is_stimulation and not self.is_inhibition:
            return Direction.ACTIVATE
        if self.is_inhibition and not self.is_stimulation:
            return Direction.INHIBIT
        if self.is_stimulation and self.is_inhibition:
            return Direction.AMBIGUOUS
        return Direction.UNSIGNED


@dataclasses.dataclass(frozen=True)
class TargetNomination:
    gene: str
    direction: Direction
    rationale: str
    edges: tuple[SignedEdge, ...]
    evidence: tuple[Evidence, ...]

    @property
    def is_nominable(self) -> bool:
        """A target may be nominated only with a resolved direction.

        Plan section 7.2: "For each target state inhibit or activate, citing the
        signed edge. Unsigned KG proximity does not qualify as a direction."
        """
        return self.direction.is_actionable

    @property
    def blocking_reason(self) -> str | None:
        if self.is_nominable:
            return None
        return {
            Direction.AMBIGUOUS: (
                "signed edges exist but contradict each other; the literature "
                "does not agree on whether this interaction is stimulatory or "
                "inhibitory, so no therapeutic direction can be stated"),
            Direction.UNSIGNED: (
                "the interaction is recorded but unsigned; proximity alone "
                "cannot say whether to inhibit or activate"),
            Direction.NO_EDGE: (
                "no causal edge to the seed genes in any signed resource"),
        }[self.direction]


def fetch_signed_edges(
    genes: list[str],
    cache_dir: str | pathlib.Path = "results/track2/omnipath",
    timeout: int = 90,
) -> dict[str, list[SignedEdge]]:
    """Fetch signed causal interactions for each gene, with an on-disk cache.

    Cached so that a rerun does not depend on the API being up, which matters
    for the reproducibility artefact: a judge running ``make reproduce`` months
    later should not silently get a different answer because OmniPath changed.
    """
    cache = pathlib.Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out: dict[str, list[SignedEdge]] = {}

    for gene in genes:
        cached = cache / f"{gene}.tsv"
        if cached.exists():
            text = cached.read_text()
        else:
            q = urllib.parse.urlencode({
                "genesymbols": "yes",
                "partners": gene,
                "fields": "sources,references,curation_effort",
                "format": "tsv",
            })
            req = urllib.request.Request(f"{OMNIPATH}?{q}",
                                         headers={"User-Agent": "mva-hackathon-2026"})
            try:
                with urllib.request.urlopen(req, timeout=timeout) as fh:
                    text = fh.read().decode()
            except Exception:
                out[gene] = []
                continue
            cached.write_text(text)
            time.sleep(0.4)

        lines = text.strip().splitlines()
        if len(lines) < 2:
            out[gene] = []
            continue
        head = lines[0].split("\t")
        idx = {c: i for i, c in enumerate(head)}
        edges = []
        for line in lines[1:]:
            f = line.split("\t")
            if len(f) < len(head):
                continue

            def flag(col: str) -> bool:
                return f[idx[col]].strip().lower() in ("true", "1")

            refs = f[idx["references"]] if "references" in idx else ""
            pmids = tuple(sorted({p.split(":")[-1] for p in refs.split(";")
                                  if p and p.split(":")[-1].isdigit()}))
            try:
                effort = int(f[idx["curation_effort"]]) if "curation_effort" in idx else 0
            except ValueError:
                effort = 0
            edges.append(SignedEdge(
                source=f[idx["source_genesymbol"]],
                target=f[idx["target_genesymbol"]],
                is_stimulation=flag("is_stimulation"),
                is_inhibition=flag("is_inhibition"),
                consensus_stimulation=flag("consensus_stimulation"),
                consensus_inhibition=flag("consensus_inhibition"),
                sources=tuple(f[idx["sources"]].split(";")) if "sources" in idx else (),
                pmids=pmids,
                curation_effort=effort,
            ))
        out[gene] = edges
    return out


def nominate(
    gene: str,
    edges: list[SignedEdge],
    seed_genes: set[str],
    desired_effect_on_seed: Direction = Direction.ACTIVATE,
) -> TargetNomination:
    """Decide whether and how to modulate ``gene``, given signed edges.

    ``desired_effect_on_seed`` is what we want to happen to the disrupted
    pathway. For a loss-of-function recessive mechanism the goal is usually to
    restore or compensate for the seed gene's activity, so the default is
    ACTIVATE.

    The inference is deliberately shallow and stated: if ``gene`` inhibits a
    seed gene and we want the seed gene's activity up, then inhibiting ``gene``
    is the proposal. Anything deeper than one hop is not defensible from edge
    signs alone, and this function does not attempt it.
    """
    relevant = [e for e in edges
                if (e.target in seed_genes and e.source == gene)]
    if not relevant:
        return TargetNomination(gene, Direction.NO_EDGE,
                                "no signed causal edge from this gene to a seed gene",
                                (), ())

    signs = {e.sign for e in relevant}
    actionable = {s for s in signs if s.is_actionable}

    if not actionable:
        d = Direction.AMBIGUOUS if Direction.AMBIGUOUS in signs else Direction.UNSIGNED
        return TargetNomination(gene, d, "sign could not be resolved", tuple(relevant), ())
    if len(actionable) > 1:
        return TargetNomination(
            gene, Direction.AMBIGUOUS,
            "this gene both activates and inhibits seed genes across different "
            "edges; a single therapeutic direction cannot be stated",
            tuple(relevant), ())

    edge_sign = actionable.pop()
    # One hop of sign algebra, and no more.
    if desired_effect_on_seed is Direction.ACTIVATE:
        action = Direction.INHIBIT if edge_sign is Direction.INHIBIT else Direction.ACTIVATE
    else:
        action = Direction.ACTIVATE if edge_sign is Direction.INHIBIT else Direction.INHIBIT

    best = max(relevant, key=lambda e: e.curation_effort)
    evidence = []
    for pmid in best.pmids[:5]:
        evidence.append(Evidence(
            criterion="signed_causal_edge",
            statement=(f"{best.source} {edge_sign.value}s {best.target}; to "
                       f"{desired_effect_on_seed.value} {best.target}, "
                       f"{action.value} {best.source}"),
            source_type=SourceType.PUBMED,
            source_id=pmid,
            weight=0.0,   # Track 2 does not use the Track 1 additive score
            detail={"sources": ";".join(best.sources[:6]),
                    "curation_effort": best.curation_effort},
        ))

    return TargetNomination(
        gene=gene, direction=action,
        rationale=(f"{gene} {edge_sign.value}s {best.target}. To "
                   f"{desired_effect_on_seed.value} {best.target}, {action.value} {gene}. "
                   f"Curated by {len(best.sources)} sources."),
        edges=tuple(relevant), evidence=tuple(evidence),
    )


def summarise(nominations: list[TargetNomination]) -> dict:
    """Counts by direction, for the report. Non-nominable targets are counted
    and named, not silently dropped: how many candidate targets had to be
    rejected for want of a sign is itself a finding about the method."""
    by = {d: [n.gene for n in nominations if n.direction is d] for d in Direction}
    return {
        "nominable": sorted(by[Direction.INHIBIT] + by[Direction.ACTIVATE]),
        "inhibit": sorted(by[Direction.INHIBIT]),
        "activate": sorted(by[Direction.ACTIVATE]),
        "rejected_ambiguous": sorted(by[Direction.AMBIGUOUS]),
        "rejected_unsigned": sorted(by[Direction.UNSIGNED]),
        "rejected_no_edge": sorted(by[Direction.NO_EDGE]),
        "n_total": len(nominations),
        "n_nominable": len(by[Direction.INHIBIT]) + len(by[Direction.ACTIVATE]),
    }
