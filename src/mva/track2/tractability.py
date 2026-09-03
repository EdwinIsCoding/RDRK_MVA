"""Directional tractability: is the required direction pharmacologically available?

The check this module exists to force
-------------------------------------
Plan section 7.2 requires every nominated target to state inhibit or activate.
It is then possible, and common, to nominate a target whose required direction
no drug achieves. A candidate list of targets that can only be modulated the
wrong way is worse than no list, because it looks actionable.

So: for each target, ask ChEMBL what drugs actually exist and what they *do*,
then check whether any of them acts in the direction the biology requires.

What this found for the direct SAC-restoration axis
---------------------------------------------------
Running ``scripts/14_track2_direction_audit.py`` over the targets nominated by
one-hop signed-edge algebra from the MVA seed genes (PLK1, AURKA, AURKB, CDK1,
TTK, ATM, EGFR, CENPE):

    118 drug-mechanism records in ChEMBL
      0 activating (agonist, activator, positive allosteric modulator)
    118 inhibitors, antagonists or blockers

The mechanism requires **activation** of these kinases, to compensate for a
hypomorphic spindle assembly checkpoint. Every drug that exists does the
opposite. The direct axis is not merely unproven, it is pharmacologically
unavailable in the required direction.

That is an empirical confirmation of plan section 7.1: you cannot fix
constitutional aneuploidy with a small molecule. It is also a safety finding in
its own right, since activating mitotic kinases in a child with a cancer
predisposition syndrome would be contraindicated even if an activator existed.

The consequence for Track 2 is that effort belongs on the axes plan section 7.1
already prefers: proteotoxic stress mitigation, mitochondrial support, and above
all cancer chemoprevention and surveillance.
"""

from __future__ import annotations

import dataclasses
import enum
import json
import pathlib
import time
import urllib.parse
import urllib.request

from mva.track2.targets import Direction

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

#: ChEMBL ``action_type`` values that increase target activity.
ACTIVATING_ACTIONS = frozenset({
    "AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "OPENER",
    "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR",
    "STABILISER", "STABILIZER",
})

#: Values that decrease it.
INHIBITING_ACTIONS = frozenset({
    "INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR",
    "NEGATIVE MODULATOR", "DISRUPTING AGENT", "INVERSE AGONIST",
    "DOWNREGULATOR", "SEQUESTERING AGENT",
})


class LookupFailed(Exception):
    """A ChEMBL lookup did not complete.

    Distinct from a lookup that completed and found nothing. Returning None for
    both made a target whose request failed indistinguishable from a target with
    no ChEMBL entry, and both surfaced as the verdict ``unknown``. During a
    period when ChEMBL was serving intermittent HTTP 500s that turned transport
    failures into apparent biology.
    """


class Tractability(str, enum.Enum):
    AVAILABLE = "available"                  # a drug acts in the required direction
    WRONG_DIRECTION_ONLY = "wrong_direction" # drugs exist, all act the other way
    NO_DRUGS = "no_drugs"                    # no mechanism records at all
    UNKNOWN = "unknown"                      # lookup failed


@dataclasses.dataclass(frozen=True)
class TractabilityResult:
    gene: str
    chembl_target_id: str | None
    required_direction: Direction
    status: Tractability
    n_mechanisms: int
    n_activating: int
    n_inhibiting: int
    max_phase: float
    action_types: dict[str, int]
    example_drugs: tuple[str, ...]

    @property
    def is_actionable(self) -> bool:
        return self.status is Tractability.AVAILABLE

    @property
    def explanation(self) -> str:
        if self.status is Tractability.AVAILABLE:
            return (f"{self.n_activating if self.required_direction is Direction.ACTIVATE else self.n_inhibiting} "
                    f"drug mechanisms act in the required direction "
                    f"({self.required_direction.value}); max clinical phase {self.max_phase:g}")
        if self.status is Tractability.WRONG_DIRECTION_ONLY:
            return (f"{self.n_mechanisms} drug mechanisms exist for {self.gene}, none in the "
                    f"required direction ({self.required_direction.value}). "
                    f"Available actions: {', '.join(sorted(self.action_types))}. "
                    f"Proposing this target would propose modulating it the wrong way.")
        if self.status is Tractability.NO_DRUGS:
            return f"no drug mechanism records in ChEMBL for {self.gene}"
        return f"tractability lookup failed for {self.gene}"


def _get(url: str, timeout: int = 45, attempts: int = 4) -> dict:
    """Fetch with backoff. A single attempt turns a transient 500 into a
    permanent-looking absence of evidence."""
    last: Exception | None = None
    for i in range(attempts):
        req = urllib.request.Request(url,
                                     headers={"User-Agent": "mva-hackathon-2026"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.load(fh)
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (i + 1))
    raise LookupFailed(f"{url}: {last}")


def resolve_target(gene: str, cache_dir: pathlib.Path) -> str | None:
    """Find the human single-protein ChEMBL target for a gene symbol.

    Searching without filtering returns other organisms first for several of
    these genes (a PLK1 search returns Xenopus and rat before human), so the
    organism and target-type filters are load-bearing rather than tidiness.
    """
    cached = cache_dir / f"target_{gene}.json"
    if cached.exists():
        d = json.loads(cached.read_text())
    else:
        q = urllib.parse.urlencode({"q": gene, "format": "json", "limit": 25})
        # Deliberately no try/except: a failed lookup must reach the caller as
        # LookupFailed, not become an apparent absence of any ChEMBL target.
        d = _get(f"{CHEMBL}/target/search?{q}")
        cached.write_text(json.dumps(d))
        time.sleep(0.3)

    human = [t for t in d.get("targets", [])
             if t.get("organism") == "Homo sapiens"
             and t.get("target_type") == "SINGLE PROTEIN"]
    if not human:
        return None
    # Prefer a target whose preferred name mentions the symbol, so a search that
    # matches a complex or a relative does not silently substitute it.
    named = [t for t in human
             if gene.lower() in (t.get("pref_name") or "").lower().replace("-", "")]
    return (named or human)[0]["target_chembl_id"]


def assess(
    gene: str,
    required_direction: Direction,
    cache_dir: str | pathlib.Path = "results/track2/chembl",
) -> TractabilityResult:
    cache = pathlib.Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    try:
        tid = resolve_target(gene, cache)
    except LookupFailed:
        # Not "no target": "we could not find out". The caller must be able to
        # tell these apart, so it is re-raised rather than folded into UNKNOWN.
        raise
    if tid is None:
        return TractabilityResult(gene, None, required_direction, Tractability.UNKNOWN,
                                  0, 0, 0, 0.0, {}, ())

    cached = cache / f"mech_{tid}.json"
    if cached.exists():
        m = json.loads(cached.read_text())
    else:
        q = urllib.parse.urlencode({"target_chembl_id": tid, "format": "json", "limit": 500})
        m = _get(f"{CHEMBL}/mechanism?{q}")
        cached.write_text(json.dumps(m))
        time.sleep(0.3)

    mechs = m.get("mechanisms", [])
    actions: dict[str, int] = {}
    for x in mechs:
        a = (x.get("action_type") or "UNKNOWN").upper()
        actions[a] = actions.get(a, 0) + 1

    n_act = sum(n for a, n in actions.items() if a in ACTIVATING_ACTIONS)
    n_inh = sum(n for a, n in actions.items() if a in INHIBITING_ACTIONS)
    phases = [float(x["max_phase"]) for x in mechs if x.get("max_phase") is not None]
    max_phase = max(phases, default=0.0)
    drugs = tuple(sorted({x.get("molecule_chembl_id", "") for x in mechs if x.get("molecule_chembl_id")})[:5])

    if not mechs:
        status = Tractability.NO_DRUGS
    else:
        wanted = n_act if required_direction is Direction.ACTIVATE else n_inh
        status = Tractability.AVAILABLE if wanted > 0 else Tractability.WRONG_DIRECTION_ONLY

    return TractabilityResult(gene, tid, required_direction, status,
                              len(mechs), n_act, n_inh, max_phase, actions, drugs)
