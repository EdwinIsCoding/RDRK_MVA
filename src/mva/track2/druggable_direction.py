"""Is a whole therapeutic axis available in the direction its biology requires?

The generalisation this module makes
------------------------------------
``tractability.py`` asks the question one target at a time: does a drug exist
that acts on this protein in the required direction? Running it over the targets
nominated for the direct spindle-checkpoint axis produced the project's most
distinctive Track 2 result. Ten targets, every one requiring activation, 118
drug-mechanism records between them, none activating.

That result raises an obvious question about the *other* axes. Proteotoxic
stress mitigation and mitochondrial support both require increasing the activity
of something, exactly as checkpoint restoration did. If activating drugs are
scarce in general, then the direct axis is not distinctively closed and the
finding means much less than it appears to.

So this module measures the denominator. It builds the set of human proteins for
which ChEMBL records **any** activating drug mechanism, then asks what fraction
of each axis's gene set falls inside it. The direct axis becomes a control on
the method rather than a claim standing on its own.

Why it is one query rather than thousands
-----------------------------------------
Asking ChEMBL per gene costs two calls per gene and there are thousands of genes
across these axes. Inverting it is cheap: there are only about 1,200 activating
mechanism records in all of ChEMBL, so fetching them all and reading off their
targets costs a handful of calls and gives the same answer.
"""

from __future__ import annotations

import dataclasses
import json
import pathlib
import time
import urllib.parse
import urllib.request

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
QUICKGO = "https://www.ebi.ac.uk/QuickGO/services"

from mva.track2.tractability import ACTIVATING_ACTIONS, INHIBITING_ACTIONS

#: GO evidence codes that rest on an experiment rather than on electronic
#: inference. IEA annotations are unreviewed and dominate large GO terms, so
#: results are reported both ways rather than silently pooled.
EXPERIMENTAL_EVIDENCE = frozenset({"EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                   "HTP", "HDA", "HMP", "HGI", "HEP"})


def binom_zero(p: float, n: int) -> float:
    """Probability of observing zero successes in ``n`` independent draws at
    rate ``p``.

    Used to ask whether "no activating drug across ten targets" is surprising.
    At the measured genome-wide activation rate it is not, which is why this
    function exists in the library rather than inline in a script: the number it
    produces weakens the project's headline Track 2 claim, so it is tested.
    """
    return (1.0 - p) ** n


def _get(url: str, accept: str = "application/json", timeout: int = 120) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "mva-hackathon-2026",
                                               "Accept": accept})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.load(fh)


def _cached(cache: pathlib.Path, key: str, fetch):
    cache.mkdir(parents=True, exist_ok=True)
    import re as _re
    f = cache / f"{_re.sub(r'[^A-Za-z0-9_.:-]', '_', key)}.json"
    if f.exists():
        return json.loads(f.read_text())
    d = fetch()
    f.write_text(json.dumps(d))
    time.sleep(0.2)
    return d


# ---------------------------------------------------------------------------
# ChEMBL: which proteins can be pushed which way
# ---------------------------------------------------------------------------

def _mechanisms_for_action(action: str, cache: pathlib.Path) -> list[dict]:
    out: list[dict] = []
    offset = 0
    while True:
        def fetch(offset=offset):
            q = urllib.parse.urlencode({"action_type": action, "format": "json",
                                        "limit": 1000, "offset": offset})
            return _get(f"{CHEMBL}/mechanism?{q}")
        d = _cached(cache, f"mech_action_{action}_{offset}", fetch)
        got = d.get("mechanisms", [])
        out.extend(got)
        total = d.get("page_meta", {}).get("total_count", len(out))
        offset += len(got)
        if not got or offset >= total:
            break
    return out


def _target_symbols(target_ids: list[str], cache: pathlib.Path) -> dict[str, set[str]]:
    """Map ChEMBL target identifiers to their gene symbols, in bulk."""
    out: dict[str, set[str]] = {}
    for i in range(0, len(target_ids), 40):
        chunk = target_ids[i:i + 40]
        def fetch(chunk=chunk):
            q = urllib.parse.urlencode({"target_chembl_id__in": ",".join(chunk),
                                        "format": "json", "limit": 100})
            return _get(f"{CHEMBL}/target?{q}")
        d = _cached(cache, f"targets_{i}_{chunk[0]}_{len(chunk)}", fetch)
        for t in d.get("targets", []):
            syms = {s.get("component_synonym")
                    for comp in t.get("target_components", [])
                    for s in comp.get("target_component_synonyms", [])
                    if s.get("syn_type") == "GENE_SYMBOL" and s.get("component_synonym")}
            out[t["target_chembl_id"]] = syms
    return out


@dataclasses.dataclass(frozen=True)
class DirectionalProteome:
    """Gene symbols reachable by a drug in each direction, with the evidence."""

    activatable: dict[str, set[str]]      # gene -> ChEMBL molecule ids
    inhibitable: dict[str, set[str]]
    n_activating_mechanisms: int
    n_inhibiting_mechanisms: int


def build_directional_proteome(cache: str | pathlib.Path) -> DirectionalProteome:
    cache = pathlib.Path(cache)
    both: dict[str, dict[str, set[str]]] = {"act": {}, "inh": {}}
    counts = {"act": 0, "inh": 0}

    for bucket, actions in (("act", ACTIVATING_ACTIONS), ("inh", INHIBITING_ACTIONS)):
        mechs: list[dict] = []
        for action in sorted(actions):
            mechs.extend(_mechanisms_for_action(action, cache))
        counts[bucket] = len(mechs)
        tids = sorted({m["target_chembl_id"] for m in mechs if m.get("target_chembl_id")})
        sym = _target_symbols(tids, cache)
        for m in mechs:
            for g in sym.get(m.get("target_chembl_id") or "", ()):
                both[bucket].setdefault(g, set()).add(m.get("molecule_chembl_id") or "")

    return DirectionalProteome(both["act"], both["inh"],
                               counts["act"], counts["inh"])


# ---------------------------------------------------------------------------
# QuickGO: what is in an axis
# ---------------------------------------------------------------------------

def resolve_go_term(name: str, cache: str | pathlib.Path) -> str | None:
    """Resolve a GO term name to its identifier by exact name match.

    An inexact match is refused. CLAUDE.md rule 2: a GO identifier that looks
    right and names a different term is worse than a gap.
    """
    cache = pathlib.Path(cache)
    def fetch():
        q = urllib.parse.urlencode({"query": name, "limit": 10, "page": 1})
        return _get(f"{QUICKGO}/ontology/go/search?{q}")
    d = _cached(cache, f"go_search_{name}", fetch)
    for r in d.get("results", []):
        if (r.get("name") or "").strip().lower() == name.strip().lower():
            return r.get("id")
    return None


def go_gene_symbols(go_id: str, cache: str | pathlib.Path,
                    experimental_only: bool = False) -> set[str]:
    cache = pathlib.Path(cache)
    syms: set[str] = set()
    page = 1
    while True:
        def fetch(page=page):
            q = urllib.parse.urlencode({"goId": go_id, "taxonId": 9606,
                                        "geneProductType": "protein",
                                        "limit": 200, "page": page})
            return _get(f"{QUICKGO}/annotation/search?{q}")
        d = _cached(cache, f"go_ann_{go_id}_{page}", fetch)
        for r in d.get("results", []):
            if experimental_only and r.get("goEvidence") not in EXPERIMENTAL_EVIDENCE:
                continue
            if r.get("symbol"):
                syms.add(r["symbol"])
        info = d.get("pageInfo") or {}
        if page >= int(info.get("total") or 1):
            break
        page += 1
    return syms
