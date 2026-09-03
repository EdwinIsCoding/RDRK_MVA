"""Chemoprevention candidate derivation for Track 2, per plan section 7.1.

Why this axis
-------------
Plan section 7.1 ranks cancer chemoprevention and surveillance as "probably the
highest-value output of the whole track", and the proband's presenting event is
rhabdomyosarcoma (HP:0002859). The direct axis, restoring spindle assembly
checkpoint function, is closed with evidence in
``results/summaries/track2_direction_audit.md``: every target nominated by
signed-edge algebra requires activation and no activating drug exists for any of
them. What remains is preventing or detecting the tumours that the chromosomal
instability causes, rather than fixing the instability.

Where the candidates come from, and why not from a model's memory
-----------------------------------------------------------------
CLAUDE.md rule 2 forbids inventing identifiers, and a list of chemoprevention
agents recalled from training is exactly that failure wearing a lab coat. Every
candidate here is derived from a query against a public registry, and every
identifier attached to it comes back from an API call.

The derivation is anchored on the MeSH heading **"Neoplastic Syndromes,
Hereditary"** rather than on a list of syndrome names chosen by us. That matters
for reproducibility: the heading is the registry's own vocabulary, so the query
is a statement about the evidence base rather than about our recall of it.

The transfer assumption, stated once and loudly
-----------------------------------------------
There are **no chemoprevention trials in mosaic variegated aneuploidy**, and
none in rhabdomyosarcoma. Every candidate this module produces is therefore
transferred from a different hereditary cancer predisposition syndrome, with a
different driver gene, a different tumour spectrum and a different mechanism.
That transfer is an assumption and not a small one. It is the dominant
limitation of the axis and it is carried into every output.

No dosing
---------
CLAUDE.md rule 3 forbids dosing, dose ranges and safety margins. Registry
intervention names routinely embed doses ("Atorvastatin 20mg"), so
``strip_dose`` removes them before a name is emitted anywhere, and a test
asserts that no output carries a dose pattern.
"""

from __future__ import annotations

import dataclasses
import json
import pathlib
import re
import time
import urllib.parse
import urllib.request

CTGOV = "https://clinicaltrials.gov/api/v2/studies"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

#: The registry's own MeSH heading for hereditary cancer predisposition.
MESH_HEREDITARY_NEOPLASTIC = "Neoplastic Syndromes, Hereditary"

#: Dose and formulation fragments that appear inside registry intervention
#: names. Removed before any name is emitted. See CLAUDE.md rule 3.
_DOSE = re.compile(
    r"""(
        \b\d+(?:\.\d+)?\s*(?:mg|mcg|µg|ug|g|ml|mL|iu|IU|units?|%)\b(?:\s*/\s*\w+)?
      | \b\d+(?:\.\d+)?\s*(?:mg|mcg|g)\s*/\s*(?:kg|m2|m\^2|day|d)\b
      | \b(?:once|twice|three times)\s+(?:a\s+)?(?:daily|day|weekly)\b
      | \bq\.?d\b | \bb\.?i\.?d\b | \bt\.?i\.?d\b | \bp\.?o\.?\b
      | \b\d+(?:\.\d+)?\s*%
    )""",
    re.IGNORECASE | re.VERBOSE)

#: Route, formulation and trial-arm noise left behind once doses are stripped.
_NOISE = re.compile(
    r"\b(oral|tablet|tablets|capsule|capsules|cream|gel|topical|injection|"
    r"placebo|arm|group|daily|dose|low[- ]dose|high[- ]dose|standard|"
    r"intravenous|iv|subcutaneous|film[- ]coated|extended[- ]release)\b",
    re.IGNORECASE)


def strip_dose(name: str) -> str:
    """Remove dose, frequency and formulation text from an intervention name.

    Registry intervention names are free text written by trial sponsors, so
    ``"Atorvastatin 20mg AND Aspirin 325 mg"`` is a realistic input. Emitting
    that verbatim would put a dose into a project output, which CLAUDE.md rule 3
    forbids outright.
    """
    s = _DOSE.sub(" ", name)
    s = _NOISE.sub(" ", s)
    s = re.sub(r"[()\[\],;:+/]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip(" -")
    return s


def split_combination(name: str) -> list[str]:
    """Split a combination arm into its component agents.

    ``"Atorvastatin AND Aspirin"`` is two candidates, not one compound called
    "atorvastatin and aspirin", and resolving the latter against ChEMBL would
    either fail or, worse, fuzzy-match something unrelated.
    """
    parts = re.split(r"\s+(?:and|plus|\+|with|versus|vs\.?)\s+", name, flags=re.I)
    return [p.strip() for p in parts if len(p.strip()) > 2]


#: Trial-arm wrappers. A sponsor writes "Celecoxib monotherapy" or "Early
#: Vigabatrin" to name an arm, not a compound, and neither string is in any
#: chemical database.
_ARM_WRAPPER = re.compile(
    r"\b(monotherapy|combination|combined|drug|therapy|treatment|early|delayed|"
    r"maintenance|adjuvant|supplement|supplements|nutritional|suppositories|"
    r"suppository|slurry|probiotic|extract|cream|ointment|comparator|"
    r"experimental|active|control|cohort)\b", re.IGNORECASE)

#: Placebo and vehicle arms. These must be **dropped**, never resolved.
#: "no active patidegib" is the control arm of a patidegib trial; resolving it
#: to patidegib would record the placebo arm as evidence for the drug.
_CONTROL_ARM = re.compile(
    r"^\s*(no active\b|vehicle\b|placebo\b|sham\b|matching placebo\b)"
    r"|\bplacebo\b|\bvehicle comparator\b", re.IGNORECASE)

#: Sponsor product codes, e.g. TAVT-18, PTX-022, TPST-1495.
_PRODUCT_CODE = re.compile(r"\b[A-Z]{2,6}[- ]?\d{2,5}\b")

#: Words that are never a compound name, so never worth a lookup on their own.
_STOPWORDS = frozenset({
    "acid", "administration", "antibody", "arm", "black", "dose", "group",
    "high", "low", "oral", "patients", "standard", "study", "then", "with",
})


def is_control_arm(name: str) -> bool:
    """Is this a placebo or vehicle arm rather than an agent?"""
    return bool(_CONTROL_ARM.search(name))


def normalise_candidates(name: str) -> list[str]:
    """Ordered candidate strings to try resolving, most specific first.

    Generating several spellings is safe here **only** because resolution
    demands an exact match against a ChEMBL preferred name or synonym. Each
    candidate is either exactly right or it returns nothing, so a longer list
    widens recall without admitting a near miss.
    """
    out: list[str] = []

    def add(c: str) -> None:
        c = re.sub(r"\s+", " ", c).strip(" -,")
        if c and len(c) > 2 and c.lower() not in {x.lower() for x in out}:
            out.append(c)

    add(name)
    add(_ARM_WRAPPER.sub(" ", name))
    add(_PRODUCT_CODE.sub(" ", _ARM_WRAPPER.sub(" ", name)))

    # A trailing or embedded all-caps abbreviation, as in "Eicosapentanoic Acid
    # EPA" or "mesalamine 5-ASA", where the expansion is the resolvable half.
    stripped = _PRODUCT_CODE.sub(" ", _ARM_WRAPPER.sub(" ", name))
    add(re.sub(r"\b[0-9]?[A-Z]{2,6}(-[A-Z]{2,6})?\b", " ", stripped))

    # Last resort: individual words. Recovers "Aspirin" from the wreckage of
    # "100 for Aspirin 100" once the dose stripper has been through it.
    for tok in re.findall(r"[A-Za-z][A-Za-z-]{3,}", stripped):
        if tok.lower() not in _STOPWORDS:
            add(tok)
    return out


@dataclasses.dataclass(frozen=True)
class TrialEvidence:
    """What the registry says about one agent, with the trials that say it."""

    agent: str
    nct_ids: tuple[str, ...]
    conditions: tuple[str, ...]
    phases: tuple[str, ...]
    #: Primary outcome measures, verbatim. Carried so that any classification of
    #: the trial endpoint can be printed next to the text it was derived from.
    outcomes: tuple[str, ...] = ()

    @property
    def n_trials(self) -> int:
        return len(self.nct_ids)


@dataclasses.dataclass(frozen=True)
class ChemblMolecule:
    chembl_id: str
    pref_name: str
    max_phase: float | None
    atc_codes: tuple[str, ...]
    molecule_type: str | None
    withdrawn: bool | None
    #: Set when the ATC codes were inherited from the parent molecule because
    #: this record is a salt or ester form carrying none of its own.
    atc_inherited_from: str | None = None


class LookupFailed(Exception):
    """A request did not complete, as distinct from completing and finding
    nothing. Collapsing the two makes a transport failure look like a gap in the
    evidence."""


def _get(url: str, timeout: int = 45, attempts: int = 4) -> dict:
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


def _cached(cache: pathlib.Path, key: str, fetch) -> dict | None:
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"{re.sub(r'[^A-Za-z0-9_.-]', '_', key)}.json"
    if f.exists():
        return json.loads(f.read_text())
    try:
        d = fetch()
    except Exception:
        return None
    f.write_text(json.dumps(d))
    time.sleep(0.25)
    return d


def query_trials(advanced: str, cache: pathlib.Path, page_size: int = 200,
                 fields: str = "NCTId,Condition,InterventionType,InterventionName,"
                               "Phase,OverallStatus,StdAge,PrimaryOutcomeMeasure"
                 ) -> dict | None:
    """One ClinicalTrials.gov v2 query, cached by its filter expression."""
    def fetch():
        q = urllib.parse.urlencode({
            "filter.advanced": advanced, "fields": fields,
            "pageSize": page_size, "countTotal": "true", "format": "json"})
        return _get(f"{CTGOV}?{q}")
    # The cache key carries page size and field list. Without them a count
    # query, which asks for one study, poisons the cache for the full query
    # that follows it, and the pipeline silently derives candidates from a
    # single trial while reporting the true total beside it.
    return _cached(cache, f"ct_{advanced}_p{page_size}_{fields}", fetch)


def count_trials(advanced: str, cache: pathlib.Path) -> int:
    d = query_trials(advanced, cache, page_size=1)
    return int((d or {}).get("totalCount") or 0)


def extract_drug_agents(studies: dict) -> dict[str, TrialEvidence]:
    """Collect DRUG interventions from a registry response, keyed by agent name."""
    acc: dict[str, dict] = {}
    for s in studies.get("studies", []):
        ps = s.get("protocolSection", {})
        nct = ps.get("identificationModule", {}).get("nctId")
        conds = tuple(ps.get("conditionsModule", {}).get("conditions", []))
        phases = tuple(ps.get("designModule", {}).get("phases", []))
        outcomes = tuple(
            o.get("measure", "") for o in
            ps.get("outcomesModule", {}).get("primaryOutcomes", [])
            if o.get("measure"))
        for iv in ps.get("armsInterventionsModule", {}).get("interventions", []):
            if iv.get("type") != "DRUG":
                continue
            for agent in split_combination(strip_dose(iv.get("name", ""))):
                key = agent.lower()
                e = acc.setdefault(key, {"agent": agent, "nct": set(),
                                         "cond": set(), "ph": set(), "out": []})
                if nct:
                    e["nct"].add(nct)
                e["cond"].update(conds)
                e["ph"].update(phases)
                for o in outcomes:
                    if o not in e["out"]:
                        e["out"].append(o)
    return {k: TrialEvidence(v["agent"], tuple(sorted(v["nct"])),
                             tuple(sorted(v["cond"])), tuple(sorted(v["ph"])),
                             tuple(v["out"]))
            for k, v in acc.items()}


#: Endpoint classification, applied to a trial's primary outcome text.
#:
#: The distinction that matters for prioritising a chemoprevention hypothesis is
#: not whether a trial exists but what it measured. A trial reporting the number
#: of new lesions is evidence about tumour prevention. A trial reporting Ki-67
#: staining is evidence about a surrogate, and a trial reporting seizure
#: frequency is evidence about something else entirely, even when its
#: participants all have a cancer predisposition syndrome.
#:
#: These are keyword sets of ours, not registry fields. Every classification is
#: printed beside the outcome text it came from so a reader can overrule it.
ENDPOINT_TUMOUR = frozenset({
    "new bcc", "bccs", "basal cell", "polyp", "adenoma", "tumor", "tumour",
    "lesion", "cancer", "carcinoma", "neoplas", "sarcoma", "melanoma",
    "incidence", "number of new", "occurrence of new", "time to progression",
    "tumor burden", "disease progression", "recurrence",
})
ENDPOINT_SURROGATE = frozenset({
    "ki-67", "ki67", "biomarker", "immunohistochem", "staining", "expression",
    "rna-seq", "proliferation", "apoptosis", "caspase", "methylation",
    "concentration", "pharmacokinetic",
})


def classify_endpoint(outcome: str) -> str:
    """Return "tumour", "surrogate" or "other" for one primary outcome string."""
    t = outcome.lower()
    if any(k in t for k in ENDPOINT_TUMOUR):
        return "tumour"
    if any(k in t for k in ENDPOINT_SURROGATE):
        return "surrogate"
    return "other"


def endpoint_classes(ev: TrialEvidence) -> set[str]:
    return {classify_endpoint(o) for o in ev.outcomes} or {"none recorded"}


def _molecule(chembl_id: str, cache: pathlib.Path) -> dict | None:
    def fetch():
        return _get(f"{CHEMBL}/molecule/{chembl_id}?format=json")
    return _cached(cache, f"molrec_{chembl_id}", fetch)


def _exact_molecules(field: str, value: str, cache: pathlib.Path) -> list[dict]:
    """Molecules whose ``field`` equals ``value`` exactly, ignoring case.

    ChEMBL's ``molecule/search`` is a ranked text search and its recall is poor
    for common drug names. Searching it for "sirolimus" returns ten molecules,
    none of them CHEMBL413, whose preferred name **is** SIROLIMUS. Relying on
    that search cost this pipeline roughly half its candidate names.

    The filter endpoints below are exact lookups and return the one right
    answer, so they are tried first and the search is kept only as a fallback.
    """
    def fetch():
        q = urllib.parse.urlencode({f"{field}__iexact": value,
                                    "format": "json", "limit": 20})
        return _get(f"{CHEMBL}/molecule?{q}")

    d = _cached(cache, f"molx_{field}_{value.lower()}", fetch)
    return (d or {}).get("molecules", []) or []


def _to_molecule(m: dict, fallback_name: str, cache: pathlib.Path) -> ChemblMolecule:
    atc = tuple(a for a in (m.get("atc_classifications") or []) if isinstance(a, str))

    # A salt or ester form carries no ATC code of its own, so the safety screen
    # sees nothing to act on. Left alone, that silently passes ERLOTINIB
    # HYDROCHLORIDE while flagging ERLOTINIB, which is the same active molecule
    # and the same concern. Inherit from the parent.
    inherited = None
    parent = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id")
    if not atc and parent and parent != m["molecule_chembl_id"]:
        pm = _molecule(parent, cache) or {}
        patc = tuple(a for a in (pm.get("atc_classifications") or [])
                     if isinstance(a, str))
        if patc:
            atc, inherited = patc, parent

    mp = m.get("max_phase")
    return ChemblMolecule(
        chembl_id=m["molecule_chembl_id"],
        pref_name=m.get("pref_name") or fallback_name,
        max_phase=float(mp) if mp is not None else None,
        atc_codes=atc,
        molecule_type=m.get("molecule_type"),
        withdrawn=m.get("withdrawn_flag"),
        atc_inherited_from=inherited,
    )


def _names_of(m: dict) -> set[str]:
    names = {(m.get("pref_name") or "").lower()}
    for syn in (m.get("molecule_synonyms") or []):
        names.add((syn.get("molecule_synonym") or "").lower())
    return names - {""}


def resolve_molecule(name: str, cache: pathlib.Path) -> ChemblMolecule | None:
    """Resolve one exact agent name to a ChEMBL molecule, or return None.

    Three lookups are tried in order: exact preferred name, exact synonym, then
    the ranked text search. **Every path applies the same exact-name check**, so
    widening recall cannot admit a near miss. A fuzzy search that silently
    returns the nearest molecule is how a wrong ChEMBL identifier gets into a
    report, and a plausible wrong identifier is worse than a gap because a
    reader cannot tell it is wrong (CLAUDE.md rule 2).
    """
    want = name.strip().lower()

    def search():
        q = urllib.parse.urlencode({"q": name, "format": "json", "limit": 10})
        return _get(f"{CHEMBL}/molecule/search?{q}")

    # Lazily, so an exact hit costs one request rather than three. The fuzzy
    # search is the slowest and least reliable of the three and should only ever
    # run when both exact lookups have come back empty.
    sources = (
        lambda: _exact_molecules("pref_name", name, cache),
        lambda: _exact_molecules("molecule_synonyms__molecule_synonym", name, cache),
        lambda: (_cached(cache, f"mol_{want}", search) or {}).get("molecules", []) or [],
    )
    for source in sources:
        for m in source():
            if want in _names_of(m):
                return _to_molecule(m, name, cache)
    return None


def merge_evidence(a: TrialEvidence, b: TrialEvidence) -> TrialEvidence:
    """Union two arms' evidence for the same molecule.

    Once normalisation is doing its job, several intervention names resolve to
    one molecule: "Aspirin", "for Aspirin 300" and "100 for Aspirin 100" are all
    CHEMBL25. Left unmerged they appear as three separate candidates with their
    trial counts split between them, and they can be given **conflicting
    verdicts**, since a name that finds no paediatric trials scores UNKNOWN while
    its twin scores ALLOWED.
    """
    def u(x, y):
        return tuple(sorted(set(x) | set(y)))
    return TrialEvidence(
        agent=a.agent if len(a.agent) <= len(b.agent) else b.agent,
        nct_ids=u(a.nct_ids, b.nct_ids),
        conditions=u(a.conditions, b.conditions),
        phases=u(a.phases, b.phases),
        outcomes=tuple(dict.fromkeys(a.outcomes + b.outcomes)),
    )


def resolve_agent(name: str, cache: pathlib.Path
                  ) -> tuple[ChemblMolecule | None, str | None]:
    """Resolve a registry intervention name, trying normalised spellings.

    Returns the molecule and **the string that actually matched**, so the output
    can print "metformin combination -> metformin -> CHEMBL1431" and a reader can
    audit the normalisation rather than trust it.

    Control arms are refused outright. "no active patidegib" is the placebo arm
    of a patidegib trial, and resolving it would record that arm as evidence for
    the drug.
    """
    if is_control_arm(name):
        return None, None
    for candidate in normalise_candidates(name):
        m = resolve_molecule(candidate, cache)
        if m is not None:
            return m, candidate
    return None, None


def mechanism_actions(chembl_id: str, cache: pathlib.Path) -> tuple[str, ...]:
    """Mechanism-of-action strings for a molecule, for the safety screen to read."""
    def fetch():
        q = urllib.parse.urlencode({"molecule_chembl_id": chembl_id,
                                    "format": "json", "limit": 100})
        return _get(f"{CHEMBL}/mechanism?{q}")

    d = _cached(cache, f"molmech_{chembl_id}", fetch)
    if not d:
        return ()
    out = []
    for m in d.get("mechanisms", []):
        bits = [m.get("mechanism_of_action"), m.get("action_type")]
        t = " ".join(b for b in bits if b)
        if t:
            out.append(t)
    return tuple(sorted(set(out)))


def paediatric_trials(agent: str, cache: pathlib.Path, limit: int = 5) -> tuple[str, ...]:
    """NCT identifiers of trials of this agent that enrolled children.

    ``StdAge`` CHILD is the registry's own standardised age band. This gives
    paediatric *exposure* evidence with a citable identifier, which is what plan
    section 7.4 asks for. It is not evidence of paediatric safety in this
    indication, and the safety module's UNKNOWN verdict exists for that reason.
    """
    expr = (f'AREA[InterventionName]"{agent}" AND AREA[StdAge]CHILD '
            f'AND AREA[InterventionType]DRUG')
    d = query_trials(expr, cache, page_size=limit, fields="NCTId")
    if not d:
        return ()
    out = []
    for s in d.get("studies", []):
        nct = s.get("protocolSection", {}).get("identificationModule", {}).get("nctId")
        if nct:
            out.append(nct)
    return tuple(out[:limit])
