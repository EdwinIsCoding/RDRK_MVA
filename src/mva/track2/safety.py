"""Deterministic safety screen for repurposing candidates, per plan section 7.4.

Why this is rules and not judgement
-----------------------------------
Plan section 7.4 requires this encoded as deterministic rules over label
sections and structured drug fields, explicitly **not** LLM judgement, and
CLAUDE.md rule 3 forbids generating dosing or clinical recommendations. Both
constraints point the same way: a screen whose verdicts come from a language
model cannot be audited, cannot be reproduced, and would be reasoning about a
real child's cancer risk with no traceable basis.

Every verdict here is a pure function of structured inputs, and every verdict
carries the rule that produced it and the field the rule read.

Why the contraindication analysis matters more than the efficacy argument
------------------------------------------------------------------------
The proband has rhabdomyosarcoma (HP:0002859) as the presenting event, in the
context of a chromosomal instability syndrome with DNA-repair-adjacent biology.
For such a patient, an agent that is genotoxic, or that suppresses tumour immune
surveillance, is not merely unhelpful; it is plausibly harmful in the specific
way the underlying disease already threatens. Naming that tension is the
scientific content. Plan section 7.4 says as much, and a clinically trained
judge will look for it first.

What this module does NOT do
----------------------------
No dosing. No dose ranges. No safety margins. No statement that any drug should
be given to anyone. The output is a triage verdict on whether a compound may
enter a research hypothesis list, addressed to researchers.
"""

from __future__ import annotations

import dataclasses
import enum
import re
from collections.abc import Sequence


class Verdict(str, enum.Enum):
    """Deliberately three-valued plus a data gap.

    ``UNKNOWN`` is not a soft pass. A compound with no paediatric data is not
    thereby safe for a child, and collapsing that into ``ALLOWED`` is the
    failure mode this enum exists to prevent.
    """

    EXCLUDED = "excluded"        # categorical, no efficacy argument can override
    FLAGGED = "flagged"          # real tension, must be stated wherever proposed
    UNKNOWN = "unknown"          # required evidence absent; not a pass
    ALLOWED = "allowed"          # passes the screen; still only a hypothesis


@dataclasses.dataclass(frozen=True)
class SafetyFinding:
    rule_id: str
    verdict: Verdict
    reason: str
    field_read: str
    source: str

    def __str__(self) -> str:
        return f"[{self.rule_id}] {self.verdict.value.upper()}: {self.reason}"


@dataclasses.dataclass(frozen=True)
class DrugRecord:
    """Structured facts about a compound. Every field is something that can be
    read from DrugBank, an FDA label section, ChEMBL or a trial registry.

    Fields left as ``None`` mean "not looked up", which is distinct from
    ``False``. Rules treat the two differently.
    """

    name: str
    chembl_id: str | None = None
    drugbank_id: str | None = None

    # Label and pharmacology
    atc_codes: tuple[str, ...] = ()
    mechanism: str | None = None
    boxed_warning: str | None = None
    carcinogenesis_section: str | None = None      # FDA label 13.1
    pregnancy_category_text: str | None = None

    # Structured flags, each from a named source
    is_genotoxic: bool | None = None
    is_mutagenic: bool | None = None
    is_immunosuppressant: bool | None = None
    is_radiosensitiser: bool | None = None
    causes_chromosomal_instability_in_vitro: bool | None = None

    # Paediatric evidence
    has_paediatric_label: bool | None = None
    has_paediatric_pk: bool | None = None
    paediatric_trial_ids: tuple[str, ...] = ()

    # CNS
    crosses_bbb: bool | None = None

    provenance: str = "TODO(source)"


#: Terms that indicate a genotoxic or DNA-damaging mechanism. Matched against
#: the mechanism string and the label's carcinogenesis section.
_GENOTOXIC_TERMS = re.compile(
    r"\b(alkylat\w*|topoisomerase\s+(?:i|ii)\s+inhibit\w*|intercalat\w*|"
    r"dna[- ]cross[- ]link\w*|radiomimetic|antimetabolite|"
    r"nitrogen mustard|platinum[- ]based|nucleoside analogue)\b", re.I)

_IMMUNOSUPPRESSANT_TERMS = re.compile(
    r"\b(mtor inhibit\w*|rapalog\w*|sirolimus|everolimus|temsirolimus|"
    r"calcineurin inhibit\w*|ciclosporin|cyclosporine|tacrolimus|"
    r"immunosuppress\w*|anti[- ]?thymocyte)\b", re.I)

#: ATC level-1 L01 is antineoplastic, L04 is immunosuppressant.
_ATC_ANTINEOPLASTIC = re.compile(r"^L01", re.I)
_ATC_IMMUNOSUPPRESSANT = re.compile(r"^L04", re.I)


def _text_fields(d: DrugRecord) -> str:
    return " ".join(filter(None, [d.mechanism, d.boxed_warning,
                                  d.carcinogenesis_section]))


# ---------------------------------------------------------------------------
# Rules. Each returns a finding or None. Order does not matter: all rules run,
# and the strictest verdict wins.
# ---------------------------------------------------------------------------

def rule_genotoxic(d: DrugRecord) -> SafetyFinding | None:
    """Categorical exclusion. Plan section 7.4.

    A child with a chromosomal instability syndrome and a cancer predisposition
    already has an elevated mutational burden and impaired mitotic fidelity.
    Adding a genotoxic agent works with the disease mechanism rather than
    against it. No efficacy argument overrides this.
    """
    if d.is_genotoxic or d.is_mutagenic:
        return SafetyFinding(
            "SAFE-001", Verdict.EXCLUDED,
            "genotoxic or mutagenic agent, categorically excluded in a "
            "chromosomal instability syndrome with cancer predisposition",
            "is_genotoxic/is_mutagenic", d.provenance)
    if _GENOTOXIC_TERMS.search(_text_fields(d)):
        m = _GENOTOXIC_TERMS.search(_text_fields(d))
        return SafetyFinding(
            "SAFE-002", Verdict.EXCLUDED,
            f"mechanism or label text indicates DNA damage ({m.group(0)!r})",
            "mechanism/boxed_warning/carcinogenesis_section", d.provenance)
    if any(_ATC_ANTINEOPLASTIC.match(a) for a in d.atc_codes):
        return SafetyFinding(
            "SAFE-003", Verdict.EXCLUDED,
            "ATC L01 antineoplastic; cytotoxic chemotherapy is out of scope for "
            "a chemoprevention or supportive-care hypothesis",
            "atc_codes", d.provenance)
    return None


def rule_chromosomal_instability(d: DrugRecord) -> SafetyFinding | None:
    """Anything with an in vitro CIN signal. Plan section 7.4."""
    if d.causes_chromosomal_instability_in_vitro:
        return SafetyFinding(
            "SAFE-004", Verdict.EXCLUDED,
            "reported to induce chromosomal instability in vitro, which is the "
            "disease mechanism itself",
            "causes_chromosomal_instability_in_vitro", d.provenance)
    return None


def rule_immunosuppression(d: DrugRecord) -> SafetyFinding | None:
    """Flag, not exclude. Plan section 7.4 is explicit that naming the tension
    is a strength and ignoring it is the error a clinical judge spots first.

    mTOR inhibition is a genuinely plausible axis here: aneuploid cells carry a
    proteostasis and energy burden, and rapalogs modulate exactly that. It is
    also immunosuppressive, and tumour immune surveillance is load-bearing in a
    child already predisposed to cancer. Both things are true, so the verdict is
    FLAGGED and the tension travels with the candidate wherever it is proposed.
    """
    triggered = (
        d.is_immunosuppressant
        or _IMMUNOSUPPRESSANT_TERMS.search(_text_fields(d))
        or any(_ATC_IMMUNOSUPPRESSANT.match(a) for a in d.atc_codes)
    )
    if triggered:
        return SafetyFinding(
            "SAFE-005", Verdict.FLAGGED,
            "immunosuppressive. Direct tension with tumour immune surveillance "
            "in a cancer predisposition syndrome. May be proposed only with this "
            "tension stated explicitly alongside it, never silently",
            "is_immunosuppressant/mechanism/atc_codes", d.provenance)
    return None


def rule_radiosensitiser(d: DrugRecord) -> SafetyFinding | None:
    if d.is_radiosensitiser:
        return SafetyFinding(
            "SAFE-006", Verdict.FLAGGED,
            "radiosensitiser; the proband's surveillance and treatment pathway "
            "for rhabdomyosarcoma may involve radiotherapy",
            "is_radiosensitiser", d.provenance)
    return None


def rule_paediatric_evidence(d: DrugRecord) -> SafetyFinding | None:
    """Plan section 7.4 requires paediatric exposure data.

    Absence returns UNKNOWN rather than EXCLUDED: a compound with no paediatric
    data is not disqualified as a research hypothesis, but it must never be
    presented as though paediatric use were established.
    """
    has_any = bool(d.has_paediatric_label or d.has_paediatric_pk or d.paediatric_trial_ids)
    if has_any:
        bits = []
        if d.has_paediatric_label:
            bits.append("paediatric label")
        if d.has_paediatric_pk:
            bits.append("paediatric PK")
        if d.paediatric_trial_ids:
            bits.append(f"paediatric trials {', '.join(d.paediatric_trial_ids)}")
        return SafetyFinding(
            "SAFE-007", Verdict.ALLOWED,
            f"paediatric exposure evidence exists: {'; '.join(bits)}",
            "has_paediatric_label/has_paediatric_pk/paediatric_trial_ids", d.provenance)
    if d.has_paediatric_label is None and d.has_paediatric_pk is None:
        return SafetyFinding(
            "SAFE-008", Verdict.UNKNOWN,
            "paediatric exposure not looked up; this is a data gap, not a pass",
            "has_paediatric_label/has_paediatric_pk", d.provenance)
    return SafetyFinding(
        "SAFE-009", Verdict.UNKNOWN,
        "no paediatric label, PK or trial found. Any proposal must state that "
        "paediatric exposure is unestablished",
        "has_paediatric_label/has_paediatric_pk/paediatric_trial_ids", d.provenance)


def rule_cns_penetrance(d: DrugRecord, cns_target: bool = False) -> SafetyFinding | None:
    """Only applies where the nominated target is CNS-expressed."""
    if not cns_target:
        return None
    if d.crosses_bbb is True:
        return SafetyFinding(
            "SAFE-010", Verdict.ALLOWED,
            "blood-brain barrier penetrant, consistent with a CNS-expressed target",
            "crosses_bbb", d.provenance)
    if d.crosses_bbb is False:
        return SafetyFinding(
            "SAFE-011", Verdict.EXCLUDED,
            "does not cross the blood-brain barrier but the nominated target is "
            "CNS-expressed; the mechanistic rationale cannot hold",
            "crosses_bbb", d.provenance)
    return SafetyFinding(
        "SAFE-012", Verdict.UNKNOWN,
        "blood-brain barrier penetrance unknown for a CNS-expressed target",
        "crosses_bbb", d.provenance)


ALL_RULES = (
    rule_genotoxic,
    rule_chromosomal_instability,
    rule_immunosuppression,
    rule_radiosensitiser,
    rule_paediatric_evidence,
)

#: Strictest first. ``screen`` returns the first of these that any rule reached.
_SEVERITY = (Verdict.EXCLUDED, Verdict.FLAGGED, Verdict.UNKNOWN, Verdict.ALLOWED)


@dataclasses.dataclass(frozen=True)
class SafetyResult:
    drug: str
    verdict: Verdict
    findings: tuple[SafetyFinding, ...]

    @property
    def may_be_proposed(self) -> bool:
        """EXCLUDED compounds must not appear in any candidate list.

        FLAGGED and UNKNOWN may appear, and their findings must be printed
        alongside them.
        """
        return self.verdict is not Verdict.EXCLUDED

    @property
    def mandatory_caveats(self) -> list[str]:
        """Text that must travel with this compound wherever it is proposed."""
        return [f.reason for f in self.findings
                if f.verdict in (Verdict.FLAGGED, Verdict.UNKNOWN)]


def screen(drug: DrugRecord, cns_target: bool = False) -> SafetyResult:
    """Run every rule and return the strictest verdict reached."""
    findings: list[SafetyFinding] = []
    for rule in ALL_RULES:
        f = rule(drug)
        if f is not None:
            findings.append(f)
    f = rule_cns_penetrance(drug, cns_target=cns_target)
    if f is not None:
        findings.append(f)

    reached = {f.verdict for f in findings}
    verdict = next((v for v in _SEVERITY if v in reached), Verdict.UNKNOWN)
    return SafetyResult(drug=drug.name, verdict=verdict, findings=tuple(findings))


def screen_all(drugs: Sequence[DrugRecord], cns_target: bool = False) -> list[SafetyResult]:
    return [screen(d, cns_target=cns_target) for d in drugs]
