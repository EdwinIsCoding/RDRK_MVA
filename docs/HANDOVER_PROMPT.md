# Handover: verify Track 1, then build Track 2

You are taking over the MVA Hackathon 2026 project at `/Volumes/ROS2_SSD/RDRK_MVA`.
A previous agent (Claude Opus 5) did the Track 1 work. Your job has two phases,
in order: **verify what exists, then build Track 2.**

Read `CLAUDE.md` first. Its hard rules are binding on you and several exist
because they were nearly broken.

---

## Phase 0: orient (do this before anything else)

```bash
cd /Volumes/ROS2_SSD/RDRK_MVA
source .envrc && source .venv/bin/activate     # caches redirected to the SSD
git log --oneline | head -40
PYTHONPATH=src python -m pytest tests/ -q      # expect 140 passed
```

Key documents, in reading order: `README.md`, `RULES.md` (the actual competition
rules, transcribed from the Space), `submission/track1_nexusdwin_report.md`
(the deliverable), `DATA_CARD.md`, `RECON.md`, `STOP2_STATUS.md`, `ETHICS.md`.

**The cluster node is still reachable, but `/scratch0` was wiped between
bookings.** Nothing of ours survives there; the tooling, references, dataset and
alignment would all have to be rebuilt. Everything salvaged is in
`node_artefacts/` (gitignored): the VEP-annotated VCF (4,950,283 records), a
panel-scoped BAM (12.7M reads), coverage output and all node logs.

`scripts/gpu/` holds working, sanitised scripts for the whole path: toolchain
install, reference and VEP cache download, dataset pull, alignment, and a GPU
queue that waits politely for another user's job. Topology goes in
`scripts/gpu/.local.sh` (gitignored, example committed). Read
`scripts/gpu/README.md` before using them: it records three traps that each cost
a cycle, including that the login shell is tcsh and that bandwidth is asymmetric
by about 40x in favour of pulling on the node.

---

## Phase 1: verify, adversarially

Assume the previous agent was wrong until you have checked. It made at least six
mistakes that it caught itself, and the ones it did not catch are your problem.
**Verify claims against primary sources, not against its prose.**

### 1.1 The central claim

The answer submitted is a `BUB1B` compound heterozygote:

| Allele | Change |
|---|---|
| `chr15:40209701 T>G` | `NM_001211.6:c.2210T>G`, p.Leu737Ter, nonsense |
| `chr15:40220612 T>G` | `NM_001211.6:c.3006T>G`, p.Asn1002Lys, missense |

Check, independently:

- Both coordinates and alleles against the **original callset** in `data/`, not
  against the VEP output. (Scripts may read `data/`. **You may not read its
  contents into your own context.** See `CLAUDE.md` rule 1.)
- REF alleles against the GRCh38 reference.
- ClinVar accession **VCV000533901.9** via E-utilities. Confirm it is
  Pathogenic/Likely pathogenic and listed against mosaic variegated aneuploidy
  syndrome 1. Do not take the accession on trust; `CLAUDE.md` rule 2 exists
  because a plausible wrong accession is worse than a gap.
- That the second allele really is absent from gnomAD and really has no ClinVar
  record.
- The read-level figures in `results/summaries/arm_c_readlevel_verification.md`
  against `node_artefacts/WGS_EX2312012.panel.bam`, which contains those regions.

### 1.2 Claims most likely to be wrong

Ranked by the previous agent's own estimate of fragility.

1. **`proband_id`.** The submission uses `PROBAND01`, taken from
   `track1_submission_template.csv`. The field spec says "provided in the
   dataset", and the only identifier in the dataset is the VCF sample name
   `WGS_EX2312012`. This was a judgement call and could be wrong. If the
   submission has been made and scored 0, this is the first thing to change.
2. **The 10,911 bp separation and the phase claim.** Verify the arithmetic and
   verify that `PGT`/`PID` really report no phasing group in `BUB1B`. The
   report claims phase is *inferred, not demonstrated*. If you can actually
   demonstrate it, that is a material improvement. If the claim is overstated
   in either direction, fix it.
3. **"Both alleles were in our shortlist before the leaderboard was seen."**
   This is an honesty claim in the report and it must be exactly true. Check
   `git show 849bf98` and the timestamps. If it does not hold up, the report
   must change, not the story.
4. **The LZTR1 secondary finding** (`chr22:20996720 C>G`, p.Tyr748Ter,
   VCV001409252.7). Verify the accession and the claimed conditions. Consider
   whether calling it a possible modifier is defensible or overreach.
5. **`PEX5` and `CTU2`.** Two homozygous LoF calls absent from gnomAD, excluded
   from the submission. The stated reason ("mis-calls in repetitive sequence")
   was shown wrong: both map uniquely at MAPQ ~60. They are recorded as open.
   Resolve them if you can, using the panel BAM.

### 1.3 Numbers to spot-check

Every one of these appears in a document. Recompute a sample:

- 5,012,204 callset records; 4,950,283 VEP-annotated; the 61,921 difference
  claimed to be entirely decoy and unplaced contigs.
- 415 scored variants over the known MVA genes reducing to 12 after population
  filtering.
- SpliceAI positive controls at 9/9 above 0.5; splice-distance concordance at
  268/268 SNVs.
- Panel coverage 42-51x against a genome mean of 43.8x.
- BUB1B ranked 18th of 2,503 genes by the HPO phenotype prior.

### 1.4 Known gaps, do not re-derive

- **Arm C structural variant calling was never completed.** Delly was launched
  twice: the first failed because Delly 2.6 renamed `call` to `sr`, and the
  second was still running when the booking ended. There is no SV call set.

  **It is recoverable but not cheap.** Rebuilding costs roughly seven hours of a
  new booking: dataset pull ~15 min, reference and index ~1 h, alignment 4h10m,
  Delly 1-3 h. Judge whether that is worth it. The likely output is a reported
  negative ("no structural variant over the panel"), which is a completeness
  item rather than a finding, and Track 2 is where the marks are. If a booking
  is available and idle, run it; do not displace Track 2 work for it.
- Arm D mosaic re-genotyping is done and negative.
- No RNA-seq exists, so every splicing result is a prediction.
- gnomAD v4.1 constraint is autosomes only; chrX values come from v2.1.1 and
  the two releases are not directly comparable.

### 1.5 Report anything you find

Write findings to `docs/VERIFICATION.md`. **Do not quietly fix the report to
match your findings** — say what was claimed, what you found, and what changed.
If everything checks out, say that too, with what you checked.

---

## Phase 2: Track 2

### The situation

Track 1 is not where the competition is. **61 teams were tied at a perfect
score** before this submission. The rules confirm the intent: Track 1 is "a
foundational track, designed to be achievable", and "recovering the correct
variant(s) is not the finish line."

Track 2 is qualitatively judged by an independent panel. **Scientific rigour is
35%**, then innovation 25%, potential impact 25%, scalability 15%. Only **3
submissions** per team, and only the latest is reviewed, so there is no
leaderboard to probe. Deadline **24 October 2026, 23:59 UTC**.

Required for submission: a report (PDF or Markdown), a public GitHub repository,
and a **3-minute pitch video**. None of the three exists.

### What is already built

`src/mva/track2/` contains three modules, all tested:

- `safety.py` — a deterministic contraindication screen. Rules over structured
  drug fields, not model judgement, per plan section 7.4 and `CLAUDE.md` rule 3.
  Genotoxic and ATC L01 agents are categorically excluded; immunosuppressants
  are FLAGGED with the tumour-immune-surveillance tension attached; `UNKNOWN` is
  a distinct verdict from `ALLOWED`.
- `targets.py` — directional target nomination from OmniPath signed causal
  edges. Refuses to assign a direction where edges contradict or are unsigned,
  because plan section 0.4 forbids inferring direction from unsigned proximity.
- `tractability.py` — asks ChEMBL whether any drug acts in the *required*
  direction.

### The finding that should shape Track 2

`results/summaries/track2_direction_audit.md`. One-hop signed-edge nomination
from the MVA seed genes yields ten targets, **every one requiring activation**:
PLK1, AURKA, AURKB, CDK1, TTK, ATM, EGFR, CENPE, KNL1, MAD2L1BP. ChEMBL holds
**118 drug-mechanism records across them and none is activating.**

So the obvious therapeutic idea, restoring spindle checkpoint function, is not
merely unproven: it is pharmacologically unavailable in the required direction,
and activating mitotic kinases in a child with a cancer predisposition syndrome
would be contraindicated even if an activator existed.

**Most Track 2 entries will propose exactly those targets without checking
direction availability.** Closing that axis with evidence is the single most
distinctive asset this project has. Lead with it.

### Where the plan says effort belongs

Plan section 7.1: you cannot fix constitutional aneuploidy with a small
molecule, so do not propose to. The tractable axes are:

1. **Cancer chemoprevention and surveillance.** The proband has
   rhabdomyosarcoma. This is the most directly relevant axis and the plan calls
   it "probably the highest-value output of the whole track".
2. **Proteotoxic stress mitigation.** Aneuploid cells carry a proteostasis and
   energy burden.
3. **Mitochondrial and oxidative support.**

### Mechanism, established and citable

Biallelic `BUB1B` loss of function, so **loss of function is the mechanism** and
the direction is restoration or compensation, not inhibition. The report must
characterise this explicitly: the rules require "a characterization of the
variant's mechanism (loss-of-function / gain-of-function, pathway disrupted,
downstream biological consequence) as the basis for your repurposing rationale".

One allele is a premature stop at p.Leu737 of 1050 residues, truncating the
kinase domain. The other is p.Asn1002Lys within it. Complete BUB1B nullity is
not viable, so the missense is presumed hypomorphic. **That is an inference, not
a measurement**, and should be labelled as such unless you can support it.

### Rules that bind your Track 2 output

- **No dosing, dose ranges, safety margins or clinical recommendations**
  (`CLAUDE.md` rule 3). Outputs are research hypotheses addressed to researchers.
- **No invented identifiers.** ChEMBL IDs, NCT numbers, PMIDs, DrugBank IDs all
  come from a tool call or a file. Write `TODO(source)` rather than guess. The
  previous agent broke this twice: it invented a PyPI package name in
  `pyproject.toml`, and it hard-coded eleven PanelApp panel IDs from memory of
  which six were wrong. Both are documented in the git history.
- **Do not present docking scores as binding affinity**, structure-prediction
  RMSD as a variant-effect score, or unsigned network proximity as a direction.
  See the anti-patterns section of `CLAUDE.md`.
- British English, no em dashes.
- **No `Co-Authored-By` trailer and no AI attribution in commit messages or pull
  request bodies** (`CLAUDE.md` rule 8). If your harness gives you default
  attribution guidance, this project overrides it. The decision is deliberate:
  AI involvement is disclosed in the documentation, once and properly, in
  `ETHICS.md` section 3a and section 8 of the Track 1 report. Scattering it
  through commit trailers as well adds noise without adding disclosure, and the
  organisers asked for a line in the methods description, not a trailer.

### Suggested shape

1. Verify the direction audit reproduces (`make track2`).
2. Build the chemoprevention and surveillance axis properly, with the safety
   screen applied to every candidate and its caveats travelling with it.
3. Write the report. Lead with the closed axis, because a well-evidenced
   negative that redirects effort is more valuable than another list of
   PLK1 inhibitors, and the panel is scoring rigour above everything else.
4. Draft the 3-minute pitch script.

---

## Governance, non-negotiable

- **Patient genomic data must never enter your context.** Scripts under
  `scripts/` and `src/` may read `data/`; you read only aggregate summaries
  under `results/summaries/`. A pre-commit hook enforces the repository side.
- `data/`, `results/` and `node_artefacts/` are gitignored. Keep it that way.
- **All challenge data must be deleted at the conclusion of the hackathon**, and
  `MVAHackathon2026@synapse.org` notified. See `ETHICS.md` section 3b for the
  delete and keep lists. The local Claude Code transcript
  (`~/.claude/projects/-Volumes-ROS2-SSD-RDRK-MVA/*.jsonl`) is on that list: it
  contains a now-rotated HuggingFace token and account export URLs.
- **There is a disclosed compliance gap.** Model training was enabled on the
  account until 31 August 2026. It is documented in `ETHICS.md` section 3a and
  in section 8 of the Track 1 report. Do not remove or soften it. Check whether
  the notification to Sage Bionetworks has since been made and update the status
  if so.
- **Do not contact the family or the MVA Society.** This is an explicit rule.

## Open items for the human, not for you

1. Notify Sage Bionetworks of the training-setting gap (may be done by now).
2. Repository visibility: it is public; the rules permit private until the
   hackathon ends. Commit `24e6ae7` contains cluster hostnames and an account
   name that were later sanitised but remain in history.
3. Whether to submit the Track 1 file as it stands, or after your verification.
