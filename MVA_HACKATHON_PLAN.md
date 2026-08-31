# MVA Hackathon 2026 — Build Plan (executable)

**Target:** Sage Bionetworks / MVA Society / HuggingFace / BEACON — *Rare Disease, Real Kid: The MVA Hackathon 2026*
**Deadline:** 24 October 2026. Winners 25 November. **~8 weeks from 30 Aug 2026.**
**Hardware:** 1× RTX 6000 (48 GB VRAM), local. Data already downloaded (~85 GB, 11 files).
**Tracks:** T1 — predict causal variant(s). T2 — drug repurposing against disrupted biology.
**Licence:** all outputs CC BY 4.0.

> **How to use this document.** Phase 0 is a script. Run it first. It writes `DATA_CARD.md`,
> which determines which of the six Track 1 arms are possible. Everything after §4 is
> conditional on that output. There are four **STOP** checkpoints where the agent must
> halt and report rather than proceed on assumptions.

---

## 0. Strategic read

### 0.1 What MVA is

Mosaic Variegated Aneuploidy: autosomal-recessive chromosomal instability syndrome.
≥25% of cells carry an abnormal chromosome count (mostly trisomies/monosomies, varying
cell-to-cell and tissue-to-tissue). Failure of the mitotic spindle assembly checkpoint
(SAC), centrosome dynamics, or chromosome cohesion.

| Gene | Locus | Notes |
|---|---|---|
| `BUB1B` | MVA1, OMIM 257300 | BubR1, core SAC. Most common. **Often 1 coding + 1 cryptic allele.** |
| `CEP57` | MVA2 | Microtubule nucleation/stabilisation. |
| `TRIP13` | MVA3 | Strong Wilms tumour predisposition. |
| `BUB1`, `BUB3` | — | SAC components, fewer cases. |
| `CEP192` | — | Recent; MVA + tetraploidy + male infertility. |
| `SMC5` | Atelis syndrome | Near-tetraploidy/MVA overlap, DNA repair. |

Phenotype: IUGR, postnatal growth failure, microcephaly (often severe), CNS anomalies,
seizures, DD/ID, and **childhood cancer predisposition** (rhabdomyosarcoma, Wilms, ALL).
The cancer risk is a hard constraint on Track 2.

### 0.2 The load-bearing inference

**This case is almost certainly not solved by a standard coding-variant pipeline.**

The family has clinical sequencing. A clean biallelic LoF pair in `BUB1B` would have been
reported. The hackathon exists because that did not happen.

> A submission that runs VEP + CADD + AlphaMissense over coding SNVs and reports
> "no pathogenic variant in known MVA genes" reproduces the clinical result and wins nothing.

Hypothesis classes, in rough prior order:

1. **Cryptic second allele in a known gene** — deep intronic, branch point, polypyrimidine
   tract, UTR, promoter, or synonymous creating a cryptic splice site. Published MVA1
   families show exactly this pattern. **Highest prior. Weight time here.**
2. **Structural variant** — deletion, duplication, inversion, MEI, complex rearrangement.
   Short-read SV calling is poor; diagnostic labs under-call it.
3. **Novel gene** in mitotic machinery — kinetochore, SAC, centrosome, cohesin/condensin,
   SMC5/6, APC/C.
4. **Mosaic / low-VAF** — MVA *is* a mosaicism disorder; standard callers filter VAF < 0.2.
5. **UPD / imprinting / repeat expansion** — cheap to exclude, occasionally decisive.
6. **Non-coding regulatory** — enhancer disruption affecting a mitotic gene.

### 0.3 What wins

Correctness, **reproducibility**, evidence quality, honest uncertainty. Not complexity.

- Every claim traceable to a resolvable identifier (rsID, ClinVar VCV, PMID, ChEMBL, NCT).
- Container + lockfile reproducing the ranked list bit-for-bit.
- Calibrated ranked lists with stated confidence, not one confident answer.
- Negative results reported. "We excluded repeat expansions, here is the evidence" earns
  real credit and costs a day.
- T2 must be mechanistically tied to T1 output, with safety contraindications explicit.

### 0.4 Anti-goals — do not implement, even if asked

- ❌ **Structure-prediction WT-vs-mutant RMSD as a variant-effect score.** Single-residue
  substitutions barely perturb the evolutionary signal ESMFold/AF2 depend on. You get
  prediction noise and read it as biology. Use FoldX / Rosetta `cartesian_ddg` /
  ThermoMPNN / RaSP, and state their r≈0.5–0.7 accuracy honestly.
- ❌ **DiffDock confidence as binding free energy.** Pose predictor only. PoseBusters showed
  a large fraction of its poses fail physical validity, and it degrades off-distribution
  from PDBBind — which is where you'd be.
- ❌ **Docking scores as affinity rankings.** ~2 kcal/mol error. Enrichment tool only.
- ❌ **Unsigned network proximity to infer drug direction.** PrimeKG/Hetionet edges are
  largely undirected and unsigned. Proximity gives *which* protein, never *inhibit or
  activate*. Use SIGNOR / OmniPath / Reactome for sign.
- ❌ **LLM-generated dosing.** See §8.
- ❌ **KG "discoveries" that are already edges in the KG.** Time-split or it's leakage.

---

## 1. Governance — before the first commit

Real child's genome, specific consent, real family reading the outputs.

- [ ] Re-read the **Hackathon Rules** you accepted at the gate. Extract clauses on
      redistribution, derived data, third-party API use, publication. Put them in `RULES.md`.
- [ ] **Decide explicitly whether variant-level data may enter a hosted LLM API.** Claude
      Code transmits file contents to Anthropic. The scaffold below assumes it may not:
      `data/` is agent-inaccessible and the agent works from aggregate summaries in
      `results/summaries/`. If the rules permit more, relax it deliberately, not by accident.
- [ ] `.gitignore` excluding `data/`, `*.vcf*`, `*.bam`, `*.cram`, `*.fastq*` **before**
      `git init` and the first commit.
- [ ] Pre-commit hook that hard-fails on genomic extensions and files > 5 MB.
- [ ] `PROVENANCE.md`: SHA256 of every input file, HF dataset revision hash, tool versions,
      DB snapshot dates.
- [ ] `ETHICS.md`: consent scope as understood, what you did and did not do, explicit
      statement that outputs are research hypotheses and not clinical advice. Ship it.

---

## 2. PHASE 0 — Recon (Day 1). Run this before anything else.

Data is local. This phase reads it, characterises it, and writes the facts that route
the rest of the plan.

### 2.1 `scripts/00_inventory.sh`

See the script in `scripts/`. It writes a file manifest with SHA256 checksums, detects
real file types, and dumps VCF headers, sample lists, contig tables and `bcftools stats`
into `results/recon/`. BAM/CRAM headers and `idxstats` are captured when alignments exist.
Tabular and clinical files have their structure captured without their contents.

### 2.2 `scripts/01_characterise.py`

Answers these and writes each into `DATA_CARD.md`. **Every one routes a downstream decision.**

| Question | Method | Routes |
|---|---|---|
| **Trio, quad, or singleton?** | `bcftools query -l`; parental sample IDs | Inheritance filtering (§6.1). No trio → de novo and segregation arms die. |
| **Reference build?** | VCF `##contig` lengths: chr1 = 249,250,621 (GRCh37) vs 248,956,422 (GRCh38) | **Everything.** Write to `config/config.yaml`. Fail loudly on mismatch. |
| **BAM/CRAM present?** | `filetypes.tsv` | Arms C (SV) and D (mosaic). No alignments → both die; document as a limitation. |
| **RNA-seq present?** | filetypes, read groups, or a counts matrix | **Highest-value single question.** Enables FRASER2/OUTRIDER functional confirmation (§6.2) *and* T2 signature reversal (§7.3). |
| **How many tissues?** | sample IDs, read groups, clinical file | VAF–aneuploidy correlation (§6.4). |
| **Caller and pre-filtering?** | `##source`, `##FILTER`, `##FORMAT` | Tells you what was already discarded — i.e. where to look. |
| **Coverage over MVA genes?** | `mosdepth` over the 7 known genes | Any region < 10× is an immediate lead — homozygous deletions present as absence. |
| **Clinical data format?** | head of tabular/JSON | HPO available for Exomiser/LIRICAL, or must be extracted from free text. |
| **Karyotype / aneuploidy %?** | clinical files | Ground truth for §6.4 correlation. |
| **Does clinical data already name a candidate het?** | clinical summary | **If yes, the project becomes "find the second allele" and §6 restructures.** |

### 2.3 🛑 STOP #1

Write `DATA_CARD.md` and `RECON.md`. `RECON.md` states, per hypothesis class in §0.2,
whether it is testable with this data, and drops the rest with a reason.

**Halt and report.** Do not proceed to Phase 1 until the branch is confirmed.

### 2.4 Branch table

| Branch | Condition | Active arms | Emphasis |
|---|---|---|---|
| **A** | Trio + BAMs + RNA-seq | All six | Arm B confirmed functionally by FRASER2. Best case; T2 gets real signature reversal. |
| **B** | Trio + BAMs, no RNA | A, B, C, D, E, F | Splicing predictions stay predictions. Weight SV and mosaic harder. |
| **C** | Singleton + BAMs | A, B, C, D, E, F (no segregation) | Cannot phase compound hets by inheritance — use read-backed + population phasing. |
| **D** | VCF only | A, B, E, F | No SV calling, no mosaic re-genotyping. Say so prominently. Splicing arm becomes ~everything. |
| **E** | Clinical data already names a het in a known MVA gene | Restructure | Collapse to one question: find the *trans* allele. Deep-dive that gene's non-coding space, SVs, and RNA. |

---

## 3. Repo scaffold

```
mva-hackathon-2026/
├── CLAUDE.md                 # agent contract — §3.1
├── README.md  RULES.md  ETHICS.md  PROVENANCE.md
├── DATA_CARD.md  RECON.md    # written by Phase 0
├── pyproject.toml            # uv, fully pinned
├── environment.yml           # bioconda: bcftools samtools vep spliceai manta delly ...
├── Dockerfile
├── Snakefile
├── config/
│   ├── config.yaml           # build, branch, sample map — written by Phase 0
│   ├── db_versions.yaml
│   └── gene_panels/
│       ├── mva_known.tsv
│       ├── mitotic_extended.tsv
│       └── panel_provenance.md
├── scripts/
│   ├── 00_inventory.sh
│   └── 01_characterise.py
├── src/mva/
│   ├── io/  track1/  track2/  llm/  report/
├── tests/
│   ├── test_positive_controls.py   # §5
│   └── fixtures/
├── benchmarks/published_mva_variants.tsv
└── results/
    ├── recon/
    └── summaries/            # ONLY directory the agent may read freely
```

### 3.1 `CLAUDE.md`

See `CLAUDE.md` in the repository root. Hard rules: no patient genomic data in agent
context; no invented identifiers; no dosing; every coordinate build-tagged; all randomness
seeded; halt at every STOP checkpoint.

---

## 4. Phase 1 — Environment and reproducibility (Day 2)

Build the container now, not in week 7. Verify tool versions inside it. Snapshot and
record: VEP cache version, gnomAD release, ClinVar release date, SpliceAI model weights
hash, Exomiser data version → `db_versions.yaml` and `PROVENANCE.md`.

GPU check: confirm CUDA is visible inside the container and SpliceAI/Pangolin run on the
RTX 6000. The 48 GB budget is ample for these; the constraint is wall-clock, not VRAM.

---

## 5. Phase 2 — Positive control harness (Days 3–6)

**Before the real pipeline.** Without it you cannot know your method works, and validation
is a large share of the judging.

### 5.1 Benchmark

`benchmarks/published_mva_variants.tsv` — every published MVA-causal variant with gene,
HGVS c./p., build-tagged coordinates, zygosity, variant class, PMID, and **whether the
source paper reported that a standard pipeline missed it.** Cover `BUB1B` (especially the
cryptic-splice-allele families), `CEP57`, `TRIP13`, `BUB1`, `BUB3`, `CEP192`, `SMC5`.

### 5.2 Spike-in test

```python
@pytest.mark.parametrize("variant", load_benchmark())
def test_pipeline_recovers_known_mva_variant(variant, pipeline):
    spiked = spike(background_vcf, variant)     # background = 1000G sample, NOT the proband
    ranked = pipeline.run(spiked, hpo=proband_hpo)
    assert variant.gene in [r.gene for r in ranked[:20]], (
        f"MISS: {variant.gene} {variant.hgvs} ({variant.pmid}) class={variant.variant_class}"
    )
```

Report **recall@10/@20/@50 broken down by variant class**. If you recover coding LoF and
miss every deep-intronic positive, you have learned the most important fact about your own
method with five weeks left to act on it.

### 5.3 Leakage control

ClinVar contains these variants. Run twice: enabled, and with all known-MVA-gene ClinVar
records masked. **The masked recall is the only number relevant to finding a novel cause.**
Report both.

### 5.4 🛑 STOP #2

Report the recall table. If masked recall@20 for splice-class positives is below ~50%, the
splicing arm needs rework before it touches the real data.

---

## 6. Phase 3 — Track 1 (Weeks 2–5)

Arms run in parallel, each emitting `Candidate(gene, variant, class, score, evidence[])`
into a common table. Activate per the §2.4 branch.

### 6.1 Arm A — Baseline (must exist, won't win)

VEP `--everything` + gnomAD v4 + ClinVar + AlphaMissense + CADD + REVEL.
Filter gnomAD popmax AF < 1e-4. Inheritance models: biallelic (primary — MVA is AR),
de novo (trio only), X-linked (low prior), **mosaic (do not VAF-filter; see Arm D)**.

Expected outcome: reproduces the clinical negative. **State this explicitly in the report** —
reproducing a known negative is evidence your pipeline is calibrated, not a failure.

### 6.2 Arm B — Splicing and the cryptic second allele ⭐ *highest prior*

- **SpliceAI** + **Pangolin**, `-D 500` minimum. The default ±50 bp window is precisely
  why these variants get missed clinically.
- Scope to introns/UTRs/2 kb promoters of the extended panel (§6.5). Genome-wide at ±500 bp
  is days of compute and mostly noise.
- **Interrogate `BUB1B` exhaustively**: every intronic, synonymous and UTR variant, scored
  for cryptic donor/acceptor gain, branch point disruption (**BPP**, **LaBranchoR**),
  polypyrimidine tract weakening, uORF creation (**UTRannotator**). Repeat for `CEP57`,
  `TRIP13`, `CEP192`, `SMC5`.
- If Arm A found a single strong coding allele in any MVA gene, **the whole project becomes
  "where is the other allele"** — redirect all remaining effort to that gene.
- **If RNA-seq exists (branch A):** run **FRASER2** (aberrant splicing) and **OUTRIDER**
  (expression outliers). A predicted cryptic splice site *plus* an observed aberrant
  junction in the same gene is a publishable answer and very likely wins the track.

### 6.3 Arm C — Structural variants *(branches A/B/C only)*

**Manta**, **Delly**, **CNVnator**/**cn.mops**, **ExpansionHunter**, **MELT**/**scramble**.
Merge with **SURVIVOR** or **Jasmine**; require ≥2-caller support for the confidence tier.
Annotate with **AnnotSV**. Priority: the mitotic panel, any gene where Arm A found one het
allele, and the low-coverage regions flagged in Phase 0.

### 6.4 Arm D — Mosaic / low-VAF *(branches A/B/C only)*

MVA is a mosaicism disorder; germline callers discard the signal.
Re-genotype the panel with **Mutect2** tumour-only or **DeepSomatic**, or force-call every
panel position with `bcftools mpileup` and inspect raw allele fractions.
If multiple tissues: **correlate candidate VAF against the reported aneuploidy percentage
per tissue.** A correlation is strong, distinctive evidence few teams will have.
Discriminate true mosaicism from artefact: strand bias, position bias, segdup/low-complexity
overlap, mappability.

### 6.5 Arm E — Novel gene, phenotype-driven

`config/gene_panels/mitotic_extended.tsv` (~300–500 genes) from:
GO (SAC GO:0007094, kinetochore GO:0000776, centrosome cycle GO:0007098, cohesion
GO:0007062, chromosome segregation GO:0007059); Reactome (Mitotic Spindle Checkpoint,
Cohesin Loading, Centrosome Maturation); CORUM (MCC, APC/C, cohesin, condensin, SMC5/6,
γ-TuRC, Ndc80, Ska); high-confidence STRING/BioGRID first-shell interactors of the seven
known genes; MGI aneuploidy/CIN phenotypes. Annotate gnomAD LOEUF and pLI.

Phenotype prioritisation with proband HPO: **Exomiser** (hiPHIVE), **LIRICAL**, **AMELIE**
(literature-driven; surfaces things ontology tools miss). Intersect the union of top hits
with the mitotic panel.

### 6.6 Arm F — Completeness checks (half a day each, all become reported negatives)

UPD/ROH (`bcftools roh`, AutoMap) · repeat expansions (ExpansionHunter Denovo) ·
mtDNA variants and heteroplasmy · recurrent pathogenic CNV loci.

### 6.7 Evidence integration

**Do not train a black-box ranker on n=1.** Transparent additive framework, ACMG/AMP-shaped:

```
score = w_gene_plausibility   # panel membership, constraint, tissue expression
      + w_variant_effect      # AlphaMissense / SpliceAI / LoF class
      + w_inheritance         # segregation consistency
      + w_frequency           # gnomAD rarity
      + w_functional          # RNA-seq corroboration — weight heavily if available
      - w_penalty             # artefact indicators, low mappability
```

Weights fixed *a priori* from literature and **frozen before the proband is scored**.
Tune on the §5 benchmark only. Document and justify each weight.

Output: ranked table with **ACMG-style classification naming the specific criteria invoked**
(PVS1, PM2, PP3, …), evidence per criterion, a confidence, and **the experiment that would
falsify it**.

### 6.8 🛑 STOP #3 (end of week 4)

Report the T1 ranked list v1 and arm-by-arm yield. Confirm the T2 target before week 5.

---

## 7. Phase 4 — Track 2 (Weeks 4–7, overlapping)

### 7.1 Reframe the therapeutic question

You cannot fix constitutional aneuploidy with a small molecule. Do not propose to.

| Axis | Rationale |
|---|---|
| **Proteotoxic stress mitigation** | Aneuploid cells carry unbalanced protein stoichiometry → chronic proteostasis burden. Recent work links proteostasis failure and mitochondrial dysfunction to CIN-induced microcephaly. HSF1 activators, autophagy inducers, chaperone modulators. |
| **mTOR / autophagy modulation** | Compensatory autophagy under aneuploidy-induced energy stress. But see §7.4 — immunosuppression in a cancer-predisposed child. |
| **Mitochondrial / oxidative support** | Documented mitochondrial dysfunction in CIN models. |
| **Selective clearance of aneuploid cells** | Aneuploidy-selective lethality, senolytics. Weak evidence — flag as speculative or omit. |
| **Cancer chemoprevention & surveillance** | Cancer predisposition is a leading cause of morbidity. **Probably the highest-value output of the whole track, and almost nobody will submit it.** |
| **Symptomatic** | Seizures, growth, feeding. Low novelty, high immediate patient value. |

### 7.2 Target nomination — must be directional

For each target state **inhibit or activate**, citing the signed edge (SIGNOR / OmniPath /
Reactome). Unsigned KG proximity does not qualify as a direction.

Gate every target on: Open Targets tractability bucket · Pharos development level
(Tclin/Tchem preferred) · expression in affected tissue (GTEx, HPA — brain, kidney, muscle) ·
BBB penetrance if CNS · an approved or clinical-stage ligand (ChEMBL, DrugBank, Broad
Repurposing Hub ~6,000 compounds).

### 7.3 Signature reversal — *branch A only*

If patient transcriptomics exist, LINCS/CMap L1000 connectivity is the highest-yield
repurposing method available and needs no pocket, no structure, no known mechanism. Derive
patient-vs-control signature, query for reversal.

If no RNA-seq: use published transcriptomic signatures from MVA/CIN models (BUB1B-hypomorph
mice, aneuploid lines) as a proxy and **label the proxy clearly**.

### 7.4 Safety screen — load-bearing, not a footnote

For a child with a cancer predisposition syndrome and DNA-repair-adjacent biology, the
contraindication analysis carries more scientific content than the efficacy argument.

- ❌ **Genotoxic / mutagenic agents** — categorically excluded.
- ⚠️ **Immunosuppressants incl. rapalogs** — real tension with tumour immune surveillance.
  Do not silently propose everolimus. **Naming the tension is a strength**; ignoring it is
  the error a clinically-trained judge will spot first.
- ⚠️ **Radiosensitisers**; anything with an in vitro chromosomal-instability signal.
- ✅ Require paediatric exposure data, BBB penetrance if CNS-targeted, known paediatric PK,
  an existing paediatric label or trial.

Encode as **deterministic rules over FDA label sections and DrugBank fields** — not LLM judgment.

### 7.5 Structural work — optional, tightly scoped

Only if T1 yields a specific missense in a protein with an **experimental structure**, and
only as mechanistic illustration. Stability: ThermoMPNN + FoldX + RaSP consensus, accuracy
stated. Interface effects (BubR1–Bub3, MCC) are a stronger argument than ΔΔG alone. Docking
for enrichment only, experimental structure only, Vina/Uni-Dock → MM/GBSA rescoring.
**If no such variant emerges, skip entirely.**

---

## 8. LLM usage — narrow and instrumented

**For:** extraction (paper/label → structured claims + PMID + section anchor) · normalisation
(free text → HPO/MONDO/ChEBI IDs) · literature triage · drafting prose from a
deterministically-built evidence table.

**Not for:** deciding causality · deciding inhibit-vs-activate · generating any number
presented as a result · acting as its own critic (a skeptic agent on the same base model
shares its failure modes; the evidence that adversarial debate improves factual accuracy is thin).

**Architecture:** every call returns a pydantic schema where each field carries a `source`
with a resolvable identifier. Post-validate against the local DB or a live API; drop and log
unresolvable claims. **Track and report your hallucinated-identifier rate** — including it is
a credibility marker.

Local Qwen 3 32B AWQ (or a biomedical fine-tune) on the RTX 6000 for bulk PubMed extraction;
Claude for orchestration, code, and synthesis — subject to §1.

**Hard prohibition:** no LLM-generated dosing, dose ranges, safety margins, or
"regulatory-ready" summaries. A polished dossier confers authority the evidence doesn't
support. Output ranked hypotheses with evidence and gaps, addressed to researchers.

---

## 9. Phase 5 — Submission (Weeks 7–8). Code freeze 17 October.

1. **T1 report** — ranked candidates, ACMG classification with criteria codes, evidence,
   explicit uncertainty, and the confirmatory experiment per candidate (Sanger, RT-PCR for
   splicing, western, karyotype correlation).
2. **T2 report** — directional target rationale, ranked candidates, safety contraindications,
   and the specific *in vitro* assay per candidate. Patient fibroblasts are likely available
   through MVA Society — propose concrete readouts: aneuploidy fraction by karyotype/FISH,
   micronucleus assay, SAC integrity under nocodazole challenge, proliferation.
3. **Benchmark results** — recall by variant class, ClinVar-masked and unmasked.
4. **Negative results** — what you excluded, and how.
5. **Reproducibility artefact** — Docker image, Snakemake DAG, pinned versions,
   `make reproduce`. **Test on a clean machine.**
6. `DATA_CARD.md`, `PROVENANCE.md`, `ETHICS.md`, `RULES.md`.

### 9.1 🛑 STOP #4 (17 October)

Freeze. Last week is writing and reproduction verification only. Submit ≥48 h early —
submission endpoints get congested at deadlines.

---

## 10. Timeline

| Week | Dates | Focus | Gate |
|---|---|---|---|
| 0 | Aug 30 – Sep 5 | Phase 0 recon, governance, scaffold, container | 🛑 STOP #1 — branch confirmed |
| 1 | Sep 6 – 12 | Benchmark curation, positive-control harness | 🛑 STOP #2 — recall table |
| 2 | Sep 13 – 19 | Arms A + B | Baseline reproduces clinical negative |
| 3 | Sep 20 – 26 | Arms C + D (+ RNA if branch A) | recall@20 by class measured |
| 4 | Sep 27 – Oct 3 | Arm E, evidence integration | 🛑 STOP #3 — T1 list v1 |
| 5 | Oct 4 – 10 | T2 nomination, evidence, safety screen | T2 candidate list v1 |
| 6 | Oct 11 – 17 | Refinement, structural work if warranted, repro test | 🛑 STOP #4 — code freeze |
| 7 | Oct 18 – 24 | Writing, Docker verification, submission | Submitted ≥48 h early |

**Buffer discipline:** if Arm B is producing signal by end of week 3, cut Arms C/E scope and
go deep. One well-evidenced cryptic allele beats six shallow arms.

---

## 11. Claude Code kickoff prompts

**P0 — recon (run first, alone)**
> Read `MVA_HACKATHON_PLAN.md` §1–2. Create `.gitignore` and the pre-commit hook from §1
> BEFORE `git init`. Write `scripts/00_inventory.sh` and `scripts/01_characterise.py` per
> §2.1–2.2 and run them against `./data`. Populate `DATA_CARD.md` and `config/config.yaml`
> (genome build, sample map, branch). Write `RECON.md` stating which of the six hypothesis
> classes in §0.2 are testable and why the others are not. Then HALT at STOP #1 and report
> the branch. Do not read any file under `data/` into your context — the scripts read it,
> you read their outputs in `results/recon/`.

**P1 — scaffold + container**
> Given the confirmed branch in `RECON.md`, create the §3 scaffold. Write `CLAUDE.md` from
> §3.1, filling in the genome build from `config/config.yaml`. Set up `pyproject.toml` (uv,
> fully pinned), `environment.yml` with the bioconda tools referenced in §6 for the active
> arms only, and a Dockerfile. Verify GPU visibility inside the container. Record all tool
> and DB versions in `config/db_versions.yaml` and `PROVENANCE.md`.

**P2 — benchmark + harness**
> Build `benchmarks/published_mva_variants.tsv` per §5.1. Search the literature for causal
> variants in BUB1B, CEP57, TRIP13, BUB1, BUB3, CEP192, SMC5 in MVA/Atelis patients. Record
> gene, HGVS c. and p., zygosity, class, PMID, and whether the source paper noted that a
> standard pipeline missed it. Cite a PMID per row; leave `TODO(source)` rather than
> guessing. Implement `tests/test_positive_controls.py` per §5.2 with the ClinVar-masked
> variant of §5.3. Run it and HALT at STOP #2 with the recall table.

**P3 — gene panel**
> Build `config/gene_panels/mitotic_extended.tsv` per §6.5. Per gene: symbol, Ensembl ID,
> nominating source (which GO term / Reactome pathway / CORUM complex / STRING seed),
> gnomAD LOEUF, pLI. Write `panel_provenance.md` documenting assembly so it is reproducible.

**P4 — splicing arm**
> Implement `src/mva/track1/splicing.py` per §6.2. SpliceAI and Pangolin over the extended
> panel's introns, UTRs and 2 kb promoters at ±500 bp. Add branch point scoring and uORF
> annotation. Every output row carries tool version and model weights hash. Unit tests using
> the splice-class positive controls.

---

## 12. Report back after STOP #1

- Confirmed branch (A–E) and the file manifest
- Genome build
- Whether RNA-seq is present
- Whether the clinical data already names a candidate het in a known MVA gene
- Reported aneuploidy percentages and tissues

Any of these change the plan materially; the last two change it structurally.
