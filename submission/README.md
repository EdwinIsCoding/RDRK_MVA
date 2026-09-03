# Submission artefacts

These are deliverables, not intermediate results, so they are tracked whereas
everything under `results/` is not.

Committing variant coordinates here is deliberate and consistent with the rules.
The organisers' data-handling guidance lists "candidate variant rankings" among
the things participants may keep, as distinct from VCF and BAM files, variant
tables with genotypes, and prompts containing pasted variant blocks, which must
be deleted. The submission itself is published under CC BY 4.0 and the
organisers state that all submitted results are made accessible to the research
community.

| File | Purpose |
|---|---|
| `track1_submission.csv` | The Track 1 file to upload. Format validated by `scripts/23_track1_submission.py`. |
| `track1_candidates.tsv` | Human-readable source with the reasoning per row. |
| `track1_nexusdwin_report.md` | The Track 1 methods write-up, judged separately. |
| `arm_c_readlevel_verification.md` | Read-level confirmation from our own alignment. |
| `arm_d_mosaic.md` | The mosaic arm, reported negative. |
| `track2_nexusdwin_report.md` | **The Track 2 report.** Mechanism characterisation, the closed direct axis, the chemoprevention axis and its safety screen. |
| `track2_nexusdwin_pitch.md` | **The 3-minute pitch script**, with production constraints. The video itself is not yet recorded. |

## Track 2 submission checklist

The rules require three artefacts. Two exist.

| Required | Status |
|---|---|
| Written report, PDF or Markdown, filename carrying the team name | **Done**, `track2_nexusdwin_report.md` |
| Public GitHub repository with reproducible code | Repository exists and is public. History was rewritten and force-pushed on 3 September 2026 to purge third-party cluster hostnames and an account name. Verified by cloning the public repository and scanning all history: zero remain. See `docs/VERIFICATION.md` section 6. |
| 3-minute pitch video on YouTube or Vimeo | **Not done.** The script is written and timed at 179 seconds; the recording and upload remain. |

Only 3 Track 2 submissions are permitted per team and only the latest is
reviewed, so there is no value in submitting early to probe a score.

## Scientific limitations: seven declared, four remain

Four of the seven were addressed on 3 September 2026. The three that remain are
blocked by data that does not exist, not by effort.

| Limitation | Status |
|---|---|
| Structural variant calling never completed | **Partly closed.** A calibrated breakpoint screen over the nine MVA genes is negative. Genome-wide calling is still undone. |
| Hypomorph model is an inference from viability | **Strengthened, still an inference.** Five curated MVA1 kinase-domain missense variants are all compound heterozygous with a nonsense allele, the nearest 10 residues from ours (PMID 15475955). No functional assay exists and no experimental structure resolves residue 1002. |
| Mosaic arm used a diploid model, not a somatic caller | **Closed.** GATK Mutect2 tumour-only agrees with `bcftools mpileup`: zero credible mosaicism. Both known alleles recovered at PASS as a positive control. |
| Chemoprevention rests on transferring evidence across syndromes | **Confirmed by search, not removed.** PubMed returns zero for `Bub1b AND chemoprevention`. The assumption stands and is stated. |
| Phase is inferred | **Cannot be closed.** Zero read pairs span the 10,911 bp, and statistical phasing carries no information for a private allele. |
| Every splicing result is a prediction | **Cannot be closed.** No RNA-seq exists. |
| Repeat expansion detection never run | **Declined.** The panel BAM covers *HTT* at 73 reads and *C9orf72* at none, so a result would look genome-wide and would not be. |

## Everything still outstanding, as of 3 September 2026

The authoritative list. Anything not here is done.

| # | Item | Who |
|---|---|---|
| 1 | **Upload the Track 1 predictions file and methods report.** Nothing has been submitted to either track. 6 Track 1 submissions are allowed per participant and the highest score counts. | owner |
| 2 | **Record and upload the 3-minute pitch video.** The script is written and timed at 178 seconds. Track 2 is incomplete without it. | owner |
| 3 | **Upload the Track 2 report and repository link.** 3 submissions per team, only the latest reviewed. | owner |
| 4 | **Notify Sage Bionetworks' Privacy and Compliance Office** of the training-setting deviation in `ETHICS.md` section 3a. Both reports state this is outstanding. | owner |
| 5 | Optionally ask GitHub Support to purge the now-unreachable pre-rewrite objects, which they will serve by SHA until asked. Not urgent: the redacted strings are third-party infrastructure names, not credentials. | owner |
| 6 | Delete the local `pre-sanitise-backup` tag and the pre-rewrite bundle, both of which still hold the unredacted history. | owner |
| 7 | At the conclusion of the hackathon, run `scripts/33_delete_challenge_data.py --execute` and notify `MVAHackathon2026@synapse.org`. | owner |

**The three provenance gaps are now closed**, 3 September 2026. The HuggingFace
dataset revision was recovered as `59e322d2` and is labelled an inference rather
than a log entry, because it rests on the repository's last-modified date
preceding our download. The Ensembl REST release number turned out to be moot:
those coordinates were replaced by the pinned Ensembl 115 GTF after the live API
overstated the `BUB3` span 9.9-fold, and `panel_provenance.md` had still been
naming REST as the source. Host tool versions are recorded, with the caveat that
they were read afterwards rather than at the time.

The three `TODO(source)` markers that remain in `config/db_versions.yaml` are
genuine: the SpliceAI model-weights hash, the OmniPath release and the GO release
behind QuickGO are not published through the interfaces used.

## Verification performed before submitting

- Coordinates read from the **original callset**, not from the VEP output.
- REF alleles checked against the GRCh38 reference: `15:40209701` T, `15:40220612` T,
  `22:20996720` C all match.
- Chromosomes chr-prefixed, since the submission requires it and the callset
  does not use it.
- The compound-heterozygous pair is in a **single row** using the `_2` columns,
  which is what earns the partial credit the metric offers for recovering one of
  two variants.
- Format gate passed. It had already rejected an earlier hand-written file in
  which one missing tab put `epcr` in the `alt_2` column.
