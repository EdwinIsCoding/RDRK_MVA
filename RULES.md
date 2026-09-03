# RULES.md

Transcribed from the official Hackathon Space on 31 August 2026:
https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space/

Source text preserved at `results/rules/page_text.md` (extracted from the Gradio
config endpoint, since the page is client-rendered).

---

## 1. Dates

| Event | Date |
|---|---|
| Submissions open | 25 August 2026 |
| **Submissions close** | **24 October 2026, 23:59 UTC** |
| Judging window | 24 October to 24 November 2026 |
| Winners announced | 25 November 2026 |

## 2. Judging criteria and weights

Applied by the independent expert panel:

| Criterion | Weight |
|---|---:|
| **Scientific rigor** | **35%** |
| Innovation | 25% |
| Potential impact | 25% |
| Scalability | 15% |

## 3. Prizes

$25,000 cash (AWS Imagine Grant) plus $25,000 in Anthropic Claude credits,
split first $12,000, second $7,000, third $4,000, and a $2,000
Innovation/Community award at the judges' discretion, each with matching credits.

## 4. Track 1: variant identification

**There is a confirmed clinical answer.** The leaderboard scores submissions
"against the confirmed clinical groundtruth", automatically and instantly. The
organisers describe Track 1 as "a **foundational track** - it's designed to be
achievable."

- **6 submissions per participant.** Each team member has their own quota; only
  the highest-scoring counts, and identical team names group entries.
- **Up to 10 candidate rows.** "This is one case, not a cohort, so we're asking
  for your best-ranked guesses, not an exhaustive list."
- **Metrics:** rank points (how high the true variant lands, with partial credit
  for recovering one of two compound-heterozygous variants) and F-max. Adapted
  from Stenton et al. 2024.
- **A methods write-up is also required and is judged separately.** "Recovering
  the correct variant(s) is not the finish line for Track 1."

### Submission format

One row per proposed causal variant or compound-heterozygous pair. **GRCh38.**

| Field | Notes |
|---|---|
| `proband_id` | from the dataset |
| `chrom_1`, `pos_1`, `ref_1`, `alt_1` | first variant. **Chromosome format is `chr15`, chr-prefixed** |
| `chrom_2`, `pos_2`, `ref_2`, `alt_2` | second variant, compound-het pairs only, blank otherwise |
| `epcr` | estimated probability of causal relationship, (0, 1]. Only the ranking relative to our own guesses matters |
| `finding_type` | `primary` or `secondary` |
| `notes` | optional rationale, especially for secondary findings |

**Secondary and incidental findings are welcomed** and do not hurt the automated
score. They are set aside for qualitative panel review.

## 5. Track 2: drug repurposing

Qualitative review by an independent panel. **3 submissions per team**, one
designated submitter, only the latest reviewed.

Required:
1. **Written report**, PDF or Markdown, proposing repositioned drug candidates.
   Must include "a characterization of the variant's mechanism (loss-of-function
   / gain-of-function, pathway disrupted, downstream biological consequence) as
   the basis for your repurposing rationale". Filename should carry the team name.
2. **Public GitHub repository.** May stay private during the hackathon, "must be
   made public once the Hackathon ends so the panel can review your code and
   methods". Must contain documented, reproducible code: all scripts and
   configuration needed to reproduce the results.
3. **3-minute pitch video** on YouTube or Vimeo.

## 6. Constraints on how we work

These are obligations, not preferences.

- **Redistribution is prohibited.** The dataset is gated under the Hackathon
  Rules and a Data Transfer Agreement.
- **All data must be deleted at the conclusion of the Hackathon.** This needs
  scheduling, not remembering. See `ETHICS.md`.
- **No recontact.** "You will not attempt to recontact the data subject, data
  subject family members, or any points of contact at the MVA Society." This
  independently confirms the Phase 0 decision not to look for the family's blog.
- **No onward access.** "You will not release or otherwise grant data access to
  anyone, and you will establish appropriate safeguards to prevent unauthorized
  data use."
- **Suspected disclosure** must be reported to Sage Bionetworks' Privacy and
  Compliance Office.
- **Outputs are CC BY 4.0**: predictions, code and reports.
- Participants must be 18 or over; each team member registers individually.

### LLM disclosure is mandatory

Added by the organisers on 28 August 2026:

> "If you used an LLM or AI assistant, please record the provider, the plan or
> tier, and the relevant data-handling setting in your methods description. A
> line is enough, for example: *'Anthropic API, Claude Sonnet, commercial terms,
> no training on customer content.'*"

**Resolved 3 September 2026.** The provider, plan and data-handling setting are
established and the disclosure line is written verbatim into `ETHICS.md`
section 3a, section 8 of the Track 1 report and section 9 of the Track 2 report.
It discloses a compliance gap rather than a clean summary: the account is a
Claude Max subscription on consumer terms and the training setting was enabled
until 31 August 2026.

**Still outstanding, and it is the owner's to do, not this repository's:**
notification to Sage Bionetworks' Privacy and Compliance Office. It had not been
made at the time of writing and the reports say so.

## 7. Compute

The Space itself runs on CPU-basic. "You are welcome to use your own compute for
training, as only the final submission files need to be uploaded here." There is
no compute restriction on participants.
