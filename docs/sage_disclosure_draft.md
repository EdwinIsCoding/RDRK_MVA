# Draft: disclosure to Sage Bionetworks' Privacy and Compliance Office

**Status: draft for the owner to review, edit and send. Not sent.**

This concerns the deviation recorded in `ETHICS.md` section 3a and disclosed in
section 8 of the Track 1 report and section 9 of the Track 2 report. The facts
below are taken from those documents; the judgement about how to characterise
them is the owner's, not this repository's.

**Before sending, check three things.** That the account and dates are right, since
this repository records them second-hand. That you are content with the framing,
because it concedes non-compliance rather than arguing the edges. And whether
your registration used a different email from the one you will send from, since
they will need to link the report to a participant.

**Where to send it.** The Hackathon Rules direct suspected unauthorised
disclosure to Sage Bionetworks' Privacy and Compliance Office, reachable through
Sage's Help Center. This is a separate obligation from the deletion notification
to `MVAHackathon2026@synapse.org`, which falls due at the conclusion of the
hackathon and is handled by `scripts/33_delete_challenge_data.py`.

---

**Subject:** Self-reported LLM data-handling deviation, Rare Disease Real Kid MVA Hackathon 2026

Dear Privacy and Compliance Office,

I am a participant in the Rare Disease, Real Kid: MVA Hackathon 2026, and I am
writing to self-report a deviation from the LLM data-handling conditions in the
guidance published on 28 August 2026.

**What happened.** I used Anthropic's Claude, through Claude Code, on a Claude
Max subscription, which is governed by consumer terms rather than the commercial
terms in your example disclosure. The "help improve Claude" setting, which
permits training on conversation content, was enabled from the start of my work
until 31 August 2026, when I discovered the requirement and disabled it. For
that period the arrangement did not meet your first mandatory condition, that
there be no training on inputs or outputs and no rights taken by the provider in
either.

**What was and was not exposed.** The project was built from the outset so that
patient genomic data never entered a model context. Scripts read the challenge
data and wrote aggregate summaries; the summaries are what the assistant saw.
That rule was written before any analysis, is enforced by a pre-commit hook, and
held throughout. Specifically:

- No VCF or BAM content, no variant table with genotypes, and no prompt
  containing a pasted variant block ever entered a model context. These are the
  categories your guidance places on the delete list.
- What did enter falls on your keep list: the eight coded HPO terms, gene-level
  and aggregate statistics, and code.

Two patient-derived items reached the model and I name them rather than gloss
them. The clinical phenotype document was read in full, including the gestational
age, approximate birth weight and family history; this was a deliberate decision,
recorded with its reasoning. And five runs-of-homozygosity interval coordinates
were printed while summarising an analysis. These are genomic intervals rather
than variants, carrying no genotypes or alleles.

The pseudonymous sample identifier assigned by the organisers also appears
throughout.

**What I have done.** The setting is disabled. Hackathon-related content has been
deleted from the assistant's stored memory and conversation history, and an
account data export taken on 31 August 2026 contained no hackathon-related
material. I should be precise that this speaks to what the account retains and is
not evidence about training pipelines, which are separate infrastructure and not
something I can inspect or control.

**What I am asking.** Whether you consider any further action necessary on my
part, and whether you wish this recorded against my participation. I have also
raised with Anthropic whether already-collected conversation content can be
removed from training pipelines under the consumer terms; that is a question only
they can answer.

**Why I am reporting it.** The Hackathon Rules require suspected unauthorised
disclosure to be reported, and the deviation is disclosed in my submitted methods
descriptions for both tracks. I would rather you heard it from me, and I would
rather a judge read an accurate account than discover an inconsistency.

I am happy to provide any further detail that would help.

Yours sincerely,

[name]
[registration email, if different from the sending address]
[team name: NexusDwin]
