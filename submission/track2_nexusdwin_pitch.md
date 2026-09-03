# Track 2 pitch: 3-minute video script

**Team:** NexusDwin
**Spoken length:** 418 words, 179 seconds at 140 words per minute,
leaving room for pauses at the six section breaks. Timings are cumulative.
Verify the count with `make pitch-wordcount` rather than trusting this line.

The argument is built around a negative result on purpose. Most entries will
show a list of drug names. What we have that they will not is evidence about
which of those lists cannot work, and why.

---

## 0:00 to 0:25 | The problem, and the obvious idea

> A child has two broken copies of `BUB1B`. That gene makes the brake that stops
> a cell dividing before its chromosomes line up. With it weakened, chromosomes
> land in the wrong daughter cells, differently in every cell. That is mosaic
> variegated aneuploidy, and it is why this child has already had a
> rhabdomyosarcoma.
>
> So the obvious idea is to fix the brake.

**On screen:** the two variants, `p.Leu737Ter` and `p.Asn1002Lys`. Then
chromosomes missegregating.

---

## 0:25 to 1:05 | The obvious idea fails

> We tested it, and it fails.
>
> We nominated targets using signed causal edges, which record whether one
> protein activates or inhibits another. Network proximity tells you which
> protein to look at, never what to do to it.
>
> Ten targets came back, and every one needs *activation*. ChEMBL holds one
> hundred and eighteen drug mechanisms across them. None activates. They are all
> inhibitors, and six of the ten have one ready to go, pointing exactly the wrong
> way for a loss-of-function disease.

**On screen:** the ten target names, then 118 mechanisms, 0 activating, 6 of 10
inhibitable.

---

## 1:05 to 1:40 | Then we attacked our own result

> That is striking, so we tried to break it.
>
> If activating drugs are rare everywhere, finding none across ten targets says
> nothing about this disease. So we measured it. In all of ChEMBL, three hundred
> and fifty-nine human genes have an activating drug and one thousand five
> hundred and forty-one have an inhibiting one. Activation reaches under two per
> cent of the genome.
>
> At that rate you would expect 0.19 activatable targets out of ten. The chance
> of seeing zero is 0.83.
>
> Our headline finding is weaker than it looked. We are telling you, rather than
> waiting for a reviewer to find it.

**On screen:** the base-rate table, then 0.83 held on screen.

---

## 1:40 to 2:05 | What survives

> Three things survive, and they are enough.
>
> There is nothing to repurpose today. Every target needs the scarce direction
> and none needs the plentiful one. And activating mitotic kinases in a child
> already predisposed to cancer would be contraindicated even if someone built
> the drug tomorrow.
>
> That axis is closed.

**On screen:** the three surviving claims, one line each.

---

## 2:05 to 2:35 | The second negative

> So we went to the cancer risk, because the tumour is what is actually hurting
> this child.
>
> Zero registered trials in mosaic variegated aneuploidy. Zero naming `BUB1B`.
> Zero prevention trials in rhabdomyosarcoma. Not a thin evidence base. An empty
> one.
>
> Everything we can offer is transferred from a different syndrome. We derived it
> from the trial registry's own vocabulary, not from what a model remembers, and
> screened every candidate with code rather than judgement.

**On screen:** four zeros. Then the candidate table with the flags visible.

---

## 2:35 to 2:50 | The close

> The honest answer is not a drug name.
>
> For this child, a programme aimed at their tumour risk is better spent on
> surveillance than on any candidate we can name. Everything reproduces with one
> command, and the number that weakens our own best finding has a unit test.

**On screen:** `make track2`, then the passing line from `make test` as it actually prints on the day, then the team name.

---

## Production notes

- **Do not show patient data.** No VCF, no BAM, no genotypes, no read pileups.
  The two variant coordinates are already public in the Track 1 submission and
  are the only patient-derived items permitted on screen.
- **Do not name or depict the family**, and use no image that could identify the
  child. The consent scope is in `ETHICS.md` section 1.
- **No dose may appear on screen**, including inside a screenshot of a trial
  record. Check every frame.
- Show real output from `results/summaries/`, not mocked tables. Every figure
  above is checked against the generated data by
  `tests/test_track2_report_matches_data.py`, so numbers on screen should be
  copied from those files.
- Keep the delivery plain and unhurried. The argument is that we checked our own
  work, and rushing it undercuts that.
