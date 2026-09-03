#!/usr/bin/env python3
"""Execute the deletion obligation in ETHICS.md section 3b.

The Hackathon Rules require all challenge data to be deleted at the conclusion
of the hackathon, and ETHICS.md section 3b says this "needs scheduling rather
than remembering". A prose list is a thing to remember. This is the thing to
run.

Safe by default
---------------
Nothing is deleted without ``--execute``. The default lists what would go, what
would stay, and what is missing, so the plan can be read before it is enacted.
Deletion is irreversible and the inputs are 79 GB that cannot be re-obtained
once the challenge closes, so the default is the cautious one.

What is deleted, and what is kept
---------------------------------
Straight from ETHICS.md section 3b, which in turn follows the organisers' own
delete and keep lists. The keep list is not a loophole: it is the material the
organisers explicitly permit participants to retain, being free of raw genomic
data.

This script deletes only what is on this machine. Provider logs, the assistant
transcript and anything already uploaded are out of its reach and are listed as
manual items, because a script that silently omitted them would imply the
obligation was discharged when it was not.
"""
from __future__ import annotations

import argparse
import pathlib
import shutil
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

#: (path, why). Everything here is challenge data or derived from it.
DELETE = [
    ("data", "the challenge distribution: VCF, index, eight FASTQ, phenotype document"),
    ("results/recon", "extracted phenotype text and Phase 0 recon outputs"),
    ("results/arm_a_shortlist.tsv", "variant-level shortlist"),
    ("results/arm_a_shortlist_rhabdo.tsv", "variant-level shortlist"),
    ("results/arm_a_shortlist_x.tsv", "variant-level shortlist"),
    ("results/arm_b", "variant-level splicing input and output"),
    ("results/arm_b_x", "variant-level splicing input and output"),
    ("results/summaries/arm_f_roh.raw.txt", "runs-of-homozygosity intervals"),
    ("results/summaries/arm_f_mtdna.raw.tsv", "mitochondrial variant table"),
    ("results/disease_cds_variants.vcf.gz", "proband variants"),
    ("results/disease_cds_variants.vcf.gz.tbi", "proband variants"),
    ("results/rhabdo_variants.vcf.gz", "proband variants"),
    ("results/rhabdo_variants.vcf.gz.tbi", "proband variants"),
    ("node_artefacts", "panel BAM, VEP-annotated VCF, coverage output, node logs"),
]

#: Public reference data cached for this analysis. Not patient data, and not
#: required to be deleted, but large and regenerable.
OPTIONAL = [
    ("refs", "public reference data: FASTA, GTF, gnomAD slices, HPO, HGNC"),
]

#: Retained per the organisers' keep list. None contains raw genomic data.
KEEP = [
    ("submission", "candidate variant rankings and the reports, both permitted"),
    ("results/summaries", "aggregate summaries, minus the raw tables listed above"),
    ("resources", "the published availability table, derived from ChEMBL only"),
    ("config", "configuration, gene panels, database versions"),
    ("scripts", "code"),
    ("src", "code"),
    ("tests", "code"),
    ("benchmarks", "public ClinVar benchmark, no proband data"),
]

#: Not reachable from here. Listed so the obligation is not mistaken for done.
MANUAL = [
    "Notify MVAHackathon2026@synapse.org that the data has been deleted, per the "
    "data access terms.",
    "Delete the local assistant transcript at "
    "~/.claude/projects/-Volumes-ROS2-SSD-RDRK-MVA/*.jsonl. ETHICS.md notes it "
    "contains a now-rotated HuggingFace token and account export URLs.",
    "Delete any copy of the dataset on the compute node or in cloud storage. "
    "Scratch was wiped when the booking ended, but confirm rather than assume.",
    "Confirm no copy remains on backup media or in a system snapshot.",
]


def size_of(p: pathlib.Path) -> int:
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true",
                    help="actually delete. Without this, nothing is removed.")
    ap.add_argument("--include-refs", action="store_true",
                    help="also delete refs/, which is public reference data")
    args = ap.parse_args()

    targets = list(DELETE) + (list(OPTIONAL) if args.include_refs else [])

    print("ETHICS.md section 3b, deletion at the conclusion of the hackathon")
    print("=" * 68)
    print(f"Repository: {REPO}")
    mode = ("EXECUTE, deletions are irreversible" if args.execute
            else "DRY RUN, nothing will be deleted")
    print(f"Mode: {mode}\n")

    total, present = 0, []
    print("TO DELETE")
    for rel, why in targets:
        p = REPO / rel
        if not p.exists():
            print(f"  [absent] {rel:44} {why}")
            continue
        n = size_of(p)
        total += n
        present.append(p)
        print(f"  [{human(n):>9}] {rel:44} {why}")
    print(f"\n  total to delete: {human(total)} across {len(present)} paths\n")

    print("TO KEEP, per the organisers' keep list")
    for rel, why in KEEP:
        p = REPO / rel
        mark = "ok" if p.exists() else "absent"
        print(f"  [{mark:>9}] {rel:44} {why}")

    if not args.include_refs:
        print("\n  refs/ is public reference data and is NOT deleted by default. "
              "Pass --include-refs to reclaim the space.")

    print("\nNOT REACHABLE FROM THIS SCRIPT, and still required")
    for item in MANUAL:
        print(f"  - {item}")

    if not args.execute:
        print("\nDry run. Re-run with --execute to delete.")
        return 0

    if not present:
        print("\nNothing to delete.")
        return 0

    print(f"\nAbout to permanently delete {human(total)}.")
    reply = input("Type the word DELETE to proceed: ").strip()
    if reply != "DELETE":
        print("Aborted. Nothing was deleted.")
        return 1

    for p in present:
        try:
            if p.is_dir():
                shutil.rmtree(p)
            else:
                p.unlink()
            print(f"  deleted {p.relative_to(REPO)}")
        except OSError as exc:
            print(f"  FAILED  {p.relative_to(REPO)}: {exc}", file=sys.stderr)

    print("\nLocal deletion complete. The manual items above are NOT done and "
          "the obligation is not discharged until they are.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
