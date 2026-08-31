#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# MVA Hackathon 2026 - Phase 0 recon, step 3: clinical phenotype extraction.
#
# GOVERNANCE NOTE. This is a separate, deliberately named script because the
# clinical phenotype document is patient data, even though it is not genomic.
# Plan section 2.2 and section 12 both require the phenotype content to route
# the analysis (HPO availability, reported aneuploidy percentages and tissues,
# and whether a candidate heterozygous variant has already been named). This
# script converts the document to plain text under results/recon/ so that the
# extraction is auditable and version-pinned rather than ad hoc.
#
# Output stays gitignored. Only the derived HPO term list and a de-identified
# aggregate summary are promoted to results/summaries/ for the report.
# ---------------------------------------------------------------------------
set -euo pipefail
DATA=${1:-./data}
OUT=${2:-results/recon}
mkdir -p "$OUT"

command -v pandoc >/dev/null || { echo "FATAL: pandoc not on PATH" >&2; exit 1; }

for d in "$DATA"/*.docx "$DATA"/*.pdf; do
  [ -e "$d" ] || continue
  b=$(basename "$d")
  echo "extracting $b"
  pandoc -f docx -t markdown --wrap=none "$d" > "$OUT/${b}.md" 2>/dev/null \
    || pandoc -t markdown --wrap=none "$d" > "$OUT/${b}.md"
  wc -w "$OUT/${b}.md"
done
