#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# MVA Hackathon 2026 - Phase 0 recon, step 1.
#
# Reads the patient data directory and writes structural descriptions of it into
# results/recon/. Nothing written here contains variant-level genotypes: headers,
# sample identifiers, contig tables and aggregate statistics only.
#
# Portable across GNU and BSD userland (macOS ships BSD find, which has no
# -printf, and shasum rather than sha256sum).
#
# Usage: scripts/00_inventory.sh [DATA_DIR] [OUT_DIR]
# ---------------------------------------------------------------------------
set -euo pipefail

DATA=${1:-./data}
OUT=${2:-results/recon}
mkdir -p "$OUT"

if [ ! -d "$DATA" ]; then
  echo "FATAL: data directory not found: $DATA" >&2
  exit 1
fi

# Portable checksum helper.
if command -v sha256sum >/dev/null 2>&1; then
  SHA() { sha256sum "$1"; }
elif command -v shasum >/dev/null 2>&1; then
  SHA() { shasum -a 256 "$1"; }
else
  echo "FATAL: no sha256sum or shasum on PATH" >&2
  exit 1
fi

# Portable byte-size helper.
if stat -f%z . >/dev/null 2>&1; then
  SIZE() { stat -f%z "$1"; }          # BSD / macOS
else
  SIZE() { stat -c%s "$1"; }          # GNU
fi

echo "== 1. File manifest =="
printf 'path\tbytes\n' > "$OUT/manifest.tsv"
find "$DATA" -type f ! -name '.DS_Store' | sort | while IFS= read -r f; do
  printf '%s\t%s\n' "$f" "$(SIZE "$f")"
done >> "$OUT/manifest.tsv"

echo "== 2. SHA256 checksums (slow on large FASTQ; skip with SKIP_SHA=1) =="
if [ "${SKIP_SHA:-0}" = "1" ]; then
  echo "SKIPPED (SKIP_SHA=1)" > "$OUT/sha256.txt"
else
  : > "$OUT/sha256.txt"
  find "$DATA" -type f ! -name '.DS_Store' | sort | while IFS= read -r f; do
    echo "  hashing $(basename "$f")" >&2
    SHA "$f" >> "$OUT/sha256.txt"
  done
fi

echo "== 3. Real file types (do not trust extensions) =="
printf 'path\tfile_type\n' > "$OUT/filetypes.tsv"
tail -n +2 "$OUT/manifest.tsv" | cut -f1 | while IFS= read -r f; do
  printf '%s\t%s\n' "$f" "$(file -b "$f")"
done >> "$OUT/filetypes.tsv"

echo "== 4. VCF characterisation =="
shopt -s nullglob 2>/dev/null || true
for v in "$DATA"/*.vcf.gz "$DATA"/*.vcf "$DATA"/*.bcf; do
  [ -e "$v" ] || continue
  b=$(basename "$v")
  echo "  $b"
  bcftools view -h  "$v" > "$OUT/${b}.header.txt"       2>"$OUT/${b}.header.err"  || true
  bcftools query -l "$v" > "$OUT/${b}.samples.txt"      2>/dev/null              || true
  # Index if absent; needed for contig statistics.
  if [ ! -e "${v}.tbi" ] && [ ! -e "${v}.csi" ]; then
    bcftools index -t "$v" 2>/dev/null || bcftools index "$v" 2>/dev/null || true
  fi
  bcftools index -s "$v" > "$OUT/${b}.contigs.tsv"      2>/dev/null              || true
  bcftools stats   "$v" > "$OUT/${b}.stats.txt"         2>"$OUT/${b}.stats.err"  || true
done

echo "== 5. Alignment characterisation =="
for a in "$DATA"/*.bam "$DATA"/*.cram; do
  [ -e "$a" ] || continue
  b=$(basename "$a")
  echo "  $b"
  samtools view -H  "$a" > "$OUT/${b}.header.txt"   2>/dev/null || true
  samtools idxstats "$a" > "$OUT/${b}.idxstats.tsv" 2>/dev/null || true
done

echo "== 6. FASTQ characterisation (headers only, no sequence) =="
for q in "$DATA"/*.fastq.gz "$DATA"/*.fq.gz "$DATA"/*.fastq "$DATA"/*.fq; do
  [ -e "$q" ] || continue
  b=$(basename "$q")
  echo "  $b"
  # First four read identifier lines only. Read names carry instrument, run,
  # flowcell, lane and index, which is what routes lane merging and read groups.
  # No sequence, no quality strings are written.
  {
    echo "# first 4 read identifier lines (@ lines only)"
    gzip -cd "$q" 2>/dev/null | head -16 | awk 'NR % 4 == 1'
    echo "# read length of first record"
    gzip -cd "$q" 2>/dev/null | head -2 | awk 'NR == 2 { print length($0) }'
  } > "$OUT/${b}.readinfo.txt" || true
done

echo "== 7. Tabular and document structure (no contents) =="
for t in "$DATA"/*.csv "$DATA"/*.tsv "$DATA"/*.txt "$DATA"/*.json; do
  [ -e "$t" ] || continue
  head -c 2000 "$t" > "$OUT/$(basename "$t").head.txt" 2>/dev/null || true
done
for d in "$DATA"/*.docx "$DATA"/*.xlsx "$DATA"/*.pdf; do
  [ -e "$d" ] || continue
  b=$(basename "$d")
  # Structure only at this stage. Clinical text extraction is a separate,
  # explicitly governed step (scripts/02_extract_clinical.sh), because the
  # phenotype document is patient data even though it is not genomic.
  printf 'file\t%s\ntype\t%s\nbytes\t%s\n' "$b" "$(file -b "$d")" "$(SIZE "$d")" \
    > "$OUT/${b}.structure.txt"
done

echo
echo "Recon written to $OUT"
