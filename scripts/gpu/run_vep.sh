#!/usr/bin/env bash
# MVA Arm A: consequence annotation. CPU-only, so it does not contend with the
# GPU. Scoped run first (disease-panel CDS, ~10k variants, minutes), then the
# whole callset in the background.
#
# The conda perl must precede /opt/ucl/bin/perl or VEP cannot find BioPerl.
set -uo pipefail
S=/REDACTED; B=$S/mva/mamba/envs/mva
export PATH=$B/bin:/usr/bin:/bin
unset PERL5LIB
export TMPDIR=$S/mva/tmp
ST=$S/mva/logs/STATUS_vep
say(){ echo "$(date -u +%H:%M) $*" | tee -a "$ST"; }

VCF=$S/mva/data/WGS_EX2312012_HGWCNDSX7.vcf.gz
CACHE=$S/mva/refs/vep
FASTA=$S/mva/refs/GRCh38.fa.gz
OUT=$S/mva/results
mkdir -p $OUT

# Wait for the cache and reference, both fetched on the node's own network.
for i in $(seq 1 240); do
  [ -d $CACHE/homo_sapiens_merged ] && [ -s $FASTA.fai ] && break
  sleep 30
done
[ -d $CACHE/homo_sapiens_merged ] || { say "FATAL cache absent"; exit 1; }
say "cache and reference ready"

common=( --cache --merged --offline --dir_cache $CACHE --fasta $FASTA
         --assembly GRCh38 --species homo_sapiens
         --everything --pick_allele_gene --check_existing
         --fork 12 --vcf --compress_output bgzip --no_stats --force_overwrite )

# 1. Disease-panel CDS. Fast, and it is where an achievable answer sits.
if [ ! -s $OUT/panel.vep.vcf.gz ]; then
  say "VEP over the disease-panel CDS"
  bcftools view -R $S/mva/refs/panels/disease_cds.nochr.bed -Oz \
    -o $S/mva/tmp/panel.vcf.gz "$VCF" && bcftools index -t -f $S/mva/tmp/panel.vcf.gz
  say "panel variants: $(bcftools index -n $S/mva/tmp/panel.vcf.gz)"
  vep "${common[@]}" -i $S/mva/tmp/panel.vcf.gz -o $OUT/panel.vep.vcf.gz \
    && say "panel VEP done" || say "panel VEP FAILED"
fi

# 2. Whole callset. Slower; the panel result is already usable by then.
if [ ! -s $OUT/all.vep.vcf.gz ]; then
  say "VEP over the whole callset (5.0M records)"
  vep "${common[@]}" -i "$VCF" -o $OUT/all.vep.vcf.gz \
    && say "full VEP done" || say "full VEP FAILED"
fi
say "ALL DONE"
