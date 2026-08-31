#!/usr/bin/env bash
# VEP cache plus the reference FASTA, downloaded on the node's own network
# rather than pushed from the laptop, which measures ~1 MB/s.
set -uo pipefail
S=/REDACTED; B=$S/mva/mamba/envs/mva
export PATH=$B/bin:/usr/bin:/bin
unset PERL5LIB
ST=$S/mva/logs/STATUS_vepcache
say(){ echo "$(date -u +%H:%M) $*" | tee -a "$ST"; }
mkdir -p $S/mva/refs/vep

say "start"
# Ensembl 116 merged cache: RefSeq plus Ensembl transcripts, GRCh38.
if [ ! -d $S/mva/refs/vep/homo_sapiens_merged ]; then
  say "downloading VEP cache (~26 GB)"
  cd $S/mva/refs/vep
  curl -sSL --retry 8 --retry-all-errors -C - \
    -O https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_merged_vep_116_GRCh38.tar.gz \
    && say "cache downloaded, extracting" \
    && tar xzf homo_sapiens_merged_vep_116_GRCh38.tar.gz \
    && rm -f homo_sapiens_merged_vep_116_GRCh38.tar.gz \
    && say "cache extracted"
fi
# Reference FASTA. Ensembl primary assembly uses no-chr naming, matching the callset.
if [ ! -s $S/mva/refs/GRCh38.fa.gz ]; then
  say "downloading reference FASTA"
  curl -sSL --retry 8 --retry-all-errors -C - -o $S/mva/refs/GRCh38.fa.gz.raw \
    https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    && gzip -dc $S/mva/refs/GRCh38.fa.gz.raw | bgzip -@ 8 -c > $S/mva/refs/GRCh38.fa.gz \
    && rm -f $S/mva/refs/GRCh38.fa.gz.raw \
    && samtools faidx $S/mva/refs/GRCh38.fa.gz \
    && say "reference ready"
fi
say "ALL DONE"
