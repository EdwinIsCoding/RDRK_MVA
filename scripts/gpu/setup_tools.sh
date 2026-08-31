#!/usr/bin/env bash
# MVA hackathon: bioinformatics toolchain on scratch. Idempotent, resumable.
set -uo pipefail
S=/REDACTED
export XDG_CACHE_HOME=$S/.cache PIP_CACHE_DIR=$S/.cache/pip TMPDIR=$S/mva/tmp
export MAMBA_ROOT_PREFIX=$S/mva/mamba
MM=$S/mva/bin/micromamba
ST=$S/mva/logs/STATUS_setup
say(){ echo "$(date -u +%H:%M) $*" | tee -a "$ST"; }

say "start"
# Core VCF handling first: everything else depends on it and it is quick.
if [ ! -x $S/mva/mamba/envs/mva/bin/bcftools ]; then
  say "installing core (bcftools samtools htslib bedtools mosdepth)"
  $MM create -y -p $S/mva/mamba/envs/mva -c conda-forge -c bioconda --no-rc \
      python=3.11 bcftools=1.24 samtools=1.24 htslib=1.24 bedtools mosdepth || say "CORE FAILED"
fi
say "core done: $($S/mva/mamba/envs/mva/bin/bcftools --version 2>&1 | head -1)"

# VEP is the high-value piece: consequence plus AlphaMissense/CADD/REVEL plugins.
if [ ! -x $S/mva/mamba/envs/mva/bin/vep ]; then
  say "installing ensembl-vep"
  $MM install -y -p $S/mva/mamba/envs/mva -c conda-forge -c bioconda --no-rc \
      ensembl-vep=116.1 || say "VEP FAILED"
fi
say "vep: $($S/mva/mamba/envs/mva/bin/vep --help 2>&1 | grep -m1 -i 'ensembl-vep\|versions' || echo absent)"

# Alignment and SV, for Arms C and D if the FASTQ is transferred.
if [ ! -x $S/mva/mamba/envs/mva/bin/bwa-mem2 ]; then
  say "installing bwa-mem2 fastp"
  $MM install -y -p $S/mva/mamba/envs/mva -c conda-forge -c bioconda --no-rc \
      bwa-mem2=2.3 fastp || say "ALIGN FAILED"
fi
say "align done"
say "ALL DONE"
