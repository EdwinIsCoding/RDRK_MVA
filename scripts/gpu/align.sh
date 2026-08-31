#!/usr/bin/env bash
# FASTQ to CRAM, unblocking Arms C (structural) and D (mosaic).
#
# Sequencing note: waits for VEP to finish before starting, because bwa-mem2
# holds the ~30 GB human index resident and VEP is already using 12 forks. The
# node has 60 GB total, so running both would swap.
#
# What this will and will not deliver, so the write-up does not overclaim:
#   WILL  read-level verification of candidate variants, coverage over the
#         panel, structural variant calling, mosaic re-genotyping.
#   WILL NOT phase the BUB1B compound heterozygote. Read-backed phasing only
#         resolves variants within a fragment, and HaplotypeCaller's PGT/PID
#         tags already found zero phasing groups in BUB1B, which means the two
#         alleles are further apart than a read pair. No amount of short-read
#         depth changes that; it needs long reads or parental samples.
set -uo pipefail
S=${REMOTE_SCRATCH:-/scratch0/$GPU_USER}
H=$S/mva/mamba/envs/htslib/bin
A=$S/mva/mamba/envs/mva/bin
export PATH=$H:$A:/usr/bin:/bin
export TMPDIR=$S/mva/tmp
ST=$S/mva/logs/STATUS_align
say(){ echo "$(date -u +%H:%M) $*" | tee -a "$ST"; }
REF=$S/mva/refs/GRCh38.fa
SM=WGS_EX2312012

say "queued; waiting for VEP and the FASTQ pull"
for i in $(seq 1 480); do
  vep_done=$(grep -c "full VEP done\|full VEP FAILED" $S/mva/logs/STATUS_vep 2>/dev/null || echo 0)
  dl_done=$(grep -c "download done" $S/mva/logs/STATUS_hfpull 2>/dev/null || echo 0)
  [ "$vep_done" -ge 1 ] && [ "$dl_done" -ge 1 ] && break
  sleep 60
done
say "prerequisites met (vep=$vep_done download=$dl_done)"

if [ ! -s $REF ]; then
  say "decompressing reference"
  gzip -dc $S/mva/refs/GRCh38.fa.gz > $REF && samtools faidx $REF
fi
if [ ! -s ${REF}.bwt.2bit.64 ]; then
  say "building bwa-mem2 index (~30 GB RAM, ~1h)"
  bwa-mem2 index $REF 2>&1 | tail -2 | tee -a "$ST"
fi
say "index ready"

D=$S/mva/data
for L in L001 L002 L003 L004; do
  OUT=$S/mva/results/${SM}_${L}.bam
  [ -s "$OUT" ] && { say "$L cached"; continue; }
  R1=$D/${SM}_HGWCNDSX7_S16_${L}_R1_001.fastq.gz
  R2=$D/${SM}_HGWCNDSX7_S16_${L}_R2_001.fastq.gz
  [ -s "$R1" ] || { say "$L FASTQ missing, skipping"; continue; }
  say "aligning $L"
  RG="@RG\tID:${SM}_${L}\tSM:${SM}\tLB:${SM}\tPL:ILLUMINA\tPU:HGWCNDSX7.${L}"
  bwa-mem2 mem -t 14 -R "$RG" -K 100000000 -Y $REF "$R1" "$R2" 2>>$S/mva/logs/align_${L}.err \
    | samtools sort -@ 4 -m 1500M -T $S/mva/tmp/sort_${L} -o "$OUT" - \
    && samtools index -@ 4 "$OUT" && say "$L done ($(du -h "$OUT" | cut -f1))"
done

say "merging and marking duplicates"
samtools merge -@ 8 -f -o $S/mva/tmp/${SM}.merged.bam $S/mva/results/${SM}_L00*.bam \
  && samtools collate -@ 8 -O -u $S/mva/tmp/${SM}.merged.bam \
   | samtools fixmate -@ 4 -m -u - - \
   | samtools sort -@ 8 -m 1500M -T $S/mva/tmp/md -u - \
   | samtools markdup -@ 8 --write-index -f $S/mva/results/${SM}.dupmetrics.txt \
       - $S/mva/results/${SM}.bam \
  && say "BAM ready: $(du -h $S/mva/results/${SM}.bam | cut -f1)"

say "coverage over the panels"
mosdepth -t 4 --by $S/mva/refs/panels/mva_known.nochr.bed --fast-mode \
  $S/mva/results/${SM}.mva $S/mva/results/${SM}.bam && say "mosdepth done"
say "ALL DONE"
