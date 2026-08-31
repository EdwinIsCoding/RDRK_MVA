#!/usr/bin/env bash
# MVA GPU queue. Waits politely for the SailSwarm labelling job to finish, then
# runs the GPU-dependent Track 1 work in sequence. Nothing here kills or
# interferes with the SailSwarm processes.
set -uo pipefail
S=${REMOTE_SCRATCH:-/scratch0/$GPU_USER}; B=$S/mva/mamba/envs/mva
export PATH=$B/bin:/usr/bin:/bin
export XDG_CACHE_HOME=$S/.cache PIP_CACHE_DIR=$S/.cache/pip TMPDIR=$S/mva/tmp
export TORCH_HOME=$S/.cache/torch HF_HOME=$S/.cache/huggingface
ST=$S/mva/logs/STATUS_gpu
say(){ echo "$(date -u +%H:%M) $*" | tee -a "$ST"; }

say "queued, waiting for the GPU"
while pgrep -f "label_bundle_dart|remote_queue|express_queue|train_student|night_queue" >/dev/null; do
  sleep 120
done
# Confirm the card is actually idle, not merely that the pattern stopped matching.
for i in $(seq 1 30); do
  used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits)
  [ "$used" -lt 1000 ] && break
  say "card still holds ${used} MiB, waiting"
  sleep 120
done
say "GPU free ($(nvidia-smi --query-gpu=memory.used --format=csv,noheader))"

# --- SpliceAI over the widened candidate set -------------------------------
VENV=$S/mva/venv
if [ ! -x $VENV/bin/python ]; then
  say "creating venv on scratch"
  $S/bin/uv venv --python 3.11 $VENV >/dev/null 2>&1 || python3 -m venv $VENV
fi
if ! $VENV/bin/python -c "import spliceai" 2>/dev/null; then
  say "installing spliceai + tensorflow"
  $S/bin/uv pip install --python $VENV/bin/python --no-cache \
      spliceai tensorflow "setuptools<81" biopython pysam >/dev/null 2>&1 \
      || say "spliceai install FAILED"
fi
say "spliceai: $($VENV/bin/python -c 'import spliceai;print("ok")' 2>&1 | tail -1)"
say "gpu visible to tf: $($VENV/bin/python -c 'import tensorflow as tf;print(len(tf.config.list_physical_devices("GPU")))' 2>&1 | tail -1)"
say "ALL DONE"
