#!/usr/bin/env bash
# Pull the gated dataset on the node. Token on stdin only: never written to
# disk, never in argv (world-readable via /proc), never echoed.
set -uo pipefail
read -r TOK
S=${REMOTE_SCRATCH:-/scratch0/$GPU_USER}
export HF_HOME=$S/mva/.hf XDG_CACHE_HOME=$S/.cache TMPDIR=$S/mva/tmp
ST=$S/mva/logs/STATUS_hfpull
say(){ echo "$(date -u +%H:%M) $*" | tee -a "$ST"; }
say "download start"
env HF_TOKEN="$TOK" $S/mva/venv/bin/python - <<'PY'
import os, time
from huggingface_hub import hf_hub_download
repo = "SageBio/mva-hackathon-2026-data"
tok = os.environ["HF_TOKEN"]
dest = "${REMOTE_SCRATCH:-/scratch0/$GPU_USER}/mva/data"
files = [f"WGS_EX2312012_HGWCNDSX7_S16_L00{l}_R{r}_001.fastq.gz"
         for l in (1,2,3,4) for r in (1,2)]
for f in files:
    t0 = time.time()
    p = hf_hub_download(repo, f, repo_type="dataset", token=tok,
                        local_dir=dest, resume_download=True)
    sz = os.path.getsize(p) / 1e9
    print(f"{time.strftime('%H:%M')} {f}: {sz:.1f} GB in {time.time()-t0:.0f}s "
          f"({sz*1000/max(1,time.time()-t0):.0f} MB/s)", flush=True)
PY
say "download done"
