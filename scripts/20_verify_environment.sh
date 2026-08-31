#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Assert that the environment matches what PROVENANCE.md and environment.yml
# claim, and report GPU visibility. Non-zero exit means the reproducibility
# claim in the submission is false, so this runs first in `make reproduce`.
# ---------------------------------------------------------------------------
set -uo pipefail
FAIL=0
ok()   { printf '  \033[32mOK\033[0m    %-26s %s\n' "$1" "$2"; }
bad()  { printf '  \033[31mFAIL\033[0m  %-26s %s\n' "$1" "$2"; FAIL=1; }
warn() { printf '  \033[33mWARN\033[0m  %-26s %s\n' "$1" "$2"; }

check() {  # name  command  expected-substring  required(yes|no)
  local name=$1 cmd=$2 want=$3 req=${4:-yes} got
  if ! got=$(eval "$cmd" 2>&1 | head -1); then
    if [ "$req" = yes ]; then bad "$name" "not installed"; else warn "$name" "not installed (optional)"; fi
    return
  fi
  if [ -z "$want" ] || [[ "$got" == *"$want"* ]]; then ok "$name" "$got"
  elif [ "$req" = yes ]; then bad "$name" "got: $got (expected to contain '$want')"
  else warn "$name" "got: $got (expected '$want')"; fi
}

echo "== core =="
check "python"    "python3 --version"                  "Python 3.1"
check "bcftools"  "bcftools --version | head -1"       "1.24"
check "samtools"  "samtools --version | head -1"       "1.24"
check "bedtools"  "bedtools --version"                 ""        no
check "pandoc"    "pandoc --version | head -1"         ""        no

echo
echo "== Arm A and B =="
check "vep"       "vep --help 2>&1 | grep -m1 -i version" "" no
check "spliceai"  "spliceai --help 2>&1 | head -1"        "" no

echo
echo "== alignment (needed for Arms C, D, F) =="
check "bwa-mem2"  "bwa-mem2 version"                   ""        no
check "fastp"     "fastp --version"                    ""        no

echo
echo "== Arm C =="
for t in configManta.py delly SURVIVOR AnnotSV ExpansionHunter; do
  check "$t" "command -v $t" "" no
done

echo
echo "== GPU =="
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader \
    | while IFS= read -r line; do ok "nvidia-smi" "$line"; done
  if python3 -c "import torch" 2>/dev/null; then
    python3 - <<'PY'
import torch
print(f"  torch {torch.__version__} cuda_available={torch.cuda.is_available()} "
      f"devices={torch.cuda.device_count()}")
PY
  else
    warn "torch" "not installed; splicing arm will run on CPU"
  fi
else
  warn "nvidia-smi" "no GPU visible; alignment and splicing will be slow or impossible"
fi

echo
echo "== config sanity =="
if [ -f config/config.yaml ]; then
  python3 - <<'PY'
import sys
try:
    import yaml
except ImportError:
    print("  WARN  pyyaml not installed, skipping config check"); sys.exit(0)
c = yaml.safe_load(open("config/config.yaml"))
build = c["genome"]["build"]; naming = c["genome"]["contig_naming"]
assert build == "GRCh38", f"unexpected build {build}"
assert naming == "ensembl_nochr", f"unexpected naming {naming}"
print(f"  OK    genome                     {build} / {naming}")
print(f"  OK    branch                     {c['branch']['id']} ({c['branch']['qualifier']})")
print(f"  OK    seed                       {c['analysis']['random_seed']}")
PY
else
  bad "config/config.yaml" "missing"
fi

echo
echo "== guards =="
[ "$(git config core.hooksPath)" = ".githooks" ] \
  && ok "pre-commit hook" "core.hooksPath=.githooks" \
  || bad "pre-commit hook" "core.hooksPath not set; run: git config core.hooksPath .githooks"
git check-ignore -q data && ok "data/ ignored" "data/ is gitignored" \
  || bad "data/ ignored" "data/ is NOT gitignored"

echo
[ "$FAIL" -eq 0 ] && echo "Environment verified." || echo "Environment does NOT match the reproducibility claim."
exit "$FAIL"
