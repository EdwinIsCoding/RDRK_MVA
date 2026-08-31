# Running on the UCL GPU node

Node `${GPU_NODE:?set in scripts/gpu/.local.sh}` via jump host `${JUMP_HOST:?set in scripts/gpu/.local.sh}`, user `${GPU_USER}`.

```bash
SSHCMD="ssh -i ${SSH_KEY:-~/.ssh/id_ed25519} -o IdentitiesOnly=yes -o IdentityAgent=none \
  -o BatchMode=yes -J ${GPU_USER}@${JUMP_HOST:?set in scripts/gpu/.local.sh}"
$SSHCMD ${GPU_USER}@${GPU_NODE:?set in scripts/gpu/.local.sh} 'bash -s' <<'REMOTE' | grep -v VBoxManage
  ...
REMOTE
```

The login shell is tcsh and prints a VirtualBox error on every connection, so
every command goes through `bash -s` and the output is filtered.

## Layout

Everything lives under `${REMOTE_SCRATCH:-/scratch0/$GPU_USER}/mva`. **Scratch is wiped when the
booking ends**, so results must be pulled before then. AFS home holds code only.

## Findings that shaped these scripts

**Bandwidth is asymmetric by roughly 40x.** The node pulls from Ensembl at about
40 MB/s; pushing from the laptop over the jump host runs at about 1 MB/s. So the
VEP cache (26 GB) and reference (762 MB) are downloaded on the node, and only the
301 MB VCF is pushed. The same arithmetic makes pushing the 85 GB FASTQ
impractical at roughly 24 hours.

**VEP needs the conda perl first on PATH.** The node's `/opt/ucl/bin/perl`
shadows it and VEP then cannot find BioPerl. Every script exports
`PATH=$ENV/bin:/usr/bin:/bin` and unsets `PERL5LIB`.

**htslib tools need their own environment.** Installing bcftools alongside VEP
let the solver satisfy VEP's perl tree by downgrading bcftools to 1.6 from 2017,
which links `libcrypto.so.1.0.0` and does not run. Worse for reproducibility, it
silently differs from the 1.24 used locally. They now live in a separate env.

**The GPU is shared.** `gpu_queue.sh` waits on the SailSwarm process patterns and
then confirms the card actually reports under 1 GB used, because the pattern
ceasing to match is not the same as the card being free. It never kills anything.
