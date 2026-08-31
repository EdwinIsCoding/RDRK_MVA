# ---------------------------------------------------------------------------
# MVA Hackathon 2026 reproducibility container.
#
# Base is the CUDA 12.4 runtime so that SpliceAI and Pangolin can use the
# RTX 6000. The image runs fine without a GPU: the splicing arm falls back to
# CPU, slowly. Build and run on the linux-64 host, not the darwin-arm64 recon
# machine (DATA_CARD.md section 7).
#
#   docker build -t mva2026:0.1.0 .
#   docker run --rm --gpus all -v "$PWD:/work" -w /work mva2026:0.1.0 make verify
#
# data/ is deliberately NOT copied into the image. Mount it read-only at run
# time so that a published image can never contain patient data:
#   -v "$PWD/data:/work/data:ro"
# ---------------------------------------------------------------------------
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/opt/conda/envs/mva2026/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 ca-certificates curl git make procps tini \
    && rm -rf /var/lib/apt/lists/*

# Micromamba, pinned. Verify the checksum before changing the version.
ARG MICROMAMBA_VERSION=2.0.5
RUN curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/${MICROMAMBA_VERSION}" \
      | tar -xvj -C /usr/local bin/micromamba

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml -p /opt/conda/envs/mva2026 \
    && micromamba clean --all --yes

# Python package layer, separate so dependency changes do not rebuild conda.
COPY pyproject.toml /tmp/pyproject.toml
COPY src /tmp/src
RUN /opt/conda/envs/mva2026/bin/pip install --no-cache-dir /tmp

WORKDIR /work
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]
