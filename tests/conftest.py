"""Shared fixtures. Tests that need data or a GPU host skip themselves rather
than failing, so the suite stays green on the recon machine and the skips are
visible in the report."""
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _complete_gzip(path: pathlib.Path) -> bool:
    """True if the file is a complete BGZF or gzip stream.

    Existence is not enough: a download still in flight leaves a truncated file
    on disk, and a fixture that only checks existence turns that into a
    confusing EOFError deep inside a test rather than a skip.
    """
    if not path.exists() or path.stat().st_size < 1024:
        return False
    try:
        import subprocess
        # gzip -t reads the whole stream and verifies the trailer.
        subprocess.run(["gzip", "-t", str(path)], check=True,
                       capture_output=True, timeout=300)
        return True
    except Exception:
        return False


@pytest.fixture(scope="session")
def repo_root() -> pathlib.Path:
    return ROOT


@pytest.fixture(scope="session")
def gtf_path() -> pathlib.Path:
    p = ROOT / "refs/Homo_sapiens.GRCh38.115.gtf.gz"
    if not _complete_gzip(p):
        pytest.skip("Ensembl GTF missing or incompletely downloaded")
    return p


@pytest.fixture(scope="session")
def background_vcf() -> pathlib.Path:
    p = ROOT / "benchmarks/background/HG002_GRCh38_v4.2.1.vcf.gz"
    if not _complete_gzip(p):
        pytest.skip("GIAB HG002 background VCF missing or incompletely downloaded")
    return p


@pytest.fixture(scope="session")
def panel_tsv() -> pathlib.Path:
    p = ROOT / "config/gene_panels/mitotic_extended.tsv"
    if not p.exists():
        pytest.skip("mitotic panel not built; run scripts/11_build_mitotic_panel.py")
    return p


@pytest.fixture(scope="session")
def gene_model(gtf_path, panel_tsv):
    """Gene model restricted to the panel, so the GTF parse stays fast."""
    import csv as _csv
    from mva.track1.regions import GeneModel
    with panel_tsv.open(newline="") as fh:
        symbols = {r["symbol"] for r in _csv.DictReader(fh, delimiter="\t")}
    return GeneModel.from_gtf(gtf_path, symbols=symbols)
