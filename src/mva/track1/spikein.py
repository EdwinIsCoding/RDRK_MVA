"""Spike benchmark variants into a background genome, per plan section 5.2.

The background is GIAB HG002, a well-characterised consented reference sample.
**Never the proband.** Spiking into the proband's own callset would contaminate
the thing under investigation, and any recall figure derived from it would be
measuring the pipeline's ability to find a variant we put there ourselves in a
genome that already contains the answer we are looking for.

What a recall figure from this harness means, and what it does not
------------------------------------------------------------------
It measures whether the ranking machinery surfaces a known-causal variant when
that variant is present in a realistic genomic background of several million
competing variants. That is a real and necessary test.

It does not measure whether the pipeline would find a *novel* variant, because
every benchmark variant is in ClinVar and most annotation resources know about
them. Plan section 5.3 requires the masked re-run for exactly this reason, and
the masked number is the only one relevant to the proband.
"""

from __future__ import annotations

import csv
import dataclasses
import gzip
import pathlib
import re
import subprocess
import tempfile

from mva.evidence import GenomicPosition, VariantClass


@dataclasses.dataclass(frozen=True)
class BenchmarkVariant:
    """One row of benchmarks/published_mva_variants.tsv or
    benchmarks/splice_mechanism_controls.tsv."""

    gene: str
    clinvar_vcv: str
    hgvs_c: str
    variant_class: VariantClass
    position: GenomicPosition
    clinical_significance: str
    review_status: str
    pmids: tuple[str, ...]
    tier: int
    offset_band: str | None = None
    trait_names: str = ""

    @property
    def label(self) -> str:
        return f"{self.gene}:{self.hgvs_c} [{self.variant_class.value}] {self.clinvar_vcv}"


_CLASS_MAP = {
    "nonsense": VariantClass.NONSENSE,
    "frameshift_or_indel": VariantClass.FRAMESHIFT,
    "missense": VariantClass.MISSENSE,
    "synonymous": VariantClass.SYNONYMOUS,
    "splice_site_canonical": VariantClass.SPLICE_CANONICAL,
    "splice_region_near": VariantClass.SPLICE_REGION,
    "deep_intronic": VariantClass.DEEP_INTRONIC,
    "utr_or_promoter": VariantClass.UTR_PROMOTER,
    "deletion": VariantClass.STRUCTURAL,
    "duplication": VariantClass.STRUCTURAL,
    "microsatellite": VariantClass.REPEAT,
}


def load_benchmark(
    path: str | pathlib.Path = "benchmarks/published_mva_variants.tsv",
    tiers: tuple[int, ...] = (1,),
    naming: str = "ensembl_nochr",
) -> list[BenchmarkVariant]:
    """Load benchmark rows that can actually be spiked.

    Rows without a resolvable GRCh38 position and explicit ref/alt alleles are
    skipped rather than guessed at: a spike-in built from an inferred allele
    would silently test the wrong thing.
    """
    out: list[BenchmarkVariant] = []
    p = pathlib.Path(path)
    if not p.exists():
        return out

    with p.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tier = int(row.get("benchmark_tier", 1) or 1)
            if tiers and tier not in tiers:
                continue
            chrom, pos = row.get("chrom_nochr", ""), row.get("pos_grch38", "")
            ref, alt = row.get("ref", ""), row.get("alt", "")
            if not (chrom and pos and ref and alt):
                continue
            # SPDI alleles are given on the plus strand with a 0-based
            # interbase start; ClinVar's variation_loc start is 1-based, which
            # is the convention a VCF wants.
            try:
                gp = GenomicPosition(
                    build="GRCh38", naming=naming, contig=chrom,
                    pos=int(pos), ref=ref, alt=alt,
                )
            except Exception:
                continue

            out.append(BenchmarkVariant(
                gene=row["gene"],
                clinvar_vcv=row.get("clinvar_vcv", ""),
                hgvs_c=row.get("hgvs_c", ""),
                variant_class=_CLASS_MAP.get(row.get("variant_class", ""),
                                             VariantClass.UNCLASSIFIED),
                position=gp,
                clinical_significance=row.get("clinical_significance", ""),
                review_status=row.get("review_status", ""),
                pmids=tuple(x for x in row.get("pmids", "").split(";") if x),
                tier=tier,
                offset_band=row.get("offset_band"),
                trait_names=row.get("trait_names", ""),
            ))
    return out


def detect_naming(vcf: str | pathlib.Path) -> str:
    """Read the contig naming convention out of a VCF header, rather than
    assuming it. Mixing the two conventions is the single most likely way to
    produce a spike-in that lands nowhere."""
    header = subprocess.run(["bcftools", "view", "-h", str(vcf)],
                            capture_output=True, text=True, check=True).stdout
    return "ucsc_chr" if re.search(r"^##contig=<ID=chr", header, re.M) else "ensembl_nochr"


def spike(
    background_vcf: str | pathlib.Path,
    variants: list[BenchmarkVariant],
    out_vcf: str | pathlib.Path,
    sample: str = "SPIKED",
    genotype: str = "0/1",
) -> pathlib.Path:
    """Insert variants into a background VCF as heterozygous calls.

    Heterozygous by default: MVA is recessive and the case of interest is a
    compound heterozygote, so a spiked homozygote would be an easier problem
    than the real one.
    """
    background_vcf = pathlib.Path(background_vcf)
    out_vcf = pathlib.Path(out_vcf)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    naming = detect_naming(background_vcf)
    records = []
    for v in variants:
        p = v.position.to_naming(naming)  # type: ignore[arg-type]
        records.append(
            f"{p.contig}\t{p.pos}\t{v.clinvar_vcv or '.'}\t{p.ref}\t{p.alt}\t"
            f"500\tPASS\tSPIKED=1;SPIKEGENE={v.gene};SPIKECLASS={v.variant_class.value}\t"
            f"GT:AD:DP:GQ\t{genotype}:20,20:40:99"
        )

    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False,
                                     dir=out_vcf.parent) as tmp:
        # A minimal header; bcftools will reconcile contigs against the
        # background when the two are concatenated.
        tmp.write("##fileformat=VCFv4.2\n")
        tmp.write('##INFO=<ID=SPIKED,Number=0,Type=Flag,Description="spiked-in positive control">\n')
        tmp.write('##INFO=<ID=SPIKEGENE,Number=1,Type=String,Description="gene of the spiked variant">\n')
        tmp.write('##INFO=<ID=SPIKECLASS,Number=1,Type=String,Description="variant class of the spiked variant">\n')
        tmp.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        tmp.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
        tmp.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        tmp.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        header = subprocess.run(["bcftools", "view", "-h", str(background_vcf)],
                                capture_output=True, text=True, check=True).stdout
        for line in header.splitlines():
            if line.startswith("##contig"):
                tmp.write(line + "\n")
        tmp.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for r in sorted(records, key=lambda r: (r.split("\t")[0], int(r.split("\t")[1]))):
            tmp.write(r + "\n")
        spike_path = pathlib.Path(tmp.name)

    try:
        subprocess.run(["bcftools", "sort", "-Oz", "-o", str(out_vcf) + ".spikes.vcf.gz",
                        str(spike_path)], check=True, capture_output=True)
        subprocess.run(["bcftools", "index", "-t", str(out_vcf) + ".spikes.vcf.gz"],
                       check=True, capture_output=True)
        # Reheader the background to the same sample name so the two can merge.
        subprocess.run(
            ["bcftools", "concat", "-a", "-Oz", "-o", str(out_vcf),
             str(background_vcf), str(out_vcf) + ".spikes.vcf.gz"],
            check=True, capture_output=True)
        subprocess.run(["bcftools", "index", "-t", str(out_vcf)],
                       check=True, capture_output=True)
    finally:
        spike_path.unlink(missing_ok=True)
    return out_vcf


def mask_clinvar(
    variants: list[BenchmarkVariant],
    clinvar_ids: set[str] | None = None,
) -> set[str]:
    """Return the set of ClinVar accessions to hide, for the leakage control of
    plan section 5.3.

    Every benchmark variant is in ClinVar, so a pipeline that consults ClinVar
    is being scored on its ability to look up the answer. The masked recall is
    the only number that speaks to finding something novel, which is what the
    proband needs. Both numbers are reported; only the masked one is used to
    make claims.
    """
    ids = {v.clinvar_vcv.split(".")[0] for v in variants if v.clinvar_vcv}
    if clinvar_ids:
        ids |= {c.split(".")[0] for c in clinvar_ids}
    return ids


def count_records(vcf: str | pathlib.Path) -> int:
    """Record count via the index where possible, falling back to a stream."""
    try:
        out = subprocess.run(["bcftools", "index", "-n", str(vcf)],
                             capture_output=True, text=True, check=True).stdout
        return int(out.strip())
    except (subprocess.CalledProcessError, ValueError):
        opener = gzip.open if str(vcf).endswith(".gz") else open
        with opener(vcf, "rt") as fh:  # type: ignore[operator]
            return sum(1 for line in fh if not line.startswith("#"))
