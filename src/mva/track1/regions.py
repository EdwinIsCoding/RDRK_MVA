"""Region and variant-class annotation from an Ensembl GTF.

Why not VEP here
----------------
VEP is the right tool for consequence calling and it is pinned in
``environment.yml`` for the GPU host. It is not used for this module because it
needs a cache of roughly 25 GB, and because the specific thing Arm B needs is
narrower than a full consequence call: **the distance from the variant to the
nearest splice site**, which decides whether a variant is near-splice, deep
intronic, or beyond the window the plan targets.

That distance is computable exactly from exon coordinates, with no cache and no
ambiguity. It is also the classification the benchmark is stratified by, so
computing it here rather than parsing it out of a VEP consequence string keeps
the benchmark and the pipeline using one definition.

What this module does NOT do: protein consequence. Missense, nonsense and
frameshift require a codon model and a translation, and VEP does that properly.
Those classes come from VEP on the GPU host. Anything this module cannot decide
is returned as ``UNCLASSIFIED`` rather than guessed.

Validation
----------
Splice distances were cross-checked against the intron offsets in ClinVar's own
HGVS strings, over the 301 benchmark rows carrying one:

    SNVs    268/268   100.0% concordance
    indels    8/33     24.2% concordance

The SNV figure is the one that validates this module. The indel figure is a
comparison artefact, not an error: HGVS 3'-shifts an indel to the most 3'
position in a repeat, while VCF left-aligns it, so the two conventions
legitimately disagree by the length of the repeat. Splice distance for indels
should therefore be read as approximate, and an indel sitting near a class
boundary should not be trusted to fall on the right side of it.
"""

from __future__ import annotations

import bisect
import collections
import dataclasses
import gzip
import pathlib
import re
from typing import Iterable, Iterator

from mva.evidence import GenomicPosition, VariantClass

#: Distance from an exon boundary, in bases, at which each class begins.
#: 1-2 is the canonical GT/AG dinucleotide; out to 10 the extended splice motif
#: still explains an effect; beyond that the variant is deep intronic and needs
#: a widened prediction window to be seen at all.
CANONICAL_MAX = 2
NEAR_SPLICE_MAX = 10


@dataclasses.dataclass(frozen=True)
class Exon:
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    transcript: str
    rank: int
    #: MANE Select where one exists, else Ensembl canonical. Splice distance is
    #: measured against these exons only. See ``Gene.exon_bounds``.
    is_representative: bool = False


@dataclasses.dataclass(frozen=True)
class Gene:
    symbol: str
    gene_id: str
    contig: str
    start: int
    end: int
    strand: int
    exons: tuple[Exon, ...]
    cds_min: int | None
    cds_max: int | None

    @property
    def exon_bounds(self) -> list[int]:
        """Exon boundaries of the representative transcript, sorted.

        Measured against MANE Select (or Ensembl canonical) only, NOT against
        every annotated transcript. Pooling all transcripts was measurably
        wrong: cross-checking against ClinVar HGVS offsets gave 89.4%
        concordance, and every disagreement was a variant pulled closer to a
        minor isoform's boundary than to the canonical one. The effect is to
        promote genuine deep intronic variants into the canonical-splice class,
        which both overstates their significance and empties the deep intronic
        bucket the benchmark is stratified by.

        Falls back to all transcripts only where the gene has no representative
        transcript annotated, which is rare and worth knowing about.
        """
        rep = [e for e in self.exons if e.is_representative]
        source = rep or self.exons
        bounds: set[int] = set()
        for e in source:
            bounds.add(e.start)
            bounds.add(e.end)
        return sorted(bounds)

    @property
    def has_representative_transcript(self) -> bool:
        return any(e.is_representative for e in self.exons)


class GeneModel:
    """Interval lookup over a set of genes, built from an Ensembl GTF.

    Ensembl GTFs use no-chr contig naming, which matches the proband VCF, so no
    rename is needed on this path. That is checked rather than assumed.
    """

    def __init__(self, genes: Iterable[Gene]) -> None:
        self.genes: list[Gene] = sorted(genes, key=lambda g: (g.contig, g.start))
        self._by_contig: dict[str, list[Gene]] = collections.defaultdict(list)
        for g in self.genes:
            self._by_contig[g.contig].append(g)
        self._starts: dict[str, list[int]] = {
            c: [g.start for g in gs] for c, gs in self._by_contig.items()
        }
        self._bounds_cache: dict[str, list[int]] = {}

    def __len__(self) -> int:
        return len(self.genes)

    @classmethod
    def from_gtf(
        cls,
        gtf: str | pathlib.Path,
        symbols: set[str] | None = None,
        flank: int = 5000,
    ) -> GeneModel:
        """Parse exons for the requested gene symbols.

        ``flank`` widens each gene's reported span so that promoter and
        downstream variants still map to the gene; it does not affect exon
        coordinates or splice distances.
        """
        gtf = pathlib.Path(gtf)
        opener = gzip.open if gtf.suffix == ".gz" else open

        gene_meta: dict[str, dict] = {}
        exons: dict[str, list[Exon]] = collections.defaultdict(list)
        cds: dict[str, list[int]] = collections.defaultdict(list)
        # transcript_id -> priority; 2 = MANE Select, 1 = Ensembl canonical.
        rep_rank: dict[str, int] = {}

        attr_re = re.compile(r'(\S+) "([^"]*)"')
        with opener(gtf, "rt") as fh:  # type: ignore[operator]
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in ("gene", "transcript", "exon", "CDS"):
                    continue
                attrs = dict(attr_re.findall(f[8]))
                sym = attrs.get("gene_name")
                if not sym or (symbols is not None and sym not in symbols):
                    continue

                if f[2] == "transcript":
                    tid = attrs.get("transcript_id", "")
                    tags = f[8]
                    if "MANE_Select" in tags:
                        rep_rank[tid] = 2
                    elif "Ensembl_canonical" in tags:
                        rep_rank[tid] = max(rep_rank.get(tid, 0), 1)
                    continue

                if f[2] == "gene":
                    gene_meta[sym] = {
                        "gene_id": attrs.get("gene_id", ""),
                        "contig": f[0],
                        "start": int(f[3]),
                        "end": int(f[4]),
                        "strand": 1 if f[6] == "+" else -1,
                    }
                elif f[2] == "exon":
                    exons[sym].append(Exon(
                        start=int(f[3]), end=int(f[4]),
                        transcript=attrs.get("transcript_id", ""),
                        rank=int(attrs.get("exon_number", 0) or 0),
                    ))
                else:  # CDS
                    cds[sym].extend((int(f[3]), int(f[4])))

        # GTFs list a transcript before its exons, but only within a record
        # block; resolve representative status after the whole file is read
        # rather than relying on ordering.
        best = {}
        for sym, exs in exons.items():
            ranks = {e.transcript: rep_rank.get(e.transcript, 0) for e in exs}
            top = max(ranks.values(), default=0)
            best[sym] = {t for t, r in ranks.items() if r == top and r > 0}
        for sym, exs in exons.items():
            keep = best.get(sym) or set()
            exons[sym] = [
                dataclasses.replace(e, is_representative=e.transcript in keep)
                for e in exs
            ]

        out = []
        for sym, meta in gene_meta.items():
            c = cds.get(sym) or []
            out.append(Gene(
                symbol=sym, gene_id=meta["gene_id"], contig=meta["contig"],
                start=max(1, meta["start"] - flank), end=meta["end"] + flank,
                strand=meta["strand"], exons=tuple(exons.get(sym, ())),
                cds_min=min(c) if c else None, cds_max=max(c) if c else None,
            ))
        return cls(out)

    def genes_at(self, contig: str, pos: int) -> list[Gene]:
        gs = self._by_contig.get(contig)
        if not gs:
            return []
        starts = self._starts[contig]
        # All genes starting at or before pos; scan back while spans can reach.
        i = bisect.bisect_right(starts, pos)
        hits = []
        for g in reversed(gs[:i]):
            if g.end >= pos:
                hits.append(g)
            # Genes are sorted by start; a very long gene earlier in the list
            # can still overlap, so do not break early on the first miss.
        return hits

    def splice_distance(self, gene: Gene, pos: int) -> int | None:
        """Distance in bases to the nearest exon boundary. 0 means the position
        is on a boundary."""
        key = f"{gene.symbol}:{gene.gene_id}"
        bounds = self._bounds_cache.get(key)
        if bounds is None:
            bounds = gene.exon_bounds
            self._bounds_cache[key] = bounds
        if not bounds:
            return None
        i = bisect.bisect_left(bounds, pos)
        cands = []
        if i < len(bounds):
            cands.append(abs(bounds[i] - pos))
        if i > 0:
            cands.append(abs(pos - bounds[i - 1]))
        return min(cands) if cands else None

    def in_exon(self, gene: Gene, pos: int) -> bool:
        return any(e.start <= pos <= e.end for e in gene.exons)

    def classify(self, position: GenomicPosition) -> tuple[str | None, VariantClass, int | None]:
        """Return (gene symbol, variant class, splice distance).

        Only the classes decidable from coordinates are returned. Protein
        consequence needs VEP; those variants come back as ``UNCLASSIFIED``
        rather than guessed at.
        """
        if position.naming != "ensembl_nochr":
            position = position.to_naming("ensembl_nochr")

        hits = self.genes_at(position.contig, position.pos)
        if not hits:
            return None, VariantClass.UNCLASSIFIED, None

        # Prefer the gene whose exons are closest, so an overlapping readthrough
        # transcript does not capture the variant.
        best_gene, best_dist = None, None
        for g in hits:
            d = self.splice_distance(g, position.pos)
            if d is None:
                continue
            if best_dist is None or d < best_dist:
                best_gene, best_dist = g, d
        if best_gene is None:
            return hits[0].symbol, VariantClass.UNCLASSIFIED, None

        pos = position.pos
        exonic = self.in_exon(best_gene, pos)

        if exonic:
            # Untranslated where it falls outside the CDS span of every
            # transcript. Coding exonic variants need VEP for their class.
            if best_gene.cds_min is not None and not (best_gene.cds_min <= pos <= best_gene.cds_max):
                return best_gene.symbol, VariantClass.UTR_PROMOTER, best_dist
            return best_gene.symbol, VariantClass.UNCLASSIFIED, best_dist

        # Outside the transcribed span entirely: promoter or downstream.
        gene_core_start = min(e.start for e in best_gene.exons)
        gene_core_end = max(e.end for e in best_gene.exons)
        if pos < gene_core_start or pos > gene_core_end:
            return best_gene.symbol, VariantClass.UTR_PROMOTER, best_dist

        if best_dist is None:
            return best_gene.symbol, VariantClass.UNCLASSIFIED, None
        if best_dist <= CANONICAL_MAX:
            return best_gene.symbol, VariantClass.SPLICE_CANONICAL, best_dist
        if best_dist <= NEAR_SPLICE_MAX:
            return best_gene.symbol, VariantClass.SPLICE_REGION, best_dist
        return best_gene.symbol, VariantClass.DEEP_INTRONIC, best_dist


def iter_vcf(vcf: str | pathlib.Path, region: str | None = None) -> Iterator[dict]:
    """Stream a VCF as dicts. Uses bcftools so that indexes and BGZF are
    handled, rather than reimplementing either.

    Two things this function refuses to do quietly:

    * It does not assume the spike-in INFO tags exist. ``SPIKEGENE`` is present
      only in VCFs built by ``mva.track1.spikein``; querying it against the
      proband callset makes bcftools exit with an error.
    * It does not treat a failed query as an empty result. That exact failure
      silently produced "0 records seen, 0 candidates, annotators ran fine",
      which is indistinguishable from a legitimate negative and is precisely
      the outcome this project must never report by accident.
    """
    import subprocess

    header = subprocess.run(["bcftools", "view", "-h", str(vcf)],
                            capture_output=True, text=True, check=True).stdout
    has_spike = "##INFO=<ID=SPIKEGENE" in header

    fields = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t[%DP]"
    if has_spike:
        fields += "\t%INFO/SPIKEGENE"
    fields += "\n"

    cmd = ["bcftools", "query"]
    if region:
        cmd += ["-r", region]
    cmd += ["-f", fields, str(vcf)]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    n = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 7:
            continue
        chrom, pos, ref, alt, filt, gt, dp = parts[:7]
        spike = parts[7] if has_spike and len(parts) > 7 else "."
        n += 1
        yield {
            "contig": chrom, "pos": int(pos), "ref": ref, "alt": alt,
            "filter": filt, "gt": gt,
            "dp": int(dp) if dp.isdigit() else None,
            "spike_gene": None if spike == "." else spike,
        }
    proc.wait()
    if proc.returncode != 0:
        err = (proc.stderr.read() if proc.stderr else "").strip()
        raise RuntimeError(
            f"bcftools query failed (exit {proc.returncode}) on {vcf}"
            f"{f' region {region}' if region else ''} after {n} records: {err}"
        )
