"""MUC1 gene model projected onto a hybrid IGV display contig.

A display contig is ``left flank + allele (motif 1..9) + right flank`` in MUC1
transcript orientation. The bundled flanks are GRCh38 sequence: the left flank is
the reverse complement of the bases just above the bundled VNTR region, the right
flank that of the bases just below it (MUC1 lies on the minus strand). A GRCh38
position ``g`` therefore maps to contig index

* ``vntr_end + left_bp - g`` when ``g > vntr_end`` (left flank), and
* ``left_bp + allele_len + (vntr_start - 1 - g)`` when ``g < vntr_start`` (right flank),

with ``left_bp`` the displayed left-flank width. The transcript holds fewer VNTR
units than GRCh38, so its exon 2 aligns as blocks around the VNTR; the exons that
overlap the VNTR region (or flank a gap containing it) are merged into one exon 2
spanning the allele's repeat array.
"""

from __future__ import annotations

import importlib.resources
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeneModel:
    """One transcript in GRCh38 plus-strand coordinates (1-based, inclusive)."""

    gene: str
    transcript: str
    label: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    cds: tuple[int, int]
    source: str


def bundled_gene_path() -> Path:
    """Path of the bundled MUC1 gene model."""
    ref = importlib.resources.files("muc_one_span.data.annotation").joinpath("muc1_gene.json")
    return Path(str(ref))


def load_gene_model(path: Path | None = None) -> GeneModel:
    """Load a gene model JSON (default: the bundled MUC1 transcript)."""
    data = json.loads((path or bundled_gene_path()).read_text(encoding="utf-8"))
    if data.get("strand") != "-":
        raise ValueError("gene model: the display flanks are in minus-strand orientation")
    exons = tuple(sorted((int(a), int(b)) for a, b in data["exons"]))
    if not exons or any(a > b for a, b in exons):
        raise ValueError("gene model: exons must be nonempty 1-based inclusive intervals")
    cds = (int(data["cds"][0]), int(data["cds"][1]))
    return GeneModel(
        data["gene"], data["transcript"], data["transcript_label"], "-", exons, cds, data["source"]
    )


def parse_region(region: str) -> tuple[int, int]:
    """``chr:start-end`` (1-based inclusive) -> (start, end)."""
    start, end = region.rsplit(":", 1)[1].replace(",", "").split("-")
    return int(start), int(end)


def merged_exons(model: GeneModel, vntr: tuple[int, int]) -> list[tuple[int, int]]:
    """Exons with the blocks around the VNTR merged into one spanning exon."""
    start, end = vntr
    exons = list(model.exons)
    touching = [i for i, (a, b) in enumerate(exons) if a <= end and b >= start]
    if not touching:
        touching = [
            i for i in range(len(exons) - 1) if exons[i][1] < start and exons[i + 1][0] > end
        ]
        touching = [*touching, touching[0] + 1] if touching else []
    if not touching:
        raise ValueError(f"{model.transcript}: no exon spans the VNTR region {start}-{end}")
    first, last = min(touching), max(touching)
    merged = (exons[first][0], exons[last][1])
    if merged[0] >= start or merged[1] <= end:
        raise ValueError(f"{model.transcript}: merged exon does not span the VNTR region")
    return [*exons[:first], merged, *exons[last + 1 :]]


def display_flanks(
    model: GeneModel, vntr: tuple[int, int], margin_bp: int, available: tuple[int, int]
) -> tuple[int, int]:
    """Left/right flank widths that hold the whole transcript plus ``margin_bp``."""
    start, end = vntr
    left = max(b for _, b in model.exons) - end + margin_bp
    right = start - min(a for a, _ in model.exons) + margin_bp
    return min(max(left, 0), available[0]), min(max(right, 0), available[1])


def gene_bed_line(
    model: GeneModel, vntr: tuple[int, int], left_bp: int, allele_len: int, contig: str
) -> str:
    """BED12 transcript feature on one display contig (``+``: transcript orientation)."""
    start, end = vntr

    def project(g: int) -> int:
        if g > end:
            return end + left_bp - g
        if g < start:
            return left_bp + allele_len + (start - 1 - g)
        raise ValueError(f"position {g} lies inside the VNTR region")

    blocks = sorted((project(b), project(a) + 1) for a, b in merged_exons(model, vntr))
    tx_start, tx_end = blocks[0][0], blocks[-1][1]
    thick_start, thick_end = project(model.cds[1]), project(model.cds[0]) + 1
    sizes = ",".join(str(b - a) for a, b in blocks)
    starts = ",".join(str(a - tx_start) for a, _ in blocks)
    return (
        f"{contig}\t{tx_start}\t{tx_end}\t{model.gene}_{model.transcript}\t0\t+\t"
        f"{thick_start}\t{thick_end}\t0\t{len(blocks)}\t{sizes},\t{starts},\n"
    )
