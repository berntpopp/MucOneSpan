"""The bundled MUC1 gene model projected onto display contigs of any allele length."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import pytest

from muc_one_span.hybrid.igv_gene import (
    bundled_gene_path,
    display_flanks,
    gene_bed_line,
    load_gene_model,
    merged_exons,
    parse_region,
)
from tests.unit.hybrid import synth

RD = synth.RD
MODEL = load_gene_model()
VNTR = parse_region(RD.vntr_region)
MARGIN = 500
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {"".join(c): AMINO[i] for i, c in enumerate(itertools.product(BASES, repeat=3))}
STOPS = {codon for codon, amino in CODON.items() if amino == "*"}


def _contig(inner: list[str]) -> tuple[str, list[str]]:
    left, right = display_flanks(
        MODEL, VNTR, MARGIN, (len(RD.flanking_left), len(RD.flanking_right))
    )
    allele = synth.allele(inner)
    contig = RD.flanking_left[-left:] + allele + RD.flanking_right[:right]
    return contig, gene_bed_line(MODEL, VNTR, left, len(allele), "c").rstrip("\n").split("\t")


def _blocks(fields: list[str]) -> list[tuple[int, int]]:
    start = int(fields[1])
    sizes = [int(v) for v in fields[10].rstrip(",").split(",")]
    offsets = [int(v) for v in fields[11].rstrip(",").split(",")]
    return [(start + o, start + o + s) for o, s in zip(offsets, sizes, strict=True)]


@pytest.mark.parametrize("inner", [["X"] * 30, ["X"] * 80, ["A", "B", "X", "C"] * 9])
def test_projected_transcript_is_spliced_and_translates_to_muc1(inner: list[str]) -> None:
    contig, fields = _contig(inner)
    blocks = _blocks(fields)
    assert len(blocks) == len(MODEL.exons) - 1  # exon 2 blocks merged across the VNTR
    for (_, donor), (acceptor, _) in itertools.pairwise(blocks):
        assert contig[donor : donor + 2] == "GT" and contig[acceptor - 2 : acceptor] == "AG"
    thick_start, thick_end = int(fields[6]), int(fields[7])
    cds = "".join(contig[max(a, thick_start) : min(b, thick_end)] for a, b in blocks)
    assert len(cds) % 3 == 0 and cds[:3] == "ATG" and cds[-3:] in STOPS
    protein = "".join(CODON[cds[i : i + 3]] for i in range(0, len(cds), 3))
    assert "*" not in protein[:-1]
    # MUC1 signal peptide start and cytoplasmic tail end (UniProt P15941).
    assert protein.startswith("MTPGTQSPFFLLLLLTVLTVVT") and protein.endswith("AATSANL*")


def test_the_merged_exon_spans_the_vntr_and_the_flanks_hold_the_gene() -> None:
    exons = merged_exons(MODEL, VNTR)
    assert any(a < VNTR[0] and b > VNTR[1] for a, b in exons)
    left, right = display_flanks(MODEL, VNTR, MARGIN, (10**6, 10**6))
    assert left == max(b for _, b in MODEL.exons) - VNTR[1] + MARGIN
    assert right == VNTR[0] - min(a for a, _ in MODEL.exons) + MARGIN
    assert display_flanks(MODEL, VNTR, MARGIN, (100, 200)) == (100, 200)


def test_a_plus_strand_model_is_rejected(tmp_path: Path) -> None:
    data = json.loads(bundled_gene_path().read_text())
    data["strand"] = "+"
    path = tmp_path / "gene.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="minus-strand"):
        load_gene_model(path)
