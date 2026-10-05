"""IGV display reference: mutated units restored, mutations marked on reference bases."""

from __future__ import annotations

from pathlib import Path

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.hybrid.igv_gene import load_gene_model
from muc_one_span.hybrid.igv_tracks import (
    UNIT_COLORS,
    display_allele,
    variant_site,
    write_igv_inputs,
)
from tests.unit.hybrid import synth

FLANK = 500
X = synth.RD.repeats["X"]
UNIT = len(X)
# Mutated unit position: after the 5 fixed pre-repeats and 3 inner X units.
INDEX = len(synth.PRE) + 4
UNIT_START = (INDEX - 1) * UNIT
WILD_TYPE = synth.allele(["X"] * 7)


def _isolated_base() -> int:
    """First unit position whose base differs from both neighbours (unambiguous indel)."""
    return next(i for i in range(1, UNIT - 1) if X[i - 1] != X[i] != X[i + 1] != X[i - 1])


def _classified(unit: str) -> tuple[str, dict]:
    seq = synth.allele(["X"] * 3 + [unit] + ["X"] * 3)
    return seq, classify_sequence(seq, synth.RD)


def test_dupc_insertion_is_left_normalized_to_the_start_of_the_c_run() -> None:
    seq, result = _classified(synth.dupc())
    shown = display_allele("allele_1", seq, result, synth.RD)
    (feature,) = shown.mutations
    assert shown.sequence == WILD_TYPE
    assert feature.name == f"allele_1:repeat_{INDEX}:X:dupC"
    assert (feature.unit_start, feature.unit_end) == (UNIT_START, UNIT_START + UNIT)
    run = UNIT_START + X.rindex("C" * 7)
    # The +C sits before the first C of the run: aligners record it there, so the
    # mark covers the base before and the first C, and the reads are sorted at that C.
    assert (feature.start, feature.end, feature.site) == (run - 1, run + 1, run)
    assert shown.sequence[run - 1] != "C" and shown.sequence[run] == "C"


def test_deletion_marks_the_deleted_reference_base() -> None:
    i = _isolated_base()
    seq, result = _classified(X[:i] + X[i + 1 :])
    shown = display_allele("allele_2", seq, result, synth.RD)
    (feature,) = shown.mutations
    assert shown.sequence == WILD_TYPE
    assert (feature.start, feature.end) == (UNIT_START + i, UNIT_START + i + 1)
    assert feature.site == UNIT_START + i
    assert feature.name == f"allele_2:repeat_{INDEX}:X:del{i + 1}{X[i]}"


def test_unmutated_allele_is_shown_as_its_consensus() -> None:
    seq = synth.allele(["X"] * 4)
    shown = display_allele("allele_1", seq, classify_sequence(seq, synth.RD), synth.RD)
    assert shown.sequence == seq and shown.mutations == []
    assert [u.label for u in shown.units] == classify_sequence(seq, synth.RD)["structure"].split()


def test_variant_site_semantics_and_left_normalization() -> None:
    unit = "ACCCGT"
    substitution = [{"pos": 5, "ref": "G", "alt": "A", "type": "substitution"}]
    assert variant_site(substitution, unit) == (4, 5, 4)
    # An extra C anywhere in the CCC run is recorded before its first C (index 1).
    insertion = [{"pos": 4, "ref": "", "alt": "C", "type": "insertion"}]
    assert variant_site(insertion, unit) == (0, 2, 1)
    # A deleted C of the run is its first C.
    deletion = [{"pos": 4, "ref": "C", "alt": "", "type": "deletion"}]
    assert variant_site(deletion, unit) == (1, 2, 1)
    other = [{"pos": 6, "ref": "", "alt": "A", "type": "insertion"}]
    assert variant_site(other, unit) == (4, 6, 5)
    assert variant_site([], unit) == (0, len(unit), 0)


def test_mutation_without_a_classified_window_is_an_error() -> None:
    result = {"repeats": [], "mutations_detected": [{"repeat_index": 3, "closest_type": "X"}]}
    with pytest.raises(ValueError, match="repeat 3"):
        display_allele("allele_1", "", result, synth.RD)


def test_written_inputs_hold_flanked_display_alleles_rows_and_tracks(tmp_path: Path) -> None:
    normal = synth.allele(["X"] * 4)
    carrier, carrier_result = _classified(synth.dupc())
    files = write_igv_inputs(
        tmp_path,
        {"allele_2": carrier, "allele_1": normal},
        {"allele_1": classify_sequence(normal, synth.RD), "allele_2": carrier_result},
        synth.RD,
        (FLANK, FLANK),
        load_gene_model(),
    )
    lines = files.fasta.read_text().splitlines()
    contigs = dict(zip(lines[::2], lines[1::2], strict=True))
    left, right = synth.RD.flanking_left[-FLANK:], synth.RD.flanking_right[:FLANK]
    assert contigs == {
        ">hybrid_allele_1": left + normal + right,
        ">hybrid_allele_2": left + WILD_TYPE + right,
    }
    assert files.contig_lengths["hybrid_allele_2"] == len(left + WILD_TYPE + right)
    assert [name for name, _ in files.tracks] == ["MUC1 gene", "Repeat units", "Detected mutations"]
    rows = [line.split("\t") for line in files.loci.read_text().splitlines()]
    name = f"allele_2:repeat_{INDEX}:X:dupC"
    assert [r[3] for r in rows] == ["MUC1_VNTR_allele_1", "MUC1_VNTR_allele_2", name]
    assert files.mutation_names == [name]
    assert rows[0][:3] == ["hybrid_allele_1", str(FLANK), str(FLANK + len(normal))]
    # A mutation row is its 1-bp site (the report adds the view context around it).
    site = FLANK + UNIT_START + X.rindex("C" * 7)
    assert rows[2][:3] == ["hybrid_allele_2", str(site), str(site + 1)]
    tracks = dict(files.tracks)
    (bed,) = tracks["Detected mutations"].read_text().splitlines()
    contig, start, end, label = bed.split("\t")
    assert (contig, int(start), int(end)) == ("hybrid_allele_2", site - 1, site + 1)
    assert ">" not in label
    units = [u.split("\t") for u in tracks["Repeat units"].read_text().splitlines()]
    carrier_units = [u for u in units if u[0] == "hybrid_allele_2"]
    assert len(carrier_units) == len(synth.PRE) + 7 + len(synth.POST)
    mutated = carrier_units[INDEX - 1]
    assert mutated[3] == f"{INDEX}:X:dupC" and mutated[8] == UNIT_COLORS["mutation"]
    assert (
        carrier_units[0][8] == UNIT_COLORS["pre"] and carrier_units[-1][8] == UNIT_COLORS["after"]
    )
    assert int(mutated[1]) == FLANK + UNIT_START and int(mutated[2]) == FLANK + UNIT_START + UNIT
    genes = tracks["MUC1 gene"].read_text().splitlines()
    assert [g.split("\t")[0] for g in genes] == ["hybrid_allele_1", "hybrid_allele_2"]
