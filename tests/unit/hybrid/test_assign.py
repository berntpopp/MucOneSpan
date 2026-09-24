"""S5/S6 edit-distance competition against ladder-flanked allele drafts."""

from __future__ import annotations

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.assign import OFF_TARGET, assign_read, assign_reads, hybrid_references
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")

# S = HybridSettings() supplies margin/max_error_rate/flank_bp everywhere below, so no
# assign tunable is duplicated as a bare literal in this test module (owner directive:
# fully config-driven, no magic numbers).
S = HybridSettings()
UNIT = synth.RD.repeat_length_bp  # the bundled dictionary's repeat-unit length (60 bp)
A = ["X"] * 10 + ["A", "A", "B"] + ["X"] * 20 + ["G", "A", "B"] + ["X"] * 10
B = ["X"] * 10 + ["A", "A", "B"] + ["X"] * 21 + ["G", "A", "B"] + ["X"] * 10  # one unit longer
REFS = hybrid_references(
    {"allele_1": synth.allele(A), "allele_2": synth.allele(B)}, synth.RD, S.assign_flank_bp
)
TRACT_START = (5 + 13) * UNIT  # pre-repeats 5, then 10 X, then A A B
TRACT_END_B = (5 + 13 + 21) * UNIT  # the 21-X tract in B


def test_fragments_spanning_the_difference_are_assigned_correctly() -> None:
    wrong = decided = 0
    for truth, inner in (("allele_1", A), ("allele_2", B)):
        seq = synth.allele(inner)
        for i, start in enumerate(range(TRACT_START - 400, TRACT_START - 100, 20)):
            piece = seq[start : TRACT_END_B + 400]
            frag = synth.reads(piece, 1, err=0.02, seed=100 + i, flank_bp=0)[0].seq
            got = assign_read(
                rc(frag) if i % 2 else frag, REFS, S.assign_margin, S.assign_max_error_rate
            )
            decided += got.allele is not None
            wrong += got.allele not in (None, truth)
    assert wrong == 0 and decided == 30


def test_fragment_inside_shared_tract_is_undecided() -> None:
    frag = synth.allele(A)[TRACT_START : TRACT_START + 15 * UNIT]
    got = assign_read(frag, REFS, S.assign_margin, S.assign_max_error_rate)
    assert got.allele is None


def test_off_target_read_is_not_assigned() -> None:
    reads = [ReadRecord("junk", "ACGT" * 400, "5" * 1600)]
    out = assign_reads(reads, REFS, S)
    assert len(out[OFF_TARGET]) == 1 and not out["allele_1"] and not out["allele_2"]


def test_single_reference_assigns_on_target_and_guards_off_target() -> None:
    single = {"allele_1": REFS["allele_1"]}
    frag = synth.allele(A)[600:2400]
    assert assign_read(frag, single, S.assign_margin, S.assign_max_error_rate).allele == "allele_1"
    off = assign_read("ACGT" * 400, single, S.assign_margin, S.assign_max_error_rate)
    assert off.allele == OFF_TARGET


def test_short_reads_are_ignored() -> None:
    out = assign_reads([ReadRecord("s", synth.allele(A)[:500], "5" * 500)], REFS, S)
    assert sum(len(v) for v in out.values()) == 0


def test_trim_to_draft_drops_ladder_flanks() -> None:
    from muc_one_span.hybrid.assign import trim_to_draft

    flank_bp = S.assign_flank_bp
    draft = synth.allele(A)
    ref = REFS["allele_1"]
    read = ref[flank_bp - 200 : flank_bp + 1200]  # 200 bp flank + 1200 bp of the draft
    assert trim_to_draft(read, ref, flank_bp, len(draft)) == draft[:1200]
    assert trim_to_draft(ref[:400], ref, flank_bp, len(draft)) == ""
