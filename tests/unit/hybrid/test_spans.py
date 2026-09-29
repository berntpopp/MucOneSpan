"""S1 anchor search, orientation and read categories on synthetic reads."""

from __future__ import annotations

from muc_one_span.classify import classify_sequence
from muc_one_span.hybrid.spans import Anchors, ReadRecord, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid.synth import LAYOUT

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S, LAYOUT)


def test_dupc_template_is_classified() -> None:
    seq = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    names = {
        (m["mutation_name"], m.get("template_match"))
        for m in classify_sequence(seq, synth.RD)["mutations_detected"]
    }
    assert ("dupC", True) in names


def test_spanning_reads_oriented_and_trimmed_to_motif1_motif9() -> None:
    seq = synth.allele(["X"] * 30)
    cats = categorize_reads(synth.reads(seq, 20, err=0.02, seed=1), ANCH, S)
    assert len(cats.spanning) == 20
    for sp in cats.spanning:
        assert abs(sp.length - len(seq)) <= 40
        assert sp.anchor_basis == "motif"
    assert {sp.strand for sp in cats.spanning} == {"+", "-"}


def test_fragments_and_offtarget_are_not_spanning() -> None:
    seq = synth.allele(["X"] * 30)
    frag = ReadRecord("frag", seq[:900], "5" * 900)
    junk = ReadRecord("junk", "ACGT" * 300, "5" * 1200)
    cats = categorize_reads([frag, junk], ANCH, S)
    assert cats.counts() == {
        "spanning": 0,
        "left_anchored": 1,
        "right_anchored": 0,
        "internal_or_offtarget": 1,
    }


def test_left_anchored_fragment_ending_at_motif_boundary_does_not_crash() -> None:
    """A genomic fragment read that ends exactly where motif 1 ends collapses the
    remaining search window (for motif 9 / the right flank) to an empty string --
    the real-world trigger for edlib's ``locations`` None-start quirk (see
    ``hybrid/align.py::infix_hit``, reproduced directly in ``test_align.py``). It
    must be categorised as left-anchored, not crash.
    """
    read_seq = synth.RD.flanking_left[-40:] + synth.RD.repeats["1"]
    cats = categorize_reads([ReadRecord("frag_left", read_seq, "5" * len(read_seq))], ANCH, S)
    assert cats.counts() == {
        "spanning": 0,
        "left_anchored": 1,
        "right_anchored": 0,
        "internal_or_offtarget": 0,
    }


def test_right_anchored_fragment_does_not_crash() -> None:
    """A fragment carrying only motif 9 (and its downstream flank), with no motif 1
    or left-flank content, must be categorised as right-anchored, not crash.
    """
    read_seq = synth.RD.repeats["9"] + synth.RD.flanking_right[:40]
    cats = categorize_reads([ReadRecord("frag_right", read_seq, "5" * len(read_seq))], ANCH, S)
    assert cats.counts() == {
        "spanning": 0,
        "left_anchored": 0,
        "right_anchored": 1,
        "internal_or_offtarget": 0,
    }


def test_internal_fragment_with_neither_anchor_does_not_crash() -> None:
    """A fragment drawn purely from mid-array repeat content, with neither motif 1
    nor motif 9 nor either flank, must fall to internal/off-target, not crash.
    """
    read_seq = synth.RD.repeats["X"] * 5
    cats = categorize_reads([ReadRecord("frag_internal", read_seq, "5" * len(read_seq))], ANCH, S)
    assert cats.counts() == {
        "spanning": 0,
        "left_anchored": 0,
        "right_anchored": 0,
        "internal_or_offtarget": 1,
    }


def test_mutated_motif1_falls_back_to_flank_anchor() -> None:
    seq = synth.allele(["X"] * 30)
    broken = seq[:5] + "T" * 40 + seq[45:]  # destroy most of motif 1
    read_seq = synth.RD.flanking_left[-40:] + broken + synth.RD.flanking_right[:40]
    cats = categorize_reads(
        [ReadRecord("m", read_seq, "5" * len(read_seq))], ANCH, HybridSettings(anchor_max_edits=6)
    )
    assert len(cats.spanning) == 1 and cats.spanning[0].anchor_basis == "flank"
