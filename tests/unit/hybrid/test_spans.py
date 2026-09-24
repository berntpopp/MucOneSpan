"""S1 anchor search, orientation and read categories on synthetic reads."""

from __future__ import annotations

import sys

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.hybrid.align import infix_hit
from muc_one_span.hybrid.spans import Anchors, ReadRecord, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD, S)


def test_dupc_template_is_classified() -> None:
    seq = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    names = {
        (m["mutation_name"], m.get("template_match"))
        for m in classify_sequence(seq, synth.RD)["mutations_detected"]
    }
    assert ("dupC", True) in names


def test_missing_edlib_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "edlib", None)
    with pytest.raises(ImportError, match="hybrid"):
        infix_hit("ACGT", "ACGT", 0)


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


def test_mutated_motif1_falls_back_to_flank_anchor() -> None:
    seq = synth.allele(["X"] * 30)
    broken = seq[:5] + "T" * 40 + seq[45:]  # destroy most of motif 1
    read_seq = synth.RD.flanking_left[-40:] + broken + synth.RD.flanking_right[:40]
    cats = categorize_reads(
        [ReadRecord("m", read_seq, "5" * len(read_seq))], ANCH, HybridSettings(anchor_max_edits=6)
    )
    assert len(cats.spanning) == 1 and cats.spanning[0].anchor_basis == "flank"
