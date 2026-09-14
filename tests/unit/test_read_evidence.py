"""Deterministic evidence tests with known record source and anchor coordinates."""

import pytest

from muc_one_span.read_evidence import assign_span, spanning_evidence

LEFT = "ACGTTGCACTGA"
RIGHT = "GTCCATAGTGAC"


def test_spanning_counts_distinct_records_despite_colliding_names():
    records = [("same", LEFT + "A" * 96 + RIGHT), ("same", LEFT + "A" * 156 + RIGHT)]
    evidence = spanning_evidence(records, LEFT, RIGHT, max_edits=0)
    assert [row.record_index for row in evidence] == [0, 1]
    assert [row.read_name for row in evidence] == ["same", "same"]
    assert [(row.span_min, row.span_max) for row in evidence] == [(120, 120), (180, 180)]
    assert [assign_span(row, {"hap1": 120, "hap2": 180}, 5) for row in evidence] == [
        ("hap1",),
        ("hap2",),
    ]


@pytest.mark.parametrize("left", ["ACGTCGCACTGA", "ACGTTTGCACTGA", "ACGTGCACTGA"])
def test_substitution_insertion_deletion_anchor_errors_are_retained(left):
    sequence = left + "A" * 96 + RIGHT
    exact = spanning_evidence([("r", sequence)], LEFT, RIGHT, max_edits=0)[0]
    tolerant = spanning_evidence([("r", sequence)], LEFT, RIGHT, max_edits=1)[0]
    assert exact.status == "partial_anchor_pair"
    assert tolerant.status == "spanning"
    assert tolerant.anchor_edits == 1
    assert tolerant.span_min <= len(sequence) <= tolerant.span_max


def test_reverse_orientation_has_same_span():
    sequence = LEFT + "A" * 96 + RIGHT
    reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    row = spanning_evidence([("r", reverse)], LEFT, RIGHT, max_edits=0)[0]
    assert (row.orientation, row.span_min, row.span_max) == ("reverse", 120, 120)


def test_partial_and_no_anchor_reads_are_not_assigned():
    rows = spanning_evidence([("r", LEFT + "A" * 80), ("r", "A" * 80)], LEFT, RIGHT)
    assert [row.status for row in rows] == ["partial_anchor_pair", "no_anchor_pair"]
    assert all(assign_span(row, {"hap1": 120}, 100) == () for row in rows)


def test_repeated_anchor_loci_are_ambiguous():
    row = spanning_evidence([("r", LEFT + LEFT + "A" * 96 + RIGHT)], LEFT, RIGHT)[0]
    assert row.status == "ambiguous_anchor_pair"
    assert row.span_min is None


def test_equal_length_and_midpoint_assignment_remain_ambiguous():
    row = spanning_evidence([("r", LEFT + "A" * 126 + RIGHT)], LEFT, RIGHT)[0]
    assert assign_span(row, {"hap1": 150, "hap2": 150}, 5) == ("hap1", "hap2")
    assert assign_span(row, {"hap1": 120, "hap2": 180}, 30) == ("hap1", "hap2")
    assert assign_span(row, {"hap1": 120, "hap2": 180}, 29) == ()


@pytest.mark.parametrize(
    "left,right,edits", [("", RIGHT, 1), (LEFT, "", 1), (LEFT, RIGHT, -1), (LEFT, RIGHT, 12)]
)
def test_invalid_anchor_configuration_rejected(left, right, edits):
    with pytest.raises(ValueError):
        spanning_evidence([], left, right, max_edits=edits)


@pytest.mark.parametrize("max_edits", [1, 2])
def test_tolerant_anchor_hits_agree_with_exhaustive_substring_distance(max_edits):
    from itertools import product

    from muc_one_span.read_evidence import _anchor_hits

    def distance(a, b):
        row = list(range(len(b) + 1))
        for i, base in enumerate(a, 1):
            new = [i]
            for j, other in enumerate(b, 1):
                new.append(min(row[j] + 1, new[-1] + 1, row[j - 1] + (base != other)))
            row = new
        return row[-1]

    for letters in product("AC", repeat=6):
        sequence = "".join(letters)
        for anchor in ("ACA", "AAC", "CCC"):
            candidates = [
                (start, end, distance(anchor, sequence[start:end]))
                for start in range(len(sequence))
                for end in range(start + 1, len(sequence) + 1)
            ]
            best = min(item[2] for item in candidates)
            expected = {item for item in candidates if item[2] == best and best <= max_edits}
            assert set(_anchor_hits(sequence, anchor, max_edits)) == expected


def test_negative_assignment_distance_rejected():
    row = spanning_evidence([("r", LEFT + RIGHT)], LEFT, RIGHT)[0]
    with pytest.raises(ValueError):
        assign_span(row, {"hap1": 24}, -1)


@pytest.mark.parametrize("gap", [0, 1, 2, 3])
def test_known_source_gap_reads_keep_length_assignment_limits(gap):
    records = [
        ("collision", LEFT + "A" * 96 + RIGHT),
        ("collision", LEFT + "C" * (96 + 60 * gap) + RIGHT),
        ("collision", LEFT + "A" * 97 + RIGHT),
        ("collision", "A" * 80 + RIGHT),
    ]
    rows = spanning_evidence(records, LEFT, RIGHT)
    candidates = {"hap1": 120, "hap2": 120 + 60 * gap}
    if gap == 0:
        assert [assign_span(r, candidates, 2) for r in rows] == [
            ("hap1", "hap2"),
            ("hap1", "hap2"),
            ("hap1", "hap2"),
            (),
        ]
    else:
        assert [assign_span(r, candidates, 2) for r in rows] == [
            ("hap1",),
            ("hap2",),
            ("hap1",),
            (),
        ]


def test_two_valid_orientations_are_rejected():
    forward = LEFT + "A" * 96 + RIGHT
    reverse = forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    row = spanning_evidence([("r", forward + reverse)], LEFT, RIGHT)[0]
    assert row.status == "ambiguous_anchor_pair"
