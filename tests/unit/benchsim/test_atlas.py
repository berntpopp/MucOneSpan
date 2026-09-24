"""benchsim.atlas: reason keys, expected-INCONCLUSIVE split and reason x stratum tables."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import pytest

from muc_one_span.benchsim.atlas import (
    ANY_REASON,
    UNRECORDED,
    build_atlas,
    case_reason_keys,
    expected_conditions,
    reason_key,
    render_atlas,
)
from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG

A = DEFAULT_BENCH_CONFIG.atlas
GATE = A.min_resolvable_depth


@pytest.mark.parametrize(
    ("text", "key"),
    [
        (
            "Allele 1: 7 primary alignments, below the per-allele depth gate (30).",
            "# primary alignments, below the per-allele depth gate (#)",
        ),
        (
            "Allele 2: 12 primary alignments, below the per-allele depth gate (30).",
            "# primary alignments, below the per-allele depth gate (#)",
        ),
        (
            "Allele 2: Observed sequence variant (dupC at repeat 23) is inconclusive "
            "(frameshift not established; localization ambiguous).",
            "observed sequence variant (<variant>) is inconclusive "
            "(frameshift not established; localization ambiguous)",
        ),
        (
            "Allele 1: allele selection unresolved (unresolved_bimodal; secondary mode "
            "fraction 0.31).",
            "allele selection unresolved (unresolved_bimodal; secondary mode fraction #)",
        ),
        (
            "Allele 1: allele selection unresolved (unresolved_bimodal; secondary mode "
            "fraction None).",
            "allele selection unresolved (unresolved_bimodal; secondary mode fraction #)",
        ),
        (
            "High number of ambiguous consensus bases (57) detected.",
            "high number of ambiguous consensus bases (#) detected",
        ),
        (
            "  Allele_2:   Genotype  phase is unphased or conflicting. ",
            "genotype phase is unphased or conflicting",
        ),
        ("model r1041 kept", "model r1041 kept"),
    ],
)
def test_reason_key_normalises_labels_numbers_and_variants(text: str, key: str) -> None:
    assert reason_key(text) == key


def _row(i: int, **kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "sample": f"s{i}",
        "profile": "hifi_amplicon",
        "decision": "INCONCLUSIVE",
        "depth": GATE,
        "realized_min_allele_depth": GATE,
        "smear": 0.05,
        "chimera": 0.01,
        "delta_class": "1",
        "event_position": None,
        "clinical_reasons": [],
        "reconstruction_flags": [],
        "reasons_recorded": True,
    }
    return row | kw


def test_case_reason_keys_are_unique_prefixed_and_unrecorded_is_explicit() -> None:
    row = _row(
        0,
        clinical_reasons=[
            "Allele 1: 7 primary alignments, below the per-allele depth gate (30).",
            "Allele 2: 9 primary alignments, below the per-allele depth gate (30).",
        ],
        reconstruction_flags=["ambiguous_reconstruction", "iupac_bases"],
    )
    assert case_reason_keys(row) == [
        "evaluator: ambiguous_reconstruction",
        "evaluator: iupac_bases",
        "gate: # primary alignments, below the per-allele depth gate (#)",
    ]
    assert case_reason_keys(_row(1, reasons_recorded=False)) == [UNRECORDED]
    assert case_reason_keys(_row(2)) == [UNRECORDED]


def test_expected_conditions_follow_the_config() -> None:
    low = _row(0, realized_min_allele_depth=GATE - 1, depth=GATE)
    assert expected_conditions(low, "dev", A) == ["depth"]
    assert expected_conditions(low, "dev", replace(A, depth_basis="design")) == []
    assert expected_conditions(_row(1), "stress", A) == ["split"]
    both = expected_conditions(low, "stress", A)
    assert both == ["split", "depth"]
    assert expected_conditions(_row(2, realized_min_allele_depth=None), "dev", A) == []
    none = replace(A, expected_inconclusive_splits=())
    assert expected_conditions(_row(3), "stress", none) == []


def _rows() -> list[dict[str, Any]]:
    iupac = ["ambiguous_reconstruction", "iupac_bases"]
    return [
        _row(0, reconstruction_flags=iupac),
        _row(1, reconstruction_flags=iupac, realized_min_allele_depth=GATE - 1, depth=GATE // 2),
        _row(2, reconstruction_flags=["missing_allele"], profile="ont_amplicon_r10"),
        _row(3, decision="NO_PATHOGENIC_VARIANT_DETECTED"),
        _row(4, decision="NO_CALL", reconstruction_flags=["insufficient_evidence"]),
    ]


def test_build_atlas_counts_reasons_split_and_strata() -> None:
    atlas = build_atlas(_rows(), "dev", A)
    assert (atlas["n_cases"], atlas["n_atlas"]) == (5, 3)
    reasons = {r["reason"]: r for r in atlas["reasons"]}
    top = reasons["evaluator: iupac_bases"]
    assert (top["k"], top["expected"], top["resolvable"]) == (2, 1, 1)
    assert top["share"] == pytest.approx(2 / 3) and top["rate"] == pytest.approx(2 / 5)
    assert atlas["reasons"][0]["reason"] == "evaluator: ambiguous_reconstruction"
    assert "evaluator: insufficient_evidence" not in reasons  # NO_CALL is not in decisions
    summary = {s["profile"]: s for s in atlas["split_summary"]}
    assert summary["all"] | {"conditions": None} == {
        "profile": "all",
        "cases": 5,
        "atlas": 3,
        "expected": 1,
        "resolvable": 2,
        "conditions": None,
    }
    assert summary["all"]["conditions"] == {"depth": 1}
    assert summary["hifi_amplicon"]["atlas"] == 2
    depth = atlas["by_stratum"]["depth"]
    cell = next(
        c
        for c in depth
        if c["profile"] == "hifi_amplicon" and c["value"] == GATE and c["reason"] == ANY_REASON
    )
    assert (cell["k"], cell["n"], cell["rate"]) == (1, 3, pytest.approx(1 / 3))
    assert set(atlas["by_stratum"]) == set(A.strata)
    wide = build_atlas(_rows(), "dev", replace(A, decisions=("INCONCLUSIVE", "NO_CALL")))
    assert wide["n_atlas"] == 4


def test_render_atlas_lists_the_split_and_top_reasons() -> None:
    cfg = replace(A, top_reasons=1)
    text = render_atlas(build_atlas(_rows(), "dev", cfg), cfg)
    assert "Reason atlas" in text and "expected" in text and "resolvable" in text
    assert "evaluator: ambiguous_reconstruction" in text
    assert "| R1 |" in text and "R2" not in text.split("### By", 1)[1]
    empty = render_atlas(build_atlas([], "dev", A), A)
    assert "no INCONCLUSIVE cases" in empty
