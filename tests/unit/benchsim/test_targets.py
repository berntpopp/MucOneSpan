"""benchsim.targets: absolute-target pass/fail evaluation (task 12e)."""

from dataclasses import replace
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, TargetsConfig
from muc_one_span.benchsim.targets import _cohort, evaluate_targets, render_targets, targets_text

CFG = DEFAULT_BENCH_CONFIG.targets
ALPHA = DEFAULT_BENCH_CONFIG.report.alpha


def _row(profile: str, truth: str, decision: str) -> dict[str, Any]:
    return {
        "sample": f"{profile}-{truth}-{decision}",
        "profile": profile,
        "pathogenic": truth == "pathogenic",
        "normal": truth == "normal",
        "benign": truth == "benign",
        "decision": decision,
        "inconclusive": int(decision == "INCONCLUSIVE"),
        "false_positive": int(decision == "PATHOGENIC" and truth in ("normal", "benign")),
    }


def _rows(spec: list[tuple[str, str, str]]) -> list[dict[str, Any]]:
    return [_row(*s) for s in spec]


PASSING_STANDARD = [
    ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
    ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
    ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
    ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
    ("ont_amplicon_r10", "normal", "NO_PATHOGENIC_VARIANT_DETECTED"),
    ("hifi_amplicon", "pathogenic", "PATHOGENIC"),
    ("hifi_amplicon", "normal", "NO_PATHOGENIC_VARIANT_DETECTED"),
]


def test_evaluate_targets_returns_none_for_an_untargeted_set() -> None:
    assert evaluate_targets(_rows(PASSING_STANDARD), "stress", CFG, ALPHA) is None


def test_evaluate_targets_pools_and_splits_by_profile() -> None:
    result = evaluate_targets(_rows(PASSING_STANDARD), "standard", CFG, ALPHA)
    assert result is not None
    assert result["bench_set"] == "standard" and result["basis"] == "point"
    assert result["profiles"] == ["hifi_amplicon", "ont_amplicon_r10"]
    groupings = {(row["grouping"], row["metric"]) for row in result["table"]}
    assert ("pooled", "pathogenic_rate") in groupings
    assert ("ont_amplicon_r10", "pathogenic_rate") in groupings
    assert ("hifi_amplicon", "inconclusive_rate") in groupings
    pooled_path = next(
        r for r in result["table"] if r["grouping"] == "pooled" and r["metric"] == "pathogenic_rate"
    )
    assert (pooled_path["k"], pooled_path["n"]) == (5, 5) and pooled_path["rate"] == 1.0
    assert pooled_path["pass"] is True and result["pass"] is True


def test_evaluate_targets_fails_when_a_metric_misses_its_threshold() -> None:
    rows = _rows(
        [
            ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
            ("ont_amplicon_r10", "pathogenic", "INCONCLUSIVE"),  # 1/2 = 0.5 < 0.80
            ("ont_amplicon_r10", "normal", "NO_PATHOGENIC_VARIANT_DETECTED"),
        ]
    )
    result = evaluate_targets(rows, "standard", CFG, ALPHA)
    assert result is not None and result["pass"] is False
    pooled_path = next(
        r for r in result["table"] if r["grouping"] == "pooled" and r["metric"] == "pathogenic_rate"
    )
    assert pooled_path["pass"] is False and pooled_path["rate"] == 0.5


def test_evaluate_targets_false_positive_must_be_exactly_zero() -> None:
    rows = _rows(
        [
            ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
            ("ont_amplicon_r10", "normal", "PATHOGENIC"),  # false positive
        ]
    )
    result = evaluate_targets(rows, "standard", CFG, ALPHA)
    assert result is not None
    fp_pooled = next(
        r
        for r in result["table"]
        if r["grouping"] == "pooled" and r["metric"] == "false_positive_rate"
    )
    assert (fp_pooled["k"], fp_pooled["n"], fp_pooled["pass"]) == (1, 1, False)


def test_evaluate_targets_empty_cohort_fails_with_a_reason() -> None:
    # No normal or benign truths at all: the false-positive cohort is empty.
    rows = _rows([("ont_amplicon_r10", "pathogenic", "PATHOGENIC")])
    result = evaluate_targets(rows, "standard", CFG, ALPHA)
    assert result is not None and result["pass"] is False
    fp_pooled = next(
        r
        for r in result["table"]
        if r["grouping"] == "pooled" and r["metric"] == "false_positive_rate"
    )
    assert fp_pooled == {
        "grouping": "pooled",
        "metric": "false_positive_rate",
        "comparator": "le",
        "threshold": 0.0,
        "k": 0,
        "n": 0,
        "rate": None,
        "ci_low": None,
        "ci_high": None,
        "judged": None,
        "pass": False,
        "reason": "no cases in this cohort",
    }


def test_evaluate_targets_ci_bound_is_more_conservative_than_point() -> None:
    # 4/5 = 0.80 exactly meets the point-estimate target but not the CI lower bound.
    rows = _rows(
        [
            ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
            ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
            ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
            ("ont_amplicon_r10", "pathogenic", "PATHOGENIC"),
            ("ont_amplicon_r10", "pathogenic", "INCONCLUSIVE"),
        ]
    )
    point = evaluate_targets(rows, "standard", CFG, ALPHA)
    assert point is not None
    point_path = next(
        r for r in point["table"] if r["grouping"] == "pooled" and r["metric"] == "pathogenic_rate"
    )
    assert point_path["judged"] == 0.8 and point_path["pass"] is True

    bound_cfg = replace(CFG, basis="ci_bound")
    bound = evaluate_targets(rows, "standard", bound_cfg, ALPHA)
    assert bound is not None
    bound_path = next(
        r for r in bound["table"] if r["grouping"] == "pooled" and r["metric"] == "pathogenic_rate"
    )
    assert bound_path["judged"] == bound_path["ci_low"] < 0.8
    assert bound_path["pass"] is False and bound["pass"] is False


def test_cohort_rejects_an_unknown_metric() -> None:
    with pytest.raises(ValueError, match="unknown target metric"):
        _cohort(_rows(PASSING_STANDARD), "bogus_rate")


def test_targets_text_lists_every_configured_set_and_metric() -> None:
    text = targets_text(CFG)
    assert "point estimate" in text
    assert "`clean` requires" in text and "`standard` requires" in text
    assert "pathogenic_rate >= 0.9" in text and "inconclusive_rate <= 0.2" in text
    assert "false_positive_rate <= 0" in text


def test_targets_text_reports_the_ci_bound_basis() -> None:
    text = targets_text(replace(CFG, basis="ci_bound"))
    assert "Clopper-Pearson" in text


def test_targets_text_handles_no_configured_targets() -> None:
    text = targets_text(TargetsConfig(by_set={}))
    assert "no set" in text


def test_render_targets_builds_a_markdown_pass_fail_table() -> None:
    result = evaluate_targets(_rows(PASSING_STANDARD), "standard", CFG, ALPHA)
    assert result is not None
    text = render_targets({"standard": result})
    assert "Absolute targets" in text
    assert "| standard | pooled | pathogenic_rate |" in text
    assert "Set verdict: standard=True" in text


def test_render_targets_skips_untargeted_sets() -> None:
    assert render_targets({"stress": None}) == ""


def test_targets_by_set_marks_an_absent_set_not_present() -> None:
    """Review I4: a targeted set with no cases in this root is not reported as FAIL."""
    from muc_one_span.benchsim.targets import targets_by_set

    rows = [r | {"bench_set": "standard"} for r in _rows(PASSING_STANDARD)]
    result = targets_by_set(rows, CFG, ALPHA)
    assert set(result) == set(CFG.by_set)
    assert result["standard"] is not None and result["standard"]["pass"] is True
    assert result["standard"]["present"] is True
    absent = result["clean"]
    assert absent is not None and absent["present"] is False
    assert absent["pass"] is None and absent["table"] == []
    text = render_targets(result)
    assert "| clean | - | - | - | 0 cases | n/a | n/a | not present |" in text
    assert "Set verdict: standard=True, clean=not present" in text
    assert "FAIL" not in text and "clean=False" not in text
