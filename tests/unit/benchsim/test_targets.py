"""benchsim.targets: absolute-target pass/fail evaluation (task 12e)."""

from dataclasses import replace
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, Target, TargetsConfig
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
        "scope": "pooled_and_profiles",
        "binding": True,
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


# Task 15o (decision rule v6): a "pooled"-scope target gates on the pooled row only.
# One profile at 2/3 INCONCLUSIVE (> 0.20) while the pooled rate is 2/12 (<= 0.20).
SCOPE_ROWS = [
    ("ont_genomic_targeted", "pathogenic", "INCONCLUSIVE"),
    ("ont_genomic_targeted", "normal", "INCONCLUSIVE"),
    ("ont_genomic_targeted", "normal", "NO_PATHOGENIC_VARIANT_DETECTED"),
    *[("ont_amplicon_r10", "pathogenic", "PATHOGENIC")] * 5,
    *[("ont_amplicon_r10", "normal", "NO_PATHOGENIC_VARIANT_DETECTED")] * 4,
]


def _scope_cfg() -> TargetsConfig:
    """Only the INCONCLUSIVE target, so the verdict isolates its scope."""
    return TargetsConfig(by_set={"standard": {"inconclusive_rate": Target("le", 0.20, "pooled")}})


def _row_of(result: dict[str, Any], grouping: str, metric: str) -> dict[str, Any]:
    return next(r for r in result["table"] if r["grouping"] == grouping and r["metric"] == metric)


def test_pooled_scope_target_ignores_a_failing_profile() -> None:
    result = evaluate_targets(_rows(SCOPE_ROWS), "standard", _scope_cfg(), ALPHA)
    assert result is not None and result["pass"] is True
    pooled = _row_of(result, "pooled", "inconclusive_rate")
    assert (pooled["k"], pooled["n"], pooled["pass"], pooled["binding"]) == (2, 12, True, True)
    profile = _row_of(result, "ont_genomic_targeted", "inconclusive_rate")
    # Still computed and shown, but informational: it does not fail the verdict.
    assert (profile["k"], profile["n"], profile["pass"]) == (2, 3, False)
    assert profile["binding"] is False and profile["scope"] == "pooled"


def test_pooled_scope_target_still_fails_on_the_pooled_rate() -> None:
    rows = _rows([*SCOPE_ROWS, *[("ont_amplicon_r10", "normal", "INCONCLUSIVE")] * 2])
    result = evaluate_targets(rows, "standard", _scope_cfg(), ALPHA)
    assert result is not None and result["pass"] is False
    assert _row_of(result, "pooled", "inconclusive_rate")["pass"] is False


def test_pooled_and_profiles_scope_still_gates_every_profile() -> None:
    by_set = {"standard": {"inconclusive_rate": Target("le", 0.20, "pooled_and_profiles")}}
    result = evaluate_targets(_rows(SCOPE_ROWS), "standard", TargetsConfig(by_set=by_set), ALPHA)
    assert result is not None and result["pass"] is False
    assert _row_of(result, "ont_genomic_targeted", "inconclusive_rate")["binding"] is True


def test_default_targets_gate_fp_and_pathogenic_per_profile_but_not_inconclusive() -> None:
    result = evaluate_targets(_rows(SCOPE_ROWS), "standard", CFG, ALPHA)
    assert result is not None
    binding = {(r["grouping"], r["metric"]): r["binding"] for r in result["table"]}
    assert binding[("ont_genomic_targeted", "inconclusive_rate")] is False
    assert binding[("ont_genomic_targeted", "pathogenic_rate")] is True
    assert binding[("ont_genomic_targeted", "false_positive_rate")] is True
    assert all(binding[("pooled", m)] for m in CFG.by_set["standard"])


def test_render_targets_marks_informational_rows() -> None:
    result = evaluate_targets(_rows(SCOPE_ROWS), "standard", _scope_cfg(), ALPHA)
    assert result is not None
    text = render_targets({"standard": result})
    assert "| standard | ont_genomic_targeted | inconclusive_rate |" in text
    assert "| False | info only |" in text and "| True | binding |" in text
    assert "Set verdict: standard=True" in text


def test_targets_text_states_the_scope_of_every_target() -> None:
    text = targets_text(CFG)
    assert "inconclusive_rate <= 0.2 (pooled only; per-profile rates for information)" in text
    assert "pathogenic_rate >= 0.8 (pooled and per profile)" in text
    assert "false_positive_rate <= 0 (pooled and per profile)" in text


def test_render_targets_omits_the_info_only_note_without_informational_rows() -> None:
    # Task 15o minor: the note explains `info only` rows, so it appears only with one.
    per_profile = TargetsConfig(
        by_set={"standard": {"inconclusive_rate": Target("le", 0.20, "pooled_and_profiles")}}
    )
    result = evaluate_targets(_rows(SCOPE_ROWS), "standard", per_profile, ALPHA)
    assert result is not None
    text = render_targets({"standard": result})
    assert "info only" not in text
    scoped = evaluate_targets(_rows(SCOPE_ROWS), "standard", _scope_cfg(), ALPHA)
    assert scoped is not None
    assert "Rows marked `info only`" in render_targets({"standard": scoped})


def test_targets_text_names_an_unknown_scope_instead_of_raising() -> None:
    # Task 15o minor: like the basis text, an unlisted scope is rendered verbatim (the
    # loader rejects it; a dataclass mutated after validation must not raise KeyError).
    cfg = _scope_cfg()
    cfg.by_set["standard"]["inconclusive_rate"] = Target("le", 0.20, "custom_scope")
    assert "inconclusive_rate <= 0.2 (custom_scope)" in targets_text(cfg)


def test_targets_text_states_the_v7_clean_inconclusive_ceiling() -> None:
    # Owner decision 2026-09-28 (task 15n, rule v7): pooled clean INCONCLUSIVE <= 0.15.
    text = targets_text(CFG)
    assert "`clean` requires false_positive_rate <= 0 (pooled and per profile), " in text
    assert "inconclusive_rate <= 0.15 (pooled only; per-profile rates for information)" in text
