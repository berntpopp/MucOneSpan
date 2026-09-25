"""benchsim.report.decide: the task 12e absolute targets and the task C1 FP ruling.

Split out of ``test_report.py`` (file-size gate): this covers Part 2 of the decision
rule (``targets.by_set``) and, since task C1, the false-positive criterion that now
lives there instead of in Part 1's relative comparison.
"""

from dataclasses import replace
from typing import Any

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, BenchConfig, TargetsConfig
from muc_one_span.benchsim.report import decide, render_markdown

# Bypasses the task 12e absolute targets (part 2) so a fixture that exercises only the
# relative rule (part 1) does not also need pathogenic/decision/inconclusive fields.
NO_TARGETS = BenchConfig(targets=TargetsConfig(by_set={}))


def _target_row(
    sample: str, profile: str, bench_set: str, exact: int, truth: str, decision: str
) -> dict[str, Any]:
    """A normalized case row exercising both the relative rule and the absolute targets."""
    fp = int(decision == "PATHOGENIC" and truth in ("normal", "benign"))
    false_negative = int(truth == "pathogenic" and decision != "PATHOGENIC")
    return {
        "sample": sample,
        "profile": profile,
        "bench_set": bench_set,
        "alleles": [{"allele": "h1", "allele_exact": exact}],
        "normal": truth == "normal",
        "benign": truth == "benign",
        "pathogenic": truth == "pathogenic",
        "decision": decision,
        "inconclusive": int(decision == "INCONCLUSIVE"),
        "false_positive": fp,
        "critical_false_negative": false_negative,
    }


def _standard_fixture(
    n_path: int, n_normal: int
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """(baseline, candidate) ``standard``-set rows: candidate superior (exact=1, fewer
    critical FN) and meeting the default ``standard`` targets (pathogenic_rate >= 0.80,
    inconclusive_rate <= 0.20, false_positive_rate <= 0)."""
    base, cand = [], []
    for i in range(n_path):
        sample = f"p{i}"
        base.append(_target_row(sample, "ont_amplicon_r10", "standard", 0, "pathogenic", "NO_CALL"))
        cand.append(
            _target_row(sample, "ont_amplicon_r10", "standard", 1, "pathogenic", "PATHOGENIC")
        )
    for i in range(n_normal):
        sample, decision = f"n{i}", "NO_PATHOGENIC_VARIANT_DETECTED"
        base.append(_target_row(sample, "ont_amplicon_r10", "standard", 0, "normal", decision))
        cand.append(_target_row(sample, "ont_amplicon_r10", "standard", 1, "normal", decision))
    return base, cand


def _clean_fixture(n_path: int, n_normal: int) -> list[dict[str, Any]]:
    """Candidate-only ``clean``-set rows meeting its targets (>= 0.90 / <= 0.10 / <= 0)."""
    cand = [
        _target_row(f"cp{i}", "ont_amplicon_r10", "clean", 1, "pathogenic", "PATHOGENIC")
        for i in range(n_path)
    ]
    cand += [
        _target_row(
            f"cn{i}", "ont_amplicon_r10", "clean", 1, "normal", "NO_PATHOGENIC_VARIANT_DETECTED"
        )
        for i in range(n_normal)
    ]
    return cand


def test_decide_adopts_only_when_headline_and_clean_targets_also_pass() -> None:
    base, cand_standard = _standard_fixture(2700, 1000)
    clean_ok = _clean_fixture(9, 1)
    result = decide({"ladder": base, "hybrid": cand_standard + clean_ok}, "ladder", "hybrid")
    assert set(result["targets"]) == {"standard", "clean"}
    assert result["targets"]["standard"]["pass"] is True
    assert result["targets"]["clean"]["pass"] is True
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True  # the relative rule too
    assert result["adopt"] is True

    # Targets are not restricted to the headline set: a `clean`-only regression (one
    # false positive) fails adoption even though the relative rule and the `standard`
    # targets are untouched.
    clean_bad = _clean_fixture(9, 1)
    clean_bad[-1] = clean_bad[-1] | {"decision": "PATHOGENIC", "false_positive": 1}
    result2 = decide({"ladder": base, "hybrid": cand_standard + clean_bad}, "ladder", "hybrid")
    assert result2["targets"]["standard"]["pass"] is True
    assert result2["targets"]["clean"]["pass"] is False
    assert result2["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result2["adopt"] is False

    # Review I4: a targeted set absent from this root is "not present", never FAIL,
    # and still blocks adoption (fail closed).
    result3 = decide({"ladder": base, "hybrid": cand_standard}, "ladder", "hybrid")
    assert result3["targets"]["clean"]["present"] is False
    assert result3["targets"]["clean"]["pass"] is None
    assert result3["targets_not_present"] == ["clean"]
    assert result3["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result3["adopt"] is False
    assert "not present" in render_markdown(result3)


_STANDARD_ONLY = BenchConfig(
    targets=TargetsConfig(by_set={"standard": DEFAULT_BENCH_CONFIG.targets.by_set["standard"]})
)


def test_decide_adopts_at_planned_test_size_with_zero_fp() -> None:
    # Task C1, the exact final-review scenario: 280 normals/profile, 0 FP in both
    # engines. Before the owner's ruling this could never ADOPT (the Newcombe
    # non-inferiority upper bound, 0.00957, exceeded the 0.005 margin even at 0
    # observed FP). FP is now judged only by the absolute `false_positive_rate <= 0`
    # target, which 0 FP clears, so this now adopts.
    base, cand = _standard_fixture(n_path=30, n_normal=280)
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", _STANDARD_ONLY)
    fp = result["profiles"]["ont_amplicon_r10"]["false_positive"]
    assert fp["fp"] == 0 and fp["fp_baseline"] == 0 and fp["n"] == 280
    assert result["targets"]["standard"]["pass"] is True
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result["adopt"] is True


def test_decide_candidate_fp_on_a_target_set_still_fails_adoption() -> None:
    # The relative rule no longer looks at FP at all (task C1); a single candidate
    # false positive on a targeted, headline-set profile still fails adoption
    # through part 2's absolute `false_positive_rate` target.
    base, cand = _standard_fixture(n_path=30, n_normal=280)
    cand[-1] = cand[-1] | {"decision": "PATHOGENIC", "false_positive": 1}
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", _STANDARD_ONLY)
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True  # relative rule ignores FP
    assert result["targets"]["standard"]["pass"] is False  # absolute target: FP > 0
    assert result["adopt"] is False


def test_decide_targets_use_point_estimate_by_default_and_ci_bound_when_configured() -> None:
    base, cand = _standard_fixture(4, 1)  # 4/4 pathogenic = 1.0, point estimate passes easily
    cand = cand + _clean_fixture(9, 1)
    point = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    assert point["targets"] == {}  # NO_TARGETS configures no sets at all

    with_targets = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    assert with_targets["targets"]["standard"]["basis"] == "point"
    ci_cfg = replace(
        DEFAULT_BENCH_CONFIG, targets=replace(DEFAULT_BENCH_CONFIG.targets, basis="ci_bound")
    )
    with_ci = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", ci_cfg)
    assert with_ci["targets"]["standard"]["basis"] == "ci_bound"
    pooled_path = next(
        r
        for r in with_ci["targets"]["standard"]["table"]
        if r["grouping"] == "pooled" and r["metric"] == "pathogenic_rate"
    )
    assert pooled_path["judged"] == pooled_path["ci_low"]


def test_render_markdown_shows_the_absolute_targets_table() -> None:
    base, cand = _standard_fixture(2700, 1000)
    cand = cand + _clean_fixture(9, 1)
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    text = render_markdown(result)
    assert "Part 2: absolute targets" in text and "Absolute targets (task 12e)" in text
    assert "| standard | pooled | pathogenic_rate |" in text
    assert "Set verdict: standard=True, clean=True" in text
    assert "**ADOPT**" in text
