"""Calibration objectives: strict loading, per-point metrics, constraints and ranking."""

import json
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG
from muc_one_span.benchsim.calibration_objective import (
    METRICS,
    judge,
    load_objective,
    point_metrics,
    rank_points,
)

REPORT = DEFAULT_BENCH_CONFIG.report
HEADLINE = DEFAULT_BENCH_CONFIG.sets.headline
OBJECTIVE = {
    "schema_version": 1,
    "constraints": {"clinical_false_negative": {"max": 0}, "smear_ambiguous_rate": {"max": 0.5}},
    "rank": ["-per_allele_exact", "inconclusive_rate"],
    "reason_metrics": {"smear_ambiguous_rate": "smear_ambiguous"},
}


def _write(tmp_path: Path, data: object) -> Path:
    path = tmp_path / "objective.json"
    path.write_text(json.dumps(data))
    return path


def _row(
    name: str,
    truth: str,
    decision: str,
    exact: int,
    reasons: tuple[str, ...] = (),
    bench_set: str = HEADLINE,
) -> dict[str, Any]:
    cfn = int(truth == "pathogenic" and decision == "NO_PATHOGENIC_VARIANT_DETECTED")
    fp = int(truth != "pathogenic" and decision == "PATHOGENIC")
    return {
        "sample": name,
        "bench_set": bench_set,
        "failed": False,
        "normal": truth == "normal",
        "benign": truth == "benign",
        "pathogenic": truth == "pathogenic",
        "alleles": [{"allele": a, "allele_exact": exact} for a in ("h1", "h2")],
        "case_exact": exact,
        "inconclusive": int(decision == "INCONCLUSIVE"),
        "no_call": int(decision == "NO_CALL"),
        "false_positive": fp,
        "critical_false_negative": cfn,
        "failure": int(not exact or fp or cfn),
        "clinical_reasons": list(reasons),
        "reconstruction_flags": [],
    }


def test_objective_round_trip_and_hash(tmp_path: Path) -> None:
    objective = load_objective(_write(tmp_path, OBJECTIVE))
    assert objective.raw == OBJECTIVE
    assert [(t.metric, t.descending) for t in objective.rank] == [
        ("per_allele_exact", True),
        ("inconclusive_rate", False),
    ]
    assert objective.sha256 == load_objective(_write(tmp_path, OBJECTIVE)).sha256
    assert objective.metric_names() == [
        "clinical_false_negative",
        "inconclusive_rate",
        "per_allele_exact",
        "smear_ambiguous_rate",
    ]


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"schema_version": 2}, "schema_version"),
        ({"extra": 1}, "unknown objective fields"),
        ({"rank": []}, "rank must be a non-empty list"),
        ({"rank": ["-per_allele_exact", "per_allele_exact"]}, "more than once"),
        ({"rank": ["nope"]}, "unknown metric 'nope'"),
        ({"constraints": {"nope": {"max": 0}}}, "unknown metric 'nope'"),
        ({"constraints": {"cases": {}}}, "min and/or max"),
        ({"constraints": {"cases": {"max": "0"}}}, "must be a number"),
        ({"constraints": {"cases": {"max": 1, "on": "ci_high"}}}, "rate metrics"),
        ({"constraints": {"case_exact": {"min": 1, "on": "median"}}}, "on must be one of"),
        ({"reason_metrics": {"cases": "x"}}, "shadows a built-in"),
        ({"reason_metrics": {"x_rate": ""}}, "non-empty reason token"),
        ({"bench_sets": []}, "bench_sets must be a non-empty list"),
        ({"reason_metrics": ["x"]}, "reason_metrics must be a JSON object"),
        ({"constraints": {"cases": 0}}, "must be an object"),
        ({"constraints": []}, "constraints must be a JSON object"),
    ],
)
def test_invalid_objectives_fail(tmp_path: Path, change: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        load_objective(_write(tmp_path, {**OBJECTIVE, **change}))


def test_objective_must_be_an_object(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="JSON object"):
        load_objective(_write(tmp_path, [OBJECTIVE]))


def test_objective_without_constraints_is_allowed(tmp_path: Path) -> None:
    data = {"schema_version": 1, "rank": ["case_exact"]}
    assert load_objective(_write(tmp_path, data)).constraints == {}


def test_point_metrics_rates_counts_and_reason_metrics(tmp_path: Path) -> None:
    objective = load_objective(_write(tmp_path, OBJECTIVE))
    rows = [
        _row("c1", "pathogenic", "PATHOGENIC", 1),
        _row("c2", "normal", "INCONCLUSIVE", 0, ("A: smear_ambiguous peak",)),
        _row("c3", "pathogenic", "NO_PATHOGENIC_VARIANT_DETECTED", 0),
        _row("c4", "normal", "PATHOGENIC", 1),
    ]
    metrics = point_metrics(rows, objective, REPORT)
    assert set(metrics) == set(METRICS) | {"smear_ambiguous_rate"}
    assert metrics["cases"]["value"] == len(rows)
    assert metrics["clinical_false_negative"]["value"] == 1
    assert metrics["false_positive"]["value"] == 1
    assert metrics["per_allele_exact"]["value"] == pytest.approx(0.5)
    assert metrics["per_allele_exact"]["n"] == 2 * len(rows)
    assert metrics["inconclusive_rate"]["value"] == pytest.approx(1 / len(rows))
    assert metrics["smear_ambiguous_rate"]["value"] == pytest.approx(1 / len(rows))
    assert metrics["false_positive_rate"]["n"] == 2
    assert metrics["critical_false_negative_rate"]["value"] == pytest.approx(0.5)
    low, high = metrics["case_exact"]["ci_low"], metrics["case_exact"]["ci_high"]
    assert low <= metrics["case_exact"]["value"] <= high
    assert metrics["cases"]["ci_low"] is None


def test_bench_sets_filter_rows(tmp_path: Path) -> None:
    objective = load_objective(_write(tmp_path, {**OBJECTIVE, "bench_sets": [HEADLINE]}))
    rows = [_row("c1", "normal", "INCONCLUSIVE", 1), _row("c2", "normal", "NO_CALL", 0, (), "x")]
    assert point_metrics(rows, objective, REPORT)["cases"]["value"] == 1


def test_empty_cohort_rate_is_none_and_fails_its_constraint(tmp_path: Path) -> None:
    data = {**OBJECTIVE, "constraints": {"false_positive_rate": {"max": 0}}}
    objective = load_objective(_write(tmp_path, data))
    metrics = point_metrics([_row("c1", "pathogenic", "PATHOGENIC", 1)], objective, REPORT)
    assert metrics["false_positive_rate"]["value"] is None
    assert judge(metrics, objective) == ["false_positive_rate: no cases in its cohort"]


def test_judge_min_max_and_ci_bounds(tmp_path: Path) -> None:
    data = {
        **OBJECTIVE,
        "constraints": {
            "clinical_false_negative": {"max": 0},
            "case_exact": {"min": 0.5, "on": "ci_low"},
        },
    }
    objective = load_objective(_write(tmp_path, data))
    metrics = {
        "clinical_false_negative": {"value": 1},
        "case_exact": {"value": 0.9, "ci_low": 0.4, "ci_high": 1.0},
    }
    assert judge(metrics, objective) == [
        "case_exact: ci_low 0.4 < min 0.5",
        "clinical_false_negative: value 1 > max 0",
    ]


def test_rank_is_feasible_first_then_lexicographic_then_hash(tmp_path: Path) -> None:
    objective = load_objective(_write(tmp_path, OBJECTIVE))

    def point(sha: str, exact: float, inconclusive: float, cfn: int) -> dict[str, Any]:
        return {
            "sha256": sha,
            "metrics": {
                "per_allele_exact": {"value": exact},
                "inconclusive_rate": {"value": inconclusive},
                "clinical_false_negative": {"value": cfn},
                "smear_ambiguous_rate": {"value": 0.0},
            },
        }

    ranked = rank_points(
        [
            point("e", 0.99, 0.0, 1),  # best accuracy but infeasible
            point("d", 0.9, 0.2, 0),
            point("c", 0.9, 0.1, 0),
            point("b", 0.8, 0.0, 0),
            point("a", 0.9, 0.1, 0),  # ties with c: hash order
        ],
        objective,
    )
    assert [p["sha256"] for p in ranked] == ["a", "c", "d", "b", "e"]
    assert [p["rank"] for p in ranked] == [1, 2, 3, 4, 5]
    assert ranked[-1]["feasible"] is False
    assert ranked[-1]["violations"] == ["clinical_false_negative: value 1 > max 0"]


LENGTHS_RATES = {"allele_count_exact": (lambda rows: list(rows), "allele_count_exact")}
LENGTHS_COUNTS: dict[str, str | None] = {"cases": None, "false_alleles": "false_alleles"}
LENGTHS_METRICS = {"allele_count_exact": "rate", "cases": "count", "false_alleles": "count"}


def test_load_objective_accepts_an_injected_metric_registry(tmp_path: Path) -> None:
    data = {"schema_version": 1, "rank": ["allele_count_exact"]}
    objective = load_objective(_write(tmp_path, data), known_metrics=LENGTHS_METRICS)
    assert objective.rank[0].metric == "allele_count_exact"
    with pytest.raises(ValueError, match="unknown metric 'per_allele_exact'"):
        load_objective(
            _write(tmp_path, {"schema_version": 1, "rank": ["per_allele_exact"]}),
            known_metrics=LENGTHS_METRICS,
        )


def test_reason_metric_collision_is_checked_against_the_injected_registry(
    tmp_path: Path,
) -> None:
    # "allele_count_exact" is not a full-pipeline metric (METRICS), but it IS a
    # lengths-stage built-in (LENGTHS_METRICS): a reason_metrics entry reusing that
    # name must be rejected for the lengths registry, or it silently overwrites the
    # built-in rate in point_metrics.
    data = {
        "schema_version": 1,
        "rank": ["allele_count_exact"],
        "reason_metrics": {"allele_count_exact": "smear_ambiguous"},
    }
    with pytest.raises(ValueError, match="shadows a built-in"):
        load_objective(_write(tmp_path, data), known_metrics=LENGTHS_METRICS)
    # Not a full-pipeline built-in, so it is fine under the default (full) registry.
    load_objective(_write(tmp_path, {**data, "rank": ["cases"]}))


def test_reason_metric_collision_is_still_checked_for_the_full_stage(tmp_path: Path) -> None:
    data = {
        "schema_version": 1,
        "rank": ["case_exact"],
        "reason_metrics": {"case_exact": "smear_ambiguous"},
    }
    with pytest.raises(ValueError, match="shadows a built-in"):
        load_objective(_write(tmp_path, data))


def test_point_metrics_accepts_injected_rates_and_counts(tmp_path: Path) -> None:
    data = {"schema_version": 1, "rank": ["allele_count_exact"]}
    objective = load_objective(_write(tmp_path, data), known_metrics=LENGTHS_METRICS)
    rows = [
        {"sample": "c1", "bench_set": HEADLINE, "allele_count_exact": 1, "false_alleles": 0},
        {"sample": "c2", "bench_set": HEADLINE, "allele_count_exact": 0, "false_alleles": 2},
    ]
    metrics = point_metrics(rows, objective, REPORT, rates=LENGTHS_RATES, counts=LENGTHS_COUNTS)
    assert set(metrics) == set(LENGTHS_METRICS)
    assert metrics["allele_count_exact"]["value"] == pytest.approx(0.5)
    assert metrics["cases"]["value"] == 2
    assert metrics["false_alleles"]["value"] == 2


def test_missing_rank_value_sorts_last(tmp_path: Path) -> None:
    objective = load_objective(_write(tmp_path, {"schema_version": 1, "rank": ["-case_exact"]}))
    ranked = rank_points(
        [
            {"sha256": "a", "metrics": {"case_exact": {"value": None}}},
            {"sha256": "b", "metrics": {"case_exact": {"value": 0.1}}},
        ],
        objective,
    )
    assert [p["sha256"] for p in ranked] == ["b", "a"]
