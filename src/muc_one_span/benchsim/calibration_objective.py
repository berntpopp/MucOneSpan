"""Calibration objectives: the selection rule as data, per-point metrics and ranking.

An objective file (``schema_version`` 1) declares:

- ``constraints``: ``{metric: {"min": x, "max": y, "on": "value"|"ci_low"|"ci_high"}}``;
  at least one of ``min``/``max``. ``on`` (rate metrics only) chooses the judged
  quantity; without it the point value is judged.
- ``rank``: a non-empty lexicographic list of metric names; a ``-`` prefix ranks
  that metric descending (higher is better).
- ``reason_metrics`` (optional): ``{name: token}`` defines a rate metric, the share
  of cases whose clinical reasons or reconstruction flags contain ``token`` (for
  example ``{"smear_ambiguous_rate": "smear_ambiguous"}``).
- ``bench_sets`` (optional): restrict the rows to these benchmark sets.

Nothing is defaulted in code: a point is feasible when it violates no declared
constraint, and points rank feasible first, then by the ``rank`` terms, then by
their content hash (deterministic). Metrics come from `report.normalize_rows`
case rows. Rates carry ``1 - alpha`` cluster-bootstrap intervals over the case
(``sample``), reusing `stats.cluster_bootstrap` with the bench ``report``
settings (Task 8); counts carry no interval. A rate with an empty cohort is
``None``: it fails any constraint on it and ranks last.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from muc_one_span.benchsim.bench_config import ReportConfig
from muc_one_span.benchsim.calibration_grid import canonical_sha256, read_strict_json
from muc_one_span.benchsim.report import allele_rows
from muc_one_span.benchsim.stats import cluster_bootstrap

SCHEMA_VERSION = 1
DESCENDING = "-"
RATE, COUNT = "rate", "count"
JUDGED = ("value", "ci_low", "ci_high")
OBJECTIVE_FIELDS = frozenset(
    {"schema_version", "constraints", "rank", "reason_metrics", "bench_sets"}
)
CONSTRAINT_FIELDS = frozenset({"min", "max", "on"})
Row = dict[str, Any]


def _all(rows: Sequence[Row]) -> list[Row]:
    return list(rows)


def _negatives(rows: Sequence[Row]) -> list[Row]:
    return [r for r in rows if r["normal"] or r["benign"]]


def _pathogenic(rows: Sequence[Row]) -> list[Row]:
    return [r for r in rows if r["pathogenic"]]


# Rate metrics: (cohort, 0/1 row field). Counts: summed 0/1 row field (None: rows).
_RATES: dict[str, tuple[Callable[[Sequence[Row]], list[Row]], str]] = {
    "per_allele_exact": (allele_rows, "allele_exact"),
    "case_exact": (_all, "case_exact"),
    "inconclusive_rate": (_all, "inconclusive"),
    "no_call_rate": (_all, "no_call"),
    "failure_rate": (_all, "failure"),
    "false_positive_rate": (_negatives, "false_positive"),
    "critical_false_negative_rate": (_pathogenic, "critical_false_negative"),
}
_COUNTS: dict[str, str | None] = {
    "cases": None,
    "clinical_false_negative": "critical_false_negative",
    "false_positive": "false_positive",
    "not_completed": "failed",
}
METRICS = {**dict.fromkeys(_RATES, RATE), **dict.fromkeys(_COUNTS, COUNT)}


@dataclass(frozen=True)
class Constraint:
    """Bounds on one metric's judged quantity."""

    minimum: float | None
    maximum: float | None
    on: str


@dataclass(frozen=True)
class RankTerm:
    """One lexicographic ranking term."""

    metric: str
    descending: bool


@dataclass(frozen=True)
class Objective:
    """A validated selection rule and its canonical hash."""

    raw: dict[str, Any]
    sha256: str
    constraints: dict[str, Constraint]
    rank: tuple[RankTerm, ...]
    reason_metrics: dict[str, str]
    bench_sets: tuple[str, ...] | None

    def metric_names(self) -> list[str]:
        """Metrics the objective uses (constraints and rank), sorted."""
        return sorted({*self.constraints, *(t.metric for t in self.rank)})


def _number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _reasons(data: Any) -> dict[str, str]:
    if not isinstance(data, dict):
        raise ValueError("reason_metrics must be a JSON object")
    for name, token in data.items():
        if name in METRICS:
            raise ValueError(f"reason metric {name!r} shadows a built-in metric")
        if not isinstance(token, str) or not token:
            raise ValueError(f"reason metric {name!r} needs a non-empty reason token")
    return dict(data)


def _constraint(name: str, spec: Any, kind: str) -> Constraint:
    if not isinstance(spec, dict) or not spec.keys() <= CONSTRAINT_FIELDS:
        raise ValueError(f"constraint {name!r} must be an object with min, max and/or on")
    if not spec.keys() & {"min", "max"}:
        raise ValueError(f"constraint {name!r} needs min and/or max")
    for bound in ("min", "max"):
        if bound in spec and not _number(spec[bound]):
            raise ValueError(f"constraint {name!r}: {bound} must be a number")
    on = spec.get("on", JUDGED[0])
    if on not in JUDGED:
        raise ValueError(f"constraint {name!r}: on must be one of {', '.join(JUDGED)}")
    if on != JUDGED[0] and kind != RATE:
        raise ValueError(f"constraint {name!r}: CI bounds exist only for rate metrics")
    return Constraint(spec.get("min"), spec.get("max"), on)


def _check_metric(name: str, kinds: dict[str, str]) -> str:
    if name not in kinds:
        raise ValueError(f"unknown metric {name!r} (known: {', '.join(sorted(kinds))})")
    return kinds[name]


def _rank(data: Any, kinds: dict[str, str]) -> tuple[RankTerm, ...]:
    if not isinstance(data, list) or not data or not all(isinstance(t, str) for t in data):
        raise ValueError("rank must be a non-empty list of metric names")
    terms = [RankTerm(t.removeprefix(DESCENDING), t.startswith(DESCENDING)) for t in data]
    for term in terms:
        _check_metric(term.metric, kinds)
    names = [t.metric for t in terms]
    if len(set(names)) != len(names):
        raise ValueError("rank names a metric more than once")
    return tuple(terms)


def load_objective(path: Path) -> Objective:
    """Read and validate an objective file (every error names the offending entry)."""
    data = read_strict_json(path)
    if not isinstance(data, dict):
        raise ValueError("an objective must be a JSON object")
    if data.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"an objective requires schema_version={SCHEMA_VERSION}")
    unknown = data.keys() - OBJECTIVE_FIELDS
    if unknown:
        raise ValueError(f"unknown objective fields: {', '.join(sorted(unknown))}")
    reasons = _reasons(data.get("reason_metrics", {}))
    kinds = {**METRICS, **dict.fromkeys(reasons, RATE)}
    constraints_data = data.get("constraints", {})
    if not isinstance(constraints_data, dict):
        raise ValueError("constraints must be a JSON object")
    constraints = {
        name: _constraint(name, spec, _check_metric(name, kinds))
        for name, spec in constraints_data.items()
    }
    sets = data.get("bench_sets")
    if sets is not None and (
        not isinstance(sets, list) or not sets or not all(isinstance(s, str) for s in sets)
    ):
        raise ValueError("bench_sets must be a non-empty list of set names")
    return Objective(
        raw=data,
        sha256=canonical_sha256(data),
        constraints=constraints,
        rank=_rank(data.get("rank"), kinds),
        reason_metrics=reasons,
        bench_sets=tuple(sets) if sets is not None else None,
    )


def _rate(rows: list[Row], field: str, report: ReportConfig) -> dict[str, Any]:
    k = sum(int(r[field]) for r in rows)
    entry: dict[str, Any] = {"kind": RATE, "k": k, "n": len(rows)}
    if not rows:
        return entry | {"value": None, "ci_low": None, "ci_high": None}
    point, low, high = cluster_bootstrap(
        rows,
        "sample",
        lambda r: float(r[field]),
        n=report.bootstrap_replicates,
        seed=report.bootstrap_seed,
        alpha=report.alpha,
    )
    return entry | {"value": point, "ci_low": low, "ci_high": high}


def _has_reason(row: Row, token: str) -> int:
    texts = [*row.get("clinical_reasons", []), *row.get("reconstruction_flags", [])]
    return int(any(token in text for text in texts))


def point_metrics(
    rows: Sequence[Row], objective: Objective, report: ReportConfig
) -> dict[str, dict[str, Any]]:
    """Every built-in and declared reason metric over the (set-filtered) case rows."""
    if objective.bench_sets is not None:
        rows = [r for r in rows if r["bench_set"] in objective.bench_sets]
    out = {name: _rate(cohort(rows), field, report) for name, (cohort, field) in _RATES.items()}
    for name, count_field in _COUNTS.items():
        value = len(rows) if count_field is None else sum(int(r[count_field]) for r in rows)
        out[name] = {"kind": COUNT, "value": value, "ci_low": None, "ci_high": None}
    for name, token in objective.reason_metrics.items():
        flagged = [r | {name: _has_reason(r, token)} for r in rows]
        out[name] = _rate(flagged, name, report)
    return out


def judge(metrics: dict[str, dict[str, Any]], objective: Objective) -> list[str]:
    """Constraint violations of one point (empty: feasible), sorted by metric."""
    violations = []
    for name in sorted(objective.constraints):
        rule = objective.constraints[name]
        value = metrics[name].get(rule.on)
        if value is None:
            violations.append(f"{name}: no cases in its cohort")
        elif rule.minimum is not None and value < rule.minimum:
            violations.append(f"{name}: {rule.on} {value:g} < min {rule.minimum:g}")
        elif rule.maximum is not None and value > rule.maximum:
            violations.append(f"{name}: {rule.on} {value:g} > max {rule.maximum:g}")
    return violations


def _sort_key(point: dict[str, Any], objective: Objective) -> tuple[Any, ...]:
    terms: list[tuple[int, float]] = []
    for term in objective.rank:
        value = point["metrics"][term.metric]["value"]
        if value is None:
            terms.append((1, 0.0))
        else:
            terms.append((0, -value if term.descending else value))
    return (bool(point["violations"]), *terms, point["sha256"])


def rank_points(points: Sequence[dict[str, Any]], objective: Objective) -> list[dict[str, Any]]:
    """Points with ``violations``/``feasible``/``rank`` added, best first."""
    judged = [dict(p) | {"violations": judge(p["metrics"], objective)} for p in points]
    ordered = sorted(judged, key=lambda p: _sort_key(p, objective))
    return [
        p | {"feasible": not p["violations"], "rank": i} for i, p in enumerate(ordered, start=1)
    ]
