"""Absolute-target evaluation for the decision rule (task 12e, spec section 6).

Owner-approved floors/ceilings per bench set (`bench_config.TargetsConfig`), judged on
the candidate engine alone: unlike the relative rule (`report.decide`'s per-profile
McNemar / non-inferiority / critical-FN tests, decided against a baseline), a target is
a fixed number the candidate's own rate must clear, whatever the baseline does.

Input is `report.normalize_rows` output (case rows) of one engine and one bench set.
`evaluate_targets` reduces it to a pass/fail table, pooled over every profile present
plus one grouping per profile, using the same Clopper-Pearson interval as
`report.stratified_table`. `targets_text` renders the configured targets into
`report.rule_text`, so a pre-registered rule's SHA-256 changes whenever a target
(threshold, comparator, judging basis or bench-set membership) does.

Metric cohorts (`_cohort`, `TARGET_METRIC_NAMES`):

- ``pathogenic_rate``: PATHOGENIC decisions over truth ``pathogenic`` rows.
- ``inconclusive_rate``: INCONCLUSIVE decisions over every row (the assay's overall
  non-definitive rate, matching ``evaluation.clinical_confusion.confusion``).
- ``false_positive_rate``: PATHOGENIC decisions over truth ``normal`` or ``benign``
  rows (the same cohort as `report._fp_test`).

A grouping with no cases in its cohort fails its target (`_judge`) rather than being
silently skipped, so a bench set the candidate was never run on cannot pass by default.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from muc_one_span.benchsim.bench_config import Target, TargetsConfig
from muc_one_span.benchsim.stats import clopper_pearson

_COMPARATOR_SYMBOL = {"ge": ">=", "le": "<="}
_BASIS_TEXT = {
    "point": "the pooled/per-profile point estimate",
    "ci_bound": (
        "the Clopper-Pearson CI bound on the threshold's side (the lower bound for a "
        ">= target, the upper bound for a <= target)"
    ),
}


def _cohort(rows: Sequence[dict[str, Any]], metric: str) -> tuple[int, int]:
    """(successes, cases) of the metric's denominator cohort within `rows`."""
    if metric == "pathogenic_rate":
        cohort = [r for r in rows if r.get("pathogenic")]
        return sum(r["decision"] == "PATHOGENIC" for r in cohort), len(cohort)
    if metric == "inconclusive_rate":
        return sum(int(r["inconclusive"]) for r in rows), len(rows)
    if metric == "false_positive_rate":
        cohort = [r for r in rows if r.get("normal") or r.get("benign")]
        return sum(int(r["false_positive"]) for r in cohort), len(cohort)
    raise ValueError(f"unknown target metric {metric!r}")


def _judge(
    rows: Sequence[dict[str, Any]], metric: str, target: Target, basis: str, alpha: float
) -> dict[str, Any]:
    """Pass/fail for one (grouping, metric): cohort counts, rate, CI and judged value."""
    k, n = _cohort(rows, metric)
    base = {"metric": metric, "comparator": target.comparator, "threshold": target.threshold}
    if n == 0:
        return base | {
            "k": k,
            "n": n,
            "rate": None,
            "ci_low": None,
            "ci_high": None,
            "judged": None,
            "pass": False,
            "reason": "no cases in this cohort",
        }
    rate = k / n
    low, high = clopper_pearson(k, n, alpha=alpha)
    judged = rate if basis == "point" else (low if target.comparator == "ge" else high)
    passed = judged >= target.threshold if target.comparator == "ge" else judged <= target.threshold
    return base | {
        "k": k,
        "n": n,
        "rate": rate,
        "ci_low": low,
        "ci_high": high,
        "judged": judged,
        "pass": passed,
    }


def evaluate_targets(
    rows: Sequence[dict[str, Any]], set_name: str, config: TargetsConfig, alpha: float
) -> dict[str, Any] | None:
    """Pass/fail table for one engine's `rows` of one bench set (``None`` if untargeted).

    ``None`` when `set_name` is absent from ``config.by_set`` or maps to no metrics
    (for example `stress`): reported descriptively, never gated on a target. Otherwise
    one table row per (grouping, metric): ``"pooled"`` (every row) plus one grouping per
    profile present in `rows`.
    """
    metrics = config.by_set.get(set_name) or {}
    if not metrics:
        return None
    profiles = sorted({str(r.get("profile")) for r in rows})
    groups = {"pooled": list(rows)}
    groups.update({p: [r for r in rows if str(r.get("profile")) == p] for p in profiles})
    table = [
        {"grouping": grouping, **_judge(group, metric, target, config.basis, alpha)}
        for metric, target in sorted(metrics.items())
        for grouping, group in groups.items()
    ]
    return {
        "bench_set": set_name,
        "basis": config.basis,
        "profiles": profiles,
        "table": table,
        "pass": all(row["pass"] for row in table),
    }


def targets_text(targets: TargetsConfig) -> str:
    """Deterministic prose for `report.rule_text`: what must pass, and how it's judged."""
    basis = _BASIS_TEXT.get(targets.basis, targets.basis)
    sets_with_targets = {name: m for name, m in sorted(targets.by_set.items()) if m}
    if not sets_with_targets:
        return f"judged on {basis}, but no set currently has a configured target"
    clauses = [
        f"`{name}` requires "
        + ", ".join(
            f"{metric} {_COMPARATOR_SYMBOL[t.comparator]} {t.threshold:g}"
            for metric, t in sorted(metrics.items())
        )
        for name, metrics in sets_with_targets.items()
    ]
    return f"judged on {basis}, " + " and ".join(clauses)


def render_targets(by_set: dict[str, dict[str, Any] | None]) -> str:
    """Markdown pass/fail table for `report.decide`'s ``result["targets"]``."""
    entries = [(name, res) for name, res in by_set.items() if res is not None]
    if not entries:
        return ""
    lines = [
        "### Absolute targets (task 12e)",
        "",
        "| set | grouping | metric | target | k/n | rate | judged | pass |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for set_name, res in entries:
        for row in res["table"]:
            symbol = _COMPARATOR_SYMBOL[row["comparator"]]
            rate = "n/a" if row["rate"] is None else f"{row['rate']:.4g}"
            judged = "n/a" if row["judged"] is None else f"{row['judged']:.4g}"
            lines.append(
                f"| {set_name} | {row['grouping']} | {row['metric']} | "
                f"{symbol} {row['threshold']:g} | {row['k']}/{row['n']} | {rate} | "
                f"{judged} | {row['pass']} |"
            )
    verdict = ", ".join(f"{name}={res['pass']}" for name, res in entries)
    return "\n".join([*lines, "", f"Set verdict: {verdict}", ""])
