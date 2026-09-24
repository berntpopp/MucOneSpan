"""Per-engine report tables for spec section 6 metrics 1-5 and the per-factor failure atlas.

Input is ``report.normalize_rows`` output (case rows). Metric 1 (per-allele
exact) is tabulated on ``report.allele_rows``; metric 2 (case exact) on the
case rows. Pooled per-allele and case-exact estimates carry 95%
cluster-bootstrap intervals over ``design_id`` (``sample``). Event recall and
precision (metric 3) use the conservative per-case bounds; ``dupC`` is
reported separately from other events. Metric 4 reuses
``evaluation.clinical_confusion.confusion`` per profile. The failure atlas
tabulates ``failure`` (case not exact, a false positive, or a critical false
negative) by every design factor in ``report.STRATA``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from muc_one_span.benchsim.report import (
    STRATA,
    allele_rows,
    fmt_value,
    render_tables,
    stratified_table,
)
from muc_one_span.benchsim.stats import clopper_pearson, cluster_bootstrap
from muc_one_span.evaluation.clinical_confusion import confusion

GROUPINGS = (("profile",), ("profile", "delta_class"), ("profile", "depth"))
N_BOOT = 2000


def _by_profile(rows: Sequence[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = {"all": list(rows)}
    for row in rows:
        groups.setdefault(str(row.get("profile")), []).append(row)
    return groups


def _pooled(rows: list[dict[str, Any]], metric: str, n_boot: int, seed: int) -> dict[str, Any]:
    if not rows:
        return {"n": 0, "clusters": 0, "point": None, "ci_low": None, "ci_high": None}
    point, low, high = cluster_bootstrap(rows, "sample", lambda r: r[metric], n_boot, seed)
    clusters = len({r["sample"] for r in rows})
    return {"n": len(rows), "clusters": clusters, "point": point, "ci_low": low, "ci_high": high}


def pooled_estimates(
    rows: Sequence[dict[str, Any]], n_boot: int = N_BOOT, seed: int = 0
) -> dict[str, dict[str, Any]]:
    """Per profile (and ``all``): per-allele and case-exact rates, cluster-bootstrap CIs."""
    out = {}
    for profile, group in _by_profile(rows).items():
        out[profile] = {
            "allele_exact": _pooled(allele_rows(group), "allele_exact", n_boot, seed),
            "case_exact": _pooled(group, "case_exact", n_boot, seed),
        }
    return out


def _ratio(k: int, n: int) -> dict[str, Any]:
    if n == 0:
        return {"k": k, "n": n, "rate": None, "ci_low": None, "ci_high": None}
    low, high = clopper_pearson(k, n)
    return {"k": k, "n": n, "rate": k / n, "ci_low": low, "ci_high": high}


def event_table(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Event recall (tp / truth events) and precision (tp / (tp + fp)) per profile x event class."""
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        cls = "dupC" if row.get("event") == "dupC" else "other"
        groups.setdefault((str(row.get("profile")), cls), []).append(row)
    table = []
    for (profile, cls), group in sorted(groups.items()):
        tp = sum(r["event_tp"] for r in group)
        fp = sum(r["event_fp"] for r in group)
        truth = sum(r["truth_events"] for r in group)
        table.append(
            {
                "profile": profile,
                "event_class": cls,
                "cases": len(group),
                "recall": _ratio(tp, truth),
                "precision": _ratio(tp, tp + fp),
            }
        )
    return table


def confusion_by_profile(rows: Sequence[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """``clinical_confusion.confusion`` per profile (and ``all``) on normalized decisions."""
    return {
        profile: confusion(
            [{"clinical": {"truth": r["truth"], "decision": r["decision"]}} for r in group]
        )
        for profile, group in _by_profile(rows).items()
    }


def build_tables(
    rows: Sequence[dict[str, Any]], n_boot: int = N_BOOT, seed: int = 0
) -> dict[str, Any]:
    """All per-engine tables: stratified CP tables, pooled bootstrap, events, confusion."""
    alleles = allele_rows(rows)
    normals = [r for r in rows if r["normal"] or r["benign"]]
    stratified: dict[str, list[dict[str, Any]]] = {}
    for by in GROUPINGS:
        label = " x ".join(by)
        stratified[f"allele_exact (metric 1, per allele) by {label}"] = stratified_table(
            alleles, by, "allele_exact"
        )
        stratified[f"case_exact (metric 2) by {label}"] = stratified_table(rows, by, "case_exact")
    for metric in ("critical_false_negative", "inconclusive", "no_call"):
        stratified[f"{metric} by profile"] = stratified_table(rows, ("profile",), metric)
    for metric in ("false_positive", "no_call"):
        stratified[f"{metric} on normal + benign truths by profile"] = stratified_table(
            normals, ("profile",), metric
        )
    for factor in STRATA:
        stratified[f"failure atlas: {factor}"] = stratified_table(rows, (factor,), "failure")
    return {
        "stratified": stratified,
        "pooled": pooled_estimates(rows, n_boot, seed),
        "events": event_table(rows),
        "confusion": confusion_by_profile(rows),
    }


def _ci(entry: dict[str, Any], rate: str = "rate") -> str:
    return f"{fmt_value(entry[rate])} [{fmt_value(entry['ci_low'])}, {fmt_value(entry['ci_high'])}]"


def render_engine_tables(tables: dict[str, Any]) -> str:
    """Markdown for ``build_tables`` output."""
    lines = [
        "### Pooled estimates (95% cluster bootstrap over design_id)",
        "",
        "| profile | per-allele exact | alleles | case exact | cases |",
        "|---|---|---|---|---|",
    ]
    for profile, est in tables["pooled"].items():
        al, ca = est["allele_exact"], est["case_exact"]
        lines.append(
            f"| {profile} | {_ci(al, 'point')} | {al['n']} | {_ci(ca, 'point')} | {ca['n']} |"
        )
    lines += [
        "",
        "### Event recall / precision (metric 3)",
        "",
        "| profile | events | cases | recall | precision |",
        "|---|---|---|---|---|",
    ]
    for row in tables["events"]:
        lines.append(
            f"| {row['profile']} | {row['event_class']} | {row['cases']} | "
            f"{_ci(row['recall'])} ({row['recall']['k']}/{row['recall']['n']}) | "
            f"{_ci(row['precision'])} |"
        )
    lines += ["", "### Clinical decision confusion (metric 4)", ""]
    for profile, conf in tables["confusion"].items():
        lines.append(
            f"- {profile}: matrix {conf['matrix']}; critical FN {conf['critical_false_negative']}; "
            f"FP normal {conf['false_positive_normal']}; FP benign {conf['false_positive_benign']}; "
            f"inconclusive rate {fmt_value(conf['inconclusive_rate'])}"
        )
    return "\n".join(lines) + "\n\n" + render_tables(tables["stratified"])
