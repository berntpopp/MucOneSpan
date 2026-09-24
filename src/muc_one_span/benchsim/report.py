"""Stratified benchmark report, paired engine comparison and the pre-registered decision rule.

Consumes ``scripts/evaluate.py`` reports (one per engine; per-sample rows under
``report["samples"]``) and generated ``case.json`` dicts keyed by
``design_id``. ``normalize_rows`` reduces each sample row to the binary
endpoints of spec section 6 (one row per case, with one entry per truth allele
under ``alleles``); ``allele_rows`` flattens those to the per-allele unit of
metric 1; ``decide`` applies the decision rule (``RULE_TEXT``) per profile.

Key mapping from ``evaluation.scoring.evaluate_sample`` rows:

- Truth alleles: the truth haplotype names in ``alternatives[0]`` (assigned
  ``pairs[*]["truth"]`` plus ``missing_truth``); their number must equal
  ``truth_haplotypes``.
- ``alleles[*]["allele_exact"]`` (metric 1): ``pairs[*]["sequence_exact"]`` of
  the least favourable optimal assignment (fewest
  ``metrics["independent_sequence_exact"]``). Exact pairs beyond that
  alternative's independent count (unproven duplicate observations) are
  demoted in descending truth-name order, so the per-allele total always
  equals ``metrics["independent_sequence_exact"]["min"]``.
- ``case_exact`` (metric 2): ``metrics["all_sequences_exact"]["min"]``.
- ``event_tp`` / ``event_fp`` / ``truth_events`` (metric 3):
  ``metrics["event_tp"]["min"]``, ``metrics["event_fp"]["max"]`` and
  ``truth_events`` (the conservative bounds ``scoring.aggregate`` uses).
- ``truth`` / ``decision`` (metrics 4-5): ``row["clinical"]``.

Runs with status ``execution_failed``, ``not_attempted`` or
``invalid_artifacts`` are kept: every truth allele scores 0, events score 0
and the decision is ``NO_CALL`` (so a pathogenic truth is a critical false
negative), whatever stale files the result directory holds. Rows without
valid truth (``invalid_truth``) raise: their truth alleles are unknown and are
never guessed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from muc_one_span.benchsim.preregistration import (
    first_evaluation,
    mark_first_evaluation,
    preregister,
    require_preregistered,
    rule_sha256,
)
from muc_one_span.benchsim.stats import clopper_pearson, holm, mcnemar_exact, noninferior

__all__ = [
    "RULE_TEXT",
    "allele_rows",
    "decide",
    "first_evaluation",
    "mark_first_evaluation",
    "normalize_rows",
    "paired",
    "preregister",
    "render_markdown",
    "render_tables",
    "require_preregistered",
    "rule_sha256",
    "stratified_table",
]

FAILED_STATUSES = frozenset({"execution_failed", "not_attempted", "invalid_artifacts"})
FALSE_NEGATIVE_DECISIONS = frozenset({"NO_PATHOGENIC_VARIANT_DETECTED", "NO_CALL"})
NI_MARGIN = 0.005
ALPHA = 0.05
HOLM_FAMILY = ("allele_exact", "false_positive", "critical_false_negative")
STRATA = (
    "profile",
    "event",
    "delta_class",
    "depth",
    "composition",
    "pcr",
    "smear",
    "chimera",
    "error",
)
ALLELE_UNIT = ("sample", "allele")

RULE_TEXT = (
    "MucSim-Bench decision rule v2 (spec section 6). Adopt the candidate engine over the "
    "baseline only if, for every profile: (1) it is superior on per-allele exact sequence "
    "(unit: each truth allele of each case, exact under the least favourable optimal "
    "assignment with independent haplotype evidence), by exact two-sided McNemar on "
    "allele pairs (same (design_id, truth allele) keys for both engines), Holm-adjusted "
    "across the primary family {allele_exact, false_positive, critical_false_negative} at "
    "alpha 0.05, with more candidate-only than baseline-only exact alleles; (2) it is "
    "non-inferior on the false-positive PATHOGENIC rate among normal and benign truths "
    "(Newcombe hybrid-score one-sided 95% upper bound of candidate minus baseline below "
    "0.005); and (3) its count of pathogenic truths called NO_PATHOGENIC_VARIANT_DETECTED "
    "or NO_CALL does not exceed the baseline's. Pooled per-allele and case-exact rates are "
    "reported with 95% cluster-bootstrap intervals over design_id (2000 replicates, seed "
    "0). Failed or unattempted runs count as NO_CALL with every truth allele not exact; no "
    "case or allele is dropped."
)


def _alleles(sample: dict[str, Any], failed: bool) -> list[dict[str, Any]]:
    """Per-truth-allele exactness from the least favourable optimal assignment."""
    name = sample["sample"]
    alternatives = sample.get("alternatives") or []
    if sample.get("truth_status") != "valid" or not alternatives:
        raise ValueError(f"{name}: no valid truth; truth alleles cannot be enumerated")
    first = alternatives[0]
    names = sorted({d["truth"] for d in first["pairs"]} | set(first["missing_truth"]))
    if len(names) != int(sample["truth_haplotypes"]):
        raise ValueError(
            f"{name}: {len(names)} truth allele names for {sample['truth_haplotypes']}"
        )
    exact: set[str] = set()
    if not failed:
        worst = min(
            alternatives,
            key=lambda a: (
                a["metrics"]["independent_sequence_exact"],
                a["metrics"]["sequence_exact"],
            ),
        )
        hits = sorted(d["truth"] for d in worst["pairs"] if d["sequence_exact"])
        exact = set(hits[: int(worst["metrics"]["independent_sequence_exact"])])
    return [{"allele": n, "allele_exact": int(n in exact)} for n in names]


def normalize_rows(
    report: dict[str, Any], cases: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    """One case row per evaluated sample (with per-allele entries), joined to its design.

    Raises:
        KeyError: If a sample has no case (strata would silently be lost).
        ValueError: If a sample has no valid truth or inconsistent truth alleles.
    """
    rows = []
    for sample in report["samples"]:
        name = sample["sample"]
        case = cases[name]
        design = case.get("design") or {}
        failed = sample["status"] in FAILED_STATUSES
        alleles = _alleles(sample, failed)
        clinical = sample.get("clinical") or {}
        truth = clinical.get("truth")
        if truth is None:
            raise ValueError(f"{name}: no clinical truth class")
        decision = "NO_CALL" if failed else str(clinical.get("decision") or "NO_CALL")
        metrics = sample["metrics"]
        fp = int(decision == "PATHOGENIC" and truth in ("normal", "benign"))
        cfn = int(truth == "pathogenic" and decision in FALSE_NEGATIVE_DECISIONS)
        case_exact = 0 if failed else int(metrics["all_sequences_exact"]["min"])
        row: dict[str, Any] = {
            "sample": name,
            "status": sample["status"],
            "failed": failed,
            "truth": truth,
            "decision": decision,
            "normal": truth == "normal",
            "benign": truth == "benign",
            "pathogenic": truth == "pathogenic",
            "alleles": alleles,
            "alleles_exact": sum(a["allele_exact"] for a in alleles),
            "truth_haplotypes": len(alleles),
            "case_exact": case_exact,
            "truth_events": int(sample["truth_events"]),
            "event_tp": 0 if failed else int(metrics["event_tp"]["min"]),
            "event_fp": 0 if failed else int(metrics["event_fp"]["max"]),
            "false_positive": fp,
            "critical_false_negative": cfn,
            "failure": int(not case_exact or fp or cfn),
            "inconclusive": int(decision == "INCONCLUSIVE"),
            "no_call": int(decision == "NO_CALL"),
            "event": design.get("event"),
            "realized_depth": _depth_sum(case.get("realized_depth")),
        }
        for key in STRATA:
            row[key] = case.get(key, design.get(key)) if key == "profile" else design.get(key)
        rows.append(row)
    return rows


def _depth_sum(realized: Any) -> int | None:
    if isinstance(realized, dict):
        return int(sum(realized.values()))
    return None if realized is None else int(realized)


def allele_rows(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Flatten case rows to one row per (sample, truth allele) with the case strata."""
    out = []
    for row in rows:
        base = {k: row.get(k) for k in ("sample", "event", "truth", *STRATA)}
        out += [base | a for a in row["alleles"]]
    return out


def _key(value: Any) -> str:
    return "" if value is None else str(value)


def stratified_table(
    rows: Sequence[dict[str, Any]], by: Sequence[str], metric: str
) -> list[dict[str, Any]]:
    """Successes ``k`` of ``n`` rows, rate and 95% Clopper-Pearson interval per stratum."""
    groups: dict[tuple[Any, ...], list[int]] = {}
    for row in rows:
        stratum = tuple(row.get(k) for k in by)
        groups.setdefault(stratum, []).append(int(row[metric]))
    table = []
    for stratum in sorted(groups, key=lambda s: tuple(_key(v) for v in s)):
        values = groups[stratum]
        k, n = sum(values), len(values)
        low, high = clopper_pearson(k, n)
        table.append(
            {
                "stratum": dict(zip(by, stratum, strict=True)),
                "metric": metric,
                "k": k,
                "n": n,
                "rate": k / n,
                "ci_low": low,
                "ci_high": high,
            }
        )
    return table


def paired(
    rows_a: Sequence[dict[str, Any]],
    rows_b: Sequence[dict[str, Any]],
    metric: str,
    unit: Sequence[str] = ("sample",),
) -> dict[str, Any]:
    """Exact McNemar on paired units: ``b`` = a-only successes, ``c`` = b-only successes.

    Units are keyed by the ``unit`` fields (``ALLELE_UNIT`` for allele pairs).

    Raises:
        ValueError: If the two engines were not scored on the same unique unit keys.
    """
    a = {tuple(r[u] for u in unit): int(r[metric]) for r in rows_a}
    b_map = {tuple(r[u] for u in unit): int(r[metric]) for r in rows_b}
    if len(a) != len(rows_a) or len(b_map) != len(rows_b) or a.keys() != b_map.keys():
        raise ValueError(f"paired comparison needs the same unique {'/'.join(unit)} set ({metric})")
    b = sum(a[s] == 1 and b_map[s] == 0 for s in a)
    c = sum(a[s] == 0 and b_map[s] == 1 for s in a)
    return {
        "metric": metric,
        "unit": list(unit),
        "n": len(a),
        "k_a": sum(a.values()),
        "k_b": sum(b_map.values()),
        "b": b,
        "c": c,
        "p": mcnemar_exact(b, c),
    }


def _fp_test(base: list[dict[str, Any]], cand: list[dict[str, Any]]) -> dict[str, Any]:
    def eligible(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        return [r for r in rows if r.get("normal") or r.get("benign")]

    b_rows, c_rows = eligible(base), eligible(cand)
    out: dict[str, Any] = {
        "denominator": "normal + benign truths",
        "n": len(c_rows),
        "n_baseline": len(b_rows),
        "fp": sum(r["false_positive"] for r in c_rows),
        "fp_baseline": sum(r["false_positive"] for r in b_rows),
        "no_call": sum(int(r.get("no_call", 0)) for r in c_rows),
        "no_call_baseline": sum(int(r.get("no_call", 0)) for r in b_rows),
        "margin": NI_MARGIN,
    }
    if not b_rows or not c_rows:
        return out | {
            "diff": None,
            "upper": None,
            "noninferior": False,
            "reason": "no normal or benign truths in this profile",
        }
    test = noninferior(out["fp"], len(c_rows), out["fp_baseline"], len(b_rows), NI_MARGIN, ALPHA)
    return out | test


def _profile(base: list[dict[str, Any]], cand: list[dict[str, Any]]) -> dict[str, Any]:
    tests = {
        "allele_exact": paired(allele_rows(base), allele_rows(cand), "allele_exact", ALLELE_UNIT)
    }
    for metric in HOLM_FAMILY[1:]:
        tests[metric] = paired(base, cand, metric)
    adjusted = holm({m: t["p"] for m, t in tests.items()})
    for metric, test in tests.items():
        test["p_holm"] = adjusted[metric]
    exact = tests["allele_exact"]
    exact["superior"] = exact["c"] > exact["b"] and exact["p_holm"] < ALPHA
    fp = _fp_test(base, cand) | {"paired": tests["false_positive"]}
    cfn = tests["critical_false_negative"]
    cfn["pass"] = cfn["k_b"] <= cfn["k_a"]
    return {
        "n": len(cand),
        "n_alleles": exact["n"],
        "allele_exact": exact,
        "false_positive": fp,
        "critical_false_negative": cfn,
        "pass": bool(exact["superior"] and fp["noninferior"] and cfn["pass"]),
    }


def decide(
    reports: dict[str, list[dict[str, Any]]], baseline: str, candidate: str
) -> dict[str, Any]:
    """Apply ``RULE_TEXT`` per profile to case rows; adopt only if every profile passes."""
    base, cand = reports[baseline], reports[candidate]
    profiles = sorted({_key(r.get("profile")) for r in (*base, *cand)})
    result: dict[str, Any] = {
        "baseline": baseline,
        "candidate": candidate,
        "rule_sha256": rule_sha256(RULE_TEXT),
        "holm_family": list(HOLM_FAMILY),
        "alpha": ALPHA,
        "margin": NI_MARGIN,
        "profiles": {},
    }
    for profile in profiles:
        b = [r for r in base if _key(r.get("profile")) == profile]
        c = [r for r in cand if _key(r.get("profile")) == profile]
        result["profiles"][profile] = _profile(b, c)
    result["adopt"] = bool(profiles) and all(p["pass"] for p in result["profiles"].values())
    return result


def fmt_value(value: Any) -> str:
    """Compact text for a table cell (``n/a`` for missing, 4 significant digits)."""
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def _table_md(title: str, table: Sequence[dict[str, Any]]) -> list[str]:
    if not table:
        return [f"### {title}", "", "(no rows)", ""]
    keys = list(table[0]["stratum"])
    lines = [
        f"### {title}",
        "",
        "| " + " | ".join([*keys, "k", "n", "rate", "95% CI"]) + " |",
        "|" + "---|" * (len(keys) + 4),
    ]
    for row in table:
        cells = [fmt_value(row["stratum"][k]) for k in keys]
        ci = f"[{fmt_value(row['ci_low'])}, {fmt_value(row['ci_high'])}]"
        lines.append(
            "| "
            + " | ".join([*cells, str(row["k"]), str(row["n"]), fmt_value(row["rate"]), ci])
            + " |"
        )
    return [*lines, ""]


def render_markdown(result: dict[str, Any]) -> str:
    """Markdown summary of a ``decide`` result plus optional ``result["tables"]``."""
    verdict = "ADOPT" if result.get("adopt") else "NOT ADOPTED"
    lines = ["# MucSim-Bench report", ""]
    if "candidate" in result:
        lines += [
            f"Candidate `{result['candidate']}` vs baseline `{result['baseline']}`: **{verdict}**",
            "",
            f"Rule sha256: `{result.get('rule_sha256')}`",
            "",
        ]
    else:
        lines += [f"Decision: **{verdict}**", ""]
    if result.get("profiles"):
        lines += [
            "## Decision rule per profile",
            "",
            "Metric 1 is per-allele exact sequence (McNemar on allele pairs); FP is the "
            "PATHOGENIC rate over normal + benign truths, shown next to their no-call rate.",
            "",
            "| profile | cases | alleles | exact alleles base/cand | b/c | p (Holm) | superior | "
            "FP/(normal+benign) cand vs base | no-call on normal+benign cand vs base | "
            "FP diff upper | non-inferior | crit. FN cand/base | pass |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        for name, prof in result["profiles"].items():
            ex, fp, cfn = (
                prof["allele_exact"],
                prof["false_positive"],
                prof["critical_false_negative"],
            )
            lines.append(
                f"| {name} | {prof['n']} | {prof['n_alleles']} | {ex['k_a']}/{ex['k_b']} | "
                f"{ex['b']}/{ex['c']} | {fmt_value(ex['p_holm'])} | {ex['superior']} | "
                f"{fp['fp']}/{fp['n']} vs {fp['fp_baseline']}/{fp['n_baseline']} | "
                f"{fp['no_call']}/{fp['n']} vs {fp['no_call_baseline']}/{fp['n_baseline']} | "
                f"{fmt_value(fp['upper'])} | {fp['noninferior']} | {cfn['k_b']}/{cfn['k_a']} | "
                f"{prof['pass']} |"
            )
        lines.append("")
    return "\n".join(lines) + "\n" + render_tables(result.get("tables") or {})


def render_tables(tables: Mapping[str, Sequence[dict[str, Any]]]) -> str:
    """Markdown for named ``stratified_table`` outputs."""
    lines: list[str] = []
    for title, table in tables.items():
        lines += _table_md(title, table)
    return "\n".join(lines) + ("\n" if lines else "")
