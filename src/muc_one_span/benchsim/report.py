"""Stratified benchmark report, paired engine comparison and the pre-registered decision rule.

Consumes ``scripts/evaluate.py`` reports (one per engine; per-sample rows under
``report["samples"]``) and generated ``case.json`` dicts keyed by
``design_id``. ``normalize_rows`` reduces each sample row to the binary
endpoints of spec section 6; ``decide`` applies the decision rule
(``RULE_TEXT``) per profile.

Key mapping from ``evaluation.scoring.evaluate_sample`` rows:

- ``allele_exact``: ``metrics["all_sequences_exact"]["min"]`` -- every truth
  allele sequence-exact with the right allele cardinality and independent
  haplotype evidence, under the least favourable optimal assignment (the same
  numerator as ``totals["exact_sample_reconstruction"]``).
- ``alleles_exact``: ``metrics["independent_sequence_exact"]["min"]`` (count of
  exact alleles, for per-allele rates over ``truth_haplotypes``).
- ``truth`` / ``decision``: ``row["clinical"]`` (``clinical_confusion``); a
  missing truth (invalid truth) falls back to the design's event class.

Runs with status ``execution_failed``, ``not_attempted``, ``invalid_truth`` or
``invalid_artifacts`` are kept: they score 0 on every exactness endpoint and
their decision is ``NO_CALL`` (so a pathogenic truth becomes a critical false
negative), whatever stale files the result directory holds.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from muc_one_span.benchsim.stats import clopper_pearson, holm, mcnemar_exact, noninferior

FAILED_STATUSES = frozenset(
    {"execution_failed", "not_attempted", "invalid_truth", "invalid_artifacts"}
)
FALSE_NEGATIVE_DECISIONS = frozenset({"NO_PATHOGENIC_VARIANT_DETECTED", "NO_CALL"})
NI_MARGIN = 0.005
ALPHA = 0.05
HOLM_FAMILY = ("allele_exact", "false_positive", "critical_false_negative")
STRATA = ("profile", "delta_class", "depth", "composition", "pcr", "smear", "chimera", "error")

RULE_TEXT = (
    "MucSim-Bench decision rule v1 (spec section 6). Adopt the candidate engine over the "
    "baseline only if, for every profile: (1) it is superior on per-allele exact sequence "
    "(sample-level: all truth alleles sequence-exact), by exact two-sided McNemar on paired "
    "samples, Holm-adjusted across the primary family {allele_exact, false_positive, "
    "critical_false_negative} at alpha 0.05, with more candidate-only than baseline-only "
    "successes; (2) it is non-inferior on the false-positive PATHOGENIC rate among normal "
    "and benign truths (Newcombe hybrid-score one-sided 95% upper bound of candidate minus "
    "baseline below 0.005); and (3) its count of pathogenic truths called "
    "NO_PATHOGENIC_VARIANT_DETECTED or NO_CALL does not exceed the baseline's. Failed or "
    "unattempted runs count as NO_CALL and not exact; no sample is dropped."
)


def _design_class(design: dict[str, Any]) -> str:
    """Clinical truth class from the design event (for rows without valid truth)."""
    event = design.get("event")
    if not event:
        return "normal"
    from muc_one_span.config import load_repeat_dictionary
    from muc_one_span.evaluation.clinical_confusion import net_length_change

    definition = load_repeat_dictionary().mutations[event]
    return "pathogenic" if net_length_change(definition) % 3 else "benign"


def _depth_sum(realized: Any) -> int | None:
    if isinstance(realized, dict):
        return int(sum(realized.values()))
    return None if realized is None else int(realized)


def normalize_rows(
    report: dict[str, Any], cases: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    """One binary-endpoint row per evaluated sample, joined to its case design.

    Raises:
        KeyError: If a sample has no case (strata would silently be lost).
    """
    rows = []
    for sample in report["samples"]:
        name = sample["sample"]
        case = cases[name]
        design = case.get("design") or {}
        failed = sample["status"] in FAILED_STATUSES or sample.get("truth_status") != "valid"
        clinical = sample.get("clinical") or {}
        truth = clinical.get("truth") or _design_class(design)
        decision = "NO_CALL" if failed else str(clinical.get("decision") or "NO_CALL")
        metrics = sample.get("metrics") or {}
        exact = 0 if failed else int(metrics["all_sequences_exact"]["min"])
        alleles = 0 if failed else int(metrics["independent_sequence_exact"]["min"])
        row: dict[str, Any] = {
            "sample": name,
            "status": sample["status"],
            "failed": failed,
            "truth": truth,
            "decision": decision,
            "normal": truth == "normal",
            "benign": truth == "benign",
            "pathogenic": truth == "pathogenic",
            "allele_exact": exact,
            "alleles_exact": alleles,
            "truth_haplotypes": int(
                sample.get("truth_haplotypes") or len(design.get("lengths") or ()) or 2
            ),
            "false_positive": int(decision == "PATHOGENIC" and truth in ("normal", "benign")),
            "critical_false_negative": int(
                truth == "pathogenic" and decision in FALSE_NEGATIVE_DECISIONS
            ),
            "inconclusive": int(decision == "INCONCLUSIVE"),
            "no_call": int(decision == "NO_CALL"),
            "event": design.get("event"),
            "realized_depth": _depth_sum(case.get("realized_depth")),
        }
        for key in STRATA:
            row[key] = case.get(key, design.get(key)) if key == "profile" else design.get(key)
        rows.append(row)
    return rows


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
    rows_a: Sequence[dict[str, Any]], rows_b: Sequence[dict[str, Any]], metric: str
) -> dict[str, Any]:
    """Exact McNemar on paired samples: ``b`` = a-only successes, ``c`` = b-only successes.

    Raises:
        ValueError: If the two engines were not scored on the same samples.
    """
    a = {r["sample"]: int(r[metric]) for r in rows_a}
    b_map = {r["sample"]: int(r[metric]) for r in rows_b}
    if len(a) != len(rows_a) or len(b_map) != len(rows_b) or a.keys() != b_map.keys():
        raise ValueError(f"paired comparison needs the same unique sample set ({metric})")
    b = sum(a[s] == 1 and b_map[s] == 0 for s in a)
    c = sum(a[s] == 0 and b_map[s] == 1 for s in a)
    return {
        "metric": metric,
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
    fp_b, fp_c = sum(r["false_positive"] for r in b_rows), sum(r["false_positive"] for r in c_rows)
    out: dict[str, Any] = {
        "n": len(c_rows),
        "n_baseline": len(b_rows),
        "fp": fp_c,
        "fp_baseline": fp_b,
        "margin": NI_MARGIN,
    }
    if not b_rows or not c_rows:
        return out | {
            "diff": None,
            "upper": None,
            "noninferior": False,
            "reason": "no normal or benign truths in this profile",
        }
    return out | noninferior(fp_c, len(c_rows), fp_b, len(b_rows), margin=NI_MARGIN, alpha=ALPHA)


def _profile(base: list[dict[str, Any]], cand: list[dict[str, Any]]) -> dict[str, Any]:
    tests = {metric: paired(base, cand, metric) for metric in HOLM_FAMILY}
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
        "allele_exact": exact,
        "false_positive": fp,
        "critical_false_negative": cfn,
        "pass": bool(exact["superior"] and fp["noninferior"] and cfn["pass"]),
    }


def decide(
    reports: dict[str, list[dict[str, Any]]], baseline: str, candidate: str
) -> dict[str, Any]:
    """Apply ``RULE_TEXT`` per profile; adopt only if every profile passes (and one exists)."""
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


def rule_sha256(rule_text: str) -> str:
    """SHA-256 of the exact rule text (UTF-8)."""
    return hashlib.sha256(rule_text.encode("utf-8")).hexdigest()


def preregister(rule_text: str, path: Path) -> str:
    """Append (never rewrite) a pre-registration entry for ``rule_text``; return its sha256."""
    digest = rule_sha256(rule_text)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "sha256": digest,
        "rule_text": rule_text,
        "registered_at": datetime.now(timezone.utc).isoformat(),
    }
    with path.open("a", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return digest


def require_preregistered(path: Path, rule_text: str) -> None:
    """Raise ``PermissionError`` unless ``rule_text`` is pre-registered in ``path``."""
    digest = rule_sha256(rule_text)
    path = Path(path)
    if not path.is_file():
        raise PermissionError(f"no pre-registration ledger at {path}; run `benchsim preregister`")
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            entry = json.loads(line)
        except ValueError as exc:
            raise PermissionError(f"corrupt pre-registration ledger {path}: {line!r}") from exc
        if entry.get("sha256") == digest and rule_sha256(str(entry.get("rule_text"))) == digest:
            return
    raise PermissionError(f"decision rule sha256 {digest} is not pre-registered in {path}")


def _fmt(value: Any) -> str:
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
        cells = [_fmt(row["stratum"][k]) for k in keys]
        ci = f"[{_fmt(row['ci_low'])}, {_fmt(row['ci_high'])}]"
        lines.append(
            "| " + " | ".join([*cells, str(row["k"]), str(row["n"]), _fmt(row["rate"]), ci]) + " |"
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
            "| profile | n | exact b/c | exact p (Holm) | superior | FP cand/base | "
            "FP diff upper | non-inferior | crit. FN cand/base | pass |",
            "|---|---|---|---|---|---|---|---|---|---|",
        ]
        for name, prof in result["profiles"].items():
            ex, fp, cfn = (
                prof["allele_exact"],
                prof["false_positive"],
                prof["critical_false_negative"],
            )
            lines.append(
                f"| {name} | {prof['n']} | {ex['b']}/{ex['c']} | {_fmt(ex['p_holm'])} | "
                f"{ex['superior']} | {fp['fp']}/{fp['n']} vs {fp['fp_baseline']}/{fp['n_baseline']} | "
                f"{_fmt(fp['upper'])} | {fp['noninferior']} | {cfn['k_b']}/{cfn['k_a']} | {prof['pass']} |"
            )
        lines.append("")
    return "\n".join(lines) + "\n" + render_tables(result.get("tables") or {})


def render_tables(tables: Mapping[str, Sequence[dict[str, Any]]]) -> str:
    """Markdown for named ``stratified_table`` outputs."""
    lines: list[str] = []
    for title, table in tables.items():
        lines += _table_md(title, table)
    return "\n".join(lines) + ("\n" if lines else "")
