"""Reason atlas: why cases end without a definitive clinical decision.

Input is ``report.normalize_rows`` output. Atlas cases are those whose
decision is in ``AtlasConfig.decisions`` (default INCONCLUSIVE). Each atlas
case contributes each of its reason keys once:

- ``gate: <key>``: one per caller reason in ``clinical_reasons`` (the
  INCONCLUSIVE banner details of ``compute_clinical_decision``), normalised
  by `reason_key`;
- ``evaluator: <flag>``: one per ``reconstruction_flags`` entry
  (``evaluation.reasons.reconstruction_flags``);
- `UNRECORDED` when the case has neither (or its evaluation predates reasons).

`reason_key` rule, applied in order: drop a leading allele label (``Allele 1:``
or ``allele_1:``); replace a variant descriptor ``(<name> at repeat <index>)``
with ``(<variant>)``; replace every standalone number and ``None`` (an unset
numeric field) with ``#``; lowercase, collapse whitespace and drop a trailing
period. Numbers inside identifiers (``r1041``) are kept.

A case is *expected* non-definitive when `expected_conditions` is nonempty:
``split`` (its split is in ``expected_inconclusive_splits``) and/or ``depth``
(its ``depth_basis`` depth is below ``min_resolvable_depth``). Other atlas
cases are *resolvable*: an engine could turn them definitive without relaxing
any clinical gate.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Sequence
from typing import Any

from muc_one_span.benchsim.bench_config import AtlasConfig
from muc_one_span.benchsim.report import fmt_value

ANY_REASON = "(any)"  # stratum row counting every atlas case
UNRECORDED = "unrecorded"
_GATE, _EVALUATOR = "gate: ", "evaluator: "
_ALLELE_LABEL = re.compile(r"^allele[ _]\d+:\s*", re.IGNORECASE)
_VARIANT = re.compile(r"\([^()]* at repeat [^()]*\)")
_NUMBER = re.compile(r"(?<![\w.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![\w.]*\w)|\bNone\b")
_DEPTH_FIELD = {"design": "depth", "realized_min_allele": "realized_min_allele_depth"}


def reason_key(text: str) -> str:
    """Normalise a caller reason string to its atlas key (module rule)."""
    key = _ALLELE_LABEL.sub("", " ".join(text.split()))
    key = _VARIANT.sub("(<variant>)", key)
    key = _NUMBER.sub("#", key)
    return key.lower().rstrip(".").strip()


def case_reason_keys(row: dict[str, Any]) -> list[str]:
    """Sorted unique reason keys of one case row (`UNRECORDED` when there are none)."""
    keys = {_GATE + reason_key(r) for r in row.get("clinical_reasons") or []}
    keys |= {_EVALUATOR + f for f in row.get("reconstruction_flags") or []}
    return sorted(keys) or [UNRECORDED]


def expected_conditions(row: dict[str, Any], split: str, cfg: AtlasConfig) -> list[str]:
    """The configured conditions that make this case expected non-definitive."""
    out = []
    if split in cfg.expected_inconclusive_splits:
        out.append("split")
    depth = row.get(_DEPTH_FIELD[cfg.depth_basis])
    if depth is not None and depth < cfg.min_resolvable_depth:
        out.append("depth")
    return out


def _summary(
    profile: str, rows: Sequence[dict[str, Any]], atlas: Sequence[dict[str, Any]]
) -> dict[str, Any]:
    conditions: Counter[str] = Counter(c for r in atlas for c in r["_expected"])
    expected = sum(bool(r["_expected"]) for r in atlas)
    return {
        "profile": profile,
        "cases": len(rows),
        "atlas": len(atlas),
        "expected": expected,
        "resolvable": len(atlas) - expected,
        "conditions": dict(sorted(conditions.items())),
    }


def _stratum(
    rows: Sequence[dict[str, Any]], atlas: Sequence[dict[str, Any]], factor: str
) -> list[dict[str, Any]]:
    cells = []
    groups = sorted(
        {(str(r.get("profile")), r.get(factor)) for r in rows},
        key=lambda g: (g[0], str(type(g[1])), "" if g[1] is None else g[1]),
    )
    for profile, value in groups:
        n = sum(str(r.get("profile")) == profile and r.get(factor) == value for r in rows)
        hits = [a for a in atlas if str(a.get("profile")) == profile and a.get(factor) == value]
        counts = Counter(k for a in hits for k in a["_keys"])
        counts[ANY_REASON] = len(hits)
        cells += [
            {"profile": profile, "value": value, "reason": reason, "k": k, "n": n, "rate": k / n}
            for reason, k in sorted(counts.items())
            if k or reason == ANY_REASON
        ]
    return cells


def build_atlas(rows: Sequence[dict[str, Any]], split: str, cfg: AtlasConfig) -> dict[str, Any]:
    """Reason counts, the expected/resolvable split and reason x profile x stratum cells."""
    atlas = [
        r | {"_keys": case_reason_keys(r), "_expected": expected_conditions(r, split, cfg)}
        for r in rows
        if r.get("decision") in cfg.decisions
    ]
    n_cases, n_atlas = len(rows), len(atlas)
    counts: Counter[str] = Counter(k for a in atlas for k in a["_keys"])
    expected: Counter[str] = Counter(k for a in atlas if a["_expected"] for k in a["_keys"])
    reasons = [
        {
            "reason": key,
            "k": k,
            "share": k / n_atlas,
            "rate": k / n_cases,
            "expected": expected[key],
            "resolvable": k - expected[key],
        }
        for key, k in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    ]
    profiles = sorted({str(r.get("profile")) for r in rows})
    summary = [
        _summary(
            p,
            [r for r in rows if str(r.get("profile")) == p],
            [a for a in atlas if str(a.get("profile")) == p],
        )
        for p in profiles
    ]
    return {
        "split": split,
        "decisions": list(cfg.decisions),
        "expected_inconclusive_splits": list(cfg.expected_inconclusive_splits),
        "depth_basis": cfg.depth_basis,
        "min_resolvable_depth": cfg.min_resolvable_depth,
        "n_cases": n_cases,
        "n_atlas": n_atlas,
        "split_summary": [*summary, _summary("all", rows, atlas)],
        "reasons": reasons,
        "by_stratum": {factor: _stratum(rows, atlas, factor) for factor in cfg.strata},
    }


def _cell(k: int, n: int) -> str:
    return f"{k} ({fmt_value(k / n)})" if k else "0"


def render_atlas(atlas: dict[str, Any], cfg: AtlasConfig) -> str:
    """Markdown for `build_atlas`: split summary, reason legend and top-reason strata."""
    label = "/".join(atlas["decisions"])
    lines = [f"### Reason atlas ({label})", ""]
    if not atlas["n_atlas"]:
        return "\n".join([*lines, f"(no {label} cases)", ""]) + "\n"
    splits = ", ".join(atlas["expected_inconclusive_splits"]) or "none"
    lines += [
        f"Expected {label}: split in [{splits}] or {atlas['depth_basis']} depth below "
        f"{atlas['min_resolvable_depth']}; every other {label} case is resolvable.",
        "",
        f"| profile | cases | {label} | expected | resolvable | expected by condition |",
        "|---|---|---|---|---|---|",
    ]
    for s in atlas["split_summary"]:
        lines.append(
            f"| {s['profile']} | {s['cases']} | {s['atlas']} | {s['expected']} | "
            f"{s['resolvable']} | {s['conditions'] or '-'} |"
        )
    top = atlas["reasons"][: cfg.top_reasons]
    ids = {r["reason"]: f"R{i}" for i, r in enumerate(top, start=1)}
    lines += [
        "",
        f"Reasons (a case counts once per reason; share of {label} cases, rate of all cases):",
        "",
        "| id | reason | cases | share | rate | expected | resolvable |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in atlas["reasons"]:
        lines.append(
            f"| {ids.get(r['reason'], '')} | {r['reason']} | {r['k']} | {fmt_value(r['share'])} | "
            f"{fmt_value(r['rate'])} | {r['expected']} | {r['resolvable']} |"
        )
    for factor, cells in atlas["by_stratum"].items():
        grid: dict[tuple[str, Any], dict[str, dict[str, Any]]] = {}
        for c in cells:
            grid.setdefault((c["profile"], c["value"]), {})[c["reason"]] = c
        head = ["profile", factor, "cases", ANY_REASON, *ids.values()]
        lines += [
            "",
            f"### By profile x {factor} (top {len(ids)} reasons, k (rate))",
            "",
            "| " + " | ".join(head) + " |",
            "|" + "---|" * len(head),
        ]
        for (profile, value), by in grid.items():
            n = by[ANY_REASON]["n"]
            row = [_cell(by[k]["k"], n) if k in by else "0" for k in (ANY_REASON, *ids)]
            lines.append(f"| {profile} | {fmt_value(value)} | {n} | " + " | ".join(row) + " |")
    return "\n".join(lines) + "\n\n"
