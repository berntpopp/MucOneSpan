"""Reason atlas: why cases end without a definitive clinical decision.

Input is ``report.normalize_rows`` output. Atlas cases are those whose
decision is in ``AtlasConfig.decisions`` (default INCONCLUSIVE). Each atlas
case contributes each of its reason keys once:

- ``gate: <key>``: one per caller reason in ``clinical_reasons`` (the banner
  details of ``compute_clinical_decision``), normalised by `reason_keys`;
- ``evaluator: <flag>``: one per ``reconstruction_flags`` entry
  (``evaluation.reasons.reconstruction_flags``);
- `NO_REASONS` when reasons were recorded (``reasons_recorded``) but both lists
  are empty, and `UNRECORDED` when the evaluation predates reason recording.

`reason_key` rule, applied in order: drop a leading allele label (``Allele 1:``
or ``allele_1:``); replace a variant descriptor ``(<name> at repeat <index>)``
with ``(<variant>)``; replace every standalone number and ``None`` (an unset
numeric field) with ``#``; lowercase, collapse whitespace and drop a trailing
period. Numbers inside identifiers (``r1041``) are kept. `reason_keys` then
splits an uncertain-variant reason (``... is inconclusive (<b1>; <b2>)``) into
one key per blocker, ``... is inconclusive: <b1>``, so each blocker is counted
on its own rather than once per blocker combination. Blockers are split on
``; `` at parenthesis depth 0 only, so a blocker's own ``(... ; ...)`` stays whole.

Each atlas case gets one class (`atlas_class`):

- *expected* when `expected_conditions` is nonempty: ``split`` (its split is in
  ``expected_inconclusive_splits``) and/or ``depth`` (its ``depth_basis`` depth
  is below ``min_resolvable_depth``);
- *depth_unknown* when no split condition holds and the ``depth_basis`` depth is
  not recorded, so the depth condition cannot be decided;
- *resolvable* otherwise: the recorded depth reaches the gate. With the default
  basis this is simulator (truth) spanning depth per allele; the caller's own
  primary-record count can still fall below the gate when reads are lost while
  alleles are split.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Sequence
from typing import Any

from muc_one_span.benchsim.bench_config import AtlasConfig
from muc_one_span.benchsim.report import fmt_value

ANY_REASON = "(any)"  # stratum row counting every atlas case
UNRECORDED = "unrecorded"  # evaluation written before reasons were recorded
NO_REASONS = "no reasons recorded"  # reasons recorded, but none given
EXPECTED, RESOLVABLE, DEPTH_UNKNOWN = "expected", "resolvable", "depth_unknown"
CLASSES = (EXPECTED, RESOLVABLE, DEPTH_UNKNOWN)
_GATE, _EVALUATOR = "gate: ", "evaluator: "
_ALLELE_LABEL = re.compile(r"^allele[ _]\d+:\s*", re.IGNORECASE)
_VARIANT = re.compile(r"\([^()]* at repeat [^()]*\)")
_NUMBER = re.compile(r"(?<![\w.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![\w.]*\w)|\bNone\b")
# Uncertain-variant reason after `reason_key`: "<head> is inconclusive (<b1>; <b2>)".
_BLOCKERS = re.compile(r"^(?P<head>.*? is inconclusive) \((?P<blockers>.*)\)$")
_BLOCKER_SEP = "; "  # separator of decision blockers in the caller's reason text
_DEPTH_FIELD = {"design": "depth", "realized_min_allele": "realized_min_allele_depth"}


def reason_key(text: str) -> str:
    """Normalise a caller reason string to its atlas key (module rule)."""
    key = _ALLELE_LABEL.sub("", " ".join(text.split()))
    key = _VARIANT.sub("(<variant>)", key)
    key = _NUMBER.sub("#", key)
    return key.lower().rstrip(".").strip()


def reason_keys(text: str) -> list[str]:
    """`reason_key`, split into one key per blocker for an uncertain-variant reason."""
    key = reason_key(text)
    match = _BLOCKERS.match(key)
    if match is None:
        return [key]
    return [f"{match['head']}: {b}" for b in _split_top_level(match["blockers"])]


def _split_top_level(text: str) -> list[str]:
    """Split on `_BLOCKER_SEP` outside parentheses."""
    parts, depth, start, i = [], 0, 0, 0
    while i < len(text):
        char = text[i]
        depth += (char == "(") - (char == ")")
        if depth == 0 and text.startswith(_BLOCKER_SEP, i):
            parts.append(text[start:i])
            i += len(_BLOCKER_SEP)
            start = i
            continue
        i += 1
    return [*parts, text[start:]]


def case_reason_keys(row: dict[str, Any]) -> list[str]:
    """Sorted unique reason keys of one case row (`NO_REASONS` / `UNRECORDED` if none)."""
    keys = {_GATE + k for r in row.get("clinical_reasons") or [] for k in reason_keys(r)}
    keys |= {_EVALUATOR + f for f in row.get("reconstruction_flags") or []}
    return sorted(keys) or [NO_REASONS if row.get("reasons_recorded") else UNRECORDED]


def expected_conditions(row: dict[str, Any], split: str, cfg: AtlasConfig) -> list[str]:
    """The configured conditions that make this case expected non-definitive."""
    out = []
    if split in cfg.expected_inconclusive_splits:
        out.append("split")
    depth = row.get(_DEPTH_FIELD[cfg.depth_basis])
    if depth is not None and depth < cfg.min_resolvable_depth:
        out.append("depth")
    return out


def atlas_class(row: dict[str, Any], split: str, cfg: AtlasConfig) -> str:
    """`EXPECTED`, `DEPTH_UNKNOWN` or `RESOLVABLE` (module doc)."""
    if expected_conditions(row, split, cfg):
        return EXPECTED
    if row.get(_DEPTH_FIELD[cfg.depth_basis]) is None:
        return DEPTH_UNKNOWN
    return RESOLVABLE


def _summary(
    profile: str, rows: Sequence[dict[str, Any]], atlas: Sequence[dict[str, Any]]
) -> dict[str, Any]:
    conditions: Counter[str] = Counter(c for r in atlas for c in r["_expected"])
    classes = Counter(r["_class"] for r in atlas)
    return {
        "profile": profile,
        "cases": len(rows),
        "atlas": len(atlas),
        **{c: classes[c] for c in CLASSES},
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
        r
        | {
            "_keys": case_reason_keys(r),
            "_expected": expected_conditions(r, split, cfg),
            "_class": atlas_class(r, split, cfg),
        }
        for r in rows
        if r.get("decision") in cfg.decisions
    ]
    n_cases, n_atlas = len(rows), len(atlas)
    counts: Counter[str] = Counter(k for a in atlas for k in a["_keys"])
    by_class = {
        c: Counter(k for a in atlas if a["_class"] == c for k in a["_keys"]) for c in CLASSES
    }
    reasons = [
        {
            "reason": key,
            "k": k,
            "share": k / n_atlas,
            "rate": k / n_cases,
            **{c: by_class[c][key] for c in CLASSES},
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
        f"{atlas['min_resolvable_depth']}. Depth unknown: no split condition and no "
        f"{atlas['depth_basis']} depth recorded. Every other {label} case is resolvable.",
        "",
        f"| profile | cases | {label} | expected | resolvable | depth unknown | "
        "expected by condition |",
        "|---|---|---|---|---|---|---|",
    ]
    for s in atlas["split_summary"]:
        lines.append(
            f"| {s['profile']} | {s['cases']} | {s['atlas']} | {s[EXPECTED]} | "
            f"{s[RESOLVABLE]} | {s[DEPTH_UNKNOWN]} | {s['conditions'] or '-'} |"
        )
    top = atlas["reasons"][: cfg.top_reasons]
    ids = {r["reason"]: f"R{i}" for i, r in enumerate(top, start=1)}
    lines += [
        "",
        f"Reasons (a case counts once per reason; share of {label} cases, rate of all cases):",
        "",
        "| id | reason | cases | share | rate | expected | resolvable | depth unknown |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in atlas["reasons"]:
        lines.append(
            f"| {ids.get(r['reason'], '')} | {r['reason']} | {r['k']} | {fmt_value(r['share'])} | "
            f"{fmt_value(r['rate'])} | {r[EXPECTED]} | {r[RESOLVABLE]} | {r[DEPTH_UNKNOWN]} |"
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
