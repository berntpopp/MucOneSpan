"""`benchsim calibrate-report`: rank a calibration's points under a declared objective.

Reads ``calibration.json`` and each evaluated point's ``evaluation.json``,
normalizes the rows (`report.normalize_rows`, or, for a Task 15d ``lengths``-stage
calibration, `calibration_lengths.normalize_lengths_rows`), computes the objective's
metrics (`calibration_objective.point_metrics`, cluster-bootstrap CIs, over the
matching metric registry) and ranks the points (`calibration_objective.rank_points`,
unchanged either way). Writes, next to ``calibration.json``:

- ``calibration-report.json`` and ``calibration-report.md``: every point's
  values, metrics and CIs, feasibility, rank and the recommendation;
- ``recommended-config.json``: the best feasible point's overlay, byte-identical
  to its ``config.json`` and loadable by ``muconespan --config``;
- ``recommended-config.provenance.json``: grid, objective, split, engine and
  hashes behind the recommendation (a sidecar, because the strict settings
  loader rejects unknown fields in the config itself);
- ``recommended-grid.json``: the recommended point as a single-point grid, for
  the confirmation run (``calibrate --split val --grid recommended-grid.json``).

With ``--shift-from NAME`` on the confirmation split, each point is compared
with the same overlay (same content hash) in the calibration-split
calibration ``NAME``: the dev -> val shift of every objective metric.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from muc_one_span.benchsim.bench_config import BenchConfig
from muc_one_span.benchsim.calibration import (
    CALIBRATION_SPLIT,
    CONFIRMATION_SPLIT,
    EVALUATED,
    EVALUATION_FILE,
    MANIFEST_FILE,
    RESULTS_DIR,
    SCHEMA_VERSION,
    calibration_dir,
    check_split,
    dump,
    file_sha256,
    write_atomic,
)
from muc_one_span.benchsim.calibration_grid import LENGTHS_STAGE, OVERLAY_FILE, canonical_sha256
from muc_one_span.benchsim.calibration_lengths import (
    LENGTHS_COUNTS,
    LENGTHS_METRICS,
    LENGTHS_RATES,
    normalize_lengths_rows,
)
from muc_one_span.benchsim.calibration_objective import (
    METRICS,
    Objective,
    load_objective,
    point_metrics,
    rank_points,
)
from muc_one_span.benchsim.report import fmt_value, normalize_rows

REPORT_JSON = "calibration-report.json"
REPORT_MD = "calibration-report.md"
RECOMMENDED_CONFIG = "recommended-config.json"
RECOMMENDED_PROVENANCE = "recommended-config.provenance.json"
RECOMMENDED_GRID = "recommended-grid.json"


def load_calibration(cal_dir: Path) -> dict[str, Any]:
    """``calibration.json`` of a calibration directory (exit if it was never run)."""
    path = cal_dir / MANIFEST_FILE
    if not path.is_file():
        raise SystemExit(f"missing {path}; run `benchsim calibrate` first")
    data: dict[str, Any] = json.loads(path.read_text())
    return data


def _objective_for(stage: str, objective_path: Path) -> Objective:
    """The full-pipeline or lengths-stage metric registry, by the calibration's own stage."""
    known = LENGTHS_METRICS if stage == LENGTHS_STAGE else METRICS
    return load_objective(objective_path, known_metrics=known)


def scored_points(
    cal_dir: Path,
    calibration: dict[str, Any],
    cases: dict[str, dict[str, Any]],
    objective: Objective,
    bench: BenchConfig,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """(evaluated points with metrics, points without an evaluation).

    Dispatches on ``calibration["inputs"]["stage"]`` (absent: the full pipeline, 15b/15c
    behaviour unchanged): a ``lengths`` calibration's ``evaluation.json`` is normalized
    and scored by `calibration_lengths`'s row schema and metric registries instead of
    `report.normalize_rows`/the built-in `calibration_objective` registries, through the
    same (now stage-agnostic) `calibration_objective.point_metrics`.
    """
    engine = calibration["inputs"]["engine"]
    stage = calibration["inputs"].get("stage", "full")
    scored, missing = [], []
    for entry in calibration["points"]:
        path = cal_dir / entry["sha256"] / RESULTS_DIR / engine / EVALUATION_FILE
        if entry["status"] != EVALUATED or not path.is_file():
            missing.append({k: entry[k] for k in ("sha256", "values", "status", "error")})
            continue
        try:
            raw = json.loads(path.read_text())
            if stage == LENGTHS_STAGE:
                rows = normalize_lengths_rows(raw, cases, bench.sets.legacy)
                metrics = point_metrics(
                    rows, objective, bench.report, LENGTHS_RATES, LENGTHS_COUNTS
                )
            else:
                rows = normalize_rows(raw, cases, bench.sets.legacy)
                metrics = point_metrics(rows, objective, bench.report)
        except (KeyError, ValueError) as exc:
            raise SystemExit(f"cannot normalize {path}: {exc}") from exc
        scored.append({"sha256": entry["sha256"], "values": entry["values"], "metrics": metrics})
    return scored, missing


def shift(
    points: list[dict[str, Any]], reference: list[dict[str, Any]], objective: Objective
) -> dict[str, Any]:
    """Per point: each objective metric on the reference split, here, and the difference."""
    by_hash = {p["sha256"]: p for p in reference}
    out: dict[str, Any] = {}
    for point in points:
        ref = by_hash.get(point["sha256"])
        if ref is None:
            out[point["sha256"]] = None
            continue
        entry = {}
        for name in objective.metric_names():
            before, after = ref["metrics"][name]["value"], point["metrics"][name]["value"]
            delta = None if before is None or after is None else after - before
            entry[name] = {CALIBRATION_SPLIT: before, CONFIRMATION_SPLIT: after, "delta": delta}
        out[point["sha256"]] = entry
    return out


def _cell(metric: dict[str, Any]) -> str:
    text = fmt_value(metric["value"])
    if metric.get("ci_low") is None:
        return text
    return f"{text} [{fmt_value(metric['ci_low'])}, {fmt_value(metric['ci_high'])}]"


def render(report: dict[str, Any], objective: Objective) -> str:
    """Markdown for ``calibration-report.json``."""
    names = objective.metric_names()
    lines = [
        f"# Calibration `{report['name']}` ({report['split']}, engine `{report['engine']}`, "
        f"stage `{report['stage']}`)",
        "",
        f"- Objective sha256: `{report['objective_sha256']}`",
        f"- Calibration inputs sha256: `{report['calibration_inputs_sha256']}`",
        f"- Constraints: `{json.dumps(objective.raw.get('constraints', {}), sort_keys=True)}`",
        f"- Rank: `{', '.join(objective.raw['rank'])}` (then point hash)",
        f"- Recommended: `{report['recommended'] or 'none (no feasible point)'}`",
        "",
        "| rank | point | values | feasible | " + " | ".join(names) + " | violations |",
        "|---|---|---|---|" + "---|" * len(names) + "---|",
    ]
    for point in report["points"]:
        cells = [_cell(point["metrics"][name]) for name in names]
        lines.append(
            f"| {point['rank']} | `{point['sha256']}` | "
            f"`{json.dumps(point['values'], sort_keys=True)}` | {point['feasible']} | "
            + " | ".join(cells)
            + f" | {'; '.join(point['violations'])} |"
        )
    if report["not_evaluated"]:
        lines += ["", "## Not evaluated", ""]
        lines += [f"- `{p['sha256']}` {p['status']}: {p['error']}" for p in report["not_evaluated"]]
    if report["shift"] is not None:
        source = report["shift"]["from"]
        lines += ["", f"## {CALIBRATION_SPLIT} -> {CONFIRMATION_SPLIT} shift (from `{source}`)"]
        lines += ["", "| point | metric | dev | val | delta |", "|---|---|---|---|---|"]
        for sha, entry in report["shift"]["points"].items():
            if entry is None:
                lines.append(f"| `{sha}` | no {CALIBRATION_SPLIT} point | | | |")
                continue
            lines += [
                f"| `{sha}` | {name} | {fmt_value(v[CALIBRATION_SPLIT])} | "
                f"{fmt_value(v[CONFIRMATION_SPLIT])} | {fmt_value(v['delta'])} |"
                for name, v in entry.items()
            ]
    return "\n".join(lines) + "\n"


def write_recommendation(
    cal_dir: Path, report: dict[str, Any], calibration: dict[str, Any], objective: Objective
) -> None:
    """Recommended config, its provenance sidecar and single-point grid (or remove stale ones)."""
    paths = [cal_dir / n for n in (RECOMMENDED_CONFIG, RECOMMENDED_PROVENANCE, RECOMMENDED_GRID)]
    best = report["recommended"]
    if best is None:
        for path in paths:
            path.unlink(missing_ok=True)
        return
    point = next(p for p in report["points"] if p["sha256"] == best)
    config_text = (cal_dir / best / OVERLAY_FILE).read_text()
    write_atomic(paths[0], config_text)
    inputs = calibration["inputs"]
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "recommended_config_sha256": file_sha256(paths[0]),
        "point_sha256": best,
        "values": point["values"],
        "split": inputs["split"],
        "name": inputs["name"],
        "engine": inputs["engine"],
        "grid": inputs["grid"],
        "grid_sha256": inputs["grid_sha256"],
        "objective": objective.raw,
        "objective_sha256": objective.sha256,
        "base_config": inputs["base_config"],
        "base_config_sha256": inputs["base_config_sha256"],
        "manifest_sha256": inputs["manifest_sha256"],
        "calibration_inputs_sha256": report["calibration_inputs_sha256"],
        "versions": inputs["versions"],
        "metrics": point["metrics"],
    }
    write_atomic(paths[1], dump(provenance))
    write_atomic(paths[2], dump({key: [value] for key, value in point["values"].items()}))


def build_report(
    out_root: Path,
    split: str,
    name: str,
    objective_path: Path,
    cases: dict[str, dict[str, Any]],
    bench: BenchConfig,
    shift_from: tuple[str, dict[str, dict[str, Any]]] | None = None,
) -> int:
    """Rank, write the report files and the recommendation; 1 when nothing is feasible.

    ``shift_from`` is (calibration name, its split's cases) on the calibration split.
    The objective is loaded here (after ``calibration.json``, hence by path, not
    pre-loaded) because a ``lengths``-stage calibration validates it against a
    different metric registry (`_objective_for`).
    """
    check_split(split)
    if shift_from is not None and split != CONFIRMATION_SPLIT:
        raise SystemExit(f"--shift-from compares with {CALIBRATION_SPLIT}; use --split val")
    cal_dir = calibration_dir(out_root, split, name)
    calibration = load_calibration(cal_dir)
    stage = calibration["inputs"].get("stage", "full")
    try:
        objective = _objective_for(stage, objective_path)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"--objective: {exc}") from exc
    scored, missing = scored_points(cal_dir, calibration, cases, objective, bench)
    ranked = rank_points(scored, objective)
    feasible = [p for p in ranked if p["feasible"]]
    shift_section = None
    if shift_from is not None:
        ref_dir = calibration_dir(out_root, CALIBRATION_SPLIT, shift_from[0])
        reference, _ = scored_points(
            ref_dir, load_calibration(ref_dir), shift_from[1], objective, bench
        )
        shift_section = {"from": shift_from[0], "points": shift(ranked, reference, objective)}
    inputs = calibration["inputs"]
    report = {
        "schema_version": SCHEMA_VERSION,
        "split": split,
        "name": name,
        "engine": inputs["engine"],
        "stage": stage,
        "calibration_inputs_sha256": canonical_sha256(inputs),
        "objective": objective.raw,
        "objective_sha256": objective.sha256,
        "bench_config_sha256": bench.sha256(),
        "points": ranked,
        "not_evaluated": missing,
        "recommended": feasible[0]["sha256"] if feasible else None,
        "shift": shift_section,
    }
    write_atomic(cal_dir / REPORT_JSON, dump(report))
    write_atomic(cal_dir / REPORT_MD, render(report, objective))
    write_recommendation(cal_dir, report, calibration, objective)
    print(f"calibration report: {cal_dir / REPORT_JSON}")
    if report["recommended"] is None:
        print("no grid point satisfies the objective's constraints")
        return 1
    print(f"recommended: {cal_dir / RECOMMENDED_CONFIG}")
    return 0
