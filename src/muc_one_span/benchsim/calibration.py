"""`benchsim calibrate`: run and score every grid point of a settings calibration.

Layout (outside Git, under the benchsim ``--out-root``)::

    calibration/<split>/<name>/
        calibration.json           inputs (grid, split, engine, hashes, seeds, versions)
                                   and per-point status
        <point-sha256>/config.json the point's settings overlay (``muconespan --config``)
        <point-sha256>/results/<engine>/...        `run_split` output
        <point-sha256>/results/<engine>/evaluation.json

Every point is validated before the first run (`calibration_grid.build_points`).
Points are content-addressed and resumable: a point whose status is
``evaluated`` and whose ``evaluation.json`` exists is not run again; a
``failed`` point is retried. A calibration name is bound to its inputs: a
rerun with a different grid, base config, engine, manifest or software version
is refused rather than mixing results. The sealed ``test`` split is refused
outright; calibration uses ``dev`` and confirms on ``val``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from muc_one_span.benchsim.calibration_grid import (
    OVERLAY_FILE,
    GridPoint,
    base_settings,
    build_points,
    canonical_sha256,
    load_grid,
)
from muc_one_span.version import __version__

SCHEMA_VERSION = 1
SEALED_SPLIT = "test"
CALIBRATION_SPLIT = "dev"
CONFIRMATION_SPLIT = "val"
CALIBRATION_DIR = "calibration"
MANIFEST_FILE = "calibration.json"
RESULTS_DIR = "results"
EVALUATION_FILE = "evaluation.json"
EVALUATED, FAILED, PENDING = "evaluated", "failed", "pending"
JSON_INDENT = 2  # output format of every JSON file written here

RunFn = Callable[[Path, Path, Path], Any]  # (manifest, results_root, overlay config)
EvaluateFn = Callable[[argparse.Namespace], tuple[dict[str, Any], int]]


@dataclass(frozen=True)
class CalibrationRequest:
    """What `run_calibration` needs from the command line."""

    split: str
    name: str
    engine: str
    grid: Path
    base_config: Path | None
    out_root: Path


def check_split(split: str) -> None:
    """Calibrate on ``dev``, confirm on ``val``; never touch the sealed ``test`` split."""
    if split == SEALED_SPLIT:
        raise SystemExit(
            f"the {SEALED_SPLIT} split is sealed: calibrate on {CALIBRATION_SPLIT} "
            f"and confirm on {CONFIRMATION_SPLIT}"
        )
    if split not in (CALIBRATION_SPLIT, CONFIRMATION_SPLIT):
        raise SystemExit(
            f"calibration uses {CALIBRATION_SPLIT} (calibrate) and {CONFIRMATION_SPLIT} "
            f"(confirm), not {split!r}"
        )


def calibration_dir(out_root: Path, split: str, name: str) -> Path:
    """``<out_root>/calibration/<split>/<name>`` (``name`` must be a plain directory name)."""
    if not name or Path(name).name != name or name in (".", ".."):
        raise SystemExit(f"--name must be a plain directory name, got {name!r}")
    return out_root / CALIBRATION_DIR / split / name


def dump(data: Any) -> str:
    """The JSON text written for every calibration file (stable key order)."""
    return json.dumps(data, indent=JSON_INDENT, sort_keys=True, allow_nan=False) + "\n"


def write_atomic(path: Path, text: str) -> None:
    """Replace ``path`` in one step, so an interrupted run never leaves half a file."""
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)


def file_sha256(path: Path) -> str:
    """SHA-256 of a file's bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest_rows(manifest: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in manifest.read_text().splitlines() if line.strip()]


def _inputs(
    request: CalibrationRequest, grid: dict[str, list[Any]], base: dict[str, Any], manifest: Path
) -> dict[str, Any]:
    rows = _manifest_rows(manifest)
    seeds = {
        row["design_id"]: {k: (row.get("design") or {}).get(k) for k in ("bio_seed", "read_seed")}
        for row in rows
    }
    muconeup = sorted({str(row["muconeup_version"]) for row in rows if row.get("muconeup_version")})
    return {
        "split": request.split,
        "name": request.name,
        "engine": request.engine,
        "grid": grid,
        "grid_sha256": canonical_sha256(grid),
        "base_config": str(request.base_config.resolve()) if request.base_config else None,
        "base_config_sha256": canonical_sha256(base),
        "manifest": str(manifest),
        "manifest_sha256": file_sha256(manifest),
        "design_seeds": seeds,
        "versions": {
            "muconespan": __version__,
            "muconeup": muconeup,
            "python": sys.version.split()[0],
        },
    }


def _point_entry(point: GridPoint) -> dict[str, Any]:
    return {
        "sha256": point.sha256,
        "values": point.values,
        "status": PENDING,
        "evaluate_exit": None,
        "error": None,
    }


def _load_or_create(cal_dir: Path, inputs: dict[str, Any], points: list[GridPoint]) -> dict:
    path = cal_dir / MANIFEST_FILE
    if path.is_file():
        saved: dict[str, Any] = json.loads(path.read_text())
        changed = sorted(k for k in inputs if saved["inputs"].get(k) != inputs[k])
        if changed:
            raise SystemExit(
                f"calibration {cal_dir} was made with different {', '.join(changed)}; "
                "use a new --name"
            )
        return saved
    cal_dir.mkdir(parents=True, exist_ok=True)
    data = {
        "schema_version": SCHEMA_VERSION,
        "inputs": inputs,
        "points": [_point_entry(p) for p in points],
    }
    write_atomic(path, dump(data))
    return data


def _run_point(
    point: GridPoint,
    point_dir: Path,
    request: CalibrationRequest,
    manifest: Path,
    run_fn: RunFn,
    evaluate_fn: EvaluateFn,
) -> dict[str, Any]:
    config_path = point_dir / OVERLAY_FILE
    text = dump(point.config)
    if config_path.is_file() and config_path.read_text() != text:
        raise SystemExit(f"{config_path} does not match its content address; remove it")
    point_dir.mkdir(parents=True, exist_ok=True)
    write_atomic(config_path, text)
    results = point_dir / RESULTS_DIR
    engine_dir = results / request.engine
    try:
        run_fn(manifest, results, config_path)
        ns = argparse.Namespace(
            result_root=engine_dir,
            truth_root=manifest.parent,
            expected_samples=engine_dir / f"inventory_{request.split}.json",
        )
        report, code = evaluate_fn(ns)
        write_atomic(engine_dir / EVALUATION_FILE, dump(report))
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        return {"status": FAILED, "evaluate_exit": None, "error": f"{type(exc).__name__}: {exc}"}
    return {"status": EVALUATED, "evaluate_exit": code, "error": None}


def run_calibration(request: CalibrationRequest, run_fn: RunFn, evaluate_fn: EvaluateFn) -> int:
    """Validate every point, then run and score the ones not yet evaluated.

    Returns 0 when every point is evaluated with a zero evaluator exit, else 1.
    """
    check_split(request.split)
    cal_dir = calibration_dir(request.out_root, request.split, request.name)
    manifest = (request.out_root / request.split / "manifest.jsonl").resolve()
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")
    try:
        grid = load_grid(request.grid)
        base = base_settings(request.base_config)
        points = build_points(grid, base, request.engine)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"invalid calibration grid or base config: {exc}") from exc
    data = _load_or_create(cal_dir, _inputs(request, grid, base, manifest), points)
    entries = {entry["sha256"]: entry for entry in data["points"]}
    for point in points:
        entry = entries[point.sha256]
        point_dir = cal_dir / point.sha256
        evaluation = point_dir / RESULTS_DIR / request.engine / EVALUATION_FILE
        if entry["status"] == EVALUATED and evaluation.is_file():
            print(f"{point.sha256}: reused ({json.dumps(point.values, sort_keys=True)})")
            continue
        entry.update(_run_point(point, point_dir, request, manifest, run_fn, evaluate_fn))
        write_atomic(cal_dir / MANIFEST_FILE, dump(data))
        print(f"{point.sha256}: {entry['status']} ({json.dumps(point.values, sort_keys=True)})")
    ok = all(e["status"] == EVALUATED and e["evaluate_exit"] == 0 for e in data["points"])
    print(f"calibration: {cal_dir / MANIFEST_FILE}")
    return 0 if ok else 1
