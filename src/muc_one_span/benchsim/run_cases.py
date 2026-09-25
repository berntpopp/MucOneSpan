"""Run caller engines over a benchsim split's manifest, keeping every denominator.

Consumes ``<out_root>/<split>/manifest.jsonl`` and the case layout it describes
(``<out_root>/<split>/<design_id>/{truth/, reads/, case.json}``, see
``muc_one_span.benchsim.generate``). Manifest rows carry at least
``design_id``, ``profile`` and ``status`` (``ok``, ``design_invalid`` or
``generation_failed``); only ``ok`` rows have a reads FASTQ to run.

Produces, per engine, under ``<results_root>/<engine>/``:

- ``<design_id>/``: the caller's own output directory (``run_pipeline``'s
  ``output_dir``).
- ``measurements.json``: the run records for every manifest row, in the array
  format ``scripts/evaluate.py`` already reads (keyed by ``sample``).
- ``inventory_<split>.json``: one entry per manifest row with ``sample``,
  ``truth_dir`` and ``result_dir``, for ``scripts/evaluate.py
  --expected-samples``. It is written per engine, since ``result_dir`` differs
  by engine; ``truth_dir`` is the same case truth across engines.

Non-``ok`` manifest rows (``design_invalid``, ``generation_failed``, or any
other non-``ok`` status) are never run; they become ``not_attempted`` records
so the denominator is kept without invoking the caller on invalid input. A
caller that raises once invoked becomes an ``execution_failed`` record (it was
attempted); a missing platform model stays ``not_attempted``.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from muc_one_span.benchmarking import run_pipeline
from muc_one_span.benchsim.generate import FASTQ


def platform_for_profile(profile: Any) -> str:
    """Map a benchsim profile to a caller platform: ``ont_*`` -> ont, ``hifi_amplicon`` -> hifi."""
    if not isinstance(profile, str):
        raise TypeError(f"profile must be a string, got {profile!r}")
    if profile.startswith("ont_"):
        return "ont"
    if profile == "hifi_amplicon":
        return "hifi"
    raise ValueError(f"no platform mapping for profile {profile!r}")


def _case_dirs(split_dir: Path, design_id: str) -> tuple[Path, Path]:
    case_dir = split_dir / design_id
    return case_dir / "truth", case_dir / "reads"


def _read_manifest(path: Path) -> list[dict[str, Any]]:
    lines = path.read_text().splitlines()
    return [json.loads(line) for line in lines if line.strip()]


def _not_attempted(
    sample: str,
    engine: str,
    result_dir: Path,
    error: str,
    status: str = "not_attempted",
    **extra: Any,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "sample": sample,
        "status": status,
        "engine": engine,
        "result_dir": str(result_dir),
        "exit_code": None,
        "error": error,
    }
    record.update({k: v for k, v in extra.items() if v is not None})
    return record


def _run_one(
    row: dict[str, Any],
    split_dir: Path,
    engine: str,
    results_root: Path,
    model_for: Callable[[str], str],
    threads: int,
    config: Path | None = None,
) -> dict[str, Any]:
    """Run one manifest row on one engine, or record why it was not attempted."""
    design_id = row["design_id"]
    profile = row.get("profile")
    status = row.get("status")
    output_dir = results_root / engine / design_id
    try:
        platform = platform_for_profile(profile)
    except (TypeError, ValueError) as exc:
        return _not_attempted(design_id, engine, output_dir, str(exc))
    assert isinstance(profile, str)  # platform_for_profile already rejected non-strings
    if status != "ok":
        reason = row.get("error") or f"manifest status: {status}"
        return _not_attempted(
            design_id, engine, output_dir, reason, platform=platform, profile=profile
        )
    _, reads_dir = _case_dirs(split_dir, design_id)
    pattern = FASTQ.get(profile)
    if pattern is None:
        return _not_attempted(
            design_id,
            engine,
            output_dir,
            f"no reads filename pattern for profile {profile!r}",
            platform=platform,
            profile=profile,
        )
    fastq = reads_dir / pattern.format(design_id)
    if not fastq.is_file():
        return _not_attempted(
            design_id,
            engine,
            output_dir,
            f"reads file not found: {fastq}",
            platform=platform,
            profile=profile,
        )
    try:
        model = model_for(platform)
    except (ValueError, KeyError) as exc:
        return _not_attempted(
            design_id, engine, output_dir, str(exc), platform=platform, profile=profile
        )
    try:
        if config is None:  # keep the call unchanged for runs without a config
            record = run_pipeline(
                design_id, fastq, output_dir, platform, model, threads, engine=engine
            )
        else:
            record = run_pipeline(
                design_id, fastq, output_dir, platform, model, threads, engine=engine, config=config
            )
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        return _not_attempted(
            design_id,
            engine,
            output_dir,
            str(exc),
            status="execution_failed",
            platform=platform,
            profile=profile,
        )
    record = dict(record)
    record["engine"] = engine
    record.setdefault("profile", profile)
    return record


def _write_engine_outputs(
    records: list[dict[str, Any]],
    engines: Sequence[str],
    rows: list[dict[str, Any]],
    split_dir: Path,
    split: str,
    results_root: Path,
) -> None:
    by_engine: dict[str, list[dict[str, Any]]] = {engine: [] for engine in engines}
    for record in records:
        by_engine.setdefault(record["engine"], []).append(record)
    for engine in engines:
        engine_dir = results_root / engine
        engine_dir.mkdir(parents=True, exist_ok=True)
        (engine_dir / "measurements.json").write_text(
            json.dumps(by_engine.get(engine, []), indent=2) + "\n"
        )
        inventory = [
            {
                "sample": row["design_id"],
                "truth_dir": str(_case_dirs(split_dir, row["design_id"])[0]),
                "result_dir": str(engine_dir / row["design_id"]),
            }
            for row in rows
        ]
        (engine_dir / f"inventory_{split}.json").write_text(json.dumps(inventory, indent=2) + "\n")


def run_split(
    manifest: Path,
    engines: Sequence[str],
    results_root: Path,
    model_for: Callable[[str], str],
    threads: int,
    jobs: int,
    config: Path | None = None,
) -> list[dict[str, Any]]:
    """Run every engine over every manifest row without dropping the denominator.

    Each (row, engine) pair is independent; ``jobs`` selects a process pool
    when greater than 1. Non-``ok`` rows never reach the caller. ``config`` is
    a runtime settings JSON forwarded to every run (``muconespan --config``).
    """
    manifest = Path(manifest).resolve()
    split_dir = manifest.parent
    split = split_dir.name
    rows = _read_manifest(manifest)
    results_root = Path(results_root)
    engines = list(engines)
    pairs = [(row, engine) for engine in engines for row in rows]
    if jobs <= 1:
        records = [
            _run_one(row, split_dir, engine, results_root, model_for, threads, config)
            for row, engine in pairs
        ]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            records = list(
                pool.map(
                    _run_one,
                    [pair[0] for pair in pairs],
                    [split_dir] * len(pairs),
                    [pair[1] for pair in pairs],
                    [results_root] * len(pairs),
                    [model_for] * len(pairs),
                    [threads] * len(pairs),
                    [config] * len(pairs),
                )
            )
    _write_engine_outputs(records, engines, rows, split_dir, split, results_root)
    return records
