#!/usr/bin/env python3
"""Evaluate caller outputs offline against strict MucOneUp truth and an inventory.

Inventory is a JSON array of sample IDs or objects. Objects may specify
truth_sample (for perturbations), truth_dir, result_dir, input, platform and
run_record. All additional inventory provenance is retained verbatim.
Historical outputs without run records have unknown execution provenance.
Accuracy of zero is valid; missing/failed/invalid inputs return nonzero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation import (
    TruthValidationError,
    aggregate,
    evaluate_sample,
    load_observation,
    load_truth,
)
from muc_one_span.evaluation.artifacts import discover_input, read_inventory


def _run_records(root: Path) -> dict[str, dict[str, Any]]:
    path = root / "measurements.json"
    if not path.exists():
        return {}
    raw = json.loads(path.read_text())
    if not isinstance(raw, list):
        raise ValueError("measurements.json must contain an array")
    records = {}
    for row in raw:
        name = row["sample"]
        if name in records:
            raise ValueError(f"duplicate run record: {name}")
        records[name] = row
    return records


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(args: argparse.Namespace) -> tuple[dict[str, Any], int]:
    """Evaluate every explicit sample and preserve per-sample adapter errors."""
    inventory = read_inventory(args.expected_samples)
    if not args.result_root.is_dir() or not args.truth_root.is_dir():
        raise ValueError("result root and truth root must be existing directories")
    runs = _run_records(args.result_root)
    rd = load_repeat_dictionary()
    rows = []
    failed = False
    for entry in inventory:
        name = entry["sample"]
        truth_dir = Path(entry.get("truth_dir", args.truth_root / entry.get("truth_sample", name)))
        result_dir = Path(entry.get("result_dir", args.result_root / name))
        run_record = entry.get("run_record", runs.get(name))
        if (
            entry.get("run_record") is not None
            and name in runs
            and entry["run_record"] != runs[name]
        ):
            raise ValueError(f"conflicting run records: {name}")
        observation = load_observation(result_dir, run_record)
        try:
            truth = replace(load_truth(truth_dir, rd), name=name)
            row = evaluate_sample(truth, observation)
        except TruthValidationError as exc:
            row = {
                "sample": name,
                "status": "invalid_truth",
                "truth_status": "invalid",
                "observation_status": observation.status,
                "error": str(exc),
                "run_record": observation.run_record,
            }
        row["inventory"] = entry
        row["truth_dir"] = str(truth_dir)
        row["result_dir"] = str(result_dir)
        if entry.get("input"):
            try:
                path = discover_input(Path(entry.get("input_dir", truth_dir)), entry["input"])
                row["input_sha256"] = _hash(path)
                row["input_path"] = str(path)
            except (ValueError, OSError) as exc:
                row["input_error"] = str(exc)
                failed = True
        rows.append(row)
        failed |= row["status"] in (
            "invalid_truth",
            "invalid_artifacts",
            "execution_failed",
            "not_attempted",
        )
    report = aggregate(rows)
    report["inventory_mode"] = "explicit"
    report["inventory_sha256"] = _hash(args.expected_samples)
    report["result_root"] = str(args.result_root)
    report["truth_root"] = str(args.truth_root)
    return report, int(failed)


def main(argv: list[str] | None = None) -> int:
    """Parse CLI options, always write valid JSON when an output path is usable."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_root", type=Path)
    parser.add_argument(
        "--truth-root", "--truth-dir", type=Path, default=Path("tests/data/generated")
    )
    parser.add_argument(
        "--expected-samples",
        type=Path,
        required=True,
        help="Explicit JSON inventory of all expected samples, including failures",
    )
    parser.add_argument("--output", type=Path, default=Path("results/evaluation/report.json"))
    args = parser.parse_args(argv)
    try:
        report, code = run(args)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        report, code = aggregate([]), 2
        report["input_error"] = str(exc)
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except OSError as exc:
        parser.exit(2, f"cannot write evaluation report: {exc}\n")
    print(json.dumps(report["totals"], indent=2, allow_nan=False))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
