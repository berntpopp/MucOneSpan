#!/usr/bin/env python3
"""Run a sample batch and score full artifacts with the strict offline evaluator."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from muc_one_span.benchmarking import (
    evaluate_inventory,
    read_sample_inventory,
    run_inventory,
)
from muc_one_span.evaluation import aggregate

SAMPLES_DIR = Path("tests/data/generated")
OUTPUT_DIR = Path("tests/results")

SCHEMA_MIGRATION = {
    "from": "unversioned array of TP/TP_partial/FN/FP/TN/error rows",
    "to": "evaluation schema_version 1",
    "removed_status": "TP_partial",
    "replacement": "samples plus totals.metrics.event_tp/event_fn/event_fp and explicit denominators",
}


def parser() -> argparse.ArgumentParser:
    """Keep legacy positionals and add explicit execution controls."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("samples_dir", nargs="?", type=Path, default=SAMPLES_DIR)
    result.add_argument("output_dir", nargs="?", type=Path, default=OUTPUT_DIR)
    result.add_argument("--platform", choices=("hifi", "ont"), default=None)
    result.add_argument(
        "--clair3-model", default=None, help="Clair3 model (defaults to CLAIR3_MODEL)"
    )
    result.add_argument("--threads", type=int, default=None, help="Threads for every sample")
    result.add_argument(
        "--expected-samples", type=Path, help="JSON inventory; every entry stays in denominators"
    )
    return result


def _write_error(output_dir: Path, error: Exception) -> None:
    report = aggregate([])
    report["input_error"] = str(error)
    report["schema_migration"] = SCHEMA_MIGRATION
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "batch_results.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n"
        )
    except OSError:
        pass


def main(argv: list[str] | None = None) -> int:
    """Run and strictly evaluate every expected sample."""
    args = parser().parse_args(argv)
    try:
        entries = read_sample_inventory(args.samples_dir, args.expected_samples)
        records = run_inventory(
            entries,
            args.samples_dir,
            args.output_dir,
            platform=args.platform,
            model=args.clair3_model if args.clair3_model is not None else os.getenv("CLAIR3_MODEL"),
            threads=args.threads,
        )
        report, failed = evaluate_inventory(
            entries,
            args.samples_dir,
            args.output_dir,
            records,
            inventory_mode="explicit" if args.expected_samples else "discovered",
        )
    except (OSError, ValueError, KeyError, TypeError) as exc:
        _write_error(args.output_dir, exc)
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    report["schema_migration"] = SCHEMA_MIGRATION
    path = args.output_dir / "batch_results.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report["totals"], indent=2, allow_nan=False))
    print(f"Full results: {path}")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
