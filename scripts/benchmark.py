#!/usr/bin/env python3
"""Benchmark the full pipeline while retaining stage timings and artifacts."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from muc_one_span.benchmarking import read_sample_inventory, run_inventory


def parser() -> argparse.ArgumentParser:
    """Build the maintained benchmark command parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--data-dir", type=Path, default=Path("tests/data/generated"), help="Sample root"
    )
    result.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/results/benchmarks"),
        help="Benchmark output root",
    )
    result.add_argument(
        "--expected-samples", type=Path, help="JSON inventory with explicit sample inputs"
    )
    result.add_argument("--platform", choices=("hifi", "ont"), help="Platform for all samples")
    result.add_argument(
        "--clair3-model", default=None, help="Clair3 model (defaults to CLAIR3_MODEL)"
    )
    result.add_argument("--threads", type=int, default=None, help="Threads for every sample")
    return result


def main(argv: list[str] | None = None) -> int:
    """Run inventoried samples and return nonzero if any sample did not complete."""
    args = parser().parse_args(argv)
    try:
        entries = read_sample_inventory(args.data_dir, args.expected_samples)
        records = run_inventory(
            entries,
            args.data_dir,
            args.output_dir,
            platform=args.platform,
            model=args.clair3_model if args.clair3_model is not None else os.getenv("CLAIR3_MODEL"),
            threads=args.threads,
        )
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    path = args.output_dir / "benchmark_results.json"
    path.write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")
    for record in records:
        timing = {key: round(value, 2) for key, value in record.get("timings", {}).items()}
        print(record["sample"], record["status"], timing)
    print(f"Results written to {path}")
    return int(any(record["status"] != "completed" for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
