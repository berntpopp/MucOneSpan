#!/usr/bin/env python3
"""Generate inventoried HiFi/ONT experiments from a strict JSON design."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from muc_one_span.experiments import run_experiment


def parser() -> argparse.ArgumentParser:
    """Build the simulator experiment CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--design",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "examples/development-experiment.json",
    )
    result.add_argument(
        "--config", type=Path, help="MucOneUp config for all platforms (overrides design)"
    )
    result.add_argument(
        "--muconeup",
        default="muconeup",
        help="MucOneUp executable; absolute path recommended for virtualenv installs",
    )
    result.add_argument("--output-dir", type=Path, default=Path("tests/data/generated"))
    result.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands; create no files and execute no tools",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    """Run every case and return nonzero if any generation fails."""
    args = parser().parse_args(argv)
    fallback = os.environ.get("MUCONEUP_CONFIG")
    try:
        manifest = run_experiment(
            args.design,
            args.output_dir,
            config=args.config,
            executable=args.muconeup,
            dry_run=args.dry_run,
            fallback_config=Path(fallback) if fallback else None,
        )
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    if args.dry_run:
        print(json.dumps(manifest, indent=2, allow_nan=False))
    else:
        for case in manifest["cases"]:
            print(case["sample"], case["status"], case.get("error", ""))
        print(f"Manifest: {args.output_dir / 'generation_manifest.json'}")
    return int(any(c["status"] not in {"completed", "planned"} for c in manifest["cases"]))


if __name__ == "__main__":
    raise SystemExit(main())
