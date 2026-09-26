#!/usr/bin/env python3
"""Inventory, prepare, freeze, run and score the PRJEB92208 clinical benchmark."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

from muc_one_span.clinical_data import fetch_inventory, prepare_run, validate_inventory
from muc_one_span.clinical_provenance import (
    freeze_environment,
    frozen_model_path,
    object_hash,
    sha256_file,
    verify_environment,
    write_json,
)
from muc_one_span.clinical_runner import LEGACY_UNHASHED_ENGINE, run_case
from muc_one_span.settings import DEFAULT_SETTINGS


def read(path: Path) -> dict[str, Any]:
    """Read a JSON object."""
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def load_records(root: Path) -> list[dict[str, Any]]:
    """Combine terminal attempts with latest invocation failures across selections."""
    records = {p.parent.name: read(p) for p in sorted(root.glob("*/benchmark_record.json"))}
    journal = root / "cohort_attempts.json"
    if journal.exists():
        for record in read(journal)["records"]:
            records[record["run_accession"]] = record
    return list(records.values())


def parser() -> argparse.ArgumentParser:
    """Build additive benchmark command without modifying the caller CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="action", required=True)
    inventory = commands.add_parser("inventory")
    inventory.add_argument("--study", default="PRJEB92208")
    inventory.add_argument("--output", type=Path, required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--manifest", type=Path, required=True)
    prepare.add_argument("--data-root", type=Path, required=True)
    prepare.add_argument("--run", action="append", default=[])
    freeze = commands.add_parser("freeze")
    freeze.add_argument("--checkout", type=Path, default=Path.cwd())
    freeze.add_argument(
        "--model", type=Path, default=None, help="Clair3 model (ladder engine only)"
    )
    freeze.add_argument("--output", type=Path, required=True)
    truth = commands.add_parser("truth")
    truth.add_argument("--ledger", type=Path, required=True)
    truth.add_argument("--cache-root", type=Path, required=True)
    truth.add_argument("--output", type=Path, required=True)
    run = commands.add_parser("run")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--data-root", type=Path, required=True)
    run.add_argument("--output-root", type=Path, required=True)
    run.add_argument("--environment", type=Path, required=True)
    run.add_argument("--threads", type=int, default=2)
    run.add_argument("--timeout", type=float, default=3600)
    run.add_argument("--run", action="append", default=[])
    run.add_argument("--resume", action="store_true")
    run.add_argument("--engine", choices=("ladder", "hybrid"), default=DEFAULT_SETTINGS.run.engine)
    run.add_argument("--assay", choices=("amplicon", "genomic"), default=None)
    score = commands.add_parser("score")
    score.add_argument("--manifest", type=Path, required=True)
    score.add_argument("--truth", type=Path, required=True)
    score.add_argument("--output-root", type=Path, required=True)
    score.add_argument("--output", type=Path, required=True)
    return result


def main(argv: list[str] | None = None) -> int:
    """Run selected operation; preserve failures with nonzero command exit."""
    args = parser().parse_args(argv)
    try:
        if args.action == "inventory":
            write_json(args.output, fetch_inventory(args.study))
            return 0
        if args.action == "freeze":
            write_json(args.output, freeze_environment(args.checkout, args.model))
            return 0
        if args.action == "truth":
            from muc_one_span.clinical_truth import hydrate_truth

            write_json(args.output, hydrate_truth(read(args.ledger), args.cache_root))
            return 0
        inventory = validate_inventory(read(args.manifest))
        if args.action == "score":
            from muc_one_span.clinical_scoring import score_cohort

            records = load_records(args.output_root)
            write_json(args.output, score_cohort(inventory, read(args.truth), records))
            return 0
        selected = set(args.run)
        accessions = {run["run_accession"] for run in inventory["runs"]}
        if selected - accessions:
            raise ValueError(f"unknown requested accessions: {sorted(selected - accessions)}")
        runs = [
            run
            for run in inventory["runs"]
            if run["arm"] != "excluded" and (not selected or run["run_accession"] in selected)
        ]
        failures = 0
        journal_root = args.data_root if args.action == "prepare" else args.output_root
        attempt_records = {r["run_accession"]: r for r in load_records(journal_root)}
        for run in runs:
            accession = run["run_accession"]
            preparation_path = args.data_root / accession / "preparation.json"
            print(f"{args.action}: {accession}", flush=True)
            try:
                if args.action == "prepare":
                    started = time.monotonic()
                    preparation = prepare_run(run, args.data_root.resolve())
                    preparation["download_preprocessing_wall_seconds"] = time.monotonic() - started
                    write_json(preparation_path, preparation)
                else:
                    environment = read(args.environment)
                    checkout = Path(__file__).resolve().parents[1]
                    verify_environment(environment, checkout)
                    model = frozen_model_path(environment)
                    if model is None and args.engine == LEGACY_UNHASHED_ENGINE:
                        raise SystemExit(
                            "the ladder engine needs a frozen Clair3 model; freeze with --model"
                        )
                    settings = {
                        "threads": args.threads,
                        "timeout": args.timeout,
                        "model": model,
                        "environment": environment,
                        "environment_sha256": object_hash(environment),
                        "report_igv": "off",
                        # The ladder hashes without an engine key (the pre-0.17 hash),
                        # so existing ladder attempts stay resumable.
                        **(
                            {"engine": args.engine} if args.engine != LEGACY_UNHASHED_ENGINE else {}
                        ),
                        **({"assay": args.assay} if args.assay else {}),
                        "repeat_policy": "one observed execution per library",
                        "harness_sha256": {
                            str(p.relative_to(checkout)): sha256_file(p)
                            for p in sorted((checkout / "src/muc_one_span").glob("clinical_*.py"))
                        },
                    }
                    record = run_case(
                        run, read(preparation_path), args.output_root, settings, resume=args.resume
                    )
                    failures += record["status"] != "completed"
                    print(
                        f"{accession}: {record['status']} ({record['analysis_state']})", flush=True
                    )
                if args.action == "run":
                    attempt_records[accession] = record
            except (OSError, ValueError, KeyError) as exc:
                failures += 1
                failed = {
                    "schema_version": 1,
                    "run_accession": accession,
                    "arm": run["arm"],
                    "biological_sample": run.get("biological_sample"),
                    "status": "execution_failed",
                    "exit_code": None,
                    "analysis_state": "unattempted",
                    "events": None,
                    "error": str(exc),
                }
                attempt_records[accession] = failed
                print(f"{accession}: {exc}", file=sys.stderr, flush=True)
            journal_root = args.data_root if args.action == "prepare" else args.output_root
            journal_root.mkdir(parents=True, exist_ok=True)
            write_json(
                journal_root / "cohort_attempts.json",
                {"schema_version": 1, "records": list(attempt_records.values())},
            )
        return int(failures > 0)
    except (OSError, ValueError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
