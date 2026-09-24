"""Auditable in-process runner and strict evaluation support for benchmark scripts."""

from __future__ import annotations

import hashlib
import json
import time
from collections.abc import Callable, Sequence
from contextlib import ExitStack
from pathlib import Path
from typing import Any, Protocol
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_span.config import RepeatDictionary, load_repeat_dictionary
from muc_one_span.evaluation import (
    TruthValidationError,
    aggregate,
    evaluate_sample,
    load_observation,
    load_truth,
)
from muc_one_span.evaluation.artifacts import discover_input, read_inventory


class Runner(Protocol):
    """Small Click runner protocol used to isolate real pipeline execution in tests."""

    def invoke(self, command: Any, args: list[str]) -> Any:
        """Invoke a Click command and return its result."""


def read_sample_inventory(data_root: Path, inventory: Path | None = None) -> list[dict[str, Any]]:
    """Read an explicit inventory or include every sample directory in sorted order."""
    if not data_root.is_dir():
        raise ValueError(f"data root is not an existing directory: {data_root}")
    if inventory is not None:
        return read_inventory(inventory)
    rows = [{"sample": path.name} for path in sorted(data_root.iterdir()) if path.is_dir()]
    if not rows:
        raise ValueError(f"no sample directories in {data_root}")
    return rows


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _timed(
    label: str, function: Callable[..., Any], timings: dict[str, float]
) -> Callable[..., Any]:
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        started = time.perf_counter()
        try:
            return function(*args, **kwargs)
        finally:
            timings[label] = timings.get(label, 0.0) + time.perf_counter() - started

    return wrapper


def _stage_patches(timings: dict[str, float]) -> list[Any]:
    from muc_one_span import alleles, calling, classify, consensus, mapping

    stages = (
        (mapping, "map_reads", "mapping"),
        (alleles, "detect_alleles", "alleles"),
        (calling, "call_variants_per_allele", "calling"),
        (consensus, "build_consensus_per_allele", "consensus"),
        (classify, "classify_sequence", "classify"),
    )
    return [
        patch.object(module, name, _timed(label, getattr(module, name), timings))
        for module, name, label in stages
    ]


def _status(output_dir: Path, exit_code: int) -> tuple[str, str | None]:
    sidecar = output_dir / "run_status.json"
    if sidecar.exists():
        try:
            data = json.loads(sidecar.read_text())
            status = data.get("status")
            if status in ("execution_failed", "insufficient_evidence"):
                return str(status), data.get("error")
            if status == "completed" and exit_code == 0:
                return "completed", None
        except (OSError, ValueError, TypeError):
            pass
    if exit_code:
        return "execution_failed", "pipeline command failed"
    if not (output_dir / "summary.json").is_file():
        return "invalid_artifacts", "successful command did not produce summary.json"
    return "completed", None


def run_pipeline(
    sample: str,
    input_path: Path,
    output_dir: Path,
    platform: str,
    model: str,
    threads: int,
    *,
    runner: Runner | None = None,
    engine: str = "ladder",
) -> dict[str, Any]:
    """Run the real full CLI while timing its five scientific stages.

    ``engine`` is appended to ``cli_args`` as ``--engine <engine>`` only when it
    is not the default ``"ladder"`` (the current CLI has no ``--engine`` flag;
    it arrives with the hybrid engine), and is always recorded in the returned
    record and ``measurement.json``.
    """
    from muc_one_span.cli import main

    output_dir.mkdir(parents=True, exist_ok=True)
    cli_args = [
        "run",
        "--input",
        str(input_path.resolve()),
        "--output-dir",
        str(output_dir.resolve()),
        "--clair3-model",
        model,
        "--threads",
        str(threads),
        "--platform",
        platform,
    ]
    if engine != "ladder":
        cli_args += ["--engine", engine]
    timings: dict[str, float] = {}
    started = time.perf_counter()
    with ExitStack() as stack:
        for stage_patch in _stage_patches(timings):
            stack.enter_context(stage_patch)
        result = (runner or CliRunner()).invoke(main, cli_args)
    timings["wall"] = time.perf_counter() - started
    status, status_error = _status(output_dir, int(result.exit_code))
    exception = str(result.exception) if result.exception is not None else None
    error = status_error or exception
    log = result.output
    if exception:
        log += f"\n{type(result.exception).__name__}: {exception}\n"
    (output_dir / "cli.log").write_text(log)
    record: dict[str, Any] = {
        "sample": sample,
        "status": status,
        "input": str(input_path),
        "input_sha256": _hash(input_path),
        "result_dir": str(output_dir),
        "platform": platform,
        "model": model,
        "threads": threads,
        "cli_args": cli_args,
        "exit_code": int(result.exit_code),
        "timings": timings,
        "engine": engine,
    }
    if error:
        record["error"] = error
    (output_dir / "measurement.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def _setting(
    entry: dict[str, Any], key: str, command_value: Any, fallback: Any, allowed: type
) -> Any:
    entry_value = entry.get(key)
    if entry_value is not None and not isinstance(entry_value, allowed):
        raise ValueError(f"{entry['sample']}: {key} has an invalid type")
    if command_value is not None and entry_value is not None and command_value != entry_value:
        raise ValueError(f"{entry['sample']}: inventory {key} conflicts with command line")
    if command_value is not None:
        return command_value
    return entry_value if entry_value is not None else fallback


def _metadata_platform(sample_dir: Path) -> str | None:
    platforms: set[str] = set()
    for path in sorted(sample_dir.glob("*_metadata.tsv")):
        fields = dict(line.split("\t", 1) for line in path.read_text().splitlines() if "\t" in line)
        technology = fields.get("Read_simulation_technology", "").strip().lower()
        if technology in ("ont", "nanopore"):
            platforms.add("ont")
        elif technology in ("pacbio", "hifi", "pacbio_hifi"):
            platforms.add("hifi")
    if len(platforms) > 1:
        raise ValueError(f"{sample_dir.name}: conflicting platform metadata")
    return next(iter(platforms), None)


def platform_for_sample(
    sample_dir: Path, entry: dict[str, Any], command_platform: str | None
) -> str:
    """Resolve platform controls and reject conflicts with simulator provenance."""
    chosen = _setting(entry, "platform", command_platform, None, str)
    metadata = _metadata_platform(sample_dir)
    if chosen is not None and metadata is not None and chosen != metadata:
        raise ValueError(
            f"{entry.get('sample', sample_dir.name)}: platform {chosen} conflicts with metadata {metadata}"
        )
    result = chosen or metadata or "hifi"
    if result not in ("hifi", "ont"):
        raise ValueError(f"{entry.get('sample', sample_dir.name)}: unsupported platform {result!r}")
    return result


def run_inventory(
    entries: Sequence[dict[str, Any]],
    data_root: Path,
    output_root: Path,
    *,
    platform: str | None = None,
    model: str | None = None,
    threads: int | None = None,
) -> list[dict[str, Any]]:
    """Run or explicitly fail every inventory entry without dropping denominators."""
    records: list[dict[str, Any]] = []
    platforms: set[str] = set()
    if model:
        for entry in entries:
            try:
                platforms.add(
                    platform_for_sample(
                        Path(entry.get("input_dir", data_root / entry["sample"])), entry, platform
                    )
                )
            except (OSError, ValueError):
                # Preserve the per-entry error below, without executing invalid input.
                continue
    shared_model_error = bool(model and len(platforms) > 1)
    output_root.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        sample = entry["sample"]
        sample_dir = Path(entry.get("input_dir", data_root / sample))
        result_dir = Path(entry.get("result_dir", output_root / sample))
        try:
            if shared_model_error:
                raise ValueError(
                    "Mixed HiFi/ONT inventory requires explicit per-sample models; omit the shared model"
                )
            selected = entry.get("input")
            if selected is not None and not isinstance(selected, str):
                raise ValueError(f"{sample}: input must be a filename")
            input_path = discover_input(sample_dir, selected)
            run_platform = platform_for_sample(sample_dir, entry, platform)
            run_model = _setting(entry, "model", model, "", str)
            run_threads = _setting(entry, "threads", threads, 4, int)
            if isinstance(run_threads, bool) or run_threads < 1:
                raise ValueError(f"{sample}: threads must be a positive integer")
            record = run_pipeline(
                sample, input_path, result_dir, run_platform, run_model, run_threads
            )
        except (OSError, ValueError) as exc:
            record = {
                "sample": sample,
                "status": "not_attempted",
                "result_dir": str(result_dir),
                "exit_code": None,
                "error": str(exc),
            }
            for key in ("platform", "model", "threads", "input"):
                if key in entry:
                    record[key] = entry[key]
        records.append(record)
        (output_root / "measurements.json").write_text(json.dumps(records, indent=2) + "\n")
    return records


def evaluate_inventory(
    entries: Sequence[dict[str, Any]],
    truth_root: Path,
    result_root: Path,
    records: Sequence[dict[str, Any]],
    repeat_dict: RepeatDictionary | None = None,
    inventory_mode: str = "explicit",
) -> tuple[dict[str, Any], bool]:
    """Score batch artifacts with the same strict APIs as scripts/evaluate.py."""
    by_sample = {record["sample"]: record for record in records}
    rd = repeat_dict or load_repeat_dictionary()
    rows: list[dict[str, Any]] = []
    failed = False
    for entry in entries:
        sample = entry["sample"]
        truth_dir = Path(entry.get("truth_dir", truth_root / entry.get("truth_sample", sample)))
        result_dir = Path(entry.get("result_dir", result_root / sample))
        record = by_sample.get(sample, {"sample": sample, "status": "not_attempted"})
        observation = load_observation(result_dir, record)
        try:
            truth = load_truth(truth_dir, rd)
            if truth.name != sample:
                from dataclasses import replace

                truth = replace(truth, name=sample)
            row = evaluate_sample(truth, observation)
        except TruthValidationError as exc:
            row = {
                "sample": sample,
                "status": "invalid_truth",
                "truth_status": "invalid",
                "observation_status": observation.status,
                "error": str(exc),
                "run_record": observation.run_record,
            }
        row["inventory"] = entry
        row["truth_dir"] = str(truth_dir)
        row["result_dir"] = str(result_dir)
        rows.append(row)
        failed |= record.get("exit_code") not in (None, 0)
        failed |= row["status"] in (
            "invalid_truth",
            "invalid_artifacts",
            "execution_failed",
            "not_attempted",
        )
    report = aggregate(rows)
    report["inventory_mode"] = inventory_mode
    report["truth_root"] = str(truth_root)
    report["result_root"] = str(result_root)
    return report, failed
