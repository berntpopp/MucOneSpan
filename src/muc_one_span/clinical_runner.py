"""Provenance-checked clinical runs without changing caller semantics."""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

from muc_one_span.clinical_provenance import now, object_hash, sha256_file, write_json
from muc_one_span.evaluation.artifacts import load_observation
from muc_one_span.tools import run_tool_pipeline

# The engine of hashed settings that carry no ``engine`` key (attempts made before engine
# selection existed, and ladder attempts, which keep that hash so they stay resumable).
LEGACY_UNHASHED_ENGINE = "ladder"


def _evidence(root: Path, exit_code: int | None) -> dict[str, Any]:
    observation = load_observation(root, {"exit_code": exit_code})
    events = [
        {
            "event": event.name,
            "allele": allele.name,
            "repeat_index": event.repeat_index,
            "parent": event.parent,
            "supported": event.supported,
            "support_status": event.support_status,
        }
        for allele in observation.predictions
        for event in allele.events
    ]
    return {
        "analysis_state": observation.status,
        "events": events if observation.predictions else None,
        "sequences": [a.sequence for a in observation.predictions],
        "sequence_independent": [a.independent_haplotype_evidence for a in observation.predictions],
        "sequence_sources": [a.sequence_source for a in observation.predictions],
        "warnings": list(observation.warnings),
    }


def validate_record(
    record: dict[str, Any], *, expected_provenance: dict[str, Any] | None = None
) -> bool:
    """Re-read bound artifacts before accepting a record for scoring or resume."""
    try:
        provenance = dict(record["provenance"])
        digest = provenance.pop("sha256")
        if digest != object_hash(provenance):
            return False
        if any(
            record.get(key) != provenance["run"].get(key)
            for key in ("run_accession", "arm", "biological_sample")
        ):
            return False
        if expected_provenance is not None and record["provenance"] != expected_provenance:
            return False
        root = Path(record["result_dir"]).resolve()
        outputs = record["outputs"]
        if not outputs or "worker.json" not in outputs:
            return False
        actual_files = {
            str(p.relative_to(root))
            for p in root.rglob("*")
            if p.is_file() and p.name != "benchmark_record.json"
        }
        if set(outputs) != actual_files:
            return False
        preparation = provenance["preparation"]
        if sha256_file(Path(preparation["output_path"])) != preparation["output_sha256"]:
            return False
        for name, metadata in outputs.items():
            path = Path(metadata["path"]).resolve()
            if path != (root / name).resolve() or not path.is_relative_to(root):
                return False
            if sha256_file(path) != metadata["sha256"]:
                return False
        worker = json.loads((root / "worker.json").read_text())
        if worker["exit_code"] != record["exit_code"]:
            return False
        if record["status"] == "completed":
            if (
                record["exit_code"] != 0
                or not {"summary.json", "run_status.json", "run_configuration.json"}
                <= outputs.keys()
            ):
                return False
            sidecar = json.loads((root / "run_status.json").read_text())
            if sidecar.get("status") != "completed":
                return False
            configuration = json.loads((root / "run_configuration.json").read_text())
            if configuration.get("input_sha256") != preparation["output_sha256"]:
                return False
        evidence = _evidence(root, record["exit_code"])
        return not any(record.get(key) != value for key, value in evidence.items())
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return False


def run_case(
    run: dict[str, Any],
    preparation: dict[str, Any],
    output_root: Path,
    settings: dict[str, Any],
    *,
    resume: bool = False,
) -> dict[str, Any]:
    """Run one accession in a new directory; resume only identical validated attempts."""
    input_path = Path(preparation["output_path"]).resolve()
    if sha256_file(input_path) != preparation["output_sha256"]:
        raise ValueError("prepared input hash mismatch")
    accession = run["run_accession"]
    if not isinstance(accession, str) or not accession.isalnum():
        raise ValueError("invalid run accession")
    if preparation.get("run_accession") != accession:
        raise ValueError("prepared input accession mismatch")
    expected_files = [
        {key: item[key] for key in ("url", "bytes", "md5")} for item in run.get("files", [])
    ]
    prepared_files = [
        {key: item[key] for key in ("url", "bytes", "md5")} for item in preparation.get("files", [])
    ]
    if expected_files != prepared_files:
        raise ValueError("prepared input source inventory mismatch")
    root = output_root.resolve() / accession
    provenance: dict[str, Any] = {"run": run, "preparation": preparation, "settings": settings}
    provenance["sha256"] = object_hash(provenance)
    if root.exists():
        if resume and (root / "benchmark_record.json").is_file():
            old = json.loads((root / "benchmark_record.json").read_text())
            if validate_record(old, expected_provenance=provenance):
                return dict(old)
        raise ValueError(f"existing attempt cannot be reused: {root}; use a new output root")
    timeout = settings["timeout"]
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout)
        or timeout <= 0
    ):
        raise ValueError("timeout must be finite and positive")
    root.mkdir(parents=True)
    argv = [
        "run",
        "--input",
        str(input_path),
        "--output-dir",
        str(root),
        "--platform",
        "ont",
        "--clair3-model",
        settings["model"],
        "--threads",
        str(settings["threads"]),
        "--report",
        "--report-igv",
        "off",
    ]
    # Hashed settings: --resume cannot mix engines. Settings without an engine key were
    # hashed before engine selection existed and ran the ladder, which is no longer the
    # CLI default, so the engine is always passed explicitly.
    argv += ["--engine", str(settings.get("engine", LEGACY_UNHASHED_ENGINE))]
    if settings.get("assay") is not None:
        argv += ["--assay", str(settings["assay"])]
    invocation = root / "invocation.json"
    write_json(invocation, {"argv": argv, "output": str(root)})
    command = [sys.executable, "-m", "muc_one_span.clinical_worker", str(invocation)]
    record: dict[str, Any] = {
        "schema_version": 1,
        "run_accession": accession,
        "biological_sample": run.get("biological_sample"),
        "arm": run["arm"],
        "status": "running",
        "exit_code": None,
        "started_at": now(),
        "ended_at": None,
        "argv": argv,
        "supervisor_argv": command,
        "provenance": provenance,
        "result_dir": str(root),
        "outputs": {},
        "events": None,
    }
    write_json(root / "benchmark_record.json", record)
    started = time.monotonic()
    error = None
    interrupted = False
    try:
        with (root / "supervisor.log").open("wb") as log:
            run_tool_pipeline([command], timeout=timeout, stdout=log)
    except (RuntimeError, TimeoutError, FileNotFoundError, KeyboardInterrupt) as exc:
        error = f"{type(exc).__name__}: {exc}"
        interrupted = isinstance(exc, KeyboardInterrupt)
    worker_path = root / "worker.json"
    try:
        worker = json.loads(worker_path.read_text()) if worker_path.is_file() else {}
        if not isinstance(worker, dict):
            raise ValueError("worker record must be an object")
    except (OSError, ValueError) as exc:
        worker = {}
        error = f"invalid worker record: {exc}"
    record.update(
        exit_code=worker.get("exit_code"), ended_at=now(), error=error or worker.get("error")
    )
    record["measurements"] = {
        "supervised_wall_seconds": time.monotonic() - started,
        "requested_threads": settings["threads"],
        "timeout_seconds": timeout,
        **worker,
    }
    evidence = _evidence(root, record["exit_code"])
    record.update(evidence)
    state = evidence["analysis_state"]
    record["status"] = (
        "execution_failed"
        if error
        else "insufficient_evidence"
        if state == "insufficient_evidence" and record["exit_code"] == 1
        else "execution_failed"
        if record["exit_code"] != 0
        else "invalid_artifacts"
        if state == "invalid_artifacts"
        else "completed"
    )
    record["callable"] = record["status"] == "completed" and state == "completed"
    record["report_state"] = "completed" if (root / "report.html").is_file() else "not_produced"
    record["outputs"] = {
        str(p.relative_to(root)): {"path": str(p), "sha256": sha256_file(p)}
        for p in sorted(root.rglob("*"))
        if p.is_file() and p.name != "benchmark_record.json"
    }
    write_json(root / "benchmark_record.json", record)
    if interrupted:
        raise KeyboardInterrupt
    return record
