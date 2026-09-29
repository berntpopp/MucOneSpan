"""Content hashes and frozen environment for clinical benchmark replay."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from muc_one_span.tools import get_tool_versions, run_tool


def sha256_file(path: Path) -> str:
    """Hash a file without retaining sequencing data in memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def object_hash(value: Any) -> str:
    """Hash canonical finite JSON."""
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def now() -> str:
    """Return a UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: Any) -> None:
    """Atomically publish finite JSON in its destination directory."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def file_manifest(root: Path) -> dict[str, str]:
    """Hash every regular file under a resource/model root."""
    return {
        str(p.relative_to(root)): sha256_file(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"
    }


def installed_packages() -> dict[str, str]:
    """Capture distributions actually visible to the caller interpreter."""
    return dict(sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions()))


LADDER_TOOLS = ["minimap2", "samtools", "bcftools", "run_clair3.sh"]
HYBRID_TOOLS = ["samtools"]  # BAM input only; FASTQ input runs no external tool


def frozen_model_path(environment: dict[str, Any]) -> str | None:
    """The frozen Clair3 model path, or ``None`` for a hybrid-only (model-less) freeze."""
    model = environment.get("model")
    return None if model is None else str(model["path"])


def freeze_environment(checkout: Path, model: Path | None) -> dict[str, Any]:
    """Record caller sources, resources, locked packages, model and executable identity.

    Without ``model`` the environment can run only the hybrid engine: no Clair3 model is
    frozen, and only the tools the hybrid path may call are recorded.
    """
    tools = HYBRID_TOOLS if model is None else LADDER_TOOLS
    executables = {}
    for tool in tools:
        path = shutil.which(tool)
        if path is None:
            raise ValueError(f"missing executable: {tool}")
        executables[tool] = {"path": path, "sha256": sha256_file(Path(path))}
    if model is not None and (not model.is_dir() or not any(model.iterdir())):
        raise ValueError("model directory missing or empty")
    harness_root = Path(__file__).resolve().parents[2]
    return {
        "schema_version": 1,
        "created_at": now(),
        "caller_commit": run_tool(
            ["git", "-C", str(checkout), "rev-parse", "HEAD"], timeout=10
        ).strip(),
        "caller_sources": {
            str(Path(name).relative_to("src/muc_one_span")): sha256_file(checkout / name)
            for name in run_tool(
                ["git", "-C", str(checkout), "ls-files", "src/muc_one_span"], timeout=10
            ).splitlines()
        },
        "lock_sha256": sha256_file(checkout / "uv.lock"),
        "harness_commit": run_tool(
            ["git", "-C", str(harness_root), "rev-parse", "HEAD"], timeout=10
        ).strip(),
        "harness_sources": {
            p.name: sha256_file(p) for p in sorted(Path(__file__).parent.glob("clinical_*.py"))
        },
        "harness_entrypoint_sha256": sha256_file(harness_root / "scripts/clinical_benchmark.py"),
        "model": None
        if model is None
        else {"path": str(model.resolve()), "files": file_manifest(model)},
        "tools": executables,
        "tool_versions": get_tool_versions(tools),
        "python": sys.version,
        "python_executable": sys.executable,
        "packages": installed_packages(),
        "execution_environment": {"CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES")},
        "hardware": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpuinfo": Path("/proc/cpuinfo").read_text()
            if Path("/proc/cpuinfo").exists()
            else None,
            "meminfo": Path("/proc/meminfo").read_text()
            if Path("/proc/meminfo").exists()
            else None,
        },
    }


def verify_environment(environment: dict[str, Any], checkout: Path) -> None:
    """Reject changes to frozen caller, resources, lock, tools or model before running."""
    package = checkout / "src" / "muc_one_span"
    for name, digest in environment["caller_sources"].items():
        if sha256_file(package / name) != digest:
            raise ValueError(f"frozen caller source/resource mismatch: {name}")
    from muc_one_span import cli

    if Path(cli.__file__).resolve() != (package / "cli.py").resolve():
        raise ValueError("imported caller is not the verified checkout")
    if environment["execution_environment"] != {
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES")
    }:
        raise ValueError("frozen execution environment mismatch")
    if environment["python"] != sys.version:
        raise ValueError("frozen Python interpreter version mismatch")
    if environment["packages"] != installed_packages():
        raise ValueError("frozen installed packages mismatch")
    if sha256_file(checkout / "uv.lock") != environment["lock_sha256"]:
        raise ValueError("frozen dependency lock mismatch")
    for name, digest in environment["harness_sources"].items():
        if sha256_file(package / name) != digest:
            raise ValueError(f"frozen benchmark harness mismatch: {name}")
    if (
        sha256_file(checkout / "scripts/clinical_benchmark.py")
        != environment["harness_entrypoint_sha256"]
    ):
        raise ValueError("frozen benchmark entrypoint mismatch")
    model = environment["model"]
    if model is not None and file_manifest(Path(model["path"])) != model["files"]:
        raise ValueError("frozen model content mismatch")
    for name, metadata in environment["tools"].items():
        path = shutil.which(name)
        if path is None or sha256_file(Path(path)) != metadata["sha256"]:
            raise ValueError(f"frozen tool mismatch: {name}")
