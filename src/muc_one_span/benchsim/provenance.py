"""Provenance of benchmark results: harness, caller and simulator versions.

``run`` writes ``<results_root>/<engine>/caller.json`` (`caller_record`): the caller
package version and the Git commit of the checkout it ran from. The caller runs
in-process (`muc_one_span.benchmarking.run_pipeline`), so its version is this
package's. ``report`` copies those records into ``report.json["provenance"]``
(`report_provenance`), together with the harness version and commit at report time
and a count of the MucOneUp versions recorded in the split manifest.
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

from muc_one_span.version import __version__

CALLER_FILE = "caller.json"


def git_commit(cwd: Path, runner: Callable[..., str]) -> str | None:
    """``git rev-parse HEAD`` in ``cwd`` through ``runner`` (``None`` outside Git)."""
    try:
        out = runner(["git", "rev-parse", "HEAD"], cwd=str(cwd)).strip()
    except (FileNotFoundError, RuntimeError):
        return None
    return out or None


def caller_record(engine: str, commit: str | None) -> dict[str, Any]:
    """What ran ``engine``: the in-process caller's version and checkout commit."""
    return {"engine": engine, "caller_version": __version__, "caller_commit": commit}


def write_caller(engine_dir: Path, record: Mapping[str, Any]) -> Path:
    """Write ``caller.json`` into ``engine_dir``; return its path."""
    engine_dir.mkdir(parents=True, exist_ok=True)
    path = engine_dir / CALLER_FILE
    path.write_text(json.dumps(dict(record), indent=2, sort_keys=True) + "\n")
    return path


def read_caller(engine_dir: Path) -> dict[str, Any] | None:
    """The ``caller.json`` of ``engine_dir`` (``None`` for results run before it existed)."""
    path = engine_dir / CALLER_FILE
    if not path.is_file():
        return None
    data: dict[str, Any] = json.loads(path.read_text())
    return data


def report_provenance(
    harness_commit: str | None,
    engine_dirs: Mapping[str, Path],
    cases: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """``report.json["provenance"]``: harness, per-engine caller and MucOneUp versions."""
    versions = Counter(str(case.get("muconeup_version")) for case in cases.values())
    return {
        "harness_version": __version__,
        "harness_commit": harness_commit,
        "callers": {engine: read_caller(path) for engine, path in engine_dirs.items()},
        "muconeup_versions": dict(sorted(versions.items())),
    }
