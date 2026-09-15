"""Persist execution status without conflating scientific no-calls with crashes."""

from __future__ import annotations

import inspect
import json
import logging
from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import ParamSpec, TypeVar

P = ParamSpec("P")
R = TypeVar("R")


class InsufficientEvidenceError(ValueError):
    """The command executed but available evidence cannot support an allele call."""


def record_run_status(function: Callable[P, R]) -> Callable[P, R]:
    """Record all pipeline outcomes in its output_dir; retain exception/exit behavior."""

    @wraps(function)
    def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
        bound = inspect.signature(function).bind(*args, **kwargs)
        output = Path(bound.arguments["output_dir"])
        output.mkdir(parents=True, exist_ok=True)
        status_path = output / "run_status.json"

        def write(status: str, error: BaseException | None = None) -> None:
            data: dict[str, str | int] = {"schema_version": 1, "status": status}
            if error is not None:
                data.update(error_type=type(error).__name__, error=str(error))
            status_path.write_text(json.dumps(data, indent=2) + "\n")
            if status != "running":
                _update_summary_status(output / "summary.json", data)

        write("running")
        try:
            result = function(*args, **kwargs)
        except InsufficientEvidenceError as error:
            write("insufficient_evidence", error)
            raise
        except KeyboardInterrupt as error:
            write("interrupted", error)
            raise
        except BaseException as error:
            write("execution_failed", error)
            raise
        write("completed")
        return result

    return wrapped


def _update_summary_status(path: Path, status: dict[str, str | int]) -> None:
    """Keep portable summaries honest; malformed old artifacts cannot mask errors.

    The sidecar remains authoritative when a summary cannot be read or written.
    Updating an existing summary after a failed rerun marks its evidence as stale.
    """
    if not path.exists():
        return
    try:
        summary = json.loads(path.read_text())
        if not isinstance(summary, dict):
            raise ValueError("summary must be a JSON object")
        summary["run_status"] = status
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(summary, indent=2) + "\n")
        temporary.replace(path)
    except (OSError, ValueError) as error:
        logging.getLogger(__name__).warning("Could not update summary execution status: %s", error)
