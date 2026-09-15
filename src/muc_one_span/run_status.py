"""Persist execution status without conflating scientific no-calls with crashes."""

from __future__ import annotations

import inspect
import json
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
            data = {"schema_version": 1, "status": status}
            if error is not None:
                data.update(error_type=type(error).__name__, error=str(error))
            status_path.write_text(json.dumps(data, indent=2) + "\n")

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
