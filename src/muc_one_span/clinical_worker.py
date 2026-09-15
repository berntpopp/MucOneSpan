"""Isolated instrumented caller, launched only through the existing tool supervisor."""

from __future__ import annotations

import json
import resource
import sys
import time
from contextlib import ExitStack, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

from muc_one_span.benchmarking import _stage_patches, _timed
from muc_one_span.clinical_provenance import write_json


def execute(argv: list[str], output: Path) -> int:
    """Execute unchanged CLI with stage timing and isolated process resource counters."""
    from muc_one_span import report
    from muc_one_span.cli import main

    timings: dict[str, float] = {}
    started = time.monotonic()
    code = 0
    error = None
    with (output / "cli.log").open("w") as log, redirect_stdout(log), redirect_stderr(log):
        try:
            with ExitStack() as stack:
                for item in _stage_patches(timings):
                    stack.enter_context(item)
                stack.enter_context(
                    patch.object(
                        report, "generate_report", _timed("report", report.generate_report, timings)
                    )
                )
                main.main(args=argv, standalone_mode=False)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
            error = str(exc)
        except Exception as exc:
            import traceback

            traceback.print_exc()
            code, error = int(getattr(exc, "exit_code", 1)), str(exc)
        finally:
            timings["wall"] = time.monotonic() - started
            own = resource.getrusage(resource.RUSAGE_SELF)
            child = resource.getrusage(resource.RUSAGE_CHILDREN)
            write_json(
                output / "worker.json",
                {
                    "exit_code": code,
                    "error": error,
                    "timings": timings,
                    "peak_memory": {
                        "value": max(own.ru_maxrss, child.ru_maxrss),
                        "units": "KiB" if sys.platform == "linux" else "bytes",
                        "method": "max(RUSAGE_SELF.ru_maxrss, RUSAGE_CHILDREN.ru_maxrss)",
                        "scope": "largest individual process high-water RSS; lower bound on concurrent process-tree peak",
                    },
                },
            )
    return code


if __name__ == "__main__":
    config = json.loads(Path(sys.argv[1]).read_text())
    # Worker transport success is separate from the actual CLI exit in worker.json.
    # This lets the supervisor reserve failure for lifecycle/transport errors.
    execute(config["argv"], Path(config["output"]))
