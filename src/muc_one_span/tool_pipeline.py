"""Bounded POSIX streaming process lifecycle, used through :mod:`tools`."""

from __future__ import annotations

import contextlib
import math
import os
import select
import signal
import subprocess
import threading
import time
from collections.abc import Callable
from typing import IO


class _StderrDrain:
    """Drain concurrently, retaining the last MiB without blocking shutdown."""

    def __init__(self, stream: IO[bytes]) -> None:
        self.stream = stream
        self.tail = bytearray()
        self.error: OSError | None = None
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._read, daemon=True)
        self.thread.start()

    def _read(self) -> None:
        try:
            fd = self.stream.fileno()
            os.set_blocking(fd, False)
            while not self.stop.is_set():
                if not select.select([fd], [], [], 0.01)[0]:
                    continue
                chunk = os.read(fd, 65536)
                if not chunk:
                    break
                self.tail.extend(chunk)
                del self.tail[:-1048576]
        except OSError as exc:
            self.error = exc
        finally:
            self.stream.close()

    def join(self, deadline: float) -> None:
        self.thread.join(timeout=max(0.0, deadline - time.monotonic()))

    def text(self) -> str:
        return self.tail.decode(errors="replace")


def validate_timeout(timeout: float) -> None:
    """Reject missing, infinite, or nonpositive lifecycle budgets."""
    if isinstance(timeout, bool) or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be finite and positive")


def _signal_groups(processes: list[subprocess.Popen[bytes]], sig: signal.Signals) -> None:
    for proc in processes:
        # start_new_session makes PID the group ID. getpgid fails once a parent
        # exits, even when its descendants still hold the pipeline pipes open.
        with contextlib.suppress(ProcessLookupError):
            os.killpg(proc.pid, sig)


def _cleanup(
    processes: list[subprocess.Popen[bytes]], drains: list[_StderrDrain], deadline: float
) -> None:
    _signal_groups(processes, signal.SIGTERM)
    grace = min(max(0.0, deadline - time.monotonic()) / 3, 0.1)
    if grace:
        time.sleep(grace)
    # Signal even groups whose direct parent has exited, to kill descendants.
    _signal_groups(processes, signal.SIGKILL)
    for proc in processes:
        # An uninterruptible kernel operation can defeat even SIGKILL. Report
        # incomplete cleanup below instead of adding an unbounded wait.
        with contextlib.suppress(subprocess.TimeoutExpired):
            proc.wait(timeout=max(0.0, deadline - time.monotonic()))
    # Let readers consume final diagnostics through EOF before cancellation.
    # Keep a short portion of the same deadline for their nonblocking stop loop.
    drain_deadline = deadline - min(0.02, max(0.0, deadline - time.monotonic()) / 2)
    for drain in drains:
        drain.join(drain_deadline)
    for drain in drains:
        drain.stop.set()
    for drain in drains:
        drain.join(deadline)
    if any(proc.poll() is None for proc in processes) or any(
        drain.thread.is_alive() for drain in drains
    ):
        raise RuntimeError("Process cleanup did not finish within the pipeline timeout")


def run_pipeline(
    commands: list[list[str]],
    *,
    timeout: float,
    env: dict[str, str],
    stdout: IO[bytes] | None = None,
    deadline: float | None = None,
    on_launch: Callable[[int], None] | None = None,
) -> None:
    """Stream commands with one deadline for execution, cleanup, and drain joins.

    Reserve up to 0.5 seconds (10% for short budgets) inside the total timeout
    for termination. Each child starts an isolated session; descendants must
    remain in that process group. Stderr diagnostics retain each tool's last MiB.
    The final stdout is discarded unless a file handle is supplied.
    """
    validate_timeout(timeout)
    deadline = time.monotonic() + timeout if deadline is None else deadline
    work_deadline = deadline - min(0.5, timeout * 0.1)
    processes: list[subprocess.Popen[bytes]] = []
    drains: list[_StderrDrain] = []
    upstream: IO[bytes] | None = None
    try:
        for index, command in enumerate(commands):
            if time.monotonic() >= work_deadline:
                raise TimeoutError("Pipeline timed out before tool launch")
            try:
                proc = subprocess.Popen(
                    command,
                    stdin=upstream,
                    stdout=subprocess.PIPE
                    if index < len(commands) - 1
                    else stdout or subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    env=env,
                    start_new_session=True,
                )
            except FileNotFoundError as exc:
                raise FileNotFoundError(f"Tool not found: {command[0]}") from exc
            processes.append(proc)
            if on_launch is not None:
                on_launch(proc.pid)
            if upstream is not None:
                upstream.close()
            upstream = proc.stdout
            if index < len(commands) - 1 and upstream is None:
                raise RuntimeError(f"{command[0]} process stdout was not captured")
            assert proc.stderr is not None
            drains.append(_StderrDrain(proc.stderr))
        while True:
            codes = [proc.poll() for proc in processes]
            if any(code is not None and code != 0 for code in codes):
                raise RuntimeError("Pipeline command failed")
            if all(code == 0 for code in codes) and all(
                not drain.thread.is_alive() for drain in drains
            ):
                break
            remaining = work_deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f"Pipeline timed out after {timeout}s (including cleanup)")
            time.sleep(min(0.01, remaining))
        if any(drain.error is not None for drain in drains):
            raise RuntimeError("Could not read pipeline stderr")
    except BaseException as exc:
        if upstream is not None:
            upstream.close()
        cleanup_error: BaseException | None = None
        try:
            _cleanup(processes, drains, deadline)
        except BaseException as failure:
            cleanup_error = failure
        diagnostics = "\n".join(
            f"{' '.join(command[:2])} failed with exit code {proc.returncode}.\n"
            f"stderr: {drain.text()}"
            for command, proc, drain in zip(commands, processes, drains, strict=False)
        )
        if cleanup_error is not None:
            raise RuntimeError(f"{exc}\n{diagnostics}\nCleanup failed: {cleanup_error}") from exc
        if isinstance(exc, (RuntimeError, TimeoutError, FileNotFoundError)):
            raise type(exc)(f"{exc}\n{diagnostics}") from exc
        raise
    finally:
        if upstream is not None:
            upstream.close()
