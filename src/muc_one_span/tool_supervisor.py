"""Isolated Linux subreaper and caller control channel for tool pipelines."""

from __future__ import annotations

import contextlib
import ctypes
import json
import logging
import os
import select
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import BinaryIO

from muc_one_span.tool_pipeline import run_pipeline, validate_timeout

logger = logging.getLogger(__name__)


def _event(value: dict[str, object]) -> None:
    print(json.dumps(value), flush=True)


def _reap_adopted(deadline: float) -> None:
    """Kill and reap adopted descendants without changing the caller's state."""
    children_path = Path(f"/proc/self/task/{os.getpid()}/children")
    while True:
        try:
            pid, _ = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if pid:
            continue
        for child in children_path.read_text().split():
            with contextlib.suppress(ProcessLookupError):
                os.kill(int(child), signal.SIGKILL)
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RuntimeError("Descendant reaping exceeded the pipeline deadline")
        time.sleep(min(0.005, remaining))


def _child_main() -> None:
    commands, timeout, deadline, output_fd = json.loads(sys.argv[1])
    # This flag affects only the isolated supervisor, never the application or
    # unrelated subprocesses. Orphans in its process tree are adopted here.
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Could not enable child subreaper")
    error: BaseException | None = None
    output = os.fdopen(output_fd, "wb", closefd=False) if output_fd is not None else None
    try:
        run_pipeline(
            commands,
            timeout=timeout,
            env=os.environ.copy(),
            stdout=output,
            deadline=deadline,
            on_launch=lambda pid: _event({"pid": pid}),
        )
    except BaseException as exc:
        error = exc
    finally:
        # A retry from the caller must not interrupt adopted-child reaping.
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        if output is not None:
            output.close()
        try:
            _reap_adopted(deadline)
        except BaseException as cleanup_error:
            error = RuntimeError(f"{error or ''}\nCleanup failed: {cleanup_error}")
    _event(
        {"result": "success" if error is None else type(error).__name__, "error": str(error or "")}
    )


def _raise_result(result: dict[str, object]) -> None:
    status = result["result"]
    if status == "success":
        return
    errors: dict[str, type[BaseException]] = {
        "TimeoutError": TimeoutError,
        "FileNotFoundError": FileNotFoundError,
        "PermissionError": PermissionError,
        "KeyboardInterrupt": KeyboardInterrupt,
    }
    raise errors.get(str(status), RuntimeError)(str(result.get("error", "")))


class _Collector:
    """Keep protocol and pipe state intact when caller interruption is retried."""

    def __init__(self, proc: subprocess.Popen[bytes]) -> None:
        assert proc.stdout is not None
        assert proc.stderr is not None
        self.streams = [proc.stdout, proc.stderr]
        self.pending = b""
        self.stderr = bytearray()
        self.result: dict[str, object] | None = None


def _collect(
    proc: subprocess.Popen[bytes],
    deadline: float,
    groups: list[int],
    state: _Collector | None = None,
) -> dict[str, object]:
    state = _Collector(proc) if state is None else state
    while state.streams:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise TimeoutError("Pipeline supervisor exceeded the lifecycle deadline")
        ready, _, _ = select.select(state.streams, [], [], min(0.01, remaining))
        for stream in ready:
            chunk = os.read(stream.fileno(), 65536)
            if not chunk:
                state.streams.remove(stream)
            elif stream is proc.stderr:
                state.stderr.extend(chunk)
                del state.stderr[:-1048576]
            else:
                state.pending += chunk
                while b"\n" in state.pending:
                    line, state.pending = state.pending.split(b"\n", 1)
                    event = json.loads(line)
                    if "pid" in event:
                        groups.append(event["pid"])
                    if "result" in event:
                        state.result = event
    proc.wait(timeout=max(0.0, deadline - time.monotonic()))
    if state.result is None:
        raise RuntimeError(
            f"Pipeline supervisor exited with {proc.returncode}: "
            f"{state.stderr.decode(errors='replace')}"
        )
    return state.result


def run_supervised_pipeline(
    commands: list[list[str]],
    *,
    timeout: float,
    env: dict[str, str],
    stdout: BinaryIO | None = None,
) -> None:
    """Run in an isolated Linux subreaper; retain POSIX fallback elsewhere.

    The parent and supervisor share one monotonic deadline. Linux supervises
    and reaps every descendant, including children orphaned before shutdown.
    Other POSIX systems retain process-group termination; their init process
    owns orphan reaping because there is no portable subreaper interface.
    """
    validate_timeout(timeout)
    if sys.platform != "linux":
        logger.warning("Orphan descendant reaping is delegated to the OS on non-Linux platforms")
        run_pipeline(commands, timeout=timeout, env=env, stdout=stdout)
        return
    deadline = time.monotonic() + timeout
    output_fd = stdout.fileno() if stdout is not None else None
    # Keep the last fraction of the caller's deadline for reaping the supervisor.
    reap_reserve = min(0.05, timeout * 0.01)
    supervisor_deadline = deadline - 2 * reap_reserve
    configuration = json.dumps([commands, timeout, supervisor_deadline, output_fd])
    proc = subprocess.Popen(
        [sys.executable, "-m", "muc_one_span.tool_supervisor", configuration],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        start_new_session=True,
        pass_fds=() if output_fd is None else (output_fd,),
    )
    groups: list[int] = []
    collector = _Collector(proc)
    try:
        result = _collect(proc, supervisor_deadline, groups, collector)
    except BaseException as original:
        # Give the supervisor the remaining budget to clean and reap its tree.
        with contextlib.suppress(ProcessLookupError):
            os.kill(proc.pid, signal.SIGINT)
        cleanup_deadline = min(deadline - reap_reserve, time.monotonic() + min(0.5, timeout * 0.1))
        try:
            while True:
                try:
                    result = _collect(proc, cleanup_deadline, groups, collector)
                    break
                except KeyboardInterrupt:
                    if time.monotonic() >= cleanup_deadline:
                        raise
        except BaseException as cleanup_error:
            for group in groups:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(group, signal.SIGKILL)
            with contextlib.suppress(ProcessLookupError):
                os.killpg(proc.pid, signal.SIGKILL)
            with contextlib.suppress(subprocess.TimeoutExpired):
                proc.wait(timeout=max(0.0, deadline - time.monotonic()))
            incomplete = " Supervisor remains unreaped." if proc.poll() is None else ""
            raise RuntimeError(
                f"{original}\nSupervisor cleanup failed: {cleanup_error}{incomplete}"
            ) from original
        raise original
    finally:
        if proc.stdout is not None:
            proc.stdout.close()
        if proc.stderr is not None:
            proc.stderr.close()
    _raise_result(result)


if __name__ == "__main__":
    _child_main()
