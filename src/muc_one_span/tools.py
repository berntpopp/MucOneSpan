# src/muc_one_span/tools.py
"""Subprocess helpers for running external bioinformatics tools."""

from __future__ import annotations

import contextlib
import logging
import os
import shutil
import signal
import subprocess
from collections.abc import Generator
from pathlib import Path

logger = logging.getLogger(__name__)


def _clean_path_for_externals(path_str: str) -> str:
    """Remove Python virtualenv bin directories from PATH.

    External tools like Clair3 internally call ``python3`` and expect the
    system or conda Python, not the project's virtualenv.  When running
    via ``uv run``, the ``.venv/bin`` directory is prepended to PATH,
    causing Clair3 to pick up the wrong Python (e.g. 3.11 instead of
    the conda env's 3.9 where tensorflow is installed).

    Uses two detection strategies (see CPython bugs.python.org/issue42041):

    1. **VIRTUAL_ENV env var:** If set, strip that exact ``bin`` directory.
    2. **Pattern matching:** Strip any PATH entry containing ``.venv/bin``
       or ``venv/bin`` as a fallback (catches cases where VIRTUAL_ENV is
       not set but a venv is on PATH).

    Conda environments are preserved — they use ``envs/`` paths, not
    ``venv`` or ``.venv``.
    """
    venv_dir = os.environ.get("VIRTUAL_ENV", "")
    venv_bin = str(Path(venv_dir) / "bin") if venv_dir else ""

    parts = path_str.split(os.pathsep)
    cleaned = [
        p
        for p in parts
        if not (
            (venv_bin and os.path.normpath(p) == os.path.normpath(venv_bin))
            or os.sep + ".venv" + os.sep in p + os.sep
            or os.sep + "venv" + os.sep in p + os.sep
        )
    ]
    return os.pathsep.join(cleaned)


def run_tool(
    cmd: list[str],
    cwd: str | None = None,
    timeout: float | None = None,
) -> str:
    """Run an external tool and return its stdout.

    Strips virtualenv ``bin`` directories from PATH so that external
    tools (Clair3, bcftools, samtools, etc.) use the system or conda
    Python rather than the project's virtualenv Python.

    Args:
        cmd: Command and arguments as a list.
        cwd: Optional working directory.
        timeout: Optional maximum seconds before terminating process group.

    Returns:
        Captured stdout as a string.

    Raises:
        FileNotFoundError: If the command is not found.
        TimeoutError: If execution exceeds timeout seconds.
        RuntimeError: If the command exits with non-zero status.
    """
    env = os.environ.copy()
    env["PATH"] = _clean_path_for_externals(env.get("PATH", ""))

    logger.debug("Running: %s", " ".join(str(c) for c in cmd))

    proc: subprocess.Popen[str] | None = None
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
            start_new_session=True,
        )
        stdout, stderr = proc.communicate(timeout=timeout)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Tool not found: {cmd[0]}") from exc
    except subprocess.TimeoutExpired as exc:
        if proc is not None:
            with contextlib.suppress(OSError):
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            proc.communicate()
        raise TimeoutError(f"Command {' '.join(cmd)} timed out after {timeout}s") from exc
    except BaseException:
        if proc is not None and proc.poll() is None:
            with contextlib.suppress(OSError):
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            proc.communicate()
        raise

    if proc.returncode != 0:
        raise RuntimeError(
            f"Command {' '.join(cmd)} failed with exit code {proc.returncode}.\nstderr: {stderr}"
        )

    return stdout


def run_tool_iter(
    cmd: list[str],
    cwd: str | None = None,
) -> Generator[str, None, None]:
    """Run an external tool and yield stdout lines without buffering.

    Unlike :func:`run_tool`, this does not capture all stdout in memory.
    Use for commands that may produce large output (e.g., ``samtools view``).

    Args:
        cmd: Command and arguments.
        cwd: Working directory for the command.

    Yields:
        Lines of stdout as strings.

    Raises:
        FileNotFoundError: If the tool is not found on PATH.
        RuntimeError: If the command exits with non-zero status.
    """
    env = os.environ.copy()
    env["PATH"] = _clean_path_for_externals(env.get("PATH", ""))

    logger.debug("Running (streaming): %s", " ".join(str(c) for c in cmd))

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        cwd=cwd,
        start_new_session=True,
    )

    # Drain stderr in a background thread to prevent blocking if the
    # subprocess writes enough to fill the stderr pipe buffer.
    import threading

    stderr_chunks: list[str] = []

    def _drain_stderr() -> None:
        if proc.stderr:
            stderr_chunks.append(proc.stderr.read())

    stderr_thread = threading.Thread(target=_drain_stderr)
    stderr_thread.start()

    try:
        if proc.stdout is not None:
            yield from proc.stdout
    finally:
        if proc.stdout is not None and hasattr(proc.stdout, "close"):
            proc.stdout.close()
        if hasattr(proc, "poll"):
            with contextlib.suppress(subprocess.TimeoutExpired, AttributeError):
                proc.wait(timeout=0.2)
            if proc.poll() is None:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                    proc.wait(timeout=1.0)
                except (OSError, subprocess.TimeoutExpired):
                    with contextlib.suppress(OSError):
                        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        stderr_thread.join(timeout=1.0)
        proc.wait()

    if proc.returncode != 0:
        stderr = stderr_chunks[0] if stderr_chunks else ""
        raise RuntimeError(
            f"Command failed: {' '.join(str(c) for c in cmd)}\n"
            f"Exit code: {proc.returncode}\nstderr: {stderr}"
        )


def check_tools(tools: list[str]) -> bool:
    """Verify that all required tools are available on PATH.

    Args:
        tools: List of tool names to check.

    Returns:
        True if all tools are found.

    Raises:
        RuntimeError: If any tools are missing, listing them all.
    """
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"Required tools not found: {', '.join(missing)}. "
            f"Install them or activate the conda environment."
        )
    logger.info("All required tools found: %s", ", ".join(tools))
    return True


def get_tool_versions(tools: list[str]) -> dict[str, str]:
    """Capture installed tool versions for reproducibility metadata.

    Attempts to run ``<tool> --version`` for each tool and captures the
    first line of output as the version string.

    Args:
        tools: List of tool names to query.

    Returns:
        Dictionary mapping tool name to version string, or "not found"
        if the tool is not available.
    """
    versions: dict[str, str] = {}
    for tool in tools:
        if not shutil.which(tool):
            versions[tool] = "not found"
            continue
        try:
            env = os.environ.copy()
            env["PATH"] = _clean_path_for_externals(env.get("PATH", ""))
            result = subprocess.run(
                [tool, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
                env=env,
            )
            first_line = result.stdout.strip().splitlines()[0] if result.stdout.strip() else ""
            if not first_line and result.stderr.strip():
                first_line = result.stderr.strip().splitlines()[0]
            versions[tool] = first_line or "unknown"
        except (subprocess.TimeoutExpired, FileNotFoundError, IndexError):
            versions[tool] = "unknown"
    return versions
