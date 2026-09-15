"""Bounded real-process lifecycle tests; no bioinformatics tools required."""

from __future__ import annotations

import contextlib
import json
import os
import signal
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

CHILD = """
import os, pathlib, signal, subprocess, sys, time
role, mode, root = sys.argv[1:]
root = pathlib.Path(root)
(root / (role + ".pid")).write_text(str(os.getpid()))
(root / (role + ".ready")).touch()
if role == "descendant":
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    while True: time.sleep(1)
if mode.startswith("descendant") and role == "producer":
    options = {} if mode == "descendant" else dict(
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.Popen([sys.executable, __file__, "descendant", mode, str(root)], **options)
    while not (root / "descendant.ready").exists(): time.sleep(0.01)
    sys.exit(7 if mode.endswith("failure") else 0)
if mode == role + "_failure":
    other = "consumer" if role == "producer" else "producer"
    while not (root / (other + ".ready")).exists(): time.sleep(0.01)
    sys.stderr.write(role + " diagnostic marker\\n")
    sys.stderr.flush()
    sys.exit(7)
if mode in (role + "_hang", "interrupt", "missing") or (
    mode.endswith("_failure") and mode != role + "_failure"
):
    sys.stderr.write(role + " waiting marker\\n")
    sys.stderr.flush()
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    while True: time.sleep(1)
if mode == "heavy":
    sys.stderr.write("x" * (2 * 1024 * 1024) + role + " diagnostic tail\\n")
    sys.stderr.flush()
if role == "producer":
    sys.stdout.buffer.write(b"SAM stream\\n" * 10000)
else:
    while chunk := sys.stdin.buffer.read(8192):
        sys.stdout.buffer.write(chunk)
"""

HARNESS = """
import json, os, pathlib, signal, subprocess, sys, threading, time
from muc_one_span.tools import run_tool_pipeline
root = pathlib.Path(sys.argv[1])
mode = sys.argv[2]
original = subprocess.Popen
children = []
def launch(*args, **kwargs):
    if children:
        end = time.monotonic() + 2
        while not (root / "producer.ready").exists():
            assert time.monotonic() < end, "producer handshake expired"
            time.sleep(0.01)
    proc = original(*args, **kwargs)
    children.append(proc)
    (root / "launched.json").write_text(json.dumps([p.pid for p in children]))
    return proc
subprocess.Popen = launch
def interrupt():
    while not (root / "consumer.ready").exists(): time.sleep(0.01)
    os.kill(os.getpid(), signal.SIGINT)
if mode == "interrupt": threading.Thread(target=interrupt, daemon=True).start()
commands = [[sys.executable, str(root / "child.py"), role, mode, str(root)]
            for role in ("producer", "consumer")]
if mode == "missing": commands[1] = [str(root / "does-not-exist")]
started = time.monotonic()
try:
    with (root / "output").open("wb") as output:
        run_tool_pipeline(commands, timeout=2, stdout=output)
    result = {"status": "success"}
except BaseException as exc:
    result = {"status": type(exc).__name__, "error": str(exc)}
result["elapsed"] = time.monotonic() - started
result["returncodes"] = [p.returncode for p in children]
(root / "result.json").write_text(json.dumps(result))
"""


def _is_running(pid: int) -> bool:
    """A zombie still exists and fails the descendant-reaping contract."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


@pytest.mark.parametrize(
    ("mode", "status"),
    [
        ("success", "success"),
        ("heavy", "success"),
        ("producer_failure", "RuntimeError"),
        ("consumer_failure", "RuntimeError"),
        ("missing", "FileNotFoundError"),
        ("producer_hang", "TimeoutError"),
        ("consumer_hang", "TimeoutError"),
        ("interrupt", "KeyboardInterrupt"),
        ("descendant", "TimeoutError"),
        ("descendant_closed_failure", "RuntimeError"),
    ],
)
def test_pipeline_lifecycle_real_processes(tmp_path: Path, mode: str, status: str) -> None:
    (tmp_path / "child.py").write_text(CHILD)
    (tmp_path / "harness.py").write_text(HARNESS)
    try:
        subprocess.run(
            [sys.executable, str(tmp_path / "harness.py"), str(tmp_path), mode],
            check=True,
            timeout=8,
            capture_output=True,
            text=True,
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["status"] == status, result
        assert result["elapsed"] < 2.5, result
        assert all(code is not None for code in result["returncodes"]), result
        for pid_file in tmp_path.glob("*.pid"):
            assert not _is_running(int(pid_file.read_text())), pid_file.name
        if mode in ("success", "heavy"):
            assert (tmp_path / "output").read_bytes() == b"SAM stream\n" * 10000
        if mode in ("producer_failure", "consumer_failure"):
            failing = mode.split("_")[0]
            waiting = "consumer" if failing == "producer" else "producer"
            assert failing + " diagnostic marker" in result["error"]
            assert waiting + " waiting marker" in result["error"]
    finally:
        # Protect the test runner even if a regression defeats the inner timeout.
        launched = tmp_path / "launched.json"
        if launched.exists():
            for pid in json.loads(launched.read_text()):
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(pid, signal.SIGKILL)
