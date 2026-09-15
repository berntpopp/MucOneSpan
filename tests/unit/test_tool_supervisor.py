"""Supervisor protocol and adoption tests without launching external tools."""

from __future__ import annotations

import json
import signal
import subprocess
from unittest.mock import MagicMock

import pytest

from muc_one_span import tool_supervisor as supervisor
from muc_one_span.tools import run_tool_pipeline


def test_tool_pipeline_cleans_path_before_supervision(mocker, monkeypatch):
    monkeypatch.setenv("PATH", "/project/.venv/bin:/usr/bin")
    run = mocker.patch("muc_one_span.tool_supervisor.run_supervised_pipeline")
    run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert run.call_args.kwargs["env"]["PATH"] == "/usr/bin"
    assert run.call_args.kwargs["timeout"] == 10


def test_nonlinux_retains_posix_cleanup_with_explicit_limit(mocker, caplog):
    mocker.patch.object(supervisor.sys, "platform", "darwin")
    run = mocker.patch.object(supervisor, "run_pipeline")
    supervisor.run_supervised_pipeline([["producer"]], timeout=10, env={})
    run.assert_called_once_with([["producer"]], timeout=10, env={}, stdout=None)
    assert "reaping is delegated to the OS" in caplog.text


@pytest.fixture
def parent(mocker):
    mocker.patch.object(supervisor.sys, "platform", "linux")
    mocker.patch.object(supervisor.time, "monotonic", return_value=100)
    proc = MagicMock(pid=12345, returncode=0)
    proc.poll.return_value = 0
    launch = mocker.patch.object(supervisor.subprocess, "Popen", return_value=proc)
    collect = mocker.patch.object(supervisor, "_collect", return_value={"result": "success"})
    kill = mocker.patch.object(supervisor.os, "kill")
    group_kill = mocker.patch.object(supervisor.os, "killpg")
    return proc, launch, collect, kill, group_kill


def test_supervisor_isolated_with_same_absolute_budget(parent):
    proc, launch, collect, _, _ = parent
    supervisor.run_supervised_pipeline([["producer"]], timeout=10, env={"PATH": "/usr/bin"})
    command = launch.call_args.args[0]
    assert command[:3] == [supervisor.sys.executable, "-m", "muc_one_span.tool_supervisor"]
    assert json.loads(command[3]) == [[["producer"]], 10, 109.9, None]
    assert launch.call_args.kwargs["start_new_session"] is True
    assert launch.call_args.kwargs["env"] == {"PATH": "/usr/bin"}
    assert collect.call_args.args[1] == 109.9
    assert proc.stdout.close.called and proc.stderr.close.called


def test_caller_interruption_waits_for_supervisor_cleanup(parent):
    _, _, collect, kill, _ = parent
    collect.side_effect = [KeyboardInterrupt(), {"result": "KeyboardInterrupt"}]
    with pytest.raises(KeyboardInterrupt):
        supervisor.run_supervised_pipeline([["producer"]], timeout=10, env={})
    kill.assert_called_once_with(12345, signal.SIGINT)
    assert collect.call_args_list[1].args[1] == 100.5


def test_supervisor_emergency_cleanup_reports_unreaped_child(parent):
    proc, _, collect, _, group_kill = parent

    def fail(process, deadline, groups, collector):
        groups.append(67890)
        raise TimeoutError("supervisor stalled")

    collect.side_effect = fail
    proc.wait.side_effect = subprocess.TimeoutExpired("supervisor", 10)
    proc.poll.return_value = None
    with pytest.raises(RuntimeError, match="Supervisor remains unreaped"):
        supervisor.run_supervised_pipeline([["producer"]], timeout=10, env={})
    assert (67890, signal.SIGKILL) in [call.args for call in group_kill.call_args_list]
    assert (12345, signal.SIGKILL) in [call.args for call in group_kill.call_args_list]
    assert proc.wait.call_args.kwargs["timeout"] == 10


@pytest.mark.parametrize("error_type", ["TimeoutError", "FileNotFoundError", "KeyboardInterrupt"])
def test_supervisor_returns_original_failure_type_and_diagnostics(parent, error_type):
    _, _, collect, _, _ = parent
    collect.return_value = {"result": error_type, "error": "both tool diagnostics"}
    expected = {
        "TimeoutError": TimeoutError,
        "FileNotFoundError": FileNotFoundError,
        "KeyboardInterrupt": KeyboardInterrupt,
    }[error_type]
    with pytest.raises(expected, match="both tool diagnostics"):
        supervisor.run_supervised_pipeline([["producer"]], timeout=10, env={})


def test_collector_handles_split_events_and_reaps_supervisor(mocker):
    proc = MagicMock()
    proc.stdout.fileno.return_value = 10
    proc.stderr.fileno.return_value = 11
    mocker.patch.object(supervisor.time, "monotonic", return_value=9)
    mocker.patch.object(
        supervisor.select,
        "select",
        side_effect=[
            ([proc.stdout], [], []),
            ([proc.stdout, proc.stderr], [], []),
            ([proc.stdout, proc.stderr], [], []),
        ],
    )
    mocker.patch.object(
        supervisor.os,
        "read",
        side_effect=[
            b'{"pid":',
            b' 101}\n{"result": "success"}\n',
            b"diagnostic",
            b"",
            b"",
        ],
    )
    groups = []
    assert supervisor._collect(proc, 10, groups) == {"result": "success"}
    assert groups == [101]
    proc.wait.assert_called_once_with(timeout=1)


def test_collector_reports_supervisor_crash_stderr(mocker):
    proc = MagicMock(returncode=1)
    mocker.patch.object(supervisor.time, "monotonic", return_value=9)
    mocker.patch.object(
        supervisor.select,
        "select",
        side_effect=[
            ([proc.stdout, proc.stderr], [], []),
            ([proc.stderr], [], []),
        ],
    )
    mocker.patch.object(supervisor.os, "read", side_effect=[b"", b"startup failure", b""])
    with pytest.raises(RuntimeError, match="startup failure"):
        supervisor._collect(proc, 10, [])


def test_adopted_descendants_are_killed_and_reaped(mocker):
    wait = mocker.patch.object(
        supervisor.os, "waitpid", side_effect=[(0, 0), (101, 0), ChildProcessError]
    )
    mocker.patch.object(supervisor.Path, "read_text", return_value="101")
    kill = mocker.patch.object(supervisor.os, "kill")
    mocker.patch.object(supervisor.time, "monotonic", return_value=9)
    mocker.patch.object(supervisor.time, "sleep")
    supervisor._reap_adopted(10)
    kill.assert_called_once_with(101, signal.SIGKILL)
    assert wait.call_count == 3


def test_adopted_reaping_has_no_wait_beyond_deadline(mocker):
    mocker.patch.object(supervisor.os, "waitpid", return_value=(0, 0))
    mocker.patch.object(supervisor.Path, "read_text", return_value="101")
    mocker.patch.object(supervisor.os, "kill")
    mocker.patch.object(supervisor.time, "monotonic", return_value=10)
    sleep = mocker.patch.object(supervisor.time, "sleep")
    with pytest.raises(RuntimeError, match="reaping exceeded"):
        supervisor._reap_adopted(10)
    sleep.assert_not_called()


def test_child_enables_subreaper_only_in_supervisor_and_reaps_on_failure(mocker):
    mocker.patch.object(
        supervisor.sys, "argv", ["supervisor", json.dumps([[["producer"]], 10, 20, None])]
    )
    mocker.patch.object(supervisor.signal, "signal")
    libc = MagicMock()
    libc.prctl.return_value = 0
    mocker.patch.object(supervisor.ctypes, "CDLL", return_value=libc)
    mocker.patch.object(
        supervisor, "run_pipeline", side_effect=RuntimeError("producer diagnostics")
    )
    reap = mocker.patch.object(supervisor, "_reap_adopted")
    event = mocker.patch.object(supervisor, "_event")
    supervisor._child_main()
    libc.prctl.assert_called_once_with(36, 1, 0, 0, 0)
    reap.assert_called_once_with(20)
    assert event.call_args.args[0] == {"result": "RuntimeError", "error": "producer diagnostics"}


def test_collector_retains_partial_event_across_interrupt(mocker):
    proc = MagicMock()
    mocker.patch.object(supervisor.time, "monotonic", return_value=9)
    mocker.patch.object(
        supervisor.select,
        "select",
        side_effect=[
            ([proc.stdout], [], []),
            KeyboardInterrupt(),
            ([proc.stdout], [], []),
            ([proc.stdout, proc.stderr], [], []),
        ],
    )
    mocker.patch.object(
        supervisor.os,
        "read",
        side_effect=[
            b'{"result":"RuntimeError","error":"diagnostic ',
            b'tail"}\n',
            b"",
            b"",
        ],
    )
    collector = supervisor._Collector(proc)
    with pytest.raises(KeyboardInterrupt):
        supervisor._collect(proc, 10, [], collector)
    assert supervisor._collect(proc, 10, [], collector) == {
        "result": "RuntimeError",
        "error": "diagnostic tail",
    }


def test_supervisor_stall_keeps_time_for_final_kill_and_reap(parent, mocker):
    proc, _, collect, _, _ = parent
    clock = [100.0]
    mocker.patch.object(supervisor.time, "monotonic", side_effect=lambda: clock[0])

    def expire(process, deadline, groups, collector):
        clock[0] = deadline
        raise TimeoutError("supervisor stalled")

    collect.side_effect = expire
    with pytest.raises(RuntimeError, match="Supervisor cleanup failed"):
        supervisor.run_supervised_pipeline([["producer"]], timeout=10, env={})
    assert 0 < proc.wait.call_args.kwargs["timeout"] <= 0.05 + 1e-9
    assert clock[0] < 110
