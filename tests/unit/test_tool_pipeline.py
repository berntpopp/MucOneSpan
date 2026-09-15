"""Deterministic process lifecycle tests with external commands mocked."""

from __future__ import annotations

import io
import signal
from unittest.mock import MagicMock

import pytest

from muc_one_span.tool_pipeline import run_pipeline


def run_tool_pipeline(commands, *, timeout):
    return run_pipeline(commands, timeout=timeout, env={"PATH": "/usr/bin"})


@pytest.fixture
def runtime(mocker):
    clock = [0.0]
    children = []

    def sleep(seconds):
        clock[0] += seconds

    def launch(*args, **kwargs):
        child = MagicMock(pid=910001 + len(children), returncode=0)
        child.stdout = io.BytesIO() if kwargs["stdout"] == -1 else None
        child.stderr = io.BytesIO(f"diagnostic {len(children)}".encode())
        child.poll.side_effect = lambda: child.returncode
        child.wait.side_effect = lambda timeout: child.returncode
        children.append(child)
        return child

    class Drain:
        def __init__(self, stream):
            self.tail = stream.read()
            stream.close()
            self.error = None
            self.thread = MagicMock()
            self.thread.is_alive.return_value = False
            self.stop = MagicMock()

        def join(self, deadline):
            assert deadline <= 10

        def text(self):
            return self.tail.decode()

    def kill(pid, sig):
        for child in children:
            if child.pid == pid and child.returncode is None:
                child.returncode = -sig

    mocker.patch("muc_one_span.tool_pipeline.time.monotonic", side_effect=lambda: clock[0])
    sleeper = mocker.patch("muc_one_span.tool_pipeline.time.sleep", side_effect=sleep)
    mocker.patch("muc_one_span.tool_pipeline._StderrDrain", Drain)
    launcher = mocker.patch("muc_one_span.tool_pipeline.subprocess.Popen", side_effect=launch)
    killer = mocker.patch("muc_one_span.tool_pipeline.os.killpg", side_effect=kill)
    return children, launcher, killer, sleeper, clock, launch


def test_success_streams_isolates_and_passes_environment(runtime):
    children, launch, killer, _, _, _ = runtime
    run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert launch.call_args_list[1].kwargs["stdin"] is children[0].stdout
    assert children[0].stdout.closed
    assert launch.call_args_list[1].kwargs["stdout"] == -3  # DEVNULL
    for call in launch.call_args_list:
        assert call.kwargs["start_new_session"] is True
        assert call.kwargs["env"]["PATH"] == "/usr/bin"
    killer.assert_not_called()


@pytest.mark.parametrize("failing_index", [0, 1])
def test_either_failure_terminates_peer_with_both_diagnostics(runtime, failing_index):
    children, launcher, killer, _, _, launch = runtime

    def failed_launch(*args, **kwargs):
        child = launch(*args, **kwargs)
        child.returncode = 7 if len(children) - 1 == failing_index else None
        return child

    launcher.side_effect = failed_launch
    with pytest.raises(RuntimeError, match=r"diagnostic 0[\s\S]*diagnostic 1"):
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert {c.args[0] for c in killer.call_args_list} == {910001, 910002}
    assert all(c.returncode is not None for c in children)
    assert all(c.wait.call_args.kwargs["timeout"] < 10 for c in children)


@pytest.mark.parametrize("exception", [FileNotFoundError(), PermissionError(), KeyboardInterrupt()])
def test_second_launch_failure_cleans_up_producer(runtime, exception):
    children, launcher, killer, _, _, launch = runtime

    def fail_second(*args, **kwargs):
        if children:
            raise exception
        child = launch(*args, **kwargs)
        child.returncode = None
        return child

    launcher.side_effect = fail_second
    with pytest.raises(type(exception)):
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert children[0].stdout.closed
    assert children[0].wait.called
    assert killer.call_args_list[0].args == (910001, signal.SIGTERM)
    assert killer.call_args_list[-1].args == (910001, signal.SIGKILL)


@pytest.mark.parametrize("hanging_index", [0, 1])
def test_one_deadline_limits_waits_and_escalates(runtime, hanging_index):
    children, launcher, killer, _, clock, launch = runtime

    def hanging_launch(*args, **kwargs):
        child = launch(*args, **kwargs)
        child.returncode = None if len(children) - 1 == hanging_index else 0
        return child

    launcher.side_effect = hanging_launch
    with pytest.raises(TimeoutError, match="timed out"):
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert clock[0] <= 10
    assert all(c.wait.call_args.kwargs["timeout"] < 0.5 for c in children)
    assert any(call.args[1] == signal.SIGKILL for call in killer.call_args_list)


def test_interruption_during_monitoring_cleans_both_groups(runtime):
    children, launcher, killer, sleeper, _, launch = runtime

    def hanging_launch(*args, **kwargs):
        child = launch(*args, **kwargs)
        child.returncode = None
        return child

    launcher.side_effect = hanging_launch
    sleeper.side_effect = [KeyboardInterrupt(), None]
    with pytest.raises(KeyboardInterrupt):
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert {call.args[0] for call in killer.call_args_list} == {910001, 910002}
    assert all(c.wait.called for c in children)


def test_missing_first_executable_is_visible(runtime):
    _, launch, killer, _, _, _ = runtime
    launch.side_effect = FileNotFoundError()
    with pytest.raises(FileNotFoundError, match="Tool not found: producer"):
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    killer.assert_not_called()


def test_uncaptured_upstream_is_cleaned(runtime):
    children, launcher, _, _, _, launch = runtime

    def bad_launch(*args, **kwargs):
        child = launch(*args, **kwargs)
        child.stdout = None
        return child

    launcher.side_effect = bad_launch
    with pytest.raises(RuntimeError, match="stdout was not captured"):
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    assert children[0].wait.called


def test_stderr_reader_bounds_tail_and_closes_stream(mocker):
    from muc_one_span.tool_pipeline import _StderrDrain

    stream = MagicMock()
    stream.fileno.return_value = 19
    mocker.patch("muc_one_span.tool_pipeline.threading.Thread")
    mocker.patch("muc_one_span.tool_pipeline.os.set_blocking")
    mocker.patch("muc_one_span.tool_pipeline.select.select", return_value=([19], [], []))
    mocker.patch(
        "muc_one_span.tool_pipeline.os.read",
        side_effect=[b"x" * 1048576, b"useful diagnostic\xff", b""],
    )
    drain = _StderrDrain(stream)
    drain._read()
    assert len(drain.tail) == 1048576
    assert drain.text().endswith("useful diagnostic�")
    assert drain.error is None
    stream.close.assert_called_once()


def test_stderr_reader_surfaces_io_errors(mocker):
    from muc_one_span.tool_pipeline import _StderrDrain

    stream = MagicMock()
    stream.fileno.side_effect = OSError("broken fd")
    mocker.patch("muc_one_span.tool_pipeline.threading.Thread")
    drain = _StderrDrain(stream)
    drain._read()
    assert str(drain.error) == "broken fd"
    stream.close.assert_called_once()


def test_stderr_join_uses_remaining_deadline(mocker):
    from muc_one_span.tool_pipeline import _StderrDrain

    mocker.patch("muc_one_span.tool_pipeline.threading.Thread")
    mocker.patch("muc_one_span.tool_pipeline.time.monotonic", return_value=9)
    drain = _StderrDrain(MagicMock())
    drain.join(10)
    drain.thread.join.assert_called_once_with(timeout=1)


def test_cleanup_reports_child_that_cannot_be_reaped_before_deadline(mocker):
    import subprocess

    from muc_one_span.tool_pipeline import _cleanup

    child = MagicMock(pid=910001)
    child.wait.side_effect = subprocess.TimeoutExpired("producer", 0)
    child.poll.return_value = None
    mocker.patch("muc_one_span.tool_pipeline.time.monotonic", return_value=10)
    mocker.patch("muc_one_span.tool_pipeline.os.killpg")
    with pytest.raises(RuntimeError, match="cleanup did not finish"):
        _cleanup([child], [], 10)
    child.wait.assert_called_once_with(timeout=0)


def test_cleanup_drains_final_diagnostics_before_requesting_reader_stop(mocker):
    from muc_one_span.tool_pipeline import _cleanup

    child = MagicMock(pid=910001)
    child.poll.return_value = 0
    drain = MagicMock()
    drain.thread.is_alive.return_value = False
    stop_states = []
    drain.join.side_effect = lambda deadline: stop_states.append(drain.stop.set.called)
    mocker.patch("muc_one_span.tool_pipeline.time.monotonic", return_value=9)
    mocker.patch("muc_one_span.tool_pipeline.time.sleep")
    mocker.patch("muc_one_span.tool_pipeline.os.killpg")
    _cleanup([child], [drain], 10)
    assert stop_states[0] is False


def test_cleanup_failure_retains_primary_failure_and_both_tool_diagnostics(runtime, mocker):
    children, launcher, _, _, _, launch = runtime

    def failed_launch(*args, **kwargs):
        child = launch(*args, **kwargs)
        child.returncode = 7
        return child

    launcher.side_effect = failed_launch
    mocker.patch("muc_one_span.tool_pipeline._cleanup", side_effect=RuntimeError("cleanup expired"))
    with pytest.raises(RuntimeError) as failure:
        run_tool_pipeline([["producer"], ["consumer"]], timeout=10)
    message = str(failure.value)
    assert "Pipeline command failed" in message
    assert "diagnostic 0" in message
    assert "diagnostic 1" in message
    assert "Cleanup failed: cleanup expired" in message
    assert len(children) == 2
