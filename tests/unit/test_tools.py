# tests/unit/test_tools.py
"""Tests for subprocess helper utilities."""

from __future__ import annotations

import logging

import pytest

from muc_one_span.tools import check_tools, run_tool


class TestRunTool:
    """Tests for the run_tool function."""

    def test_run_tool_success(self):
        """run_tool returns stdout for a successful command."""
        result = run_tool(["echo", "hello"])
        assert result.strip() == "hello"

    def test_run_tool_failure_raises(self):
        """run_tool raises RuntimeError on non-zero exit."""
        with pytest.raises(RuntimeError, match="failed with exit code"):
            run_tool(["false"])

    def test_run_tool_not_found_raises(self):
        """run_tool raises FileNotFoundError for missing commands."""
        with pytest.raises(FileNotFoundError):
            run_tool(["nonexistent_tool_xyz"])

    def test_run_tool_captures_stderr(self):
        """run_tool includes stderr in error message."""
        with pytest.raises(RuntimeError, match="No such file"):
            run_tool(["ls", "/nonexistent_path_xyz"])


class TestRunToolPathSanitization:
    """Tests for PATH sanitization in run_tool."""

    def test_venv_bin_stripped_from_path(self):
        """run_tool strips virtualenv bin dirs from PATH for subprocesses."""
        from muc_one_span.tools import _clean_path_for_externals

        # Simulate a PATH with .venv/bin entries
        fake_path = "/home/user/project/.venv/bin:/usr/local/bin:/usr/bin"
        cleaned = _clean_path_for_externals(fake_path)
        assert ".venv/bin" not in cleaned
        assert "/usr/local/bin" in cleaned
        assert "/usr/bin" in cleaned

    def test_non_venv_paths_preserved(self):
        """run_tool keeps conda and system paths intact."""
        from muc_one_span.tools import _clean_path_for_externals

        fake_path = "/home/user/miniforge3/envs/env_clair3/bin:/usr/bin:/home/user/.venv/bin"
        cleaned = _clean_path_for_externals(fake_path)
        assert "env_clair3/bin" in cleaned
        assert "/usr/bin" in cleaned
        assert ".venv/bin" not in cleaned

    def test_run_tool_uses_cleaned_path(self):
        """run_tool subprocess does not see .venv/bin in PATH."""
        import os

        # Get the PATH that run_tool would pass to subprocess
        result = run_tool(["printenv", "PATH"])
        # If we're in a venv, the .venv/bin should NOT appear
        venv_dir = os.environ.get("VIRTUAL_ENV", "")
        if venv_dir:
            assert f"{venv_dir}/bin" not in result


class TestCheckTools:
    """Tests for the check_tools function."""

    def test_check_tools_all_present(self):
        """check_tools returns True when all tools are found."""
        assert check_tools(["echo", "ls"]) is True

    def test_check_tools_missing_tool(self):
        """check_tools raises RuntimeError listing missing tools."""
        with pytest.raises(RuntimeError, match="nonexistent_tool_xyz"):
            check_tools(["echo", "nonexistent_tool_xyz"])

    def test_check_tools_empty_list(self):
        """check_tools with empty list returns True."""
        assert check_tools([]) is True


def test_get_tool_versions_returns_versions(mocker):
    """get_tool_versions captures version strings from external tools."""
    mocker.patch(
        "muc_one_span.tools.subprocess.run",
        return_value=mocker.MagicMock(returncode=0, stdout="minimap2 2.28-r1209\n"),
    )
    mocker.patch("muc_one_span.tools.shutil.which", return_value="/usr/bin/minimap2")

    from muc_one_span.tools import get_tool_versions

    versions = get_tool_versions(["minimap2"])
    assert "minimap2" in versions
    assert "2.28" in versions["minimap2"]


def test_get_tool_versions_handles_missing_tool(mocker):
    """get_tool_versions returns 'not found' for missing tools."""
    mocker.patch("muc_one_span.tools.shutil.which", return_value=None)

    from muc_one_span.tools import get_tool_versions

    versions = get_tool_versions(["nonexistent_tool"])
    assert versions["nonexistent_tool"] == "not found"


def test_run_tool_logs_command(mocker, caplog):
    """run_tool logs the command at DEBUG level."""
    mocker.patch(
        "muc_one_span.tools.subprocess.run",
        return_value=mocker.MagicMock(returncode=0, stdout="output"),
    )
    mocker.patch("muc_one_span.tools.shutil.which", return_value="/usr/bin/echo")

    with caplog.at_level(logging.DEBUG, logger="muc_one_span.tools"):
        from muc_one_span.tools import run_tool

        run_tool(["echo", "hello"])

    assert any("echo hello" in record.message for record in caplog.records)


def test_run_tool_iter_yields_lines(mocker):
    """run_tool_iter yields stdout lines without buffering."""
    from unittest.mock import MagicMock

    from muc_one_span.tools import run_tool_iter

    mock_proc = MagicMock()
    mock_proc.stdout = iter(["line1\n", "line2\n", "line3\n"])
    mock_proc.wait.return_value = 0
    mock_proc.returncode = 0

    mocker.patch("muc_one_span.tools.subprocess.Popen", return_value=mock_proc)

    lines = list(run_tool_iter(["echo", "hello"]))
    assert lines == ["line1\n", "line2\n", "line3\n"]


def test_run_tool_iter_raises_on_failure(mocker):
    """run_tool_iter raises RuntimeError on non-zero exit."""
    from unittest.mock import MagicMock

    from muc_one_span.tools import run_tool_iter

    mock_proc = MagicMock()
    mock_proc.stdout = iter([])
    mock_proc.wait.return_value = 1
    mock_proc.returncode = 1
    mock_proc.stderr = MagicMock()
    mock_proc.stderr.read.return_value = "error output"

    mocker.patch("muc_one_span.tools.subprocess.Popen", return_value=mock_proc)

    with pytest.raises(RuntimeError, match="failed"):
        list(run_tool_iter(["failing_tool"]))


def test_run_tool_timeout_raises_and_cleans_up(mocker):
    """run_tool raises TimeoutError on timeout and kills process group."""
    import subprocess
    from unittest.mock import MagicMock

    from muc_one_span.tools import run_tool

    mock_proc = MagicMock()
    mock_proc.pid = 12345
    mock_proc.communicate.side_effect = [
        subprocess.TimeoutExpired(cmd=["sleep", "10"], timeout=0.1),
        ("out", "err"),
    ]

    mocker.patch("muc_one_span.tools.subprocess.Popen", return_value=mock_proc)
    mocker.patch("os.getpgid", return_value=12345)
    killpg_mock = mocker.patch("os.killpg")

    with pytest.raises(TimeoutError, match="timed out"):
        run_tool(["sleep", "10"], timeout=0.1)

    assert killpg_mock.called
