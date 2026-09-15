"""Mapping lifecycle regressions at the command boundary."""

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.mapping import map_reads


@pytest.mark.parametrize("failure", [RuntimeError("failed"), KeyboardInterrupt()])
def test_failed_mapping_removes_partial_bam_and_indexes(tmp_path, failure):
    artifacts = [tmp_path / name for name in ("mapping.bam", "mapping.bam.bai", "mapping.bam.csi")]
    for artifact in artifacts:
        artifact.write_text("stale or partial output")
    with (
        patch("muc_one_span.mapping._run_mapping_pipeline", side_effect=failure),
        pytest.raises(type(failure)),
    ):
        map_reads(Path("reads.fq"), Path("reference.fa"), tmp_path)
    assert not any(artifact.exists() for artifact in artifacts)


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf"), True, False])
def test_mapping_rejects_invalid_timeout_before_launch(tmp_path, timeout):
    with (
        patch("muc_one_span.mapping._run_mapping_pipeline") as pipeline,
        patch("muc_one_span.mapping.run_tool_pipeline"),
    ):
        with pytest.raises(ValueError, match=r"finite.*positive"):
            map_reads(Path("reads.fq"), Path("reference.fa"), tmp_path, timeout=timeout)
        pipeline.assert_not_called()


@pytest.mark.parametrize("suffix", ["fq", "bam"])
def test_mapping_stages_share_remaining_budget(tmp_path, mocker, suffix):
    ticks = [0, 1, 3] if suffix == "fq" else [0, 1, 3, 5]
    mocker.patch("muc_one_span.mapping.time.monotonic", side_effect=ticks)
    pipeline = mocker.patch("muc_one_span.mapping._run_mapping_pipeline")
    tools = mocker.patch("muc_one_span.mapping.run_tool_pipeline")
    map_reads(Path(f"reads.{suffix}"), Path("reference.fa"), tmp_path, timeout=10)
    assert pipeline.call_args.kwargs["timeout"] == (9 if suffix == "fq" else 7)
    assert tools.call_args.kwargs["timeout"] == (7 if suffix == "fq" else 5)
    if suffix == "bam":
        assert tools.call_args_list[0].kwargs["timeout"] == 9
        assert tools.call_args_list[0].kwargs["stdout"].closed


def test_index_failure_removes_new_bam_and_partial_index(tmp_path, mocker):
    bam = tmp_path / "mapping.bam"
    index = tmp_path / "mapping.bam.bai"
    mocker.patch(
        "muc_one_span.mapping._run_mapping_pipeline", side_effect=lambda *a, **k: bam.touch()
    )

    def fail_index(*args, **kwargs):
        index.touch()
        raise RuntimeError("index failed")

    mocker.patch("muc_one_span.mapping.run_tool_pipeline", side_effect=fail_index)
    with pytest.raises(RuntimeError, match="index failed"):
        map_reads(Path("reads.fq"), Path("reference.fa"), tmp_path)
    assert not bam.exists()
    assert not index.exists()


def test_exhausted_budget_does_not_launch_index(tmp_path, mocker):
    mocker.patch("muc_one_span.mapping.time.monotonic", side_effect=[0, 1, 11])
    mocker.patch("muc_one_span.mapping._run_mapping_pipeline")
    tool = mocker.patch("muc_one_span.mapping.run_tool_pipeline")
    with pytest.raises(TimeoutError, match="before the next stage"):
        map_reads(Path("reads.fq"), Path("reference.fa"), tmp_path, timeout=10)
    tool.assert_not_called()


def test_mapping_rejects_input_output_collision_without_deleting_source(tmp_path, mocker):
    source = tmp_path / "mapping.bam"
    source.write_bytes(b"original reads")
    mocker.patch("muc_one_span.mapping.run_tool_pipeline")
    mocker.patch("muc_one_span.mapping._run_mapping_pipeline")
    with pytest.raises(ValueError, match=r"input.*output"):
        map_reads(source, Path("reference.fa"), tmp_path)
    assert source.read_bytes() == b"original reads"
