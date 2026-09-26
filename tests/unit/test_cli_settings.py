"""Configuration precedence and provenance through real Click parsing."""

import json
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main


def config_file(tmp_path, **sections):
    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"schema_version": 1, **sections}))
    return path


@pytest.mark.parametrize("extra,expected", [([], 2), (["--threads", "7"], 7)])
def test_cli_overrides_file_which_overrides_defaults(tmp_path, extra, expected):
    config = config_file(tmp_path, run={"threads": 2, "min_qual": 9.0})
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nAC\n+\nII\n")
    with patch("muc_one_span.pipeline.execute_pipeline") as execute:
        result = CliRunner().invoke(
            main, ["--config", str(config), "run", "--input", str(reads), *extra]
        )
    assert result.exit_code == 0, result.output
    assert execute.call_args.args[4] == expected
    assert execute.call_args.args[6] == 9.0


@pytest.mark.parametrize(
    "text",
    [
        '{"schema_version":1,"unknown":{}}',
        '{"schema_version":1,"run":{"threads":true}}',
        '{"schema_version":1,"run":{"min_qual":NaN}}',
        '{"schema_version":1,"run":{"threads":2,"threads":3}}',
    ],
)
def test_invalid_configuration_never_enters_pipeline(tmp_path, text):
    config = tmp_path / "invalid.json"
    config.write_text(text)
    with patch("muc_one_span.pipeline.execute_pipeline") as execute:
        result = CliRunner().invoke(main, ["--config", str(config), "run", "--input", "unused"])
    assert result.exit_code == 2
    assert "configuration" in result.output.lower()
    execute.assert_not_called()


def test_effective_settings_survive_tool_failure(tmp_path):
    config = config_file(tmp_path, run={"threads": 2, "min_qual": 9.0})
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nAC\n+\nII\n")
    out = tmp_path / "result"
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("missing tool")):
        result = CliRunner().invoke(
            main,
            [
                "--config",
                str(config),
                "run",
                "--input",
                str(reads),
                "--output-dir",
                str(out),
                "--threads",
                "7",
            ],
        )
    assert result.exit_code != 0
    recorded = json.loads((out / "run_configuration.json").read_text())
    assert recorded["settings"]["run"]["threads"] == 7
    assert recorded["settings"]["run"]["min_qual"] == 9.0
    assert recorded["settings"]["calling"]["read_phase"] is False
    assert len(recorded["configuration_sha256"]) == 64
    assert len(recorded["input_sha256"]) == 64
    assert len(recorded["repeat_dictionary_sha256"]) == 64
    assert json.loads((out / "run_status.json").read_text())["status"] == "execution_failed"


def test_classify_uses_configured_probe_settings(tmp_path):
    config = config_file(tmp_path, classification={"max_indel_probe": 12})
    sequence = tmp_path / "sequence.fa"
    sequence.write_text(">s\nAC\n")
    with patch(
        "muc_one_span.classify.classify_sequence",
        return_value={
            "structure": "X",
            "mutations_detected": [],
        },
    ) as classify:
        result = CliRunner().invoke(
            main,
            [
                "--config",
                str(config),
                "classify",
                "--input",
                str(sequence),
                "--output-dir",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    assert classify.call_args.kwargs["settings"].max_indel_probe == 12


@pytest.mark.parametrize(
    "command,options,target",
    [
        ("map", ["--threads", "0"], "muc_one_span.mapping.map_reads"),
        (
            "ladder",
            ["--min-units", "20", "--max-units", "10"],
            "muc_one_span.ladder.generate_ladder_fasta",
        ),
    ],
)
def test_explicit_invalid_stage_values_are_rejected(tmp_path, command, options, target):
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nAC\n+\nII\n")
    args = [command, *options]
    if command == "map":
        args += ["--input", str(reads)]
    with patch(target) as operation, patch("muc_one_span.tools.check_tools"):
        result = CliRunner().invoke(main, args)
    assert result.exit_code == 2, result.output
    operation.assert_not_called()


def test_failed_rerun_does_not_retain_old_configuration(tmp_path):
    old = tmp_path / "run_configuration.json"
    old.write_text('{"input_path":"previous-input"}')
    result = CliRunner().invoke(
        main, ["run", "--input", str(tmp_path / "missing"), "-o", str(tmp_path)]
    )
    assert result.exit_code == 2
    assert not old.exists()


def test_programmatic_defaults_survive_without_config(tmp_path):
    with patch("muc_one_span.pipeline.execute_pipeline") as execute:
        result = CliRunner().invoke(
            main,
            ["run", "--input", "unused", "-o", str(tmp_path)],
            default_map={"run": {"threads": 11}},
        )
    assert result.exit_code == 0, result.output
    assert execute.call_args.args[4] == 11


def test_ladder_range_alone_does_not_require_custom_run_reference(tmp_path):
    config = config_file(tmp_path, reference_layout={"max_units": 100})
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nAC\n+\nII\n")
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("tool sentinel")):
        result = CliRunner().invoke(
            main, ["--config", str(config), "run", "--input", str(reads), "-o", str(tmp_path)]
        )
    assert str(result.exception) == "tool sentinel"


def test_standalone_consensus_passes_configured_dictionary(tmp_path):
    from muc_one_span.config import load_repeat_dictionary

    dictionary = load_repeat_dictionary()
    custom = tmp_path / "custom.json"
    custom.write_text("{}")
    config = config_file(tmp_path, repeat_dictionary=str(custom))
    source = tmp_path / "source"
    source.write_text("{}")
    with (
        patch("muc_one_span.config.load_repeat_dictionary", return_value=dictionary) as load,
        patch("muc_one_span.tools.check_tools"),
        patch("muc_one_span.consensus.build_consensus_per_allele", return_value={}) as build,
    ):
        result = CliRunner().invoke(
            main,
            [
                "--config",
                str(config),
                "consensus",
                "--input",
                str(source),
                "--reference",
                str(source),
                "--alleles-json",
                str(source),
                "-o",
                str(tmp_path),
            ],
        )
    assert result.exit_code == 0, result.output
    load.assert_called_once_with(custom)
    assert build.call_args.kwargs["repeat_dict"] is dictionary


@pytest.mark.parametrize(
    "option,value", [("--threads", "0"), ("--min-coverage", "0"), ("--min-qual", "nan")]
)
def test_invalid_run_override_is_usage_error(tmp_path, option, value):
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nAC\n+\nII\n")
    with patch("muc_one_span.mapping.map_reads") as mapping:
        result = CliRunner().invoke(
            main, ["run", "-i", str(reads), "-o", str(tmp_path), option, value]
        )
    assert result.exit_code == 2, result.output
    assert "Invalid value" in result.output
    mapping.assert_not_called()
    assert json.loads((tmp_path / "run_status.json").read_text())["error_type"] == "BadParameter"


def test_recorded_auto_preset_remains_auto_when_reused(tmp_path):
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nAC\n+\nII\n")
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("tool sentinel")):
        CliRunner().invoke(
            main, ["run", "-i", str(reads), "-o", str(tmp_path), "--engine", "ladder"]
        )
    first = json.loads((tmp_path / "run_configuration.json").read_text())
    config = tmp_path / "reused.json"
    config.write_text(json.dumps(first["settings"]))
    out = tmp_path / "ont"
    with patch("muc_one_span.tools.check_tools", side_effect=RuntimeError("tool sentinel")):
        CliRunner().invoke(
            main,
            [
                *("--config", str(config), "run", "-i", str(reads), "-o", str(out)),
                *("--platform", "ont", "--engine", "ladder"),
            ],
        )
    second = json.loads((out / "run_configuration.json").read_text())
    assert second["settings"]["run"]["minimap2_preset"] is None
    assert second["resolved_minimap2_preset"] == "lr:hq"


def test_invalid_global_config_leaves_previous_attempt_untouched(tmp_path):
    config = tmp_path / "invalid.json"
    config.write_text('{"bad":true}')
    status = tmp_path / "run_status.json"
    previous = '{"schema_version":1,"status":"completed"}'
    status.write_text(previous)
    result = CliRunner().invoke(
        main, ["--config", str(config), "run", "-i", "missing", "-o", str(tmp_path)]
    )
    assert result.exit_code == 2
    assert status.read_text() == previous


def test_dictionary_flank_mismatch_is_usage_error_before_tools(tmp_path):
    from dataclasses import replace

    from muc_one_span.config import load_repeat_dictionary

    rd = replace(load_repeat_dictionary(), flanking_left="AAAA")
    reads = tmp_path / "reads.fq"
    reads.write_text("@r\nAC\n+\nII\n")
    with (
        patch("muc_one_span.config.load_repeat_dictionary", return_value=rd),
        patch("muc_one_span.tools.check_tools") as check,
    ):
        result = CliRunner().invoke(main, ["run", "-i", str(reads), "-o", str(tmp_path)])
    assert result.exit_code == 2, result.output
    assert "flank_length" in result.output
    check.assert_not_called()
