"""Wave 1 execution provenance and reporting interfaces through real Click parsing."""

import json
from contextlib import ExitStack
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main
from muc_one_span.settings import RunSettings


@pytest.mark.parametrize("value", [0, -1, float("inf"), float("nan"), True, "30"])
def test_mapping_timeout_rejects_invalid_values(value):
    with pytest.raises(ValueError, match="mapping_timeout"):
        RunSettings(mapping_timeout=value)


@pytest.mark.parametrize("command", ["map", "run"])
@pytest.mark.parametrize("override,expected", [([], 45.0), (["--mapping-timeout", "90"], 90.0)])
def test_mapping_timeout_cli_precedence(tmp_path, command, override, expected):
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"schema_version": 1, "run": {"mapping_timeout": 45}}))
    reads = tmp_path / "reads.fq"
    reads.touch()
    target = (
        "muc_one_span.mapping.map_reads"
        if command == "map"
        else "muc_one_span.pipeline.execute_pipeline"
    )
    with patch(target) as operation, patch("muc_one_span.tools.check_tools"):
        result = CliRunner().invoke(
            main,
            ["--config", str(config), command, "-i", str(reads), "-o", str(tmp_path), *override],
        )
    assert result.exit_code == 0, result.output
    key = "timeout" if command == "map" else "mapping_timeout"
    assert operation.call_args.kwargs[key] == expected


def summary_file(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    path = source / "summary.json"
    path.write_text(
        json.dumps(
            {
                "alleles": {"allele_1": {"length": 60, "reads": 100}},
                "classifications": {"allele_1": {"mutations": []}},
                "run_status": {"status": "completed"},
            }
        )
    )
    return path


@pytest.mark.parametrize(
    "status", ["insufficient_evidence", "execution_failed", "interrupted", "running"]
)
def test_standalone_report_uses_source_sidecar_over_stale_summary(tmp_path, status):
    source = summary_file(tmp_path)
    (source.parent / "run_status.json").write_text(json.dumps({"status": status}))
    output = tmp_path / "relocated.html"
    result = CliRunner().invoke(main, ["report", "-i", str(source), "-o", str(output)])
    assert result.exit_code == 0, result.output
    html = output.read_text()
    assert "decision-inconclusive" in html.split("<body", 1)[1]
    assert 'class="decision-banner decision-negative"' not in html


def test_standalone_never_infers_status_from_report_destination(tmp_path):
    source = summary_file(tmp_path)
    (tmp_path / "run_status.json").write_text('{"status":"execution_failed"}')
    output = tmp_path / "relocated.html"
    result = CliRunner().invoke(main, ["report", "-i", str(source), "-o", str(output)])
    assert result.exit_code == 0, result.output
    assert "decision-banner decision-negative" in output.read_text()


def test_standalone_explicit_status_and_all_vcfs(tmp_path):
    source = summary_file(tmp_path)
    status = tmp_path / "failed.json"
    status.write_text('{"status":"interrupted"}')
    paths = [tmp_path / f"{label}.vcf" for label in ("a", "b", "legacy")]
    for path in paths:
        path.touch()
    with patch("muc_one_span.report.generate_report") as generate:
        result = CliRunner().invoke(
            main,
            [
                "report",
                "-i",
                str(source),
                "--run-status",
                str(status),
                "--allele-vcf",
                "allele_1",
                str(paths[0]),
                "--allele-vcf",
                "allele_2",
                str(paths[1]),
                "--vcf",
                str(paths[2]),
            ],
        )
    assert result.exit_code == 0, result.output
    assert generate.call_args.kwargs["vcf_paths"] == dict(
        zip(("allele_1", "allele_2"), paths[:2], strict=True)
    )
    assert generate.call_args.kwargs["execution_status"]["status"] == "interrupted"


@pytest.fixture
def pipeline_stages(tmp_path):
    reads = tmp_path / "reads.fq"
    reads.touch()
    fa = tmp_path / "consensus.fa"
    fa.write_text(">a\nACGT\n")
    vcfs = {key: tmp_path / f"{key}.vcf" for key in ("allele_1", "allele_2")}
    alleles = {
        key: {"length": 60, "reads": 100, "contig_name": f"contig_{index}"}
        for index, key in enumerate(vcfs, 1)
    }
    classification = {
        "structure": "X",
        "mutations_detected": [],
        "reconstruction_status": "resolved",
    }
    mocks = {
        "tools.check_tools": None,
        "tools.get_tool_versions": {},
        "mapping.map_reads": tmp_path / "mapping.bam",
        "mapping.get_idxstats": "",
        "alleles.parse_idxstats": {},
        "alleles.detect_alleles": alleles,
        "calling.call_variants_per_allele": vcfs,
        "consensus.build_consensus_per_allele": dict.fromkeys(vcfs, fa),
        "classify.classify_sequence": classification,
        "classify.validate_mutations_against_vcf": classification,
        "vcf.parse_vcf_variants": [],
    }
    with ExitStack() as stack:
        for target, value in mocks.items():
            stack.enter_context(patch(f"muc_one_span.{target}", return_value=value))
        yield reads, vcfs


@pytest.mark.parametrize(
    "failure", [None, RuntimeError("render failed"), ImportError("report unavailable")]
)
def test_full_run_status_commits_only_after_report_success(tmp_path, pipeline_stages, failure):
    reads, vcfs = pipeline_stages
    output = tmp_path / "out"
    from muc_one_span.report import generate_report

    def render(summary, path, **kwargs):
        assert json.loads((output / "run_status.json").read_text())["status"] == "running"
        assert kwargs["execution_status"]["status"] == "analysis_completed"
        assert kwargs["vcf_paths"] == vcfs
        if failure:
            raise failure
        return generate_report(summary, path, **kwargs)

    with patch("muc_one_span.report.generate_report", side_effect=render):
        result = CliRunner().invoke(main, ["run", "-i", str(reads), "-o", str(output), "--report"])
    status = json.loads((output / "run_status.json").read_text())
    if failure:
        assert result.exit_code != 0
        assert status["status"] == "execution_failed"
        assert status["error_type"] == type(failure).__name__
        assert not (output / "report.html").exists()
    else:
        assert result.exit_code == 0, str(result.exception)
        assert status["status"] == "completed"
        assert (
            json.loads((output / "summary.json").read_text())["run_status"]["status"] == "completed"
        )
        assert "decision-banner decision-negative" in (output / "report.html").read_text()


def test_standalone_uses_absolute_consensus_vcf_provenance(tmp_path):
    source = summary_file(tmp_path)
    summary = json.loads(source.read_text())
    vcf = tmp_path / "actual.vcf"
    vcf.touch()
    summary["alleles"]["allele_1"].update(
        vcf_path="old-working-directory/results/allele_1/variants.vcf.gz",
        consensus_context={"vcf_path": str(vcf)},
    )
    source.write_text(json.dumps(summary))
    with patch("muc_one_span.report.generate_report") as generate:
        result = CliRunner().invoke(main, ["report", "-i", str(source)])
    assert result.exit_code == 0, result.output
    assert generate.call_args.kwargs["vcf_paths"] == {"allele_1": vcf}


@pytest.mark.parametrize("value", ["{bad", "null", "[]"])
def test_invalid_execution_sidecar_cannot_fall_back_to_negative(tmp_path, value):
    source = summary_file(tmp_path)
    (source.parent / "run_status.json").write_text(value)
    output = tmp_path / "report.html"
    result = CliRunner().invoke(main, ["report", "-i", str(source), "-o", str(output)])
    assert result.exit_code != 0
    assert "execution status" in result.output
    assert not output.exists()


def test_invalid_recorded_clinical_decision_section_is_a_clean_error(tmp_path):
    """M2: a bad recorded clinical_decision section must not raise a raw traceback."""
    source = summary_file(tmp_path)
    summary = json.loads(source.read_text())
    summary["configuration"] = {"settings": {"clinical_decision": {"bogus_field": 1}}}
    source.write_text(json.dumps(summary))
    output = tmp_path / "report.html"
    result = CliRunner().invoke(main, ["report", "-i", str(source), "-o", str(output)])
    assert result.exit_code != 0
    assert not isinstance(result.exception, ValueError)
    assert "clinical_decision" in result.output
    assert "bogus_field" in result.output
    assert "Traceback" not in result.output
    assert not output.exists()
