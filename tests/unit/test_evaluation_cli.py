"""CLI keeps explicit inventoried failures visible and returns nonzero."""

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

from muc_one_span.evaluation.models import TruthHaplotype, TruthSample


def module():
    path = Path(__file__).parents[2] / "scripts" / "evaluate.py"
    spec = importlib.util.spec_from_file_location("evaluation_script", path)
    assert spec and spec.loader
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_cli_accounting_and_output_parent(tmp_path):
    inventory = tmp_path / "expected.json"
    inventory.write_text('["good", "missing"]')
    root = tmp_path / "results"
    root.mkdir()
    good = root / "good"
    good.mkdir()
    (good / "summary.json").write_text(
        json.dumps(
            {
                "alleles": {"p": {"length": 1, "canonical_repeats": -8}},
                "classifications": {"p": {"structure": "X", "mutations": []}},
            }
        )
    )
    (good / "consensus_p.fa").write_text(">p\nA\n")
    truth = TruthSample("good", (TruthHaplotype("h", "C", ("X",)),))
    output = tmp_path / "new" / "report.json"
    script = module()
    with patch.object(script, "load_truth", return_value=truth):
        exit_code = script.main(
            [
                str(root),
                "--truth-dir",
                str(tmp_path),
                "--expected-samples",
                str(inventory),
                "--output",
                str(output),
            ]
        )
    report = json.loads(output.read_text())
    assert exit_code == 1
    assert report["schema_version"] == 1
    assert report["totals"]["samples"] == 2
    assert report["totals"]["truth_haplotypes"] == 2
    assert report["totals"]["sequence_accuracy"]["value"] == 0
    assert [s["sample"] for s in report["samples"]] == ["good", "missing"]
    good_row, missing_row = report["samples"]
    assert good_row["clinical"]["reasons"] == []
    assert good_row["reconstruction_flags"] == []
    assert missing_row["status"] == "not_attempted"
    assert missing_row["reconstruction_flags"] == ["not_attempted", "missing_allele"]


def test_cli_invalid_truth_does_not_hide_other_rows(tmp_path):
    inventory = tmp_path / "expected.json"
    inventory.write_text('["a", "b"]')
    output = tmp_path / "report.json"
    code = module().main(
        [
            str(tmp_path),
            "--truth-root",
            str(tmp_path),
            "--expected-samples",
            str(inventory),
            "--output",
            str(output),
        ]
    )
    report = json.loads(output.read_text())
    assert code == 1
    assert report["totals"]["invalid_truth_samples"] == 2


def test_cli_bad_inventory_writes_error_report(tmp_path):
    inventory = tmp_path / "expected.json"
    inventory.write_text("[]")
    output = tmp_path / "report.json"
    assert (
        module().main(
            [str(tmp_path), "--expected-samples", str(inventory), "--output", str(output)]
        )
        == 2
    )
    assert json.loads(output.read_text())["input_error"]
