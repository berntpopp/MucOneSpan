"""Compatibility tests for rendering legacy and strict batch reports."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def module():
    path = Path(__file__).parents[2] / "scripts" / "validate_catalog.py"
    spec = importlib.util.spec_from_file_location("validate_catalog_script", path)
    assert spec and spec.loader
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_versioned_batch_report_renders_strict_metrics(tmp_path: Path) -> None:
    report = {
        "schema_version": 1,
        "totals": {
            "samples": 1,
            "statuses": {"completed": 1},
            "metrics": {
                "event_tp": {"min": 0, "max": 0},
                "event_fn": {"min": 1, "max": 1},
                "event_fp": {"min": 1, "max": 1},
            },
            "event_recall": {"numerator": 0, "denominator": 1, "value": 0.0},
            "event_precision": {"numerator": 0, "denominator": 1, "value": 0.0},
            "normal_specificity": {"numerator": 0, "denominator": 0, "value": None},
        },
        "samples": [
            {
                "sample": "wrong-site",
                "status": "completed",
                "truth_events": 1,
                "predicted_alleles": 1,
                "truth_haplotypes": 1,
                "missing_alleles": 0,
                "extra_alleles": 0,
                "prediction_events": {
                    "allele_1": [{"repeat_index": 999, "parent": "X", "name": "dupC"}]
                },
                "metrics": {
                    "event_tp": {"min": 0, "max": 0},
                    "event_fn": {"min": 1, "max": 1},
                    "event_fp": {"min": 1, "max": 1},
                },
                "alternatives": [
                    {
                        "pairs": [
                            {
                                "prediction": "allele_1",
                                "called_repeats": 60,
                            }
                        ]
                    }
                ],
            }
        ],
    }
    (tmp_path / "batch_results.json").write_text(json.dumps(report))
    output = tmp_path / "catalog.tsv"
    assert module().main(["--results-dir", str(tmp_path), "--output", str(output)]) == 0
    text = output.read_text()
    assert "wrong-site\tcompleted\t1\t1\t0..0\t1..1\t1..1" in text
    assert "allele_1:999:X:dupC" in text
    assert "# Event annotation: 0 TP, 1 FN, 1 FP" in text
    assert "# Event recall: 0/1 (0.0%)" in text
    assert "# Normal specificity: 0/0 (undefined)" in text


def test_legacy_batch_array_remains_readable(tmp_path: Path) -> None:
    legacy = [
        {
            "sample": "old",
            "expected_mutation": "dupC",
            "status": "TP",
            "details": {
                "detected_mutations": [{"mutation_name": "dupC"}],
                "allele_1_structure_len": 60,
            },
        }
    ]
    (tmp_path / "batch_results.json").write_text(json.dumps(legacy))
    output = tmp_path / "catalog.tsv"
    assert module().main(["--results-dir", str(tmp_path), "--output", str(output)]) == 0
    text = output.read_text()
    assert "old\tdupC\tTP\tdupC\t60\t" in text
    assert "# Summary: 1 samples, 1 TP, 0 TN, 0 FP, 0 FN" in text


def test_unknown_version_returns_nonzero(tmp_path: Path) -> None:
    (tmp_path / "batch_results.json").write_text('{"schema_version":99,"samples":[]}')
    assert module().main(["--results-dir", str(tmp_path)]) == 2


def test_invalid_truth_is_not_rendered_as_zero_error_metrics(tmp_path: Path) -> None:
    report = {
        "schema_version": 1,
        "totals": {"samples": 1, "statuses": {"invalid_truth": 1}, "metrics": {}},
        "samples": [{"sample": "invalid", "status": "invalid_truth"}],
    }
    (tmp_path / "batch_results.json").write_text(json.dumps(report))
    output = tmp_path / "catalog.tsv"
    assert module().main(["--results-dir", str(tmp_path), "--output", str(output)]) == 0
    row = output.read_text().splitlines()[1].split("\t")
    assert row[4:7] == ["not_assessable"] * 3
