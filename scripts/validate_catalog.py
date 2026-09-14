#!/usr/bin/env python3
"""Render legacy or strict batch_results.json as a TSV validation report."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, TextIO


def parser() -> argparse.ArgumentParser:
    """Build the catalog validation parser while preserving existing flags."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--results-dir",
        type=Path,
        default=Path("tests/results"),
        help="Directory containing batch_results.json",
    )
    result.add_argument("--output", type=Path, default=None, help="TSV path (default: stdout)")
    return result


def _legacy(results: list[dict[str, Any]], out: TextIO) -> None:
    header = (
        "sample",
        "expected_mutation",
        "status",
        "detected_mutations",
        "allele_1_len",
        "allele_2_len",
    )
    rows = []
    for entry in results:
        details = entry.get("details", {})
        detected = ", ".join(
            str(mutation.get("mutation_name", "unknown"))
            for mutation in details.get("detected_mutations", [])
        )
        rows.append(
            (
                str(entry["sample"]),
                str(entry.get("expected_mutation", "none")),
                str(entry["status"]),
                detected or "none",
                str(details.get("allele_1_structure_len", "")),
                str(details.get("allele_2_structure_len", "")),
            )
        )
    counts = {
        status: sum(entry.get("status") == status for entry in results)
        for status in ("TP", "TN", "FP", "FN")
    }
    tp, tn, fp, fn = (counts[key] for key in ("TP", "TN", "FP", "FN"))
    sensitivity = tp / (tp + fn) if tp + fn else 0
    specificity = tn / (tn + fp) if tn + fp else 0
    _write_rows(out, header, rows)
    out.write(f"\n# Summary: {len(results)} samples, {tp} TP, {tn} TN, {fp} FP, {fn} FN\n")
    out.write(f"# Sensitivity: {sensitivity:.1%}, Specificity: {specificity:.1%}\n")


def _bound(row: dict[str, Any], metric: str) -> str:
    value = row.get("metrics", {}).get(metric)
    if not isinstance(value, dict) or "min" not in value or "max" not in value:
        return "not_assessable"
    return f"{value['min']}..{value['max']}"


def _events(row: dict[str, Any]) -> tuple[int, str]:
    labels = []
    for allele, events in row.get("prediction_events", {}).items():
        for event in events:
            labels.append(
                ":".join(
                    str(value)
                    for value in (
                        allele,
                        event.get("repeat_index", "?"),
                        event.get("parent", "?"),
                        event.get("name", "?"),
                    )
                )
            )
    return len(labels), ", ".join(labels) or "none"


def _called_counts(row: dict[str, Any]) -> str:
    alternatives = row.get("alternatives", [])
    if not alternatives:
        return ""
    return ",".join(str(pair["called_repeats"]) for pair in alternatives[0].get("pairs", []))


def _ratio(label: str, ratio: dict[str, Any]) -> str:
    numerator, denominator = ratio.get("numerator", 0), ratio.get("denominator", 0)
    value = ratio.get("value")
    rendered = "undefined" if value is None else f"{value:.1%}"
    return f"# {label}: {numerator}/{denominator} ({rendered})\n"


def _strict(report: dict[str, Any], out: TextIO) -> None:
    samples = report.get("samples")
    totals = report.get("totals")
    if not isinstance(samples, list) or not isinstance(totals, dict):
        raise ValueError("schema version 1 requires samples and totals")
    header = (
        "sample",
        "status",
        "truth_events",
        "predicted_events",
        "event_tp",
        "event_fn",
        "event_fp",
        "truth_haplotypes",
        "predicted_alleles",
        "missing_alleles",
        "extra_alleles",
        "detected_events",
        "called_repeat_counts",
    )
    rows = []
    for row in samples:
        if not isinstance(row, dict):
            raise ValueError("schema version 1 sample rows must be objects")
        predicted_count, events = _events(row)
        rows.append(
            (
                str(row["sample"]),
                str(row["status"]),
                str(row.get("truth_events", "")),
                str(predicted_count),
                _bound(row, "event_tp"),
                _bound(row, "event_fn"),
                _bound(row, "event_fp"),
                str(row.get("truth_haplotypes", "")),
                str(row.get("predicted_alleles", "")),
                str(row.get("missing_alleles", "")),
                str(row.get("extra_alleles", "")),
                events,
                _called_counts(row),
            )
        )
    _write_rows(out, header, rows)
    metrics = totals.get("metrics", {})
    tp = metrics.get("event_tp", {}).get("min", 0)
    fn = metrics.get("event_fn", {}).get("max", 0)
    fp = metrics.get("event_fp", {}).get("max", 0)
    statuses = ", ".join(
        f"{key}={value}" for key, value in sorted(totals.get("statuses", {}).items())
    )
    out.write(f"\n# Strict summary: {totals.get('samples', len(samples))} samples; {statuses}\n")
    out.write(f"# Event annotation: {tp} TP, {fn} FN, {fp} FP\n")
    out.write(_ratio("Event recall", totals.get("event_recall", {})))
    out.write(_ratio("Event precision", totals.get("event_precision", {})))
    out.write(_ratio("Normal specificity", totals.get("normal_specificity", {})))


def _write_rows(out: TextIO, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(row) + "\n")


def main(argv: list[str] | None = None) -> int:
    """Render the detected schema and return nonzero for missing or invalid input."""
    args = parser().parse_args(argv)
    results_path = args.results_dir / "batch_results.json"
    if not results_path.is_file():
        print(f"Error: {results_path} not found. Run batch_analyze.py first.", file=sys.stderr)
        return 1
    try:
        document = json.loads(results_path.read_text())
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
        out = args.output.open("w") if args.output else sys.stdout
        try:
            if isinstance(document, list):
                _legacy(document, out)
            elif isinstance(document, dict) and document.get("schema_version") == 1:
                _strict(document, out)
            else:
                raise ValueError("unsupported batch_results.json schema")
        finally:
            if args.output:
                out.close()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
