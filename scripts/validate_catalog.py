#!/usr/bin/env python3
"""Validate repeat classification against MucOneUp test samples.

Reads batch_results.json (produced by batch_analyze.py) and generates
a TSV validation report with sensitivity/specificity summary.

Usage:
    python scripts/validate_catalog.py [--results-dir tests/results]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate repeat classification catalog")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("tests/results"),
        help="Directory containing batch_results.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV path (default: stdout)",
    )
    args = parser.parse_args()

    results_path = args.results_dir / "batch_results.json"
    if not results_path.exists():
        print(f"Error: {results_path} not found. Run batch_analyze.py first.", file=sys.stderr)
        sys.exit(1)

    results = json.loads(results_path.read_text())

    header = [
        "sample",
        "expected_mutation",
        "status",
        "detected_mutations",
        "allele_1_len",
        "allele_2_len",
    ]

    rows = []
    for entry in results:
        detected = ", ".join(
            m.get("mutation_name", "unknown")
            for m in entry.get("details", {}).get("detected_mutations", [])
        )
        rows.append(
            [
                entry["sample"],
                entry.get("expected_mutation", "none"),
                entry["status"],
                detected or "none",
                str(entry.get("details", {}).get("allele_1_structure_len", "")),
                str(entry.get("details", {}).get("allele_2_structure_len", "")),
            ]
        )

    tp = sum(1 for r in results if r["status"] == "TP")
    tn = sum(1 for r in results if r["status"] == "TN")
    fp = sum(1 for r in results if r["status"] == "FP")
    fn = sum(1 for r in results if r["status"] == "FN")
    total = len(results)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

    out = args.output.open("w") if args.output else sys.stdout
    try:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
        out.write(f"\n# Summary: {total} samples, {tp} TP, {tn} TN, {fp} FP, {fn} FN\n")
        out.write(f"# Sensitivity: {sensitivity:.1%}, Specificity: {specificity:.1%}\n")
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
