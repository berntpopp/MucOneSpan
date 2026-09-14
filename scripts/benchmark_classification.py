#!/usr/bin/env python3
"""Compare scalar and bit-vector scoring inside the sequence classifier."""

from __future__ import annotations

import argparse
import ast
import hashlib
import inspect
import json
import platform
import sys
import textwrap
from collections.abc import Callable
from pathlib import Path
from statistics import median
from time import perf_counter
from typing import Any, cast

import muc_one_span.classify as classifier
from muc_one_span.config import RepeatDictionary, load_repeat_dictionary
from muc_one_span.repeat_alignment import _edit_distance_bitvector

Scorer = Callable[[str, str], int]
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "tests/results/deep_validation_20260914/hifi" / "sample_bench_5003/consensus_allele_2.fa"
)
DEFAULT_PANEL_ROOT = ROOT / "tests/results/deep_validation_20260914"
DEFAULT_OUTPUT = ROOT / "results/production_validation_20260914/scoring/benchmark.json"
BASELINE_COMMIT = "d8390b3c244ef8f3240af74b92db12b50dfc77d1"
BASELINE_REPEAT_ALIGNMENT_SHA256 = (
    "4e9fcdfdea4e2ed6a0c0ab7ed69c7e434a29fa050c6f4d36cd51f6575b222ccf"
)
BASELINE_SCALAR_AST_SHA256 = "2c2eb12bd7c52c46fc62b7e19a20b4bf1db8f7edf4760fc3ec591307b1ce5ac2"


def _scalar_edit_distance(s1: str, s2: str) -> int:
    """Compute Levenshtein edit distance between two sequences.

    Args:
        s1: First sequence.
        s2: Second sequence.

    Returns:
        Minimum number of single-character edits (insert, delete, substitute).
    """
    m, n = len(s1), len(s2)

    # Use single-row optimization for memory efficiency
    prev = list(range(n + 1))
    curr = [0] * (n + 1)

    for i in range(1, m + 1):
        curr[0] = i
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                curr[j] = prev[j - 1]
            else:
                curr[j] = 1 + min(prev[j], curr[j - 1], prev[j - 1])
        prev, curr = curr, prev

    return prev[n]


def _normalized_function_hash(function: Scorer) -> str:
    """Hash a function AST after normalizing its local symbol name."""
    tree = ast.parse(textwrap.dedent(inspect.getsource(function)))
    node = tree.body[0]
    if not isinstance(node, ast.FunctionDef):
        raise TypeError("scalar baseline source is not a function definition")
    node.name = "edit_distance"
    return hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()


def _classify_with_scorer(
    sequence: str,
    repeat_dict: RepeatDictionary,
    scorer: Scorer,
) -> dict[str, Any]:
    """Classify once with a process-local scorer and always restore the original."""
    original = classifier.edit_distance
    mutable_classifier = cast(Any, classifier)
    try:
        mutable_classifier.edit_distance = scorer
        return classifier.classify_sequence(sequence, repeat_dict)
    finally:
        mutable_classifier.edit_distance = original


def _run_batch(
    sequence: str,
    repeat_dict: RepeatDictionary,
    scorer: Scorer,
    batch_size: int,
) -> tuple[float, dict[str, Any], bool]:
    """Time repeated classifications with one scorer patch per batch."""
    original_scorer = classifier.edit_distance
    mutable_classifier = cast(Any, classifier)
    first_result: dict[str, Any] | None = None
    equal_within_batch = True
    try:
        mutable_classifier.edit_distance = scorer
        start = perf_counter()
        for _ in range(batch_size):
            result = classifier.classify_sequence(sequence, repeat_dict)
            if first_result is None:
                first_result = result
            else:
                equal_within_batch = equal_within_batch and result == first_result
        elapsed = perf_counter() - start
    finally:
        mutable_classifier.edit_distance = original_scorer
    if first_result is None:
        raise ValueError("batch_size must be positive")
    return elapsed, first_result, equal_within_batch


def _read_fasta(path: Path) -> str:
    """Read one FASTA sequence, rejecting missing or empty benchmark input."""
    if not path.is_file():
        raise FileNotFoundError(f"benchmark FASTA not found: {path}")
    sequence = "".join(
        line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line[:1] != ">"
    )
    if not sequence:
        raise ValueError(f"benchmark FASTA contains no sequence: {path}")
    return sequence


def _find_panel_fastas(root: Path) -> list[Path]:
    """Find trimmed allele consensus FASTAs across cached platform directories."""
    return sorted(root.glob("*/*/consensus_allele_?.fa"))


def _near_exact_control(repeat_dict: RepeatDictionary) -> str:
    """Build a long classification input with one novel substitution."""
    repeat = repeat_dict.repeats[repeat_dict.canonical_repeat]
    for offset in range(len(repeat) // 2, len(repeat)):
        for replacement in "ACGT":
            candidate = repeat[:offset] + replacement + repeat[offset + 1 :]
            if (
                replacement != repeat[offset]
                and candidate not in repeat_dict.seq_to_id
                and candidate not in repeat_dict.mutated_sequences
            ):
                return repeat * 60 + candidate + repeat * 59
    raise ValueError("could not construct a novel near-exact control")


def _sha256_file(path: Path) -> str:
    """Return a reproducibility hash for one input or source file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _result_hash(result: dict[str, Any]) -> str:
    """Hash a complete classification dictionary in canonical JSON form."""
    payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _benchmark_case(
    name: str,
    sequence: str,
    repeat_dict: RepeatDictionary,
    repetitions: int,
    batch_size: int,
) -> dict[str, Any]:
    """Run alternating paired scalar/bit-vector classifications."""
    scorers = {
        "scalar": _scalar_edit_distance,
        "bitvector": _edit_distance_bitvector,
    }
    rows: list[dict[str, Any]] = []
    reference: dict[str, Any] | None = None
    for repetition in range(repetitions):
        order = ("scalar", "bitvector") if repetition % 2 == 0 else ("bitvector", "scalar")
        for mode in order:
            elapsed, result, stable = _run_batch(sequence, repeat_dict, scorers[mode], batch_size)
            if reference is None:
                reference = result
            equal = stable and result == reference
            row = {
                "repetition": repetition + 1,
                "order": order.index(mode) + 1,
                "mode": mode,
                "batch_size": batch_size,
                "elapsed_seconds": elapsed,
                "seconds_per_classification": elapsed / batch_size,
                "complete_dictionary_equal": equal,
                "result_sha256": _result_hash(result),
            }
            rows.append(row)
            print(
                f"{name} repetition={repetition + 1} mode={mode} "
                f"seconds_per_classification={elapsed / batch_size:.9f} equal={equal}",
                flush=True,
            )

    medians = {
        mode: median(row["seconds_per_classification"] for row in rows if row["mode"] == mode)
        for mode in scorers
    }
    return {
        "name": name,
        "sequence_length": len(sequence),
        "batch_size": batch_size,
        "rows": rows,
        "median_seconds_per_classification": medians,
        "scalar_over_bitvector_speedup": medians["scalar"] / medians["bitvector"],
        "complete_dictionary_equal": all(row["complete_dictionary_equal"] for row in rows),
    }


def _compare_panel(root: Path, repeat_dict: RepeatDictionary) -> dict[str, Any]:
    """Compare complete classifier dictionaries for every cached consensus."""
    paths = _find_panel_fastas(root)
    if not paths:
        raise ValueError(f"no cached consensus FASTAs found under: {root}")
    scorers = {
        "scalar": _scalar_edit_distance,
        "bitvector": _edit_distance_bitvector,
    }
    rows: list[dict[str, Any]] = []
    input_hashes: dict[str, str] = {}
    for index, path in enumerate(paths):
        relative_path = path.relative_to(root).as_posix()
        input_hashes[relative_path] = _sha256_file(path)
        sequence = _read_fasta(path)
        order = ("scalar", "bitvector") if index % 2 == 0 else ("bitvector", "scalar")
        results: dict[str, dict[str, Any]] = {}
        elapsed: dict[str, float] = {}
        for mode in order:
            start = perf_counter()
            results[mode] = _classify_with_scorer(sequence, repeat_dict, scorers[mode])
            elapsed[mode] = perf_counter() - start
        equal = results["scalar"] == results["bitvector"]
        rows.append(
            {
                "path": relative_path,
                "sequence_length": len(sequence),
                "order": list(order),
                "scalar_seconds": elapsed["scalar"],
                "bitvector_seconds": elapsed["bitvector"],
                "complete_dictionary_equal": equal,
                "scalar_result_sha256": _result_hash(results["scalar"]),
                "bitvector_result_sha256": _result_hash(results["bitvector"]),
            }
        )
        print(f"panel {relative_path} equal={equal}", flush=True)
    return {
        "root": str(root),
        "file_count": len(paths),
        "rows": rows,
        "input_hashes": input_hashes,
        "complete_dictionary_equal": all(row["complete_dictionary_equal"] for row in rows),
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--panel-root", type=Path, default=DEFAULT_PANEL_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--difficult-batch", type=int, default=1)
    parser.add_argument("--control-batch", type=int, default=1000)
    parser.add_argument("--near-exact-batch", type=int, default=100)
    return parser.parse_args()


def main() -> int:
    """Run paired benchmarks, write evidence, and enforce fixed performance gates."""
    args = _parse_args()
    if args.repetitions < 3:
        raise ValueError("at least three paired repetitions are required")
    if args.difficult_batch < 1 or args.control_batch < 1 or args.near_exact_batch < 1:
        raise ValueError("batch sizes must be positive")

    scalar_hash = _normalized_function_hash(_scalar_edit_distance)
    if scalar_hash != BASELINE_SCALAR_AST_SHA256:
        raise RuntimeError("embedded scalar scorer differs from the baseline recurrence")

    repeat_dict = load_repeat_dictionary()
    difficult_sequence = _read_fasta(args.input)
    exact_control = repeat_dict.repeats[repeat_dict.canonical_repeat] * 120
    near_exact_control = _near_exact_control(repeat_dict)
    cases = [
        _benchmark_case(
            "difficult_bench5003",
            difficult_sequence,
            repeat_dict,
            args.repetitions,
            args.difficult_batch,
        ),
        _benchmark_case(
            "exact_X_120",
            exact_control,
            repeat_dict,
            args.repetitions,
            args.control_batch,
        ),
        _benchmark_case(
            "near_exact_X_120",
            near_exact_control,
            repeat_dict,
            args.repetitions,
            args.near_exact_batch,
        ),
    ]
    panel = _compare_panel(args.panel_root, repeat_dict)
    cases_by_name = {case["name"]: case for case in cases}
    difficult = cases_by_name["difficult_bench5003"]
    control = cases_by_name["exact_X_120"]
    near_exact = cases_by_name["near_exact_X_120"]
    control_medians = control["median_seconds_per_classification"]
    near_exact_medians = near_exact["median_seconds_per_classification"]
    acceptance = {
        "zero_complete_dictionary_differences": all(
            case["complete_dictionary_equal"] for case in cases
        )
        and panel["complete_dictionary_equal"],
        "difficult_speedup_at_least_2x": difficult["scalar_over_bitvector_speedup"] >= 2.0,
        "exact_control_regression_at_most_20_percent": (
            control_medians["bitvector"] <= control_medians["scalar"] * 1.2
        ),
        "near_exact_control_regression_at_most_20_percent": (
            near_exact_medians["bitvector"] <= near_exact_medians["scalar"] * 1.2
        ),
    }
    repeat_source = ROOT / "src/muc_one_span/data/repeats/repeats.json"
    source_paths = {
        "repeat_alignment": ROOT / "src/muc_one_span/repeat_alignment.py",
        "classify": ROOT / "src/muc_one_span/classify.py",
        "benchmark": Path(__file__).resolve(),
    }
    report = {
        "schema_version": 1,
        "conditions": {
            "repetitions": args.repetitions,
            "order": "alternating scalar/bitvector in one Python process",
            "python": sys.version,
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "hashes": {
            "sources": {name: _sha256_file(path) for name, path in source_paths.items()},
            "inputs": {
                "difficult_fasta": _sha256_file(args.input),
                "repeat_dictionary": _sha256_file(repeat_source),
            },
        },
        "scalar_baseline": {
            "source_commit": BASELINE_COMMIT,
            "repeat_alignment_sha256": BASELINE_REPEAT_ALIGNMENT_SHA256,
            "embedded_function_ast_sha256": scalar_hash,
            "expected_function_ast_sha256": BASELINE_SCALAR_AST_SHA256,
            "hash_verified": scalar_hash == BASELINE_SCALAR_AST_SHA256,
        },
        "cases": cases,
        "cached_consensus_panel": panel,
        "acceptance": acceptance,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(acceptance, indent=2), flush=True)
    print(f"wrote {args.output}", flush=True)
    return 0 if all(acceptance.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
