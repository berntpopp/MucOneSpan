"""Exactness and traceback tests for repeat alignment."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest

import muc_one_span.repeat_alignment as repeat_alignment
from muc_one_span import classify
from muc_one_span.config import load_repeat_dictionary


def _scalar_edit_distance(left: str, right: str) -> int:
    """Independent full-matrix Levenshtein oracle used only by these tests."""
    rows = [[0] * (len(right) + 1) for _ in range(len(left) + 1)]
    for row in range(len(left) + 1):
        rows[row][0] = row
    for column in range(len(right) + 1):
        rows[0][column] = column
    for row, left_symbol in enumerate(left, start=1):
        for column, right_symbol in enumerate(right, start=1):
            rows[row][column] = min(
                rows[row - 1][column] + 1,
                rows[row][column - 1] + 1,
                rows[row - 1][column - 1] + (left_symbol != right_symbol),
            )
    return rows[-1][-1]


def _mutate(sequence: str, rng: random.Random, alphabet: str) -> str:
    """Apply deterministic substitutions and indels to a source sequence."""
    result = list(sequence)
    for _ in range(rng.randrange(4)):
        operation = rng.choice(("delete", "insert", "substitute"))
        if operation == "delete" and result:
            del result[rng.randrange(len(result))]
        elif operation == "insert":
            result.insert(rng.randrange(len(result) + 1), rng.choice(alphabet))
        elif result:
            result[rng.randrange(len(result))] = rng.choice(alphabet)
    return "".join(result)


def _seeded_pairs(count: int = 10_000) -> list[tuple[str, str]]:
    """Build a reproducible mix of random and adversarial sequence pairs."""
    rng = random.Random(20260914)
    alphabet = "ACGTNRYKMSWBDHVacgti"
    boundaries = (0, 1, 31, 32, 63, 64, 65, 127, 128, 129)
    pairs: list[tuple[str, str]] = [
        ("", ""),
        ("", "ACGT"),
        ("ACGT", ""),
        ("A" * 129, "A" * 64),
        ("NRYKMSWBDHV", "nrykmswbdhv"),
        ("ACGT" * 32, "TGCA" * 32),
    ]
    while len(pairs) < count:
        index = len(pairs)
        left_length = boundaries[index % len(boundaries)] if index < 406 else rng.randrange(73)
        if index % 7 == 0:
            left = rng.choice(alphabet) * left_length
        else:
            left = "".join(rng.choices(alphabet, k=left_length))
        right = _mutate(left, rng, alphabet)
        if index % 11 == 0:
            right = "".join(rng.choices(alphabet, k=rng.choice(boundaries)))
        pairs.append((left, right))
    return pairs


def test_public_edit_distance_matches_scalar_oracle_on_10k_seeded_pairs() -> None:
    """Literal edit distance agrees with an independent scalar recurrence."""
    for left, right in _seeded_pairs():
        expected = _scalar_edit_distance(left, right)
        assert repeat_alignment.edit_distance(left, right) == expected, (left, right)
        assert repeat_alignment.edit_distance(right, left) == expected, (right, left)


def test_bitvector_entry_point_handles_machine_word_boundaries() -> None:
    """The optimized recurrence remains exact across integer word boundaries."""
    optimized = repeat_alignment._edit_distance_bitvector
    for left, right in _seeded_pairs(406):
        assert optimized(left, right) == _scalar_edit_distance(left, right)


def test_edit_distance_uses_literal_symbol_equality() -> None:
    """Case and ambiguity-code symbols do not receive wildcard treatment."""
    assert repeat_alignment.edit_distance("ACGTNI", "acgtnA") == 6


def test_edit_distance_is_symmetric_above_five_kilobases() -> None:
    """Unbounded Python bit vectors remain exact for whole-VNTR-sized inputs."""
    left = "A" * 6000
    right = left[:3000] + "C" + left[3000:]
    assert repeat_alignment.edit_distance(left, right) == 1
    assert repeat_alignment.edit_distance(right, left) == 1


def test_traceback_preserves_first_hit_ties_and_insertion_position() -> None:
    """Traceback retains substitution-first ties and existing coordinates."""
    assert repeat_alignment.characterize_differences("ACGT", "ACAGT") == [
        {"pos": 3, "ref": "", "alt": "A", "type": "insertion"}
    ]
    assert repeat_alignment.characterize_differences("AA", "A") == [
        {"pos": 1, "ref": "A", "alt": "", "type": "deletion"}
    ]
    assert repeat_alignment.characterize_differences("A", "AA") == [
        {"pos": 1, "ref": "", "alt": "A", "type": "insertion"}
    ]
    assert repeat_alignment.characterize_differences("AG", "GA") == [
        {"pos": 1, "ref": "A", "alt": "G", "type": "substitution"},
        {"pos": 2, "ref": "G", "alt": "A", "type": "substitution"},
    ]


def test_benchmark_restores_classifier_scorer_after_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The paired harness cannot leak a patched scorer after classification fails."""
    script_path = Path(__file__).parents[2] / "scripts" / "benchmark_classification.py"
    assert script_path.exists(), "classification benchmark harness is missing"
    spec = importlib.util.spec_from_file_location("benchmark_classification", script_path)
    assert spec is not None and spec.loader is not None
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)

    original = classify.edit_distance

    def fail_classification(sequence: str, repeat_dict: object) -> dict:
        raise RuntimeError("controlled failure")

    monkeypatch.setattr(classify, "classify_sequence", fail_classification)
    with pytest.raises(RuntimeError, match="controlled failure"):
        benchmark._classify_with_scorer("ACGT", object(), lambda left, right: 0)
    assert classify.edit_distance is original


def test_benchmark_discovers_complete_cached_consensus_panel(tmp_path: Path) -> None:
    """Panel discovery includes both platforms and only trimmed allele FASTAs."""
    script_path = Path(__file__).parents[2] / "scripts" / "benchmark_classification.py"
    spec = importlib.util.spec_from_file_location("benchmark_classification", script_path)
    assert spec is not None and spec.loader is not None
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)

    expected = []
    for platform_name in ("hifi", "ont"):
        sample = tmp_path / platform_name / f"sample_{platform_name}"
        sample.mkdir(parents=True)
        fasta = sample / "consensus_allele_1.fa"
        fasta.write_text(">allele\nACGT\n", encoding="utf-8")
        expected.append(fasta)
        (sample / "consensus_allele_1_full.fa").write_text(">full\nACGT\n", encoding="utf-8")

    assert benchmark._find_panel_fastas(tmp_path) == sorted(expected)


def test_benchmark_near_exact_control_exercises_scorer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The near-exact timing control reaches classifier edit-distance scoring."""
    script_path = Path(__file__).parents[2] / "scripts" / "benchmark_classification.py"
    spec = importlib.util.spec_from_file_location("benchmark_classification", script_path)
    assert spec is not None and spec.loader is not None
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    repeat_dict = load_repeat_dictionary()
    sequence = benchmark._near_exact_control(repeat_dict)

    calls = 0
    original = classify.edit_distance

    def counting_scorer(left: str, right: str) -> int:
        nonlocal calls
        calls += 1
        return original(left, right)

    monkeypatch.setattr(classify, "edit_distance", counting_scorer)
    classify.classify_sequence(sequence, repeat_dict)
    assert calls > 0
