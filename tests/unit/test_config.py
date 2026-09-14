# tests/unit/test_config.py
"""Tests for the config module (repeat dictionary loading and classification)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from muc_one_span.config import (
    RepeatDictionary,
    classify_repeat_id,
    load_repeat_dictionary,
)


class TestLoadRepeatDictionary:
    """Tests for load_repeat_dictionary."""

    def test_loads_bundled_dictionary(self):
        """Bundled dictionary loads without error and contains >= 34 repeats."""
        rd = load_repeat_dictionary()
        assert isinstance(rd, RepeatDictionary)
        assert len(rd.repeats) >= 34
        assert rd.repeat_length_bp == 60

    def test_canonical_repeat_x_exists(self):
        """Repeat 'X' is present in the bundled dictionary with length 60."""
        rd = load_repeat_dictionary()
        assert "X" in rd.repeats
        assert len(rd.repeats["X"]) == 60

    def test_pre_repeats_present(self):
        """Repeats '1' through '5' are present in the bundled dictionary."""
        rd = load_repeat_dictionary()
        for rid in ["1", "2", "3", "4", "5"]:
            assert rid in rd.repeats, f"Pre-repeat '{rid}' missing from dictionary"

    def test_after_repeats_present(self):
        """Repeats '6' through '9' are present in the bundled dictionary."""
        rd = load_repeat_dictionary()
        for rid in ["6", "7", "8", "9"]:
            assert rid in rd.repeats, f"After-repeat '{rid}' missing from dictionary"

    def test_all_repeats_are_60bp(self):
        """Every repeat sequence in the bundled dictionary is exactly 60 bp."""
        rd = load_repeat_dictionary()
        for rid, seq in rd.repeats.items():
            assert len(seq) == 60, f"Repeat '{rid}' has length {len(seq)}, expected 60"

    def test_flanking_sequences_present(self):
        """Left and right flanking sequences are non-empty strings."""
        rd = load_repeat_dictionary()
        assert isinstance(rd.flanking_left, str) and len(rd.flanking_left) > 0
        assert isinstance(rd.flanking_right, str) and len(rd.flanking_right) > 0

    def test_loads_from_custom_path(self, tmp_path: Path):
        """load_repeat_dictionary accepts a custom Path and loads it correctly."""
        custom_data = {
            "source": "test",
            "repeat_length_bp": 60,
            "repeats": {"X": "A" * 60, "A": "C" * 60},
            "flanking_hg38": {
                "left": "LEFTSEQ",
                "right": "RIGHTSEQ",
                "vntr_region": "chr1:100-200",
            },
            "pre_repeat_ids": [],
            "after_repeat_ids": [],
            "canonical_repeat": "X",
        }
        custom_file = tmp_path / "custom_repeats.json"
        custom_file.write_text(json.dumps(custom_data), encoding="utf-8")

        rd = load_repeat_dictionary(path=custom_file)
        assert len(rd.repeats) == 2
        assert rd.repeats["X"] == "A" * 60
        assert rd.flanking_left == "LEFTSEQ"
        assert rd.flanking_right == "RIGHTSEQ"
        assert rd.vntr_region == "chr1:100-200"
        assert rd.source == "test"

    def test_missing_file_raises(self, tmp_path: Path):
        """load_repeat_dictionary raises FileNotFoundError for missing paths."""
        with pytest.raises(FileNotFoundError):
            load_repeat_dictionary(path=tmp_path / "does_not_exist.json")


class TestClassifyRepeatId:
    """Tests for classify_repeat_id."""

    @pytest.fixture
    def dictionary(self) -> RepeatDictionary:
        return load_repeat_dictionary()

    def test_classify_pre_repeat(self, dictionary: RepeatDictionary):
        """Pre-repeat IDs are classified as 'pre'."""
        for rid in ["1", "2", "3", "4", "4p", "5", "5C"]:
            assert classify_repeat_id(rid, dictionary) == "pre", f"Expected 'pre' for '{rid}'"

    def test_classify_after_repeat(self, dictionary: RepeatDictionary):
        """After-repeat IDs are classified as 'after'."""
        for rid in ["6", "6p", "7", "8", "9"]:
            assert classify_repeat_id(rid, dictionary) == "after", f"Expected 'after' for '{rid}'"

    def test_classify_canonical_repeat(self, dictionary: RepeatDictionary):
        """Canonical repeat ID 'X' is classified as 'canonical'."""
        assert classify_repeat_id("X", dictionary) == "canonical"

    def test_classify_variable_repeat(self, dictionary: RepeatDictionary):
        """Known variant repeat IDs (e.g. 'A', 'B') are classified as 'variable'."""
        for rid in ["A", "B", "C", "D"]:
            assert classify_repeat_id(rid, dictionary) == "variable", (
                f"Expected 'variable' for '{rid}'"
            )

    def test_classify_unknown_id_is_variable(self, dictionary: RepeatDictionary):
        """An unrecognised repeat ID falls back to 'variable'."""
        assert classify_repeat_id("UNKNOWN_ZZZ", dictionary) == "variable"

    def test_classify_uses_bundled_dictionary_by_default(self):
        """classify_repeat_id loads the bundled dictionary when none is provided."""
        assert classify_repeat_id("X") == "canonical"
        assert classify_repeat_id("1") == "pre"
        assert classify_repeat_id("6") == "after"
        assert classify_repeat_id("A") == "variable"


class TestMutationCatalog:
    """Tests for mutation catalog loading and sequence pre-computation."""

    def test_mutations_loaded(self):
        """Repeat dictionary loads mutation definitions."""
        rd = load_repeat_dictionary()
        assert hasattr(rd, "mutations")
        assert "dupC" in rd.mutations

    def test_mutation_has_required_fields(self):
        """Each mutation has allowed_repeats and changes."""
        rd = load_repeat_dictionary()
        dupc = rd.mutations["dupC"]
        assert "allowed_repeats" in dupc
        assert "changes" in dupc
        assert "X" in dupc["allowed_repeats"]

    def test_mutated_sequences_precomputed(self):
        """Pre-computed mutated sequences map sequence -> (repeat_id, mutation_name)."""
        rd = load_repeat_dictionary()
        assert hasattr(rd, "mutated_sequences")
        assert len(rd.mutated_sequences) > 0
        # dupC on X should produce a 61bp sequence
        found_dupc = any(
            name == "dupC" and repeat == "X" for seq, (repeat, name) in rd.mutated_sequences.items()
        )
        assert found_dupc

    def test_dupc_sequence_is_61bp(self):
        """dupC on X inserts C before position 60 (1-based), producing 61bp."""
        rd = load_repeat_dictionary()
        for seq, (repeat_id, mut_name) in rd.mutated_sequences.items():
            if repeat_id == "X" and mut_name == "dupC":
                assert len(seq) == 61
                # dupC inserts C before 1-based pos 60 = 0-based index 59
                x_seq = rd.repeats["X"]
                assert seq == x_seq[:59] + "C" + x_seq[59:]
                break
        else:
            pytest.fail("dupC on X not found in mutated_sequences")
