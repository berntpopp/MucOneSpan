"""Exact seven-C-tract duplication across supported repeat backgrounds."""

import pytest

from muc_one_span.classify import classify_repeat, classify_sequence, validate_mutations_against_vcf
from muc_one_span.config import load_repeat_dictionary


def insert_c(sequence: str) -> str:
    """Insert before base60, within the terminal seven-C homopolymer."""
    return sequence[:59] + "C" + sequence[59:]


def test_b_tract_duplication_keeps_61bp_boundary_and_parent():
    dictionary = load_repeat_dictionary()
    sequence = insert_c(dictionary.repeats["B"]) + dictionary.repeats["X"] * 2
    result = classify_sequence(sequence, dictionary)
    assert result["structure"] == "B:dupC X X"
    assert [(repeat["start"], repeat["end"]) for repeat in result["repeats"]] == [
        (0, 61),
        (61, 121),
        (121, 181),
    ]
    assert len(result["mutations_detected"]) == 1
    mutation = result["mutations_detected"][0]
    assert mutation["closest_type"] == "B"
    assert mutation["mutation_name"] == "dupC"
    assert mutation["repeat_index"] == 1
    assert mutation["frameshift"] is True


@pytest.mark.parametrize("parent", ["5", "X", "B", "D", "K", "W", "M", "N", "O", "P", "Q", "R"])
def test_exact_terminal_seven_c_backgrounds(parent):
    dictionary = load_repeat_dictionary()
    sequence = dictionary.repeats[parent]
    assert sequence[51] == "G" and sequence[52:59] == "C" * 7 and sequence[59] != "C"
    result = classify_repeat(insert_c(sequence), dictionary)
    assert result["type"] == f"{parent}:dupC"
    assert result["parent_repeat"] == parent


@pytest.mark.parametrize("parent", ["B", "X", "C", "F", "H", "L"])
def test_normal_and_ineligible_tracts_never_gain_dupc(parent):
    dictionary = load_repeat_dictionary()
    unchanged = classify_sequence(dictionary.repeats[parent] + dictionary.repeats["X"], dictionary)
    assert unchanged["mutations_detected"] == []
    if parent not in ("B", "X"):
        result = classify_sequence(
            insert_c(dictionary.repeats[parent]) + dictionary.repeats["X"], dictionary
        )
        assert all(
            mutation.get("mutation_name") != "dupC" for mutation in result["mutations_detected"]
        )


@pytest.mark.parametrize("position,base", [(10, "C"), (30, "C"), (59, "A"), (59, "G"), (59, "T")])
def test_other_insertions_do_not_become_dupc(position, base):
    dictionary = load_repeat_dictionary()
    parent = dictionary.repeats["B"]
    sequence = parent[:position] + base + parent[position:] + dictionary.repeats["X"] * 2
    result = classify_sequence(sequence, dictionary)
    assert all(mutation.get("mutation_name") != "dupC" for mutation in result["mutations_detected"])


def test_b_tract_duplication_has_exact_vcf_support_at_equivalent_anchor(tmp_path):
    dictionary = load_repeat_dictionary()
    reference = dictionary.repeats["B"] + dictionary.repeats["X"] * 2
    sequence = insert_c(dictionary.repeats["B"]) + dictionary.repeats["X"] * 2
    reference_path = tmp_path / "reference.fa"
    full_path = tmp_path / "consensus.fa"
    reference_path.write_text(">test\n" + reference + "\n")
    full_path.write_text(">test\n" + sequence + "\n")
    context = {
        "chrom": "test",
        "reference_path": str(reference_path),
        "full_consensus_path": str(full_path),
        "haplotype": 1,
        "trim_start": 0,
        "trim_end": len(sequence),
    }
    variants = [
        {"chrom": "test", "pos": 52, "ref": "G", "alt": "GC", "genotype": "1/1", "qual": 35.0}
    ]
    result = validate_mutations_against_vcf(
        classify_sequence(sequence, dictionary),
        variants,
        sequence=sequence,
        repeat_dict=dictionary,
        consensus_context=context,
    )
    mutation = result["mutations_detected"][0]
    assert mutation["mutation_name"] == "dupC"
    assert mutation["closest_type"] == "B"
    assert mutation["vcf_support"] is True
    assert mutation["vcf_support_status"] == "exact_sequence_concordance"


def test_fixed_precursor_five_uses_same_exact_tract_without_shifting_next_unit():
    """Structural placement does not alter the identical seven-C event context."""
    dictionary = load_repeat_dictionary()
    prefix = "".join(dictionary.repeats[parent] for parent in ("1", "2", "3", "4"))
    normal = prefix + dictionary.repeats["5"] + dictionary.repeats["X"]
    assert classify_sequence(normal, dictionary)["mutations_detected"] == []
    mutated = prefix + insert_c(dictionary.repeats["5"]) + dictionary.repeats["X"]
    result = classify_sequence(mutated, dictionary)
    assert result["structure"] == "1 2 3 4 5:dupC X"
    assert result["mutations_detected"][0]["repeat_index"] == 5
    assert result["mutations_detected"][0]["closest_type"] == "5"
    assert result["repeats"][5]["start"] == 301
