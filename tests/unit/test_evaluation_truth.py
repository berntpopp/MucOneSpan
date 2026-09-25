"""Strict MucOneUp adapter fixtures, independent of generated data."""

import gzip
import json
from dataclasses import replace

import pytest

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.truth import TruthValidationError, load_truth


def fixture(tmp_path, mutant=True):
    rd = replace(
        load_repeat_dictionary(),
        repeats={"X": "AC"},
        flanking_left="TT",
        flanking_right="GG",
        mutations={
            "dupC": {
                "allowed_repeats": ["X"],
                "changes": [{"type": "insert", "start": 2, "sequence": "C"}],
            }
        },
    )
    seq = "ACC" if mutant else "AC"
    (tmp_path / "s.simulated.fa").write_text(f">haplotype_1\nTT{seq}GG\n")
    (tmp_path / "s.vntr_structure.txt").write_text(
        "haplotype_1\tXm\n" if mutant else "haplotype_1\tX\n"
    )
    stats = {
        "haplotype_statistics": [
            {
                "repeat_count": 1,
                "vntr_length": len(seq),
                "repeat_lengths": [2],
                "mutant_repeat_count": int(mutant),
                "mutation_details": [{"position": 1, "repeat": "X"}] if mutant else [],
            }
        ],
        "mutation_info": {"mutation_name": "dupC", "mutation_targets": ["1,1"]} if mutant else {},
        "provenance": {"seed": None},
    }
    (tmp_path / "s.simulation_stats.json").write_text(json.dumps(stats))
    if mutant:
        (tmp_path / "s.mutated_unit.fa").write_text(">haplotype_1_repeat_1\nACC\n")
    return rd


def test_valid_actual_mutated_unit_and_historical_warning(tmp_path):
    rd = fixture(tmp_path)
    result = load_truth(tmp_path, rd)
    assert result.haplotypes[0].sequence == "ACC"
    assert result.haplotypes[0].structure == ("X:dupC",)
    assert "stale_repeat_lengths:haplotype_1" in result.warnings
    assert "missing_seed" in result.warnings
    assert result.hashes


@pytest.mark.parametrize(
    "file,content",
    [
        ("s.simulated.fa", ">haplotype_1\nTTACGG\n"),
        ("s.simulated.fa", ">haplotype_1\nTTACCGG\n>haplotype_1\nTTACCGG\n"),
        ("s.vntr_structure.txt", "haplotype_1\tY\n"),
        ("s.mutated_unit.fa", ">haplotype_1_repeat_1\nAAA\n"),
    ],
)
def test_reject_inconsistent_truth(tmp_path, file, content):
    rd = fixture(tmp_path)
    (tmp_path / file).write_text(content)
    with pytest.raises(TruthValidationError):
        load_truth(tmp_path, rd)


def test_missing_mutated_unit_and_ambiguous_files(tmp_path):
    rd = fixture(tmp_path)
    (tmp_path / "s.mutated_unit.fa").unlink()
    with pytest.raises(TruthValidationError):
        load_truth(tmp_path, rd)
    fixture(tmp_path, False)
    (tmp_path / "other.simulated.fa").write_text(">h\nA\n")
    with pytest.raises(TruthValidationError):
        load_truth(tmp_path, rd)


def test_empty_flanks_are_not_negative_zero_slice(tmp_path):
    rd = fixture(tmp_path, False)
    rd = replace(rd, flanking_left="", flanking_right="")
    (tmp_path / "s.simulated.fa").write_text(">haplotype_1\nAC\n")
    assert load_truth(tmp_path, rd).haplotypes[0].sequence == "AC"


@pytest.mark.parametrize(
    "change",
    ["count", "length", "events", "target", "snp", "missing_haplotype", "unknown_mutation"],
)
def test_invalid_metadata_does_not_create_truth(tmp_path, change):
    rd = fixture(tmp_path)
    path = tmp_path / "s.simulation_stats.json"
    stats = json.loads(path.read_text())
    row = stats["haplotype_statistics"][0]
    if change == "count":
        row["repeat_count"] = 2
    elif change == "length":
        row["vntr_length"] = 2
    elif change == "events":
        row["mutation_details"] = []
    elif change == "target":
        stats["mutation_info"]["mutation_targets"] = ["1,999"]
    elif change == "snp":
        row["snp_count"] = 1
    elif change == "missing_haplotype":
        stats["haplotype_statistics"].append(dict(row))
    else:
        stats["mutation_info"]["mutation_name"] = "unknown"
    path.write_text(json.dumps(stats))
    with pytest.raises(TruthValidationError):
        load_truth(tmp_path, rd)


def test_metadata_nominal_molecules_are_separate_from_stale_coverage(tmp_path):
    rd = fixture(tmp_path, False)
    (tmp_path / "reads_metadata.tsv").write_text(
        "Parameter\tValue\nCoverage\t30\nCommand\tmuconeup reads amplicon input.fa --coverage 200 --seed 9101\n"
    )
    result = load_truth(tmp_path, rd)
    assert result.provenance["nominal_molecules"] == 200
    assert result.provenance["retained_records"] is None
    assert result.provenance["read_simulation_seed"] == 9101
    assert result.provenance["read_metadata"]["Coverage"] == "30"


_READ_TRUTH_HEADER = (
    "read_id\thap\tmolecule\tkind\tstrand\tsrc_start\tsrc_end\tn_hp_edits\thp_edits\tdetail\n"
)


def test_read_truth_manifest_marked_available(tmp_path):
    rd = fixture(tmp_path)
    with gzip.open(tmp_path / "x_read_truth.tsv.gz", "wt") as fh:
        fh.write(_READ_TRUTH_HEADER)
    truth = load_truth(tmp_path, rd)
    assert truth.provenance["read_source_truth"] == "available"
    assert "x_read_truth.tsv.gz" in truth.hashes


def test_read_truth_manifest_missing_or_ambiguous_is_unavailable(tmp_path):
    rd = fixture(tmp_path)
    truth = load_truth(tmp_path, rd)
    assert truth.provenance["read_source_truth"] == "unavailable"
    assert not any(name.endswith("_read_truth.tsv.gz") for name in truth.hashes)

    for name in ("x_read_truth.tsv.gz", "y_read_truth.tsv.gz"):
        with gzip.open(tmp_path / name, "wt") as fh:
            fh.write(_READ_TRUTH_HEADER)
    truth = load_truth(tmp_path, rd)
    assert truth.provenance["read_source_truth"] == "unavailable"
    assert not any(name.endswith("_read_truth.tsv.gz") for name in truth.hashes)
