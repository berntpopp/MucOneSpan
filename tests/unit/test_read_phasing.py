"""Guard optional read-backed phase and preserve independent primary records."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from muc_one_span import read_phasing

READ_LIST_HEADER = (
    "#readname\tsource_id\tsample\tphaseset\thaplotype\tcovered_variants"
    "\tfirst_variant_pos\tlast_variant_pos\n"
)


def records(*genotypes: str, phase_sets: list[str] | None = None) -> list[dict]:
    return [
        {
            "chrom": "c",
            "pos": index + 1,
            "ref": "A",
            "alt": "C",
            "qual": 60.0,
            "genotype": genotype,
            "sample": "sample",
            "phase_set": phase_sets[index] if phase_sets else None,
        }
        for index, genotype in enumerate(genotypes)
    ]


@pytest.fixture
def boundary(monkeypatch):
    monkeypatch.setattr(read_phasing, "select_vcf_sample", lambda *a, **k: "sample")
    monkeypatch.setattr(read_phasing.shutil, "which", lambda *a, **k: "/tools/whatshap")
    monkeypatch.setattr(read_phasing, "parse_vcf_genotypes", lambda *a, **k: records("0/1", "0/1"))


@pytest.mark.parametrize(
    "variants,status",
    [
        ([], "no_informative_heterozygosity"),
        (records("1/1"), "no_informative_heterozygosity"),
        (records("0/1"), "single_heterozygous_unordered"),
        (records("0|1", "1|0", phase_sets=["1", "1"]), "phased"),
        (records("0/.", "0/1"), "missing_genotype"),
        (records("1", "0/1"), "non_diploid"),
    ],
)
def test_existing_evidence_is_unchanged_without_read_tool(tmp_path, monkeypatch, variants, status):
    monkeypatch.setattr(read_phasing, "parse_vcf_genotypes", lambda *a, **k: variants)
    monkeypatch.setattr(read_phasing, "select_vcf_sample", lambda *a, **k: "sample")
    monkeypatch.setattr(
        read_phasing, "run_tool", lambda *a: pytest.fail("Unexpected external call")
    )
    vcf = tmp_path / "input.vcf.gz"
    selected, metadata = read_phasing.phase_same_length_reads(
        vcf, tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == vcf
    assert metadata["status"] == "not_needed"
    assert metadata["input_phase_status"] == status
    assert metadata["output_phase_status"] == status
    assert not (tmp_path / "phase").exists()


def test_overlapping_variants_never_run_read_phaser(tmp_path, monkeypatch):
    variants = records("0|1", "1|0", phase_sets=["1", "1"])
    variants[0]["ref"] = "AAA"
    monkeypatch.setattr(read_phasing, "parse_vcf_genotypes", lambda *a, **k: variants)
    monkeypatch.setattr(read_phasing, "select_vcf_sample", lambda *a, **k: "sample")
    selected, metadata = read_phasing.phase_same_length_reads(
        tmp_path / "vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == tmp_path / "vcf"
    assert metadata["input_phase_status"] == "conflicting_variant_records"
    assert metadata["status"] == "not_needed"


def test_missing_optional_tool_is_explicit_and_retains_input(tmp_path, monkeypatch, boundary):
    monkeypatch.setattr(read_phasing.shutil, "which", lambda *a, **k: None)
    selected, metadata = read_phasing.phase_same_length_reads(
        tmp_path / "vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == tmp_path / "vcf"
    assert metadata["status"] == "unavailable"
    assert metadata["version"] is None
    assert metadata["phase_command"] is None
    assert metadata["output_phase_status"] == "unphased"


def alignment(name="collision", flag=0, tail="NM:i:0"):
    return f"{name}\t{flag}\tc\t1\t60\t2M\t*\t0\t0\tAC\tII\t{tail}\n"


def mock_tools(
    monkeypatch, source_lines, output_records=None, fail_phase=False, read_list_text=None
):
    """Mock only external executables; persist their relevant filesystem outputs."""
    monkeypatch.setattr(read_phasing, "run_tool_iter", lambda *a: iter(source_lines))
    captured = {}

    def run(command):
        if command == ["whatshap", "--version"]:
            return "1.7\n"
        if command[:2] == ["samtools", "view"]:
            captured["sam"] = Path(command[-1]).read_text()
            Path(command[command.index("-o") + 1]).write_bytes(b"BAM")
        elif command[:2] == ["whatshap", "phase"]:
            if fail_phase:
                raise RuntimeError("Phasing execution failed")
            captured["phase_command"] = command
            Path(command[command.index("-o") + 1]).write_text("phased output")
            Path(command[command.index("--output-read-list") + 1]).write_text(
                read_list_text
                if read_list_text is not None
                else (
                    READ_LIST_HEADER
                    + "phase_record_000000000000\t0\tsample\t1\t0\t2\t1\t2\n"
                    + "phase_record_000000000001\t0\tsample\t1\t1\t2\t1\t2\n"
                )
            )
        elif command[:2] == ["bcftools", "view"]:
            Path(command[command.index("-o") + 1]).write_bytes(b"VCF")
        return ""

    monkeypatch.setattr(read_phasing, "run_tool", run)
    if output_records is not None:
        monkeypatch.setattr(
            read_phasing,
            "parse_vcf_genotypes",
            lambda path, **kw: (
                records("0/1", "0/1") if path.name == "input.vcf" else output_records
            ),
        )
    return captured


def test_collision_free_phase_retains_distinct_records_and_all_other_fields(
    tmp_path, monkeypatch, boundary
):
    lines = [
        "@HD\tVN:1.6\tSO:coordinate\n",
        alignment(),
        alignment(flag=256),
        alignment(tail="NM:i:1"),
        alignment(flag=2048),
        alignment(flag=4),
    ]
    captured = mock_tools(monkeypatch, lines, records("0|1", "1|0", phase_sets=["1", "1"]))
    original = tmp_path / "input.vcf"
    selected, metadata = read_phasing.phase_same_length_reads(
        original, tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected != original and selected.is_file()
    assert metadata["status"] == "phased"
    assert metadata["version"] == "1.7"
    assert metadata["primary_records"] == 2
    assert metadata["excluded_records"] == 3
    assert metadata["phasing_assigned_records"] == 2
    written = [
        line.split("\t") for line in captured["sam"].splitlines() if not line.startswith("@")
    ]
    assert len({fields[0] for fields in written}) == 2
    assert [fields[1:] for fields in written] == [
        lines[index].rstrip("\n").split("\t")[1:] for index in (1, 3)
    ]
    sidecar = [json.loads(line) for line in Path(metadata["source_map"]).read_text().splitlines()]
    assert [row["input_record_ordinal"] for row in sidecar] == [1, 3]
    assert [row["original_qname"] for row in sidecar] == ["collision", "collision"]
    command = captured["phase_command"]
    assert metadata["phase_command"] == command
    assert "--indels" in command and "--ignore-read-groups" in command
    assert command[command.index("--sample") + 1] == "sample"
    assert not Path(command[-1]).exists(), "Temporary BAM should be removed"


def test_explicit_tool_overrides_are_recorded_and_forwarded(tmp_path, monkeypatch, boundary):
    from muc_one_span.settings import ReadPhasingSettings

    captured = mock_tools(monkeypatch, [alignment(), alignment()])
    _, metadata = read_phasing.phase_same_length_reads(
        tmp_path / "input.vcf",
        tmp_path / "bam",
        tmp_path / "ref",
        tmp_path / "phase",
        settings=ReadPhasingSettings(internal_downsampling=27, mapping_quality=6),
    )
    command = captured["phase_command"]
    assert command[command.index("--internal-downsampling") + 1] == "27"
    assert command[command.index("--mapping-quality") + 1] == "6"
    assert metadata["tool_parameter_overrides"] == {
        "internal_downsampling": 27,
        "mapping_quality": 6,
    }


@pytest.mark.parametrize(
    "output", [records("0/1", "0/1"), records("0|1", "1|0", phase_sets=["1", "2"])]
)
def test_unresolved_tool_output_cannot_replace_original_vcf(
    tmp_path, monkeypatch, boundary, output
):
    mock_tools(monkeypatch, [alignment(), alignment()], output)
    selected, metadata = read_phasing.phase_same_length_reads(
        tmp_path / "input.vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == tmp_path / "input.vcf"
    assert metadata["status"] == "unresolved"
    assert metadata["output_phase_status"] == "unphased"
    assert Path(metadata["candidate_vcf"]).is_file()


def test_empty_primary_input_is_insufficient_evidence(tmp_path, monkeypatch, boundary):
    captured = mock_tools(monkeypatch, ["@HD\tVN:1.6\n", alignment(flag=256)])
    selected, metadata = read_phasing.phase_same_length_reads(
        tmp_path / "input.vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == tmp_path / "input.vcf"
    assert metadata["status"] == "insufficient_reads"
    assert metadata["primary_records"] == 0
    assert "phase_command" not in captured


def test_present_tool_execution_failure_propagates(tmp_path, monkeypatch, boundary):
    mock_tools(monkeypatch, [alignment(), alignment()], fail_phase=True)
    with pytest.raises(RuntimeError, match="Phasing execution failed"):
        read_phasing.phase_same_length_reads(
            tmp_path / "input.vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
        )
    assert not list((tmp_path / "phase").glob("**/*.bam"))


def test_phaser_must_not_change_variant_alleles_or_genotypes(tmp_path, monkeypatch, boundary):
    mock_tools(
        monkeypatch, [alignment(), alignment()], records("1|1", "0|1", phase_sets=["1", "1"])
    )
    with pytest.raises(ValueError, match="changed variant identity or genotype"):
        read_phasing.phase_same_length_reads(
            tmp_path / "input.vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
        )


@pytest.mark.parametrize(
    "read_list_text",
    [
        READ_LIST_HEADER,
        READ_LIST_HEADER + "phase_record_000000000000\t0\tOTHER\t1\t0\t2\t1\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t99\t0\t2\t1\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t2\t2\t1\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t0\t0\t1\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t0\t3\t1\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t0\t2\t0\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t0\t2\t2\t1\n",
        READ_LIST_HEADER + "phase_record_000000000000\t1\tsample\t1\t0\t2\t1\t2\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t0\n",
        READ_LIST_HEADER + "phase_record_000000000000\t0\tsample\t1\t0\tX\t1\t2\n",
        "#unexpected header\nphase_record_000000000000\t0\tsample\t1\t0\t2\t1\t2\n",
    ],
    ids=[
        "empty",
        "sample",
        "phase_set",
        "haplotype",
        "zero_variants",
        "too_many_variants",
        "position",
        "reversed_interval",
        "source",
        "truncated",
        "noninteger",
        "header",
    ],
)
def test_phased_claim_requires_consistent_read_assignments(
    tmp_path, monkeypatch, boundary, read_list_text
):
    mock_tools(
        monkeypatch,
        [alignment(), alignment()],
        records("0|1", "1|0", phase_sets=["1", "1"]),
        read_list_text=read_list_text,
    )
    with pytest.raises(ValueError, match="phasing"):
        read_phasing.phase_same_length_reads(
            tmp_path / "input.vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
        )


def test_empty_unphased_read_list_remains_unresolved(tmp_path, monkeypatch, boundary):
    mock_tools(monkeypatch, [alignment()], records("0/1", "0/1"), read_list_text=READ_LIST_HEADER)
    selected, metadata = read_phasing.phase_same_length_reads(
        tmp_path / "input.vcf", tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == tmp_path / "input.vcf"
    assert metadata["status"] == "unresolved"
    assert metadata["phasing_assigned_records"] == 0
    assert metadata["phasing_selected_records_by_phase_set"] == {}


@pytest.mark.parametrize("haplotype_counts", [(10, 1), (11, 0)])
def test_selected_haplotype_counts_do_not_claim_all_primary_read_support(
    tmp_path, monkeypatch, boundary, haplotype_counts
):
    """Expose internal selection and imbalance without introducing a support gate."""
    haplotypes = [0] * haplotype_counts[0] + [1] * haplotype_counts[1]
    assignments = READ_LIST_HEADER + "".join(
        f"phase_record_{index:012d}\t0\tsample\t1\t{haplotype}\t2\t1\t2\n"
        for index, haplotype in enumerate(haplotypes)
    )
    mock_tools(
        monkeypatch,
        [alignment() for _ in range(30)],
        records("0|1", "1|0", phase_sets=["1", "1"]),
        read_list_text=assignments,
    )
    original = tmp_path / "input.vcf"
    selected, metadata = read_phasing.phase_same_length_reads(
        original, tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected != original
    assert metadata["status"] == "phased"
    assert metadata["primary_records"] == 30
    assert metadata["phasing_assigned_records"] == 11
    assert metadata["phasing_selected_records_by_phase_set"] == {
        "1": {"0": haplotype_counts[0], "1": haplotype_counts[1]}
    }
    assert metadata["phasing_read_count_scope"] == "internal_selection_not_all_primary_records"


def test_selected_haplotype_counts_keep_disconnected_phase_sets_separate(
    tmp_path, monkeypatch, boundary
):
    mock_tools(
        monkeypatch,
        [alignment(), alignment(), alignment()],
        records("0|1", "1|0", phase_sets=["1", "2"]),
        read_list_text=(
            READ_LIST_HEADER
            + "phase_record_000000000000\t0\tsample\t1\t0\t1\t1\t1\n"
            + "phase_record_000000000001\t0\tsample\t2\t1\t1\t2\t2\n"
        ),
    )
    original = tmp_path / "input.vcf"
    selected, metadata = read_phasing.phase_same_length_reads(
        original, tmp_path / "bam", tmp_path / "ref", tmp_path / "phase"
    )
    assert selected == original
    assert metadata["status"] == "unresolved"
    assert metadata["phasing_assigned_records"] == 2
    assert metadata["phasing_selected_records_by_phase_set"] == {
        "1": {"0": 1, "1": 0},
        "2": {"0": 0, "1": 1},
    }
