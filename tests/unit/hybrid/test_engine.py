"""Engine orchestration and the end-to-end hybrid pipeline on synthetic FASTQ."""

from __future__ import annotations

import json
from dataclasses import replace
from importlib import metadata
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.evaluation import load_observation
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.phase import PhaseResult
from muc_one_span.hybrid.spans import ReadRecord, SpanRead
from muc_one_span.pipeline import execute_pipeline
from muc_one_span.report import compute_clinical_decision
from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth

A = synth.allele(["X"] * 25)
B = synth.allele(["X"] * 14 + [synth.dupc()] + ["X"] * 30)


def _fastq(path: Path, reads: list[ReadRecord]) -> Path:
    path.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in reads))
    return path


def _sample(tmp_path: Path, n_b: int = 60) -> Path:
    reads = synth.reads(A, 150, err=0.02, seed=1) + synth.reads(B, n_b, err=0.02, seed=2)
    return _fastq(tmp_path / "in.fastq", reads)


def test_heterozygous_dupc_sample_is_reconstructed(tmp_path: Path) -> None:
    result = reconstruct_alleles(_sample(tmp_path), tmp_path, synth.RD, DEFAULT_SETTINGS)
    seqs = {
        p.read_text().split("\n", 1)[1].replace("\n", "") for p in result.consensus_paths.values()
    }
    assert seqs == {A, B}
    assert "hybrid" not in result.alleles and "_members" not in result.alleles
    for key in ("allele_1", "allele_2"):
        info = result.alleles[key]
        assert info["engine"] == "hybrid" and info["depth_status"] == "adequate"
        assert info["selection_status"] == "resolved" and info["split_basis"] == "length"
        # Real hybrid read-support evidence, additive alongside (not replacing) the
        # ladder's dictionary-fit classify.py confidence -- see allele_fields.py.
        assert 0.0 <= info["consensus_concordance_fraction"] <= 1.0
        assert info["classification_confidence_status"] == "not_applicable_dictionary_fit_heuristic"
        fixed = DEFAULT_SETTINGS.reference_layout.fixed_repeat_count
        assert info["length"] == info["canonical_repeats"] + fixed
    assert result.block["read_categories"]["spanning"] == 210
    assert result.block["poa_backend"] == "pyabpoa"
    saved = json.loads((tmp_path / "hybrid_reads.json").read_text())
    assert saved["read_categories"]["spanning"] == 210


def test_dimer_products_are_recorded_and_kept_out_of_the_alleles(tmp_path: Path) -> None:
    # Task 15h: head-to-tail PCR dimers (A+A, A+B, B+B) make length peaks near L_a + L_b;
    # they are recorded as 'dimer' rejected peaks, not gate-relevant, and never polished.
    reads = (
        synth.reads(A, 150, err=0.02, seed=1)
        + synth.reads(B, 60, err=0.02, seed=2)
        + synth.concatemers(A, A, 3, err=0.02, seed=5)
        + synth.concatemers(A, B, 3, err=0.02, seed=6)
        + synth.concatemers(B, B, 3, err=0.02, seed=7)
    )
    path = _fastq(tmp_path / "dimers.fastq", reads)
    result = reconstruct_alleles(path, tmp_path, synth.RD, DEFAULT_SETTINGS)
    block = result.block
    dimers = [r for r in block["rejected_peaks"] if r["reason"] == "dimer"]
    assert dimers and block["dimer_product_reads"] == sum(r["support"] for r in dimers)
    assert block["dimer_product_fraction"] > 0
    assert block["selection_status"] == "resolved"
    seqs = {
        p.read_text().split("\n", 1)[1].replace("\n", "") for p in result.consensus_paths.values()
    }
    assert seqs == {A, B}
    member_lengths = [len(seq) for group in result.members.values() for seq, _ in group]
    assert max(member_lengths) < len(A) + len(A)  # no dimer product joined an allele


def test_low_depth_allele_is_marked_not_dropped(tmp_path: Path) -> None:
    result = reconstruct_alleles(_sample(tmp_path, n_b=15), tmp_path, synth.RD, DEFAULT_SETTINGS)
    statuses = sorted(result.alleles[k]["depth_status"] for k in ("allele_1", "allele_2"))
    assert statuses == ["adequate", "low"]


def test_single_unconfirmed_site_is_not_homozygous(tmp_path: Path) -> None:
    a = synth.allele(["X"] * 10 + ["Q"] + ["X"] * 19)
    reads = synth.reads(a, 40, err=0.02, seed=3) + synth.reads(
        synth.allele(["X"] * 30), 40, err=0.02, seed=4
    )
    result = reconstruct_alleles(
        _fastq(tmp_path / "s.fastq", reads), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert result.alleles["homozygous"] is False
    assert result.alleles["allele_1"]["selection_status"] == "unresolved_single_site"
    assert result.alleles["allele_2"]["candidate_duplicate_of"] == "allele_1"


def test_no_spanning_reads_is_insufficient_evidence(tmp_path: Path) -> None:
    fq = _fastq(tmp_path / "e.fastq", [ReadRecord("j", "ACGT" * 400, "5" * 1600)])
    with pytest.raises(InsufficientEvidenceError):
        reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)


def test_hybrid_pipeline_end_to_end(tmp_path: Path) -> None:
    out = tmp_path / "out"
    with patch("muc_one_span.tools.check_tools") as check:
        execute_pipeline(
            str(_sample(tmp_path)),
            str(out),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )
    check.assert_called_once_with([])
    summary = json.loads((out / "summary.json").read_text())
    assert summary["hybrid"]["poa_backend"] == "pyabpoa"
    assert "hybrid" not in summary["alleles"]
    mutations = [m for c in summary["classifications"].values() for m in c["mutations"]]
    dupc = [m for m in mutations if m["mutation_name"] == "dupC"]
    assert len(dupc) == 1 and dupc[0]["read_support"]["status"] == "supported"
    assert dupc[0]["vcf_support"] is False
    assert dupc[0]["vcf_support_status"] == "not_applicable_read_consensus"
    assert compute_clinical_decision(summary)["state"] == "PATHOGENIC"
    assert load_observation(out).status != "invalid_artifacts"


def test_hybrid_rejects_igv_tracks(tmp_path: Path) -> None:
    import click

    with patch("muc_one_span.tools.check_tools"), pytest.raises(click.BadParameter):
        execute_pipeline(
            str(_sample(tmp_path)),
            str(tmp_path / "o"),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            report_igv="embedded",
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )


def test_read_input_streams_fastq_gzip_and_bam(tmp_path: Path) -> None:
    import gzip

    from muc_one_span.hybrid.engine import read_input

    record = "@r1 extra\nacgt\n+\nIIII\n"
    (tmp_path / "a.fastq").write_text(record)
    with gzip.open(tmp_path / "a.fastq.gz", "wt") as handle:
        handle.write(record)
    lines = iter(["@r1", "ACGT", "+", "IIII"])
    with patch("muc_one_span.tools.run_tool_iter", return_value=lines) as tool:
        bam = list(read_input(tmp_path / "a.bam"))
    tool.assert_called_once_with(["samtools", "fastq", "-F", "0x900", str(tmp_path / "a.bam")])
    for got in (
        list(read_input(tmp_path / "a.fastq")),
        list(read_input(tmp_path / "a.fastq.gz")),
        bam,
    ):
        assert got == [ReadRecord("r1", "ACGT", "IIII")]


def test_truncated_fastq_is_an_error(tmp_path: Path) -> None:
    from muc_one_span.hybrid.engine import read_input

    (tmp_path / "t.fastq").write_text("@r1\nACGT\n")
    with pytest.raises(ValueError, match="truncated"):
        list(read_input(tmp_path / "t.fastq"))


def test_user_settings_reach_every_stage(tmp_path: Path) -> None:
    """A non-default configuration is honoured, never replaced by package defaults."""
    hybrid = replace(DEFAULT_SETTINGS.hybrid, polish_rounds=1, depth_adequate_spanning=1000)
    settings = replace(DEFAULT_SETTINGS, hybrid=hybrid)
    result = reconstruct_alleles(_sample(tmp_path), tmp_path, synth.RD, settings)
    for key in ("allele_1", "allele_2"):
        info = result.alleles[key]
        assert len(info["polish"]["rounds"]) == 1
        assert info["depth_status"] == "low" and info["depth_threshold"] == 1000


def _split_first_peak(unassigned_extra: list[SpanRead]) -> Any:
    """Wrap split_by_linked_sites: force a linked split of the first peak seen."""
    from muc_one_span.hybrid import engine

    real = engine.split_by_linked_sites
    seen: list[int] = []

    def fake(cons: str, members: list[SpanRead], settings: Any, rng: Any) -> PhaseResult:
        seen.append(1)
        if len(seen) > 1:
            return real(cons, members, settings, rng)
        half = len(members) // 2
        return PhaseResult(
            [members[:half], members[half:]],
            "linked_sites",
            unassigned=list(unassigned_extra),
        )

    return fake


def test_more_than_two_groups_is_unresolved_and_reported(tmp_path: Path) -> None:
    with patch("muc_one_span.hybrid.engine.split_by_linked_sites", _split_first_peak([])):
        result = reconstruct_alleles(_sample(tmp_path), tmp_path, synth.RD, DEFAULT_SETTINGS)
    assert result.block["selection_status"] == "unresolved_max_alleles"
    assert len(result.block["dropped_groups"]) == 1
    assert {result.alleles[k]["selection_status"] for k in ("allele_1", "allele_2")} == {
        "unresolved_max_alleles"
    }


def test_phase_unassigned_reads_are_reassigned_or_counted(tmp_path: Path) -> None:
    junk = SpanRead("junk", "ACGT" * 600, 20.0, "+", 0, "motif")
    fq = _sample(tmp_path)
    with patch("muc_one_span.hybrid.engine.split_by_linked_sites", _split_first_peak([junk])):
        result = reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)
    assert result.block["phase_unassigned_spanning_reads"] == 1
    assert result.block["unassigned_spanning_reads"] >= 1


def test_reassign_places_reads_by_edit_distance() -> None:
    from muc_one_span.hybrid.engine import _Group, _reassign

    groups = [_Group([], "linked_sites", A), _Group([], "linked_sites", B)]
    spans = [SpanRead(f"a{i}", r.seq, 20.0, "+", 0, "motif") for i, r in enumerate(_reads(A))]
    junk = SpanRead("junk", "ACGT" * 600, 20.0, "+", 0, "motif")
    leftover = _reassign([*spans, junk], groups, DEFAULT_SETTINGS.hybrid)
    assert leftover == [junk]
    assert len(groups[0].members) == len(spans) and groups[1].members == []


def _reads(seq: str) -> list[ReadRecord]:
    return synth.reads(seq, 3, err=0.02, seed=5, strand_mix=False, flank_bp=0)


def test_hybrid_pipeline_streams_bam_through_samtools(tmp_path: Path) -> None:
    fq = _sample(tmp_path)
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"")
    lines = fq.read_text().splitlines()
    with (
        patch("muc_one_span.tools.check_tools") as check,
        patch("muc_one_span.tools.run_tool_iter", return_value=iter(lines)) as tool,
        patch("muc_one_span.selection_qc.annotate_selection_qc") as qc,
    ):
        execute_pipeline(
            str(bam),
            str(tmp_path / "o"),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )
    check.assert_called_once_with(["samtools"])
    assert tool.call_args.args[0][:2] == ["samtools", "fastq"]
    qc.assert_not_called()
    summary = json.loads((tmp_path / "o" / "summary.json").read_text())
    assert summary["run_status"]["status"] == "completed"
    assert summary["tool_versions"].keys() == {"edlib", "pyabpoa"}


def test_malformed_fastq_header_is_an_error(tmp_path: Path) -> None:
    from muc_one_span.hybrid.engine import read_input

    (tmp_path / "m.fastq").write_text(">r1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="malformed"):
        list(read_input(tmp_path / "m.fastq"))


def test_extra_versions_reports_missing_packages() -> None:
    from muc_one_span.hybrid.engine import extra_versions

    with patch("importlib.metadata.version", side_effect=metadata.PackageNotFoundError):
        assert extra_versions("pyspoa") == {"edlib": "unknown", "pyspoa": "unknown"}
