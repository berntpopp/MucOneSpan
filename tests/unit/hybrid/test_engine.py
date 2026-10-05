"""Engine orchestration and the end-to-end hybrid pipeline on synthetic FASTQ."""

from __future__ import annotations

import json
from dataclasses import replace
from importlib import metadata
from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.evaluation import load_observation
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.igv_gene import display_flanks, load_gene_model, parse_region
from muc_one_span.hybrid.phase import PhaseResult
from muc_one_span.hybrid.spans import ReadRecord, SpanRead
from muc_one_span.pipeline import execute_pipeline
from muc_one_span.report import compute_clinical_decision
from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth

A = synth.allele(["X"] * 25)
RD_X = synth.RD.repeats["X"]
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
    # Head-to-tail PCR dimers (A+A, A+B, B+B) make length peaks near L_a + L_b;
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


def test_hybrid_igv_report_shows_each_alleles_reads_and_the_sorted_dupc_site(
    tmp_path: Path,
) -> None:
    out = tmp_path / "o"
    captured: dict[str, Any] = {}
    maps: list[tuple[list[str], str, str]] = []

    def fake_map(reads, reference, output_dir, threads, *, preset, timeout):
        # Per-allele inputs are removed after merging: keep their content.
        headers = reads.read_text().splitlines()[::4]
        maps.append((headers, reference.read_text().splitlines()[0], preset))
        bam = output_dir / "mapping.bam"
        bam.write_bytes(b"")
        return bam

    def fake_context(**kwargs):
        captured.update(kwargs)
        return {
            "mode": "off",
            "has_igv": False,
            "igv_content": "",
            "sidecar_path": None,
            "table_json": "[]",
            "session_dictionary": "{}",
        }

    with (
        patch("muc_one_span.tools.check_tools") as check,
        patch("muc_one_span.tools.get_tool_versions", return_value={"minimap2": "2.28"}),
        patch("muc_one_span.report_igv.preflight_igv_report") as preflight,
        patch("muc_one_span.hybrid.igv_alignment.map_reads", side_effect=fake_map),
        patch("muc_one_span.hybrid.igv_alignment.run_tool") as tool,
        patch("muc_one_span.report_igv.build_igv_context", side_effect=fake_context),
    ):
        execute_pipeline(
            str(_sample(tmp_path)),
            str(out),
            None,
            "",
            1,
            10,
            5.0,
            True,
            "ont",
            None,
            report_igv="embedded",
            igv_session=True,
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )
    h = DEFAULT_SETTINGS.hybrid
    check.assert_called_once_with(["minimap2", "samtools", "create_report"])
    preflight.assert_called_once()
    igv = out / "igv"
    reference = igv / "igv_reference.fa"
    assert captured["fasta_path"] == reference and captured["bam_path"] == igv / "mapping.bam"
    commands = [c.args[0][:2] for c in tool.call_args_list]
    assert commands == [["samtools", "merge"], ["samtools", "index"], ["samtools", "faidx"]]
    assert not (igv / "allele_1").exists() and not (igv / "allele_2").exists()
    # Each allele's own reads (A reads r1_*, B reads r2_*) go to its own contig only.
    assert [m[2] for m in maps] == [h.igv_minimap2_preset] * 2
    for (headers, header, _), allele, prefix in zip(
        maps, ("allele_1", "allele_2"), ("@r1_", "@r2_"), strict=True
    ):
        assert headers and all(name.startswith(prefix) for name in headers)
        assert header == f">hybrid_{allele}"
    session = (igv / "igv_session.xml").read_text()
    assert 'genome="igv_reference.fa"' in session and 'path="mapping.bam"' in session
    lines = reference.read_text().splitlines()
    contigs = dict(zip(lines[::2], lines[1::2], strict=True))
    # The dupC unit is shown as its canonical parent X, so dupC reads show the +C; the
    # flanks hold the MUC1 gene model plus the margin.
    left, right = display_flanks(
        load_gene_model(),
        parse_region(synth.RD.vntr_region),
        h.igv_gene_margin_bp,
        (len(synth.RD.flanking_left), len(synth.RD.flanking_right)),
    )
    shown = contigs[">hybrid_allele_2"]
    assert shown[left : len(shown) - right] == B.replace(synth.dupc(), RD_X)
    tracks = dict(captured["annotation_tracks"])
    assert list(tracks) == ["MUC1 gene", "Repeat units", "Detected mutations"]
    (feature,) = tracks["Detected mutations"].read_text().splitlines()
    contig, start, _, name = feature.split("\t")
    assert contig == "hybrid_allele_2" and name.endswith(":X:dupC")
    site = int(start) + 1
    assert contigs[">hybrid_allele_2"][site] == "C" != contigs[">hybrid_allele_2"][site - 1]
    loci = captured["bed_path"].read_text().splitlines()
    assert loci[2].split("\t") == [contig, str(site), str(site + 1), name]
    unit = synth.RD.repeat_length_bp
    assert captured["window"] == 2 * h.igv_context_units * unit
    assert captured["flanking"] == 2 * max(len(seq) for seq in contigs.values())
    summary = json.loads((out / "summary.json").read_text())
    assert summary["tool_versions"]["minimap2"] == "2.28"
    assert summary["run_status"]["status"] == "completed"


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

    def fake(
        cons: str, members: list[SpanRead], settings: Any, rng: Any, **kwargs: Any
    ) -> PhaseResult:
        seen.append(1)
        if len(seen) > 1:
            return real(cons, members, settings, rng, **kwargs)
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
    # The junk read matches no allele, so it is counted exactly once on top of the
    # length model's own unassigned reads (the same split without it).
    (tmp_path / "plain").mkdir()
    with patch("muc_one_span.hybrid.engine.split_by_linked_sites", _split_first_peak([])):
        plain = reconstruct_alleles(fq, tmp_path / "plain", synth.RD, DEFAULT_SETTINGS)
    assert plain.block["phase_unassigned_spanning_reads"] == 0
    assert result.block["unassigned_spanning_reads"] == (
        plain.block["unassigned_spanning_reads"] + 1
    )


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


def test_igv_reads_per_allele_are_capped_by_a_seeded_subsample(tmp_path: Path) -> None:
    from muc_one_span.pipeline import _hybrid_igv_inputs

    cap = 5
    settings = replace(
        DEFAULT_SETTINGS, hybrid=replace(DEFAULT_SETTINGS.hybrid, igv_max_reads_per_allele=cap)
    )
    seq = synth.allele(["X"] * 4)
    names = {"allele_1": [f"r{i}" for i in range(40)], "allele_2": ["a", "b"]}
    shown: list[dict[str, list[str]]] = []

    def fake_align(input_path, fasta, read_names, out_dir, **_):
        shown.append(read_names)
        return out_dir / "mapping.bam"

    with patch("muc_one_span.hybrid.igv_alignment.align_assigned_reads", side_effect=fake_align):
        for _ in range(2):
            _hybrid_igv_inputs(
                tmp_path,
                "in.fastq",
                synth.RD,
                settings,
                {"allele_1": seq, "allele_2": seq},
                names,
                {k: classify_sequence(seq, synth.RD) for k in names},
            )
    first, second = shown
    assert first == second  # seeded: identical across runs
    assert len(first["allele_1"]) == cap and set(first["allele_1"]) <= set(names["allele_1"])
    assert first["allele_2"] == ["a", "b"]


def _igv_run(tmp_path: Path, **options: Any) -> tuple[Path, Any, Any]:
    out = tmp_path / "o"

    def fake_align(input_path, fasta, read_names, out_dir, **_):
        bam = out_dir / "mapping.bam"
        bam.write_bytes(b"")
        return bam

    with (
        patch("muc_one_span.tools.check_tools") as check,
        patch("muc_one_span.tools.get_tool_versions", return_value={}),
        patch("muc_one_span.report_igv.preflight_igv_report") as preflight,
        patch("muc_one_span.hybrid.igv_alignment.align_assigned_reads", side_effect=fake_align),
        patch("muc_one_span.report_igv.build_igv_context", return_value={"mode": "off"}),
    ):
        execute_pipeline(
            str(_sample(tmp_path)),
            str(out),
            None,
            "",
            1,
            10,
            5.0,
            True,
            "ont",
            None,
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
            **options,
        )
    return out, check, preflight


def test_igv_session_alone_keeps_igv_without_create_report(tmp_path: Path) -> None:
    out, check, preflight = _igv_run(tmp_path, igv_session=True)
    check.assert_called_once_with(["minimap2", "samtools"])
    preflight.assert_not_called()
    kept = {p.name for p in (out / "igv").iterdir()}
    assert {"igv_session.xml", "igv_reference.fa", "mapping.bam", "mutations.bed"} <= kept


def test_report_igv_without_session_removes_the_intermediate_folder(tmp_path: Path) -> None:
    out, check, preflight = _igv_run(tmp_path, report_igv="embedded")
    check.assert_called_once_with(["minimap2", "samtools", "create_report"])
    preflight.assert_called_once()
    assert (out / "report.html").exists() and not (out / "igv").exists()
