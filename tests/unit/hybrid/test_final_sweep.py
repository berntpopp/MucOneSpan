"""Final-review sweep: safety floor, configured anchors, fail-closed refit, output hygiene."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from muc_one_span.hybrid import engine, evidence, lengths
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.reads_io import parse_fastq
from muc_one_span.hybrid.spans import Anchors, ReadRecord
from muc_one_span.settings import DEFAULT_SETTINGS, ReferenceLayoutSettings
from muc_one_span.settings_hybrid import MIN_LINKED_EVENTS
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_length_artefacts as art
from tests.unit.hybrid import test_single_event as base
from tests.unit.hybrid import test_single_event_gate as gate

S = DEFAULT_SETTINGS.hybrid


# --- I1: a "linked" split needs two events --------------------------------------------


@pytest.mark.parametrize("value", [0, 1])
def test_min_linked_sites_below_two_is_rejected(value: int) -> None:
    assert MIN_LINKED_EVENTS == 2
    with pytest.raises(ValueError, match="min_linked_sites"):
        dataclasses.replace(S, min_linked_sites=value)


def test_min_linked_sites_one_is_rejected_from_a_config_file(tmp_path: Path) -> None:
    from muc_one_span.cli import main

    config = tmp_path / "c.json"
    config.write_text(json.dumps({"schema_version": 1, "hybrid": {"min_linked_sites": 1}}))
    reads = tmp_path / "r.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    args = ["--config", str(config), "run", "-i", str(reads), "-o", str(tmp_path / "o")]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 2 and "min_linked_sites" in result.output


def test_wild_type_run_excess_at_the_minimum_linkage_is_not_pathogenic(tmp_path: Path) -> None:
    """The reviewer's probe (30% site-specific +1 excess, light stutter) at the floor."""
    hybrid = dataclasses.replace(S, min_linked_sites=MIN_LINKED_EVENTS)
    settings = dataclasses.replace(DEFAULT_SETTINGS, hybrid=hybrid)
    records = gate._stress_records(0.30, "light", S.seed)
    summary, decision = base._run(tmp_path, records, settings)
    assert "linked_sites" not in summary["hybrid"]["split_bases"]
    assert decision["state"] == "INCONCLUSIVE", decision


# --- I4: anchors come from the configured reference layout ---------------------------


def test_anchor_repeats_follow_the_configured_layout() -> None:
    default = Anchors.from_dictionary(synth.RD, S, DEFAULT_SETTINGS.reference_layout)
    layout = DEFAULT_SETTINGS.reference_layout
    assert (default.left, default.right) == (
        synth.RD.repeats[layout.left_anchor_id],
        synth.RD.repeats[layout.right_anchor_id],
    )
    other = ReferenceLayoutSettings(pre=layout.pre[1:], after=layout.after[:-1])
    moved = Anchors.from_dictionary(synth.RD, S, other)
    assert moved.left == synth.RD.repeats[other.left_anchor_id] != default.left
    assert moved.right == synth.RD.repeats[other.right_anchor_id] != default.right


# --- L210: a dimer refit that adds a peak falls back to the first fit ----------------


def test_refit_that_adds_a_new_peak_keeps_the_first_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    """Dimers form a false second allele and push a real 14-read allele into max_alleles.

    Removing the dimers would let the real allele become a peak the first fit never
    accepted, so the dimer explanation is dropped and the first fit stands (fail closed).
    """
    major = art._spans(art.SHORT_INNER, art.MAJOR_READS, 40)
    dimers = art._dimer_spans(art.SHORT_INNER, art.SHORT_INNER, 2 * art.DIMER_READS, 41)
    real = art._spans(art.LONG_INNER, 2 * art.DIMER_READS - 2, 42)
    first = art._fit(major + dimers + real, dataclasses.replace(S, dimer_recognition=False))
    assert [r["reason"] for r in first.rejected if r["reason"] == "max_alleles"] == ["max_alleles"]
    model = art._fit(major + dimers + real)
    assert art._units(model) == art._units(first)
    assert model.dimer_products == [] and art._dimers(model) == []
    assert "max_alleles" in [r["reason"] for r in model.gate_relevant_rejections]
    # Without the fail-closed branch the refit would be taken: guard the guard.
    monkeypatch.setattr(lengths, "_near", lambda *_: True)
    assert art._fit(major + dimers + real).dimer_products != []


# --- minors -----------------------------------------------------------------------------


def test_fractions_share_one_precision() -> None:
    assert engine.FRACTION_DECIMALS is evidence.FRACTION_DECIMALS


def test_resolved_homozygote_has_resolved_multiplicity(tmp_path: Path) -> None:
    reads = synth.reads(synth.allele(["X"] * 25), 150, err=0.02, seed=5)
    fq = tmp_path / "h.fastq"
    fq.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in reads))
    result = reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)
    assert result.alleles["homozygous"] is True
    assert result.alleles["allele_multiplicity_status"] == "resolved"


def test_alleles_json_is_written_once_by_the_shared_tail(tmp_path: Path) -> None:
    from muc_one_span import pipeline

    fake = engine.HybridResult({"allele_1": {}}, {}, {"poa_backend": "pyabpoa"}, {})
    with (
        patch("muc_one_span.hybrid.engine.reconstruct_alleles", return_value=fake),
        patch("muc_one_span.pipeline_tail.finish_run") as tail,
        patch("muc_one_span.tools.check_tools"),
    ):
        pipeline._run_hybrid(tmp_path, "r.fastq", synth.RD, DEFAULT_SETTINGS, {}, False)
    tail.assert_called_once()
    assert not (tmp_path / "alleles.json").exists()


@pytest.mark.parametrize(
    ("text", "match"),
    [
        ("@r1\nACGT\nIIII\n@r2\n", "separator"),
        ("@r1\nACGT\n+\nIII\n", "length"),
    ],
)
def test_malformed_fastq_record_is_a_visible_error(text: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        list(parse_fastq(text.splitlines(keepends=True)))


def test_fastq_error_excerpt_is_bounded() -> None:
    from muc_one_span.hybrid import reads_io

    header = ">" + "N" * (10 * reads_io.ERROR_EXCERPT_CHARS)
    with pytest.raises(ValueError) as err:
        list(parse_fastq([header + "\n", "A\n", "+\n", "I\n"]))
    assert len(str(err.value)) < 3 * reads_io.ERROR_EXCERPT_CHARS
    assert list(parse_fastq(["@r x\n", "acgt\n", "+r\n", "IIII\n"])) == [
        ReadRecord("r", "ACGT", "IIII")
    ]
