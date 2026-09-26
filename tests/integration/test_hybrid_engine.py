"""Hybrid engine on the repository's generated HiFi sample (samtools + generated data)."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main

pytestmark = pytest.mark.integration
GENERATED = Path(__file__).resolve().parents[1] / "data" / "generated"
DATA = GENERATED / "sample_close_51_58"


def test_hybrid_run_on_generated_sample(tmp_path: Path) -> None:
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not on PATH")
    bams = sorted(DATA.glob("*_amplicon_aligned.bam")) if DATA.is_dir() else []
    if not bams:
        pytest.skip(f"generated test data missing: {DATA} (run make generate-testdata)")
    args = ["run", "-i", str(bams[0]), "-o", str(tmp_path), "--engine", "hybrid", "--no-report"]
    res = CliRunner().invoke(main, args)
    assert res.exit_code == 0, res.output
    summary = json.loads((tmp_path / "summary.json").read_text())
    lengths = sorted(summary["alleles"][k]["length"] for k in ("allele_1", "allele_2"))
    assert lengths == [51, 58]
    assert summary["hybrid"]["poa_backend"] == "pyabpoa"


def test_default_engine_is_hybrid_and_calls_dupc(tmp_path: Path) -> None:
    """Without --engine, run uses the hybrid engine (v0.17.0 default) end to end."""
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not on PATH")
    sample = GENERATED / "sample_dupc_60_80"
    bams = sorted(sample.glob("*_amplicon_aligned.bam")) if sample.is_dir() else []
    if not bams:
        pytest.skip(f"generated test data missing: {sample} (run make generate-testdata)")
    res = CliRunner().invoke(main, ["run", "-i", str(bams[0]), "-o", str(tmp_path), "--no-report"])
    assert res.exit_code == 0, res.output
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["configuration"]["settings"]["run"]["engine"] == "hybrid"
    assert "hybrid" in summary and summary["deprecations"] == []
    lengths = sorted(summary["alleles"][k]["length"] for k in ("allele_1", "allele_2"))
    assert lengths == [60, 80]
    names = {
        m.get("mutation_name") for c in summary["classifications"].values() for m in c["mutations"]
    }
    assert "dupC" in names
