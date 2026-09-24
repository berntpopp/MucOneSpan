"""Real MucOneUp case generation: one dev case per profile at depth 10."""

from __future__ import annotations

import os
import shutil
from dataclasses import replace
from pathlib import Path

import pytest

from muc_one_span.benchsim.design import PROFILES, build_split
from muc_one_span.benchsim.generate import FASTQ, GenerateContext, generate_case
from muc_one_span.benchsim.muconeup import require_muconeup
from muc_one_span.benchsim.read_truth import load_read_truth
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.truth import load_truth

pytestmark = pytest.mark.integration


def _environment() -> tuple[str, Path, Path]:
    missing = [tool for tool in ("muconeup", "pbsim", "ccs") if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"missing tools on PATH: {missing}")
    config = os.environ.get("MUCONEUP_CONFIG")
    if not config:
        pytest.skip("MUCONEUP_CONFIG is not set")
    try:
        version = require_muconeup("muconeup")
    except RuntimeError as exc:
        pytest.skip(str(exc))
    profiles = Path(
        os.environ.get("MUCONEUP_PROFILES")
        or Path(config).parent / "muc_one_up" / "data" / "read_profiles"
    )
    if not profiles.is_dir():
        pytest.skip("MucOneUp read profiles not found; set MUCONEUP_PROFILES")
    return version, Path(config).resolve(), profiles


def test_generate_one_case_per_profile(tmp_path: Path) -> None:
    version, config, profiles = _environment()
    ctx = GenerateContext("muconeup", config, tmp_path, profiles, None, None, version)
    designs = build_split("dev", 12, "integration", ["dupC"])
    for profile in PROFILES:
        design = next(
            d
            for d in designs
            if d.profile == profile
            and d.event
            and d.composition == "markov"
            and d.delta_class != "0_identical"
        )
        design = replace(design, depth=10, lengths=(25, 30), targets=((1, 12),))
        case = generate_case(design, ctx)
        assert case["status"] == "ok", case.get("error")
        case_dir = tmp_path / "dev" / design.design_id
        fastq = case_dir / "reads" / FASTQ[profile].format(design.design_id)
        records = sum(1 for _ in fastq.open()) // 4
        assert records > 0
        manifest = next((case_dir / "truth").glob("*_read_truth.tsv.gz"))
        assert len(load_read_truth(manifest)) == records
        truth = load_truth(case_dir / "truth", load_repeat_dictionary())
        assert truth.provenance["read_source_truth"] == "available"
        assert case["actual_targets"] == [[1, 12]]
