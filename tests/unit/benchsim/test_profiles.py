"""Profile variants: artefact, error and PCR levels."""

import json
from dataclasses import replace
from pathlib import Path
from unittest import mock

import pytest

from muc_one_span.benchsim.design import build_split
from muc_one_span.benchsim.profiles import builtin_profile_dir, variant_name, write_variant

BASE = {
    "schema_version": 1,
    "name": "ont_r10_sup_amplicon_v1",
    "platform": "ont",
    "config_overrides": {"amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}}},
    "molecules": {"forward_frac": 0.5, "smear_rate": 0.24, "chimera_rate": 0.023},
    "errors": {
        "mismatch_rate": 0.007,
        "insertion_rate": 0.006,
        "deletion_rate": 0.008,
        "insertion_len_pmf": {"1": 1.0},
        "deletion_len_pmf": {"1": 1.0},
    },
}


def _design():
    d = next(x for x in build_split("dev", 30, "s", ["dupC"]) if x.profile == "ont_amplicon_r10")
    return replace(d, smear=0.5, chimera=0.05, pcr="strong", error="poor")


def test_variant_applies_all_levels(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    path, sha = write_variant(base, _design(), tmp_path / "v")
    data = json.loads(path.read_text())
    assert data["molecules"]["smear_rate"] == 0.5 and data["molecules"]["chimera_rate"] == 0.05
    assert data["errors"]["mismatch_rate"] == 0.007 * 1.5
    assert data["config_overrides"]["amplicon_params"]["pcr_bias"]["alpha"] == 2 * 9.27e-5
    assert data["name"] == variant_name(_design()) and len(sha) == 64
    assert data["provenance"]["derived_from"]["name"] == "ont_r10_sup_amplicon_v1"


def test_variant_written_once_and_stable(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    a = write_variant(base, _design(), tmp_path / "v")
    b = write_variant(base, _design(), tmp_path / "v")
    assert a == b


def test_calibrated_levels_keep_base_values(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    d = replace(_design(), smear=0.24, chimera=0.023, pcr="calibrated", error="calibrated")
    data = json.loads(write_variant(base, d, tmp_path / "v")[0].read_text())
    assert data["errors"] == BASE["errors"]
    assert data["config_overrides"] == BASE["config_overrides"]


def test_builtin_profile_dir_explicit(tmp_path: Path) -> None:
    """Explicit path returned as-is."""
    result = builtin_profile_dir(tmp_path)
    assert result == tmp_path


def test_builtin_profile_dir_from_spec(tmp_path: Path) -> None:
    """Find data/read_profiles from muc_one_up package."""
    # Create the directory structure as the code expects it
    muc_one_up_pkg = tmp_path / "muc_one_up"
    muc_one_up_pkg.mkdir()
    profiles_dir = muc_one_up_pkg / "data" / "read_profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)

    mock_spec = mock.Mock()
    mock_spec.origin = str(muc_one_up_pkg / "__init__.py")
    with mock.patch(
        "muc_one_span.benchsim.profiles.importlib.util.find_spec", return_value=mock_spec
    ):
        result = builtin_profile_dir(None)
        assert result == profiles_dir


def test_builtin_profile_dir_no_spec_raises(tmp_path: Path) -> None:
    """Raise FileNotFoundError when neither explicit nor spec found."""
    with (
        mock.patch("muc_one_span.benchsim.profiles.importlib.util.find_spec", return_value=None),
        pytest.raises(FileNotFoundError, match="--muconeup-profiles"),
    ):
        builtin_profile_dir(None)
