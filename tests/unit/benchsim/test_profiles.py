"""Profile variants: artefact, error and PCR levels."""

import json
import time
from dataclasses import replace
from pathlib import Path
from unittest import mock

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, ProfileConfig
from muc_one_span.benchsim.design import build_split
from muc_one_span.benchsim.profiles import builtin_profile_dir, variant_name, write_variant

P = DEFAULT_BENCH_CONFIG.profiles
BASE = {
    "schema_version": 1,
    "name": "ont_r10_sup_amplicon_v1",
    "platform": "ont",
    "config_overrides": {"amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}}},
    "molecules": {
        "forward_frac": 0.5,
        "smear_rate": 0.24,
        "chimera_rate": 0.023,
        "concatemer_rate": 0.024,
        "offtarget_frac": 0.3,
        "offtarget_median_bp": 367,
    },
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
    assert data["errors"]["mismatch_rate"] == 0.007 * P.poor_error_scale
    alpha = P.strong_pcr_alpha_factor * P.r10_pcr_alpha
    assert data["config_overrides"]["amplicon_params"]["pcr_bias"]["alpha"] == alpha
    assert data["name"] == variant_name(_design()) and len(sha) == 64
    assert data["provenance"]["derived_from"]["name"] == "ont_r10_sup_amplicon_v1"


def test_variant_written_once_and_stable(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    a = write_variant(base, _design(), tmp_path / "v")
    mtime_first = a[0].stat().st_mtime
    time.sleep(0.01)  # Ensure time has passed
    b = write_variant(base, _design(), tmp_path / "v")
    mtime_second = b[0].stat().st_mtime
    assert a == b
    assert mtime_first == mtime_second  # File not rewritten


def test_calibrated_levels_keep_base_values(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    d = replace(_design(), smear=0.24, chimera=0.023, pcr="calibrated", error="calibrated")
    data = json.loads(write_variant(base, d, tmp_path / "v")[0].read_text())
    assert data["errors"] == BASE["errors"]
    threads = {"ont_amplicon_params": {"threads": P.simulator_threads}}
    assert data["config_overrides"] == BASE["config_overrides"] | threads


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


def test_hifi_error_poor_sets_accuracy(tmp_path: Path) -> None:
    """HiFi profile without 'errors' dict: error=poor sets the configured accuracy_mean."""
    hifi_base = {
        "schema_version": 1,
        "name": "hifi_amplicon_v1",
        "platform": "pacbio",
        "config_overrides": {},
        "molecules": {"forward_frac": 0.5},
    }
    base = tmp_path / "hifi.json"
    base.write_text(json.dumps(hifi_base))
    d = next(x for x in build_split("dev", 30, "s", ["dupC"]) if x.profile == "hifi_amplicon")
    d = replace(d, error="poor")
    data = json.loads(write_variant(base, d, tmp_path / "v")[0].read_text())
    assert data["config_overrides"]["pacbio_params"]["accuracy_mean"] == P.hifi_poor_accuracy_mean


def test_pcr_none_sets_no_bias(tmp_path: Path) -> None:
    """PCR level 'none' sets preset to no_bias."""
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    d = next(x for x in build_split("dev", 30, "s", ["dupC"]) if x.profile == "ont_amplicon_r10")
    d = replace(d, pcr="none", error="calibrated")
    data = json.loads(write_variant(base, d, tmp_path / "v")[0].read_text())
    assert data["config_overrides"]["amplicon_params"]["pcr_bias"] == {"preset": "no_bias"}


def test_genomic_profile_excludes_artefact_and_pcr(tmp_path: Path) -> None:
    """ont_genomic_targeted: no smear/chimera in name or molecules, no PCR override."""
    genomic_base = {
        "schema_version": 1,
        "name": "ont_r10_genomic_v1",
        "platform": "ont",
        "config_overrides": {"amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}}},
        "molecules": {"forward_frac": 0.5, "smear_rate": 0.0, "chimera_rate": 0.0},
    }
    base = tmp_path / "genomic.json"
    base.write_text(json.dumps(genomic_base))
    d = next(
        x for x in build_split("dev", 30, "s", ["dupC"]) if x.profile == "ont_genomic_targeted"
    )
    d = replace(d, smear=0.5, chimera=0.05, pcr="strong", error="calibrated")
    path, _ = write_variant(base, d, tmp_path / "v")
    name = variant_name(d)
    assert "s0.5" not in name and "c0.05" not in name
    data = json.loads(path.read_text())
    assert data["molecules"]["smear_rate"] == 0.0 and data["molecules"]["chimera_rate"] == 0.0
    assert data["config_overrides"]["amplicon_params"]["pcr_bias"]["preset"] == "madritsch2025_r10"


def test_profile_levels_come_from_the_config(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    cfg = ProfileConfig(poor_error_scale=2.0, strong_pcr_alpha_factor=3.0)
    path, _ = write_variant(base, _design(), tmp_path / "v", cfg)
    data = json.loads(path.read_text())
    assert data["errors"]["mismatch_rate"] == 0.007 * 2.0
    alpha = data["config_overrides"]["amplicon_params"]["pcr_bias"]["alpha"]
    assert alpha == 3.0 * P.r10_pcr_alpha


def test_variants_are_content_addressed_across_configs(tmp_path: Path) -> None:
    # Same design levels, different variant-shaping setting: two distinct files,
    # neither overwritten nor reused; each file name carries its content hash.
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    out = tmp_path / "v"
    other = ProfileConfig(poor_error_scale=P.poor_error_scale * 2)
    path_a, sha_a = write_variant(base, _design(), out, P)
    path_b, sha_b = write_variant(base, _design(), out, other)
    assert path_a != path_b and sha_a != sha_b
    assert sha_a in path_a.name and sha_b in path_b.name
    assert json.loads(path_a.read_text())["errors"]["mismatch_rate"] == 0.007 * P.poor_error_scale
    assert json.loads(path_b.read_text())["name"] == variant_name(_design())
    assert write_variant(base, _design(), out, P) == (path_a, sha_a)


@pytest.mark.parametrize(
    ("platform", "section"), [("ont", "ont_amplicon_params"), ("pacbio", "pacbio_params")]
)
def test_variant_bounds_simulator_threads(tmp_path: Path, platform: str, section: str) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE | {"platform": platform}))
    path, sha = write_variant(base, _design(), tmp_path / "v")
    data = json.loads(path.read_text())
    assert data["config_overrides"][section]["threads"] == P.simulator_threads
    more = replace(P, simulator_threads=P.simulator_threads + 1)
    other, other_sha = write_variant(base, _design(), tmp_path / "v", more)
    assert json.loads(other.read_text())["config_overrides"][section]["threads"] == (
        P.simulator_threads + 1
    )
    assert json.loads(path.read_text())["config_overrides"][section]["threads"] == (
        P.simulator_threads + 1
    )
    assert other_sha == sha and other == path  # threads are a run-time knob, not content


def test_simulator_threads_must_be_positive() -> None:
    with pytest.raises(ValueError, match="simulator_threads"):
        replace(P, simulator_threads=0)


def test_clean_variant_has_no_molecule_artefact() -> None:
    import tempfile

    from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG as CFG

    clean = next(
        x
        for x in build_split("dev", 3, "s", ["dupC"], CFG, "clean")
        if x.profile == "ont_amplicon_r10"
    )
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "base.json"
        base.write_text(json.dumps(BASE))
        data = json.loads(write_variant(base, clean, Path(tmp) / "v")[0].read_text())
    rates = {
        k: v
        for k, v in data["molecules"].items()
        if (k.endswith("_rate") or k.endswith("_frac")) and k != "forward_frac"
    }
    assert set(rates) >= {"smear_rate", "chimera_rate", "concatemer_rate", "offtarget_frac"}
    assert all(v == 0 for v in rates.values()), rates
    assert data["config_overrides"]["amplicon_params"]["pcr_bias"] == {"preset": "no_bias"}


def test_variant_sets_concatemer_and_offtarget_levels(tmp_path: Path) -> None:
    base = tmp_path / "base.json"
    base.write_text(json.dumps(BASE))
    design = replace(_design(), concatemer=0.01, offtarget=0.2)
    path, _ = write_variant(base, design, tmp_path / "v")
    mol = json.loads(path.read_text())["molecules"]
    assert (mol["concatemer_rate"], mol["offtarget_frac"]) == (0.01, 0.2)
    assert "k0.01" in variant_name(design) and "o0.2" in variant_name(design)
    legacy = replace(_design(), concatemer=None, offtarget=None)
    mol = json.loads(write_variant(base, legacy, tmp_path / "v")[0].read_text())["molecules"]
    assert (mol["concatemer_rate"], mol["offtarget_frac"]) == (0.024, 0.3)  # base untouched
