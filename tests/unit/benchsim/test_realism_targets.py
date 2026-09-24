"""Realism target loading and spec section 4 scoring (no edlib needed).

Target shapes follow the real `targets/prjeb92208_v1.json` keys, not the
brief's illustrative shape.
"""

from pathlib import Path

import pytest

from muc_one_span.benchsim.realism_targets import (
    PROFILE_SECTIONS,
    compare,
    js_distance,
    load_targets,
)


def test_packaged_targets_are_public_only() -> None:
    targets = load_targets()
    assert set(targets) == {"_meta", "ont_amplicon_PRJEB92208", "ont_wgs_PRJEB92208"}
    assert set(PROFILE_SECTIONS.values()) <= set(targets)
    amp = targets["ont_amplicon_PRJEB92208"]
    assert amp["hp_P_obs_given_true"]["C7|+"]["p_correct"] == pytest.approx(0.5174)


def test_js_distance_bounds() -> None:
    assert js_distance([1, 2, 3], [2, 4, 6]) == pytest.approx(0.0)
    assert js_distance([1, 0], [0, 1]) == pytest.approx(1.0)
    assert js_distance([0, 0], [1, 1]) is None


TARGETS = {
    "ont_amplicon_PRJEB92208": {
        "hp_P_obs_given_true": {"C7|+": {"p_correct": 0.52}, "C7|-": {"p_correct": 0.89}},
        "error_rates_per_ref_base": {
            "all:total": {"min": 0.016, "median": 0.021, "max": 0.033},
            "+:del": {"min": 0.009, "median": 0.011, "max": 0.017},
        },
        "category_frac": {"off_target": {"min": 0.06, "median": 0.30, "max": 0.77}},
        "span_between_alleles_frac": {"min": 0.01, "median": 0.023, "max": 0.042},
        "span_offset_pmf_15bp_bins": {"lt55u": {"p": [0.0] * 12 + [1.0] + [0.0] * 4}},
        "allelic_ratio": {"through_origin_b_per_unit": -0.056},
    }
}


def test_compare_flags_out_of_tolerance() -> None:
    metrics = {
        "c7_correct": {"+": 0.60, "-": 0.88, "both": 0.7},
        "error_rate": {"all": 0.022, "+": None},
        "deletion_rate": {"+": 0.02},
        "offtarget_frac": 0.9,
        "span_between_alleles_frac": 0.02,
        "span_offset_hist": {"lt55u": [0] * 12 + [9] + [0] * 4, "ge55u": [0] * 17},
        "log_ratio_slope": -0.05,
    }
    res = compare(metrics, TARGETS, "ont_amplicon_r10")
    assert res["c7_correct_+"]["pass"] is False and res["c7_correct_-"]["pass"] is True
    assert res["error_rate_all"]["pass"] is True
    assert res["deletion_rate_+"]["pass"] is False
    assert res["offtarget_frac"]["pass"] is False
    assert res["span_between_alleles_frac"]["pass"] is True
    assert res["span_offset_jsd_lt55u"]["pass"] is True
    assert "span_offset_jsd_ge55u" not in res and "error_rate_+" not in res
    assert res["allele_ratio_slope"]["pass"] is True


def test_compare_aggregate_range_and_median() -> None:
    inside = {"span_between_alleles_frac": {"min": 0.015, "median": 0.024, "max": 0.04}}
    shifted = {"span_between_alleles_frac": {"min": 0.015, "median": 0.035, "max": 0.04}}
    wide = {"span_between_alleles_frac": {"min": 0.0, "median": 0.024, "max": 0.04}}
    key = "span_between_alleles_frac"
    assert compare(inside, TARGETS, "ont_amplicon_PRJEB92208")[key]["pass"] is True
    assert compare(shifted, TARGETS, "ont_amplicon_PRJEB92208")[key]["pass"] is False
    assert compare(wide, TARGETS, "ont_amplicon_PRJEB92208")[key]["pass"] is False


def test_compare_unknown_profile() -> None:
    with pytest.raises(KeyError, match="hifi_amplicon"):
        compare({}, TARGETS, "hifi_amplicon")


def test_load_local_targets(tmp_path: Path) -> None:
    path = tmp_path / "local.json"
    path.write_text('{"custom": {"category_frac": {"spanning": {"min": 0.1, "max": 0.2}}}}')
    res = compare({"spanning_frac": 0.15}, load_targets(path), "custom")
    assert res["spanning_frac"]["pass"] is True
