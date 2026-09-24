"""Realism targets (public PRJEB92208 aggregates) and spec section 4 scoring.

The packaged ``targets/prjeb92208_v1.json`` holds only ``_meta`` and the
``ont_amplicon_PRJEB92208`` / ``ont_wgs_PRJEB92208`` sections exported from the
local ``realprofile/targets.json``. In-house genomic targets stay local: load
them by path (`load_targets(path)`) and pass their section name explicitly as
`profile` (e.g. ``ont_genomic_inhouse``); `PROFILE_SECTIONS` maps benchmark
profiles to the public sections only.

Spec section 4 tolerance -> target key (per section):

==========================  ================================================  ==========================
metric (`realism`)          target key                                        pass rule
==========================  ================================================  ==========================
``<rate>[<strand>]``        ``error_rates_per_ref_base["<strand>:<type>"]``   median +/- 20 % relative
``c7_correct[+|-]``         ``hp_P_obs_given_true["C7|<strand>"].p_correct``  +/- 0.03 absolute
``span_off_gt1unit_frac``   ``span_off_gt1unit_frac`` (smear products)        real [min, max]
``span_between_alleles_``   ``span_between_alleles_frac``                     real [min, max]
``span_below_short_frac``   ``span_below_short_frac``                         real [min, max]
``offtarget_frac``          ``category_frac["off_target"]``                   real [min, max]
``spanning_frac``           ``category_frac["spanning"]``                     real [min, max]
``span_offset_hist``        ``span_offset_pmf_15bp_bins[lt55u|ge55u].p``      JS distance <= 0.1
``log_ratio_slope``         ``allelic_ratio.through_origin_b_per_unit``       +/- 0.01 per unit
==========================  ================================================  ==========================

``<rate>`` is mismatch/insertion/deletion/error_rate for type
mismatch/ins/del/total and ``<strand>`` is ``+``, ``-`` or ``all``; for an
`aggregate` it is the median of per-case rates, compared with the median of
per-library rates. A range
metric given as an `aggregate` spread must also keep its median within
+/- 25 % of the real median. Checks whose simulated value or target is
missing are omitted from the result.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from importlib import resources
from pathlib import Path
from typing import Any

PROFILE_SECTIONS = {
    "ont_amplicon_r10": "ont_amplicon_PRJEB92208",
    "ont_genomic_targeted": "ont_wgs_PRJEB92208",
}
PACKAGED_TARGETS = "prjeb92208_v1.json"
ERROR_REL_TOL = 0.20
C7_ABS_TOL = 0.03
RANGE_MEDIAN_REL_TOL = 0.25
JSD_MAX = 0.1
SLOPE_ABS_TOL = 0.01
ERROR_TYPES = {
    "mismatch_rate": "mismatch",
    "insertion_rate": "ins",
    "deletion_rate": "del",
    "error_rate": "total",
}
RANGE_TARGETS: dict[str, tuple[str, ...]] = {
    "span_off_gt1unit_frac": ("span_off_gt1unit_frac",),
    "span_between_alleles_frac": ("span_between_alleles_frac",),
    "span_below_short_frac": ("span_below_short_frac",),
    "offtarget_frac": ("category_frac", "off_target"),
    "spanning_frac": ("category_frac", "spanning"),
}


def load_targets(path: Path | None = None) -> dict[str, Any]:
    """Load realism targets: the packaged public file, or a local file by path."""
    if path is not None:
        data: dict[str, Any] = json.loads(path.read_text())
        return data
    ref = resources.files("muc_one_span.benchsim") / "targets" / PACKAGED_TARGETS
    packaged: dict[str, Any] = json.loads(ref.read_text())
    return packaged


def js_distance(p: Sequence[float], q: Sequence[float]) -> float | None:
    """Jensen-Shannon distance (base 2, in [0, 1]) of two unnormalized pmfs."""
    sp, sq = sum(p), sum(q)
    if not sp or not sq:
        return None
    pn = [x / sp for x in p]
    qn = [x / sq for x in q]

    def kl(a: list[float], m: list[float]) -> float:
        return sum(x * math.log2(x / y) for x, y in zip(a, m, strict=True) if x > 0)

    mid = [(x + y) / 2 for x, y in zip(pn, qn, strict=True)]
    return math.sqrt(max(0.0, (kl(pn, mid) + kl(qn, mid)) / 2))


def _section(targets: dict[str, Any], profile: str) -> dict[str, Any]:
    name = profile if profile in targets else PROFILE_SECTIONS.get(profile, "")
    if name not in targets:
        raise KeyError(f"no realism targets for profile {profile!r}")
    section: dict[str, Any] = targets[name]
    return section


def _dig(data: dict[str, Any], path: Sequence[str]) -> Any:
    for key in path:
        if not isinstance(data, dict) or key not in data:
            return None
        data = data[key]
    return data


def _check(sim: Any, target: Any, tolerance: str, ok: bool) -> dict[str, Any]:
    return {"sim": sim, "target": target, "tolerance": tolerance, "pass": bool(ok)}


def _range(sim: Any, real: dict[str, float]) -> dict[str, Any]:
    lo, hi = real["min"], real["max"]
    if not isinstance(sim, dict):
        return _check(sim, [lo, hi], "within real [min, max]", lo <= sim <= hi)
    median = real.get("median")
    in_range = lo <= sim["min"] and sim["max"] <= hi
    near = median is None or abs(sim["median"] - median) <= RANGE_MEDIAN_REL_TOL * median
    target = {"min": lo, "median": median, "max": hi}
    return _check(sim, target, "within real [min, max]; median +/- 25 %", in_range and near)


def compare(metrics: dict[str, Any], targets: dict[str, Any], profile: str) -> dict[str, Any]:
    """Score realism metrics against one target section (spec section 4).

    Args:
        metrics: `realism.read_metrics` (one case) or `realism.aggregate`.
        targets: Loaded targets (`load_targets`).
        profile: A benchmark profile (mapped by `PROFILE_SECTIONS`) or an
            explicit section name present in `targets`.

    Returns:
        Check name -> `{sim, target, tolerance, pass}`.

    Raises:
        KeyError: If no section exists for `profile`.
    """
    sec = _section(targets, profile)
    res: dict[str, Any] = {}
    rates = sec.get("error_rates_per_ref_base", {})
    for key, etype in ERROR_TYPES.items():
        for strand, sim in (metrics.get(key) or {}).items():
            real = _dig(rates, (f"{strand}:{etype}", "median"))
            if sim is not None and real:
                ok = abs(sim - real) <= ERROR_REL_TOL * real
                res[f"{key}_{strand}"] = _check(sim, real, "+/- 20 % relative", ok)
    for strand in ("+", "-"):
        sim = (metrics.get("c7_correct") or {}).get(strand)
        real = _dig(sec, ("hp_P_obs_given_true", f"C7|{strand}", "p_correct"))
        if sim is not None and real is not None:
            ok = abs(sim - real) <= C7_ABS_TOL
            res[f"c7_correct_{strand}"] = _check(sim, real, "+/- 0.03 absolute", ok)
    for key, path in RANGE_TARGETS.items():
        real = _dig(sec, path)
        if metrics.get(key) is not None and isinstance(real, dict):
            res[key] = _range(metrics[key], real)
    for size, hist in (metrics.get("span_offset_hist") or {}).items():
        real = _dig(sec, ("span_offset_pmf_15bp_bins", size, "p"))
        dist = js_distance(hist, real) if real else None
        if dist is not None:
            res[f"span_offset_jsd_{size}"] = _check(
                dist, 0.0, "JS distance <= 0.1", dist <= JSD_MAX
            )
    sim = metrics.get("log_ratio_slope")
    real = _dig(sec, ("allelic_ratio", "through_origin_b_per_unit"))
    if sim is not None and real is not None:
        ok = abs(sim - real) <= SLOPE_ABS_TOL
        res["allele_ratio_slope"] = _check(sim, real, "+/- 0.01 per unit", ok)
    return res
