"""Benchmark sets: named technical factor mixes, per profile (``sets`` config section).

A *set* fixes the technical factor levels a design may draw (depth, PCR bias,
error level, smear and chimera rates) for every profile; the biological factors
(length classes, compositions, events and positions, the normal fraction) stay
in `bench_config.DesignConfig` and are shared by every set. Splits keep their
meaning (seed stream and size); ``design --split S --set X`` crosses the two.

Default sets (docs/benchmark.md, "Benchmark sets"):

- ``standard``, the headline set: realistic read depths, calibrated error, no
  or calibrated PCR bias, and the typical (median) smear and chimera rates of
  the public PRJEB92208 amplicons. The decision rule applies to it only.
- ``clean``, a control: the top ``standard`` depth without PCR bias, smear or
  chimeras (comparable with historical, artefact-free benchmarks).
- ``stress``, reported separately and never a headline number: the harsh mix
  of the first pilots (depths down to 3 reads, poor error, strong PCR bias and
  smear above the real maximum).

``offpeak_share_cap`` bounds max smear + max chimera per profile (the expected
off-by->1-unit product share); ``None`` exempts a set.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

from .bench_checks import (
    ERROR_LEVEL_NAMES,
    PCR_LEVEL_NAMES,
    PROFILE_NAMES,
    check_levels,
    check_names,
    check_num,
)

SET_NAME = re.compile(r"[a-z][a-z0-9_]*")  # set names become part of design IDs


@dataclass(frozen=True)
class FactorLevels:
    """Technical factor levels of one set for one profile (stratified in `design`)."""

    depths: tuple[int, ...]
    pcr_levels: tuple[str, ...]
    error_levels: tuple[str, ...]
    smear_levels: tuple[float, ...]
    chimera_levels: tuple[float, ...]
    concatemer_levels: tuple[float, ...]
    offtarget_levels: tuple[float, ...]


# Molecule artefact levels (MucOneUp ``molecules`` rates; 0 = no such product).
ARTEFACT_LEVELS = ("smear_levels", "chimera_levels", "concatemer_levels", "offtarget_levels")


@dataclass(frozen=True)
class BenchSet:
    """One benchmark set: a description, an off-peak cap and levels per profile."""

    description: str
    offpeak_share_cap: float | None
    profiles: dict[str, FactorLevels]


def _uniform(
    depths: Mapping[str, tuple[int, ...]],
    pcr: Mapping[str, tuple[str, ...]],
    error: tuple[str, ...],
    artefacts: tuple[Mapping[str, tuple[float, ...]], ...],
) -> dict[str, FactorLevels]:
    """Levels per profile; ``artefacts`` maps are in `ARTEFACT_LEVELS` order."""
    smear, chimera, concatemer, offtarget = artefacts
    return {
        p: FactorLevels(depths[p], pcr[p], error, smear[p], chimera[p], concatemer[p], offtarget[p])
        for p in depths
    }


AMPLICON_PROFILES = ("ont_amplicon_r10", "hifi_amplicon")
GENOMIC_PROFILE = "ont_genomic_targeted"  # no PCR or molecule-artefact step (`profiles`)
NO_ARTEFACT = (0.0,)
# Typical amplicon molecule artefacts: the ont_r10_sup_amplicon_v1 MucOneUp profile
# rates, fitted to the PRJEB92208 median smear share (0.2361 of spanning reads).
TYPICAL_SMEAR = 0.24
TYPICAL_CHIMERA = 0.023
# Concatemer rate and off-target fraction of the base profiles, kept by standard and
# stress: ont_r10_sup_amplicon_v1 sets them; hifi_amplicon_v1 keeps MucOneUp's 0.
PROFILE_CONCATEMER = {"ont_amplicon_r10": (0.024,), "hifi_amplicon": NO_ARTEFACT}
PROFILE_OFFTARGET = {"ont_amplicon_r10": (0.3,), "hifi_amplicon": NO_ARTEFACT}
# Median PRJEB92208 span_off_gt1unit_frac (packaged realism targets).
TYPICAL_OFFPEAK_SHARE = 0.2695


def _amplicon(levels: tuple[float, ...]) -> dict[str, tuple[float, ...]]:
    """``levels`` for both amplicon profiles, none for the genomic profile."""
    return {**dict.fromkeys(AMPLICON_PROFILES, levels), GENOMIC_PROFILE: NO_ARTEFACT}


def _profile_values() -> tuple[dict[str, tuple[float, ...]], ...]:
    return (
        {**PROFILE_CONCATEMER, GENOMIC_PROFILE: NO_ARTEFACT},
        {**PROFILE_OFFTARGET, GENOMIC_PROFILE: NO_ARTEFACT},
    )


def _standard() -> BenchSet:
    return BenchSet(
        "Headline: realistic depths, calibrated error, no or calibrated PCR bias, "
        "typical (median) smear and chimera rates, profile concatemer and off-target rates",
        TYPICAL_OFFPEAK_SHARE,
        _uniform(
            {
                "ont_amplicon_r10": (500, 1000, 2000),
                "hifi_amplicon": (200, 500, 1000),
                GENOMIC_PROFILE: (30, 60, 100),
            },
            {
                **dict.fromkeys(AMPLICON_PROFILES, ("none", "calibrated")),
                GENOMIC_PROFILE: ("none",),
            },
            ("calibrated",),
            (_amplicon((TYPICAL_SMEAR,)), _amplicon((TYPICAL_CHIMERA,)), *_profile_values()),
        ),
    )


def _clean() -> BenchSet:
    top = {p: (max(levels.depths),) for p, levels in _standard().profiles.items()}
    zero = dict.fromkeys(top, NO_ARTEFACT)
    return BenchSet(
        "Control: top standard depth, calibrated error, no PCR bias, no molecule artefacts "
        "(smear, chimera, concatemer, off-target)",
        0.0,
        _uniform(top, dict.fromkeys(top, ("none",)), ("calibrated",), (zero,) * 4),
    )


def _stress() -> BenchSet:
    amplicon_depths = (5, 10, 20, 30, 60, 150, 500, 2000)
    return BenchSet(
        "Stress (never a headline number): low depths, poor error, strong PCR bias, "
        "smear above the real maximum",
        None,
        _uniform(
            {
                **dict.fromkeys(AMPLICON_PROFILES, amplicon_depths),
                GENOMIC_PROFILE: (3, 6, 10, 20, 40, 80),
            },
            {**dict.fromkeys(AMPLICON_PROFILES, PCR_LEVEL_NAMES), GENOMIC_PROFILE: ("none",)},
            ERROR_LEVEL_NAMES,
            (_amplicon((0.05, 0.25, 0.5)), _amplicon((0.01, 0.05)), *_profile_values()),
        ),
    )


def _defaults() -> dict[str, BenchSet]:
    return {"standard": _standard(), "clean": _clean(), "stress": _stress()}


@dataclass(frozen=True)
class SetsConfig:
    """Named benchmark sets; ``headline`` is decided on and reported first.

    ``default`` is the set `design` uses without ``--set``; ``legacy`` names
    the set of designs written before sets existed (no ``bench_set`` field).
    """

    default: str = "standard"
    headline: str = "standard"
    legacy: str = "stress"
    definitions: dict[str, BenchSet] = field(default_factory=_defaults)

    def __post_init__(self) -> None:
        if not isinstance(self.definitions, dict) or not self.definitions:
            raise ValueError("sets.definitions must name at least one set")
        for name, bench_set in self.definitions.items():
            _check_set(name, bench_set)
        for key in ("default", "headline", "legacy"):
            if getattr(self, key) not in self.definitions:
                raise ValueError(f"sets.{key} must name a set in sets.definitions")

    def order(self, present: Iterable[str] | None = None) -> tuple[str, ...]:
        """Report order: the headline, the other defined sets, then unknown names sorted.

        With ``present``, only those names are returned.
        """
        names = [self.headline, *(n for n in self.definitions if n != self.headline)]
        if present is None:
            return tuple(names)
        wanted = set(present)
        return (*(n for n in names if n in wanted), *sorted(wanted - set(names)))


def _check_set(name: str, bench_set: BenchSet) -> None:
    key = f"sets.definitions.{name}"
    if not SET_NAME.fullmatch(name):
        raise ValueError(f"{key}: set name must match {SET_NAME.pattern}")
    if not isinstance(bench_set.description, str):
        raise ValueError(f"{key}.description must be a string")
    cap = bench_set.offpeak_share_cap
    if cap is not None:
        check_num(f"{key}.offpeak_share_cap", cap, 0, 1)
    if not isinstance(bench_set.profiles, dict) or not bench_set.profiles:
        raise ValueError(f"{key}.profiles must name at least one profile")
    for profile, levels in bench_set.profiles.items():
        where = f"{key}.profiles.{profile}"
        if profile not in PROFILE_NAMES:
            raise ValueError(f"{where}: unknown profile (known: {', '.join(PROFILE_NAMES)})")
        depths = levels.depths
        if (
            not isinstance(depths, tuple)
            or not depths
            or any(type(d) is not int or d < 1 for d in depths)
        ):
            raise ValueError(f"{where}.depths must be positive integers")
        check_names(f"{where}.pcr_levels", levels.pcr_levels, PCR_LEVEL_NAMES)
        check_names(f"{where}.error_levels", levels.error_levels, ERROR_LEVEL_NAMES)
        check_levels(f"{where}.smear_levels", levels.smear_levels, 0, 1)
        check_levels(f"{where}.chimera_levels", levels.chimera_levels, 0, 1)
        check_levels(f"{where}.concatemer_levels", levels.concatemer_levels, 0, 1)
        check_levels(f"{where}.offtarget_levels", levels.offtarget_levels, 0, 1)
        if max(levels.offtarget_levels) >= 1:  # MucOneUp: off-target share in [0, 1)
            raise ValueError(f"{where}.offtarget_levels must be below 1")
        products = (
            max(levels.smear_levels) + max(levels.chimera_levels) + max(levels.concatemer_levels)
        )
        if products > 1:  # MucOneUp: shares of one molecule pool
            raise ValueError(f"{where}: max smear + chimera + concatemer = {products:g} exceeds 1")
        offpeak = max(levels.smear_levels) + max(levels.chimera_levels)
        if cap is not None and offpeak > cap:
            raise ValueError(
                f"{where}: max smear + max chimera = {offpeak:g} exceeds offpeak_share_cap {cap}"
            )
