"""Benchmark designs: factors from spec §5, stratified sampling, splits, sets, seeds.

A design belongs to a split (seed stream, size) and a benchmark set (technical
factor levels per profile, `bench_sets`). Its ID is
``<split>-<set>-<profile>-<NNNN>``. The biological draws (events, length classes,
alleles, positions, compositions, lengths, targets) and the MucOneUp structure
seed depend on the split and profile only, so every set of a split simulates the
same haplotypes; the technical draws and the read seed depend on the set.
"""

from __future__ import annotations

import hashlib
import math
import random
from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields
from fractions import Fraction
from typing import Any

from muc_one_span.settings import DEFAULT_LAYOUT

from .bench_config import DEFAULT_BENCH_CONFIG, BenchConfig, DesignConfig

PROFILES = ("ont_amplicon_r10", "ont_genomic_targeted", "hifi_amplicon")
ALLELE_CHOICES = ("shorter", "longer", "equal")
POSITION_CHOICES = ("first10", "middle", "last10")

# Conserved repeat positions from the bundled reference layout: the head (units
# 1-5; unit 1 holds part of the forward amplicon primer site) and the tail
# (units 6-9). Events and rare units never go there (`event_bounds`).
CONSERVED_HEAD = len(DEFAULT_LAYOUT.pre)
CONSERVED_TAIL = len(DEFAULT_LAYOUT.after)


@dataclass(frozen=True)
class Design:
    design_id: str
    split: str
    profile: str
    lengths: tuple[int, int]
    delta_class: str
    composition: str
    event: str | None
    event_allele: str | None
    event_position: str | None
    targets: tuple[tuple[int, int], ...]
    depth: int
    pcr: str
    smear: float
    chimera: float
    error: str
    bio_seed: int
    read_seed: int
    target_clamped: bool = False  # a drawn event target was moved inside `event_bounds`
    bench_set: str | None = None  # None: written before benchmark sets existed

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Design:
        known = {f.name for f in fields(cls)}
        row = {k: v for k, v in data.items() if k in known}
        row["lengths"] = tuple(row["lengths"])
        row["targets"] = tuple(tuple(t) for t in row["targets"])
        return cls(**row)


def derive_seed(salt: str, design_id: str, stream: str) -> int:
    """Stable 63-bit seed per (salt, design, stream)."""
    digest = hashlib.sha256(f"{salt}|{design_id}|{stream}".encode()).digest()
    return int.from_bytes(digest[:8], "big") >> 1


def _stratum(values: Sequence[Any], n: int, rng: random.Random) -> list[Any]:
    """Latin-hypercube style: each level ~n/len(values) times, shuffled."""
    out = [values[i % len(values)] for i in range(n)]
    rng.shuffle(out)
    return out


def _lengths(delta_class: str, rng: random.Random, cfg: DesignConfig) -> tuple[int, int]:
    lo, hi = cfg.delta_ranges[delta_class]
    delta = rng.randint(lo, hi)
    a = rng.randint(cfg.length_min, cfg.length_max - delta)
    pair = [a, a + delta]
    rng.shuffle(pair)
    return pair[0], pair[1]


def event_bounds(length: int) -> tuple[int, int]:
    """1-based ``(first, last)`` repeat an event (or rare unit) may occupy.

    Raises:
        ValueError: If the chain has no unit outside the conserved head and tail.
    """
    lo, hi = CONSERVED_HEAD + 1, length - CONSERVED_TAIL
    if hi < lo:
        raise ValueError(
            f"allele of {length} units is too short for an event "
            f"(needs > {CONSERVED_HEAD + CONSERVED_TAIL} units)"
        )
    return lo, hi


def _target(
    lengths: tuple[int, int], allele: str, position: str, rng: random.Random, fraction: float
) -> tuple[tuple[int, int], bool]:
    """Event target ``(hap, repeat)`` and whether the drawn repeat was clamped."""
    if allele == "equal" or lengths[0] == lengths[1]:
        hap = rng.choice((1, 2))
    else:
        short = 1 if lengths[0] < lengths[1] else 2
        hap = short if allele == "shorter" else 3 - short
    length = lengths[hap - 1]
    edge = max(1, int(Fraction(str(fraction)) * length))
    if position == "first10":
        repeat = rng.randint(1, edge)
    elif position == "last10":
        repeat = rng.randint(length - edge + 1, length)
    else:
        repeat = rng.randint(edge + 1, max(edge + 1, length - edge))
    # Clamp (not redraw) so every design consumes the same random draws.
    lo, hi = event_bounds(length)
    clamped = min(max(repeat, lo), hi)
    return (hap, clamped), clamped != repeat


def build_split(
    split: str,
    n_per_profile: int,
    salt: str,
    mutations: Sequence[str],
    config: BenchConfig = DEFAULT_BENCH_CONFIG,
    bench_set: str | None = None,
) -> list[Design]:
    """Designs for one split and benchmark set; every profile gets ``n_per_profile`` cases.

    Biological factor levels, the normal fraction and length ranges come from
    ``config.design``; technical levels from the set (default
    ``config.sets.default``) in ``config.sets.definitions``.

    Raises:
        ValueError: If there are no mutations, the split or set is unknown, or
            the set has no levels for a profile.
    """
    cfg = config.design
    name = config.sets.default if bench_set is None else bench_set
    if not mutations:
        raise ValueError("at least one mutation name is required")
    if split not in cfg.split_sizes:
        raise ValueError(f"unknown split {split!r}")
    if name not in config.sets.definitions:
        raise ValueError(f"unknown set {name!r}")
    levels = config.sets.definitions[name].profiles
    missing = [p for p in PROFILES if p not in levels]
    if missing:
        raise ValueError(f"set {name!r} has no levels for {', '.join(missing)}")
    designs: list[Design] = []
    for profile in PROFILES:
        rng = random.Random(derive_seed(salt, f"{split}:{profile}", "design"))
        tech = random.Random(derive_seed(salt, f"{split}:{name}:{profile}", "technical"))
        n_normal = math.ceil(Fraction(str(cfg.normal_fraction)) * n_per_profile)
        events: list[str | None] = [None] * n_normal + _stratum(
            list(mutations), n_per_profile - n_normal, rng
        )
        rng.shuffle(events)
        deltas = _stratum(tuple(cfg.delta_ranges), n_per_profile, rng)
        alleles = _stratum(ALLELE_CHOICES, n_per_profile, rng)
        positions = _stratum(POSITION_CHOICES, n_per_profile, rng)
        default_comp = next(iter(cfg.compositions))
        comps = [c for c, w in cfg.compositions.items() for _ in range(round(w * n_per_profile))]
        comps = (comps + [default_comp] * n_per_profile)[:n_per_profile]
        rng.shuffle(comps)
        factor = levels[profile]
        depths = _stratum(factor.depths, n_per_profile, tech)
        pcrs = _stratum(factor.pcr_levels, n_per_profile, tech)
        smears = _stratum(factor.smear_levels, n_per_profile, tech)
        chimeras = _stratum(factor.chimera_levels, n_per_profile, tech)
        errors = _stratum(factor.error_levels, n_per_profile, tech)
        for i in range(n_per_profile):
            bio_id = f"{split}-{profile}-{i + 1:04d}"
            design_id = f"{split}-{name}-{profile}-{i + 1:04d}"
            lengths = _lengths(deltas[i], rng, cfg)
            event = events[i]
            clamped = False
            if event:
                allele = alleles[i]
                position = positions[i]
                target, clamped = _target(lengths, allele, position, rng, cfg.position_fraction)
                targets: tuple[tuple[int, int], ...] = (target,)
            else:
                allele = None
                position = None
                targets = ()
            designs.append(
                Design(
                    design_id,
                    split,
                    profile,
                    lengths,
                    deltas[i],
                    comps[i],
                    event,
                    allele,
                    position,
                    targets,
                    depths[i],
                    pcrs[i],
                    smears[i],
                    chimeras[i],
                    errors[i],
                    derive_seed(salt, bio_id, "bio"),
                    derive_seed(salt, design_id, "reads"),
                    clamped,
                    name,
                )
            )
    return designs
