"""Benchmark designs: factors from spec §5, stratified sampling, splits, seeds."""

from __future__ import annotations

import hashlib
import random
from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields
from typing import Any

PROFILES = ("ont_amplicon_r10", "ont_genomic_targeted", "hifi_amplicon")
DEPTHS = {
    "ont_amplicon_r10": (5, 10, 20, 30, 60, 150, 500, 2000),
    "hifi_amplicon": (5, 10, 20, 30, 60, 150, 500, 2000),
    "ont_genomic_targeted": (3, 6, 10, 20, 40, 80),
}
DELTA_CLASSES = ("0_identical", "0_different", "1", "2", "3-5", "6-20", ">20")
_DELTA_RANGE = {
    "0_identical": (0, 0),
    "0_different": (0, 0),
    "1": (1, 1),
    "2": (2, 2),
    "3-5": (3, 5),
    "6-20": (6, 20),
    ">20": (21, 90),
}
COMPOSITIONS = (("markov", 0.75), ("real_derived", 0.20), ("rare_units", 0.05))
NORMAL_FRACTION = 0.35
LENGTH_MIN, LENGTH_MAX = 20, 130
PCR_LEVELS = ("calibrated", "strong", "none")
SMEAR_LEVELS = (0.05, 0.25, 0.5)
CHIMERA_LEVELS = (0.01, 0.05)
ERROR_LEVELS = ("calibrated", "poor")


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


def _lengths(delta_class: str, rng: random.Random) -> tuple[int, int]:
    lo, hi = _DELTA_RANGE[delta_class]
    delta = rng.randint(lo, hi)
    a = rng.randint(LENGTH_MIN, LENGTH_MAX - delta)
    pair = [a, a + delta]
    rng.shuffle(pair)
    return pair[0], pair[1]


def _target(
    lengths: tuple[int, int], allele: str, position: str, rng: random.Random
) -> tuple[int, int]:
    if allele == "equal" or lengths[0] == lengths[1]:
        hap = rng.choice((1, 2))
    else:
        short = 1 if lengths[0] < lengths[1] else 2
        hap = short if allele == "shorter" else 3 - short
    length = lengths[hap - 1]
    tenth = max(1, length // 10)
    if position == "first10":
        repeat = rng.randint(1, tenth)
    elif position == "last10":
        repeat = rng.randint(length - tenth + 1, length)
    else:
        repeat = rng.randint(tenth + 1, max(tenth + 1, length - tenth))
    return hap, repeat


def build_split(
    split: str, n_per_profile: int, salt: str, mutations: Sequence[str]
) -> list[Design]:
    """Designs for one split; every profile gets ``n_per_profile`` cases."""
    if not mutations:
        raise ValueError("at least one mutation name is required")
    designs: list[Design] = []
    for profile in PROFILES:
        rng = random.Random(derive_seed(salt, f"{split}:{profile}", "design"))
        n_normal = -(-int(n_per_profile * NORMAL_FRACTION * 100) // 100)
        events: list[str | None] = [None] * n_normal + _stratum(
            list(mutations), n_per_profile - n_normal, rng
        )
        rng.shuffle(events)
        deltas = _stratum(DELTA_CLASSES, n_per_profile, rng)
        depths = _stratum(DEPTHS[profile], n_per_profile, rng)
        alleles = _stratum(("shorter", "longer", "equal"), n_per_profile, rng)
        positions = _stratum(("first10", "middle", "last10"), n_per_profile, rng)
        pcrs = _stratum(PCR_LEVELS, n_per_profile, rng)
        smears = _stratum(SMEAR_LEVELS, n_per_profile, rng)
        chimeras = _stratum(CHIMERA_LEVELS, n_per_profile, rng)
        errors = _stratum(ERROR_LEVELS, n_per_profile, rng)
        comps = [c for c, w in COMPOSITIONS for _ in range(round(w * n_per_profile))]
        comps = (comps + ["markov"] * n_per_profile)[:n_per_profile]
        rng.shuffle(comps)
        for i in range(n_per_profile):
            design_id = f"{split}-{profile}-{i + 1:04d}"
            lengths = _lengths(deltas[i], rng)
            event = events[i]
            if event:
                allele = alleles[i]
                position = positions[i]
                targets: tuple[tuple[int, int], ...] = (_target(lengths, allele, position, rng),)
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
                    derive_seed(salt, design_id, "bio"),
                    derive_seed(salt, design_id, "reads"),
                )
            )
    return designs
