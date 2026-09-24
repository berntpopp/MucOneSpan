"""Predefined VNTR structures for ``0_identical``, ``rare_units`` and ``real_derived`` designs.

MucOneUp is only run as an external executable. Structures that need editing
are made by a first ``muconeup simulate`` with ``--fixed-lengths`` (a valid
Markov chain), then edited and written as an ``--input-structure`` file:

- ``0_identical``: one pre-run haplotype, written twice.
- ``rare_units``: about 10% of the interior units are replaced by dictionary
  units that the Markov model uses in < 1% of positions; never the first four
  motifs (1-4), the last five units, conserved unit types, or event targets.
- ``real_derived``: diploid entry ``bio_seed % n`` of a local structure pool
  (a MucOneUp structure file, haplotype lines in consecutive pairs). Only the
  pool's SHA-256 and the entry index are recorded, never the chain. Without a
  pool the design falls back to Markov and records ``composition_effective``.
"""

from __future__ import annotations

import hashlib
import json
import random
from collections.abc import Callable, Collection
from pathlib import Path
from typing import Any

from .design import Design, event_bounds

Runner = Callable[[list[str]], str]
CONSERVED = frozenset({"1", "2", "3", "4", "4p", "5", "5C", "6", "6p", "7", "8", "9"})
RARE_FRACTION = 0.10
RARE_USAGE = 0.01
_HEAD, _TAIL = 4, 5


def read_chains(path: Path) -> list[list[str]]:
    """Haplotype chains from a MucOneUp structure file, mutation markers removed."""
    chains = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 2 or not columns[1]:
            raise ValueError(f"{path}: malformed structure line {line!r}")
        chains.append([u[:-1] if u.endswith("m") else u for u in columns[1].split("-")])
    return chains


def write_chains(path: Path, chains: list[list[str]]) -> Path:
    """Write chains in the MucOneUp ``--input-structure`` format."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"haplotype_{i}\t{'-'.join(c)}\n" for i, c in enumerate(chains, 1)))
    return path


def rare_units(config: Path, known: Collection[str] | None = None) -> list[str]:
    """Non-conserved units used at < 1% of positions by the MucOneUp Markov model.

    Usage is the long-run visit frequency of the transition chain (``END``
    restarts at ``1``), averaged over iterations so periodic chains converge.
    """
    data = json.loads(config.read_text())
    probs: dict[str, dict[str, float]] = data["probabilities"]
    units = sorted(set(data["repeats"]) | set(probs))
    state = dict.fromkeys(units, 0.0)
    state["1"] = 1.0
    total = dict.fromkeys(units, 0.0)
    steps = 2000
    for _ in range(steps):
        nxt = dict.fromkeys(units, 0.0)
        for unit, mass in state.items():
            if mass:
                for target, p in probs.get(unit, {"1": 1.0}).items():
                    nxt["1" if target == "END" else target] += mass * p
        state = nxt
        for unit, mass in state.items():
            total[unit] += mass / steps
    return [
        u
        for u in units
        if u in data["repeats"]
        and u not in CONSERVED
        and total[u] < RARE_USAGE
        and (known is None or u in known)
    ]


def inject_rare(
    chain: list[str], rare: list[str], rng: random.Random, protect: Collection[int] = ()
) -> list[str]:
    """Replace ~10% of interior units with rare units (``protect`` is 1-based)."""
    if not rare:
        raise ValueError("no rare units available in the MucOneUp config")
    positions = [
        i
        for i in range(_HEAD, len(chain) - _TAIL)
        if i + 1 not in protect and chain[i] not in CONSERVED
    ]
    n = min(len(positions), max(1, round(RARE_FRACTION * (len(chain) - _HEAD - _TAIL))))
    out = list(chain)
    for i in rng.sample(positions, n):
        choices = [u for u in rare if u != chain[i]] or rare
        out[i] = rng.choice(choices)
    return out


def scale_targets(
    targets: tuple[tuple[int, int], ...],
    design_lengths: tuple[int, int],
    actual_lengths: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    """Keep each target's relative position when a structure has other lengths.

    Scaled targets are clamped to `event_bounds` of the actual chain.
    """
    out = []
    for hap, repeat in targets:
        actual = actual_lengths[hap - 1]
        scaled = round(repeat * actual / design_lengths[hap - 1])
        lo, hi = event_bounds(actual)
        out.append((hap, max(lo, min(hi, scaled))))
    return tuple(out)


def _prerun(
    design: Design, executable: str, config: Path, work: Path, n_hap: int, run: Runner
) -> list[list[str]]:
    out = work / "prerun"
    args = [executable, "--config", str(config), "simulate", "--out-dir", str(out)]
    args += ["--out-base", "prerun", "--num-haplotypes", str(n_hap), "--output-structure"]
    args += ["--seed", str(design.bio_seed)]
    for length in design.lengths[:n_hap]:
        args += ["--fixed-lengths", str(length)]
    run(args)
    files = sorted(out.glob("*.vntr_structure.txt"))
    if len(files) != 1:
        raise ValueError(f"structure pre-run wrote {len(files)} structure files")
    chains = read_chains(files[0])
    if len(chains) != n_hap:
        raise ValueError(f"structure pre-run returned {len(chains)} haplotypes, not {n_hap}")
    return chains


def _pool_entry(pool: Path, seed: int) -> tuple[list[list[str]], int, str]:
    chains = read_chains(pool)
    if len(chains) < 2 or len(chains) % 2:
        raise ValueError(f"structure pool {pool} needs haplotype lines in pairs")
    index = seed % (len(chains) // 2)
    sha = hashlib.sha256(pool.read_bytes()).hexdigest()
    return chains[2 * index : 2 * index + 2], index, sha


def prepare_structure(
    design: Design,
    executable: str,
    config: Path,
    pool: Path | None,
    work: Path,
    run: Runner,
    known: Collection[str] | None = None,
) -> tuple[Path | None, dict[str, Any]]:
    """Return an ``--input-structure`` file (or ``None`` for plain Markov) and provenance."""
    effective = design.composition
    if effective == "real_derived" and pool is None:
        effective = "markov"
    identical = design.delta_class == "0_identical"
    info: dict[str, Any] = {"composition_effective": effective, "structure_source": "markov"}
    if effective == "markov" and not identical:
        return None, info
    if effective == "real_derived":
        assert pool is not None
        chains, index, sha = _pool_entry(pool, design.bio_seed)
        info.update(structure_source="pool", pool_index=index, pool_sha256=sha)
    else:
        chains = _prerun(design, executable, config, work, 1 if identical else 2, run)
        info["structure_source"] = "prerun"
    if identical:
        chains = [chains[0], list(chains[0])]
    if effective == "rare_units":
        rng = random.Random(design.bio_seed)
        rare = rare_units(config, known)
        protect = {repeat for _, repeat in design.targets}
        edited = [inject_rare(c, rare, rng, protect) for c in chains[: 1 if identical else 2]]
        chains = [edited[0], list(edited[0])] if identical else edited
        info["rare_units_present"] = sorted({u for c in chains for u in c} & set(rare))
    return write_chains(work / "input_structure.txt", chains), info
