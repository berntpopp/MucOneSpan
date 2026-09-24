"""Structure files for 0_identical, rare_units and real_derived designs."""

import json
import random
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, StructureConfig
from muc_one_span.benchsim.design import build_split, event_bounds
from muc_one_span.benchsim.structures import (
    conserved_units,
    inject_rare,
    prepare_structure,
    rare_units,
    read_chains,
    scale_targets,
)

DESIGN = next(d for d in build_split("dev", 30, "s", ["dupC"]) if d.event)
CHAIN = ["1", "2", "3", "4", "5", "C", *["X"] * 20, "6p", "7", "8", "9"]
CONFIG = {
    "repeats": dict.fromkeys(
        ["1", "2", "3", "4", "5", "C", "X", "A", "L", "K", "6p", "7", "8", "9"], "N"
    ),
    "probabilities": {
        "1": {"2": 1.0},
        "2": {"3": 1.0},
        "3": {"4": 1.0},
        "4": {"5": 1.0},
        "5": {"C": 1.0},
        "C": {"X": 1.0},
        "X": {"X": 0.8, "A": 0.1995, "L": 0.0005, "6p": 0.0},
        "A": {"X": 0.9, "6p": 0.1},
        "L": {"X": 1.0},
        "6p": {"7": 1.0},
        "7": {"8": 1.0},
        "8": {"9": 1.0},
        "9": {"END": 1.0},
    },
}


def _config(tmp_path: Path) -> Path:
    path = tmp_path / "config.json"
    path.write_text(json.dumps(CONFIG))
    return path


class Prerun:
    """Fake one-shot MucOneUp simulate that writes a structure file."""

    def __init__(self) -> None:
        self.calls: list[list[str]] = []

    def __call__(self, args: list[str], **_: Any) -> str:
        self.calls.append(args)
        out = Path(args[args.index("--out-dir") + 1])
        lengths = [int(args[i + 1]) for i, a in enumerate(args) if a == "--fixed-lengths"]
        chains = [
            ["1", "2", "3", "4", "5", "C", *["X"] * (n - 10), "6p", "7", "8", "9"] for n in lengths
        ]
        out.mkdir(parents=True, exist_ok=True)
        (out / "prerun.001.vntr_structure.txt").write_text(
            "".join(f"haplotype_{i}\t{'-'.join(c)}\n" for i, c in enumerate(chains, 1))
        )
        return ""


def test_rare_units_exclude_conserved_and_common(tmp_path: Path) -> None:
    rare = rare_units(_config(tmp_path), {"X", "A", "L", "K", "1", "9"})
    assert "L" in rare and "K" in rare  # K never reached by the model
    assert "X" not in rare and "A" not in rare
    assert not set(rare) & conserved_units()


def test_conserved_units_come_from_the_repeat_dictionary() -> None:
    from muc_one_span.config import load_repeat_dictionary

    rd = load_repeat_dictionary()
    assert conserved_units() == frozenset(rd.pre_repeat_ids) | frozenset(rd.after_repeat_ids)


def test_rare_usage_threshold_is_configured(tmp_path: Path) -> None:
    strict = StructureConfig(rare_usage_max=1e-6)
    assert "L" not in rare_units(_config(tmp_path), {"X", "A", "L", "K"}, strict)


def test_inject_rare_never_touches_motifs_or_targets() -> None:
    lo, hi = event_bounds(len(CHAIN))  # the same conserved head/tail as event targets
    chain = inject_rare(CHAIN, ["L"], random.Random(1), protect={10})
    assert chain[: lo - 1] == CHAIN[: lo - 1] and chain[hi:] == CHAIN[hi:] and chain[9] == "X"
    changed = [i for i, (a, b) in enumerate(zip(CHAIN, chain, strict=True)) if a != b]
    fraction = DEFAULT_BENCH_CONFIG.structures.rare_fraction
    assert len(changed) == round(fraction * (hi - lo + 1))
    assert all(chain[i] == "L" and lo <= i + 1 <= hi for i in changed)


def test_read_chains_strips_markers_and_comments(tmp_path: Path) -> None:
    path = tmp_path / "s.txt"
    path.write_text("# Mutation Applied: dupC\nhaplotype_1\t1-Xm-9\nhaplotype_2\t1-X-9\n")
    assert read_chains(path) == [["1", "X", "9"], ["1", "X", "9"]]


def test_scale_targets_keeps_relative_position() -> None:
    assert scale_targets(((2, 50),), (100, 100), (40, 20)) == ((2, 10),)
    assert scale_targets(((1, 15),), (100, 100), (40, 20)) == ((1, 6),)


def test_scale_targets_stay_inside_event_bounds() -> None:
    assert scale_targets(((1, 6),), (100, 100), (40, 20)) == ((1, event_bounds(40)[0]),)
    assert scale_targets(((2, 95),), (100, 100), (40, 30)) == ((2, event_bounds(30)[1]),)
    with pytest.raises(ValueError, match="too short"):
        scale_targets(((1, 5),), (100, 100), (3, 3))


def test_markov_needs_no_structure(tmp_path: Path) -> None:
    design = replace(DESIGN, composition="markov", delta_class="1")
    prerun = Prerun()
    path, info = prepare_structure(design, "mu", _config(tmp_path), None, tmp_path, prerun)
    assert path is None and info["composition_effective"] == "markov" and not prerun.calls


def test_identical_uses_one_chain_for_both(tmp_path: Path) -> None:
    design = replace(DESIGN, composition="markov", delta_class="0_identical", lengths=(30, 30))
    prerun = Prerun()
    path, info = prepare_structure(design, "mu", _config(tmp_path), None, tmp_path, prerun)
    assert path is not None
    chains = read_chains(path)
    assert len(chains) == 2 and chains[0] == chains[1] and len(chains[0]) == 30
    assert prerun.calls[0][prerun.calls[0].index("--num-haplotypes") + 1] == "1"
    assert info["structure_source"] == "prerun"


def test_rare_units_structure(tmp_path: Path) -> None:
    design = replace(DESIGN, composition="rare_units", delta_class="1", lengths=(40, 41))
    path, info = prepare_structure(design, "mu", _config(tmp_path), None, tmp_path, Prerun())
    assert path is not None
    chains = read_chains(path)
    assert [len(c) for c in chains] == [40, 41]
    assert info["rare_units_present"] and any("L" in c or "K" in c for c in chains)


def test_real_derived_pool_entry_and_fallback(tmp_path: Path) -> None:
    pool = tmp_path / "pool.txt"
    pool.write_text("s1_h1\t1-X-X-9\ns1_h2\t1-X-9\ns2_h1\t1-A-9\ns2_h2\t1-A-A-9\n")
    design = replace(DESIGN, composition="real_derived", delta_class="1", bio_seed=3)
    path, info = prepare_structure(design, "mu", _config(tmp_path), pool, tmp_path, Prerun())
    assert path is not None and read_chains(path) == [["1", "A", "9"], ["1", "A", "A", "9"]]
    assert info["pool_index"] == 1 and len(info["pool_sha256"]) == 64
    assert "sequence" not in json.dumps(info)
    ident = replace(design, delta_class="0_identical")
    path, _ = prepare_structure(ident, "mu", _config(tmp_path), pool, tmp_path, Prerun())
    assert path is not None and read_chains(path) == [["1", "A", "9"], ["1", "A", "9"]]
    _, info = prepare_structure(design, "mu", _config(tmp_path), None, tmp_path, Prerun())
    assert info["composition_effective"] == "markov"


def test_bad_pool_is_rejected(tmp_path: Path) -> None:
    pool = tmp_path / "pool.txt"
    pool.write_text("only\t1-X-9\n")
    design = replace(DESIGN, composition="real_derived")
    with pytest.raises(ValueError, match="pool"):
        prepare_structure(design, "mu", _config(tmp_path), pool, tmp_path, Prerun())
