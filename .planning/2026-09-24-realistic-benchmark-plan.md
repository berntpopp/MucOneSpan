# Realistic Simulated Benchmark (MucSim-Bench) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate, seal, run and score a large stratified benchmark of realistic
simulated MUC1 VNTR cases (ONT R10 amplicon, ONT genomic, PacBio HiFi amplicon)
with per-read truth, and compare engines with pre-registered statistics.

**Architecture:** MucOneUp ≥ 0.45.0 generates everything that is sequenced: haplotypes
(`simulate`), then reads with a read profile and a per-read truth manifest
(`reads amplicon --read-profile …`, `reads ont --simulator pbsim3-fragments
--read-profile …`). MucOneSpan adds a tested library `src/muc_one_span/benchsim/`
that designs cases, derives profile variants, plans depth, drives MucOneUp through
`tools.run_tool`, validates the truth, runs engines, and reports statistics.
Scoring extends `evaluation/`; there is no second scorer.

**Tech Stack:** Python 3.10+, the MucOneUp 0.45.0 CLI (external executable), pbsim3,
ccs, edlib (already a dependency via the hybrid extra, or lazy-imported), numpy,
pytest, and the existing `evaluation/` and `durable_ledger.py`.

**Spec:** `.planning/2026-09-23-realistic-benchmark-spec.md` (read it first).
**Decision recorded 2026-09-24:** the spec's §3 "read layer is our own code" and its
open genomic-profile decision are superseded. MucOneUp 0.45.0 provides the read layer:

- `--read-profile ont_r10_sup_amplicon_v1` (PRJEB92208-calibrated amplicon)
- `--read-profile ont_r10_genomic_v1` with `reads ont --simulator pbsim3-fragments`
- `hifi_amplicon_v1` (uncalibrated)
- the `{base}_read_truth.tsv.gz` manifest

MucOneSpan therefore has no pbsim3 or stutter code of its own.

## Global Constraints

- Every authored file has **fewer than 650 physical lines** (max 649); split by responsibility.
- External commands only through `muc_one_span.tools.run_tool` with argument lists; no shell strings, no machine-specific paths.
- MucOneUp is an external executable, never imported. Require `muconeup --version` ≥ **0.45.0**.
- Unit tests are deterministic, synthetic, and mock `run_tool`. Real-tool tests are `@pytest.mark.integration` and skip with a reason when a tool is missing.
- Datasets, reads, truth and results live outside Git in `../MucOneSpan-bench-data/`. Commit only code, design definitions, and **public PRJEB92208 aggregate** targets.
- Patient reads, in-house structures, in-house aggregates and LB-level outputs stay local. They are passed by path at run time, and only their SHA-256 is recorded.
- Splits: `dev` 900 (300/profile), `val` 900, `test` 2,400 (800/profile), `stress` ~300. ≥35% normals in every split and profile.
- `test` seeds derive from a secret salt stored outside the working tree. `test` truth is not read by any scoring command until the decision rule is pre-registered in the ledger.
- Decision rule (spec §6): superior on per-allele exact sequence and non-inferior (margin 0.5 pp, one-sided 95%) on normal false-positive rate for every profile, and no increase in false NO_PATHOGENIC.
- `make ci-check` before every commit; `make docs-check` for Task 12.

## Review Focus

1. **Short allele with a large Δ in genomic mode.** Only few reads span, so the depth planner must still reach the target spanning depth on the *longer* allele. Expect realized spanning reads per allele to be recorded and to be within 30% of the target for most cases (Task 4).
2. **A mutation requested at a repeat that MucOneUp cannot mutate** (a disallowed parent unit). MucOneUp strict mode fails, or the output lacks the event. Expect the case to be rejected and logged as `design_invalid`, never scored as a normal (Task 6).
3. **A HiFi case where ccs drops molecules.** Expect realized depth to come from the truth manifest, not the requested coverage (Task 5, Task 6).
4. **An engine run that crashes or writes no `summary.json`.** Expect the case to stay in the denominator as `execution_failed`, scored as a failure for every metric and as NO_CALL in the clinical matrix. It must never be dropped (Task 8, Task 10).
5. **Opening the sealed `test` truth before pre-registration.** Expect `benchsim evaluate --split test` to refuse without a matching pre-registration entry (Task 11).

---

## File Structure

| Path | Responsibility |
| --- | --- |
| `src/muc_one_span/benchsim/__init__.py` | package docstring only |
| `src/muc_one_span/benchsim/muconeup.py` | version gate; argument builders for `simulate` / `reads` |
| `src/muc_one_span/benchsim/design.py` | `Design` model, factor levels, stratified sampler, splits, seeds |
| `src/muc_one_span/benchsim/profiles.py` | derive profile variants (artefact, error, PCR) from a built-in profile JSON |
| `src/muc_one_span/benchsim/depth.py` | template / read counts that reach a target spanning depth |
| `src/muc_one_span/benchsim/read_truth.py` | load `*_read_truth.tsv.gz`; realized depth and composition |
| `src/muc_one_span/benchsim/generate.py` | generate one case: simulate, reads, validate, `case.json`, ledger |
| `src/muc_one_span/benchsim/realism.py` | realism metrics of simulated reads vs targets |
| `src/muc_one_span/benchsim/targets/prjeb92208_v1.json` | public PRJEB92208 aggregate targets (amplicon + WGS) |
| `src/muc_one_span/benchsim/stats.py` | Clopper–Pearson, exact McNemar, Holm, non-inferiority, cluster bootstrap |
| `src/muc_one_span/benchsim/run_cases.py` | run engines per case; measurements per engine |
| `src/muc_one_span/benchsim/report.py` | stratified tables, paired tests, decision rule, pre-registration |
| `src/muc_one_span/evaluation/clinical_confusion.py` | truth clinical class; predicted decision; confusion matrix |
| `src/muc_one_span/evaluation/truth.py` | mark read-source truth available and hash the manifest |
| `scripts/benchsim.py` | CLI: `design`, `generate`, `run`, `evaluate`, `realism`, `report`, `preregister` |
| `tests/unit/benchsim/test_*.py` | one test module per library module |
| `tests/integration/test_benchsim_generate.py` | one real MucOneUp case per profile |
| `docs/benchmark.md` | user documentation |

---

### Task 1: MucOneUp version gate and argument builders

**Files:**
- Create: `src/muc_one_span/benchsim/__init__.py`, `src/muc_one_span/benchsim/muconeup.py`
- Test: `tests/unit/benchsim/__init__.py`, `tests/unit/benchsim/test_muconeup.py`

**Interfaces:**
- Produces:
  - `MIN_MUCONEUP = (0, 45, 0)`
  - `require_muconeup(executable: str) -> str`: returns the version string; raises `RuntimeError`.
  - `simulate_args(executable: str, config: Path, out_dir: Path, base: str, seed: int, lengths: tuple[int, int] | None, structure_file: Path | None, mutation: str | None, targets: tuple[tuple[int, int], ...]) -> list[str]`
  - `reads_args(executable: str, config: Path, profile: str, truth_fa: Path, out_dir: Path, base: str, seed: int, *, amount: int, profile_path: Path | None, flank_fasta: Path | None, pcr_preset: str | None) -> list[str]`, where `profile` is one of `"ont_amplicon_r10" | "ont_genomic_targeted" | "hifi_amplicon"`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/benchsim/test_muconeup.py
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.benchsim.muconeup import reads_args, require_muconeup, simulate_args

MOD = "muc_one_span.benchsim.muconeup"


@pytest.mark.parametrize("out,ok", [("muconeup, version 0.45.0\n", True),
                                    ("muconeup, version 0.44.5\n", False),
                                    ("muconeup, version 1.2.3\n", True)])
def test_version_gate(out: str, ok: bool) -> None:
    with patch(f"{MOD}.run_tool", return_value=out):
        if ok:
            assert require_muconeup("muconeup") == out.split()[-1]
        else:
            with pytest.raises(RuntimeError, match="0.45.0"):
                require_muconeup("muconeup")


def test_simulate_args_fixed_lengths_and_mutation(tmp_path: Path) -> None:
    args = simulate_args("mu", tmp_path / "c.json", tmp_path, "b", 7, (40, 61), None,
                         "dupC", ((2, 30),))
    assert args[:4] == ["mu", "--config", str(tmp_path / "c.json"), "simulate"]
    assert args.count("--fixed-lengths") == 2 and "--output-structure" in args
    assert args[args.index("--mutation-name") + 1] == "dupC"
    assert args[args.index("--mutation-targets") + 1] == "2,30"


def test_simulate_args_structure_file_excludes_lengths(tmp_path: Path) -> None:
    s = tmp_path / "s.txt"
    args = simulate_args("mu", tmp_path / "c.json", tmp_path, "b", 7, None, s, None, ())
    assert "--fixed-lengths" not in args and args[args.index("--input-structure") + 1] == str(s)


def test_reads_args_per_profile(tmp_path: Path) -> None:
    fa = tmp_path / "t.fa"
    amp = reads_args("mu", tmp_path / "c.json", "ont_amplicon_r10", fa, tmp_path, "b", 3,
                     amount=900, profile_path=None, flank_fasta=None, pcr_preset=None)
    assert amp[3:5] == ["reads", "amplicon"] and "--no-align" in amp
    assert amp[amp.index("--read-profile") + 1] == "ont_r10_sup_amplicon_v1"
    assert amp[amp.index("--coverage") + 1] == "900"
    gen = reads_args("mu", tmp_path / "c.json", "ont_genomic_targeted", fa, tmp_path, "b", 3,
                     amount=250, profile_path=tmp_path / "p.json", flank_fasta=tmp_path / "f.fa",
                     pcr_preset=None)
    assert gen[3:5] == ["reads", "ont"]
    assert gen[gen.index("--simulator") + 1] == "pbsim3-fragments"
    assert gen[gen.index("--read-profile") + 1] == str(tmp_path / "p.json")
    assert gen[gen.index("--n-reads") + 1] == "250"
    assert gen[gen.index("--flank-fasta") + 1] == str(tmp_path / "f.fa")
    hifi = reads_args("mu", tmp_path / "c.json", "hifi_amplicon", fa, tmp_path, "b", 3,
                      amount=60, profile_path=None, flank_fasta=None, pcr_preset="no_bias")
    assert hifi[hifi.index("--read-profile") + 1] == "hifi_amplicon_v1"
    assert hifi[hifi.index("--pcr-preset") + 1] == "no_bias"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --locked --all-extras pytest tests/unit/benchsim/test_muconeup.py -v`
Expected: FAIL with `ModuleNotFoundError: muc_one_span.benchsim`

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/benchsim/__init__.py
"""Realistic simulated benchmark (MucSim-Bench) built on MucOneUp >= 0.45.0."""
```

```python
# src/muc_one_span/benchsim/muconeup.py
"""Argument builders for the external MucOneUp CLI (never imported as a library)."""

from __future__ import annotations

import re
from pathlib import Path

from muc_one_span.tools import run_tool

MIN_MUCONEUP = (0, 45, 0)
BUILTIN_PROFILE = {
    "ont_amplicon_r10": "ont_r10_sup_amplicon_v1",
    "ont_genomic_targeted": "ont_r10_genomic_v1",
    "hifi_amplicon": "hifi_amplicon_v1",
}


def require_muconeup(executable: str) -> str:
    """Return the MucOneUp version; fail if it predates read profiles and read truth."""
    out = run_tool([executable, "--version"])
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    if not match or tuple(map(int, match.groups())) < MIN_MUCONEUP:
        raise RuntimeError(f"MucOneUp >= 0.45.0 is required, found: {out.strip()!r}")
    return match.group(0)


def simulate_args(
    executable: str, config: Path, out_dir: Path, base: str, seed: int,
    lengths: tuple[int, int] | None, structure_file: Path | None,
    mutation: str | None, targets: tuple[tuple[int, int], ...],
) -> list[str]:
    """`muconeup simulate` for one diploid design."""
    args = [executable, "--config", str(config), "simulate", "--out-dir", str(out_dir),
            "--out-base", base, "--num-haplotypes", "2", "--output-structure",
            "--seed", str(seed)]
    if structure_file is not None:
        args += ["--input-structure", str(structure_file)]
    elif lengths is not None:
        for length in lengths:
            args += ["--fixed-lengths", str(length)]
    else:
        raise ValueError("either lengths or structure_file is required")
    if mutation:
        args += ["--mutation-name", mutation]
        for hap, repeat in targets:
            args += ["--mutation-targets", f"{hap},{repeat}"]
    return args


def reads_args(
    executable: str, config: Path, profile: str, truth_fa: Path, out_dir: Path, base: str,
    seed: int, *, amount: int, profile_path: Path | None, flank_fasta: Path | None,
    pcr_preset: str | None,
) -> list[str]:
    """`muconeup reads …` for one benchmark profile; always FASTQ (`--no-align`)."""
    ref = str(profile_path) if profile_path else BUILTIN_PROFILE[profile]
    head = [executable, "--config", str(config), "reads"]
    tail = ["--read-profile", ref, "--out-dir", str(out_dir), "--out-base", base,
            "--seed", str(seed), "--no-align", str(truth_fa)]
    if profile == "ont_genomic_targeted":
        args = head + ["ont", "--simulator", "pbsim3-fragments", "--n-reads", str(amount)]
        if flank_fasta is not None:
            args += ["--flank-fasta", str(flank_fasta)]
        return args + tail
    platform = "ont" if profile == "ont_amplicon_r10" else "pacbio"
    args = head + ["amplicon", "--platform", platform, "--coverage", str(amount)]
    if pcr_preset:
        args += ["--pcr-preset", pcr_preset]
    return args + tail
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --locked --all-extras pytest tests/unit/benchsim/test_muconeup.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
make ci-check
git add src/muc_one_span/benchsim tests/unit/benchsim
git commit -m "feat(benchsim): MucOneUp 0.45 version gate and argument builders"
```

---

### Task 2: Design model, factor levels, stratified sampler, splits and seeds

**Files:**
- Create: `src/muc_one_span/benchsim/design.py`
- Test: `tests/unit/benchsim/test_design.py`

**Interfaces:**
- Consumes: `load_repeat_dictionary()` from `muc_one_span.config` (mutation names and `changes`).
- Produces:
  - `@dataclass(frozen=True) Design` with fields `design_id: str`, `split: str`, `profile: str`, `lengths: tuple[int, int]`, `delta_class: str`, `composition: str`, `event: str | None`, `event_allele: str | None` (`"shorter" | "longer" | "equal"`), `event_position: str | None` (`"first10" | "middle" | "last10"`), `targets: tuple[tuple[int, int], ...]`, `depth: int`, `pcr: str` (`"calibrated" | "strong" | "none"`), `smear: float`, `chimera: float`, `error: str` (`"calibrated" | "poor"`), `bio_seed: int`, `read_seed: int`. Also `to_dict()` / `from_dict()`.
  - `derive_seed(salt: str, design_id: str, stream: str) -> int` (63-bit, SHA-256)
  - `build_split(split: str, n_per_profile: int, salt: str, mutations: Sequence[str]) -> list[Design]`
  - Constants `PROFILES`, `DEPTHS`, `DELTA_CLASSES`, `NORMAL_FRACTION = 0.35`.

Factor levels are copied verbatim from spec §5. Two cases are special. For `delta_class="0_identical"`, both alleles get the same structure file, which the generator writes (Task 6). For `"0_different"`, lengths are equal and the structures come from independent draws. `composition="real_derived"` is only allowed when a local structure pool is given (Task 6), and designs record the pool index, never the sequence.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/benchsim/test_design.py
from collections import Counter

from muc_one_span.benchsim.design import (
    DEPTHS, NORMAL_FRACTION, PROFILES, Design, build_split, derive_seed,
)

MUTS = ["dupC", "dupA", "insG", "delGCCCA", "insCCC_benign"]


def test_seed_derivation_is_stable_and_stream_specific() -> None:
    a = derive_seed("s", "dev-0001", "bio")
    assert a == derive_seed("s", "dev-0001", "bio")
    assert a != derive_seed("s", "dev-0001", "reads") != derive_seed("t", "dev-0001", "bio")
    assert 0 <= a < 2**63


def test_split_sizes_normal_floor_and_profiles() -> None:
    designs = build_split("dev", 60, "salt", MUTS)
    assert len(designs) == 180 and {d.profile for d in designs} == set(PROFILES)
    for profile in PROFILES:
        rows = [d for d in designs if d.profile == profile]
        assert sum(d.event is None for d in rows) / len(rows) >= NORMAL_FRACTION
        assert all(d.depth in DEPTHS[profile] for d in rows)


def test_targets_match_event_allele_and_position() -> None:
    for d in build_split("dev", 90, "salt", MUTS):
        if d.event is None:
            assert d.targets == () and d.event_allele is None
            continue
        (hap, repeat), *_ = d.targets
        length = d.lengths[hap - 1]
        assert 1 <= repeat <= length
        if d.event_allele == "shorter" and d.lengths[0] != d.lengths[1]:
            assert length == min(d.lengths)
        if d.event_position == "first10":
            assert repeat <= max(1, length // 10) + 4  # MucOneUp allowed-parent search window


def test_designs_are_deterministic_and_round_trip() -> None:
    a, b = build_split("val", 30, "salt", MUTS), build_split("val", 30, "salt", MUTS)
    assert a == b
    assert Design.from_dict(a[0].to_dict()) == a[0]
    assert len({d.design_id for d in a}) == len(a)


def test_delta_classes_are_stratified() -> None:
    counts = Counter(d.delta_class for d in build_split("dev", 300, "salt", MUTS))
    assert min(counts.values()) >= 20
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/benchsim/test_design.py -v`
Expected: FAIL (`ModuleNotFoundError`)

- [ ] **Step 3: Implement `design.py`**

```python
# src/muc_one_span/benchsim/design.py
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
_DELTA_RANGE = {"0_identical": (0, 0), "0_different": (0, 0), "1": (1, 1), "2": (2, 2),
                "3-5": (3, 5), "6-20": (6, 20), ">20": (21, 90)}
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


def _target(lengths: tuple[int, int], allele: str, position: str,
            rng: random.Random) -> tuple[int, int]:
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


def build_split(split: str, n_per_profile: int, salt: str,
                mutations: Sequence[str]) -> list[Design]:
    """Designs for one split; every profile gets ``n_per_profile`` cases."""
    if not mutations:
        raise ValueError("at least one mutation name is required")
    designs: list[Design] = []
    for profile in PROFILES:
        rng = random.Random(derive_seed(salt, f"{split}:{profile}", "design"))
        n_normal = -(-int(n_per_profile * NORMAL_FRACTION * 100) // 100)
        events: list[str | None] = [None] * n_normal + _stratum(
            list(mutations), n_per_profile - n_normal, rng)
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
            allele = alleles[i] if event else None
            position = positions[i] if event else None
            targets = (_target(lengths, allele, position, rng),) if event else ()
            designs.append(Design(
                design_id, split, profile, lengths, deltas[i], comps[i], event, allele,
                position, targets, depths[i], pcrs[i], smears[i], chimeras[i], errors[i],
                derive_seed(salt, design_id, "bio"), derive_seed(salt, design_id, "reads"),
            ))
    return designs
```

MucOneUp picks the nearest allowed parent unit when the requested repeat is not
allowed for the mutation. The generator (Task 6) therefore records the **actual**
target from MucOneUp's output, and the test's `+ 4` window reflects that.
Designs keep the requested target; `case.json` keeps both.

- [ ] **Step 4: Run tests** — `uv run --locked --all-extras pytest tests/unit/benchsim/test_design.py -v` → PASS
- [ ] **Step 5: Commit** — `make ci-check && git add src/muc_one_span/benchsim/design.py tests/unit/benchsim/test_design.py && git commit -m "feat(benchsim): stratified design model with split-specific seeds"`

---

### Task 3: Profile variants (artefact, error and PCR levels)

**Files:**
- Create: `src/muc_one_span/benchsim/profiles.py`
- Test: `tests/unit/benchsim/test_profiles.py`

**Interfaces:**
- Consumes: `Design` (Task 2).
- Produces:
  - `variant_name(design: Design) -> str`, e.g. `ont_r10_sup_amplicon_v1__s0.25_c0.05_pcr-strong_err-poor`. The genomic profile ignores the smear and chimera levels.
  - `write_variant(base_profile: Path, design: Design, out_dir: Path) -> tuple[Path, str]`: writes the variant JSON once per name and returns `(path, sha256)`.
  - `builtin_profile_dir(explicit: Path | None) -> Path`: the explicit path, else `importlib.util.find_spec("muc_one_up")`-based `data/read_profiles`, else raise with instructions (`--muconeup-profiles DIR`).

Rules (documented in the module docstring):
- smear → `molecules.smear_rate`; chimera → `molecules.chimera_rate`
- error `poor` → multiply `errors.mismatch_rate/insertion_rate/deletion_rate` by 1.5; stutter is unchanged
- PCR `strong` → `config_overrides.amplicon_params.pcr_bias = {"preset": "madritsch2025_r10", "alpha": 2 × 9.27e-5}`
- PCR `none` → `{"preset": "no_bias"}`
- `calibrated` levels leave the base untouched
- `hifi_amplicon_v1` has no `errors`, so error `poor` for HiFi instead writes `config_overrides.pacbio_params.accuracy_mean = 0.95` (base pbsim default 0.99)

The variant `name` is set to the variant name, and `provenance.derived_from` records the base name and sha256.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/benchsim/test_profiles.py
import json
from dataclasses import replace
from pathlib import Path

from muc_one_span.benchsim.design import build_split
from muc_one_span.benchsim.profiles import variant_name, write_variant

BASE = {
    "schema_version": 1, "name": "ont_r10_sup_amplicon_v1", "platform": "ont",
    "config_overrides": {"amplicon_params": {"pcr_bias": {"preset": "madritsch2025_r10"}}},
    "molecules": {"forward_frac": 0.5, "smear_rate": 0.24, "chimera_rate": 0.023},
    "errors": {"mismatch_rate": 0.007, "insertion_rate": 0.006, "deletion_rate": 0.008,
               "insertion_len_pmf": {"1": 1.0}, "deletion_len_pmf": {"1": 1.0}},
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
```

- [ ] **Step 2: Run to verify failure** (`ModuleNotFoundError`).
- [ ] **Step 3: Implement**

```python
# src/muc_one_span/benchsim/profiles.py
"""Profile variants: artefact, error and PCR levels layered on a built-in MucOneUp profile."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

from .design import Design
from .muconeup import BUILTIN_PROFILE

ALPHA_R10 = 9.27e-5
POOR_ERROR_SCALE = 1.5


def builtin_profile_dir(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit
    spec = importlib.util.find_spec("muc_one_up")
    if spec is not None and spec.origin:
        candidate = Path(spec.origin).parent / "data" / "read_profiles"
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        "MucOneUp read profiles not found; pass --muconeup-profiles "
        "<MucOneUp checkout>/muc_one_up/data/read_profiles"
    )


def variant_name(design: Design) -> str:
    base = BUILTIN_PROFILE[design.profile]
    parts = [f"pcr-{design.pcr}", f"err-{design.error}"]
    if design.profile != "ont_genomic_targeted":
        parts = [f"s{design.smear}", f"c{design.chimera}", *parts]
    return f"{base}__{'_'.join(parts)}"


def _apply(data: dict[str, Any], design: Design) -> None:
    mol = data.setdefault("molecules", {})
    if design.profile != "ont_genomic_targeted":
        mol["smear_rate"] = design.smear
        mol["chimera_rate"] = design.chimera
    if design.error == "poor":
        if "errors" in data:
            for key in ("mismatch_rate", "insertion_rate", "deletion_rate"):
                data["errors"][key] = data["errors"][key] * POOR_ERROR_SCALE
        else:
            data.setdefault("config_overrides", {}).setdefault("pacbio_params", {})[
                "accuracy_mean"] = 0.95
    if design.profile != "ont_genomic_targeted" and design.pcr != "calibrated":
        pcr = ({"preset": "no_bias"} if design.pcr == "none"
               else {"preset": "madritsch2025_r10", "alpha": 2 * ALPHA_R10})
        data.setdefault("config_overrides", {}).setdefault("amplicon_params", {})["pcr_bias"] = pcr


def write_variant(base_profile: Path, design: Design, out_dir: Path) -> tuple[Path, str]:
    raw = base_profile.read_bytes()
    base = json.loads(raw)
    data = copy.deepcopy(base)
    _apply(data, design)  # calibrated levels still get a variant file, for provenance
    data["name"] = variant_name(design)
    data.setdefault("provenance", {})["derived_from"] = {
        "name": base["name"], "sha256": hashlib.sha256(raw).hexdigest()}
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{data['name']}.json"
    text = json.dumps(data, indent=2, sort_keys=True) + "\n"
    if not path.exists() or path.read_text() != text:
        path.write_text(text)
    return path, hashlib.sha256(text.encode()).hexdigest()
```

- [ ] **Step 4: Run tests** → PASS
- [ ] **Step 5: Commit** — `git commit -m "feat(benchsim): derive artefact, error and PCR profile variants"`

---

### Task 4: Depth planning

**Files:**
- Create: `src/muc_one_span/benchsim/depth.py`
- Test: `tests/unit/benchsim/test_depth.py`

**Interfaces:**
- Produces:
  - `amplicon_templates(target_full_per_allele: int, artefact_rate: float, minor_share: float) -> int`. It returns the `--coverage` (template molecules) needed so that the minor allele gets the target number of full-length reads. `artefact_rate = smear + chimera + concatemer`; `minor_share` is the expected PCR share of the minor allele.
  - `pcr_minor_share(lengths: tuple[int, int], pcr: str) -> float`: slope −0.056/unit for calibrated, −0.112 for strong, and 0 for none, so `share = 1 / (1 + exp(0.056 * Δ))`.
  - `genomic_reads(target_spanning: int, source_len: int, span_start: int, span_end: int, median: float, sigma: float, n_hap: int = 2, draws: int = 20000, seed: int = 0) -> int`. It uses the MucOneUp ≥ 0.45 fragment model: start uniform in `[-(L-1), S)` and clipped at both ends. A fragment spans when `start ≤ span_start` and `start + L ≥ span_end`. With `p` the Monte-Carlo spanning probability, it returns `ceil(target * n_hap / p)`.

- [ ] **Step 1: Failing tests**

```python
# tests/unit/benchsim/test_depth.py
import math

import pytest

from muc_one_span.benchsim.depth import amplicon_templates, genomic_reads, pcr_minor_share


def test_minor_share() -> None:
    assert pcr_minor_share((40, 40), "calibrated") == 0.5
    assert pcr_minor_share((40, 60), "none") == 0.5
    assert pcr_minor_share((40, 60), "calibrated") == pytest.approx(1 / (1 + math.exp(0.056 * 20)))
    assert pcr_minor_share((40, 60), "strong") < pcr_minor_share((40, 60), "calibrated")


def test_amplicon_templates_reach_target() -> None:
    n = amplicon_templates(30, 0.287, 0.25)
    assert n * 0.25 * (1 - 0.287) >= 30 > (n - 1) * 0.25 * (1 - 0.287)


def test_genomic_reads_hits_expected_spanning() -> None:
    import random
    n = genomic_reads(20, 26000, 10000, 16000, 6000, 0.5, n_hap=2, seed=1)
    rng, spans = random.Random(2), 0
    for _ in range(n):
        length = max(1, round(rng.lognormvariate(math.log(6000), 0.5)))
        start = rng.randrange(-(length - 1), 26000)
        spans += start <= 10000 and start + length >= 16000
    assert spans / 2 == pytest.approx(20, rel=0.3)


def test_genomic_reads_rejects_impossible_span() -> None:
    with pytest.raises(ValueError, match="cannot span"):
        genomic_reads(10, 30000, 0, 30000, 500, 0.1, seed=1)
```

- [ ] **Step 2: Verify failure.**
- [ ] **Step 3: Implement**

```python
# src/muc_one_span/benchsim/depth.py
"""Template and read counts that reach a target spanning depth per allele."""

from __future__ import annotations

import math
import random

SLOPE = {"calibrated": 0.056, "strong": 0.112, "none": 0.0}


def pcr_minor_share(lengths: tuple[int, int], pcr: str) -> float:
    delta = abs(lengths[0] - lengths[1])
    return 1.0 / (1.0 + math.exp(SLOPE[pcr] * delta))


def amplicon_templates(target_full_per_allele: int, artefact_rate: float,
                       minor_share: float) -> int:
    if not 0 < minor_share <= 0.5 or not 0 <= artefact_rate < 1:
        raise ValueError("minor_share in (0, 0.5] and artefact_rate in [0, 1) required")
    return math.ceil(target_full_per_allele / (minor_share * (1.0 - artefact_rate)))


def genomic_reads(target_spanning: int, source_len: int, span_start: int, span_end: int,
                  median: float, sigma: float, n_hap: int = 2, draws: int = 20000,
                  seed: int = 0) -> int:
    rng = random.Random(seed)
    hits = 0
    for _ in range(draws):
        length = max(1, round(rng.lognormvariate(math.log(median), sigma)))
        start = rng.randrange(-(length - 1), source_len)
        hits += start <= span_start and start + length >= span_end
    if hits == 0:
        raise ValueError("fragment length model cannot span the VNTR; lengthen reads or flanks")
    return math.ceil(target_spanning * n_hap * draws / hits)
```

In the generator (Task 6), `span_start`/`span_end` are the anchor-to-anchor interval
in the fragment source: the left flank-FASTA length, plus MucOneUp's left flank, plus
motif 1 start through motif 9 end. Use the longer allele's source, so both alleles
reach at least the target.

- [ ] **Step 4: Run tests** → PASS
- [ ] **Step 5: Commit** — `git commit -m "feat(benchsim): depth planning for amplicon templates and genomic fragments"`

---

### Task 5: Read-truth loader, realized depth, truth adapter flag

**Files:**
- Create: `src/muc_one_span/benchsim/read_truth.py`
- Modify: `src/muc_one_span/evaluation/truth.py`, in `_load`: replace the line `provenance["read_source_truth"] = "unavailable"`
- Test: `tests/unit/benchsim/test_read_truth.py`, `tests/unit/test_evaluation_truth.py` (append)

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) ReadTruth(read_id: str, hap: int, molecule: int, kind: str, strand: str, src_start: int, src_end: int, hp_edits: str, detail: str)`
  - `load_read_truth(path: Path) -> list[ReadTruth]`: validates the header against the MucOneUp 0.45 columns (`read_id hap molecule kind strand src_start src_end n_hp_edits hp_edits detail`); raises `ValueError` on duplicates or a bad header.
  - `realized_depth(rows, span: dict[int, tuple[int, int]] | None) -> dict[int, int]`. For amplicons (`span=None`) it counts `kind == "full"` per hap. For fragments it counts rows whose `[src_start, src_end)` covers the hap's span interval.
  - `composition(rows) -> dict[str, Any]`: kind fractions, strand fractions, reads per hap.
  - In `evaluation/truth.py`, `provenance["read_source_truth"]` is `"available"` when exactly one `*_read_truth.tsv.gz` exists in the sample dir (its hash is added to `hashes`), and `"unavailable"` otherwise.

- [ ] **Step 1: Failing tests**

```python
# tests/unit/benchsim/test_read_truth.py
import gzip
from pathlib import Path

import pytest

from muc_one_span.benchsim.read_truth import composition, load_read_truth, realized_depth

HEADER = "read_id\thap\tmolecule\tkind\tstrand\tsrc_start\tsrc_end\tn_hp_edits\thp_edits\tdetail\n"


def _write(path: Path, rows: list[str], header: str = HEADER) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write(header + "".join(r + "\n" for r in rows))
    return path


def test_load_and_amplicon_depth(tmp_path: Path) -> None:
    p = _write(tmp_path / "t.tsv.gz", [
        "b_h1_m0000001\t1\t1\tfull\t+\t0\t2500\t0\t\t",
        "b_h2_m0000002\t2\t2\tsmear\t-\t0\t2500\t0\t\tdeletion:100-900",
        "b_h2_m0000003\t2\t3\tfull\t-\t0\t2600\t1\t52:C:7>8\t",
        "b_h1_m0000004\t1\t4\tofftarget\t+\t10\t300\t0\t\t",
    ])
    rows = load_read_truth(p)
    assert realized_depth(rows, None) == {1: 1, 2: 1}
    comp = composition(rows)
    assert comp["kind_frac"]["offtarget"] == 0.25 and comp["strand_frac"]["-"] == 0.5


def test_fragment_depth_requires_covering_span(tmp_path: Path) -> None:
    p = _write(tmp_path / "t.tsv.gz", [
        "b_h1_m0000001\t1\t1\tfragment\t+\t900\t7000\t0\t\t",
        "b_h1_m0000002\t1\t2\tfragment\t+\t1100\t7000\t0\t\t",
    ])
    assert realized_depth(load_read_truth(p), {1: (1000, 5000), 2: (1000, 6000)}) == {1: 1, 2: 0}


def test_bad_header_and_duplicates(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="header"):
        load_read_truth(_write(tmp_path / "a.tsv.gz", [], header="read_id\thap\n"))
    row = "b_h1_m0000001\t1\t1\tfull\t+\t0\t10\t0\t\t"
    with pytest.raises(ValueError, match="duplicate"):
        load_read_truth(_write(tmp_path / "b.tsv.gz", [row, row]))
```

Also append to `tests/unit/test_evaluation_truth.py`. Reuse the module's existing
fixture that builds a valid MucOneUp sample dir; call it `valid_sample_dir`, or
whatever the module's name is, at implementation time.

```python
def test_read_truth_manifest_marked_available(valid_sample_dir, repeat_dict):
    import gzip
    with gzip.open(valid_sample_dir / "x_read_truth.tsv.gz", "wt") as fh:
        fh.write("read_id\thap\tmolecule\tkind\tstrand\tsrc_start\tsrc_end\tn_hp_edits\thp_edits\tdetail\n")
    truth = load_truth(valid_sample_dir, repeat_dict)
    assert truth.provenance["read_source_truth"] == "available"
    assert "x_read_truth.tsv.gz" in truth.hashes
```

- [ ] **Step 2: Verify failures.**
- [ ] **Step 3: Implement**

```python
# src/muc_one_span/benchsim/read_truth.py
"""MucOneUp >= 0.45 per-read truth manifests ({base}_read_truth.tsv.gz)."""

from __future__ import annotations

import gzip
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

COLUMNS = ("read_id", "hap", "molecule", "kind", "strand", "src_start", "src_end",
           "n_hp_edits", "hp_edits", "detail")


@dataclass(frozen=True)
class ReadTruth:
    read_id: str
    hap: int
    molecule: int
    kind: str
    strand: str
    src_start: int
    src_end: int
    hp_edits: str
    detail: str


def load_read_truth(path: Path) -> list[ReadTruth]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if tuple(header) != COLUMNS:
            raise ValueError(f"{path}: unexpected read truth header {header}")
        rows, seen = [], set()
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] in seen:
                raise ValueError(f"{path}: duplicate read id {f[0]}")
            seen.add(f[0])
            rows.append(ReadTruth(f[0], int(f[1]), int(f[2]), f[3], f[4], int(f[5]),
                                  int(f[6]), f[8], f[9]))
    return rows


def realized_depth(rows: Sequence[ReadTruth],
                   span: dict[int, tuple[int, int]] | None) -> dict[int, int]:
    haps = sorted({r.hap for r in rows} | set(span or {}))
    depth = dict.fromkeys(haps, 0)
    for r in rows:
        if span is None:
            depth[r.hap] += r.kind == "full"
        elif r.kind == "fragment" and r.hap in span:
            lo, hi = span[r.hap]
            depth[r.hap] += r.src_start <= lo and r.src_end >= hi
    return depth


def composition(rows: Sequence[ReadTruth]) -> dict[str, Any]:
    n = len(rows) or 1
    kinds, strands, haps = (Counter(r.kind for r in rows), Counter(r.strand for r in rows),
                            Counter(r.hap for r in rows))
    return {"n_reads": len(rows),
            "kind_frac": {k: v / n for k, v in sorted(kinds.items())},
            "strand_frac": {k: v / n for k, v in sorted(strands.items())},
            "reads_per_hap": dict(sorted(haps.items()))}
```

In `evaluation/truth.py` `_load`, replace `provenance["read_source_truth"] = "unavailable"` with:

```python
    manifests = sorted(root.glob("*_read_truth.tsv.gz"))
    provenance["read_source_truth"] = "available" if len(manifests) == 1 else "unavailable"
    if len(manifests) == 1:
        paths.append(manifests[0])
```

(`paths` feeds `hashes` at the end of `_load`, so the manifest is hashed automatically.)

- [ ] **Step 4: Run tests** → PASS (`make test-fast`)
- [ ] **Step 5: Commit** — `git commit -m "feat(benchsim): read truth loader and realized depth; mark read truth available"`

---

### Task 6: Case generation, validation, manifest and ledger

**Files:**
- Create: `src/muc_one_span/benchsim/generate.py`
- Create: `scripts/benchsim.py` (subcommands `design` and `generate` here; others added later)
- Test: `tests/unit/benchsim/test_generate.py`, `tests/integration/test_benchsim_generate.py`

**Interfaces:**
- Consumes: Tasks 1–5; `load_truth`, `load_repeat_dictionary`; `DurableLedger`, `LedgerEntry`, `compute_sha256` from `durable_ledger.py`.
- Produces:
  - `@dataclass GenerateContext(executable: str, config: Path, out_root: Path, profile_dir: Path, flank_fasta: Path | None, structure_pool: Path | None, muconeup_version: str)`
  - `generate_case(design: Design, ctx: GenerateContext) -> dict[str, Any]`. It returns the `case.json` dict, with `status` one of `"ok" | "design_invalid" | "generation_failed"`.
  - Case layout: `<out_root>/<split>/<design_id>/{truth/, reads/, case.json}`. `truth/` holds the MucOneUp `simulate` outputs plus a copy of the reads' `*_read_truth.tsv.gz` and `*_metadata.tsv`, so that `load_truth(truth_dir)` sees one metadata TSV and one manifest. `reads/` holds the FASTQ.
  - `write_manifest(split_dir: Path, cases: list[dict]) -> Path` → `manifest.jsonl`, sorted by `design_id`.

Behaviour:
1. `require_muconeup` once per run (in the CLI); `ctx.muconeup_version` goes into every `case.json`.
2. `simulate`:
   - `composition="real_derived"` requires `ctx.structure_pool`; the design is then built from pool entry `bio_seed % len(pool)`, and only its sha256 is recorded.
   - `rare_units`: the structure file is written by drawing ~10% units from rarely used dictionary units.
   - `0_identical` uses one structure file for both alleles.
3. Validation: `load_truth(truth_dir)` must succeed, and the truth events must equal the design event (name) on the design haplotype. Record the actual `(hap, repeat)`. On failure → `design_invalid` with the error. Never score it as normal.
4. Reads:
   - Amplicon: `amount = amplicon_templates(depth, smear + chimera + base concatemer, pcr_minor_share(...))`.
   - Genomic: `amount = genomic_reads(depth, …)`, using the source length and the span interval computed from the truth FASTA and flank FASTA.
   - Profile variant path from Task 3.
5. Realized depth and composition come from the manifest; `case.json` stores the design, actual targets, requested amount, realized depth, composition, profile variant name and sha256, the MucOneUp version, and the hashes of the truth FASTA, FASTQ and manifest.
6. Ledger: `DurableLedger(out_root / split).commit_entry(LedgerEntry(…))`, with `token=design_id`, `category=f"{profile}:{delta_class}:{event or 'normal'}"` and `platform=profile`. Rerunning a verified-complete case is a no-op (`is_verified_complete`).

- [ ] **Step 1: Failing unit tests (run_tool mocked; fake MucOneUp writes minimal outputs)**

```python
# tests/unit/benchsim/test_generate.py
import json
from pathlib import Path
from unittest.mock import patch

from muc_one_span.benchsim.design import build_split
from muc_one_span.benchsim.generate import GenerateContext, generate_case

MOD = "muc_one_span.benchsim.generate"


def _ctx(tmp_path: Path) -> GenerateContext:
    prof = tmp_path / "profiles"
    prof.mkdir()
    for name in ("ont_r10_sup_amplicon_v1", "ont_r10_genomic_v1", "hifi_amplicon_v1"):
        (prof / f"{name}.json").write_text(json.dumps(
            {"schema_version": 1, "name": name, "platform": "ont", "molecules": {}}))
    return GenerateContext("muconeup", tmp_path / "c.json", tmp_path / "out", prof, None, None,
                           "0.45.0")


def test_invalid_design_is_recorded_not_scored(tmp_path: Path) -> None:
    design = next(d for d in build_split("dev", 30, "s", ["dupC"]) if d.event)
    with patch(f"{MOD}.run_tool"), patch(f"{MOD}.load_truth", side_effect=ValueError("no event")):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "design_invalid" and "no event" in case["error"]
    saved = json.loads((tmp_path / "out/dev" / design.design_id / "case.json").read_text())
    assert saved["status"] == "design_invalid"


def test_event_mismatch_is_design_invalid(tmp_path: Path) -> None:
    from muc_one_span.evaluation.models import Event, TruthHaplotype, TruthSample
    design = next(d for d in build_split("dev", 30, "s", ["dupC"]) if d.event)
    truth = TruthSample("x", (TruthHaplotype("haplotype_1", "A", ("X",)),
                              TruthHaplotype("haplotype_2", "A", ("X",))))
    with patch(f"{MOD}.run_tool"), patch(f"{MOD}.load_truth", return_value=truth):
        case = generate_case(design, _ctx(tmp_path))
    assert case["status"] == "design_invalid" and "event" in case["error"]
```

Also add a test in which the fake `run_tool` creates:
- `truth/` with the three MucOneUp truth files, reused from the `tests/unit/test_evaluation_truth.py` fixture helper;
- `reads/b_amplicon_ont.fastq`;
- `reads/b_read_truth.tsv.gz` (2 full reads per hap).

Assert:
- `status == "ok"`, `realized_depth == {"1": 2, "2": 2}`, and `profile_sha256` has length 64;
- the argument list passed for reads contains `--read-profile` pointing into `out/profiles/`;
- `ledger_sealed.jsonl` has one line, and a second call does not call `run_tool` again.

- [ ] **Step 2: Verify failures.**
- [ ] **Step 3: Implement `generate.py` (≤ 400 lines) and the `design` / `generate` subcommands.** The `design` subcommand writes `designs_<split>.jsonl`; `--salt-file` is mandatory for `test` and must resolve outside the repository (`Path.resolve()` not under `git rev-parse --show-toplevel`). `generate` takes `--designs`, `--jobs N` (process pool; each case independent), `--muconeup`, `--muconeup-config`, `--muconeup-profiles`, `--flank-fasta`, `--structure-pool`, and `--out-root` (default `../MucOneSpan-bench-data`).
- [ ] **Step 4: Integration test** (`@pytest.mark.integration`; skip unless `muconeup` ≥ 0.45, `pbsim`, `ccs` are on PATH and `MUCONEUP_CONFIG` is set): generate one `dev` case per profile at depth 10, then assert status `ok`, a non-empty FASTQ, manifest rows equal to FASTQ records, and `load_truth` passing.
- [ ] **Step 5: Run** `make test-fast` and `make test-int` (report skips). **Commit** — `git commit -m "feat(benchsim): generate validated cases with read truth, manifest and ledger"`

---

### Task 7: Clinical-decision confusion matrix in `evaluation/`

**Files:**
- Create: `src/muc_one_span/evaluation/clinical_confusion.py`
- Modify: `src/muc_one_span/evaluation/__init__.py` (export), `scripts/evaluate.py` (add `clinical` block to each row and `clinical_confusion` to the report)
- Test: `tests/unit/test_clinical_confusion.py`

**Interfaces:**
- Consumes: `TruthSample`, `RepeatDictionary` (`rd.mutations[name]["changes"]`), `compute_clinical_decision` from `report.py`.
- Produces:
  - `truth_class(truth: TruthSample, rd: RepeatDictionary) -> str`, one of `"pathogenic" | "benign" | "normal"`. An event is pathogenic when its net length change `Σ len(insert) − Σ deleted` is not divisible by 3 (frameshift). Truth is pathogenic if any event is, benign if only in-frame events exist, and normal otherwise.
  - `predicted_decision(result_dir: Path) -> str`, one of `"PATHOGENIC" | "INCONCLUSIVE" | "NO_PATHOGENIC_VARIANT_DETECTED" | "NO_CALL"`. It calls `compute_clinical_decision(json.loads(summary.json))["state"]`; a missing or invalid `summary.json` gives `NO_CALL`.
  - `confusion(rows: Sequence[dict]) -> dict`: counts by `(truth_class, decision)`, plus:
    - `critical_false_negative` (pathogenic → NO_PATHOGENIC or NO_CALL)
    - `false_positive_normal` and `false_positive_benign` (→ PATHOGENIC)
    - `inconclusive_rate`

- [ ] **Step 1: Failing tests**

```python
# tests/unit/test_clinical_confusion.py
import json
from pathlib import Path

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.clinical_confusion import confusion, predicted_decision, truth_class
from muc_one_span.evaluation.models import Event, TruthHaplotype, TruthSample

RD = load_repeat_dictionary()


def _truth(*names: str) -> TruthSample:
    events = tuple(Event(5, "X", n) for n in names)
    return TruthSample("s", (TruthHaplotype("h1", "A", ("X",), events),
                             TruthHaplotype("h2", "A", ("X",))))


def test_truth_class_by_frameshift() -> None:
    assert truth_class(_truth(), RD) == "normal"
    assert truth_class(_truth("dupC"), RD) == "pathogenic"
    benign = next(n for n, d in RD.mutations.items()
                  if sum(len(c.get("sequence", "")) - (c.get("end", c["start"]) - c["start"])
                         for c in d["changes"]) % 3 == 0)
    assert truth_class(_truth(benign), RD) == "benign"


def test_missing_summary_is_no_call(tmp_path: Path) -> None:
    assert predicted_decision(tmp_path) == "NO_CALL"
    (tmp_path / "summary.json").write_text("{not json")
    assert predicted_decision(tmp_path) == "NO_CALL"


def test_confusion_counts_critical_errors() -> None:
    rows = [{"clinical": {"truth": "pathogenic", "decision": "NO_CALL"}},
            {"clinical": {"truth": "pathogenic", "decision": "PATHOGENIC"}},
            {"clinical": {"truth": "normal", "decision": "PATHOGENIC"}},
            {"clinical": {"truth": "normal", "decision": "INCONCLUSIVE"}}]
    c = confusion(rows)
    assert c["critical_false_negative"] == 1 and c["false_positive_normal"] == 1
    assert c["inconclusive_rate"] == 0.25
    assert c["matrix"]["pathogenic"]["NO_CALL"] == 1
```

Check the dictionary's change schema (`type`, `start`, `end`, `sequence`) in
`src/muc_one_span/config.py::_apply_mutation` before implementing, and compute the
net change the same way `_apply_mutation` applies it. The test's benign lookup must
use the same helper once it exists: replace the inline expression with
`net_length_change(d)` from the new module.

- [ ] **Step 2: Verify failure.** **Step 3: Implement** (≤ 150 lines; `net_length_change(definition) -> int` is public). In `scripts/evaluate.py`, add per row `row["clinical"] = {"truth": truth_class(truth, rd), "decision": predicted_decision(result_dir)}` for scored rows, and `{"truth": None, "decision": "NO_CALL"}` otherwise; then `report["clinical_confusion"] = confusion(rows)`.
- [ ] **Step 4: Run** `make test-fast` → PASS. **Step 5: Commit** — `git commit -m "feat(evaluation): clinical decision confusion matrix against simulated truth"`

---

### Task 8: Statistics

**Files:**
- Create: `src/muc_one_span/benchsim/stats.py`
- Test: `tests/unit/benchsim/test_stats.py`

**Interfaces:**
- Produces:
  - `clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]`
  - `mcnemar_exact(b: int, c: int) -> float` (two-sided exact binomial, `b`/`c` discordant counts)
  - `holm(pvalues: dict[str, float]) -> dict[str, float]` (adjusted, monotone, capped at 1)
  - `noninferior(fp_new: int, n_new: int, fp_ref: int, n_ref: int, margin: float = 0.005, alpha: float = 0.05) -> dict` (Newcombe hybrid-score upper bound of the difference; `{"diff", "upper", "noninferior"}`)
  - `cluster_bootstrap(rows: Sequence[dict], key: str, value: Callable[[dict], float], n: int = 2000, seed: int = 0) -> tuple[float, float, float]`

Pure Python with `math`, plus `scipy.stats.beta` if scipy is already a dependency. Otherwise implement Clopper–Pearson through a bisection on the regularized incomplete beta with `math.lgamma` (check `pyproject.toml` first; do not add scipy for this).

- [ ] **Step 1: Failing tests** (reference values from R `binom.test` / `mcnemar.exact`):

```python
# tests/unit/benchsim/test_stats.py
import pytest

from muc_one_span.benchsim.stats import (
    clopper_pearson, cluster_bootstrap, holm, mcnemar_exact, noninferior,
)


def test_clopper_pearson_reference() -> None:
    lo, hi = clopper_pearson(59, 59)
    assert lo == pytest.approx(0.9394, abs=1e-4) and hi == 1.0
    assert clopper_pearson(0, 280)[1] == pytest.approx(0.01308, abs=1e-4)
    lo, hi = clopper_pearson(7, 20)
    assert (lo, hi) == pytest.approx((0.1539, 0.5922), abs=1e-4)


def test_mcnemar_exact() -> None:
    assert mcnemar_exact(0, 0) == 1.0
    assert mcnemar_exact(1, 9) == pytest.approx(0.02148, abs=1e-5)


def test_holm_monotone() -> None:
    adj = holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert adj == pytest.approx({"a": 0.03, "b": 0.06, "c": 0.06})


def test_noninferiority() -> None:
    assert noninferior(0, 280, 0, 280)["noninferior"] is False  # CI too wide at n=280
    assert noninferior(0, 2000, 0, 2000)["noninferior"] is True


def test_cluster_bootstrap_groups() -> None:
    rows = [{"g": g, "ok": 1 if g % 2 else 0} for g in range(20) for _ in range(2)]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1)
    assert mean == 0.5 and lo < 0.5 < hi
```

(Verify the numbers against R once while implementing; the ones above are
`binom.test(59,59)`, `binom.test(0,280)`, `binom.test(7,20)`, and
`binom.test(1,10)$p.value`. If R disagrees in the 4th decimal, update the test from R.)

- [ ] **Step 2–4:** implement and pass. **Step 5: Commit** — `git commit -m "feat(benchsim): exact CIs, McNemar, Holm, non-inferiority, cluster bootstrap"`

---

### Task 9: Realism report and public targets

**Files:**
- Create: `src/muc_one_span/benchsim/realism.py`, `src/muc_one_span/benchsim/targets/prjeb92208_v1.json`
- Create (outside Git, run once): the export command below
- Test: `tests/unit/benchsim/test_realism.py`

**Interfaces:**
- Consumes: `load_read_truth`, `ReadTruth`; the truth FASTA (haplotype sequences) and the reads FASTQ of a case.
- Produces:
  - `read_metrics(fastq: Path, truth_rows: Sequence[ReadTruth], sources: dict[int, str]) -> dict`. For every `full`/`fragment` read, align it to its source interval (reverse-complemented for `-`) with edlib (`task="path"`). It returns:
    - mismatch, insertion and deletion rates per strand;
    - the C7 correct-length fraction per strand (runs of exactly 7 C in the source; the read length is measured through the alignment);
    - `smear_frac` among amplicon products;
    - `chimera_frac`;
    - `offtarget_frac`;
    - `log_ratio_slope`: the per-case `ln(n_long/n_short)/Δ`, used across cases.
  - `compare(metrics: dict, targets: dict, profile: str) -> dict`: pass/fail per spec §4 tolerance (±20% relative error rates, ±0.03 C7, smear/offtarget within the real range, allele slope ±0.01).
  - `targets/prjeb92208_v1.json`: only the `ont_amplicon_PRJEB92208` and `ont_wgs_PRJEB92208` sections, plus `_meta`, from the local `realprofile/targets.json`. **Never `ont_genomic_inhouse`.** In-house genomic targets are passed as `--targets-local PATH` at run time.

Export (run once, review the diff, then commit the JSON):

```bash
uv run --locked --all-extras python - <<'EOF'
import json, pathlib
src = json.loads(pathlib.Path("../MucOneSpan-review-20260923/realprofile/targets.json").read_text())
keep = {k: src[k] for k in ("_meta", "ont_amplicon_PRJEB92208", "ont_wgs_PRJEB92208")}
keep["_meta"] = {**keep["_meta"], "source": "ENA PRJEB92208 public libraries (aggregates only)",
                 "excluded": "in-house genomic targets are local-only"}
out = pathlib.Path("src/muc_one_span/benchsim/targets/prjeb92208_v1.json")
out.write_text(json.dumps(keep, indent=1, sort_keys=True) + "\n")
EOF
grep -c inhouse src/muc_one_span/benchsim/targets/prjeb92208_v1.json  # must print 0
```

Also add `"src/muc_one_span/benchsim/targets/*.json"` to package data. It is data,
so it falls outside the file-size gate.

- [ ] **Step 1: Failing tests** — synthetic source with a C7, reads built from it with known edits:

```python
# tests/unit/benchsim/test_realism.py
from pathlib import Path

from muc_one_span.benchsim.read_truth import ReadTruth
from muc_one_span.benchsim.realism import compare, read_metrics

SRC = "ACGTACGTAA" + "C" * 7 + "GATTACAGATTACA" * 4


def _fq(path: Path, reads: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in reads))
    return path


def test_metrics_per_strand(tmp_path: Path) -> None:
    rc = SRC.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    short_c = SRC.replace("C" * 7, "C" * 6)
    rows = [ReadTruth("r1", 1, 1, "full", "+", 0, len(SRC), "", ""),
            ReadTruth("r2", 1, 2, "full", "+", 0, len(SRC), "", ""),
            ReadTruth("r3", 1, 3, "full", "-", 0, len(SRC), "", "")]
    fq = _fq(tmp_path / "r.fq", [("r1", SRC), ("r2", short_c), ("r3", rc)])
    m = read_metrics(fq, rows, {1: SRC})
    assert m["c7_correct"]["+"] == 0.5 and m["c7_correct"]["-"] == 1.0
    assert m["error_rate"]["-"] == 0.0 and m["deletion_rate"]["+"] > 0


def test_compare_flags_out_of_tolerance() -> None:
    targets = {"ont_amplicon_PRJEB92208": {"c7_correct": {"+": 0.52, "-": 0.89},
                                           "error_rate": {"median": 0.021}}}
    res = compare({"c7_correct": {"+": 0.60, "-": 0.88}, "error_rate": {"all": 0.022}},
                  targets, "ont_amplicon_r10")
    assert res["c7_correct_+"]["pass"] is False and res["c7_correct_-"]["pass"] is True
    assert res["error_rate"]["pass"] is True
```

Before writing `compare`, open `targets.json` and map its actual key names; the
test's target shape is illustrative. Adjust the test to the real keys (for example,
`c7_correct_by_strand`), not the reverse.

- [ ] **Step 2–4:** implement (edlib lazy import with a clear error if missing) and pass.
- [ ] **Step 5: Commit** — `git commit -m "feat(benchsim): realism metrics and public PRJEB92208 targets"`

---

### Task 10: Engine runner over cases

**Files:**
- Create: `src/muc_one_span/benchsim/run_cases.py`; add `run` subcommand to `scripts/benchsim.py`
- Test: `tests/unit/benchsim/test_run_cases.py`

**Interfaces:**
- Consumes: `benchmarking.run_pipeline(sample, input_path, output_dir, platform, model, threads, *, runner=None, engine="ladder")`. The `engine` keyword comes from hybrid-engine plan Task 10. **Precondition:** that task is merged. If it is not, this task adds only `engine` forwarding (a keyword appended as `--engine <engine>` to `cli_args` and recorded in `measurement.json`), with a test, and hybrid Task 10 then reuses it.
- Produces:
  - `run_split(manifest: Path, engines: Sequence[str], results_root: Path, model_for: Callable[[str], str], threads: int, jobs: int) -> list[dict]`
  - Output layout: `<results_root>/<engine>/<design_id>/` plus `<results_root>/<engine>/measurements.json` (the format `scripts/evaluate.py` already reads)
  - Platform mapping: `ont_*` → `ont`, `hifi_amplicon` → `hifi`
  - `design_invalid` and `generation_failed` cases are written as `not_attempted` records, keeping the denominator
  - `inventory_<split>.json` for `scripts/evaluate.py`, with `truth_dir` and `result_dir` per case

- [ ] **Step 1: Failing tests** — fake runner; one ok case, one `design_invalid`, one run whose output lacks `summary.json`:

```python
# tests/unit/benchsim/test_run_cases.py
import json
from pathlib import Path
from unittest.mock import patch

from muc_one_span.benchsim.run_cases import run_split


def _manifest(tmp_path: Path) -> Path:
    rows = [{"design_id": "dev-a", "profile": "ont_amplicon_r10", "status": "ok",
             "reads": str(tmp_path / "a.fastq"), "truth_dir": str(tmp_path / "ta")},
            {"design_id": "dev-b", "profile": "hifi_amplicon", "status": "design_invalid"}]
    (tmp_path / "a.fastq").write_text("@r\nA\n+\nI\n")
    p = tmp_path / "manifest.jsonl"
    p.write_text("".join(json.dumps(r) + "\n" for r in rows))
    return p


def test_denominator_kept_and_engine_forwarded(tmp_path: Path) -> None:
    calls = []

    def fake(sample, input_path, output_dir, platform, model, threads, *, runner=None,
             engine="ladder"):
        calls.append((sample, platform, engine))
        return {"sample": sample, "status": "completed", "result_dir": str(output_dir)}

    with patch("muc_one_span.benchsim.run_cases.run_pipeline", side_effect=fake):
        records = run_split(_manifest(tmp_path), ["ladder", "hybrid"], tmp_path / "res",
                            lambda p: "r1041", 1, 1)
    assert sorted(calls) == [("dev-a", "ont", "hybrid"), ("dev-a", "ont", "ladder")]
    statuses = {(r["engine"], r["sample"]): r["status"] for r in records}
    assert statuses[("ladder", "dev-b")] == "not_attempted"
    assert json.loads((tmp_path / "res/ladder/measurements.json").read_text())
```

- [ ] **Step 2–4:** implement and pass. **Step 5: Commit** — `git commit -m "feat(benchsim): run engines over a split without dropping denominators"`

---

### Task 11: Report, paired comparison, decision rule and pre-registration

**Files:**
- Create: `src/muc_one_span/benchsim/report.py`; add `evaluate`, `realism`, `report`, `preregister` subcommands to `scripts/benchsim.py`
- Test: `tests/unit/benchsim/test_report.py`

**Interfaces:**
- Consumes:
  - Task 7 rows (from `scripts/evaluate.py` JSON, one report per engine);
  - Task 8 stats;
  - case `design` fields (joined by `design_id`);
  - `DurableLedger` (pre-registration entries live in `<out_root>/test/preregistration.jsonl`, append-only, with sha256 of the rule text).
- Produces:
  - `stratified_table(rows, by: Sequence[str], metric: str) -> list[dict]` (k, n, rate, CP interval per stratum)
  - `paired(rows_a, rows_b, metric) -> dict` (b, c, McNemar p)
  - `decide(reports: dict[str, list[dict]], baseline: str, candidate: str) -> dict`. For every profile it computes the per-allele exact sequence (superiority by McNemar, Holm across primary metrics), the normal FP (non-inferiority, margin 0.005), and the critical false negatives (candidate ≤ baseline). It returns `{"profiles": {...}, "adopt": bool}`.
  - `preregister(rule_text: str, path: Path) -> str` (sha256), and `require_preregistered(path: Path, rule_text: str) -> None`
  - `render_markdown(result: dict) -> str`

- [ ] **Step 1: Failing tests**

```python
# tests/unit/benchsim/test_report.py
from pathlib import Path

import pytest

from muc_one_span.benchsim.report import decide, paired, preregister, require_preregistered


def _rows(exact: list[int], fp: list[int], profile="ont_amplicon_r10"):
    return [{"sample": f"s{i}", "profile": profile, "allele_exact": e, "normal": True,
             "false_positive": f, "critical_false_negative": 0}
            for i, (e, f) in enumerate(zip(exact, fp))]


def test_paired_counts_discordant() -> None:
    a, b = _rows([1, 0, 0, 1], [0] * 4), _rows([1, 1, 1, 0], [0] * 4)
    res = paired(a, b, "allele_exact")
    assert (res["b"], res["c"]) == (1, 2)


def test_decide_requires_superiority_and_noninferiority() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    assert decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")["adopt"] is True
    worse = _rows([1] * n, [1] * 50 + [0] * (n - 50))
    assert decide({"ladder": base, "hybrid": worse}, "ladder", "hybrid")["adopt"] is False


def test_sealed_split_requires_preregistration(tmp_path: Path) -> None:
    path = tmp_path / "preregistration.jsonl"
    with pytest.raises(PermissionError):
        require_preregistered(path, "rule v1")
    preregister("rule v1", path)
    require_preregistered(path, "rule v1")
    with pytest.raises(PermissionError):
        require_preregistered(path, "rule v2")
```

Row fields are normalized from `scripts/evaluate.py` sample rows by a
`normalize_rows(report: dict, cases: dict[str, dict]) -> list[dict]` helper in this
module:
- `allele_exact` = all truth alleles sequence-exact;
- `normal` = `truth == "normal"`;
- `false_positive` = decision is PATHOGENIC on a normal or benign truth;
- `critical_false_negative` from Task 7.

Test it with one row from each clinical class.

- [ ] **Step 2–4:** implement. `scripts/benchsim.py evaluate --split test` calls `require_preregistered` before reading any test truth. `realism` runs Task 9 over a split and writes `realism.json` + `.md`. `report` writes `report.json` + `report.md` into `<results_root>`.
- [ ] **Step 5: Commit** — `git commit -m "feat(benchsim): stratified report, paired tests, decision rule, pre-registration"`

---

### Task 12: Documentation, changelog, first dev run

**Files:**
- Create: `docs/benchmark.md`
- Modify: `mkdocs.yml` (nav), `CHANGELOG.md` (Unreleased), `docs/development.md` (link)

- [ ] **Step 1:** Write `docs/benchmark.md` covering:
  - prerequisites: MucOneUp ≥ 0.45.0, pbsim3 (`pbsim`), ccs, and edlib;
  - the three profiles, including that HiFi is **uncalibrated**;
  - the data layout outside Git, the seeds and sealing, and the pre-registration workflow;
  - every `scripts/benchsim.py` subcommand with an example;
  - the realism report and the sim-to-real gap statement;
  - what may be committed: only public PRJEB92208 aggregates.
- [ ] **Step 2:** Run `make docs-check` and `make ci-check`.
- [ ] **Step 3: First dev run (outside Git), recorded in `.planning/`:** run `benchsim design --split dev --n 30`, then generate, realism, `run --engines ladder`, evaluate and report. Put the numbers, runtimes, realism pass/fail and failures into `.planning/2026-09-2x-benchsim-dev-pilot.md`. Keep no LB-level or patient data.
- [ ] **Step 4: Commit** — `git commit -m "docs(benchsim): benchmark guide, changelog and dev pilot"`

---

## Self-review notes

- **Spec coverage:**
  - §2 goals 1–4 → Tasks 2–6, 9.
  - §3 architecture → superseded as recorded above (MucOneUp read layer); code location and tool rules → Global Constraints.
  - §4 calibration → MucOneUp profiles, plus Task 3 variants and Task 9 realism with public targets. Stutter deconvolution is done upstream in MucOneUp (`helpers/calibrate_read_profile.py`), so it is not duplicated here.
  - §5 factors, splits, normals floor and sealing → Tasks 2, 6, 11.
  - §6 metrics and statistics → Tasks 7, 8, 11. Read-assignment accuracy needs per-read engine output, which the hybrid plan does not yet emit. It is deferred until the hybrid engine writes `read_assignments.tsv.gz`; the read truth needed for it is already stored (Task 5).
  - §7 outputs → `case.json`, manifest, `read_truth.tsv.gz` in `truth/`, evaluation extension.
  - §8 compute → `--jobs`, and `dev` is regenerable.
  - §9 risks → realism report, HiFi flagged, strict event validation (`design_invalid`).
- **Known gaps / follow-ups:**
  - A MucOneUp `reads profiles --show NAME` would remove the need for `--muconeup-profiles`. File it upstream if the importlib lookup proves fragile.
  - The in-house genomic realism check runs only locally.
- **Type consistency:** `Design` fields are used by Tasks 3, 4, 6, 10 and 11. `ReadTruth` is used by Tasks 5 and 9. `run_pipeline(..., engine=)` is shared with hybrid-plan Task 10.
