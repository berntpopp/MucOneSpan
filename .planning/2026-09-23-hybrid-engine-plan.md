# Hybrid Read-Centric Engine Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `muconespan run --engine hybrid`. It reconstructs MUC1 VNTR alleles from spanning reads (POA), assigns all reads to allele drafts wrapped in ladder flanks, polishes, and reports explicit read-level evidence. It keeps the existing output contract and leaves the ladder engine as the default.

**Architecture:** A new subpackage `muc_one_span/hybrid/` has one module per stage: anchor spans, length peaks, POA backend, polishing, phase split, read assignment, evidence, and the engine orchestrator. `pipeline.execute_pipeline` branches on `settings.run.engine` after writing the run configuration. Both engines then share one classification/summary/report tail. Clinical decisions and evaluation are changed so that support must be explicit (VCF or read-level).

**Tech Stack:** Python 3.10+, edlib (MIT), pyabpoa (MIT) with a pyspoa (MIT) fallback, numpy, pytest, Click, and the existing `classify_sequence`.

**Spec:** `.planning/2026-09-23-hybrid-engine-spec.md` (read it first). Prototype for reference only (outside Git): `../MucOneSpan-review-20260923/poa-prototype/{proto.py,hetsplit.py}`, `homopolymer/hp_model.py`, `inhouse/hybrid/assign_test.py`.

## Global Constraints

- Every authored file has **fewer than 650 physical lines** (max 649). Split by responsibility; do not compress lines.
- Keep CLI flags, output schemas, bundled resources, and scientific semantics of the ladder engine. The hybrid engine **adds fields only**.
- `--engine` default stays `ladder` until the benchmark decision rule (benchmark spec §6) passes on the sealed test split.
- New Python dependencies go in a new optional extra `hybrid` only. `pyproject.toml` and `uv.lock` change together (`make lock`).
- External commands only through `tools.py` (argument lists, no shell). The hybrid default path runs no external tools.
- Deterministic: every random choice uses `random.Random(settings.hybrid.seed)` (or a per-allele derivation of it).
- Unit tests are deterministic and synthetic, with no tool binaries. Mark integration tests `@pytest.mark.integration`.
- Never commit patient reads, sequences, or LB-level outputs. PRJEB92208 and in-house validation outputs stay outside Git.
- `make ci-check` before every commit; `make test-int`, `make docs-check`, `make build-check` where the task says so.

## Review Focus

1. **A sample whose reads carry a mutation inside motif 1 or motif 9.** Expect a flank-anchor fallback to find spans, recorded as `anchor_basis="flank"`, not zero spans and a silent no-call. Test lives in Task 3.
2. **An amplicon library with 40–50% short smear products.** Expect no phantom short allele: smear is reported as `short_product_fraction`, never as allele 2. Test lives in Task 4.
3. **Two alleles 44 bp apart** (in-frame 16 bp insertion on a Δ1 background). Expect either two alleles, or INCONCLUSIVE with `unassigned_spanning_fraction > 0.2`, never a confident homozygous call. Test lives in Task 4.
4. **A sample where the long allele has 5 spanning reads.** Expect a rejected peak to be listed and the clinical decision to be INCONCLUSIVE, not NO_PATHOGENIC. Tests live in Tasks 4 and 9.
5. **All reads on one strand** (MucOneUp artefact) versus mixed strands. Expect homopolymer evidence to report per-strand fractions and not to fail when one strand has zero reads. Test lives in Task 8.

---

## File Structure

| Path | Responsibility |
| --- | --- |
| `src/muc_one_span/hybrid/__init__.py` | public `reconstruct_alleles` re-export |
| `src/muc_one_span/hybrid/align.py` | edlib wrappers: `rc`, `cigar_ops`, `infix_hit`, `global_columns` (per-column read bases, insertions, target→query map) |
| `src/muc_one_span/hybrid/spans.py` | S1: anchors, orientation, read categories (`ReadRecord`, `SpanRead`, `categorize_reads`) |
| `src/muc_one_span/hybrid/lengths.py` | S2: length peaks with scaled windows, rejected peaks, short products (`LengthModel`) |
| `src/muc_one_span/hybrid/poa.py` | POA backend protocol with pyabpoa / pyspoa implementations |
| `src/muc_one_span/hybrid/polish.py` | S3/S7: draft consensus, pileup polish, homopolymer vote |
| `src/muc_one_span/hybrid/phase.py` | S4: candidate sites, linkage, split decision |
| `src/muc_one_span/hybrid/assign.py` | S5/S6: ladder-flanked references, edit-distance competition |
| `src/muc_one_span/hybrid/evidence.py` | S8/S10: residual QC, per-event read support, homopolymer LLR |
| `src/muc_one_span/hybrid/engine.py` | orchestration and output-contract dicts |
| `src/muc_one_span/pipeline_tail.py` | shared classify/summary/report tail (moved out of `pipeline.py`) |
| `src/muc_one_span/cli_run.py` | the `run` command, moved out of `cli.py` (cli.py is at 638/649 lines) |
| `src/muc_one_span/settings.py` | `HybridSettings`, `RunSettings.engine`, `RunSettings.assay` |
| `src/muc_one_span/report.py` | explicit-support clinical decision |
| `src/muc_one_span/evaluation/{artifacts,models}.py` | accept hybrid statuses and read support |
| `src/muc_one_span/benchmarking.py`, `clinical_runner.py` | `engine` forwarding |
| `tests/unit/hybrid/synth.py` | deterministic synthetic allele/read factory used by all hybrid tests |
| `tests/unit/hybrid/test_*.py` | one test module per stage |

---

### Task 1: Explicit-support clinical decision (safety prerequisite)

Both engines need this. Today a mutation with no `vcf_support` key counts as supported (`report.py:91-98`).

**Files:**
- Modify: `src/muc_one_span/report.py:85-100`
- Test: `tests/unit/test_clinical_decision.py`

**Interfaces:**
- Produces: `_mutation_supported(mut: dict) -> bool` in `report.py`; later tasks rely on `read_support.status == "supported"` being accepted.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/test_clinical_decision.py (append)
from muc_one_span.report import compute_clinical_decision


def _summary(mut: dict) -> dict:
    allele = {"length": 50, "canonical_repeats": 41, "reads": 200,
              "phase_status": "phased", "independent_haplotype_evidence": True}
    return {
        "alleles": {"allele_1": dict(allele), "allele_2": dict(allele, length=60, canonical_repeats=51),
                    "homozygous": False},
        "classifications": {
            "allele_1": {"mutations": [mut], "ambiguous_bases": 0,
                         "reconstruction_status": "complete_segmentation"},
            "allele_2": {"mutations": [], "ambiguous_bases": 0,
                         "reconstruction_status": "complete_segmentation"},
        },
    }


BASE = {"repeat_index": 20, "mutation_name": "dupC", "frameshift": True,
        "localization_status": "exact", "template_match": True}


def test_missing_support_fields_are_not_pathogenic() -> None:
    decision = compute_clinical_decision(_summary(dict(BASE)))
    assert decision["status"] != "PATHOGENIC"


def test_read_support_is_accepted_as_explicit_support() -> None:
    mut = dict(BASE, vcf_support=False, vcf_support_status="not_applicable_read_consensus",
               read_support={"status": "supported", "alt": 180, "ref": 10, "other": 10})
    assert compute_clinical_decision(_summary(mut))["status"] == "PATHOGENIC"


def test_vcf_exact_concordance_still_pathogenic() -> None:
    mut = dict(BASE, vcf_support=True, vcf_support_status="exact_sequence_concordance")
    assert compute_clinical_decision(_summary(mut))["status"] == "PATHOGENIC"
```

- [ ] **Step 2: Run to verify the first test fails**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py -k "missing_support or read_support or exact_concordance" --no-cov -v`
Expected: `test_missing_support_fields_are_not_pathogenic` FAILS (status PATHOGENIC); `test_read_support...` FAILS.
First check the actual key name returned by `compute_clinical_decision` and the name of its status values (e.g. `decision["status"]` or `decision["call"]`), and adjust the assertions to the existing key before continuing.

- [ ] **Step 3: Implement**

In `report.py`, replace the `supp_ok` expression with a helper:

```python
_SUPPORTED_VCF_STATUSES = {"exact_sequence_concordance"}


def _mutation_supported(mut: dict) -> bool:
    """Return True only for explicit sequence-level support (VCF or read evidence)."""
    read_support = mut.get("read_support")
    if isinstance(read_support, dict) and read_support.get("status") == "supported":
        return True
    if mut.get("vcf_support") is True:
        status = mut.get("vcf_support_status")
        # Legacy ladder artifacts without a status keep their historical meaning.
        return status is None or status in _SUPPORTED_VCF_STATUSES
    return False
```

and in the loop: `supp_ok = _mutation_supported(mut)`.

- [ ] **Step 4: Run the whole clinical and report suites**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py tests/unit/test_report.py tests/unit/test_report_wave1.py --no-cov -q`
Expected: PASS. If a legacy test relied on "missing support = pass", update it and state the reason in the commit message. This is an intended safety change (spec §4).

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/report.py tests/unit/test_clinical_decision.py
git commit -m "fix(report): require explicit VCF or read-level support for pathogenic decisions"
```

---

### Task 2: Settings, CLI engine/assay options, optional extra

**Files:**
- Modify: `src/muc_one_span/settings.py` (add `HybridSettings`, extend `RunSettings`, `RuntimeSettings`, `_SECTIONS`)
- Create: `src/muc_one_span/cli_run.py` (move the `run` command from `cli.py:437-543` unchanged, then add options)
- Modify: `src/muc_one_span/cli.py` (register `run` from `cli_run`)
- Modify: `pyproject.toml` (extra `hybrid`, mypy overrides), `Makefile` (`UV_TEST` adds `--extra hybrid`), `docker/Dockerfile` (`uv sync --extra report --extra hybrid`)
- Test: `tests/unit/test_runtime_settings.py`, `tests/unit/test_cli_run.py`

**Interfaces:**
- Produces: `RunSettings.engine: str` ("ladder"|"hybrid"); `RunSettings.assay: str` ("amplicon"|"genomic"); `RuntimeSettings.hybrid: HybridSettings` with the fields below; `execute_pipeline(..., engine=...)` does **not** change. The engine is read from `settings.run.engine`.

- [ ] **Step 1: Write failing settings tests**

```python
# tests/unit/test_runtime_settings.py (append)
import pytest

from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings, RunSettings, load_settings, settings_as_dict


def test_hybrid_defaults() -> None:
    h = DEFAULT_SETTINGS.hybrid
    assert (h.anchor_max_edits, h.n_poa, h.assign_margin, h.depth_adequate_spanning) == (12, 40, 3, 30)
    assert DEFAULT_SETTINGS.run.engine == "ladder"
    assert DEFAULT_SETTINGS.run.assay == "amplicon"


@pytest.mark.parametrize("field,value", [("n_poa", 0), ("assign_margin", -1), ("het_af_min", 0.9)])
def test_hybrid_rejects_invalid(field: str, value: object) -> None:
    with pytest.raises(ValueError):
        HybridSettings(**{field: value})


def test_engine_choice_validated() -> None:
    with pytest.raises(ValueError):
        RunSettings(engine="assembly")


def test_hybrid_roundtrip(tmp_path) -> None:
    data = settings_as_dict(DEFAULT_SETTINGS)
    data["hybrid"]["n_poa"] = 25
    data["run"]["engine"] = "hybrid"
    path = tmp_path / "s.json"
    import json
    path.write_text(json.dumps(data))
    loaded = load_settings(path)
    assert loaded.hybrid.n_poa == 25 and loaded.run.engine == "hybrid"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/test_runtime_settings.py --no-cov -q`
Expected: ImportError for `HybridSettings`.

- [ ] **Step 3: Implement the settings**

Add `engine: str = "ladder"` and `assay: str = "amplicon"` to `RunSettings`, with `_choice("engine", self.engine, ("ladder", "hybrid"))` and `_choice("assay", self.assay, ("amplicon", "genomic"))` in its `__post_init__`. Add:

```python
@dataclass(frozen=True)
class HybridSettings:
    """Read-centric engine thresholds (spec 2026-09-23 §5); tuned on dev/val splits only."""

    anchor_max_edits: int = 12
    min_span_units: int = 15
    max_span_units: int = 160
    peak_window_base_bp: float = 30.0
    peak_window_per_unit_bp: float = 0.6
    min_peak_reads: int = 8
    far_peak_min_frac: float = 0.03
    near_peak_min_frac: float = 0.20
    n_poa: int = 40
    n_poa_min: int = 10
    polish_rounds: int = 2
    hp_vote: bool = True
    het_af_min: float = 0.2
    het_min_group: float = 0.15
    link_phi_min: float = 0.5
    min_linked_sites: int = 2
    min_fragment_bp: int = 1000
    assign_margin: int = 3
    qc_residual_af: float = 0.25
    depth_adequate_spanning: int = 30
    depth_low_spanning: int = 10
    hp_llr_min: float = 10.0
    hp_min_reads: int = 20
    hp_min_alt_frac: float = 0.30
    seed: int = 1

    def __post_init__(self) -> None:
        for name in ("anchor_max_edits", "min_peak_reads", "n_poa_min", "polish_rounds",
                     "min_linked_sites", "min_fragment_bp", "assign_margin",
                     "depth_low_spanning", "hp_min_reads", "seed"):
            _integer(name, getattr(self, name))
        _integer("n_poa", self.n_poa, 1)
        _integer("min_span_units", self.min_span_units, 1)
        _integer("max_span_units", self.max_span_units, self.min_span_units + 1)
        _integer("depth_adequate_spanning", self.depth_adequate_spanning, self.depth_low_spanning)
        _number("peak_window_base_bp", self.peak_window_base_bp, 1)
        _number("peak_window_per_unit_bp", self.peak_window_per_unit_bp, 0)
        for name in ("far_peak_min_frac", "near_peak_min_frac", "het_min_group", "qc_residual_af",
                     "hp_min_alt_frac"):
            _number(name, getattr(self, name), 0, 1)
        _number("het_af_min", self.het_af_min, 0.01, 0.5)
        _number("link_phi_min", self.link_phi_min, 0, 1)
        _number("hp_llr_min", self.hp_llr_min, 0)
        _boolean("hp_vote", self.hp_vote)
```

Add `hybrid: HybridSettings = field(default_factory=HybridSettings)` to `RuntimeSettings` and `"hybrid": HybridSettings` to `_SECTIONS`.

- [ ] **Step 4: Run the settings tests**

Run: `uv run --locked --all-extras pytest tests/unit/test_runtime_settings.py tests/unit/test_cli_settings.py --no-cov -q`
Expected: PASS.

- [ ] **Step 5: Move `run` into `cli_run.py` and add the options**

Move the `run` command function and its decorators from `cli.py` into `cli_run.py`, keeping the same imports it uses. In `cli.py`, keep `from muc_one_span.cli_run import run` and `main.add_command(run)` (use the registration pattern already in `cli.py`; check how subcommands attach, `@main.command` vs `add_command`). Add two options to `run`:

```python
@click.option("--engine", type=click.Choice(["ladder", "hybrid"]), default=None,
              help="Allele reconstruction engine (default from settings: ladder).")
@click.option("--assay", type=click.Choice(["amplicon", "genomic"]), default=None,
              help="Library type; tunes hybrid-engine expectations (default: amplicon).")
```

Pass `engine=engine, assay=assay` into the existing `effective_run_settings(...)` call, following how `platform` is passed there. `configure_context` already maps `run.*` names to parameters.

- [ ] **Step 6: Write and run a CLI test**

```python
# tests/unit/test_cli_run.py (append)
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_span.cli import main


def test_engine_option_reaches_settings(tmp_path) -> None:
    fq = tmp_path / "r.fastq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, ["run", "-i", str(fq), "-o", str(tmp_path / "o"),
                                        "--engine", "hybrid", "--assay", "genomic"])
    assert res.exit_code == 0, res.output
    settings = ex.call_args.kwargs["settings"]
    assert settings.run.engine == "hybrid" and settings.run.assay == "genomic"
```

Run: `uv run --locked --all-extras pytest tests/unit/test_cli_run.py tests/unit/test_cli.py tests/unit/test_wave1_cli.py --no-cov -q`
Expected: PASS. If the patch target differs (run imports `execute_pipeline` lazily), patch where `cli_run` resolves it.

- [ ] **Step 7: Add the optional extra**

`pyproject.toml`:

```toml
[project.optional-dependencies]
hybrid = [
    "edlib>=1.3.9",
    "pyabpoa>=1.5.3",
    "pyspoa>=0.2.1",
    "numpy>=1.24",
]

[[tool.mypy.overrides]]
module = ["edlib", "pyabpoa", "spoa"]
ignore_missing_imports = true
```

Run `make lock`, review the `uv.lock` diff, then `make dev`. If `pyabpoa` fails to build from source on CI, keep `pyspoa` only in the extra and document pyabpoa as an optional speed-up (the backend is selected at runtime in Task 5). Update `Makefile` `UV_TEST` and `docker/Dockerfile` to include `--extra hybrid`.

- [ ] **Step 8: Full check and commit**

Run: `make ci-check`
Expected: PASS, with `cli.py` and `cli_run.py` both below 650 lines.

```bash
git add pyproject.toml uv.lock Makefile docker/Dockerfile src/muc_one_span/settings.py \
  src/muc_one_span/cli.py src/muc_one_span/cli_run.py tests/unit/test_runtime_settings.py tests/unit/test_cli_run.py
git commit -m "feat(settings): add hybrid engine settings, --engine/--assay options, hybrid extra"
```

---

### Task 3: Synthetic read factory, alignment helpers, and span categorisation (S1)

**Files:**
- Create: `tests/unit/hybrid/__init__.py`, `tests/unit/hybrid/synth.py`, `tests/unit/hybrid/test_spans.py`
- Create: `src/muc_one_span/hybrid/__init__.py`, `src/muc_one_span/hybrid/align.py`, `src/muc_one_span/hybrid/spans.py`

**Interfaces:**
- Produces:
  - `align.rc(seq: str) -> str`
  - `align.cigar_ops(cigar: str) -> list[tuple[int, str]]`
  - `align.infix_hit(query: str, target: str, k: int) -> tuple[int, int, int] | None` (start, end-exclusive, edits)
  - `align.global_columns(read: str, cons: str) -> Columns`, where `Columns` is a dataclass with `cols: list[str]`, `ins: dict[int, str]`, `t2q: list[int]`
  - `spans.ReadRecord(name: str, seq: str, qual: str)`
  - `spans.SpanRead(name, seq, mean_q, strand, anchor_edits, anchor_basis)` with property `length`
  - `spans.ReadCategories(spanning: list[SpanRead], left_anchored: list[ReadRecord], right_anchored: list[ReadRecord], internal_or_offtarget: list[ReadRecord])` with `counts() -> dict[str, int]`
  - `spans.categorize_reads(reads: Iterable[ReadRecord], anchors: Anchors, settings: HybridSettings) -> ReadCategories`
  - `spans.Anchors.from_dictionary(rd: RepeatDictionary) -> Anchors` with `left`, `right` (motif 1/9) and `left_flank`, `right_flank` (last/first 30 bp of the flanks)
  - `synth.allele(structure: list[str]) -> str`, `synth.mutate_unit(unit: str, name: str) -> str`, `synth.reads(allele_seq: str, n: int, *, err: float, seed: int, strand_mix: bool = True, flank_bp: int = 40, smear_frac: float = 0.0) -> list[ReadRecord]`

- [ ] **Step 1: Write the synthetic factory (test helper, deterministic)**

```python
# tests/unit/hybrid/synth.py
"""Deterministic synthetic MUC1 alleles and noisy long reads for hybrid-engine tests."""

from __future__ import annotations

import random

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.spans import ReadRecord

RD = load_repeat_dictionary()
PRE = ["1", "2", "3", "4", "5"]
POST = ["6", "7", "8", "9"]


def allele(inner: list[str]) -> str:
    """Motif 1..9 sequence for pre-repeats + inner units + after-repeats."""
    return "".join(RD.repeats[u] for u in PRE + inner + POST)


def dupc(unit: str = "X") -> str:
    """Unit with the dupC insertion (C after the 7C tract; dictionary insert at 60)."""
    return RD.repeats[unit] + "C"


def _noisy(seq: str, err: float, rng: random.Random) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < err / 3:
            continue  # deletion
        if r < 2 * err / 3:
            out.append(rng.choice("ACGT"))  # substitution
            continue
        out.append(base)
        if r < err:
            out.append(rng.choice("ACGT"))  # insertion
    return "".join(out)


def reads(allele_seq: str, n: int, *, err: float, seed: int, strand_mix: bool = True,
          flank_bp: int = 40, smear_frac: float = 0.0) -> list[ReadRecord]:
    rng = random.Random(seed)
    left = RD.flanking_left[-flank_bp:] if flank_bp else ""
    right = RD.flanking_right[:flank_bp] if flank_bp else ""
    out = []
    for i in range(n):
        template = left + allele_seq + right
        if smear_frac and rng.random() < smear_frac:
            cut_a = rng.randint(len(left) + 600, len(left) + len(allele_seq) // 2)
            cut_b = rng.randint(cut_a + 300, len(left) + len(allele_seq) - 600)
            template = template[:cut_a] + template[cut_b:]
        seq = _noisy(template, err, rng)
        if strand_mix and rng.random() < 0.5:
            seq = rc(seq)
        out.append(ReadRecord(f"r{seed}_{i}", seq, "5" * len(seq)))
    return out
```

- [ ] **Step 2: Write failing span tests**

```python
# tests/unit/hybrid/test_spans.py
from muc_one_span.hybrid.spans import Anchors, ReadRecord, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

ANCH = Anchors.from_dictionary(synth.RD)
S = HybridSettings()


def test_spanning_reads_oriented_and_trimmed_to_motif1_motif9() -> None:
    seq = synth.allele(["X"] * 30)
    cats = categorize_reads(synth.reads(seq, 20, err=0.02, seed=1), ANCH, S)
    assert len(cats.spanning) == 20
    for sp in cats.spanning:
        assert abs(sp.length - len(seq)) <= 40
        assert sp.seq.startswith(synth.RD.repeats["1"][:10]) or sp.anchor_edits > 0
    assert {sp.strand for sp in cats.spanning} == {"+", "-"}


def test_fragments_and_offtarget_are_not_spanning() -> None:
    seq = synth.allele(["X"] * 30)
    frag = ReadRecord("frag", seq[:900], "5" * 900)
    junk = ReadRecord("junk", "ACGT" * 300, "5" * 1200)
    cats = categorize_reads([frag, junk], ANCH, S)
    assert cats.counts() == {"spanning": 0, "left_anchored": 1, "right_anchored": 0,
                             "internal_or_offtarget": 1}


def test_mutated_motif1_falls_back_to_flank_anchor() -> None:
    seq = synth.allele(["X"] * 30)
    broken = seq[:5] + "TTTTTTTTTTTTTTTTTTTT" + seq[25:]  # destroy most of motif 1
    read = ReadRecord("m", synth.RD.flanking_left[-40:] + broken + synth.RD.flanking_right[:40], "5" * 2000)
    cats = categorize_reads([read], ANCH, HybridSettings(anchor_max_edits=6))
    assert len(cats.spanning) == 1 and cats.spanning[0].anchor_basis == "flank"
```

- [ ] **Step 3: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_spans.py --no-cov -q`
Expected: ModuleNotFoundError `muc_one_span.hybrid`.

- [ ] **Step 4: Implement `align.py`**

```python
# src/muc_one_span/hybrid/align.py
"""edlib-based alignment primitives shared by the hybrid engine stages."""

from __future__ import annotations

from dataclasses import dataclass

import edlib

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    """Reverse complement."""
    return seq.translate(_COMP)[::-1]


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    """Parse an extended CIGAR (=, X, I, D) into (length, op) pairs."""
    ops: list[tuple[int, str]] = []
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ""
    return ops


def infix_hit(query: str, target: str, k: int) -> tuple[int, int, int] | None:
    """Best infix location of query in target with at most k edits (start, end_excl, edits)."""
    res = edlib.align(query, target, mode="HW", task="locations", k=k)
    if res["editDistance"] < 0:
        return None
    start, end = res["locations"][0]
    return int(start), int(end) + 1, int(res["editDistance"])


def edit_distance_infix(query: str, target: str) -> int:
    """Edit distance of query aligned fully inside target (semi-global)."""
    return int(edlib.align(query, target, mode="HW", task="distance")["editDistance"])


@dataclass
class Columns:
    """A read projected onto consensus columns."""

    cols: list[str]  # read base or "-" per consensus column
    ins: dict[int, str]  # inserted read bases before consensus column i (i == len at end)
    t2q: list[int]  # consensus position -> read position (len(cons)+1 entries)


def global_columns(read: str, cons: str) -> Columns:
    """Globally align read to cons and project onto consensus columns."""
    res = edlib.align(read, cons, mode="NW", task="path")
    n = len(cons)
    cols = ["-"] * n
    ins: dict[int, str] = {}
    t2q = [0] * (n + 1)
    ti = qi = 0
    for length, op in cigar_ops(res["cigar"]):
        if op in "=X":
            for _ in range(length):
                cols[ti] = read[qi]
                t2q[ti] = qi
                ti += 1
                qi += 1
        elif op == "I":
            ins[ti] = ins.get(ti, "") + read[qi:qi + length]
            qi += length
        else:
            for _ in range(length):
                t2q[ti] = qi
                ti += 1
    t2q[n] = qi
    return Columns(cols, ins, t2q)
```

- [ ] **Step 5: Implement `spans.py`**

```python
# src/muc_one_span/hybrid/spans.py
"""S1: locate motif-1/motif-9 anchors, orient reads, and categorise them."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import infix_hit, rc
from muc_one_span.settings import HybridSettings

FLANK_ANCHOR_BP = 30


@dataclass(frozen=True)
class ReadRecord:
    name: str
    seq: str
    qual: str


@dataclass(frozen=True)
class SpanRead:
    name: str
    seq: str  # oriented motif1..motif9 (inclusive), forward strand of the VNTR
    mean_q: float
    strand: str  # "+" read was already forward, "-" read was reverse-complemented
    anchor_edits: int
    anchor_basis: str  # "motif" or "flank"

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass(frozen=True)
class Anchors:
    left: str
    right: str
    left_flank: str
    right_flank: str

    @classmethod
    def from_dictionary(cls, rd: RepeatDictionary) -> Anchors:
        return cls(rd.repeats["1"], rd.repeats["9"],
                   rd.flanking_left[-FLANK_ANCHOR_BP:], rd.flanking_right[:FLANK_ANCHOR_BP])


@dataclass
class ReadCategories:
    spanning: list[SpanRead] = field(default_factory=list)
    left_anchored: list[ReadRecord] = field(default_factory=list)
    right_anchored: list[ReadRecord] = field(default_factory=list)
    internal_or_offtarget: list[ReadRecord] = field(default_factory=list)

    def counts(self) -> dict[str, int]:
        return {"spanning": len(self.spanning), "left_anchored": len(self.left_anchored),
                "right_anchored": len(self.right_anchored),
                "internal_or_offtarget": len(self.internal_or_offtarget)}


def _mean_q(qual: str) -> float:
    return sum(ord(c) - 33 for c in qual) / len(qual) if qual else 0.0


def _span_in(target: str, anchors: Anchors, k: int) -> tuple[int, int, int, str] | None:
    """Return (start, end_excl, edits, basis) of motif1..motif9 in target, or None."""
    left, right = infix_hit(anchors.left, target, k), None
    basis = "motif"
    if left is not None:
        tail = infix_hit(anchors.right, target[left[1]:], k)
        right = None if tail is None else (tail[0] + left[1], tail[1] + left[1], tail[2])
    if left is None or right is None:
        kf = max(2, k // 4)
        lf = infix_hit(anchors.left_flank, target, kf)
        rf = infix_hit(anchors.right_flank, target[lf[1]:], kf) if lf else None
        if lf is None or rf is None:
            return None
        start, end = lf[1], rf[0] + lf[1]
        return start, end, lf[2] + rf[2], "flank"
    return left[0], right[1], left[2] + right[2], basis


def categorize_reads(reads: Iterable[ReadRecord], anchors: Anchors,
                     settings: HybridSettings) -> ReadCategories:
    """Classify reads as spanning (oriented, trimmed) / one-end anchored / other."""
    k = settings.anchor_max_edits
    lo, hi = settings.min_span_units * 60, settings.max_span_units * 60
    cats = ReadCategories()
    for rec in reads:
        seq = rec.seq.upper()
        best = None
        for strand, target, qual in (("+", seq, rec.qual), ("-", rc(seq), rec.qual[::-1])):
            hit = _span_in(target, anchors, k)
            if hit and (best is None or hit[2] < best[0][2]):
                best = (hit, strand, target, qual)
        if best is not None:
            (start, end, edits, basis), strand, target, qual = best
            if lo <= end - start <= hi:
                cats.spanning.append(SpanRead(rec.name, target[start:end],
                                              _mean_q(qual[start:end]), strand, edits, basis))
                continue
        has_left = any(infix_hit(anchors.left, t, k) for t in (seq, rc(seq)))
        has_right = any(infix_hit(anchors.right, t, k) for t in (seq, rc(seq)))
        if has_left and not has_right:
            cats.left_anchored.append(rec)
        elif has_right and not has_left:
            cats.right_anchored.append(rec)
        else:
            cats.internal_or_offtarget.append(rec)
    return cats
```

`src/muc_one_span/hybrid/__init__.py`: a one-line module docstring for now (Task 9 adds the export).

- [ ] **Step 6: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_spans.py --no-cov -v`
Expected: PASS. If the flank-fallback test fails because 20 T's still leave a motif-1 hit within 6 edits, lengthen the destroyed stretch to 40 bp, not the tolerance.

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/hybrid tests/unit/hybrid
git commit -m "feat(hybrid): anchor-based read categorisation with flank fallback"
```

---

### Task 4: Length model with scaled windows, rejected peaks, and short products (S2)

**Files:**
- Create: `src/muc_one_span/hybrid/lengths.py`
- Test: `tests/unit/hybrid/test_lengths.py`

**Interfaces:**
- Consumes: `SpanRead` (Task 3), `HybridSettings`.
- Produces:
  - `LengthPeak(center_bp: float, support: int, members: list[SpanRead])`
  - `LengthModel(peaks: list[LengthPeak], rejected: list[dict], short_products: list[SpanRead], unassigned: list[SpanRead])` with `unassigned_fraction` and `short_product_fraction` properties
  - `fit_length_model(spans: list[SpanRead], settings: HybridSettings) -> LengthModel`
  - `window_bp(length_bp: float, settings) -> float`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/hybrid/test_lengths.py
from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.spans import Anchors, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD)


def spans(inner_units: int, n: int, seed: int, **kw):
    seq = synth.allele(["X"] * inner_units)
    return categorize_reads(synth.reads(seq, n, err=0.02, seed=seed, **kw), ANCH, S).spanning


def test_two_distant_alleles_with_minor_long_allele() -> None:
    model = fit_length_model(spans(30, 200, 1) + spans(70, 12, 2), S)
    assert len(model.peaks) == 2
    assert sorted(round(p.center_bp / 60) for p in model.peaks) == [39, 79]


def test_smear_is_short_product_not_allele() -> None:
    model = fit_length_model(spans(60, 120, 3, smear_frac=0.45), S)
    assert len(model.peaks) == 1
    assert model.short_product_fraction > 0.3


def test_close_alleles_are_not_silently_dropped() -> None:
    # Alleles one unit apart must yield two peaks or a large unassigned fraction.
    model = fit_length_model(spans(40, 150, 4) + spans(41, 150, 5), S)
    assert len(model.peaks) == 2 or model.unassigned_fraction > 0.2


def test_rejected_minor_peak_is_reported() -> None:
    model = fit_length_model(spans(30, 200, 6) + spans(80, 4, 7), S)
    assert len(model.peaks) == 1
    assert model.rejected and model.rejected[0]["support"] == 4
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_lengths.py --no-cov -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/hybrid/lengths.py
"""S2: allele length peaks from spanning-read lengths with length-scaled windows."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

UNIT = 60


@dataclass
class LengthPeak:
    center_bp: float
    support: int
    members: list[SpanRead] = field(default_factory=list)


@dataclass
class LengthModel:
    peaks: list[LengthPeak]
    rejected: list[dict]
    short_products: list[SpanRead]
    unassigned: list[SpanRead]
    total: int

    @property
    def unassigned_fraction(self) -> float:
        return len(self.unassigned) / self.total if self.total else 0.0

    @property
    def short_product_fraction(self) -> float:
        return len(self.short_products) / self.total if self.total else 0.0


def window_bp(length_bp: float, settings: HybridSettings) -> float:
    """Assignment half-window: grows with allele length (ONT span noise ~0.25 bp/unit)."""
    return settings.peak_window_base_bp + settings.peak_window_per_unit_bp * length_bp / UNIT


def _density(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grid = np.arange(arr.min() - 100, arr.max() + 100, 2.0)
    dens = np.zeros_like(grid)
    for value in arr:
        bw = 8.0 + 0.004 * value
        dens += np.exp(-0.5 * ((grid - value) / bw) ** 2)
    return grid, dens


def fit_length_model(spans: list[SpanRead], settings: HybridSettings) -> LengthModel:
    """Pick up to two allele peaks; report rejected peaks, short products, unassigned reads."""
    if not spans:
        return LengthModel([], [], [], [], 0)
    arr = np.asarray([s.length for s in spans], dtype=float)
    grid, dens = _density(arr)
    maxima = [i for i in range(1, len(grid) - 1) if dens[i] >= dens[i - 1] and dens[i] > dens[i + 1]]
    maxima.sort(key=lambda i: -dens[i])
    centers: list[float] = []
    for i in maxima:  # keep maxima at least ~0.7 unit apart (Δ1 alleles stay separable)
        if all(abs(grid[i] - c) >= 0.7 * UNIT for c in centers):
            centers.append(float(grid[i]))
    support = {c: int(np.sum(np.abs(arr - c) <= window_bp(c, settings))) for c in centers}
    top = max(centers, key=lambda c: support[c])
    kept, rejected = [top], []
    for c in sorted(centers, key=lambda c: -support[c]):
        if c == top:
            continue
        far = abs(c - top) >= 2 * UNIT
        frac = settings.far_peak_min_frac if far else settings.near_peak_min_frac
        smear = c < min(kept) - 1.5 * UNIT and not far
        ok = (support[c] >= settings.min_peak_reads and support[c] >= frac * support[top]
              and not smear and len(kept) < 2)
        if ok:
            kept.append(c)
            continue
        if smear:
            reason = "smear"
        elif len(kept) < 2:
            reason = "support_below_threshold"
        else:
            reason = "max_alleles"
        rejected.append({"center_bp": round(c, 1), "units": round(c / UNIT),
                         "support": support[c], "reason": reason})
    peaks = [LengthPeak(c, support[c]) for c in sorted(kept)]
    short_products: list[SpanRead] = []
    unassigned: list[SpanRead] = []
    shortest = peaks[0].center_bp
    for sp in spans:
        dist = [abs(sp.length - p.center_bp) for p in peaks]
        j = int(np.argmin(dist))
        if dist[j] <= window_bp(peaks[j].center_bp, settings):
            peaks[j].members.append(sp)
        elif sp.length < shortest - 1.5 * UNIT:
            short_products.append(sp)
        else:
            unassigned.append(sp)
    return LengthModel(peaks, rejected, short_products, unassigned, len(spans))
```

- [ ] **Step 4: Run the tests and tune only the window formula if needed**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_lengths.py --no-cov -v`
Expected: PASS. If `test_close_alleles...` fails, do **not** loosen the assertion. Check that the Δ1 peaks (60 bp apart) survive the 0.7-unit separation and that `near_peak_min_frac` (0.20) passes with equal support.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/hybrid/lengths.py tests/unit/hybrid/test_lengths.py
git commit -m "feat(hybrid): length peaks with scaled windows, rejected-peak and smear reporting"
```

---

### Task 5: POA backend and polishing (S3, S7)

**Files:**
- Create: `src/muc_one_span/hybrid/poa.py`, `src/muc_one_span/hybrid/polish.py`
- Test: `tests/unit/hybrid/test_polish.py`

**Interfaces:**
- Produces:
  - `poa.PoaBackend` (Protocol) with `consensus(seqs: list[str]) -> str`
  - `poa.get_backend(name: str | None = None) -> PoaBackend`; tries pyabpoa, then pyspoa; raises `ImportError("... pip install muc_one_span[hybrid]")`
  - `polish.draft_consensus(members: list[SpanRead], n_poa: int, rng: random.Random, backend: PoaBackend) -> str`
  - `polish.pileup_polish(cons: str, reads: list[str]) -> tuple[str, int]`
  - `polish.homopolymer_vote(cons: str, reads: list[str], min_len: int = 4) -> tuple[str, int]`
  - `polish.polish(cons: str, reads: list[str], rounds: int, hp_vote: bool) -> tuple[str, dict]`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/hybrid/test_polish.py
import random

import pytest

from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import draft_consensus, homopolymer_vote, polish
from muc_one_span.hybrid.spans import Anchors, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib")


def _members(seq: str, n: int, seed: int):
    return categorize_reads(synth.reads(seq, n, err=0.03, seed=seed),
                            Anchors.from_dictionary(synth.RD), HybridSettings()).spanning


def test_poa_plus_polish_recovers_exact_allele_with_dupc() -> None:
    inner = ["X"] * 12 + ["A", "B"] + ["X"] * 10
    seq = "".join(synth.RD.repeats[u] for u in synth.PRE) + "".join(
        synth.RD.repeats[u] for u in inner[:5]) + synth.dupc("X") + "".join(
        synth.RD.repeats[u] for u in inner[5:] + synth.POST)
    members = _members(seq, 60, 11)
    cons = draft_consensus(members, 40, random.Random(1), get_backend())
    cons, info = polish(cons, [m.seq for m in members], rounds=2, hp_vote=True)
    assert cons == seq, info


def test_homopolymer_vote_uses_median_run_length() -> None:
    cons = "ACGT" + "C" * 7 + "ACGT"
    reads = ["ACGT" + "C" * 8 + "ACGT"] * 6 + ["ACGT" + "C" * 7 + "ACGT"] * 3
    new, changes = homopolymer_vote(cons, reads)
    assert new == "ACGT" + "C" * 8 + "ACGT" and changes == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_polish.py --no-cov -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement `poa.py`**

```python
# src/muc_one_span/hybrid/poa.py
"""Partial-order-alignment consensus backends (pyabpoa preferred, pyspoa fallback)."""

from __future__ import annotations

from typing import Protocol

_HINT = "Hybrid engine requires the 'hybrid' extra: pip install muc_one_span[hybrid]"


class PoaBackend(Protocol):
    name: str

    def consensus(self, seqs: list[str]) -> str: ...


class _Abpoa:
    name = "pyabpoa"

    def __init__(self) -> None:
        import pyabpoa

        self._aligner = pyabpoa.msa_aligner(aln_mode="g")

    def consensus(self, seqs: list[str]) -> str:
        res = self._aligner.msa(seqs, out_cons=True, out_msa=False)
        return str(res.cons_seq[0])


class _Spoa:
    name = "pyspoa"

    def __init__(self) -> None:
        import spoa

        self._poa = spoa.poa

    def consensus(self, seqs: list[str]) -> str:
        cons, _msa = self._poa(seqs, algorithm=1)  # global alignment
        return str(cons)


def get_backend(name: str | None = None) -> PoaBackend:
    """Return the requested backend or the first importable one."""
    order = [name] if name else ["pyabpoa", "pyspoa"]
    errors = []
    for candidate in order:
        try:
            return _Abpoa() if candidate == "pyabpoa" else _Spoa()
        except ImportError as exc:  # pragma: no cover - depends on installed extras
            errors.append(f"{candidate}: {exc}")
    raise ImportError(f"{_HINT} ({'; '.join(errors)})")
```

- [ ] **Step 4: Implement `polish.py`**

```python
# src/muc_one_span/hybrid/polish.py
"""S3/S7: POA draft, majority pileup polishing and homopolymer median vote."""

from __future__ import annotations

import random
import statistics
from collections import Counter

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.poa import PoaBackend
from muc_one_span.hybrid.spans import SpanRead


def draft_consensus(members: list[SpanRead], n_poa: int, rng: random.Random,
                    backend: PoaBackend) -> str:
    """POA over a random sample of near-modal members (random, not quality-ranked)."""
    median = statistics.median(m.length for m in members)
    tol = max(15.0, 0.006 * median)
    near = [m for m in members if abs(m.length - median) <= tol] or list(members)
    sample = rng.sample(near, min(n_poa, len(near)))
    return backend.consensus([m.seq for m in sample])


def pileup_polish(cons: str, reads: list[str]) -> tuple[str, int]:
    """One majority-vote round over columns and insertion slots."""
    n = len(cons)
    col = [Counter() for _ in range(n)]
    ins = [Counter() for _ in range(n + 1)]
    for read in reads:
        proj = global_columns(read, cons)
        for pos, base in enumerate(proj.cols):
            col[pos][base] += 1
        for pos in range(n + 1):
            ins[pos][proj.ins.get(pos, "")] += 1
    out: list[str] = []
    changes = 0
    for pos in range(n + 1):
        ins_vote, ins_count = ins[pos].most_common(1)[0]
        if ins_vote and ins_count > len(reads) / 2:
            out.append(ins_vote)
            changes += 1
        if pos < n:
            base = col[pos].most_common(1)[0][0]
            if base != "-":
                out.append(base)
            if base != cons[pos]:
                changes += 1
    return "".join(out), changes


def _runs(seq: str, min_len: int) -> list[tuple[int, int, str]]:
    runs, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            runs.append((i, j, seq[i]))
        i = j
    return runs


def read_run_length(read: str, t2q: list[int], start: int, end: int, base: str) -> int:
    """Longest stretch of `base` in the read around the interval mapped to [start, end)."""
    qs, qe = t2q[start], t2q[end]
    a = qs
    while a > 0 and read[a - 1] == base:
        a -= 1
    b = max(qe, qs)
    while b < len(read) and read[b] == base:
        b += 1
    seg = "".join(c if c == base else " " for c in read[a:b])
    return max((len(m) for m in seg.split()), default=0)


def homopolymer_vote(cons: str, reads: list[str], min_len: int = 4) -> tuple[str, int]:
    """Set each run (>= min_len) to the median read run length (ONT deletion bias)."""
    runs = _runs(cons, min_len)
    if not runs or not reads:
        return cons, 0
    maps = [(r, global_columns(r, cons).t2q) for r in reads]
    out, pos, changes = [], 0, 0
    for start, end, base in runs:
        lengths = [read_run_length(r, t2q, start, end, base) for r, t2q in maps]
        new_len = max(1, int(statistics.median(lengths) + 0.5))
        out.append(cons[pos:start])
        out.append(base * new_len)
        changes += int(new_len != end - start)
        pos = end
    out.append(cons[pos:])
    return "".join(out), changes


def polish(cons: str, reads: list[str], rounds: int, hp_vote: bool) -> tuple[str, dict]:
    """Run `rounds` pileup rounds, each followed by an optional homopolymer vote."""
    info: dict = {"rounds": []}
    for _ in range(rounds):
        cons, changes = pileup_polish(cons, reads)
        hp = 0
        if hp_vote:
            cons, hp = homopolymer_vote(cons, reads)
        info["rounds"].append({"changes": changes, "hp_changes": hp})
    return cons, info
```

- [ ] **Step 5: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_polish.py --no-cov -v`
Expected: PASS. At 3% synthetic error and 60 reads the prototype was exact. If the dupC test fails, inspect `info` and the diff before changing any setting.

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/hybrid/poa.py src/muc_one_span/hybrid/polish.py tests/unit/hybrid/test_polish.py
git commit -m "feat(hybrid): POA backends and pileup/homopolymer polishing"
```

---

### Task 6: Linked-site phase split for equal and close lengths (S4)

**Files:**
- Create: `src/muc_one_span/hybrid/phase.py`
- Test: `tests/unit/hybrid/test_phase.py`

**Interfaces:**
- Consumes: `global_columns`, `read_run_length`, `SpanRead`, `HybridSettings`.
- Produces:
  - `PhaseResult(groups: list[list[SpanRead]], basis: str, sites: list[dict], candidate: dict | None)`, where basis ∈ {"none", "linked_sites", "unconfirmed_single_site"}
  - `split_by_linked_sites(cons: str, members: list[SpanRead], settings: HybridSettings, rng: random.Random) -> PhaseResult`

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/hybrid/test_phase.py
import random

from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.spans import Anchors, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD)


def _spans(inner, n, seed):
    return categorize_reads(synth.reads(synth.allele(inner), n, err=0.02, seed=seed), ANCH, S).spanning


def test_equal_length_alleles_with_two_linked_differences_split() -> None:
    a = ["X"] * 10 + ["A"] + ["X"] * 10 + ["B"] + ["X"] * 8
    b = ["X"] * 10 + ["X"] + ["X"] * 10 + ["X"] + ["X"] * 8
    members = _spans(a, 40, 1) + _spans(b, 40, 2)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "linked_sites" and sorted(len(g) for g in res.groups) == [40, 40]


def test_identical_alleles_do_not_split() -> None:
    a = ["X"] * 30
    members = _spans(a, 80, 3)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "none" and len(res.groups) == 1
```

(A single-site test, where alleles differ by one unit type, asserts `basis == "unconfirmed_single_site"` and `len(res.groups) == 1`. Build it from `a` above with only the `A` difference.)

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_phase.py --no-cov -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement**

Port `site_table`, `candidate_sites`, `_indicator` and `linked_sites` from the prototype `hetsplit.py` (shown in `.planning/` review notes; the source is in `../MucOneSpan-review-20260923/poa-prototype/hetsplit.py:87-260`), with these adaptations:

```python
# src/muc_one_span/hybrid/phase.py
"""S4: split a length peak only when >=2 read-linked difference sites support it."""

from __future__ import annotations

import math
import random
import statistics
from collections import Counter
from dataclasses import dataclass, field

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

MAX_SITE_READS = 300


@dataclass
class PhaseResult:
    groups: list[list[SpanRead]]
    basis: str
    sites: list[dict] = field(default_factory=list)
    candidate: dict | None = None


def _features(cons: str, reads: list[str]) -> tuple[list[dict], dict]:
    runs = _runs(cons, 3)
    in_run = {i for s, e, _ in runs for i in range(s, e)}
    feats = []
    for read in reads:
        proj = global_columns(read, cons)
        f: dict = {}
        for pos, base in enumerate(proj.cols):
            if pos not in in_run:
                f[("col", pos)] = base
                if pos - 1 not in in_run:
                    f[("ins", pos)] = proj.ins.get(pos, "")
        for s, e, b in runs:
            f[("run", s)] = read_run_length(read, proj.t2q, s, e, b)
        feats.append(f)
    return feats, {("run", s): (b, e - s) for s, e, b in runs}


def _candidates(feats: list[dict], meta: dict, af_min: float) -> list[dict]:
    counts: dict = {}
    for f in feats:
        for site, allele in f.items():
            counts.setdefault(site, Counter())[allele] += 1
    fracs: dict = {}
    for site, (base, length) in meta.items():
        c = counts.get(site, Counter())
        tot = sum(c.values()) or 1
        for x in set(range(max(0, length - 3), length + 4)) - {length}:
            fracs.setdefault((base, length, x), []).append(c.get(x, 0) / tot)
    bg = {k: statistics.median(v) for k, v in fracs.items()}
    out = []
    for site, c in counts.items():
        tot = sum(c.values())
        major, _ = c.most_common(1)[0]
        minors = [(a, n) for a, n in c.most_common(4) if a != major and n >= 5]
        if not minors:
            continue
        if site[0] == "run":
            base, length = meta[site]
            scored = [(n / tot - 4 * bg.get((base, length, a), 0.0), a, n) for a, n in minors]
            score, minor, n_minor = max(scored)
            if n_minor / tot < max(0.1, 4 * bg.get((base, length, minor), 0.0)):
                continue
        else:
            minor, n_minor = minors[0]
            af = n_minor / tot
            if af < af_min or ("-" in (major, minor) and af < 1.5 * af_min):
                continue
        out.append({"site": site, "major": major, "minor": minor,
                    "af": round(n_minor / tot, 3), "n": tot})
    return out


def _linked(feats: list[dict], sites: list[dict], min_phi: float) -> list[dict]:
    def ind(f: dict, s: dict) -> int | None:
        a = f.get(s["site"])
        return 1 if a == s["minor"] else 0 if a == s["major"] else None

    vec = [[ind(f, s) for s in sites] for f in feats]
    n = len(sites)
    adj: dict[int, set[int]] = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            pairs = [(v[i], v[j]) for v in vec if v[i] is not None and v[j] is not None]
            if len(pairs) < 10:
                continue
            a = sum(1 for x, y in pairs if x and y)
            b = sum(1 for x, y in pairs if x and not y)
            c = sum(1 for x, y in pairs if not x and y)
            d = len(pairs) - a - b - c
            den = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
            if den and abs(a * d - b * c) / den >= min_phi:
                adj[i].add(j)
                adj[j].add(i)
    seen: set[int] = set()
    best: list[int] = []
    for i in range(n):
        if i in seen:
            continue
        comp, stack = [], [i]
        seen.add(i)
        while stack:
            k = stack.pop()
            comp.append(k)
            for m in adj[k] - seen:
                seen.add(m)
                stack.append(m)
        best = comp if len(comp) > len(best) else best
    return [sites[i] for i in sorted(best)]


def split_by_linked_sites(cons: str, members: list[SpanRead], settings: HybridSettings,
                          rng: random.Random) -> PhaseResult:
    """Return one group (no split) unless >= min_linked_sites linked sites support two groups."""
    sample = members if len(members) <= MAX_SITE_READS else rng.sample(members, MAX_SITE_READS)
    feats, meta = _features(cons, [m.seq for m in sample])
    sites = _candidates(feats, meta, settings.het_af_min)
    if not sites:
        return PhaseResult([members], "none")
    linked = _linked(feats, sites, settings.link_phi_min) if len(sites) > 1 else sites[:1]
    if len(linked) < settings.min_linked_sites:
        top = max(sites, key=lambda s: s["af"] * (1 - s["af"]))
        return PhaseResult([members], "unconfirmed_single_site", sites, candidate=top)
    all_feats, _ = _features(cons, [m.seq for m in members])
    groups: list[list[SpanRead]] = [[], []]
    for m, f in zip(members, all_feats):
        votes = [1 if f.get(s["site"]) == s["minor"] else 0 if f.get(s["site"]) == s["major"] else None
                 for s in linked]
        known = [v for v in votes if v is not None]
        if known:
            groups[int(sum(known) * 2 >= len(known) + 1)].append(m)
    if min(len(g) for g in groups) < settings.het_min_group * len(members):
        return PhaseResult([members], "unconfirmed_single_site", linked, candidate=linked[0])
    return PhaseResult(groups, "linked_sites", linked)
```

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_phase.py --no-cov -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/hybrid/phase.py tests/unit/hybrid/test_phase.py
git commit -m "feat(hybrid): linked-site phasing for equal and close allele lengths"
```

---

### Task 7: All-read assignment to ladder-flanked drafts (S5, S6)

**Files:**
- Create: `src/muc_one_span/hybrid/assign.py`
- Test: `tests/unit/hybrid/test_assign.py`

**Interfaces:**
- Consumes: `ReadRecord`, `SpanRead`, `rc`, `edit_distance_infix`, `RepeatDictionary`.
- Produces:
  - `hybrid_references(drafts: dict[str, str], rd: RepeatDictionary, flank_bp: int = 500) -> dict[str, str]`
  - `Assignment(allele: str | None, margin: int, distances: dict[str, int])`
  - `assign_read(seq: str, refs: dict[str, str], margin: int) -> Assignment`, which tries both orientations and picks the orientation with the smaller best distance
  - `assign_reads(reads: list[ReadRecord], refs: dict[str, str], settings: HybridSettings) -> dict[str, list[str]]` (allele → oriented read sequences, plus key `"undecided"`); only reads ≥ `min_fragment_bp` are considered

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/hybrid/test_assign.py
import random

from muc_one_span.hybrid.assign import assign_read, hybrid_references
from muc_one_span.hybrid.align import rc
from tests.unit.hybrid import synth

A = ["X"] * 10 + ["A", "A", "B"] + ["X"] * 20 + ["G", "A", "B"] + ["X"] * 10
B = ["X"] * 10 + ["A", "A", "B"] + ["X"] * 21 + ["G", "A", "B"] + ["X"] * 10  # one unit longer


def test_internal_fragments_assigned_to_correct_near_identical_allele() -> None:
    refs = hybrid_references({"allele_1": synth.allele(A), "allele_2": synth.allele(B)}, synth.RD)
    rng = random.Random(3)
    wrong = decided = 0
    for truth, inner in (("allele_1", A), ("allele_2", B)):
        seq = synth.allele(inner)
        for i in range(30):
            start = rng.randint(0, len(seq) - 1500)
            frag = synth.reads(seq[start:start + 1500], 1, err=0.02, seed=100 + i, flank_bp=0)[0].seq
            got = assign_read(rc(frag) if i % 2 else frag, refs, margin=3)
            if got.allele is not None:
                decided += 1
                wrong += got.allele != truth
    assert wrong == 0 and decided >= 20


def test_identical_region_is_undecided() -> None:
    refs = hybrid_references({"allele_1": synth.allele(A), "allele_2": synth.allele(B)}, synth.RD)
    frag = synth.RD.repeats["X"] * 5
    assert assign_read(frag, refs, margin=3).allele is None
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_assign.py --no-cov -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/hybrid/assign.py
"""S5/S6: ladder-flanked allele references and edit-distance competition for every read."""

from __future__ import annotations

from dataclasses import dataclass

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import edit_distance_infix, rc
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.settings import HybridSettings


def hybrid_references(drafts: dict[str, str], rd: RepeatDictionary,
                      flank_bp: int = 500) -> dict[str, str]:
    """Wrap each motif1..motif9 draft in the ladder's hg38 flanks."""
    left, right = rd.flanking_left[-flank_bp:], rd.flanking_right[:flank_bp]
    return {name: left + seq + right for name, seq in drafts.items()}


@dataclass(frozen=True)
class Assignment:
    allele: str | None
    margin: int
    distances: dict[str, int]
    oriented: str


def assign_read(seq: str, refs: dict[str, str], margin: int) -> Assignment:
    """Assign seq to the reference with the smallest infix edit distance, if clearly best."""
    best: Assignment | None = None
    for oriented in (seq, rc(seq)):
        dist = {name: edit_distance_infix(oriented, ref) for name, ref in refs.items()}
        ranked = sorted(dist.values())
        gap = ranked[1] - ranked[0] if len(ranked) > 1 else 10**9
        winner = min(dist, key=lambda k: dist[k])
        cand = Assignment(winner if gap >= margin else None, gap, dist, oriented)
        if best is None or ranked[0] < min(best.distances.values()):
            best = cand
    assert best is not None
    return best


def assign_reads(reads: list[ReadRecord], refs: dict[str, str],
                 settings: HybridSettings) -> dict[str, list[str]]:
    """Map allele -> oriented read sequences; ambiguous reads go to 'undecided'."""
    out: dict[str, list[str]] = {name: [] for name in refs}
    out["undecided"] = []
    for rec in reads:
        if len(rec.seq) < settings.min_fragment_bp:
            continue
        got = assign_read(rec.seq.upper(), refs, settings.assign_margin)
        out[got.allele or "undecided"].append(got.oriented)
    return out
```

For a single-allele (homozygous) sample `assign_read` sees one reference, so every read of at least `min_fragment_bp` is assigned. That is the intended behaviour.

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_assign.py --no-cov -v`
Expected: PASS (the in-house experiment gave 0 wrong out of 318 at 1.5–2.5 kb for 82/83-unit alleles).

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/hybrid/assign.py tests/unit/hybrid/test_assign.py
git commit -m "feat(hybrid): assign all reads to ladder-flanked allele drafts by edit-distance margin"
```

---

### Task 8: Residual QC and per-event read evidence, including homopolymer LLR (S8, S10)

**Files:**
- Create: `src/muc_one_span/hybrid/evidence.py`
- Test: `tests/unit/hybrid/test_evidence.py`

**Interfaces:**
- Consumes: `global_columns`, `read_run_length`, `_runs`, `HybridSettings`, and the classification dict from `classify_sequence` (mutations carry `repeat_index`, `mutation_name`, `differences`, …; each repeat has `index`, `start`, `end`).
- Produces:
  - `residual_sites(cons: str, reads: list[str], af: float) -> list[dict]`
  - `homopolymer_background(cons: str, reads: list[tuple[str, str]], base: str, length: int) -> dict[str, list[float]]` (per strand "+", "-", and "all": smoothed P(observed length) over 0..16)
  - `event_read_support(cons: str, classification: dict, reads: list[tuple[str, str]], settings: HybridSettings) -> dict[int, dict]` keyed by the mutation's position in `classification["mutations_detected"]`. Reads are `(sequence, strand)` pairs, already oriented to the consensus.

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/hybrid/test_evidence.py
from muc_one_span.hybrid.evidence import homopolymer_llr, residual_sites


def test_residual_sites_flag_mixture_but_not_clean_reads() -> None:
    cons = "ACGTACGTTTGACCATGCA" * 10
    alt = cons[:50] + "G" + cons[51:]
    assert residual_sites(cons, [cons] * 20, af=0.25) == []
    sites = residual_sites(cons, [cons] * 12 + [alt] * 8, af=0.25)
    assert [s["pos"] for s in sites] == [50]


def test_homopolymer_llr_separates_8c_from_7c_background() -> None:
    background = {"+": [0.02] * 6 + [0.20, 0.55, 0.20] + [0.03 / 8] * 8,
                  "-": [0.01] * 6 + [0.05, 0.88, 0.05] + [0.01 / 8] * 8}
    obs_mut = [("+", 8)] * 12 + [("+", 7)] * 8 + [("-", 8)] * 18 + [("-", 7)] * 2
    obs_wt = [("+", 7)] * 12 + [("+", 6)] * 8 + [("-", 7)] * 18 + [("-", 8)] * 2
    assert homopolymer_llr(obs_mut, background, 7) > 10
    assert homopolymer_llr(obs_wt, background, 7) < 0


def test_llr_tolerates_single_strand_input() -> None:
    background = {"+": [0.05] * 17, "-": [0.05] * 17}
    assert isinstance(homopolymer_llr([("-", 8)] * 25, background, 7), float)
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_evidence.py --no-cov -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/hybrid/evidence.py
"""S8/S10: residual heterogeneity QC and per-event read-level support."""

from __future__ import annotations

import math
from collections import Counter

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.settings import HybridSettings

MAXLEN = 16


def residual_sites(cons: str, reads: list[str], af: float) -> list[dict]:
    """Consensus columns where a non-consensus base/deletion reaches `af` (outside runs >= 3)."""
    in_run = {i for s, e, _ in _runs(cons, 3) for i in range(s, e)}
    counts = [Counter() for _ in cons]
    for read in reads:
        for pos, base in enumerate(global_columns(read, cons).cols):
            counts[pos][base] += 1
    out = []
    for pos, c in enumerate(counts):
        if pos in in_run or not c:
            continue
        tot = sum(c.values())
        for base, n in c.most_common():
            if base != cons[pos] and n / tot >= af:
                out.append({"pos": pos, "ref": cons[pos], "alt": base, "af": round(n / tot, 3), "n": tot})
                break
    return out


def _smooth(counter: Counter) -> list[float]:
    tot = sum(counter.values())
    return [(counter.get(i, 0) + 0.5) / (tot + 0.5 * (MAXLEN + 1)) for i in range(MAXLEN + 1)]


def _shift(p: list[float], d: int) -> list[float]:
    q = [p[min(max(i - d, 0), MAXLEN)] for i in range(MAXLEN + 1)]
    s = sum(q)
    return [x / s for x in q]


def homopolymer_background(cons: str, reads: list[tuple[str, str]], base: str,
                           length: int, exclude: set[int] | None = None) -> dict[str, list[float]]:
    """Per-strand observed-length profile at all consensus runs of `base` x `length`."""
    runs = [(s, e) for s, e, b in _runs(cons, length) if b == base and e - s == length
            and s not in (exclude or set())]
    per: dict[str, Counter] = {"+": Counter(), "-": Counter()}
    for seq, strand in reads:
        t2q = global_columns(seq, cons).t2q
        for s, e in runs:
            per[strand][min(read_run_length(seq, t2q, s, e, base), MAXLEN)] += 1
    return {s: _smooth(c) for s, c in per.items()}


def homopolymer_llr(obs: list[tuple[str, int]], background: dict[str, list[float]], n0: int) -> float:
    """log L(true n0+1) - log L(true n0), strand-specific; the n0+1 profile is the shifted n0 profile."""
    llr = 0.0
    for strand, observed in obs:
        p0 = background.get(strand) or background.get("+") or [1 / (MAXLEN + 1)] * (MAXLEN + 1)
        p1 = _shift(p0, 1)
        k = min(observed, MAXLEN)
        llr += math.log(p1[k] / p0[k])
    return llr


def event_read_support(cons: str, classification: dict, reads: list[tuple[str, str]],
                       settings: HybridSettings) -> dict[int, dict]:
    """Read-level support for each detected mutation on the reads assigned to this allele."""
    results: dict[int, dict] = {}
    repeats = {r["index"]: r for r in classification.get("repeats", [])}
    projections = [(seq, strand, global_columns(seq, cons)) for seq, strand in reads]
    for idx, mut in enumerate(classification.get("mutations_detected", [])):
        rep = repeats.get(mut.get("repeat_index"))
        if rep is None:
            results[idx] = {"status": "not_localized"}
            continue
        start, end = rep["start"], rep["end"]
        runs = [(s, e, b) for s, e, b in _runs(cons[start:end], 6)]
        if runs:  # homopolymer event: use stutter-aware LLR on the longest run in the unit
            s, e, b = max(runs, key=lambda r: r[1] - r[0])
            s, e = s + start, e + start
            n_event = e - s
            bg = homopolymer_background(cons, reads, b, n_event - 1, exclude={s})
            obs = [(strand, read_run_length(seq, p.t2q, s, e, b)) for seq, strand, p in projections]
            alt = sum(1 for _, k in obs if k >= n_event)
            llr = homopolymer_llr(obs, bg, n_event - 1)
            by_strand = {st: [k for sd, k in obs if sd == st] for st in "+-"}
            strand_alt = {st: (round(sum(k >= n_event for k in v) / len(v), 3) if v else None)
                          for st, v in by_strand.items()}
            n = len(obs)
            status = ("insufficient_depth" if n < settings.hp_min_reads else
                      "supported" if llr >= settings.hp_llr_min and alt / n >= settings.hp_min_alt_frac
                      else "not_supported")
            results[idx] = {"kind": "homopolymer", "n": n, "alt": alt, "alt_frac": round(alt / n, 3) if n else 0.0,
                            "strand_alt_frac": strand_alt, "llr": round(llr, 1), "status": status}
            continue
        window = cons[start:end]
        alt = ref = other = 0
        for seq, _strand, proj in projections:
            q0, q1 = proj.t2q[start], proj.t2q[end]
            piece = seq[q0:q1]
            if piece == window:
                alt += 1
            elif abs(len(piece) - len(window)) >= 1 and len(piece) % 60 == 0:
                ref += 1
            else:
                other += 1
        n = alt + ref + other
        frac = alt / n if n else 0.0
        status = ("insufficient_depth" if n < settings.hp_min_reads else
                  "supported" if frac >= settings.hp_min_alt_frac and alt > ref else "not_supported")
        results[idx] = {"kind": "exact_window", "n": n, "alt": alt, "ref": ref, "other": other,
                        "alt_frac": round(frac, 3), "status": status}
    return results
```

For non-homopolymer events the exact-window rule is deliberately conservative: an exact unit match at 2% error happens in only ~30% of reads. Its thresholds are tuned on the benchmark `dev` split (benchmark plan Task 9). Until then, those events can reach at most INCONCLUSIVE unless the VCF-free support passes.

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_evidence.py --no-cov -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/hybrid/evidence.py tests/unit/hybrid/test_evidence.py
git commit -m "feat(hybrid): residual QC and stutter-aware per-event read support"
```

---

### Task 9: Engine orchestration, output contract, and pipeline branch

**Files:**
- Create: `src/muc_one_span/hybrid/engine.py`, `src/muc_one_span/pipeline_tail.py`
- Modify: `src/muc_one_span/pipeline.py` (move lines 142–217 into `pipeline_tail.finish_run`; add the engine branch after line 89; make `check_tools` engine-dependent), `src/muc_one_span/hybrid/__init__.py`
- Test: `tests/unit/hybrid/test_engine.py`, and the existing `tests/unit/test_wave1_cli.py` and `tests/unit/test_cli_run.py` must stay green

**Interfaces:**
- Consumes: everything from Tasks 3–8, `classify_sequence`, `load_repeat_dictionary`.
- Produces:
  - `hybrid.engine.reconstruct_alleles(input_path: Path, output_dir: Path, rd: RepeatDictionary, settings: RuntimeSettings) -> tuple[dict, dict[str, Path]]`, returning the alleles dict (same top-level keys as `detect_alleles`: `allele_1`, `allele_2`, `homozygous`, `same_length`, `observed_length_candidates`, `allele_multiplicity_status`) plus `hybrid` and `consensus_paths`
  - `hybrid.engine.read_input(path: Path) -> list[ReadRecord]` (FASTQ, FASTQ.gz, or BAM via `samtools fastq -F 0x900` through `tools.run_tool`)
  - `pipeline_tail.finish_run(alleles_result, consensus_paths, vcf_paths, rd, settings, output_dir, ...) -> dict` (summary). It is exactly the moved code, parameterised.

- [ ] **Step 1: Write the failing end-to-end unit test (synthetic FASTQ, no external tools)**

```python
# tests/unit/hybrid/test_engine.py
import json
from pathlib import Path

from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth


def _fastq(path: Path, reads) -> None:
    path.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in reads))


def test_heterozygous_dupc_sample_reconstructed_with_read_support(tmp_path: Path) -> None:
    a_inner = ["X"] * 25
    b_units = ["X"] * 14 + ["DUPC"] + ["X"] * 30
    a = synth.allele(a_inner)
    b = "".join(synth.RD.repeats[u] for u in synth.PRE) + "".join(
        synth.dupc() if u == "DUPC" else synth.RD.repeats[u] for u in b_units) + "".join(
        synth.RD.repeats[u] for u in synth.POST)
    reads = synth.reads(a, 150, err=0.02, seed=1) + synth.reads(b, 60, err=0.02, seed=2)
    fq = tmp_path / "in.fastq"
    _fastq(fq, reads)
    alleles, paths = reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)
    seqs = {Path(p).read_text().split("\n", 1)[1].replace("\n", "") for p in paths.values()}
    assert seqs == {a, b}
    for key in ("allele_1", "allele_2"):
        info = alleles[key]
        assert info["engine"] == "hybrid" and info["depth_status"] == "adequate"
        assert info["length"] == info["canonical_repeats"] + 9
    assert json.loads((tmp_path / "hybrid_reads.json").read_text())["read_categories"]["spanning"] == 210
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_engine.py --no-cov -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement `engine.py`**

```python
# src/muc_one_span/hybrid/engine.py
"""Hybrid engine orchestration: spans -> lengths -> POA -> phase -> assign -> polish -> contract."""

from __future__ import annotations

import gzip
import json
import random
from pathlib import Path

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.assign import assign_reads, hybrid_references
from muc_one_span.hybrid.evidence import residual_sites
from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import draft_consensus, polish
from muc_one_span.hybrid.spans import Anchors, ReadRecord, SpanRead, categorize_reads
from muc_one_span.settings import RuntimeSettings
from muc_one_span.tools import run_tool

UNIT = 60


def read_input(path: Path) -> list[ReadRecord]:
    """Load FASTQ(.gz) or primary reads of a BAM."""
    if path.suffix == ".bam":
        text = run_tool(["samtools", "fastq", "-F", "0x900", str(path)])
        lines = text.splitlines()
    else:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as fh:  # type: ignore[operator]
            lines = fh.read().splitlines()
    return [ReadRecord(lines[i][1:].split()[0], lines[i + 1].upper(), lines[i + 3])
            for i in range(0, len(lines) - 3, 4)]


def _depth_status(spanning: int, s: RuntimeSettings) -> str:
    h = s.hybrid
    if spanning >= h.depth_adequate_spanning:
        return "adequate"
    return "low" if spanning >= h.depth_low_spanning else "insufficient"


def _allele_info(name: str, seq: str, members: list[SpanRead], assigned: int, basis: str,
                 residual: list[dict], settings: RuntimeSettings) -> dict:
    units = round(len(seq) / UNIT)
    status = _depth_status(len(members), settings)
    phase = {"length": "phased", "linked_sites": "phased", "none": "no_informative_heterozygosity",
             "unconfirmed_single_site": "unresolved_single_site"}[basis]
    return {
        "engine": "hybrid", "length": units, "canonical_repeats": units - 9,
        "fixed_repeat_count": 9, "reads": len(members), "spanning_reads": len(members),
        "assigned_reads": assigned, "support_basis": "spanning_reads",
        "depth_status": status, "length_basis": "spanning_peak", "consensus_basis": "poa",
        "split_basis": basis, "phase_status": phase, "residual_sites": residual,
        "reconstruction_status": "complete_segmentation" if status == "adequate" else "read_consensus_low_depth",
        "independent_haplotype_evidence": basis in ("length", "linked_sites"),
        "sequence_source": f"hybrid:{name}", "contig_name": f"hybrid_{name}",
        "vcf_path": None,
    }


def reconstruct_alleles(input_path: Path, output_dir: Path, rd: RepeatDictionary,
                        settings: RuntimeSettings) -> tuple[dict, dict[str, Path]]:
    h = settings.hybrid
    rng = random.Random(h.seed)
    backend = get_backend()
    reads = read_input(input_path)
    cats = categorize_reads(reads, Anchors.from_dictionary(rd), h)
    model = fit_length_model(cats.spanning, h)
    groups: list[tuple[list[SpanRead], str]] = []
    for peak in model.peaks:
        draft = draft_consensus(peak.members, h.n_poa, rng, backend)
        split = split_by_linked_sites(draft, peak.members, h, rng)
        basis = "length" if len(model.peaks) == 2 else split.basis
        groups.extend((g, basis if split.basis != "linked_sites" else "linked_sites") for g in split.groups)
    groups = sorted(groups, key=lambda g: -len(g[0]))[:2]
    drafts = {f"allele_{i + 1}": draft_consensus(g, h.n_poa, rng, backend) for i, (g, _) in enumerate(groups)}
    refs = hybrid_references(drafts, rd)
    non_spanning = cats.left_anchored + cats.right_anchored + cats.internal_or_offtarget
    extra = assign_reads(non_spanning, refs, h) if len(refs) > 1 else {k: [] for k in [*refs, "undecided"]}
    alleles: dict = {}
    paths: dict[str, Path] = {}
    for i, (members, basis) in enumerate(groups):
        name = f"allele_{i + 1}"
        pool = [m.seq for m in members] + extra.get(name, [])
        cons, info = polish(drafts[name], pool[: max(h.n_poa * 3, 120)], h.polish_rounds, h.hp_vote)
        residual = residual_sites(cons, [m.seq for m in members][:200], h.qc_residual_af)
        alleles[name] = _allele_info(name, cons, members, len(extra.get(name, [])), basis, residual, settings)
        alleles[name]["polish"] = info
        path = output_dir / f"consensus_{name}.fa"
        path.write_text(f">{name}\n{cons}\n")
        paths[name] = path
        (output_dir / f"consensus_{name}_context.json").write_text(json.dumps(
            {"engine": "hybrid", "reads": len(members), "full_consensus_path": str(path),
             "vcf_path": None}, indent=2))
    if len(alleles) == 1:
        alleles["allele_2"] = {**alleles["allele_1"], "candidate_duplicate_of": "allele_1",
                               "reconstruction_status": "not_separately_resolved"}
    hybrid_block = {
        "read_categories": cats.counts(), "rejected_peaks": model.rejected,
        "short_product_fraction": round(model.short_product_fraction, 4),
        "unassigned_spanning_fraction": round(model.unassigned_fraction, 4),
        "undecided_reads": len(extra.get("undecided", [])), "poa_backend": backend.name,
    }
    (output_dir / "hybrid_reads.json").write_text(json.dumps(hybrid_block, indent=2))
    (output_dir / "hybrid_references.fa").write_text(
        "".join(f">hybrid_{k}\n{v}\n" for k, v in refs.items()))
    lengths = [a["length"] for k, a in alleles.items() if k.startswith("allele_")]
    alleles.update({
        "homozygous": len(groups) == 1, "same_length": len(set(lengths)) == 1,
        "observed_length_candidates": [round(p.center_bp / UNIT) for p in model.peaks],
        "allele_multiplicity_status": "resolved" if len(groups) == 2 else "unresolved",
        "hybrid": hybrid_block,
    })
    return alleles, paths
```

`__init__.py`: `from muc_one_span.hybrid.engine import reconstruct_alleles` and `__all__ = ["reconstruct_alleles"]`.

Two constraints to check:
- **Ruff line length** is set in `pyproject.toml`. Wrap long lines when running `make format`.
- **File size:** keep `engine.py` below 650 lines. The code above is ~140 lines.

- [ ] **Step 4: Run the engine test**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_engine.py --no-cov -v`
Expected: PASS. If the sequences differ, print both and diff them before changing anything. Common causes are a length-window miss (Task 4) or insufficient polish reads.

- [ ] **Step 5: Extract the shared tail and branch the pipeline**

Move `pipeline.py` lines 142–217 (classification, validation, `alleles.json`/`repeats.*`/`summary.json` writing, and report) into `pipeline_tail.finish_run(...)`, keeping the same behaviour for the ladder path. Then, in `execute_pipeline`, after `write_run_configuration(...)`:

```python
    if settings.run.engine == "hybrid":
        from muc_one_span.hybrid import reconstruct_alleles

        check_tools(["samtools"] if str(input_path).endswith(".bam") else [])
        alleles_result, consensus_paths = reconstruct_alleles(
            Path(input_path), out, rd, settings)
        vcf_paths: dict[str, Path] = {}
        tool_versions = {"poa_backend": alleles_result["hybrid"]["poa_backend"]}
    else:
        # existing ladder steps 1-4, unchanged
        ...
    summary = finish_run(alleles_result, consensus_paths, vcf_paths, rd, settings, out, ...)
```

Inside `finish_run`, for hybrid runs (`alleles_result.get("hybrid")`):
- call `event_read_support(...)` per allele on its spanning member sequences with strands (store members on the side, e.g. `alleles_result[k]["_members"]`, and pop them before JSON writing);
- set each mutation's `read_support`, `vcf_support=False`, and `vcf_support_status="not_applicable_read_consensus"`;
- copy `alleles_result["hybrid"]` into `summary["hybrid"]`.

The report is called with `bam_path=None` for hybrid runs when IGV is off. When IGV is on, pass `fasta_path=out/"hybrid_references.fa"` and skip BAM tracks, with a notice.

- [ ] **Step 6: Run the full unit suite**

Run: `make test-fast`
Expected: PASS, including `test_wave1_cli.py` (the ladder path is untouched apart from the moved tail).

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/hybrid src/muc_one_span/pipeline.py src/muc_one_span/pipeline_tail.py tests/unit/hybrid/test_engine.py
git commit -m "feat(hybrid): engine orchestration, output contract, and pipeline branch"
```

---

### Task 10: Clinical gating on hybrid evidence, evaluation acceptance, harness forwarding

**Files:**
- Modify: `src/muc_one_span/report.py` (NO_PATHOGENIC gating), `src/muc_one_span/evaluation/artifacts.py:130-153`, `src/muc_one_span/evaluation/models.py`, `src/muc_one_span/benchmarking.py:102-126`, `src/muc_one_span/clinical_runner.py:144-158`, and the evaluation observation loader (read `read_support` into `Event`)
- Test: `tests/unit/test_clinical_decision.py`, `tests/unit/test_evaluation_artifacts.py`, `tests/unit/test_benchmark_tools.py`

**Interfaces:**
- Consumes: the hybrid fields from Task 9.
- Produces: `benchmarking.run_pipeline(..., engine: str | None = None)`, which appends `["--engine", engine]` only when it is not None, so the existing pinned argv test stays valid. The same applies to `clinical_runner.run_case`.

- [ ] **Step 1: Write failing tests**

```python
# tests/unit/test_clinical_decision.py (append)
def test_hybrid_rejected_peak_blocks_negative() -> None:
    s = _summary(dict(BASE, frameshift=False))
    s["classifications"]["allele_1"]["mutations"] = []
    s["hybrid"] = {"rejected_peaks": [{"units": 80, "support": 5, "reason": "support_below_threshold"}],
                   "unassigned_spanning_fraction": 0.05}
    for k in ("allele_1", "allele_2"):
        s["alleles"][k]["depth_status"] = "adequate"
    assert compute_clinical_decision(s)["status"] == "INCONCLUSIVE"


def test_hybrid_low_depth_blocks_negative() -> None:
    s = _summary(dict(BASE))
    s["classifications"]["allele_1"]["mutations"] = []
    s["hybrid"] = {"rejected_peaks": [], "unassigned_spanning_fraction": 0.0}
    s["alleles"]["allele_1"]["depth_status"] = "adequate"
    s["alleles"]["allele_2"]["depth_status"] = "low"
    assert compute_clinical_decision(s)["status"] == "INCONCLUSIVE"
```

```python
# tests/unit/test_evaluation_artifacts.py (append; reuse write_valid_artifacts helper at line 40)
def test_hybrid_statuses_count_as_completed(tmp_path) -> None:
    root = write_valid_artifacts(tmp_path)  # adapt to the helper's real signature
    summary = json.loads((root / "summary.json").read_text())
    for k in ("allele_1", "allele_2"):
        summary["alleles"][k]["phase_status"] = "no_informative_heterozygosity"
        summary["alleles"][k]["engine"] = "hybrid"
    (root / "summary.json").write_text(json.dumps(summary))
    obs = load_observation(root)
    assert obs.state != "ambiguous_reconstruction"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py tests/unit/test_evaluation_artifacts.py --no-cov -q`
Expected: the two hybrid decision tests FAIL (NO_PATHOGENIC). Adjust the evaluation test to the helper's real API first; it may already pass.

- [ ] **Step 3: Implement**

- `report.compute_clinical_decision`: before returning a negative, if `summary.get("hybrid")`, collect reasons:
  - any allele `depth_status != "adequate"`;
  - any allele `split_basis == "unconfirmed_single_site"`;
  - any rejected peak with `reason == "support_below_threshold"`;
  - `unassigned_spanning_fraction > 0.2`.
  If any reason exists, return INCONCLUSIVE with those reasons listed. Use the decision dict's existing reason field (inspect its current keys).
- `evaluation/artifacts._statuses_allow_completion`: `phase_status` stays as-is. Add `"read_consensus_low_depth"` **only** to a new `low_depth` observation flag, not to completion. Low depth scores as ambiguous by design.
- `evaluation/models.Event.supported`: add `read_support_status: str = field(default="unknown", compare=False)`. In `supported`, return `self.frameshift and self.template_match and (self.support_status == "exact_sequence_concordance" and self.vcf_support or self.read_support_status == "supported")`. Populate it in the observation loader from `mutation.get("read_support", {}).get("status")`.
- `benchmarking.run_pipeline` / `clinical_runner.run_case`: add the `engine: str | None = None` keyword and append `--engine` when set. `scripts/benchmark.py` and the clinical CLI get `--engine` passthrough options.

- [ ] **Step 4: Run the suites**

Run: `make test-fast`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/report.py src/muc_one_span/evaluation src/muc_one_span/benchmarking.py \
  src/muc_one_span/clinical_runner.py scripts tests/unit
git commit -m "feat(hybrid): gate negatives on hybrid evidence; score read support; forward --engine"
```

---

### Task 11: Integration test and frozen-panel, clinical, and in-house regression runs

**Files:**
- Create: `tests/integration/test_hybrid_engine.py`
- Create (outside Git): `../MucOneSpan-review-20260923/engine-regression/run.sh` and results

**Interfaces:**
- Consumes: CLI `muconespan run --engine hybrid`; `scripts/evaluate.py`; `simpanel/score.py`; `heldout/manifest.json`.

- [ ] **Step 1: Integration test on the repo's generated HiFi sample**

```python
# tests/integration/test_hybrid_engine.py
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main

pytestmark = [pytest.mark.integration]
DATA = Path("tests/data/generated/sample_close_51_58")


@pytest.mark.skipif(not DATA.exists(), reason="generated test data missing")
def test_hybrid_run_on_generated_sample(tmp_path: Path) -> None:
    bam = next(DATA.glob("*_amplicon_aligned.bam"))
    res = CliRunner().invoke(main, ["run", "-i", str(bam), "-o", str(tmp_path), "--engine", "hybrid",
                                    "--no-report"])
    assert res.exit_code == 0, res.output
    summary = json.loads((tmp_path / "summary.json").read_text())
    lengths = sorted(summary["alleles"][k]["length"] for k in ("allele_1", "allele_2"))
    assert lengths == [51, 58]
```

Run: `make test-int`
Expected: PASS. Report the test as skipped if samtools or the data is missing, and do not count a skip as validation.

- [ ] **Step 2: Frozen simulated panels (outside Git)**

```bash
R=../MucOneSpan-review-20260923
for panel in simpanel heldout; do
  mkdir -p $R/engine-regression/$panel
  python3 - "$R/$panel/manifest.json" <<'EOF' > $R/engine-regression/$panel/cases.tsv
import json,sys; m=json.load(open(sys.argv[1]))
cases = m["cases"] if isinstance(m, dict) else m
for c in cases: print(c["case_id"], c["reads"], c["platform"], sep="\t")
EOF
  while IFS=$'\t' read -r cid reads plat; do
    uv run --locked --all-extras muconespan run -i "$reads" -o $R/engine-regression/$panel/runs/$cid \
      --engine hybrid --platform $([ "$plat" = ont ] && echo ont || echo hifi) --no-report -t 1
  done < $R/engine-regression/$panel/cases.tsv
  python3 $R/simpanel/adapt_muconespan.py $R/engine-regression/$panel/runs $R/engine-regression/$panel/preds
  uv run --locked --all-extras python $R/simpanel/score.py --manifest $R/$panel/manifest.json \
    --pred $R/engine-regression/$panel/preds --out $R/engine-regression/$panel/scores --name hybrid
done
```

Check the manifest's actual key names (`case_id`, `reads`, `platform`) before running.
Expected: simpanel ≥ 80/80 alleles sequence-exact; held-out ≥ 76/80. The three known prototype failures must now be either correct or reported INCONCLUSIVE, not confident wrong calls.

- [ ] **Step 3: PRJEB92208 cohort (outside Git)**

Run the clinical harness with `--engine hybrid` on the 11 libraries, as for cohort-v3 but in a new output root. Expected:
- MP1–MP4 PATHOGENIC with `read_support.status == "supported"` at 78:17, 78:17, 44:7, 80:49;
- HG002 PCR and WGS 2/2 exact against Q100;
- HG001–HG004 not PATHOGENIC.

Record in `benchmarks/clinical/prjeb92208/hybrid-engine.json`: scores, hashes and versions only, with no sequences.

- [ ] **Step 4: In-house genomic (local only)**

Run `--engine hybrid --assay genomic` on `../MucOneSpan-review-20260923/inhouse/reads/*.fastq`. Expected:
- no phantom fragment alleles;
- The in-house dupC sample (ID kept local; motif 33) of the 82-unit allele reported with `depth_status` low or insufficient, and the decision INCONCLUSIVE, not PATHOGENIC and not NEGATIVE, at 7 spanning reads.

Keep the results out of Git.

- [ ] **Step 5: Commit the integration test and the public clinical record**

```bash
git add tests/integration/test_hybrid_engine.py benchmarks/clinical/prjeb92208/hybrid-engine.json
git commit -m "test(hybrid): integration test and PRJEB92208 regression record"
```

---

### Task 12: Documentation and changelog

**Files:**
- Modify: `docs/reference/cli.md`, `docs/guides/configuration.md`, `docs/reference/limitations.md`, `docs/getting-started/concepts.md`, `examples/runtime-settings.json`, `CHANGELOG.md`, `README.md` (one line on `--engine hybrid`, experimental)

- [ ] **Step 1: Write docs**

Cover:
- the algorithm stages (spec §3), the evidence fields (spec §4), and every `hybrid.*` setting with its default;
- that `hybrid` is experimental and not default until benchmark gating;
- genomic-assay depth caveats;
- the licence note (pyabpoa/pyspoa/edlib MIT; medaka/dorado not used).

Add a CHANGELOG "Unreleased" entry: explicit-support clinical decision (safety change), `--engine hybrid` (experimental), `--assay`.

- [ ] **Step 2: Verify**

Run: `make docs-check && make build-check && make ci-check`
Expected: all PASS; the distribution contains no new data files beyond the code.

- [ ] **Step 3: Commit**

```bash
git add docs examples CHANGELOG.md README.md
git commit -m "docs: document experimental hybrid engine, evidence fields and settings"
```

---

## Self-review notes

- **Spec coverage:**
  - §3 S1–S10 → Tasks 3–9. S8 and S10 are in Task 8; the optional Clair3-on-own-consensus QC (`--hybrid-qc clair3`) is **deferred** to a follow-up, and the residual Python QC covers the default path.
  - §3 "ladder-assisted length prior" and "ladder_seeded" consensus are **deferred** to a follow-up plan: they need the benchmark genomic low-depth strata to be tuned properly. The fields exist, and are always `spanning_peak` / `poa` until then.
  - §4 → Tasks 1, 9 and 10. §5 → Task 2. §6 → Task 2. §7 acceptance → Tasks 11–12, plus the benchmark plan.
- **Deferred items are listed explicitly:** Clair3-QC mode, ladder-assisted length and ladder-seeded consensus, and exact-window event thresholds (tuned in the benchmark plan).
