# Hybrid Read-Centric Engine Implementation Plan (v2, corrected)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task by task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `muconespan run --engine hybrid` (experimental, opt-in). It reconstructs MUC1 VNTR alleles from spanning reads (POA), assigns all other reads to allele drafts wrapped in ladder flanks, polishes, and reports explicit read-level evidence through the existing output contract and the P0 clinical gates. `ladder` stays the default engine.

**Architecture:** A new subpackage `muc_one_span/hybrid/` has one pure module per stage: `align` (lazy edlib), `spans` (S1), `lengths` (S2), `poa` + `polish` (S3/S7), `phase` (S4), `assign` (S5/S6), `evidence` (S8/S10) and `engine` (orchestration, output contract). `pipeline.execute_pipeline` branches on `settings.run.engine` right after `write_run_configuration` and before any tool check. Both engines share `pipeline_tail.finish_run` (classification, summary, report). Hybrid evidence reaches the clinical decision only through the existing P0 per-allele fields (`depth_status`, `selection_status`, `allele_genotype_status`, `length_status`, …) and `mutations[].read_support`. Engine-level data goes to top-level `summary["hybrid"]`, never into `alleles`.

**Tech Stack:** Python 3.10+, edlib (MIT), pyabpoa (MIT, default POA backend), pyspoa (MIT, alternative backend), pytest, Click, and the existing `classify_sequence`. No numpy and no scipy.

**Spec:** `.planning/2026-09-23-hybrid-engine-spec.md`. **Supersedes** `.planning/2026-09-23-hybrid-engine-plan.md`. It applies every finding and correction in `.planning/2026-09-24-hybrid-plan-preflight.md`; where v1 and the pre-flight disagree, v2 follows the pre-flight. The only exceptions are listed under "Decisions" below. Prototype (outside Git, reference only): `../MucOneSpan-review-20260923/poa-prototype/{proto.py,hetsplit.py}`, `homopolymer/hp_model.py`, `inhouse/hybrid/assign_test.py`.

**Base:** worktree `/home/bernt-popp/development/MucOneSpan-hybrid`, branch `feat/hybrid-engine`. It sits on `fix/p0-clinical-safety` HEAD `2675e9b` (v0.16.0).

**Validation of this plan (2026-09-24, not a substitute for executing it):**

- Every new module and test below, plus every diff, was run in a scratch copy of `2675e9b`. That copy used CPython 3.14.4 with edlib 1.3.9.post1, pyabpoa 1.5.7.1 and pyspoa 0.3.2.
- The full unit suite passed (1146 passed, `test_file_size.py` deselected). `ruff check`, `ruff format --check` and `mypy src/muc_one_span scripts` were clean.
- Branch coverage was 91% in total. **Without the hybrid tests it was 79%**, so the hybrid tests must run in the coverage job (see Decision D1).
- The integration test (Task 13) passed on the main checkout's generated `sample_close_51_58` BAM (lengths `[51, 58]`, under 1 s).
- Expected-failure claims were checked by running the new tests against the unmodified `2675e9b` sources.
- Code blocks are the validated, ruff-formatted text. The implementer still runs every step, because `uv.lock`, Docker and the other Python versions were **not** exercised.

**Smoke findings on generated MucOneUp 0.44.2 amplicon samples (not a benchmark):**

- 14 positives were run with `--engine hybrid`: 13 were PATHOGENIC with `read_support.status == "supported"`, and `dupc_60_80_cov50` was INCONCLUSIVE (low depth).
- 5 negatives were run: `gap3_50_53` gave NO_PATHOGENIC. `normal_50_55`, `normal_50_60`, `normal_100_120` and `normal_60_80` were INCONCLUSIVE, caused by `unresolved_single_site` or `unresolved_max_alleles`.
- The cause is minor read sub-populations (allele fraction 0.20–0.45) inside one length peak. In `normal_60_80` that is 15/51 reads carrying six linked substitutions.
- Polishing the draft before phasing did not change this, so it is not a draft artefact. Whether these are simulated PCR chimeras or something else is **not known**. Task 13 Step 2 must answer it before any threshold is tuned.
- The conservative direction (INCONCLUSIVE, never a false negative) is intended. Expect a low NO_PATHOGENIC yield until the question is answered.

## Global Constraints

- Every authored file has **fewer than 650 physical lines** (max 649). Split by responsibility and never compress lines. **Do not touch `alleles.py` (641) or `clinical_scoring.py` (626)**; nothing in this plan needs them. Estimated sizes after the plan:

  | File | Lines |
  | --- | --- |
  | `settings.py` | 518 |
  | `cli.py` | 531 |
  | `cli_run.py` | 134 |
  | `pipeline.py` | 219 |
  | `pipeline_tail.py` | 119 |
  | `report.py` | 433 |
  | `clinical_gates.py` | 112 |
  | `evaluation/artifacts.py` | 267 |
  | `benchmarking.py` | 334 |
  | `clinical_runner.py` | 231 |
  | `hybrid/*.py` | ≤ 292 each |
  | `test_clinical_decision.py` | 399 |

- `--engine` defaults to `ladder`. It stays the default until the benchmark decision rule (benchmark spec §6) passes on the sealed test split. The hybrid engine **adds fields only**; the ladder output and gates are unchanged, apart from the P0 `insufficient` gap fix in Task 3.
- New Python dependencies go only into the optional extra `hybrid`. `pyproject.toml` and `uv.lock` change together (`make lock`).
- External commands run only through `tools.py` (argument lists, no shell), with streamed output for reads (`run_tool_iter`). The hybrid default path runs no external tool for FASTQ input.
- Deterministic: every random choice uses `random.Random(settings.hybrid.seed)`, and synthetic tests use fixed seeds.
- Unit tests are synthetic, with no tool binaries and no patient data. Integration tests are `@pytest.mark.integration`. Report every skip with its reason, and never count a skip as validation.
- Never commit patient reads, sequences or LB-level outputs. Outputs from PRJEB92208, in-house runs and frozen panels stay outside Git.
- Run `make ci-check` before **every** commit, and add `make test-int`, `make docs-check`, `make security-check` and `make build-check` where a task says so. Every commit message ends with a blank line and `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Use `decision["state"]` (not `status`), and reuse the existing test helpers `_summary`, `BASE`, `RESOLVED_GATES` and `_gated_summary` (test_clinical_decision.py:140-221) without redefining them.

## Review Focus

1. **Minor read sub-populations inside one length peak.** The smoke runs show 4/5 generated negatives INCONCLUSIVE because of `unresolved_single_site` or `unresolved_max_alleles` (Task 8/11). Confirm the cause before tuning `het_af_min`, `het_min_group` or `link_phi_min` on the dev split. Never accept a silent drop to two alleles.
2. **Smear vs a real short allele (Task 5).** Below the major peak by more than 1.5 units, a candidate counts as an allele only with ≥`min_peak_reads` reads **and** prominence ≥`smear_min_prominence` over its shoulders. So a real short allele with fewer than 8 reads in an amplicon is called smear. Check that this trade-off is acceptable.
3. **A long allele with 5 spanning reads (Tasks 5, 11).** It must appear in `rejected_peaks` as `support_below_threshold`, set `selection_status = unresolved_rejected_peak` on both alleles, and give INCONCLUSIVE for a negative. The smoke sample `asymmetric_25_140` shows exactly this.
4. **Read-support producer contract (Task 9) vs the P0 gates (Task 3).** A homopolymer event needs n ≥ 20, LLR ≥ 10 and alt fraction ≥ 0.30, and no well-covered strand may have a negative LLR (otherwise `discordant`). A zero-read strand never fails. The event type comes from the dictionary template, not from "any 7C run". Only `supported` is support.
5. **`summary["hybrid"]` stays out of `alleles`, and `annotate_selection_qc` never runs on hybrid alleles.** Otherwise `evaluation/artifacts.py:229` reports `invalid_artifacts`, or `depth_status` becomes `not_assessed` and the depth gate is silently disabled.

## Decisions (made in v2; revisit only with evidence)

- **D1 edlib marker: remove it in BOTH extras.** Use `edlib>=1.3.9` with no `python_version` marker, in `hybrid` here and in `bench` on `feat/benchsim` (which currently has `edlib>=1.3; python_version < "3.14"`).
  - Why: coverage is measured on 3.14 only (test.yml:104-111), and without the hybrid tests total coverage is 79% (measured), below the 80% gate. The edlib sdist builds and imports on CPython 3.14.4 with system gcc (pre-flight (d), and again in this validation).
  - CI implications: on 3.14 the edlib and pyabpoa sdists compile in the `test-unit` job, using the ubuntu-24.04 gcc and the uv cache keyed on `uv.lock`.
  - Coordination: whichever of `feat/hybrid-engine` and `feat/benchsim` merges second must end with the identical line. Task 2 Step 7 records this.
- **D2 POA backend.** `hybrid.poa_backend` defaults to `pyabpoa` (the prototype evidence is pyabpoa-based) and `pyspoa` is selectable. There is **no silent fallback**. The backend name and package versions are recorded in `summary["hybrid"]` and `summary["tool_versions"]`, and a missing backend raises `ImportError` naming the extra.
- **D3 Imports.** `edlib`, `pyabpoa` and `spoa` are imported through `importlib.import_module`, so the package imports without the extra. No mypy overrides are needed; adding them only produces a "unused section" note.
- **D4 Docker.** The `app-builder` stage gets `gcc`, `libc6-dev` and `zlib1g-dev` through apt, because pyabpoa is sdist-only and the micromamba base has no compiler. Only `/opt/venv` is copied to the runtime, so the compilers do not ship.
- **D5 Single-site test unit.** Unit `Q` replaces the pre-flight's `C`. Measured: `C` differs from `X` at position 58, which also shortens the 7C run. That gives two features (`col` + `run`), so the test yields `linked_sites`, not `unconfirmed_single_site`. `Q` differs from `X` at position 6 only, outside any run. Sites are still counted per column (the pre-flight default): a unit type that differs at several bases (`A`, 4 bases) counts as linked sites, and a test documents this.
- **D6 Run-length candidate sites use `het_af_min`.** The threshold is `max(het_af_min, 4·background)`, not the prototype's `max(0.1, …)`. At 0.1, homopolymer stutter in generated HiFi amplicons produced `unconfirmed_single_site` in resolved samples (`dupa_60_80`).
- **D7 New provisional settings.**
  - `rejected_peak_noise_reads = 2`: a KDE maximum with ≤ 2 reads is `noise` and not gate-relevant. Without it, one outlier read makes every sample INCONCLUSIVE.
  - `n_poa_min` is **dropped** until ladder-seeded consensus exists; adding it later is additive.
  - Tune every hybrid default on the dev/val splits only.
- **D8 IGV.** `--engine hybrid` with `--report-igv embedded|sidecar` is rejected with `click.BadParameter`, because hybrid runs produce no BAM tracks. This is a new combination, so no existing behaviour changes.
- **D9 Harness engine field.** `benchmarking.run_pipeline` uses exactly the `feat/benchsim` signature and code (`engine: str = "ladder"`; `--engine` is appended only when the engine is not `ladder`; `engine` is always recorded). The clinical harness puts `engine`/`assay` into the **hashed** settings only when they are not the default, so existing ladder attempts stay resumable.
- **Deviations from spec §3, recorded in `engine.py`'s docstring:**
  - no ladder-assisted length prior and no ladder-seeded consensus (S2/S3);
  - depth is judged on spanning reads only (the spec's "or ≥ 40 assigned" is deferred);
  - read-level support uses spanning members only;
  - phasing uses a majority vote over linked sites (the prototype EM is not ported);
  - optional Clair3 QC is not included.

---

## File Structure

| Path | Responsibility |
| --- | --- |
| `src/muc_one_span/hybrid/__init__.py` | re-exports `reconstruct_alleles`, `HybridResult` |
| `src/muc_one_span/hybrid/align.py` | lazy edlib; `rc`, `cigar_ops`, `infix_hit`, `edit_distance(_infix)`, `Columns`, `project` (global or partial), `global_columns` |
| `src/muc_one_span/hybrid/spans.py` | S1: `ReadRecord`, `SpanRead`, `Anchors`, `ReadCategories`, `categorize_reads` |
| `src/muc_one_span/hybrid/lengths.py` | S2: `LengthPeak`, `LengthModel`, `window_bp`, `fit_length_model` (pure-Python KDE, smear rule, noise rule) |
| `src/muc_one_span/hybrid/poa.py` | `PoaBackend` protocol, `get_backend(name)` (pyabpoa, pyspoa) |
| `src/muc_one_span/hybrid/polish.py` | S3/S7: `draft_consensus`, `pileup_polish`, `homopolymer_vote`, `polish`, `_runs`, `read_run_length` |
| `src/muc_one_span/hybrid/assign.py` | S5/S6: `hybrid_references`, `assign_read`, `assign_reads`, `trim_to_draft` |
| `src/muc_one_span/hybrid/phase.py` | S4: `PhaseResult`, `split_by_linked_sites` |
| `src/muc_one_span/hybrid/evidence.py` | S8/S10: `residual_sites`, `homopolymer_llr`, `homopolymer_event_run`, `hp_status`, `event_read_support` |
| `src/muc_one_span/hybrid/engine.py` | `read_input`, `reconstruct_alleles` → `HybridResult`, `annotate_read_support`, `extra_versions` |
| `src/muc_one_span/pipeline_tail.py` | shared classify/summary/report tail (`finish_run`), moved from `pipeline.py:152-228` |
| `src/muc_one_span/cli_run.py` | the `run` command, moved from `cli.py:437-543`, plus `--engine/--assay` |
| `src/muc_one_span/{settings,clinical_gates,report,pipeline,cli,benchmarking,clinical_runner}.py`, `evaluation/{models,artifacts}.py`, `scripts/{benchmark,clinical_benchmark}.py` | modified (diffs below) |
| `tests/unit/hybrid/{__init__,synth,test_*}.py` | synthetic factory and one test module per stage |
| `tests/unit/test_hybrid_settings.py`, `tests/unit/test_cli_engine.py` | new settings and CLI tests |
| `tests/integration/test_hybrid_engine.py` | CLI on the generated HiFi BAM |

---

### Task 1: Explicit-support clinical decision — **DONE in P0 (v0.16.0)**

`clinical_gates.mutation_supported` accepts `read_support.status == "supported"` (clinical_gates.py:14-27), with tests at test_clinical_decision.py:178-200. Do not repeat it or re-announce it in the CHANGELOG.

---

### Task 2: Settings, optional extra `hybrid`, and tooling

**Files:**
- Modify: `src/muc_one_span/settings.py`, `pyproject.toml`, `uv.lock`, `Makefile` (lines 3-4), `docker/Dockerfile` (app-builder, lines 18-28), `conda/environment-dev.yml`
- Create: `tests/unit/test_hybrid_settings.py`

**Interfaces:**
- Consumes: `_integer`, `_number`, `_boolean`, `_choice` (settings.py:16-52).
- Produces:
  - `RunSettings.engine: str = "ladder"` (∈ {ladder, hybrid}) and `RunSettings.assay: str = "amplicon"` (∈ {amplicon, genomic});
  - `HybridSettings` (frozen dataclass; fields in the diff);
  - `RuntimeSettings.hybrid: HybridSettings`;
  - `_SECTIONS["hybrid"]`, so `load_settings`/`settings_as_dict` cover the new section automatically;
  - validation errors are named `hybrid.<field>` and `run.<field>`.

- [ ] **Step 1: Write the failing settings tests**

```python
# tests/unit/test_hybrid_settings.py
"""Hybrid-engine settings: defaults, validation, and configuration round trip."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from muc_one_span.settings import (
    DEFAULT_SETTINGS,
    HybridSettings,
    RunSettings,
    load_settings,
    settings_as_dict,
)


def test_hybrid_defaults() -> None:
    h = DEFAULT_SETTINGS.hybrid
    assert (h.anchor_max_edits, h.n_poa, h.assign_margin, h.depth_adequate_spanning) == (
        12,
        40,
        3,
        30,
    )
    assert (h.poa_backend, h.assign_max_error_rate, h.max_unassigned_spanning_fraction) == (
        "pyabpoa",
        0.15,
        0.2,
    )
    assert (h.smear_min_prominence, h.rejected_peak_noise_reads, h.hp_min_strand_reads) == (
        3.0,
        2,
        5,
    )
    assert (DEFAULT_SETTINGS.run.engine, DEFAULT_SETTINGS.run.assay) == ("ladder", "amplicon")


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("n_poa", 0),
        ("assign_margin", -1),
        ("het_af_min", 0.9),
        ("poa_backend", "medaka"),
        ("assign_max_error_rate", 1.5),
        ("smear_min_prominence", 0.5),
        ("max_span_units", 10),
        ("depth_adequate_spanning", 5),
        ("hp_vote", 1),
    ],
)
def test_hybrid_rejects_invalid(field: str, value: object) -> None:
    with pytest.raises(ValueError, match=f"hybrid.{field}"):
        HybridSettings(**{field: value})  # type: ignore[arg-type]


@pytest.mark.parametrize(("field", "value"), [("engine", "assembly"), ("assay", "wgs")])
def test_run_engine_and_assay_are_validated(field: str, value: str) -> None:
    with pytest.raises(ValueError, match=f"run.{field}"):
        RunSettings(**{field: value})  # type: ignore[arg-type]


def test_hybrid_roundtrip(tmp_path: Path) -> None:
    data = settings_as_dict(DEFAULT_SETTINGS)
    data["hybrid"]["n_poa"] = 25
    data["run"]["engine"] = "hybrid"
    path = tmp_path / "s.json"
    path.write_text(json.dumps(data))
    loaded = load_settings(path)
    assert loaded.hybrid.n_poa == 25 and loaded.run.engine == "hybrid"


def test_partial_hybrid_section_keeps_other_defaults(tmp_path: Path) -> None:
    path = tmp_path / "s.json"
    path.write_text('{"schema_version": 1, "hybrid": {"poa_backend": "pyspoa"}}')
    loaded = load_settings(path)
    assert loaded.hybrid.poa_backend == "pyspoa" and loaded.hybrid.n_poa == 40
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/test_hybrid_settings.py --no-cov -q`
Expected: collection error, `ImportError: cannot import name 'HybridSettings' from 'muc_one_span.settings'`.

- [ ] **Step 3: Implement the settings**

```diff
--- a/src/muc_one_span/settings.py
+++ b/src/muc_one_span/settings.py
@@ -65,6 +65,8 @@
     mapping_timeout: float = 3600.0
     report: bool = False
     report_igv: str = "off"
+    engine: str = "ladder"
+    assay: str = "amplicon"
 
     def __post_init__(self) -> None:
         _integer("run.threads", self.threads, 1)
@@ -81,6 +83,8 @@
         _boolean("run.report", self.report)
         if self.report_igv not in ("off", "embedded", "sidecar"):
             raise ValueError("run.report_igv must be off, embedded, or sidecar")
+        _choice("run.engine", self.engine, ("ladder", "hybrid"))
+        _choice("run.assay", self.assay, ("amplicon", "genomic"))
 
 
 @dataclass(frozen=True)
@@ -306,6 +310,86 @@
 
 
 @dataclass(frozen=True)
+class HybridSettings:
+    """Read-centric engine thresholds (spec 2026-09-23 §5).
+
+    Every default is provisional (prototype-derived) and is tuned on the benchmark
+    dev/validation splits only; the sealed test split never informs a default.
+    """
+
+    anchor_max_edits: int = 12
+    min_span_units: int = 15
+    max_span_units: int = 160
+    peak_window_base_bp: float = 30.0
+    peak_window_per_unit_bp: float = 0.6
+    min_peak_reads: int = 8
+    far_peak_min_frac: float = 0.03
+    near_peak_min_frac: float = 0.20
+    smear_min_prominence: float = 3.0
+    rejected_peak_noise_reads: int = 2
+    n_poa: int = 40
+    poa_backend: str = "pyabpoa"
+    polish_rounds: int = 2
+    hp_vote: bool = True
+    het_af_min: float = 0.2
+    het_min_group: float = 0.15
+    link_phi_min: float = 0.5
+    min_linked_sites: int = 2
+    min_fragment_bp: int = 1000
+    assign_margin: int = 3
+    assign_max_error_rate: float = 0.15
+    qc_residual_af: float = 0.25
+    max_unassigned_spanning_fraction: float = 0.2
+    depth_adequate_spanning: int = 30
+    depth_low_spanning: int = 10
+    hp_llr_min: float = 10.0
+    hp_min_reads: int = 20
+    hp_min_alt_frac: float = 0.30
+    hp_min_strand_reads: int = 5
+    seed: int = 1
+
+    def __post_init__(self) -> None:
+        for name in (
+            "anchor_max_edits",
+            "min_peak_reads",
+            "rejected_peak_noise_reads",
+            "polish_rounds",
+            "min_linked_sites",
+            "min_fragment_bp",
+            "assign_margin",
+            "depth_low_spanning",
+            "hp_min_reads",
+            "hp_min_strand_reads",
+            "seed",
+        ):
+            _integer(f"hybrid.{name}", getattr(self, name))
+        _integer("hybrid.n_poa", self.n_poa, 1)
+        _integer("hybrid.min_span_units", self.min_span_units, 1)
+        _integer("hybrid.max_span_units", self.max_span_units, self.min_span_units + 1)
+        _integer(
+            "hybrid.depth_adequate_spanning", self.depth_adequate_spanning, self.depth_low_spanning
+        )
+        _number("hybrid.peak_window_base_bp", self.peak_window_base_bp, 1)
+        _number("hybrid.peak_window_per_unit_bp", self.peak_window_per_unit_bp)
+        _number("hybrid.smear_min_prominence", self.smear_min_prominence, 1)
+        for name in (
+            "far_peak_min_frac",
+            "near_peak_min_frac",
+            "het_min_group",
+            "qc_residual_af",
+            "hp_min_alt_frac",
+            "link_phi_min",
+            "assign_max_error_rate",
+            "max_unassigned_spanning_fraction",
+        ):
+            _number(f"hybrid.{name}", getattr(self, name), 0, 1)
+        _number("hybrid.het_af_min", self.het_af_min, 0.01, 0.5)
+        _number("hybrid.hp_llr_min", self.hp_llr_min)
+        _boolean("hybrid.hp_vote", self.hp_vote)
+        _choice("hybrid.poa_backend", self.poa_backend, ("pyabpoa", "pyspoa"))
+
+
+@dataclass(frozen=True)
 class RuntimeSettings:
     """Complete schema-one settings; sections remain immutable when passed to workers."""
 
@@ -318,6 +402,7 @@
     calling: CallingSettings = field(default_factory=CallingSettings)
     read_phasing: ReadPhasingSettings = field(default_factory=ReadPhasingSettings)
     reference_layout: ReferenceLayoutSettings = field(default_factory=ReferenceLayoutSettings)
+    hybrid: HybridSettings = field(default_factory=HybridSettings)
     repeat_dictionary: str | None = None
 
     def __post_init__(self) -> None:
@@ -338,6 +423,7 @@
     "calling": CallingSettings,
     "read_phasing": ReadPhasingSettings,
     "reference_layout": ReferenceLayoutSettings,
+    "hybrid": HybridSettings,
 }
 DEFAULT_SETTINGS = RuntimeSettings()
 DEFAULT_LAYOUT = DEFAULT_SETTINGS.reference_layout
```

- [ ] **Step 4: Run the settings tests and the existing settings suites**

Run: `uv run --locked --all-extras pytest tests/unit/test_hybrid_settings.py tests/unit/test_runtime_settings.py tests/unit/test_cli_settings.py tests/unit/test_api_defaults.py --no-cov -q`
Expected: PASS.

- [ ] **Step 5: Add the optional extra (D1–D3)**

In `pyproject.toml` `[project.optional-dependencies]`, after `report`:

```toml
hybrid = [
    "edlib>=1.3.9",
    "pyabpoa>=1.5.7",
    "pyspoa>=0.3.2",
]
```

Do **not** add numpy or mypy overrides (D3). In the `Makefile`:

```make
UV_TEST = $(UV_RUN) --group test --extra report --extra hybrid
UV_QUALITY = $(UV_RUN) --group quality --extra report --extra hybrid
```

(If `feat/benchsim` is already merged, the lines read `--extra report --extra bench --extra hybrid`.)

Run: `make lock && make dev`
Expected: `uv.lock` gains `edlib`, `pyabpoa` and `pyspoa` only; review the lock diff. pyabpoa and edlib (on 3.14) build from source. If a build fails, stop and report the compiler error. Do not drop pyabpoa silently, because the prototype evidence depends on it (D2).

- [ ] **Step 6: Verify the extra builds on the CI interpreters**

Run:

```bash
for v in 3.10 3.13 3.14; do
  uv run --locked --no-default-groups --group test --extra report --extra hybrid --python "$v" \
    python -c "import edlib, pyabpoa, spoa; print('$v ok')"
done
```

Expected: three `ok` lines. pyabpoa on 3.10–3.13 is **unverified** by the pre-flight, so this step is the verification. If an interpreter fails, record it and ask the user before adding any marker. A marker on 3.14 would bring the coverage problem back (79%).

- [ ] **Step 7: Docker, conda, and the benchsim note**

In `docker/Dockerfile`, stage `app-builder`, right after `USER root` (line 18):

```dockerfile
# pyabpoa is sdist-only; compilers stay in the builder, only /opt/venv is copied.
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
```

Change both `uv sync` lines (25 and 28) from `--extra report` to `--extra report --extra hybrid`. In `conda/environment-dev.yml`, add `edlib`, `pyabpoa` and `pyspoa` under `dependencies` (bioconda/conda-forge binaries) for conda-based development.

Append this line to `.planning/2026-09-24-session-handoff.md`: "feat/benchsim `bench` extra must use `edlib>=1.3.9` (no marker) to match `hybrid` (plan v2 D1)."

Run: `make docker-test`
Expected: PASS. If Docker is unavailable, report the step as **untested**.

- [ ] **Step 8: Full check and commit**

Run: `make ci-check && make security-check && make build-check`
Expected: PASS. `settings.py` is ≈518 lines.

```bash
git add pyproject.toml uv.lock Makefile docker/Dockerfile conda/environment-dev.yml \
  src/muc_one_span/settings.py tests/unit/test_hybrid_settings.py
git commit -m "feat(settings): add hybrid engine settings, run.engine/run.assay and the hybrid extra

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: P0 gate reconciliation and the read-support producer contract

This task is pure summary-dict work. It fixes the latent P0 gap: `depth_status == "insufficient"` was recognised nowhere. It also fixes the gate vocabulary that the engine will produce.

**Files:**
- Modify: `src/muc_one_span/clinical_gates.py`, `src/muc_one_span/report.py` (lines 22-26, 98, 116)
- Test: `tests/unit/test_clinical_decision.py`, `tests/unit/test_clinical_gates.py`

**Interfaces:**
- Produces:
  - `clinical_gates.LOW_DEPTH_STATUSES = frozenset({"low", "insufficient"})`;
  - `clinical_gates.READ_SUPPORT_STATUSES` = {supported, insufficient_depth, discordant, not_supported, not_localized};
  - `_GENOTYPE_REASONS["residual_heterogeneity"]`;
  - the depth message uses `depth_basis`. The default `primary_alignment_records` keeps the text "12 primary alignments" that test_clinical_decision.py:307,328 assert; `spanning_reads` gives "12 spanning reads";
  - the selection message omits "secondary mode fraction" when it is `None`, and appends `selection_detail`;
  - `mutation_blockers` gives `read-level support <status>` when a `read_support` dict is present but not supported.
- Producer contract (implemented in Task 9, consumed unchanged by `mutation_supported`): only `read_support.status == "supported"` is support.
  - Homopolymer events (dictionary template = single-base indel inside a run ≥ 4): n ≥ `hp_min_reads` (20), stutter-aware LLR ≥ `hp_llr_min` (10), alt fraction ≥ `hp_min_alt_frac` (0.30), and no strand with ≥ `hp_min_strand_reads` (5) reads has a negative strand LLR (otherwise `discordant`). A zero-read strand never fails the event.
  - Other templated events: parent-vs-template competition, `supported` iff n ≥ 20, alt/n ≥ 0.30 and alt > ref.
  - n < 20 gives `insufficient_depth`, which is blocked, never negative.
  - A `supported` event on an allele whose `depth_status` is in `LOW_DEPTH_STATUSES` is still blocked by the carrier gate.

- [ ] **Step 1: Write the failing tests**

Add `import pytest` to the imports of `tests/unit/test_clinical_decision.py`, and append:

```python
HYBRID_GATES = dict(
    RESOLVED_GATES,
    depth_basis="spanning_reads",
    spanning_reads=120,
    secondary_mode_fraction=None,
    allele_genotype_status="not_applicable_read_consensus",
    engine="hybrid",
)


def _hybrid_summary(**allele2: object) -> dict:
    summary = _gated_summary()
    for key in ("allele_1", "allele_2"):
        summary["alleles"][key].update(HYBRID_GATES)
    summary["alleles"]["allele_2"].update(allele2)
    return summary


def test_hybrid_resolved_is_negative() -> None:
    assert compute_clinical_decision(_hybrid_summary())["state"] == (
        "NO_PATHOGENIC_VARIANT_DETECTED"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("depth_status", "insufficient"),
        ("selection_status", "unresolved_rejected_peak"),
        ("selection_status", "unresolved_unassigned_spanning"),
        ("selection_status", "unresolved_single_site"),
        ("selection_status", "unresolved_max_alleles"),
        ("allele_genotype_status", "residual_heterogeneity"),
    ],
)
def test_hybrid_gate_blocks_negative(field: str, value: str) -> None:
    decision = compute_clinical_decision(_hybrid_summary(**{field: value}))
    assert decision["state"] == "INCONCLUSIVE"


def test_insufficient_depth_carrier_blocks_pathogenic() -> None:
    mutation = dict(
        BASE,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "supported"},
    )
    summary = _hybrid_summary()
    summary["classifications"]["allele_1"]["mutations"] = [mutation]
    assert compute_clinical_decision(summary)["state"] == "PATHOGENIC"
    summary["alleles"]["allele_1"]["depth_status"] = "insufficient"
    assert compute_clinical_decision(summary)["state"] == "INCONCLUSIVE"


def test_hybrid_depth_message_names_spanning_reads() -> None:
    decision = compute_clinical_decision(_hybrid_summary(depth_status="low", spanning_reads=12))
    assert any("Allele 2: 12 spanning reads" in detail for detail in decision["details"])


def test_selection_message_omits_missing_secondary_fraction() -> None:
    decision = compute_clinical_decision(
        _hybrid_summary(selection_status="unresolved_rejected_peak", selection_detail="peak 80u")
    )
    reason = next(d for d in decision["details"] if "allele selection unresolved" in d)
    assert "secondary mode fraction" not in reason and "peak 80u" in reason
```

Append to `tests/unit/test_clinical_gates.py` (it already imports `allele_gate_reasons` and `mutation_blockers`):

```python
_TEMPLATED = {
    "frameshift": True,
    "template_match": True,
    "mutation_name": "dupC",
    "localization_status": "exact",
    "vcf_support": False,
    "vcf_support_status": "not_applicable_read_consensus",
}


def test_read_support_statuses_name_the_blocker() -> None:
    for status in ("insufficient_depth", "discordant", "not_supported", "not_localized"):
        blockers = mutation_blockers({**_TEMPLATED, "read_support": {"status": status}})
        assert blockers == [f"read-level support {status}"]
    assert mutation_blockers({**_TEMPLATED, "read_support": {"status": "supported"}}) == []


def test_insufficient_depth_is_gated_like_low() -> None:
    info = {
        "depth_status": "insufficient",
        "depth_basis": "spanning_reads",
        "spanning_reads": 7,
        "depth_threshold": 30,
    }
    assert allele_gate_reasons(info, "Allele 1") == [
        "Allele 1: 7 spanning reads, below the per-allele depth gate (30)."
    ]
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py tests/unit/test_clinical_gates.py --no-cov -q`
Expected: 7 failures:
- `test_hybrid_gate_blocks_negative[depth_status-insufficient]`
- `test_hybrid_gate_blocks_negative[allele_genotype_status-residual_heterogeneity]`
- `test_insufficient_depth_carrier_blocks_pathogenic`
- `test_hybrid_depth_message_names_spanning_reads`
- `test_selection_message_omits_missing_secondary_fraction`
- `test_read_support_statuses_name_the_blocker`
- `test_insufficient_depth_is_gated_like_low`

The four `unresolved_*` parametrizations already pass on P0 and act as regression guards.

- [ ] **Step 3: Implement**

```diff
--- a/src/muc_one_span/clinical_gates.py
+++ b/src/muc_one_span/clinical_gates.py
@@ -9,6 +9,16 @@
 from typing import Any
 
 SUPPORTED_VCF_STATUSES = frozenset({"exact_sequence_concordance"})
+# Per-allele depth statuses that block a negative call and a PATHOGENIC carrier.
+LOW_DEPTH_STATUSES = frozenset({"low", "insufficient"})
+# read_support.status values a producer may emit; only "supported" is support.
+READ_SUPPORT_STATUSES = frozenset(
+    {"supported", "insufficient_depth", "discordant", "not_supported", "not_localized"}
+)
+_DEPTH_BASIS_LABELS = {
+    "primary_alignment_records": "primary alignments",
+    "spanning_reads": "spanning reads",
+}
 
 
 def mutation_supported(mutation: dict[str, Any]) -> bool:
@@ -35,6 +45,10 @@
         "consensus uses unresolved (IUPAC) selection"
     ),
     "unresolved_genotype_records": "conflicting or incomplete genotype records",
+    "residual_heterogeneity": (
+        "residual read heterogeneity on this allele "
+        "(possible unresolved mixture, chimera or mosaicism)"
+    ),
 }
 
 
@@ -48,12 +62,17 @@
     if mutation.get("localization_status") == "ambiguous":
         blockers.append("localization ambiguous")
     if not mutation_supported(mutation):
+        read_support = mutation.get("read_support")
         status = mutation.get("vcf_support_status")
-        blockers.append(
-            "heterozygous genotype not resolved to one allele"
-            if status == "heterozygous_genotype_unresolved"
-            else f"no explicit sequence-level support ({status or 'status unavailable'})"
-        )
+        if isinstance(read_support, dict):
+            state = read_support.get("status") or "status unavailable"
+            blockers.append(f"read-level support {state}")
+        elif status == "heterozygous_genotype_unresolved":
+            blockers.append("heterozygous genotype not resolved to one allele")
+        else:
+            blockers.append(
+                f"no explicit sequence-level support ({status or 'status unavailable'})"
+            )
     return blockers
 
 
@@ -64,10 +83,10 @@
     reasons: list[str] = []
     selection = info.get("selection_status")
     if isinstance(selection, str) and selection.startswith("unresolved"):
-        reasons.append(
-            f"{label}: allele selection unresolved ({selection}; secondary mode fraction "
-            f"{info.get('secondary_mode_fraction')})."
-        )
+        fraction = info.get("secondary_mode_fraction")
+        detail = "" if fraction is None else f"; secondary mode fraction {fraction}"
+        extra = f" {info['selection_detail']}" if info.get("selection_detail") else ""
+        reasons.append(f"{label}: allele selection unresolved ({selection}{detail}).{extra}")
     length, reference_length = info.get("length"), info.get("reference_length")
     length_status = info.get("length_status")
     if length_status is None:  # Legacy summary without selection_qc: compare conservatively.
@@ -81,9 +100,10 @@
             f"{label}: reported length {length} differs from the consensus contig length "
             f"{reference_length}."
         )
-    if info.get("depth_status") == "low":
+    if info.get("depth_status") in LOW_DEPTH_STATUSES:
+        basis = info.get("depth_basis") or "primary_alignment_records"
         reasons.append(
-            f"{label}: {info.get('primary_alignment_records')} primary alignments, below the "
+            f"{label}: {info.get(basis)} {_DEPTH_BASIS_LABELS.get(basis, basis)}, below the "
             f"per-allele depth gate ({info.get('depth_threshold')})."
         )
     genotype = info.get("allele_genotype_status")
```

```diff
--- a/src/muc_one_span/report.py
+++ b/src/muc_one_span/report.py
@@ -21,6 +21,7 @@
 
 from muc_one_span.clinical_gates import (
     LEGACY_MIN_TOTAL_READS,
+    LOW_DEPTH_STATUSES,
     allele_gate_reasons,
     mutation_blockers,
 )
@@ -95,7 +96,9 @@
     a1 = alleles.get("allele_1", {}) if isinstance(alleles, dict) else {}
     a2 = alleles.get("allele_2", {}) if isinstance(alleles, dict) else {}
     carriers = {"allele_1": a1, "allele_2": a2}
-    depth_assessed = any(a.get("depth_status") in ("adequate", "low") for a in (a1, a2))
+    depth_assessed = any(
+        a.get("depth_status") in ("adequate", *LOW_DEPTH_STATUSES) for a in (a1, a2)
+    )
     total_reads = (a1.get("reads", 0) or 0) + (a2.get("reads", 0) or 0)
     low_coverage = (
         not depth_assessed and total_reads < LEGACY_MIN_TOTAL_READS and (bool(a1) or bool(a2))
@@ -113,7 +116,7 @@
             mut_copy = dict(mut)
             mut_copy["allele"] = allele_key
             blockers = mutation_blockers(mut)
-            if carrier.get("depth_status") == "low":
+            if carrier.get("depth_status") in LOW_DEPTH_STATUSES:
                 blockers.append("carrying allele is below the per-allele depth gate")
             if low_coverage:
                 blockers.append("total read depth is below the diagnostic threshold")
```

Do **not** add a `summary["hybrid"]` branch to `report.py`. Hybrid sample-level blockers reach the decision only through `selection_status` and `allele_genotype_status`.

- [ ] **Step 4: Run the decision, gate and report suites**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py tests/unit/test_clinical_gates.py tests/unit/test_report.py tests/unit/test_report_wave1.py --no-cov -q`
Expected: PASS.

- [ ] **Step 5: Full check and commit**

Run: `make ci-check`
Expected: PASS.

```bash
git add src/muc_one_span/clinical_gates.py src/muc_one_span/report.py \
  tests/unit/test_clinical_decision.py tests/unit/test_clinical_gates.py
git commit -m "fix(gates): gate insufficient depth, residual heterogeneity and read-level support

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Synthetic factory, lazy alignment primitives, and span categorisation (S1)

**Files:**
- Create: `src/muc_one_span/hybrid/__init__.py` (docstring only for now), `src/muc_one_span/hybrid/align.py`, `src/muc_one_span/hybrid/spans.py`
- Create: `tests/unit/hybrid/__init__.py` (empty), `tests/unit/hybrid/synth.py`, `tests/unit/hybrid/test_spans.py`

**Interfaces:**
- Produces (align):
  - `HYBRID_HINT: str`
  - `rc(seq: str) -> str`
  - `cigar_ops(cigar: str) -> list[tuple[int, str]]`
  - `infix_hit(query: str, target: str, k: int) -> tuple[int, int, int] | None`
  - `edit_distance_infix(query: str, target: str) -> int`
  - `edit_distance(a: str, b: str) -> int`
  - `Columns(cols: list[str | None], ins: dict[int, str], t2q: list[int], first: int, last: int)` with `covers(start: int, end: int) -> bool`
  - `project(read: str, cons: str, *, partial: bool = False) -> Columns`
  - `global_columns(read: str, cons: str) -> Columns`
- Produces (spans):
  - `ReadRecord(name, seq, qual)`
  - `SpanRead(name, seq, mean_q, strand, anchor_edits, anchor_basis)` with `.length`
  - `Anchors.from_dictionary(rd) -> Anchors`
  - `ReadCategories(...)` with `.counts() -> dict[str, int]`
  - `categorize_reads(reads: Iterable[ReadRecord], anchors: Anchors, settings: HybridSettings) -> ReadCategories`
- Produces (tests):
  - `synth.RD`, `synth.PRE`, `synth.POST`
  - `synth.dupc(unit: str = "X") -> str` (the dictionary dupC template)
  - `synth.allele(inner: list[str]) -> str` (a token is a repeat ID or a raw unit sequence)
  - `synth.reads(allele_seq, n, *, err, seed, strand_mix=True, flank_bp=40, smear_frac=0.0) -> list[ReadRecord]`

The `edlib` import is lazy, so `synth`, `rc` and `cigar_ops` import without the extra (pre-flight C3.2). Tests that align call `pytest.importorskip("edlib")` at module level. Because D1 installs the extra everywhere, these skips should never happen in CI; a skip there is a configuration error.

- [ ] **Step 1: Write the synthetic factory**

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


def dupc(unit: str = "X") -> str:
    """Parent unit carrying the dictionary dupC template (7C->8C; insert at 1-based 60)."""
    return next(
        seq
        for seq, (parent, name) in RD.mutated_sequences.items()
        if parent == unit and name == "dupC"
    )


def allele(inner: list[str]) -> str:
    """Motif 1..9 sequence: pre-repeats + inner tokens + after-repeats.

    A token is a repeat ID from the bundled dictionary or, if it is not an ID, a raw
    unit sequence (for example ``dupc()``).
    """
    return "".join(RD.repeats.get(u, u) for u in PRE + inner + POST)


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


def reads(
    allele_seq: str,
    n: int,
    *,
    err: float,
    seed: int,
    strand_mix: bool = True,
    flank_bp: int = 40,
    smear_frac: float = 0.0,
) -> list[ReadRecord]:
    """``n`` noisy reads of flank + allele + flank; smear reads lose an internal block."""
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

- [ ] **Step 2: Write the failing tests**

```python
# tests/unit/hybrid/test_spans.py
"""S1 anchor search, orientation and read categories on synthetic reads."""

from __future__ import annotations

import sys

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.hybrid.align import infix_hit
from muc_one_span.hybrid.spans import Anchors, ReadRecord, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

ANCH = Anchors.from_dictionary(synth.RD)
S = HybridSettings()


def test_dupc_template_is_classified() -> None:
    seq = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    names = {
        (m["mutation_name"], m.get("template_match"))
        for m in classify_sequence(seq, synth.RD)["mutations_detected"]
    }
    assert ("dupC", True) in names


def test_missing_edlib_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "edlib", None)
    with pytest.raises(ImportError, match="hybrid"):
        infix_hit("ACGT", "ACGT", 0)


def test_spanning_reads_oriented_and_trimmed_to_motif1_motif9() -> None:
    seq = synth.allele(["X"] * 30)
    cats = categorize_reads(synth.reads(seq, 20, err=0.02, seed=1), ANCH, S)
    assert len(cats.spanning) == 20
    for sp in cats.spanning:
        assert abs(sp.length - len(seq)) <= 40
        assert sp.anchor_basis == "motif"
    assert {sp.strand for sp in cats.spanning} == {"+", "-"}


def test_fragments_and_offtarget_are_not_spanning() -> None:
    seq = synth.allele(["X"] * 30)
    frag = ReadRecord("frag", seq[:900], "5" * 900)
    junk = ReadRecord("junk", "ACGT" * 300, "5" * 1200)
    cats = categorize_reads([frag, junk], ANCH, S)
    assert cats.counts() == {
        "spanning": 0,
        "left_anchored": 1,
        "right_anchored": 0,
        "internal_or_offtarget": 1,
    }


def test_mutated_motif1_falls_back_to_flank_anchor() -> None:
    seq = synth.allele(["X"] * 30)
    broken = seq[:5] + "T" * 40 + seq[45:]  # destroy most of motif 1
    read_seq = synth.RD.flanking_left[-40:] + broken + synth.RD.flanking_right[:40]
    cats = categorize_reads(
        [ReadRecord("m", read_seq, "5" * len(read_seq))], ANCH, HybridSettings(anchor_max_edits=6)
    )
    assert len(cats.spanning) == 1 and cats.spanning[0].anchor_basis == "flank"
```

- [ ] **Step 3: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_spans.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid'`.

- [ ] **Step 4: Implement `align.py`**

```python
# src/muc_one_span/hybrid/align.py
"""edlib-based alignment primitives shared by the hybrid engine stages.

edlib is imported lazily: ``rc`` and ``cigar_ops`` (and the synthetic test factory)
work without the optional ``hybrid`` extra, and alignment calls fail with a clear hint.
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from types import ModuleType

HYBRID_HINT = "The hybrid engine needs the 'hybrid' extra: pip install 'muc_one_span[hybrid]'"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _edlib() -> ModuleType:
    try:
        return importlib.import_module("edlib")
    except ImportError as exc:
        raise ImportError(HYBRID_HINT) from exc


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
    res = _edlib().align(query, target, mode="HW", task="locations", k=k)
    if res["editDistance"] < 0:
        return None
    start, end = res["locations"][0]
    return int(start), int(end) + 1, int(res["editDistance"])


def edit_distance_infix(query: str, target: str) -> int:
    """Edit distance of query aligned fully inside target (semi-global)."""
    return int(_edlib().align(query, target, mode="HW", task="distance")["editDistance"])


def edit_distance(a: str, b: str) -> int:
    """Global edit distance."""
    return int(_edlib().align(a, b, mode="NW", task="distance")["editDistance"])


@dataclass
class Columns:
    """A read projected onto consensus columns.

    ``cols[i]`` is the read base, ``"-"`` for a deletion, or ``None`` where a partial
    read does not cover column ``i``. ``ins[i]`` holds read bases inserted before
    column ``i`` (``i == len(cons)`` at the end). ``t2q`` maps consensus positions
    (``len(cons) + 1`` entries) to read positions, ``-1`` where uncovered.
    """

    cols: list[str | None]
    ins: dict[int, str]
    t2q: list[int]
    first: int
    last: int  # covered consensus interval [first, last)

    def covers(self, start: int, end: int) -> bool:
        """True when the read covers [start - 1, end] (one anchoring column each side)."""
        return self.first <= max(start - 1, 0) and min(end + 1, len(self.cols)) <= self.last


def project(read: str, cons: str, *, partial: bool = False) -> Columns:
    """Project read onto consensus columns (global, or infix for partial reads)."""
    res = _edlib().align(read, cons, mode="HW" if partial else "NW", task="path")
    n = len(cons)
    first = int(res["locations"][0][0]) if partial else 0
    cols: list[str | None] = [None] * n
    t2q = [-1] * (n + 1)
    ins: dict[int, str] = {}
    ti, qi = first, 0
    for length, op in cigar_ops(res["cigar"]):
        if op in "=X":
            for _ in range(length):
                cols[ti] = read[qi]
                t2q[ti] = qi
                ti += 1
                qi += 1
        elif op == "I":
            ins[ti] = ins.get(ti, "") + read[qi : qi + length]
            qi += length
        else:
            for _ in range(length):
                cols[ti] = "-"
                t2q[ti] = qi
                ti += 1
    t2q[ti] = qi
    return Columns(cols, ins, t2q, first, ti)


def global_columns(read: str, cons: str) -> Columns:
    """Globally align a spanning read to cons and project onto consensus columns."""
    return project(read, cons)
```

`src/muc_one_span/hybrid/__init__.py` for now:

```python
"""Experimental read-centric (hybrid) MUC1 VNTR reconstruction engine."""
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
UNIT = 60


@dataclass(frozen=True)
class ReadRecord:
    """One input read as stored in FASTQ."""

    name: str
    seq: str
    qual: str


@dataclass(frozen=True)
class SpanRead:
    """A spanning read trimmed to motif 1..motif 9 and oriented to the VNTR forward strand."""

    name: str
    seq: str
    mean_q: float
    strand: str  # "+" read was already forward, "-" read was reverse-complemented
    anchor_edits: int
    anchor_basis: str  # "motif" or "flank"

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass(frozen=True)
class Anchors:
    """Motif 1/9 anchors plus flank anchors used when a motif carries a mutation."""

    left: str
    right: str
    left_flank: str
    right_flank: str

    @classmethod
    def from_dictionary(cls, rd: RepeatDictionary) -> Anchors:
        return cls(
            rd.repeats["1"],
            rd.repeats["9"],
            rd.flanking_left[-FLANK_ANCHOR_BP:],
            rd.flanking_right[:FLANK_ANCHOR_BP],
        )


@dataclass
class ReadCategories:
    """Reads grouped by anchor evidence."""

    spanning: list[SpanRead] = field(default_factory=list)
    left_anchored: list[ReadRecord] = field(default_factory=list)
    right_anchored: list[ReadRecord] = field(default_factory=list)
    internal_or_offtarget: list[ReadRecord] = field(default_factory=list)

    def counts(self) -> dict[str, int]:
        return {
            "spanning": len(self.spanning),
            "left_anchored": len(self.left_anchored),
            "right_anchored": len(self.right_anchored),
            "internal_or_offtarget": len(self.internal_or_offtarget),
        }


def _mean_q(qual: str) -> float:
    return sum(ord(c) - 33 for c in qual) / len(qual) if qual else 0.0


def _span_in(target: str, anchors: Anchors, k: int) -> tuple[int, int, int, str] | None:
    """Return (start, end_excl, edits, basis) of motif1..motif9 in target, or None."""
    left = infix_hit(anchors.left, target, k)
    right = None
    if left is not None:
        tail = infix_hit(anchors.right, target[left[1] :], k)
        right = None if tail is None else (tail[0] + left[1], tail[1] + left[1], tail[2])
    if left is not None and right is not None:
        return left[0], right[1], left[2] + right[2], "motif"
    kf = max(2, k // 4)
    lf = infix_hit(anchors.left_flank, target, kf)
    rf = infix_hit(anchors.right_flank, target[lf[1] :], kf) if lf else None
    if lf is None or rf is None:
        return None
    return lf[1], rf[0] + lf[1], lf[2] + rf[2], "flank"


def categorize_reads(
    reads: Iterable[ReadRecord], anchors: Anchors, settings: HybridSettings
) -> ReadCategories:
    """Classify reads as spanning (oriented, trimmed) / one-end anchored / other."""
    k = settings.anchor_max_edits
    lo, hi = settings.min_span_units * UNIT, settings.max_span_units * UNIT
    cats = ReadCategories()
    for rec in reads:
        seq = rec.seq.upper()
        best: tuple[tuple[int, int, int, str], str, str, str] | None = None
        for strand, target, qual in (("+", seq, rec.qual), ("-", rc(seq), rec.qual[::-1])):
            hit = _span_in(target, anchors, k)
            if hit and (best is None or hit[2] < best[0][2]):
                best = (hit, strand, target, qual)
        if best is not None:
            (start, end, edits, basis), strand, target, qual = best
            if lo <= end - start <= hi:
                cats.spanning.append(
                    SpanRead(
                        rec.name, target[start:end], _mean_q(qual[start:end]), strand, edits, basis
                    )
                )
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

- [ ] **Step 6: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_spans.py --no-cov -v`
Expected: 5 passed. `test_dupc_template_is_classified` shows that the synthetic dupC is the dictionary template (7C→8C), so classify reports `dupC` with `template_match=True` (pre-flight C3.1).

- [ ] **Step 7: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/hybrid tests/unit/hybrid
git commit -m "feat(hybrid): lazy edlib primitives, synthetic reads and anchor-based categorisation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Length model with scaled windows, smear, noise and rejected peaks (S2)

**Files:**
- Create: `src/muc_one_span/hybrid/lengths.py`
- Test: `tests/unit/hybrid/test_lengths.py`

**Interfaces:**
- Consumes: `SpanRead`, `HybridSettings` (`peak_window_*`, `min_peak_reads`, `far/near_peak_min_frac`, `smear_min_prominence`, `rejected_peak_noise_reads`).
- Produces:
  - `UNIT = 60`
  - `GATE_RELEVANT_REJECTIONS = frozenset({"support_below_threshold", "max_alleles"})`
  - `LengthPeak(center_bp, support, members)`
  - `LengthModel(peaks, rejected, short_products, unassigned, total)`, with properties `unassigned_fraction`, `short_product_fraction` and `gate_relevant_rejections`
  - `window_bp(length_bp: float, settings: HybridSettings) -> float`
  - `fit_length_model(spans: list[SpanRead], settings: HybridSettings) -> LengthModel`
- Each rejected entry is `{"center_bp", "units", "support", "reason"}`, with reason ∈ {smear, noise, support_below_threshold, max_alleles}.
- The smear rule (C4.2, tightened, see Review Focus #2): a candidate more than 1.5 units below the top peak is smear unless it has ≥ `min_peak_reads` reads and inside ≥ `smear_min_prominence` × max(shoulders, 1). Both constants are provisional and tuned on the benchmark dev split.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/hybrid/test_lengths.py
"""S2 length peaks: scaled windows, rejected peaks, smear and close alleles."""

from __future__ import annotations

from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD)


def spans(inner_units: int, n: int, seed: int, **kw: float) -> list[SpanRead]:
    seq = synth.allele(["X"] * inner_units)
    return categorize_reads(synth.reads(seq, n, err=0.02, seed=seed, **kw), ANCH, S).spanning


def test_two_distant_alleles_with_minor_long_allele() -> None:
    model = fit_length_model(spans(30, 200, 1) + spans(70, 12, 2), S)
    assert sorted(round(p.center_bp / 60) for p in model.peaks) == [39, 79]
    assert model.gate_relevant_rejections == []


def test_smear_is_short_product_not_allele() -> None:
    model = fit_length_model(spans(60, 120, 3, smear_frac=0.45), S)
    assert len(model.peaks) == 1
    assert model.short_product_fraction > 0.3
    assert model.gate_relevant_rejections == []


def test_real_short_allele_amid_smear_is_kept() -> None:
    model = fit_length_model(spans(30, 40, 8) + spans(60, 120, 9, smear_frac=0.3), S)
    assert sorted(round(p.center_bp / 60) for p in model.peaks) == [39, 69]


def test_close_alleles_are_not_silently_dropped() -> None:
    # Alleles one unit apart must yield two peaks or a large unassigned fraction.
    model = fit_length_model(spans(40, 150, 4) + spans(41, 150, 5), S)
    assert len(model.peaks) == 2 or model.unassigned_fraction > 0.2


def test_rejected_minor_peak_is_reported() -> None:
    model = fit_length_model(spans(30, 200, 6) + spans(80, 5, 7), S)
    assert len(model.peaks) == 1
    assert [(r["support"], r["reason"]) for r in model.gate_relevant_rejections] == [
        (5, "support_below_threshold")
    ]


def test_third_real_peak_is_max_alleles() -> None:
    model = fit_length_model(spans(30, 60, 10) + spans(50, 60, 11) + spans(70, 60, 12), S)
    assert len(model.peaks) == 2
    assert [r["reason"] for r in model.gate_relevant_rejections] == ["max_alleles"]


def test_single_outlier_read_is_noise_not_gate_relevant() -> None:
    model = fit_length_model(spans(30, 100, 13) + spans(90, 1, 14), S)
    assert [r["reason"] for r in model.rejected] == ["noise"]
    assert model.gate_relevant_rejections == []
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_lengths.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid.lengths'`.

- [ ] **Step 3: Implement (pure Python, no numpy; C4.1)**

```python
# src/muc_one_span/hybrid/lengths.py
"""S2: allele length peaks from spanning-read lengths with length-scaled windows."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

UNIT = 60
GATE_RELEVANT_REJECTIONS = frozenset({"support_below_threshold", "max_alleles"})


@dataclass
class LengthPeak:
    """One accepted allele length peak and its spanning members."""

    center_bp: float
    support: int
    members: list[SpanRead] = field(default_factory=list)


@dataclass
class LengthModel:
    """Accepted peaks plus everything that was not used, so nothing is dropped silently."""

    peaks: list[LengthPeak]
    rejected: list[dict[str, Any]]
    short_products: list[SpanRead]
    unassigned: list[SpanRead]
    total: int

    @property
    def unassigned_fraction(self) -> float:
        return len(self.unassigned) / self.total if self.total else 0.0

    @property
    def short_product_fraction(self) -> float:
        return len(self.short_products) / self.total if self.total else 0.0

    @property
    def gate_relevant_rejections(self) -> list[dict[str, Any]]:
        """Rejected candidates that may be a real allele (not smear or single-read noise)."""
        return [r for r in self.rejected if r["reason"] in GATE_RELEVANT_REJECTIONS]


def window_bp(length_bp: float, settings: HybridSettings) -> float:
    """Assignment half-window: grows with allele length (ONT span noise ~0.25 bp/unit)."""
    return settings.peak_window_base_bp + settings.peak_window_per_unit_bp * length_bp / UNIT


def _density(lengths: list[float]) -> tuple[list[float], list[float]]:
    """Gaussian kernel density on a 2 bp grid; bandwidth grows with length."""
    lo, hi = min(lengths) - 100.0, max(lengths) + 100.0
    grid = [lo + 2.0 * i for i in range(int((hi - lo) / 2.0))]
    dens = [0.0] * len(grid)
    for value in lengths:
        bw = 8.0 + 0.004 * value
        first = max(0, int((value - 4 * bw - lo) / 2.0))
        last = min(len(grid), int((value + 4 * bw - lo) / 2.0) + 1)
        for i in range(first, last):
            dens[i] += math.exp(-0.5 * ((grid[i] - value) / bw) ** 2)
    return grid, dens


def _is_smear(c: float, top: float, lengths: list[float], settings: HybridSettings) -> bool:
    """Below the major peak by >1.5 units, a candidate is a short product (smear) unless it
    has at least ``min_peak_reads`` reads and stands out from its shoulders."""
    if c >= top - 1.5 * UNIT:
        return False
    w = window_bp(c, settings)
    inside = sum(abs(x - c) <= w for x in lengths)
    shoulders = sum(w < abs(x - c) <= 3 * w for x in lengths) / 2  # same width as inside
    return inside < settings.min_peak_reads or inside < settings.smear_min_prominence * max(
        shoulders, 1.0
    )


def _reason(
    c: float, top: float, support: int, n_kept: int, lengths: list[float], settings: HybridSettings
) -> str | None:
    """None when the candidate is accepted, else the rejection reason."""
    if _is_smear(c, top, lengths, settings):
        return "smear"
    if support <= settings.rejected_peak_noise_reads:
        return "noise"
    far = abs(c - top) >= 2 * UNIT
    frac = settings.far_peak_min_frac if far else settings.near_peak_min_frac
    top_support = sum(abs(x - top) <= window_bp(top, settings) for x in lengths)
    if support < settings.min_peak_reads or support < frac * top_support:
        return "support_below_threshold"
    return "max_alleles" if n_kept >= 2 else None


def fit_length_model(spans: list[SpanRead], settings: HybridSettings) -> LengthModel:
    """Pick up to two allele peaks; report rejected peaks, short products, unassigned reads."""
    if not spans:
        return LengthModel([], [], [], [], 0)
    lengths = [float(s.length) for s in spans]
    grid, dens = _density(lengths)
    maxima = [
        i for i in range(1, len(grid) - 1) if dens[i] >= dens[i - 1] and dens[i] > dens[i + 1]
    ]
    maxima.sort(key=lambda i: -dens[i])
    centers: list[float] = []
    for i in maxima:  # keep maxima at least ~0.7 unit apart (Δ1 alleles stay separable)
        if all(abs(grid[i] - c) >= 0.7 * UNIT for c in centers):
            centers.append(grid[i])
    support = {c: sum(abs(x - c) <= window_bp(c, settings) for x in lengths) for c in centers}
    top = max(centers, key=lambda c: support[c])
    kept, rejected = [top], []
    for c in sorted(centers, key=lambda c: -support[c]):
        if c == top:
            continue
        reason = _reason(c, top, support[c], len(kept), lengths, settings)
        if reason is None:
            kept.append(c)
            continue
        rejected.append(
            {
                "center_bp": round(c, 1),
                "units": round(c / UNIT),
                "support": support[c],
                "reason": reason,
            }
        )
    peaks = [LengthPeak(c, support[c]) for c in sorted(kept)]
    short_products: list[SpanRead] = []
    unassigned: list[SpanRead] = []
    shortest = peaks[0].center_bp
    for sp in spans:
        dist = [abs(sp.length - p.center_bp) for p in peaks]
        j = min(range(len(peaks)), key=lambda k: dist[k])
        if dist[j] <= window_bp(peaks[j].center_bp, settings):
            peaks[j].members.append(sp)
        elif sp.length < shortest - 1.5 * UNIT:
            short_products.append(sp)
        else:
            unassigned.append(sp)
    return LengthModel(peaks, rejected, short_products, unassigned, len(spans))
```

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_lengths.py --no-cov -v`
Expected: 7 passed. In validation, the smear test also passed for seeds 21–25. If a test fails, do not loosen its assertion. Inspect `model.rejected` and report the failure.

- [ ] **Step 5: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/hybrid/lengths.py tests/unit/hybrid/test_lengths.py
git commit -m "feat(hybrid): length peaks with scaled windows, smear, noise and rejected-peak reporting

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: POA backend and polishing with partial-read projection (S3, S7)

**Files:**
- Create: `src/muc_one_span/hybrid/poa.py`, `src/muc_one_span/hybrid/polish.py`
- Test: `tests/unit/hybrid/test_polish.py`

**Interfaces:**
- Consumes: `project`, `Columns`, `HYBRID_HINT` (Task 4); `SpanRead`.
- Produces:
  - `poa.PoaBackend` (Protocol with `name` and `consensus(seqs) -> str`)
  - `poa.get_backend(name: str) -> PoaBackend`: `ValueError` for an unknown name, `ImportError` naming the extra for a missing package, no fallback (C5.1)
  - `polish.draft_consensus(members: list[SpanRead], n_poa: int, rng: random.Random, backend: PoaBackend) -> str`
  - `polish.pileup_polish(cons: str, full: list[str], partial: list[str] | None = None) -> tuple[str, int]`
  - `polish._runs(seq: str, min_len: int) -> list[tuple[int, int, str]]`
  - `polish.read_run_length(read: str, t2q: list[int], start: int, end: int, base: str) -> int`
  - `polish.homopolymer_vote(cons: str, full: list[str], partial: list[str] | None = None, min_len: int = 4) -> tuple[str, int]`
  - `polish.polish(cons: str, full: list[str], rounds: int, hp_vote: bool, partial: list[str] | None = None) -> tuple[str, dict[str, Any]]`
- Partial reads (C5.2) align as infixes (`project(..., partial=True)`). Uncovered columns are `None` and do not vote. Insertion slots vote only strictly inside the covered interval. Majority denominators are per-column covered counts. A homopolymer is voted only by reads that cover it plus one column on each side.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/hybrid/test_polish.py
"""S3/S7 POA backend selection and polishing, including partial-read projection."""

from __future__ import annotations

import random
import sys

import pytest

from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import draft_consensus, homopolymer_vote, pileup_polish, polish
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
BACKEND = DEFAULT_SETTINGS.hybrid.poa_backend
DUPC_ALLELE = synth.allele(["X"] * 5 + [synth.dupc()] + ["X"] * 6 + ["A", "B"] + ["X"] * 10)


def _members(seq: str, n: int, seed: int) -> list[SpanRead]:
    return categorize_reads(
        synth.reads(seq, n, err=0.03, seed=seed),
        Anchors.from_dictionary(synth.RD),
        HybridSettings(),
    ).spanning


def test_poa_plus_polish_recovers_exact_allele_with_dupc() -> None:
    pytest.importorskip("pyabpoa", reason="pyabpoa (extra 'hybrid') is not installed")
    members = _members(DUPC_ALLELE, 60, 11)
    cons = draft_consensus(members, 40, random.Random(1), get_backend(BACKEND))
    cons, info = polish(cons, [m.seq for m in members], rounds=2, hp_vote=True)
    assert cons == DUPC_ALLELE, info


def test_pyspoa_backend_is_selectable() -> None:
    pytest.importorskip("spoa", reason="pyspoa (extra 'hybrid') is not installed")
    backend = get_backend("pyspoa")
    assert backend.name == "pyspoa"
    assert backend.consensus(["ACGTACGT", "ACGTACGT", "ACGAACGT"]) == "ACGTACGT"


def test_missing_backend_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "pyabpoa", None)
    with pytest.raises(ImportError, match="hybrid"):
        get_backend("pyabpoa")
    with pytest.raises(ValueError, match="unknown"):
        get_backend("medaka")


def test_partial_reads_do_not_truncate_the_consensus() -> None:
    truth = DUPC_ALLELE
    noisy = truth[:700] + truth[701:1400] + "G" + truth[1400:]  # one deletion, one insertion
    full = [r.seq for r in _members(truth, 20, 12)]
    half = len(truth) // 2
    fragments = [
        r.seq
        for r in synth.reads(truth[:half], 60, err=0.03, seed=13, strand_mix=False, flank_bp=0)
    ]
    cons = noisy
    for _ in range(2):
        cons, _changes = pileup_polish(cons, full, fragments)
        cons, _hp = homopolymer_vote(cons, full, fragments)
    assert cons == truth


def test_homopolymer_vote_uses_median_run_length() -> None:
    cons = "ACGTTGCA" + "C" * 7 + "AGTTGCAT"
    reads = [cons.replace("C" * 7, "C" * 8)] * 6 + [cons] * 3
    new, changes = homopolymer_vote(cons, reads)
    assert new == cons.replace("C" * 7, "C" * 8) and changes == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_polish.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid.poa'`.

- [ ] **Step 3: Implement `poa.py`**

```python
# src/muc_one_span/hybrid/poa.py
"""Partial-order-alignment consensus backends, selected explicitly by setting.

There is no silent fallback: results must be reproducible from the recorded
``hybrid.poa_backend`` setting. The prototype evidence was produced with pyabpoa.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable
from typing import Protocol

from muc_one_span.hybrid.align import HYBRID_HINT


class PoaBackend(Protocol):
    """Minimal consensus interface used by the engine."""

    name: str

    def consensus(self, seqs: list[str]) -> str: ...


class _Abpoa:
    name = "pyabpoa"

    def __init__(self) -> None:
        self._aligner = importlib.import_module("pyabpoa").msa_aligner(aln_mode="g")

    def consensus(self, seqs: list[str]) -> str:
        res = self._aligner.msa(seqs, out_cons=True, out_msa=False)
        return str(res.cons_seq[0])


class _Spoa:
    name = "pyspoa"

    def __init__(self) -> None:
        self._poa = importlib.import_module("spoa").poa

    def consensus(self, seqs: list[str]) -> str:
        cons, _msa = self._poa(seqs, algorithm=1)  # global alignment
        return str(cons)


def get_backend(name: str) -> PoaBackend:
    """Return the named backend; raise ImportError naming the extra when unavailable."""
    factories: dict[str, Callable[[], PoaBackend]] = {"pyabpoa": _Abpoa, "pyspoa": _Spoa}
    if name not in factories:
        raise ValueError(f"unknown POA backend {name!r}")
    try:
        backend = factories[name]()
    except ImportError as exc:
        raise ImportError(f"{HYBRID_HINT} ({name}: {exc})") from exc
    return backend
```

- [ ] **Step 4: Implement `polish.py`**

```python
# src/muc_one_span/hybrid/polish.py
"""S3/S7: POA draft, majority pileup polishing and homopolymer median vote."""

from __future__ import annotations

import random
import statistics
from collections import Counter
from typing import Any

from muc_one_span.hybrid.align import Columns, project
from muc_one_span.hybrid.poa import PoaBackend
from muc_one_span.hybrid.spans import SpanRead


def draft_consensus(
    members: list[SpanRead], n_poa: int, rng: random.Random, backend: PoaBackend
) -> str:
    """POA over a random sample of near-modal members (random, not quality-ranked)."""
    median = statistics.median(m.length for m in members)
    tol = max(15.0, 0.006 * median)
    near = [m for m in members if abs(m.length - median) <= tol] or list(members)
    sample = rng.sample(near, min(n_poa, len(near)))
    return backend.consensus([m.seq for m in sample])


def _projections(cons: str, full: list[str], partial: list[str]) -> list[tuple[str, Columns, bool]]:
    """Spanning reads align globally; fragments align as infixes and cover only part."""
    return [(r, project(r, cons), False) for r in full] + [
        (r, project(r, cons, partial=True), True) for r in partial
    ]


def pileup_polish(cons: str, full: list[str], partial: list[str] | None = None) -> tuple[str, int]:
    """One majority-vote round; each column/insertion slot is voted only by covering reads."""
    n = len(cons)
    col: list[Counter[str | None]] = [Counter() for _ in range(n)]
    ins: list[Counter[str]] = [Counter() for _ in range(n + 1)]
    for _read, proj, is_partial in _projections(cons, full, partial or []):
        for pos in range(proj.first, proj.last):
            col[pos][proj.cols[pos]] += 1
        slots = range(proj.first + 1, proj.last) if is_partial else range(n + 1)
        for pos in slots:
            ins[pos][proj.ins.get(pos, "")] += 1
    out: list[str] = []
    changes = 0
    for pos in range(n + 1):
        if ins[pos]:
            ins_vote, ins_count = ins[pos].most_common(1)[0]
            if ins_vote and ins_count > sum(ins[pos].values()) / 2:
                out.append(ins_vote)
                changes += 1
        if pos < n:
            base = col[pos].most_common(1)[0][0] if col[pos] else cons[pos]
            if base != "-":
                out.append(str(base))
            changes += int(base != cons[pos])
    return "".join(out), changes


def _runs(seq: str, min_len: int) -> list[tuple[int, int, str]]:
    """Homopolymer runs of at least ``min_len`` as (start, end_excl, base)."""
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
    """Longest stretch of ``base`` in the read around the interval mapped to [start, end)."""
    qs, qe = t2q[start], t2q[end]
    a = qs
    while a > 0 and read[a - 1] == base:
        a -= 1
    b = max(qe, qs)
    while b < len(read) and read[b] == base:
        b += 1
    seg = "".join(c if c == base else " " for c in read[a:b])
    return max((len(m) for m in seg.split()), default=0)


def homopolymer_vote(
    cons: str, full: list[str], partial: list[str] | None = None, min_len: int = 4
) -> tuple[str, int]:
    """Set each run (>= min_len) to the median run length of the reads covering it."""
    runs = _runs(cons, min_len)
    if not runs or not (full or partial):
        return cons, 0
    projections = _projections(cons, full, partial or [])
    out, pos, changes = [], 0, 0
    for start, end, base in runs:
        lengths = [
            read_run_length(r, p.t2q, start, end, base)
            for r, p, _ in projections
            if p.covers(start, end)
        ]
        new_len = max(1, int(statistics.median(lengths) + 0.5)) if lengths else end - start
        out.append(cons[pos:start])
        out.append(base * new_len)
        changes += int(new_len != end - start)
        pos = end
    out.append(cons[pos:])
    return "".join(out), changes


def polish(
    cons: str,
    full: list[str],
    rounds: int,
    hp_vote: bool,
    partial: list[str] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Run ``rounds`` pileup rounds, each followed by an optional homopolymer vote."""
    info: dict[str, Any] = {"rounds": []}
    for _ in range(rounds):
        cons, changes = pileup_polish(cons, full, partial)
        hp = 0
        if hp_vote:
            cons, hp = homopolymer_vote(cons, full, partial)
        info["rounds"].append({"changes": changes, "hp_changes": hp})
    return cons, info
```

- [ ] **Step 5: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_polish.py --no-cov -v`
Expected: 5 passed. The first test proves that pyabpoa followed by polishing reproduces a dictionary-dupC allele exactly. `test_partial_reads_do_not_truncate_the_consensus` covers 20 spanning reads plus 60 left-half fragments, and must reproduce the full-length truth.

- [ ] **Step 6: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/hybrid/poa.py src/muc_one_span/hybrid/polish.py tests/unit/hybrid/test_polish.py
git commit -m "feat(hybrid): explicit POA backend and coverage-aware pileup/homopolymer polishing

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: All-read assignment with an off-target guard (S5, S6)

**Files:**
- Create: `src/muc_one_span/hybrid/assign.py`
- Test: `tests/unit/hybrid/test_assign.py`

**Interfaces:**
- Consumes: `edit_distance_infix`, `project`, `rc`, `ReadRecord`, `RepeatDictionary`, `HybridSettings` (`min_fragment_bp`, `assign_margin`, `assign_max_error_rate`).
- Produces:
  - `OFF_TARGET = "off_target"`, `UNDECIDED = "undecided"`, `FLANK_BP = 500`
  - `hybrid_references(drafts: dict[str, str], rd: RepeatDictionary, flank_bp: int = FLANK_BP) -> dict[str, str]`
  - `Assignment(allele: str | None, margin: int, distances: dict[str, int], oriented: str)`
  - `assign_read(seq: str, refs: dict[str, str], margin: int, max_error_rate: float = 1.0) -> Assignment`
  - `assign_reads(reads: list[ReadRecord], refs: dict[str, str], settings: HybridSettings) -> dict[str, list[str]]`, with keys = allele names + `undecided` + `off_target`
  - `trim_to_draft(oriented: str, ref: str, flank_bp: int, draft_len: int) -> str`
- Homozygous behaviour, stated once: with a single reference every read of at least `min_fragment_bp` is assigned to it, subject to the same off-target guard (C7.2).

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/hybrid/test_assign.py
"""S5/S6 edit-distance competition against ladder-flanked allele drafts."""

from __future__ import annotations

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.assign import assign_read, assign_reads, hybrid_references
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
A = ["X"] * 10 + ["A", "A", "B"] + ["X"] * 20 + ["G", "A", "B"] + ["X"] * 10
B = ["X"] * 10 + ["A", "A", "B"] + ["X"] * 21 + ["G", "A", "B"] + ["X"] * 10  # one unit longer
REFS = hybrid_references({"allele_1": synth.allele(A), "allele_2": synth.allele(B)}, synth.RD)
TRACT_START = (5 + 13) * 60  # pre-repeats 5, then 10 X, then A A B
TRACT_END_B = (5 + 13 + 21) * 60  # the 21-X tract in B


def test_fragments_spanning_the_difference_are_assigned_correctly() -> None:
    wrong = decided = 0
    for truth, inner in (("allele_1", A), ("allele_2", B)):
        seq = synth.allele(inner)
        for i, start in enumerate(range(TRACT_START - 400, TRACT_START - 100, 20)):
            piece = seq[start : TRACT_END_B + 400]
            frag = synth.reads(piece, 1, err=0.02, seed=100 + i, flank_bp=0)[0].seq
            got = assign_read(rc(frag) if i % 2 else frag, REFS, margin=3, max_error_rate=0.15)
            decided += got.allele is not None
            wrong += got.allele not in (None, truth)
    assert wrong == 0 and decided == 30


def test_fragment_inside_shared_tract_is_undecided() -> None:
    frag = synth.allele(A)[TRACT_START : TRACT_START + 15 * 60]
    assert assign_read(frag, REFS, margin=3).allele is None


def test_off_target_read_is_not_assigned() -> None:
    reads = [ReadRecord("junk", "ACGT" * 400, "5" * 1600)]
    out = assign_reads(reads, REFS, HybridSettings())
    assert len(out["off_target"]) == 1 and not out["allele_1"] and not out["allele_2"]


def test_single_reference_assigns_on_target_and_guards_off_target() -> None:
    single = {"allele_1": REFS["allele_1"]}
    frag = synth.allele(A)[600:2400]
    assert assign_read(frag, single, margin=3, max_error_rate=0.15).allele == "allele_1"
    assert assign_read("ACGT" * 400, single, margin=3, max_error_rate=0.15).allele == "off_target"


def test_short_reads_are_ignored() -> None:
    out = assign_reads([ReadRecord("s", synth.allele(A)[:500], "5" * 500)], REFS, HybridSettings())
    assert sum(len(v) for v in out.values()) == 0


def test_trim_to_draft_drops_ladder_flanks() -> None:
    from muc_one_span.hybrid.assign import FLANK_BP, trim_to_draft

    draft = synth.allele(A)
    ref = REFS["allele_1"]
    read = ref[FLANK_BP - 200 : FLANK_BP + 1200]  # 200 bp flank + 1200 bp of the draft
    assert trim_to_draft(read, ref, FLANK_BP, len(draft)) == draft[:1200]
    assert trim_to_draft(ref[:400], ref, FLANK_BP, len(draft)) == ""
```

Tract offsets (C7.1, verified): pre-repeats 5, then 10 X, then `A A B`, so the X tract starts at unit 18 (1 080 bp). It holds 20 X in `A` and 21 X in `B`. Every fragment spans the whole tract plus unique `AAB`/`GAB` context on both sides, so every fragment is decidable.

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_assign.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid.assign'`.

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/hybrid/assign.py
"""S5/S6: ladder-flanked allele references and edit-distance competition for every read."""

from __future__ import annotations

from dataclasses import dataclass

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import edit_distance_infix, project, rc
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.settings import HybridSettings

OFF_TARGET = "off_target"
UNDECIDED = "undecided"
FLANK_BP = 500


def hybrid_references(
    drafts: dict[str, str], rd: RepeatDictionary, flank_bp: int = FLANK_BP
) -> dict[str, str]:
    """Wrap each motif1..motif9 draft in the ladder's hg38 flanks."""
    left, right = rd.flanking_left[-flank_bp:], rd.flanking_right[:flank_bp]
    return {name: left + seq + right for name, seq in drafts.items()}


@dataclass(frozen=True)
class Assignment:
    """Outcome for one read: an allele name, ``undecided``, or ``off_target``."""

    allele: str | None
    margin: int
    distances: dict[str, int]
    oriented: str


def assign_read(
    seq: str, refs: dict[str, str], margin: int, max_error_rate: float = 1.0
) -> Assignment:
    """Assign seq to the reference with the smallest infix edit distance, if clearly best.

    Returns ``allele=None`` when the margin is too small, and ``allele="off_target"`` when
    even the best reference needs more than ``max_error_rate * len(seq)`` edits. A single
    reference (homozygous sample) is assigned subject to the same off-target guard.
    """
    best: Assignment | None = None
    for oriented in (seq, rc(seq)):
        dist = {name: edit_distance_infix(oriented, ref) for name, ref in refs.items()}
        ranked = sorted(dist.values())
        gap = ranked[1] - ranked[0] if len(ranked) > 1 else len(seq)
        winner = min(dist, key=lambda k: dist[k])
        cand = Assignment(winner if gap >= margin else None, gap, dist, oriented)
        if best is None or ranked[0] < min(best.distances.values()):
            best = cand
    if best is None:
        raise ValueError("no reference to assign against")
    if min(best.distances.values()) > max_error_rate * len(seq):
        return Assignment(OFF_TARGET, best.margin, best.distances, best.oriented)
    return best


def assign_reads(
    reads: list[ReadRecord], refs: dict[str, str], settings: HybridSettings
) -> dict[str, list[str]]:
    """Map allele -> oriented read sequences, plus ``undecided`` and ``off_target``.

    Reads shorter than ``min_fragment_bp`` are not considered.
    """
    out: dict[str, list[str]] = {name: [] for name in refs}
    out[UNDECIDED] = []
    out[OFF_TARGET] = []
    for rec in reads:
        if len(rec.seq) < settings.min_fragment_bp:
            continue
        got = assign_read(
            rec.seq.upper(), refs, settings.assign_margin, settings.assign_max_error_rate
        )
        out[got.allele or UNDECIDED].append(got.oriented)
    return out


def trim_to_draft(oriented: str, ref: str, flank_bp: int, draft_len: int) -> str:
    """Cut an assigned read to the part aligned inside the draft (drop ladder flanks)."""
    proj = project(oriented, ref, partial=True)
    lo, hi = max(proj.first, flank_bp), min(proj.last, flank_bp + draft_len)
    if hi - lo <= 0:
        return ""
    return oriented[proj.t2q[lo] : proj.t2q[hi]]
```

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_assign.py --no-cov -v`
Expected: 6 passed (30/30 decidable fragments assigned, 0 wrong).

- [ ] **Step 5: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/hybrid/assign.py tests/unit/hybrid/test_assign.py
git commit -m "feat(hybrid): edit-distance read assignment with off-target guard and flank trimming

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Linked-site phase split for equal and close lengths (S4)

**Files:**
- Create: `src/muc_one_span/hybrid/phase.py`
- Test: `tests/unit/hybrid/test_phase.py`

**Interfaces:**
- Consumes: `global_columns`, `_runs`, `read_run_length`, `SpanRead`, `HybridSettings` (`het_af_min`, `link_phi_min`, `min_linked_sites`, `het_min_group`).
- Produces:
  - `PhaseResult(groups: list[list[SpanRead]], basis: str, sites: list[dict], candidate: dict | None)`, with basis ∈ {none, linked_sites, unconfirmed_single_site}
  - `split_by_linked_sites(cons: str, members: list[SpanRead], settings: HybridSettings, rng: random.Random) -> PhaseResult`
- Port of `hetsplit.py:87-251` (the prototype file has 251 lines). The EM `phase_reads` (164-203) is **not** ported; majority vote over linked sites is used instead (recorded deviation). The run-site threshold is `max(het_af_min, 4·bg)` (D6).

- [ ] **Step 1: Write the failing tests (C6.1 with unit Q, D5)**

```python
# tests/unit/hybrid/test_phase.py
"""S4 linked-site phase split for equal and close allele lengths."""

from __future__ import annotations

import random

import pytest

from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
S = HybridSettings()
ANCH = Anchors.from_dictionary(synth.RD)


def _spans(inner: list[str], n: int, seed: int) -> list[SpanRead]:
    reads = synth.reads(synth.allele(inner), n, err=0.02, seed=seed)
    return categorize_reads(reads, ANCH, S).spanning


def test_equal_length_alleles_with_two_linked_differences_split() -> None:
    a = ["X"] * 10 + ["A"] + ["X"] * 10 + ["B"] + ["X"] * 8
    b = ["X"] * 30
    members = _spans(a, 40, 1) + _spans(b, 40, 2)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "linked_sites" and sorted(len(g) for g in res.groups) == [40, 40]


def test_identical_alleles_do_not_split() -> None:
    members = _spans(["X"] * 30, 80, 3)
    res = split_by_linked_sites(synth.allele(["X"] * 30), members, S, random.Random(1))
    assert res.basis == "none" and len(res.groups) == 1


def test_single_base_unit_difference_is_unconfirmed() -> None:
    a = ["X"] * 10 + ["Q"] + ["X"] * 19  # Q differs from X at 1 base, outside any run
    members = _spans(a, 40, 1) + _spans(["X"] * 30, 40, 2)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "unconfirmed_single_site" and len(res.groups) == 1
    assert res.candidate is not None


def test_multi_base_unit_difference_counts_as_linked_sites() -> None:
    a = ["X"] * 10 + ["A"] + ["X"] * 19  # A differs from X at 4 linked bases (intended)
    members = _spans(a, 40, 1) + _spans(["X"] * 30, 40, 2)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "linked_sites"
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_phase.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid.phase'`.

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/hybrid/phase.py
"""S4: split a length peak only when >=2 read-linked difference sites support it.

Ported from the prototype ``hetsplit.py`` (site table, candidate sites, linkage).
Deviation: groups are formed by a majority vote over linked sites; the prototype's EM
read phasing is not ported. Sites are counted per consensus column, so one unit type
that differs from its neighbour at several bases yields several perfectly linked sites.
"""

from __future__ import annotations

import math
import random
import statistics
from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import HybridSettings

MAX_SITE_READS = 300
Site = tuple[str, int]


@dataclass
class PhaseResult:
    """Groups of members and the basis of the (non-)split."""

    groups: list[list[SpanRead]]
    basis: str  # "none", "linked_sites" or "unconfirmed_single_site"
    sites: list[dict[str, Any]] = field(default_factory=list)
    candidate: dict[str, Any] | None = None


def _features(cons: str, reads: list[str]) -> tuple[list[dict[Site, Any]], dict[Site, Any]]:
    runs = _runs(cons, 3)
    in_run = {i for s, e, _ in runs for i in range(s, e)}
    feats = []
    for read in reads:
        proj = global_columns(read, cons)
        f: dict[Site, Any] = {}
        for pos, base in enumerate(proj.cols):
            if pos not in in_run:
                f[("col", pos)] = base
                if pos - 1 not in in_run:
                    f[("ins", pos)] = proj.ins.get(pos, "")
        for s, e, b in runs:
            f[("run", s)] = read_run_length(read, proj.t2q, s, e, b)
        feats.append(f)
    return feats, {("run", s): (b, e - s) for s, e, b in runs}


def _candidates(
    feats: list[dict[Site, Any]], meta: dict[Site, Any], af_min: float
) -> list[dict[str, Any]]:
    counts: dict[Site, Counter[Any]] = {}
    for f in feats:
        for site, allele in f.items():
            counts.setdefault(site, Counter())[allele] += 1
    fracs: dict[tuple[str, int, int], list[float]] = {}
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
            _score, minor, n_minor = max(
                (n / tot - 4 * bg.get((base, length, a), 0.0), a, n) for a, n in minors
            )
            if n_minor / tot < max(af_min, 4 * bg.get((base, length, minor), 0.0)):
                continue
        else:
            minor, n_minor = minors[0]
            af = n_minor / tot
            if af < af_min or ("-" in (major, minor) and af < 1.5 * af_min):
                continue
        out.append(
            {"site": site, "major": major, "minor": minor, "af": round(n_minor / tot, 3), "n": tot}
        )
    return out


def _phi(pairs: list[tuple[int, int]]) -> float:
    a = sum(1 for x, y in pairs if x and y)
    b = sum(1 for x, y in pairs if x and not y)
    c = sum(1 for x, y in pairs if not x and y)
    d = len(pairs) - a - b - c
    den = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    return abs(a * d - b * c) / den if den else 0.0


def _linked(
    feats: list[dict[Site, Any]], sites: list[dict[str, Any]], min_phi: float
) -> list[dict[str, Any]]:
    """Largest connected component of sites whose read indicators correlate (|phi|)."""

    def ind(f: dict[Site, Any], s: dict[str, Any]) -> int | None:
        a = f.get(s["site"])
        return 1 if a == s["minor"] else 0 if a == s["major"] else None

    vec = [[ind(f, s) for s in sites] for f in feats]
    n = len(sites)
    adj: dict[int, set[int]] = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            pairs = [(v[i], v[j]) for v in vec if v[i] is not None and v[j] is not None]
            if len(pairs) >= 10 and _phi(pairs) >= min_phi:  # type: ignore[arg-type]
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


def split_by_linked_sites(
    cons: str, members: list[SpanRead], settings: HybridSettings, rng: random.Random
) -> PhaseResult:
    """Return one group (no split) unless >= min_linked_sites linked sites support two."""
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
    for m, f in zip(members, all_feats, strict=True):
        votes = [
            1 if f.get(s["site"]) == s["minor"] else 0 if f.get(s["site"]) == s["major"] else None
            for s in linked
        ]
        known = [v for v in votes if v is not None]
        if known:
            groups[int(sum(known) * 2 >= len(known) + 1)].append(m)
    if min(len(g) for g in groups) < settings.het_min_group * len(members):
        return PhaseResult([members], "unconfirmed_single_site", linked, candidate=linked[0])
    return PhaseResult(groups, "linked_sites", linked)
```

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_phase.py --no-cov -v`
Expected: 4 passed (validated with seed pairs 1/2, 5/6, 7/8 and 9/10).

- [ ] **Step 5: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/hybrid/phase.py tests/unit/hybrid/test_phase.py
git commit -m "feat(hybrid): linked-site phase split with unconfirmed single-site flag

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Residual QC and per-event read evidence (S8, S10)

**Files:**
- Create: `src/muc_one_span/hybrid/evidence.py`
- Test: `tests/unit/hybrid/test_evidence.py`

**Interfaces:**
- Consumes:
  - `global_columns`, `edit_distance`, `Columns`, `_runs`, `read_run_length`;
  - `RepeatDictionary` (`rd.mutations[name]["changes"]`, `rd.repeats`);
  - the `classify_sequence` output (`repeats[].index/start/end`, `mutations_detected[].repeat_index/closest_type/mutation_name`);
  - `HybridSettings` (`hp_*`, `qc_residual_af`).
- Produces:
  - `residual_sites(cons: str, reads: list[str], af: float) -> list[dict[str, Any]]`
  - `homopolymer_llr(obs: list[tuple[str, int]], background: dict[str, list[float]], shift: int = 1) -> float`
  - `homopolymer_event_run(mutation: dict, rd: RepeatDictionary, cons: str, start: int, end: int) -> tuple[int, int, str, int] | None` (C8.1: typed from the template)
  - `hp_status(n, llr, alt_frac, strand_llr, strand_n, s) -> str` (C8.2)
  - `event_read_support(cons: str, classification: dict, reads: list[tuple[str, str]], rd: RepeatDictionary, settings: HybridSettings) -> dict[int, dict[str, Any]]`
- Each `read_support` is `{kind: homopolymer|competition|none, n, alt, ref, other, alt_frac, strand_alt_frac, llr?, strand_llr?, status}` with status ∈ `READ_SUPPORT_STATUSES` (Task 3). The competition rule is C8.3.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/hybrid/test_evidence.py
"""S8/S10 residual QC and read-level event support (homopolymer LLR and competition)."""

from __future__ import annotations

import pytest

from muc_one_span.classify import classify_sequence
from muc_one_span.hybrid.evidence import (
    event_read_support,
    homopolymer_event_run,
    homopolymer_llr,
    hp_status,
    residual_sites,
)
from muc_one_span.hybrid.spans import Anchors, categorize_reads
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
S = HybridSettings()
BACKGROUND = {
    "+": [0.02] * 6 + [0.20, 0.55, 0.20] + [0.03 / 8] * 8,
    "-": [0.01] * 6 + [0.05, 0.88, 0.05] + [0.01 / 8] * 8,
}


def _oriented(allele: str, n: int, seed: int) -> list[tuple[str, str]]:
    cats = categorize_reads(
        synth.reads(allele, n, err=0.02, seed=seed), Anchors.from_dictionary(synth.RD), S
    )
    return [(sp.seq, sp.strand) for sp in cats.spanning]


def test_residual_sites_flag_mixture_but_not_clean_reads() -> None:
    cons = "ACGTACGTTTGACCATGCA" * 10
    alt = cons[:50] + "G" + cons[51:]
    assert residual_sites(cons, [cons] * 20, af=0.25) == []
    sites = residual_sites(cons, [cons] * 12 + [alt] * 8, af=0.25)
    assert [s["pos"] for s in sites] == [50]


def test_homopolymer_llr_separates_8c_from_7c_background() -> None:
    obs_mut = [("+", 8)] * 12 + [("+", 7)] * 8 + [("-", 8)] * 18 + [("-", 7)] * 2
    obs_wt = [("+", 7)] * 12 + [("+", 6)] * 8 + [("-", 7)] * 18 + [("-", 8)] * 2
    assert homopolymer_llr(obs_mut, BACKGROUND) > 10
    assert homopolymer_llr(obs_wt, BACKGROUND) < 0


def test_single_strand_input_is_supported_not_failed() -> None:
    obs = [("-", 8)] * 25
    llr = homopolymer_llr(obs, BACKGROUND)
    strand = {"+": homopolymer_llr([], BACKGROUND), "-": llr}
    assert hp_status(25, llr, 1.0, strand, {"+": 0, "-": 25}, S) == "supported"


def test_negative_well_covered_strand_is_discordant() -> None:
    plus = [("+", 7)] * 6
    minus = [("-", 8)] * 24
    total = homopolymer_llr(plus + minus, BACKGROUND)
    strand = {"+": homopolymer_llr(plus, BACKGROUND), "-": homopolymer_llr(minus, BACKGROUND)}
    assert total > S.hp_llr_min and strand["+"] < 0
    assert hp_status(30, total, 0.8, strand, {"+": 6, "-": 24}, S) == "discordant"


def test_low_read_count_is_insufficient_depth() -> None:
    assert hp_status(19, 50.0, 1.0, {"+": 5.0, "-": 5.0}, {"+": 10, "-": 9}, S) == (
        "insufficient_depth"
    )


def test_event_type_comes_from_the_dictionary_template() -> None:
    cons = synth.allele(["X"] * 3 + [synth.dupc()] + ["X"] * 3)
    start = (5 + 3) * 60
    run = homopolymer_event_run({"mutation_name": "dupC"}, synth.RD, cons, start, start + 61)
    assert run == (start + 52, start + 60, "C", 1)
    for name in ("dupA", "insG_pos58", "del18_31", "insCCCC"):
        assert (
            homopolymer_event_run({"mutation_name": name}, synth.RD, cons, start, start + 61)
            is None
        )


def _support(allele: str, reads: list[tuple[str, str]], name: str) -> dict:
    classification = classify_sequence(allele, synth.RD)
    idx = next(
        i for i, m in enumerate(classification["mutations_detected"]) if m["mutation_name"] == name
    )
    return event_read_support(allele, classification, reads, synth.RD, S)[idx]


def test_dupc_carrier_reads_support_the_event() -> None:
    allele = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    got = _support(allele, _oriented(allele, 40, 21), "dupC")
    assert got["kind"] == "homopolymer" and got["status"] == "supported", got


def test_wild_type_reads_do_not_support_a_consensus_dupc() -> None:
    allele = synth.allele(["X"] * 10 + [synth.dupc()] + ["X"] * 10)
    wild = _oriented(synth.allele(["X"] * 21), 40, 22)
    assert _support(allele, wild, "dupC")["status"] == "not_supported"


def test_non_homopolymer_event_uses_competition() -> None:
    dupa = next(s for s, (p, n) in synth.RD.mutated_sequences.items() if p == "X" and n == "dupA")
    allele = synth.allele(["X"] * 10 + [dupa] + ["X"] * 10)
    got = _support(allele, _oriented(allele, 40, 23), "dupA")
    assert got["kind"] == "competition" and got["status"] == "supported", got
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_evidence.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid.evidence'`.

- [ ] **Step 3: Implement**

```python
# src/muc_one_span/hybrid/evidence.py
"""S8/S10: residual heterogeneity QC and per-event read-level support.

Producer contract for ``read_support.status == "supported"`` (consumed unchanged by
``clinical_gates.mutation_supported``):

* Homopolymer event (dictionary template = single-base indel inside a run >= 4):
  ``n >= hp_min_reads``, stutter-aware LLR >= ``hp_llr_min``, alt fraction >=
  ``hp_min_alt_frac``, and no strand with >= ``hp_min_strand_reads`` reads has a
  negative strand LLR (else ``discordant``). A strand with zero reads never fails.
* Other templated events: parent-vs-template edit-distance competition per read;
  ``n >= hp_min_reads``, ``alt / n >= hp_min_alt_frac`` and ``alt > ref``.
* ``n < hp_min_reads`` is ``insufficient_depth``: blocked, never negative.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.align import Columns, edit_distance, global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.settings import HybridSettings

MAXLEN = 16
Profile = dict[str, list[float]]


def residual_sites(cons: str, reads: list[str], af: float) -> list[dict[str, Any]]:
    """Consensus columns where a non-consensus base/deletion reaches ``af`` (outside runs >= 3)."""
    in_run = {i for s, e, _ in _runs(cons, 3) for i in range(s, e)}
    counts: list[Counter[str | None]] = [Counter() for _ in cons]
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
                out.append(
                    {"pos": pos, "ref": cons[pos], "alt": base, "af": round(n / tot, 3), "n": tot}
                )
                break
    return out


def _smooth(counter: Counter[int]) -> list[float]:
    tot = sum(counter.values())
    return [(counter.get(i, 0) + 0.5) / (tot + 0.5 * (MAXLEN + 1)) for i in range(MAXLEN + 1)]


def _shift(p: list[float], d: int) -> list[float]:
    q = [p[min(max(i - d, 0), MAXLEN)] for i in range(MAXLEN + 1)]
    s = sum(q)
    return [x / s for x in q]


def homopolymer_background(
    cons: str, reads: list[tuple[str, str, Columns]], base: str, length: int, exclude: int
) -> Profile:
    """Per-strand observed-length profile at consensus runs of ``base`` x ``length``."""
    runs = [
        (s, e)
        for s, e, b in _runs(cons, length)
        if b == base and e - s == length and not s <= exclude < e
    ]
    per: dict[str, Counter[int]] = {"+": Counter(), "-": Counter()}
    for seq, strand, proj in reads:
        for s, e in runs:
            per[strand][min(read_run_length(seq, proj.t2q, s, e, base), MAXLEN)] += 1
    return {s: _smooth(c) for s, c in per.items()}


def homopolymer_llr(obs: list[tuple[str, int]], background: Profile, shift: int = 1) -> float:
    """log L(event length) - log L(reference length); the event profile is the shifted one."""
    llr = 0.0
    for strand, observed in obs:
        p0 = background.get(strand) or [1 / (MAXLEN + 1)] * (MAXLEN + 1)
        p1 = _shift(p0, shift)
        k = min(observed, MAXLEN)
        llr += math.log(p1[k] / p0[k])
    return llr


def homopolymer_event_run(
    mutation: dict[str, Any], rd: RepeatDictionary, cons: str, start: int, end: int
) -> tuple[int, int, str, int] | None:
    """(run_start, run_end, base, shift) for a templated single-base indel in a run >= 4."""
    template = rd.mutations.get(mutation.get("mutation_name") or "") or {}
    changes = template.get("changes", [])
    if len(changes) != 1 or changes[0].get("type") not in ("insert", "delete"):
        return None
    change = changes[0]
    if change["type"] == "insert" and len(change.get("sequence", "")) != 1:
        return None
    if change["type"] == "delete" and int(change["end"]) != int(change["start"]):
        return None
    pos = start + int(change["start"]) - 1
    base = change.get("sequence") or cons[max(pos - 1, start)]
    for s, e, b in _runs(cons[start:end], 4):
        if b == base and s + start <= pos <= e + start:
            return s + start, e + start, b, 1 if change["type"] == "insert" else -1
    return None


def hp_status(
    n: int,
    llr: float,
    alt_frac: float,
    strand_llr: dict[str, float],
    strand_n: dict[str, int],
    s: HybridSettings,
) -> str:
    """Status of a homopolymer event (spec §5 thresholds plus strand consistency)."""
    if n < s.hp_min_reads:
        return "insufficient_depth"
    if llr < s.hp_llr_min or alt_frac < s.hp_min_alt_frac:
        return "not_supported"
    if any(strand_n[st] >= s.hp_min_strand_reads and strand_llr[st] < 0 for st in "+-"):
        return "discordant"
    return "supported"


def _homopolymer_support(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    run: tuple[int, int, str, int],
    s: HybridSettings,
) -> dict[str, Any]:
    start, end, base, shift = run
    length = end - start
    background = homopolymer_background(cons, reads, base, length - shift, exclude=start)
    obs = [(st, read_run_length(seq, p.t2q, start, end, base)) for seq, st, p in reads]
    is_alt = (lambda k: k >= length) if shift > 0 else (lambda k: k <= length)
    alt = sum(1 for _, k in obs if is_alt(k))
    n = len(obs)
    strand_obs = {st: [o for o in obs if o[0] == st] for st in "+-"}
    strand_llr = {st: homopolymer_llr(v, background, shift) for st, v in strand_obs.items()}
    strand_n = {st: len(v) for st, v in strand_obs.items()}
    llr = homopolymer_llr(obs, background, shift)
    alt_frac = alt / n if n else 0.0
    return {
        "kind": "homopolymer",
        "n": n,
        "alt": alt,
        "ref": n - alt,
        "other": 0,
        "alt_frac": round(alt_frac, 3),
        "strand_alt_frac": {
            st: round(sum(1 for _, k in v if is_alt(k)) / len(v), 3) if v else None
            for st, v in strand_obs.items()
        },
        "llr": round(llr, 1),
        "strand_llr": {k: round(v, 1) for k, v in strand_llr.items()},
        "status": hp_status(n, llr, alt_frac, strand_llr, strand_n, s),
    }


def _competition_support(
    cons: str,
    reads: list[tuple[str, str, Columns]],
    start: int,
    end: int,
    parent: str,
    s: HybridSettings,
) -> dict[str, Any]:
    window = cons[start:end]
    alt = ref = other = 0
    per_strand: dict[str, list[int]] = {"+": [], "-": []}
    for seq, strand, proj in reads:
        piece = seq[proj.t2q[start] : proj.t2q[end]]
        d_alt, d_ref = edit_distance(piece, window), edit_distance(piece, parent)
        alt += d_alt < d_ref
        ref += d_alt > d_ref
        other += d_alt == d_ref
        per_strand[strand].append(int(d_alt < d_ref))
    n = alt + ref + other
    frac = alt / n if n else 0.0
    status = (
        "insufficient_depth"
        if n < s.hp_min_reads
        else "supported"
        if frac >= s.hp_min_alt_frac and alt > ref
        else "not_supported"
    )
    return {
        "kind": "competition",
        "n": n,
        "alt": alt,
        "ref": ref,
        "other": other,
        "alt_frac": round(frac, 3),
        "strand_alt_frac": {
            st: round(sum(v) / len(v), 3) if v else None for st, v in per_strand.items()
        },
        "status": status,
    }


def event_read_support(
    cons: str,
    classification: dict[str, Any],
    reads: list[tuple[str, str]],
    rd: RepeatDictionary,
    settings: HybridSettings,
) -> dict[int, dict[str, Any]]:
    """Read-level support per detected mutation (keyed by its index in the classification).

    ``reads`` are (sequence, strand) pairs of spanning reads assigned to this allele,
    already oriented to the consensus.
    """
    repeats = {r.get("index"): r for r in classification.get("repeats", [])}
    projected = [(seq, strand, global_columns(seq, cons)) for seq, strand in reads]
    results: dict[int, dict[str, Any]] = {}
    for idx, mut in enumerate(classification.get("mutations_detected", [])):
        rep = repeats.get(mut.get("repeat_index"))
        parent = rd.repeats.get(mut.get("closest_type") or "")
        if rep is None or parent is None:
            results[idx] = {"kind": "none", "status": "not_localized"}
            continue
        start, end = int(rep["start"]), int(rep["end"])
        run = homopolymer_event_run(mut, rd, cons, start, end)
        results[idx] = (
            _homopolymer_support(cons, projected, run, settings)
            if run is not None
            else _competition_support(cons, projected, start, end, parent, settings)
        )
    return results
```

Do not port `scipy.stats.binom` from `hp_model.py`; the LLR does not need it.

- [ ] **Step 4: Run the tests**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_evidence.py --no-cov -v`
Expected: 9 passed.
- dupC is typed `homopolymer`; `dupA`, `insG_pos58`, `del18_31` and `insCCCC` are not.
- dupA support is `competition` / `supported`.
- Wild-type reads against a dupC consensus give `not_supported`.

- [ ] **Step 5: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/hybrid/evidence.py tests/unit/hybrid/test_evidence.py
git commit -m "feat(hybrid): residual QC and template-typed read-level event support

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 10: Move `run` to `cli_run.py`; `--engine/--assay`; pipeline kwargs; temporary guard

**Files:**
- Create: `src/muc_one_span/cli_run.py`, `tests/unit/test_cli_engine.py`
- Modify: `src/muc_one_span/cli.py`, `src/muc_one_span/pipeline.py`

**Interfaces:**
- Produces:
  - `muc_one_span.cli_run.run` (Click command), still importable as `muc_one_span.cli.run`;
  - `execute_pipeline(..., *, engine: str | None = None, assay: str | None = None)`, forwarded into `effective_run_settings` only when not None;
  - a temporary guard: `settings.run.engine == "hybrid"` raises `click.BadParameter` (exit code 2). Task 11 removes it (pre-flight (a) Task 2 item 6).
- `--config` defaults for `run.engine`/`run.assay` work through `configure_context`, which maps `run.*` fields to same-named Click parameters.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/test_cli_engine.py
"""The run command's --engine/--assay options and their configuration defaults."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_span.cli import main


def _fastq(tmp_path: Path) -> Path:
    fq = tmp_path / "r.fastq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    return fq


def test_engine_option_reaches_pipeline(tmp_path: Path) -> None:
    args = ["run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path / "o")]
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, [*args, "--engine", "hybrid", "--assay", "genomic"])
    assert res.exit_code == 0, res.output
    assert (ex.call_args.kwargs["engine"], ex.call_args.kwargs["assay"]) == ("hybrid", "genomic")


def test_engine_defaults_to_ladder(tmp_path: Path) -> None:
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, ["run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path)])
    assert res.exit_code == 0, res.output
    assert (ex.call_args.kwargs["engine"], ex.call_args.kwargs["assay"]) == ("ladder", "amplicon")


def test_engine_default_comes_from_config(tmp_path: Path) -> None:
    cfg = tmp_path / "c.json"
    cfg.write_text('{"schema_version": 1, "run": {"engine": "hybrid"}}')
    args = ["--config", str(cfg), "run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path / "o")]
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, args)
    assert res.exit_code == 0, res.output
    assert ex.call_args.kwargs["engine"] == "hybrid"


def test_run_is_still_importable_from_cli() -> None:
    from muc_one_span.cli import run
    from muc_one_span.cli_run import run as moved

    assert run is moved and main.commands["run"] is moved


def test_hybrid_engine_is_rejected_until_available(tmp_path: Path) -> None:
    """Temporary (removed in Task 11): never silently run the ladder as 'hybrid'."""
    res = CliRunner().invoke(
        main, ["run", "-i", str(_fastq(tmp_path)), "-o", str(tmp_path / "o"), "--engine", "hybrid"]
    )
    assert res.exit_code == 2, res.output
    assert "hybrid engine is not available" in res.output
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/test_cli_engine.py --no-cov -q`
Expected: `No such option: --engine` failures (exit code 2 ≠ 0), plus `ImportError: cannot import name 'run'` from `cli_run` in `test_run_is_still_importable_from_cli`.

- [ ] **Step 3: Move `run` and add the options**

Create `cli_run.py`. Its body is cli.py:437-543 moved **verbatim**, with `@main.command()` changed to `@click.command()`, two new options before `@record_run_status`, two new parameters, and two new keyword arguments:

```python
# src/muc_one_span/cli_run.py
"""The ``run`` command: full pipeline execution (moved from ``cli.py``)."""

from __future__ import annotations

import click

from muc_one_span.cli_settings import current_configuration_path, current_settings
from muc_one_span.run_status import record_run_status
from muc_one_span.settings import DEFAULT_SETTINGS


@click.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(),
    help="Input FASTQ or BAM file.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(),
    default="results",
    help="Output directory.",
)
@click.option(
    "--reference",
    "-r",
    type=click.Path(),
    default=None,
    help="Reference FASTA (defaults to bundled ladder).",
)
@click.option(
    "--clair3-model",
    type=str,
    default=DEFAULT_SETTINGS.run.clair3_model,
    help="Path to Clair3 model.",
)
@click.option(
    "--threads", "-t", type=int, default=DEFAULT_SETTINGS.run.threads, help="Number of threads."
)
@click.option(
    "--min-coverage",
    type=int,
    default=DEFAULT_SETTINGS.run.min_coverage,
    help="Minimum read coverage.",
)
@click.option(
    "--min-qual",
    type=float,
    default=DEFAULT_SETTINGS.run.min_qual,
    help="Minimum VCF QUAL (default 5.0); see calling.haploid_min_qual for length-split calls.",
)
@click.option(
    "--report/--no-report",
    default=False,
    help="Generate HTML report (requires jinja2: pip install muc_one_span[report]).",
)
@click.option(
    "--report-igv",
    type=click.Choice(["embedded", "sidecar", "off"], case_sensitive=False),
    default="off",
    help="IGV alignment browser mode in HTML report (default: off).",
)
@click.option(
    "--platform",
    type=click.Choice(["hifi", "ont"], case_sensitive=False),
    default=DEFAULT_SETTINGS.run.platform,
    help="Sequencing platform (default: hifi).",
)
@click.option(
    "--mapping-timeout",
    type=float,
    default=DEFAULT_SETTINGS.run.mapping_timeout,
    help="Total mapping timeout in seconds (finite and positive; default 3600).",
)
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (auto-selected from --platform if not set).",
)
@click.option(
    "--engine",
    type=click.Choice(["ladder", "hybrid"]),
    default=DEFAULT_SETTINGS.run.engine,
    help="Allele reconstruction engine (hybrid is experimental; default: ladder).",
)
@click.option(
    "--assay",
    type=click.Choice(["amplicon", "genomic"]),
    default=DEFAULT_SETTINGS.run.assay,
    help="Library type used by the hybrid engine (default: amplicon).",
)
@record_run_status
def run(
    input_path: str,
    output_dir: str,
    reference: str | None,
    clair3_model: str,
    threads: int,
    min_coverage: int,
    min_qual: float,
    report: bool,
    report_igv: str,
    platform: str,
    minimap2_preset: str | None,
    mapping_timeout: float = DEFAULT_SETTINGS.run.mapping_timeout,
    engine: str = DEFAULT_SETTINGS.run.engine,
    assay: str = DEFAULT_SETTINGS.run.assay,
) -> None:
    """Run the full MucOneSpan pipeline."""
    from muc_one_span.pipeline import execute_pipeline

    execute_pipeline(
        input_path,
        output_dir,
        reference,
        clair3_model,
        threads,
        min_coverage,
        min_qual,
        report,
        platform,
        minimap2_preset,
        report_igv=report_igv,
        mapping_timeout=mapping_timeout,
        settings=current_settings(),
        configuration=current_configuration_path(),
        engine=engine,
        assay=assay,
    )
```

`cli_run.py` must not import `muc_one_span.cli`, because that would be a cycle. In `cli.py`, delete lines 437-545 (the old `run` command and its trailing blank lines), then:

```diff
--- a/src/muc_one_span/cli.py
+++ b/src/muc_one_span/cli.py
@@ -9,14 +9,13 @@
 
 import click
 
+from muc_one_span.cli_run import run
 from muc_one_span.cli_settings import (
     configure_context,
-    current_configuration_path,
     current_settings,
     validate_stage_options,
 )
 from muc_one_span.mapping import PLATFORM_PRESETS
-from muc_one_span.run_status import record_run_status
 from muc_one_span.settings import DEFAULT_SETTINGS
 from muc_one_span.version import __version__
 
@@ -53,6 +52,9 @@
     )
 
 
+main.add_command(run)
+
+
 @main.command()
 @click.option(
     "--output",
```

(`current_configuration_path` and `record_run_status` are no longer used in `cli.py`, and nothing imports them from there. Checked with a grep over `src`, `tests`, `scripts` and `docs`.)

- [ ] **Step 4: Pipeline kwargs and the temporary guard**

```diff
--- a/src/muc_one_span/pipeline.py
+++ b/src/muc_one_span/pipeline.py
@@ -28,6 +28,8 @@
     mapping_timeout: float | None = None,
     settings: RuntimeSettings | None = None,
     configuration: Path | None = None,
+    engine: str | None = None,
+    assay: str | None = None,
 ) -> None:
     """Run the full MucOneSpan pipeline."""
     from muc_one_span.alleles import detect_alleles, parse_idxstats
@@ -66,6 +68,7 @@
         mapping_timeout=mapping_timeout
         if mapping_timeout is not None
         else (settings or DEFAULT_SETTINGS).run.mapping_timeout,
+        **{k: v for k, v in (("engine", engine), ("assay", assay)) if v is not None},
     )
     if reference is None and (
         settings.repeat_dictionary is not None
@@ -88,6 +91,10 @@
     configuration_record = write_run_configuration(
         settings, configuration, Path(input_path), ref, out
     )
+    if settings.run.engine == "hybrid":  # Temporary: removed when the engine lands.
+        raise click.BadParameter(
+            "the hybrid engine is not available in this build", param_hint="--engine"
+        )
     igv_requested = settings.run.report_igv != "off"
     check_tools(
         ["minimap2", "samtools", "bcftools", "run_clair3.sh"]
```

- [ ] **Step 5: Run the CLI suites**

Run: `uv run --locked --all-extras pytest tests/unit/test_cli_engine.py tests/unit/test_cli_run.py tests/unit/test_cli.py tests/unit/test_wave1_cli.py tests/unit/test_pipeline_gates.py tests/unit/test_cli_settings.py --no-cov -q`
Expected: PASS. `cli.py` is ≈531 lines and `cli_run.py` 134. `make docs-check` still renders the CLI reference, because mkdocs-click walks `main` and `run` is registered.

- [ ] **Step 6: Full check and commit**

Run: `make ci-check && make docs-check`

```bash
git add src/muc_one_span/cli.py src/muc_one_span/cli_run.py src/muc_one_span/pipeline.py \
  tests/unit/test_cli_engine.py
git commit -m "refactor(cli): move run to cli_run; add --engine/--assay behind a temporary guard

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 11: Engine orchestration, shared tail, and pipeline branch

**Files:**
- Create: `src/muc_one_span/hybrid/engine.py`, `src/muc_one_span/pipeline_tail.py`, `tests/unit/hybrid/test_engine.py`
- Modify: `src/muc_one_span/hybrid/__init__.py`, `src/muc_one_span/pipeline.py`, `tests/unit/test_cli_engine.py` (delete the guard test)

**Interfaces:**
- Consumes: everything from Tasks 4–9; `InsufficientEvidenceError` (run_status.py:17); `tools.run_tool_iter` (tools.py:118); `write_run_configuration`.
- Produces:
  - `HybridResult(alleles, consensus_paths, block, members)`
  - `reconstruct_alleles(input_path: Path, output_dir: Path, rd: RepeatDictionary, settings: RuntimeSettings) -> HybridResult`
  - `read_input(path: Path) -> Iterator[ReadRecord]`, which streams FASTQ, FASTQ.gz, or BAM via `samtools fastq -F 0x900` (C9.5)
  - `extra_versions(backend: str) -> dict[str, str]`
  - `annotate_read_support(allele_key, sequence, result, *, rd, members, settings) -> dict`
  - `pipeline_tail.finish_run(*, out, input_path, rd, settings, alleles_result, consensus_paths, vcf_paths, tool_versions, configuration_record, report, bam_path, fasta_path, annotate=None, extra_summary=None) -> dict` (C9.3). It imports `classify_sequence`, `validate_mutations_against_vcf`, `parse_vcf_variants` and `generate_report` **inside** the function, so that the patch targets in test_pipeline_gates.py:66-88 and test_wave1_cli.py keep working.
- Output contract (pre-flight (c)), per allele:
  - `engine`, `length`, `canonical_repeats`, `reference_length`, `length_status = consistent_with_consensus_contig`, `length_basis = spanning_peak`, `consensus_basis = poa`;
  - `reads`, `spanning_reads`, `assigned_reads`;
  - `depth_status` ∈ {adequate, low, insufficient}, `depth_basis = spanning_reads`, `depth_threshold`;
  - `selection_status`, stamped on **both** alleles: resolved, unresolved_max_alleles, unresolved_single_site, unresolved_rejected_peak or unresolved_unassigned_spanning;
  - `secondary_mode_fraction = None`, `selection_detail`, `split_basis`, `phase_status`, `independent_haplotype_evidence`;
  - `residual_sites`, `heterozygous_sites`, `variant_filter = None`, `allele_genotype_status` (residual_heterogeneity or not_applicable_read_consensus);
  - `reconstruction_status` (complete_segmentation or read_consensus_low_depth), `sequence_source`, `contig_name`, `vcf_path = None`.
- Sample level:
  - `homozygous` is True only for a resolved single group with no residual sites. Otherwise it is False, `sequence_identity_status` is `unresolved`, and an `allele_2` alias is written (`candidate_duplicate_of`, `not_separately_resolved`).
  - More than 2 groups never silently drops an allele: `unresolved_max_alleles` is set.
  - Zero peaks raises `InsufficientEvidenceError`.
- Pipeline: the hybrid branch sits right after `write_run_configuration` and **before** `igv_requested`/`check_tools`, so Clair3 is never checked. `annotate_selection_qc` is **never** called on hybrid alleles. `summary["hybrid"]` is top level.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/hybrid/test_engine.py
"""Engine orchestration and the end-to-end hybrid pipeline on synthetic FASTQ."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.evaluation import load_observation
from muc_one_span.hybrid.engine import reconstruct_alleles
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.pipeline import execute_pipeline
from muc_one_span.report import compute_clinical_decision
from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import DEFAULT_SETTINGS
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
pytest.importorskip("pyabpoa", reason="pyabpoa (extra 'hybrid') is not installed")
A = synth.allele(["X"] * 25)
B = synth.allele(["X"] * 14 + [synth.dupc()] + ["X"] * 30)


def _fastq(path: Path, reads: list[ReadRecord]) -> Path:
    path.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in reads))
    return path


def _sample(tmp_path: Path, n_b: int = 60) -> Path:
    reads = synth.reads(A, 150, err=0.02, seed=1) + synth.reads(B, n_b, err=0.02, seed=2)
    return _fastq(tmp_path / "in.fastq", reads)


def test_heterozygous_dupc_sample_is_reconstructed(tmp_path: Path) -> None:
    result = reconstruct_alleles(_sample(tmp_path), tmp_path, synth.RD, DEFAULT_SETTINGS)
    seqs = {
        p.read_text().split("\n", 1)[1].replace("\n", "") for p in result.consensus_paths.values()
    }
    assert seqs == {A, B}
    assert "hybrid" not in result.alleles and "_members" not in result.alleles
    for key in ("allele_1", "allele_2"):
        info = result.alleles[key]
        assert info["engine"] == "hybrid" and info["depth_status"] == "adequate"
        assert info["selection_status"] == "resolved" and info["split_basis"] == "length"
        assert info["length"] == info["canonical_repeats"] + 9
    assert result.block["read_categories"]["spanning"] == 210
    assert result.block["poa_backend"] == "pyabpoa"
    saved = json.loads((tmp_path / "hybrid_reads.json").read_text())
    assert saved["read_categories"]["spanning"] == 210


def test_low_depth_allele_is_marked_not_dropped(tmp_path: Path) -> None:
    result = reconstruct_alleles(_sample(tmp_path, n_b=15), tmp_path, synth.RD, DEFAULT_SETTINGS)
    statuses = sorted(result.alleles[k]["depth_status"] for k in ("allele_1", "allele_2"))
    assert statuses == ["adequate", "low"]


def test_single_unconfirmed_site_is_not_homozygous(tmp_path: Path) -> None:
    a = synth.allele(["X"] * 10 + ["Q"] + ["X"] * 19)
    reads = synth.reads(a, 40, err=0.02, seed=3) + synth.reads(
        synth.allele(["X"] * 30), 40, err=0.02, seed=4
    )
    result = reconstruct_alleles(
        _fastq(tmp_path / "s.fastq", reads), tmp_path, synth.RD, DEFAULT_SETTINGS
    )
    assert result.alleles["homozygous"] is False
    assert result.alleles["allele_1"]["selection_status"] == "unresolved_single_site"
    assert result.alleles["allele_2"]["candidate_duplicate_of"] == "allele_1"


def test_no_spanning_reads_is_insufficient_evidence(tmp_path: Path) -> None:
    fq = _fastq(tmp_path / "e.fastq", [ReadRecord("j", "ACGT" * 400, "5" * 1600)])
    with pytest.raises(InsufficientEvidenceError):
        reconstruct_alleles(fq, tmp_path, synth.RD, DEFAULT_SETTINGS)


def test_hybrid_pipeline_end_to_end(tmp_path: Path) -> None:
    out = tmp_path / "out"
    with patch("muc_one_span.tools.check_tools") as check:
        execute_pipeline(
            str(_sample(tmp_path)),
            str(out),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )
    check.assert_called_once_with([])
    summary = json.loads((out / "summary.json").read_text())
    assert summary["hybrid"]["poa_backend"] == "pyabpoa"
    assert "hybrid" not in summary["alleles"]
    mutations = [m for c in summary["classifications"].values() for m in c["mutations"]]
    dupc = [m for m in mutations if m["mutation_name"] == "dupC"]
    assert len(dupc) == 1 and dupc[0]["read_support"]["status"] == "supported"
    assert dupc[0]["vcf_support"] is False
    assert dupc[0]["vcf_support_status"] == "not_applicable_read_consensus"
    assert compute_clinical_decision(summary)["state"] == "PATHOGENIC"
    assert load_observation(out).status != "invalid_artifacts"


def test_hybrid_rejects_igv_tracks(tmp_path: Path) -> None:
    import click

    with patch("muc_one_span.tools.check_tools"), pytest.raises(click.BadParameter):
        execute_pipeline(
            str(_sample(tmp_path)),
            str(tmp_path / "o"),
            None,
            "",
            1,
            10,
            5.0,
            False,
            "ont",
            None,
            report_igv="embedded",
            engine="hybrid",
            settings=DEFAULT_SETTINGS,
        )


def test_read_input_streams_fastq_gzip_and_bam(tmp_path: Path) -> None:
    import gzip

    from muc_one_span.hybrid.engine import read_input

    record = "@r1 extra\nacgt\n+\nIIII\n"
    (tmp_path / "a.fastq").write_text(record)
    with gzip.open(tmp_path / "a.fastq.gz", "wt") as handle:
        handle.write(record)
    lines = iter(["@r1", "ACGT", "+", "IIII"])
    with patch("muc_one_span.tools.run_tool_iter", return_value=lines) as tool:
        bam = list(read_input(tmp_path / "a.bam"))
    tool.assert_called_once_with(["samtools", "fastq", "-F", "0x900", str(tmp_path / "a.bam")])
    for got in (
        list(read_input(tmp_path / "a.fastq")),
        list(read_input(tmp_path / "a.fastq.gz")),
        bam,
    ):
        assert got == [ReadRecord("r1", "ACGT", "IIII")]


def test_truncated_fastq_is_an_error(tmp_path: Path) -> None:
    from muc_one_span.hybrid.engine import read_input

    (tmp_path / "t.fastq").write_text("@r1\nACGT\n")
    with pytest.raises(ValueError, match="truncated"):
        list(read_input(tmp_path / "t.fastq"))
```

Delete `test_hybrid_engine_is_rejected_until_available` from `tests/unit/test_cli_engine.py`.

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid/test_engine.py --no-cov -q`
Expected: `ModuleNotFoundError: No module named 'muc_one_span.hybrid.engine'`.

- [ ] **Step 3: Implement `engine.py` and the package export**

```python
# src/muc_one_span/hybrid/engine.py
"""Hybrid engine orchestration: spans -> lengths -> POA -> phase -> assign -> polish.

Returns the existing allele contract (plus added fields) and a separate ``block`` that
the pipeline stores as ``summary["hybrid"]``; nothing engine-specific is stored inside
``alleles`` except per-allele fields. Deviations from spec §3 in this version: no
ladder-assisted length prior or ladder-seeded consensus (S2/S3), depth is judged on
spanning reads only, and read-level support uses the spanning members only.
"""

from __future__ import annotations

import gzip
import json
import random
import statistics
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

from muc_one_span.config import RepeatDictionary
from muc_one_span.hybrid.assign import (
    FLANK_BP,
    OFF_TARGET,
    UNDECIDED,
    assign_reads,
    hybrid_references,
    trim_to_draft,
)
from muc_one_span.hybrid.evidence import event_read_support, residual_sites
from muc_one_span.hybrid.lengths import fit_length_model
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import draft_consensus, polish
from muc_one_span.hybrid.spans import Anchors, ReadRecord, SpanRead, categorize_reads
from muc_one_span.run_status import InsufficientEvidenceError
from muc_one_span.settings import HybridSettings, RuntimeSettings

UNIT = 60
FIXED_UNITS = 9
MAX_QC_READS = 200
PHASE_STATUS = {
    "length": "phased",
    "linked_sites": "phased",
    "none": "no_informative_heterozygosity",
    "unconfirmed_single_site": "unresolved_single_site",
}
_PACKAGES = {"pyabpoa": "pyabpoa", "pyspoa": "pyspoa"}


@dataclass
class HybridResult:
    """Engine output; ``block`` becomes ``summary["hybrid"]`` and is never put in alleles."""

    alleles: dict[str, Any]
    consensus_paths: dict[str, Path]
    block: dict[str, Any]
    members: dict[str, list[tuple[str, str]]]


def _fastq(lines: Iterable[str]) -> Iterator[ReadRecord]:
    it = iter(lines)
    for header in it:
        if not header.strip():
            continue
        try:
            seq, _plus, qual = next(it), next(it), next(it)
        except StopIteration:
            raise ValueError("truncated FASTQ record") from None
        yield ReadRecord(header[1:].split()[0], seq.strip().upper(), qual.strip())


def read_input(path: Path) -> Iterator[ReadRecord]:
    """Stream FASTQ(.gz) records, or the primary reads of a BAM via ``samtools fastq``."""
    if path.suffix == ".bam":
        from muc_one_span.tools import run_tool_iter

        yield from _fastq(run_tool_iter(["samtools", "fastq", "-F", "0x900", str(path)]))
    elif path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield from _fastq(handle)
    else:
        with path.open() as handle:
            yield from _fastq(handle)


def extra_versions(backend: str) -> dict[str, str]:
    """Versions of the optional packages that produced the consensus."""
    out = {}
    for pkg in ("edlib", _PACKAGES[backend]):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = "unknown"
    return out


def _depth_status(spanning: int, h: HybridSettings) -> str:
    if spanning >= h.depth_adequate_spanning:
        return "adequate"
    return "low" if spanning >= h.depth_low_spanning else "insufficient"


def _selection(model: Any, split_bases: list[str], n_groups: int, h: HybridSettings) -> str:
    """Sample-level selection status; any unresolved state blocks a negative call."""
    if n_groups > 2:
        return "unresolved_max_alleles"
    if "unconfirmed_single_site" in split_bases:
        return "unresolved_single_site"
    if model.gate_relevant_rejections:
        return "unresolved_rejected_peak"
    if model.unassigned_fraction > h.max_unassigned_spanning_fraction:
        return "unresolved_unassigned_spanning"
    return "resolved"


def _allele_info(
    name: str,
    seq: str,
    members: list[SpanRead],
    assigned: int,
    basis: str,
    residual: list[dict[str, Any]],
    selection: str,
    detail: str,
    h: HybridSettings,
) -> dict[str, Any]:
    units = round(len(seq) / UNIT)
    depth = _depth_status(len(members), h)
    return {
        "engine": "hybrid",
        "length": units,
        "canonical_repeats": units - FIXED_UNITS,
        "reference_length": units,
        "length_status": "consistent_with_consensus_contig",
        "length_basis": "spanning_peak",
        "consensus_basis": "poa",
        "reads": len(members),
        "spanning_reads": len(members),
        "assigned_reads": assigned,
        "depth_status": depth,
        "depth_basis": "spanning_reads",
        "depth_threshold": h.depth_adequate_spanning,
        "selection_status": selection,
        "secondary_mode_fraction": None,
        "selection_detail": detail,
        "split_basis": basis,
        "phase_status": PHASE_STATUS[basis],
        "independent_haplotype_evidence": basis in ("length", "linked_sites"),
        "residual_sites": residual,
        "heterozygous_sites": residual,
        "variant_filter": None,
        "allele_genotype_status": (
            "residual_heterogeneity" if residual else "not_applicable_read_consensus"
        ),
        "reconstruction_status": (
            "complete_segmentation" if depth == "adequate" else "read_consensus_low_depth"
        ),
        "sequence_source": f"hybrid:{name}",
        "contig_name": f"hybrid_{name}",
        "vcf_path": None,
    }


def reconstruct_alleles(
    input_path: Path, output_dir: Path, rd: RepeatDictionary, settings: RuntimeSettings
) -> HybridResult:
    """Reconstruct up to two allele sequences with explicit evidence and uncertainty."""
    h = settings.hybrid
    rng = random.Random(h.seed)
    backend = get_backend(h.poa_backend)
    cats = categorize_reads(read_input(input_path), Anchors.from_dictionary(rd), h)
    model = fit_length_model(cats.spanning, h)
    if not model.peaks:
        raise InsufficientEvidenceError("hybrid: no allele length peak passed the thresholds")
    groups: list[tuple[list[SpanRead], str, str]] = []
    split_bases = []
    for peak in model.peaks:
        draft = draft_consensus(peak.members, h.n_poa, rng, backend)
        split = split_by_linked_sites(draft, peak.members, h, rng)
        split_bases.append(split.basis)
        basis = split.basis
        if basis == "none" and len(model.peaks) == 2:
            basis = "length"
        for g in split.groups:
            groups.append((g, basis, draft if len(split.groups) == 1 else ""))
    selection = _selection(model, split_bases, len(groups), h)
    detail = f"{len(groups)} group(s); rejected peaks {model.gate_relevant_rejections}"
    groups = sorted(groups, key=lambda g: -len(g[0]))[:2]
    groups.sort(key=lambda g: statistics.median(m.length for m in g[0]))  # allele_1 shorter
    names = [f"allele_{i + 1}" for i in range(len(groups))]
    drafts = {
        n: d or draft_consensus(g, h.n_poa, rng, backend)
        for n, (g, _b, d) in zip(names, groups, strict=True)
    }
    refs = hybrid_references(drafts, rd)
    extra = assign_reads(
        cats.left_anchored + cats.right_anchored + cats.internal_or_offtarget, refs, h
    )
    alleles: dict[str, Any] = {}
    paths: dict[str, Path] = {}
    members: dict[str, list[tuple[str, str]]] = {}
    cap = h.n_poa * 3
    for name, (group, basis, _d) in zip(names, groups, strict=True):
        partial = [
            trim_to_draft(r, refs[name], FLANK_BP, len(drafts[name])) for r in extra[name][:cap]
        ]
        cons, info = polish(
            drafts[name],
            [m.seq for m in group][:cap],
            h.polish_rounds,
            h.hp_vote,
            partial=[p for p in partial if len(p) >= UNIT],
        )
        residual = residual_sites(cons, [m.seq for m in group][:MAX_QC_READS], h.qc_residual_af)
        alleles[name] = _allele_info(
            name, cons, group, len(extra[name]), basis, residual, selection, detail, h
        )
        alleles[name]["polish"] = info
        members[name] = [(m.seq, m.strand) for m in group]
        paths[name] = output_dir / f"consensus_{name}.fa"
        paths[name].write_text(f">{name}\n{cons}\n")
        (output_dir / f"consensus_{name}_context.json").write_text(
            json.dumps(
                {
                    "engine": "hybrid",
                    "reads": len(group),
                    "full_consensus_path": str(paths[name]),
                    "vcf_path": None,
                },
                indent=2,
            )
            + "\n"
        )
    residual_any = any(alleles[n]["residual_sites"] for n in names)
    homozygous = len(groups) == 1 and selection == "resolved" and not residual_any
    if len(groups) == 1:
        alleles["allele_2"] = {
            **alleles["allele_1"],
            "candidate_duplicate_of": "allele_1",
            "reconstruction_status": "not_separately_resolved",
            "independent_haplotype_evidence": homozygous,
        }
        alleles["sequence_identity_status"] = "resolved" if homozygous else "unresolved"
    lengths = [alleles[n]["length"] for n in names]
    alleles.update(
        {
            "homozygous": homozygous,
            "same_length": len(set(lengths)) == 1,
            "observed_length_candidates": [round(p.center_bp / UNIT) for p in model.peaks],
            "allele_multiplicity_status": "resolved" if len(groups) == 2 else "unresolved",
        }
    )
    block = {
        "engine": "hybrid",
        "assay": settings.run.assay,
        "poa_backend": backend.name,
        "read_categories": cats.counts(),
        "rejected_peaks": model.rejected,
        "short_product_fraction": round(model.short_product_fraction, 4),
        "unassigned_spanning_fraction": round(model.unassigned_fraction, 4),
        "undecided_reads": len(extra[UNDECIDED]),
        "off_target_reads": len(extra[OFF_TARGET]),
        "selection_status": selection,
        "split_bases": split_bases,
    }
    (output_dir / "hybrid_reads.json").write_text(json.dumps(block, indent=2) + "\n")
    (output_dir / "hybrid_references.fa").write_text(
        "".join(f">hybrid_{k}\n{v}\n" for k, v in refs.items())
    )
    return HybridResult(alleles, paths, block, members)


def annotate_read_support(
    allele_key: str,
    sequence: str,
    result: dict[str, Any],
    *,
    rd: RepeatDictionary,
    members: dict[str, list[tuple[str, str]]],
    settings: HybridSettings,
) -> dict[str, Any]:
    """Attach read-level support; VCF support does not apply to a read consensus."""
    support = event_read_support(sequence, result, members.get(allele_key, []), rd, settings)
    for idx, mutation in enumerate(result.get("mutations_detected", [])):
        read_support = support.get(idx, {"kind": "none", "status": "not_localized"})
        mutation["read_support"] = read_support
        mutation["vcf_support"] = False
        mutation["vcf_support_status"] = "not_applicable_read_consensus"
        mutation["support_status"] = f"read_level_{read_support['status']}"
    return result
```

```python
# src/muc_one_span/hybrid/__init__.py
"""Experimental read-centric (hybrid) MUC1 VNTR reconstruction engine."""

from muc_one_span.hybrid.engine import HybridResult, reconstruct_alleles

__all__ = ["HybridResult", "reconstruct_alleles"]
```

- [ ] **Step 4: Extract the shared tail**

Move pipeline.py:152-228 (`# Step 5` through `Pipeline complete`) verbatim into `finish_run`. Only the parameters and the two additions are new: the `annotate` hook after VCF validation, and `extra_summary` merged at the top level.

```python
# src/muc_one_span/pipeline_tail.py
"""Shared classification, summary and report tail for both reconstruction engines."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import click

from muc_one_span.config import RepeatDictionary
from muc_one_span.settings import RuntimeSettings
from muc_one_span.version import __version__

Annotate = Callable[[str, str, dict[str, Any]], dict[str, Any]]


def finish_run(
    *,
    out: Path,
    input_path: str,
    rd: RepeatDictionary,
    settings: RuntimeSettings,
    alleles_result: dict[str, Any],
    consensus_paths: dict[str, Path],
    vcf_paths: dict[str, Path],
    tool_versions: dict[str, str],
    configuration_record: dict[str, Any],
    report: bool,
    bam_path: Path | None,
    fasta_path: Path | None,
    annotate: Annotate | None = None,
    extra_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Classify each consensus, write repeats/summary files and the optional report."""
    # Imported here so that tests patching these module attributes keep working.
    from muc_one_span.classify import classify_sequence, validate_mutations_against_vcf
    from muc_one_span.vcf import parse_vcf_variants

    # Step 5: Classify repeats
    click.echo("Step 5/5: Classifying repeats...")
    all_results: dict[str, dict] = {}
    for allele_key, fa_path in consensus_paths.items():
        fa_lines = fa_path.read_text().strip().splitlines()
        sequence = "".join(line for line in fa_lines if not line.startswith(">"))
        result = classify_sequence(sequence, rd, settings=settings.classification)

        # VCF-backed validation if VCF available
        if allele_key in vcf_paths:
            vcf_variants = parse_vcf_variants(vcf_paths[allele_key])
            result = validate_mutations_against_vcf(
                result,
                vcf_variants=vcf_variants,
                sequence=sequence,
                repeat_dict=rd,
                consensus_context=alleles_result[allele_key].get("consensus_context"),
                settings=settings.confidence,
            )
        if annotate is not None:
            result = annotate(allele_key, sequence, result)

        all_results[allele_key] = result
        click.echo(f"  {allele_key}: {result['structure']}")
        if result.get("allele_confidence") is not None:
            click.echo(f"    confidence: {result['allele_confidence']:.2f}")

    (out / "alleles.json").write_text(json.dumps(alleles_result, indent=2) + "\n")

    # Write combined outputs
    (out / "repeats.json").write_text(json.dumps(all_results, indent=2) + "\n")
    structures = {k: v["structure"] for k, v in all_results.items()}
    (out / "repeats.txt").write_text("\n".join(f"{k}: {v}" for k, v in structures.items()) + "\n")

    # Summary
    summary = {
        "run_status": {"status": "analysis_completed"},
        "alleles": alleles_result,
        "classifications": {
            k: {
                "structure": v["structure"],
                "mutations": v["mutations_detected"],
                "reconstruction_status": v.get("reconstruction_status", "unverified"),
                "ambiguous_bases": v.get("ambiguous_bases", 0),
                "classification_coverage": v.get("classification_coverage"),
                "vcf_projection": v.get("vcf_projection"),
                "confidence_semantics": "heuristic_dictionary_fit_not_probability",
            }
            for k, v in all_results.items()
        },
        "tool_versions": tool_versions,
        "pipeline_version": __version__,
        "configuration": configuration_record,
        **(extra_summary or {}),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    effective_igv = settings.run.report_igv
    if report or effective_igv != "off":
        from muc_one_span.report import generate_report

        report_path = out / "report.html"
        generate_report(
            summary,
            report_path,
            sample_name=Path(input_path).stem,
            detailed_repeats=all_results,
            report_igv=effective_igv,
            bam_path=bam_path,
            vcf_paths=vcf_paths,
            fasta_path=fasta_path,
            execution_status={"status": "analysis_completed"},
        )
        click.echo(f"Report: {report_path}")

    summary["run_status"] = {"status": "completed"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    click.echo("Pipeline complete.")
    return summary
```

- [ ] **Step 5: Branch the pipeline (removes the Task 10 guard)**

```diff
--- a/src/muc_one_span/pipeline.py
+++ b/src/muc_one_span/pipeline.py
@@ -3,13 +3,17 @@
 from __future__ import annotations
 
 import json
+from functools import partial
 from pathlib import Path
+from typing import TYPE_CHECKING, Any
 
 import click
 
 from muc_one_span.cli_settings import effective_run_settings, write_run_configuration
 from muc_one_span.settings import DEFAULT_SETTINGS, RuntimeSettings
-from muc_one_span.version import __version__
+
+if TYPE_CHECKING:
+    from muc_one_span.config import RepeatDictionary
 
 
 def execute_pipeline(
@@ -34,14 +38,13 @@
     """Run the full MucOneSpan pipeline."""
     from muc_one_span.alleles import detect_alleles, parse_idxstats
     from muc_one_span.calling import call_variants_per_allele
-    from muc_one_span.classify import classify_sequence
     from muc_one_span.cli import PLATFORM_PRESETS, _bundled_reference
     from muc_one_span.config import load_repeat_dictionary
     from muc_one_span.consensus import build_consensus_per_allele
     from muc_one_span.mapping import get_idxstats, map_reads
+    from muc_one_span.pipeline_tail import finish_run
     from muc_one_span.selection_qc import annotate_selection_qc
     from muc_one_span.tools import check_tools, get_tool_versions
-    from muc_one_span.vcf import parse_vcf_variants
 
     out = Path(output_dir)
     out.mkdir(parents=True, exist_ok=True)
@@ -91,10 +94,9 @@
     configuration_record = write_run_configuration(
         settings, configuration, Path(input_path), ref, out
     )
-    if settings.run.engine == "hybrid":  # Temporary: removed when the engine lands.
-        raise click.BadParameter(
-            "the hybrid engine is not available in this build", param_hint="--engine"
-        )
+    if settings.run.engine == "hybrid":
+        _run_hybrid(out, input_path, rd, settings, configuration_record, report)
+        return
     igv_requested = settings.run.report_igv != "off"
     check_tools(
         ["minimap2", "samtools", "bcftools", "run_clair3.sh"]
@@ -156,80 +158,62 @@
         reference_layout=settings.reference_layout,
     )
 
-    # Step 5: Classify repeats
-    click.echo("Step 5/5: Classifying repeats...")
-    from muc_one_span.classify import validate_mutations_against_vcf
-
-    all_results: dict[str, dict] = {}
-    for allele_key, fa_path in consensus_paths.items():
-        fa_lines = fa_path.read_text().strip().splitlines()
-        sequence = "".join(line for line in fa_lines if not line.startswith(">"))
-        result = classify_sequence(sequence, rd, settings=settings.classification)
-
-        # VCF-backed validation if VCF available
-        if allele_key in vcf_paths:
-            vcf_variants = parse_vcf_variants(vcf_paths[allele_key])
-            result = validate_mutations_against_vcf(
-                result,
-                vcf_variants=vcf_variants,
-                sequence=sequence,
-                repeat_dict=rd,
-                consensus_context=alleles_result[allele_key].get("consensus_context"),
-                settings=settings.confidence,
-            )
-
-        all_results[allele_key] = result
-        click.echo(f"  {allele_key}: {result['structure']}")
-        if result.get("allele_confidence") is not None:
-            click.echo(f"    confidence: {result['allele_confidence']:.2f}")
+    finish_run(
+        out=out,
+        input_path=input_path,
+        rd=rd,
+        settings=settings,
+        alleles_result=alleles_result,
+        consensus_paths=consensus_paths,
+        vcf_paths=vcf_paths,
+        tool_versions=tool_versions,
+        configuration_record=configuration_record,
+        report=report,
+        bam_path=bam,
+        fasta_path=ref,
+    )
 
-    (out / "alleles.json").write_text(json.dumps(alleles_result, indent=2) + "\n")
 
-    # Write combined outputs
-    (out / "repeats.json").write_text(json.dumps(all_results, indent=2) + "\n")
-    structures = {k: v["structure"] for k, v in all_results.items()}
-    (out / "repeats.txt").write_text("\n".join(f"{k}: {v}" for k, v in structures.items()) + "\n")
-
-    # Summary
-    summary = {
-        "run_status": {"status": "analysis_completed"},
-        "alleles": alleles_result,
-        "classifications": {
-            k: {
-                "structure": v["structure"],
-                "mutations": v["mutations_detected"],
-                "reconstruction_status": v.get("reconstruction_status", "unverified"),
-                "ambiguous_bases": v.get("ambiguous_bases", 0),
-                "classification_coverage": v.get("classification_coverage"),
-                "vcf_projection": v.get("vcf_projection"),
-                "confidence_semantics": "heuristic_dictionary_fit_not_probability",
-            }
-            for k, v in all_results.items()
-        },
-        "tool_versions": tool_versions,
-        "pipeline_version": __version__,
-        "configuration": configuration_record,
-    }
-    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
-
-    effective_igv = settings.run.report_igv
-    if report or effective_igv != "off":
-        from muc_one_span.report import generate_report
-
-        report_path = out / "report.html"
-        generate_report(
-            summary,
-            report_path,
-            sample_name=Path(input_path).stem,
-            detailed_repeats=all_results,
-            report_igv=effective_igv,
-            bam_path=bam,
-            vcf_paths=vcf_paths,
-            fasta_path=ref,
-            execution_status={"status": "analysis_completed"},
-        )
-        click.echo(f"Report: {report_path}")
+def _run_hybrid(
+    out: Path,
+    input_path: str,
+    rd: RepeatDictionary,
+    settings: RuntimeSettings,
+    configuration_record: dict[str, Any],
+    report: bool,
+) -> None:
+    """Hybrid engine: no mapping, Clair3 or VCF; evidence comes from the reads."""
+    from muc_one_span.hybrid.engine import (
+        annotate_read_support,
+        extra_versions,
+        reconstruct_alleles,
+    )
+    from muc_one_span.pipeline_tail import finish_run
+    from muc_one_span.tools import check_tools
 
-    summary["run_status"] = {"status": "completed"}
-    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
-    click.echo("Pipeline complete.")
+    if settings.run.report_igv != "off":
+        raise click.BadParameter(
+            "IGV tracks are not available for the hybrid engine", param_hint="--report-igv"
+        )
+    check_tools(["samtools"] if str(input_path).endswith(".bam") else [])
+    click.echo("Hybrid engine: reconstructing alleles from reads...")
+    hybrid = reconstruct_alleles(Path(input_path), out, rd, settings)
+    (out / "alleles.json").write_text(json.dumps(hybrid.alleles, indent=2) + "\n")
+    finish_run(
+        out=out,
+        input_path=input_path,
+        rd=rd,
+        settings=settings,
+        alleles_result=hybrid.alleles,
+        consensus_paths=hybrid.consensus_paths,
+        vcf_paths={},
+        tool_versions=extra_versions(hybrid.block["poa_backend"]),
+        configuration_record=configuration_record,
+        report=report,
+        bam_path=None,
+        fasta_path=out / "hybrid_references.fa",
+        annotate=partial(
+            annotate_read_support, rd=rd, members=hybrid.members, settings=settings.hybrid
+        ),
+        extra_summary={"hybrid": hybrid.block},
+    )
```

- [ ] **Step 6: Run the engine tests and the ladder regression suites**

Run: `uv run --locked --all-extras pytest tests/unit/hybrid tests/unit/test_cli_engine.py tests/unit/test_pipeline_gates.py tests/unit/test_wave1_cli.py tests/unit/test_cli_run.py tests/unit/test_run_status.py --no-cov -q`
Expected: PASS.
- `test_hybrid_pipeline_end_to_end` reaches PATHOGENIC through read support, with `vcf_support False` and `not_applicable_read_consensus`, and `load_observation` is not `invalid_artifacts`.
- The ladder suites are unchanged (the patch targets are still honoured).

If sequences differ, print both and diff them before changing any setting.

- [ ] **Step 7: Full check and commit**

Run: `make ci-check`
Expected: PASS; total coverage ≈ 91%.

```bash
git add src/muc_one_span/hybrid src/muc_one_span/pipeline.py src/muc_one_span/pipeline_tail.py \
  tests/unit/hybrid/test_engine.py tests/unit/test_cli_engine.py
git commit -m "feat(hybrid): engine orchestration, P0 output contract and pipeline branch

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 12: Evaluation acceptance and harness forwarding

**Files:**
- Modify: `src/muc_one_span/evaluation/models.py`, `src/muc_one_span/evaluation/artifacts.py` (loader at :80-99), `src/muc_one_span/benchmarking.py`, `src/muc_one_span/clinical_runner.py` (argv at :142-157), `scripts/benchmark.py`, `scripts/clinical_benchmark.py`
- Test: `tests/unit/test_evaluation_artifacts.py`, `tests/unit/test_benchmark_tools.py`, `tests/unit/test_clinical_runner.py`

**Interfaces:**
- Produces:
  - `Event.read_support_status: str = "unknown"` (compare=False). `Event.supported` accepts `frameshift ∧ template_match ∧ read_support_status == "supported"`, and `legacy_supported` is unchanged (scoring.py:21,126).
  - The loader raises `"<name>: malformed read support"` unless `read_support` is a dict with a string `status`.
  - `_statuses_allow_completion` is **unchanged**: `unresolved_single_site` and `read_consensus_low_depth` fall to `ambiguous_reconstruction` by design.
  - `benchmarking.run_pipeline(sample, input_path, output_dir, platform, model, threads, *, runner=None, engine: str = "ladder")` is identical to `feat/benchsim` (D9).
  - `run_inventory(..., engine: str = "ladder")`.
  - `scripts/benchmark.py --engine`.
  - `clinical_runner.run_case` reads `settings.get("engine")`/`settings.get("assay")` from the hashed settings.
  - `scripts/clinical_benchmark.py run --engine/--assay` adds them to the hashed settings only when not default.
  - `benchmarking._stage_patches` times ladder stages only. Hybrid records carry wall time and empty stage timings, which is recorded as expected rather than a failure.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_evaluation_artifacts.py` (it uses the real helper `write_valid_artifacts(tmp_path)`, which returns None, and `observation.status`):

```python
def _edit_first_mutation(tmp_path, **changes):
    path = tmp_path / "summary.json"
    summary = json.loads(path.read_text())
    summary["classifications"]["p"]["mutations"][0].update(changes)
    path.write_text(json.dumps(summary))


def test_read_support_is_loaded_and_counts_as_support(tmp_path):
    write_valid_artifacts(tmp_path)
    _edit_first_mutation(
        tmp_path,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "supported", "n": 40},
    )
    event = load_observation(tmp_path).predictions[0].events[0]
    assert event.read_support_status == "supported"
    assert event.supported and not event.legacy_supported


def test_unsupported_read_evidence_is_not_support(tmp_path):
    write_valid_artifacts(tmp_path)
    _edit_first_mutation(
        tmp_path,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "discordant"},
    )
    assert not load_observation(tmp_path).predictions[0].events[0].supported


def test_malformed_read_support_is_invalid(tmp_path):
    write_valid_artifacts(tmp_path)
    _edit_first_mutation(tmp_path, read_support="supported")
    assert load_observation(tmp_path).status == "invalid_artifacts"


def test_unconfirmed_single_site_allele_is_ambiguous(tmp_path):
    write_valid_artifacts(tmp_path)
    path = tmp_path / "summary.json"
    summary = json.loads(path.read_text())
    summary["alleles"]["p"]["phase_status"] = "unresolved_single_site"
    path.write_text(json.dumps(summary))
    assert load_observation(tmp_path).status == "ambiguous_reconstruction"
```

Append to `tests/unit/test_benchmark_tools.py`, verbatim from `feat/benchsim` so the merge is trivial:

```python
def test_engine_is_forwarded_only_when_not_ladder(tmp_path: Path) -> None:
    from muc_one_span.benchmarking import run_pipeline

    class Result:
        exit_code = 0
        output = "ok"
        exception = None

    class Runner:
        def invoke(self, command, args):
            return Result()

    reads = tmp_path / "reads.fastq"
    reads.touch()
    default_record = run_pipeline(
        "sample", reads, tmp_path / "ladder", "ont", "model", 1, runner=Runner()
    )
    assert "--engine" not in default_record["cli_args"]
    assert default_record["engine"] == "ladder"

    hybrid_record = run_pipeline(
        "sample", reads, tmp_path / "hybrid", "ont", "model", 1, runner=Runner(), engine="hybrid"
    )
    assert hybrid_record["cli_args"][-2:] == ["--engine", "hybrid"]
    assert hybrid_record["engine"] == "hybrid"
```

Append to `tests/unit/test_clinical_runner.py`:

```python
def test_engine_and_assay_come_from_hashed_settings(tmp_path: Path, monkeypatch):
    import contextlib
    import json

    from muc_one_span.clinical_provenance import sha256_file

    reads = tmp_path / "reads.fastq"
    reads.write_text("@a\nACGT\n+\nIIII\n")
    seen: list[list[str]] = []

    def execute(commands, **kwargs):
        seen.append(json.loads(Path(commands[0][-1]).read_text())["argv"])
        raise RuntimeError("stop after recording argv")

    monkeypatch.setattr("muc_one_span.clinical_runner.run_tool_pipeline", execute)
    run = {"run_accession": "ERR1", "arm": "primary_amplicon"}
    prep = {"run_accession": "ERR1", "output_path": str(reads), "output_sha256": sha256_file(reads)}
    base = {"threads": 2, "timeout": 10, "model": "unused"}
    for root, extra in (("ladder", {}), ("hybrid", {"engine": "hybrid", "assay": "genomic"})):
        with contextlib.suppress(RuntimeError):
            run_case(run, prep, tmp_path / root, {**base, **extra})
    assert "--engine" not in seen[0]
    assert seen[1][-4:] == ["--engine", "hybrid", "--assay", "genomic"]
```

- [ ] **Step 2: Run to verify failure**

Run: `uv run --locked --all-extras pytest tests/unit/test_evaluation_artifacts.py tests/unit/test_benchmark_tools.py tests/unit/test_clinical_runner.py --no-cov -q`
Expected: 4 failures:
- `test_read_support_is_loaded_and_counts_as_support`
- `test_malformed_read_support_is_invalid`
- `test_engine_is_forwarded_only_when_not_ladder`
- `test_engine_and_assay_come_from_hashed_settings`

`test_unsupported_read_evidence_is_not_support` and `test_unconfirmed_single_site_allele_is_ambiguous` already pass as guards.

- [ ] **Step 3: Implement**

```diff
--- a/src/muc_one_span/evaluation/models.py
+++ b/src/muc_one_span/evaluation/models.py
@@ -19,11 +19,14 @@
     template_match: bool = field(default=False, compare=False)
     vcf_support: bool = field(default=False, compare=False)
     support_status: str = field(default="unknown", compare=False)
+    read_support_status: str = field(default="unknown", compare=False)
 
     @property
     def supported(self) -> bool:
-        """Return whether exact template and sequence-projection evidence agree."""
-        return self.legacy_supported and self.support_status == "exact_sequence_concordance"
+        """Exact template plus sequence-projection (VCF) or read-level support."""
+        if self.legacy_supported and self.support_status == "exact_sequence_concordance":
+            return True
+        return self.frameshift and self.template_match and self.read_support_status == "supported"
 
     @property
     def legacy_supported(self) -> bool:
```

```diff
--- a/src/muc_one_span/evaluation/artifacts.py
+++ b/src/muc_one_span/evaluation/artifacts.py
@@ -86,6 +86,9 @@
             raise ValueError(f"{name}: malformed event annotation")
         if not isinstance(m.get("vcf_support_status", "unknown"), str):
             raise ValueError(f"{name}: malformed variant support status")
+        read_support = m.get("read_support", {"status": "unknown"})
+        if not isinstance(read_support, dict) or not isinstance(read_support.get("status"), str):
+            raise ValueError(f"{name}: malformed read support")
         events.append(
             Event(
                 index,
@@ -95,6 +98,7 @@
                 m.get("template_match") is True,
                 m.get("vcf_support") is True,
                 m.get("vcf_support_status", "unknown"),
+                read_support["status"],
             )
         )
     from muc_one_span.settings import DEFAULT_SETTINGS
```

```diff
--- a/src/muc_one_span/benchmarking.py
+++ b/src/muc_one_span/benchmarking.py
@@ -108,8 +108,15 @@
     threads: int,
     *,
     runner: Runner | None = None,
+    engine: str = "ladder",
 ) -> dict[str, Any]:
-    """Run the real full CLI while timing its five scientific stages."""
+    """Run the real full CLI while timing its five scientific stages.
+
+    ``engine`` is appended to ``cli_args`` as ``--engine <engine>`` only when it
+    is not the default ``"ladder"``, and is always recorded in the returned
+    record and ``measurement.json``. Stage timings cover the ladder stages only;
+    a hybrid run records wall time.
+    """
     from muc_one_span.cli import main
 
     output_dir.mkdir(parents=True, exist_ok=True)
@@ -126,6 +133,8 @@
         "--platform",
         platform,
     ]
+    if engine != "ladder":
+        cli_args += ["--engine", engine]
     timings: dict[str, float] = {}
     started = time.perf_counter()
     with ExitStack() as stack:
@@ -152,6 +161,7 @@
         "cli_args": cli_args,
         "exit_code": int(result.exit_code),
         "timings": timings,
+        "engine": engine,
     }
     if error:
         record["error"] = error
@@ -210,6 +220,7 @@
     platform: str | None = None,
     model: str | None = None,
     threads: int | None = None,
+    engine: str = "ladder",
 ) -> list[dict[str, Any]]:
     """Run or explicitly fail every inventory entry without dropping denominators."""
     records: list[dict[str, Any]] = []
@@ -246,7 +257,13 @@
             if isinstance(run_threads, bool) or run_threads < 1:
                 raise ValueError(f"{sample}: threads must be a positive integer")
             record = run_pipeline(
-                sample, input_path, result_dir, run_platform, run_model, run_threads
+                sample,
+                input_path,
+                result_dir,
+                run_platform,
+                run_model,
+                run_threads,
+                engine=engine,
             )
         except (OSError, ValueError) as exc:
             record = {
```

```diff
--- a/src/muc_one_span/clinical_runner.py
+++ b/src/muc_one_span/clinical_runner.py
@@ -155,6 +155,9 @@
         "--report-igv",
         "off",
     ]
+    for option in ("engine", "assay"):  # hashed settings: --resume cannot mix engines
+        if settings.get(option) is not None:
+            argv += [f"--{option}", str(settings[option])]
     invocation = root / "invocation.json"
     write_json(invocation, {"argv": argv, "output": str(root)})
     command = [sys.executable, "-m", "muc_one_span.clinical_worker", str(invocation)]
```

```diff
--- a/scripts/benchmark.py
+++ b/scripts/benchmark.py
@@ -32,6 +32,9 @@
         "--clair3-model", default=None, help="Clair3 model (defaults to CLAIR3_MODEL)"
     )
     result.add_argument("--threads", type=int, default=None, help="Threads for every sample")
+    result.add_argument(
+        "--engine", choices=("ladder", "hybrid"), default="ladder", help="Reconstruction engine"
+    )
     return result
 
 
@@ -47,6 +50,7 @@
             platform=args.platform,
             model=args.clair3_model if args.clair3_model is not None else os.getenv("CLAIR3_MODEL"),
             threads=args.threads,
+            engine=args.engine,
         )
     except (OSError, ValueError) as exc:
         print(f"Error: {exc}", file=sys.stderr)
```

```diff
--- a/scripts/clinical_benchmark.py
+++ b/scripts/clinical_benchmark.py
@@ -67,6 +67,8 @@
     run.add_argument("--timeout", type=float, default=3600)
     run.add_argument("--run", action="append", default=[])
     run.add_argument("--resume", action="store_true")
+    run.add_argument("--engine", choices=("ladder", "hybrid"), default="ladder")
+    run.add_argument("--assay", choices=("amplicon", "genomic"), default=None)
     score = commands.add_parser("score")
     score.add_argument("--manifest", type=Path, required=True)
     score.add_argument("--truth", type=Path, required=True)
@@ -130,6 +132,10 @@
                         "environment": environment,
                         "environment_sha256": object_hash(environment),
                         "report_igv": "off",
+                        # Only non-default engine choices enter the hashed settings, so
+                        # existing ladder attempts stay resumable.
+                        **({"engine": args.engine} if args.engine != "ladder" else {}),
+                        **({"assay": args.assay} if args.assay else {}),
                         "repeat_policy": "one observed execution per library",
                         "harness_sha256": {
                             str(p.relative_to(checkout)): sha256_file(p)
```

- [ ] **Step 4: Run the suites**

Run: `uv run --locked --all-extras pytest tests/unit/test_evaluation_artifacts.py tests/unit/test_evaluation_scoring.py tests/unit/test_evaluation_cli.py tests/unit/test_benchmark_tools.py tests/unit/test_clinical_runner.py tests/unit/test_clinical_benchmark_cli.py --no-cov -q`
Expected: PASS. The pinned argv tests (test_benchmark_tools.py:108-141) are unchanged.

- [ ] **Step 5: Full check and commit**

Run: `make ci-check`

```bash
git add src/muc_one_span/evaluation src/muc_one_span/benchmarking.py src/muc_one_span/clinical_runner.py \
  scripts/benchmark.py scripts/clinical_benchmark.py tests/unit/test_evaluation_artifacts.py \
  tests/unit/test_benchmark_tools.py tests/unit/test_clinical_runner.py
git commit -m "feat(eval): score read-level support; forward --engine/--assay through the harnesses

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 13: Integration test and regression runs (outputs outside Git)

**Files:**
- Create: `tests/integration/test_hybrid_engine.py`, `benchmarks/clinical/prjeb92208/hybrid-engine.json` (scores, hashes and versions only)
- Create (outside Git): `../MucOneSpan-review-20260923/engine-regression/{simpanel,heldout,prjeb92208,inhouse}/…`

- [ ] **Step 1: Integration test (C11.1)**

```python
# tests/integration/test_hybrid_engine.py
"""Hybrid engine on the repository's generated HiFi sample (samtools + generated data)."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from muc_one_span.cli import main

pytestmark = pytest.mark.integration
DATA = Path(__file__).resolve().parents[1] / "data" / "generated" / "sample_close_51_58"


def test_hybrid_run_on_generated_sample(tmp_path: Path) -> None:
    pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
    pytest.importorskip("pyabpoa", reason="pyabpoa (extra 'hybrid') is not installed")
    if shutil.which("samtools") is None:
        pytest.skip("samtools is not on PATH")
    bams = sorted(DATA.glob("*_amplicon_aligned.bam")) if DATA.is_dir() else []
    if not bams:
        pytest.skip(f"generated test data missing: {DATA} (run make generate-testdata)")
    args = ["run", "-i", str(bams[0]), "-o", str(tmp_path), "--engine", "hybrid", "--no-report"]
    res = CliRunner().invoke(main, args)
    assert res.exit_code == 0, res.output
    summary = json.loads((tmp_path / "summary.json").read_text())
    lengths = sorted(summary["alleles"][k]["length"] for k in ("allele_1", "allele_2"))
    assert lengths == [51, 58]
    assert summary["hybrid"]["poa_backend"] == "pyabpoa"
```

Run: `make test-int`
Expected: PASS where `samtools` and the generated data exist (validated on the main checkout: `[51, 58]`). The hybrid worktree has no `tests/data/generated` (untracked), so run `make generate-testdata` first or **report the skip and its reason**. A skip is not validation.

- [ ] **Step 2: Sub-population investigation (Review Focus #1; blocks any tuning)**

For the generated negatives `normal_50_55`, `normal_50_60`, `normal_60_80` and `normal_100_120`:
1. Run `--engine hybrid`.
2. Extract the linked or unconfirmed sites. The site table is `PhaseResult.sites`; get it with a short script over `split_by_linked_sites`, as in the plan validation.
3. Compare the minor-group reads against the two truth haplotypes (`*.simulated.fa`, `*.vntr_structure.txt`).

Classify each minor group as an inter-haplotype chimera, a simulator artefact, or unexplained. Write the result, with no sequences, to `../MucOneSpan-review-20260923/engine-regression/subpopulations.md`. Do not change `het_*` or `link_phi_min` in this task.

- [ ] **Step 3: Frozen simulated panels (C11.2)**

```bash
R=../MucOneSpan-review-20260923
for panel in simpanel heldout; do
  out=$R/engine-regression/$panel
  mkdir -p "$out/runs"
  python3 - "$R/$panel/manifest.json" > "$out/cases.tsv" <<'EOF'
import json, sys
m = json.load(open(sys.argv[1]))
for c in m["cases"].values():
    print(c["case_id"], c["reads"], c["muconespan_platform"], sep="\t")
EOF
  while IFS=$'\t' read -r cid reads plat; do
    uv run --locked --all-extras muconespan run -i "$reads" -o "$out/runs/$cid" \
      --engine hybrid --platform "$plat" --no-report -t 1
  done < "$out/cases.tsv"
  python3 $R/simpanel/adapt_muconespan.py "$out/runs" "$out/preds"
  uv run --locked --all-extras python $R/simpanel/score.py --manifest $R/$panel/manifest.json \
    --pred "$out/preds" --out "$out/scores" --name hybrid
done
```

Expected (spec §7.2): simpanel ≥ 80/80 alleles sequence-exact, held-out ≥ 76/80. The 3 known prototype failures must be correct or INCONCLUSIVE, never confident and wrong. Report the NO_PATHOGENIC yield on negatives next to the ladder's, and do not tune.

- [ ] **Step 4: PRJEB92208 (C11.3)**

Run `scripts/clinical_benchmark.py run ... --engine hybrid` into a **new** output root. Add `--assay genomic` for the WGS libraries, as a separate invocation. Expected (spec §7.3):
- MP1–MP4 PATHOGENIC with `read_support.status == "supported"`;
- HG002 2/2 exact against Q100;
- HG001–HG004 not PATHOGENIC.

Write `benchmarks/clinical/prjeb92208/hybrid-engine.json` with scores, hashes, `poa_backend` and the `hybrid` extra's package versions, and **no sequences**.

- [ ] **Step 5: In-house genomic (local only)**

Run `--engine hybrid --assay genomic` on `../MucOneSpan-review-20260923/inhouse/reads/*.fastq`. Expected:
- no phantom fragment alleles;
- the in-house dupC sample (82-unit allele, 7 spanning reads) has `depth_status` `insufficient` and the decision is INCONCLUSIVE, neither PATHOGENIC nor NEGATIVE.

Results stay out of Git.

- [ ] **Step 6: Commit the test and the public record**

Run: `make ci-check && make test-int`

```bash
git add tests/integration/test_hybrid_engine.py benchmarks/clinical/prjeb92208/hybrid-engine.json
git commit -m "test(hybrid): integration test and PRJEB92208 regression record

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 14: Documentation and changelog

**Files:**
- Modify: `docs/reference/cli.md` (mkdocs-click; `run` stays registered), `docs/guides/configuration.md`, `docs/reference/limitations.md`, `docs/getting-started/concepts.md`, `examples/runtime-settings.json`, `CHANGELOG.md` (`## [Unreleased]`), `README.md` (one line: `--engine hybrid`, experimental)

- [ ] **Step 1: Write the docs**

Cover the following:
- The stages (spec §3) and the deviations listed under "Decisions".
- The evidence fields (pre-flight (c)) and `summary["hybrid"]`.
- Every `hybrid.*` setting with its default, including `poa_backend`, `assign_max_error_rate`, `max_unassigned_spanning_fraction`, `smear_min_prominence`, `rejected_peak_noise_reads` and `hp_min_strand_reads`, each marked provisional and tuned on the dev/val splits only.
- That `hybrid` is experimental and not the default, and that `--report-igv` is unavailable with it.
- The genomic-assay depth caveat.
- The platform/wheel matrix: pyabpoa is sdist-only and needs a C compiler and zlib; edlib builds from source on 3.14; pyspoa has no macOS wheel.
- The licence note (MIT; medaka and dorado are not used).

CHANGELOG `Unreleased` (do **not** re-announce the explicit-support decision released in 0.16.0):
- `--engine hybrid` (experimental) and `--assay`;
- the new `hybrid` extra and settings section;
- the fix that makes `insufficient` depth, `residual_heterogeneity` and non-supported read-level evidence gate the clinical decision;
- `run` moved to `cli_run.py` (import path preserved).

- [ ] **Step 2: Verify and commit**

Run: `make docs-check && make build-check && make ci-check`
Expected: PASS; the wheel contains no test data.

```bash
git add docs examples CHANGELOG.md README.md
git commit -m "docs: document the experimental hybrid engine, evidence fields and settings

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## Self-review

- **Pre-flight coverage:**
  - C2.1 → Task 2 (plus the D7 noise setting and dropping `n_poa_min`);
  - C2.2 → Task 2 (D1–D4);
  - C2.2b/C2.3 → Task 10;
  - C3.1–C3.3 → Task 4 (C3.3's `collect_ignore_glob` is not needed under D1: module-level `importorskip` handles a missing extra, and CI always installs it);
  - C4.1/C4.2 → Task 5;
  - C5.1–C5.3 → Task 6;
  - C6.1 → Task 8 (unit Q, D5);
  - C7.1/C7.2 → Task 7;
  - C8.1–C8.3 → Task 9;
  - C9.1–C9.6 → Task 11 (the IGV-with-hybrid case is rejected, D8);
  - C10.1–C10.3 → Task 3;
  - C10.4/C10.5 → Task 12 (D9 signature);
  - C11.1–C11.3 → Task 13;
  - Task 12 row → Task 14.
- **Test-defect fixes:**
  - `synth.dupc` uses the dictionary template (Task 4);
  - decided-read expectations only count decidable fragments (Task 7);
  - `decision["state"]` everywhere (Task 3);
  - the CLI test asserts the `execute_pipeline` kwargs, not `current_settings()` (Task 10);
  - `importorskip` is at module level with lazy edlib (Tasks 4–9);
  - the manifest loop uses `.values()` (Task 13).
- **Untested in validation:** `uv.lock` resolution, pyabpoa builds on 3.10–3.13, `make docker-test`, `make docs-check`, `make build-check`, the frozen panels, PRJEB92208 and the in-house runs.
