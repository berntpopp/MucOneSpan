# Pre-flight: hybrid engine plan (Tasks 2–12) against v0.16.0 (P0)

Date: 2026-09-24. Base: `feat/hybrid-engine` = `fix/p0-clinical-safety` HEAD `2675e9b`
(v0.16.0). Checked: `.planning/2026-09-23-hybrid-engine-plan.md` (the plan) and
`.planning/2026-09-23-hybrid-engine-spec.md` (the spec). The plan file is unchanged.
Apply this addendum task by task. Where the addendum and the plan disagree, the addendum wins.

Method: every file, line, symbol, and fixture the plan names was read at `2675e9b`.
The plan's synthetic test data was checked against the bundled `repeats.json`.
Wheel availability was queried from PyPI JSON on 2026-09-24. `edlib`, `pyspoa`, and
`pyabpoa` were installed into a scratch CPython 3.14.4 venv. No code was changed and
nothing was committed.

---

## (a) Findings by task

| Task | Findings | Required correction |
| --- | --- | --- |
| 1 | Done in P0. `clinical_gates.mutation_supported` keeps the `read_support.status == "supported"` branch (clinical_gates.py:14-27, verified). The decision key is `state`, not `status`. The P0 tests `test_read_support_is_accepted_as_explicit_support` etc. are at test_clinical_decision.py:178-200. | None. Later tasks must use `decision["state"]` and reuse the existing `_summary`, `BASE`, `RESOLVED_GATES`, and `_gated_summary` helpers (test_clinical_decision.py:140-221). Do not redefine them. |
| 2 | (1) `settings.py` is 432 lines; with `HybridSettings` it becomes about 500 (OK). (2) The helpers `_integer/_number/_boolean/_choice` exist with the signatures the plan assumes. The settings convention names errors `run.<field>`. (3) The `run` command is at cli.py:437-543 (anchor valid). It calls `execute_pipeline` positionally, and **`effective_run_settings` is called in pipeline.py:55, not in the CLI**. The plan's "pass engine into the existing effective_run_settings call" points at the wrong file. (4) The plan's CLI test reads `ex.call_args.kwargs["settings"].run.engine`. That value is `current_settings()` (config file or defaults), not the CLI override, so the test fails. (5) `configure_context` maps any `run.*` field name to a same-named Click parameter, so `--config` defaults for `engine`/`assay` work once the options exist (verified). (6) A commit that accepts `--engine hybrid` before Task 9 would silently run the ladder engine and record `engine=hybrid` in `run_configuration.json`. (7) Dependency facts: see (d). `numpy` is **not** installed and not a core dependency. The quality env (`UV_QUALITY`) lacks any hybrid extra, so mypy fails on `import numpy`. (8) Coverage is measured on 3.14 only (test.yml:104-109). If hybrid tests skip there, coverage drops from about 89.6% to about 79–80% (estimated). | Split into **2a** (settings + extra) and **2b** (move `run` to `cli_run.py` + options + pipeline kwargs). Use `run.engine`/`run.assay` names in `_choice`. Add kw-only `engine`/`assay` to `execute_pipeline` and forward them into `effective_run_settings` only when not None. Replace the CLI test (C2.3). Add a temporary guard that rejects `engine == "hybrid"` until Task 9. Drop `numpy` (C4.1). Decide the edlib marker per (d), together with the coverage interpreter. Add the settings the gates need (C2.1). |
| 3 | (1) `synth.dupc` appends `C` after the unit. That is **not** the dupC template. The dictionary inserts at 1-based 60, which is 0-based 59, before the final `A`, giving 7C to 8C (config.py:63-67; `rd.mutated_sequences` maps sequence to `(parent, "dupC")`). The plan's allele therefore carries an untemplated `...CCCCCCCAC` unit. Tasks 5/9 would still "pass", because they only compare sequences, but classify would not report dupC. (2) `align.py` imports `edlib` at module import, and `synth.py` imports `align.rc`. Any environment without edlib fails at collection, and `pytest.importorskip` placed after the imports (Task 5 test) cannot help. (3) Flanks are 10 000 bp and motifs 1/9 are 60 bp. The `from_dictionary` fields exist (`rd.repeats`, `rd.flanking_left/right`). (4) `tests/__init__.py` and `tests/unit/__init__.py` exist, so `from tests.unit.hybrid import synth` resolves. (5) `categorize_reads` sends off-target reads into `internal_or_offtarget`, which Task 9 later assigns (see Task 7). | Fix `dupc` (C3.1). Load edlib lazily in `align.py` with the extra hint (C3.2). Put `rc`/`cigar_ops` before any edlib use so that `synth` imports without edlib. Add `tests/unit/hybrid/conftest.py` with a `collect_ignore_glob` only if the edlib marker is kept (C3.3). |
| 4 | (1) Uses numpy, which is not available (see Task 2). (2) The smear rule is dead logic. `smear = c < min(kept) - 1.5*UNIT and not far` is only true for candidates 1.5–2 units below the major peak. A broad smear 3–40 units below is judged by `far_peak_min_frac` (3%) and can be **accepted as allele 2**. That is Review Focus #2, and the prototype `length_peaks` has no smear logic either (proto.py:121-153). `test_smear_is_short_product_not_allele` passes only if no smear KDE maximum reaches 8 reads. (3) `"max_alleles"` rejections (a third real peak) are not treated as gate-relevant by Task 10. | Pure-Python density (C4.1). Replace the smear test with a prominence rule (C4.2); its constant is new and must be tuned on the dev split. Record `max_alleles` as gate-relevant (Task 9/10). |
| 5 | (1) `test_polish.py` calls `pytest.importorskip("edlib")` *after* importing modules that import edlib, which is ineffective. (2) `get_backend()` silently falls back pyabpoa→pyspoa, so results depend on what is installed (non-reproducible). The prototype evidence (80/80, 76/80) is **pyabpoa**-based (proto.py:229-232). (3) **Critical:** `pileup_polish` projects every read with global NW (`global_columns`). Task 9 feeds non-spanning fragments into the polish pool. Their uncovered columns vote `-`, and in genomic data (fragments ≫ spanning) this truncates the consensus ends or drops bases. (4) pyabpoa/pyspoa call signatures used by the plan were verified on 3.14 (`msa_aligner(aln_mode="g").msa(..., out_cons=True, out_msa=False).cons_seq[0]`; `spoa.poa(seqs, algorithm=1)` returns `(cons, msa)`). | Make the backend an explicit setting with no silent fallback, and record it (C2.1, C5.1). Add an infix projection for partial reads, and make pileup ignore uncovered columns (C5.2). Until that lands, polish with spanning members only. |
| 6 | (1) Prototype reference: `hetsplit.py` is 251 lines (functions at 87-251), not 87-260. The EM `phase_reads` (164-203) is not ported; majority vote is used instead. That is acceptable, but record the deviation. (2) **Test mismatch:** the described single-site test ("only the `A` difference") does not produce `unconfirmed_single_site`. Unit `A` differs from `X` at **4** bases (B: 2, G: 2, C: 1; measured), so one A-vs-X unit gives 4 perfectly linked sites and results in `linked_sites`. | Single-site test uses unit `C` (1 base) (C6.1). Add an explicit test that a multi-base single-unit difference (`A`) yields `linked_sites`, and document this as intended. Alternatively decide to count loci rather than sites; this is a scientific decision for the user. |
| 7 | (1) **Test will fail:** A/B differ only by one extra `X` inside a 20/21-X tract (1 200–1 260 bp). A 1 500 bp fragment is decided only if it spans the whole tract: start ∈ [~780, ~1 080] of [0, 1 800], about 17%, or about 10/60 decided. `decided >= 20` fails. (2) Off-target reads ≥ `min_fragment_bp` are assigned whenever the distance gap ≥ 3, even at 40% error. They then enter polishing. (3) The Task 7 prose says a homozygous single reference assigns all reads, but Task 9 skips assignment when `len(refs) == 1`. | Rewrite the test (C7.1). Add a maximum-error guard `hybrid.assign_max_error_rate` (default 0.15; a new setting to tune) and an `off_target` bucket (C7.2). State the homozygous behaviour once: assign to the single reference with the error guard. |
| 8 | (1) **Wrong event typing:** `_runs(cons[start:end], 6)` detects a ≥6 run in *every* X-family unit, since each has a 7C tract. Every event in such a unit, including substitutions, would be scored with the homopolymer LLR. (2) The `supported` rule ignores strand and never emits the spec's `discordant`. (3) The `exact_window` "ref" rule (`len(piece) % 60 == 0`) is not meaningful. (4) The test expectations were recomputed by hand and the plan's LLR tests pass (mut ≈ +46.6, wt < 0). `residual_sites` test position 50 is outside runs ≥3 (OK). (5) Reads are spanning members only (spec: "assigned reads"). This is acceptable for now; record it. | Type homopolymer events from the dictionary template (C8.1). Add strand discordance and the `discordant` status (C8.2). Replace exact-window with parent-vs-template competition (C8.3). Producer contract in (c). |
| 9 | (1) Anchors moved. The shared tail is **pipeline.py:152-228** (`# Step 5` through `Pipeline complete`), not 142-217. `write_run_configuration` is at 88. `check_tools` for Clair3 is at 91-100 and must **not** run for hybrid, so the branch goes at line 91, before `igv_requested`/`check_tools`. (2) P0 calls `annotate_selection_qc` (pipeline.py:120) on the ladder path. On hybrid alleles it would **overwrite** `depth_status` with `not_assessed` (there are no `fit_metrics`/`primary_alignment_records`), which disables the depth gate. (3) Test patch targets: `test_pipeline_gates.py:66-88` and `test_wave1_cli.py` patch `muc_one_span.classify.classify_sequence`, `muc_one_span.vcf.parse_vcf_variants`, `muc_one_span.report.generate_report`, and `muc_one_span.tools.check_tools`. `pipeline_tail` must import these **inside the function** or the patches stop applying. (4) **Invalid artifacts:** the plan puts `alleles["hybrid"] = {...}` inside `alleles`. `evaluation/artifacts.py:229` treats every dict-valued key of `alleles` as an allele, so `"hybrid"` becomes an unresolved alias and `load_observation` returns `invalid_artifacts`. The same applies to the `_members` side channel. (5) With 2 length peaks where one peak also splits by linked sites, there are 3 groups; the plan keeps the top 2 by support and **silently drops one allele**. `basis = "length"` also masks an `unconfirmed_single_site` result inside a peak. (6) With zero accepted peaks, `alleles` has no `allele_1` and the tail/report crash. (7) `homozygous = len(groups) == 1` makes an `unconfirmed_single_site` sample "homozygous". `report.py` only checks reconstruction reasons when not homozygous, and P0 set the principle that one length is not identity (calling.py:407-410). (8) `read_input` uses `run_tool` (captures all stdout), against package guidance "avoid capturing unbounded sequencing data"; `tools.run_tool_iter` exists. (9) The P0 gate fields are missing from `_allele_info` (see (c)). (10) `benchmarking._stage_patches` times ladder stages only. For hybrid runs the timings are empty; record that and do not treat it as a failure. | Apply C9.1–C9.6: `HybridResult` return type; branch placement; skip `annotate_selection_qc` for hybrid; top-level `summary["hybrid"]`; group/basis logic; zero-peak `InsufficientEvidenceError`; streamed input; the full output contract; stronger engine test (dupC detected, read support supported, decision PATHOGENIC). |
| 10 | (1) The Task 1 helper exists; the NO_PATHOGENIC hybrid gates are **not** implemented yet. (2) The plan's tests use `decision["status"]`; use `state`. `test_hybrid_low_depth_blocks_negative` already passes on P0 (`depth_status == "low"` is gated), so it is not a red test. (3) **P0 gap:** `depth_status == "insufficient"` (a spec value) is recognised nowhere. `report.py:98` (`depth_assessed`), `report.py:116` (carrier blocker), and `clinical_gates.allele_gate_reasons` only check `"low"`. An insufficient-depth hybrid allele would pass the depth gate, and with no assessed depth the legacy total-read fallback applies. (4) Preferred design: express hybrid sample-level blockers through the existing P0 per-allele fields (`selection_status = "unresolved_*"`, `allele_genotype_status`) rather than a new `summary["hybrid"]` branch in report.py. This keeps one gate path for both engines. (5) The depth message is hard-wired to `primary_alignment_records`. Existing tests assert `"Allele 2: 12 primary alignments"` (test_clinical_decision.py:307,328), so keep that text for the ladder basis. (6) `evaluation/artifacts.py:130-153` anchor still valid; `Event` at models.py:12-31; the loader at artifacts.py:80-99. (7) `benchmarking.run_pipeline` is at :102-128 (anchor valid). In `clinical_runner.run_case` the argv is at :142-157. The engine must go into the hashed `settings` dict (provenance, :125), not a separate kwarg, so that `--resume` cannot mix engines. | C10.1–C10.5. |
| 11 | (1) The integration test uses a cwd-relative `Path("tests/data/generated/...")`. Existing tests use `Path(__file__).resolve().parents[1] / "data" / "generated"`. The data exists in the main checkout (`sample_close_51_58_reads_amplicon_aligned.bam`) but **not** in the hybrid worktree (untracked), so the test skips there. (2) Heldout `manifest.json` is `{"cases": {case_id: {...}}}`, a **dict**. The plan's loop `for c in cases` iterates keys and fails. The fields are `case_id`, `reads`, `platform`, `muconespan_platform`. (3) `../MucOneSpan-review-20260923/engine-regression/` does not exist yet (created by the task). `simpanel/score.py`, `simpanel/adapt_muconespan.py`, and both manifests exist. (4) PRJEB92208 via `clinical_runner` hard-codes `--platform ont` (clinical_runner.py:148-149). Hybrid needs `settings["engine"]` and, for WGS libraries, `settings["assay"]` in the hashed settings. | C11.1–C11.3. |
| 12 | All named docs exist (`docs/reference/cli.md` uses mkdocs-click on `muc_one_span.cli:main`, so moving `run` is transparent once it is registered). `CHANGELOG.md` has `## [Unreleased]`. The CHANGELOG entry "explicit-support clinical decision" is **already released in 0.16.0**; do not re-announce it. | Replace that bullet with "hybrid read-level support accepted by existing gates; `insufficient` depth now gated". Document `hybrid.*` settings including the new ones (C2.1), the backend choice, and the platform/wheel matrix. |

---

## (b) Corrections addendum

### Task 2a — settings and extra

**C2.1 `settings.py`.** In `RunSettings`, add after `report_igv`:

```python
    engine: str = "ladder"
    assay: str = "amplicon"
```

and in `RunSettings.__post_init__`:

```python
        _choice("run.engine", self.engine, ("ladder", "hybrid"))
        _choice("run.assay", self.assay, ("amplicon", "genomic"))
```

`HybridSettings`: use the plan's block. Also add these fields. Each one is needed by a later correction; mark all as provisional and tune on dev/val only:

```python
    poa_backend: str = "pyabpoa"            # C5.1: explicit, recorded; no silent fallback
    max_unassigned_spanning_fraction: float = 0.2   # spec §4 NO_PATHOGENIC rule (was a literal)
    assign_max_error_rate: float = 0.15     # C7.2: off-target guard (new)
    smear_min_prominence: float = 3.0       # C4.2: smear rule (new)
    hp_min_strand_reads: int = 5            # C8.2: strand discordance (new)
```

with validation:

```python
        _choice("hybrid.poa_backend", self.poa_backend, ("pyabpoa", "pyspoa"))
        _number("max_unassigned_spanning_fraction", self.max_unassigned_spanning_fraction, 0, 1)
        _number("assign_max_error_rate", self.assign_max_error_rate, 0, 1)
        _number("smear_min_prominence", self.smear_min_prominence, 1)
        _integer("hp_min_strand_reads", self.hp_min_strand_reads)
```

Register it in `RuntimeSettings` and in `_SECTIONS` (settings.py:332-341). `load_settings` iterates `_SECTIONS`, so no other change is needed. The plan's settings tests work unchanged. Add `assert DEFAULT_SETTINGS.hybrid.poa_backend == "pyabpoa"`.

**C2.2 Extra/tooling.** See (d) for the decision. Recommended `pyproject.toml`:

```toml
hybrid = [
    "edlib>=1.3.9",          # must use the SAME requirement string as feat/benchsim's `bench`
    "pyabpoa>=1.5.7",        # sdist only
    "pyspoa>=0.3.2",
]

[[tool.mypy.overrides]]
module = ["edlib", "pyabpoa", "spoa"]
ignore_missing_imports = true
```

Drop `numpy` (C4.1). In the `Makefile`, add `--extra hybrid` to **both** `UV_TEST` and `UV_QUALITY` (lines 3-4, matching how benchsim adds `bench`). In `docker/Dockerfile` (lines 25 and 28), add `--extra hybrid`. The builder stage is micromamba `tools-runtime` and has no C compiler, so either add a toolchain to the builder, or install `pyabpoa` from bioconda in `conda/environment.yml` and exclude it from the uv extra. Validate with `make docker-test`. Add `edlib pyspoa pyabpoa` (bioconda/conda-forge) to `conda/environment-dev.yml` if developers use conda.

### Task 2b — `cli_run.py` and pipeline kwargs

**C2.2b** Create `src/muc_one_span/cli_run.py` with the decorators and body of cli.py:437-543 moved verbatim. Replace `@main.command()` with `@click.command()`, and add the two options just before `@record_run_status`:

```python
@click.option("--engine", type=click.Choice(["ladder", "hybrid"]),
              default=DEFAULT_SETTINGS.run.engine,
              help="Allele reconstruction engine (hybrid is experimental; default: ladder).")
@click.option("--assay", type=click.Choice(["amplicon", "genomic"]),
              default=DEFAULT_SETTINGS.run.assay,
              help="Library type used by the hybrid engine (default: amplicon).")
```

Add the parameters `engine: str = DEFAULT_SETTINGS.run.engine, assay: str = DEFAULT_SETTINGS.run.assay` after `mapping_timeout`, and pass `engine=engine, assay=assay` to `execute_pipeline`. `cli_run.py` imports only `click`, `cli_settings`, `run_status`, and `settings`; it must **not** import `muc_one_span.cli` (cycle). In `cli.py`, add `from muc_one_span.cli_run import run` to the top-level imports and `main.add_command(run)` right after the `main` group definition (cli.py:38-54). This keeps `muc_one_span.cli.run` importable (package rule: preserve import paths).

In `pipeline.py`, add kw-only `engine: str | None = None, assay: str | None = None` to `execute_pipeline`. In the `effective_run_settings(...)` call (pipeline.py:55), append:

```python
        **{k: v for k, v in (("engine", engine), ("assay", assay)) if v is not None},
```

Temporary guard, to be removed in Task 9, placed right after `write_run_configuration` (pipeline.py:88-90):

```python
    if settings.run.engine == "hybrid":
        raise click.BadParameter("the hybrid engine is not available in this build", param_hint="--engine")
```

**C2.3 CLI test (replaces the plan's Step 6):**

```python
def test_engine_option_reaches_pipeline(tmp_path) -> None:
    fq = tmp_path / "r.fastq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, ["run", "-i", str(fq), "-o", str(tmp_path / "o"),
                                        "--engine", "hybrid", "--assay", "genomic"])
    assert res.exit_code == 0, res.output
    assert (ex.call_args.kwargs["engine"], ex.call_args.kwargs["assay"]) == ("hybrid", "genomic")


def test_engine_default_comes_from_config(tmp_path) -> None:
    cfg = tmp_path / "c.json"
    cfg.write_text('{"schema_version": 1, "run": {"engine": "hybrid"}}')
    fq = tmp_path / "r.fastq"
    fq.write_text("@r\nACGT\n+\nIIII\n")
    with patch("muc_one_span.pipeline.execute_pipeline") as ex:
        res = CliRunner().invoke(main, ["--config", str(cfg), "run", "-i", str(fq),
                                        "-o", str(tmp_path / "o")])
    assert res.exit_code == 0, res.output
    assert ex.call_args.kwargs["engine"] == "hybrid"
```

After the move, check that `cli.py` is at most about 540 lines and `cli_run.py` about 130.

### Task 3

**C3.1 `synth.dupc`:**

```python
def dupc(unit: str = "X") -> str:
    """Parent unit carrying the dictionary dupC template (7C->8C; insert at 1-based 60)."""
    return next(seq for seq, (parent, name) in RD.mutated_sequences.items()
                if parent == unit and name == "dupC")
```

Add `test_dupc_template_is_classified` in `test_spans.py` or a new `test_synth.py`: `classify_sequence(allele_with_dupc, RD)["mutations_detected"]` contains `mutation_name == "dupC"` with `template_match is True`.

**C3.2 `align.py` lazy edlib** (so `rc`/`cigar_ops`/`synth` import without the extra; the pattern matches benchsim `realism._edlib`):

```python
import importlib
from types import ModuleType

HYBRID_HINT = "The hybrid engine needs the 'hybrid' extra: pip install 'muc_one_span[hybrid]'"


def _edlib() -> ModuleType:
    try:
        return importlib.import_module("edlib")
    except ImportError as exc:  # pragma: no cover - depends on installed extras
        raise ImportError(HYBRID_HINT) from exc
```

Replace each `edlib.align(` with `_edlib().align(`.

**C3.3** Only if the edlib marker is kept (see (d)), add `tests/unit/hybrid/conftest.py`:

```python
import importlib.util

collect_ignore_glob = [] if importlib.util.find_spec("edlib") else ["test_*.py"]
```

If you do this, the coverage job must run on an interpreter where edlib installs.

### Task 4

**C4.1 Drop numpy** (it avoids a new runtime dependency and the mypy/quality-env gap). In `lengths.py`:

```python
import math

def _density(lengths: list[float]) -> tuple[list[float], list[float]]:
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
```

Then `support = {c: sum(abs(x - c) <= window_bp(c, settings) for x in lengths) for c in centers}`, and use `min(range(len(peaks)), key=lambda j: dist[j])` instead of `np.argmin`.

**C4.2 Smear rule** (replaces `smear = ...` in `fit_length_model`):

```python
def _is_smear(c: float, top: float, lengths: list[float], settings: HybridSettings) -> bool:
    """A candidate >1.5 units below the major peak must stand out from its shoulders."""
    if c >= top - 1.5 * UNIT:
        return False
    w = window_bp(c, settings)
    inside = sum(abs(x - c) <= w for x in lengths)
    shoulders = sum(w < abs(x - c) <= 3 * w for x in lengths) / 2  # same width as inside
    return inside < settings.smear_min_prominence * max(shoulders, 1.0)
```

In the loop: `smear = _is_smear(c, top, lengths, settings)`, and do **not** combine it with `not far`. Keep `test_smear_is_short_product_not_allele`. Add `test_real_short_allele_amid_smear_is_kept`: a 30-unit allele with 40 reads plus a 60-unit allele with 120 reads at `smear_frac=0.3` should give 2 peaks.

### Task 5

**C5.1 `poa.get_backend(name: str)`** takes `settings.hybrid.poa_backend` and raises `ImportError(HYBRID_HINT + f" ({name}: {exc})")`. There is no fallback. The engine records `backend.name` in `summary["hybrid"]["poa_backend"]`, and `run_configuration.json` already records the setting. Test with `monkeypatch.setitem(sys.modules, "pyabpoa", None)`: the result must be an `ImportError` that mentions the extra.

**C5.2 Partial-read projection** (needed before fragments enter polishing):

```python
def project(read: str, cons: str, *, partial: bool = False) -> Columns:
    """Project read onto consensus columns; partial reads leave uncovered columns as None."""
    res = _edlib().align(read, cons, mode="HW" if partial else "NW", task="path")
    n = len(cons)
    ti = res["locations"][0][0] if partial else 0
    cols: list[str | None] = [None] * n if partial else ["-"] * n
    ...  # same CIGAR walk as global_columns, starting at ti
```

Keep `global_columns(read, cons)` as `project(read, cons)` so the Task 6/8 callers stay unchanged. In `pileup_polish`, skip `None` columns, and skip insertion votes at `ins[pos]` unless `cols[pos-1]` and `cols[pos]` are both covered (use `n` as the right sentinel for spanning reads). The majority denominator becomes the per-column covered count, not `len(reads)`. Test: 20 spanning reads plus 60 left-half fragments must reproduce the full-length consensus exactly.

**C5.3** In `test_polish.py`, move `pytest.importorskip("edlib")` above the `muc_one_span.hybrid.*` imports, or rely on C3.3.

### Task 6

**C6.1** Tests:

```python
def test_single_base_unit_difference_is_unconfirmed() -> None:
    a = ["X"] * 10 + ["C"] + ["X"] * 19        # C differs from X at 1 base
    b = ["X"] * 30
    members = _spans(a, 40, 1) + _spans(b, 40, 2)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "unconfirmed_single_site" and len(res.groups) == 1


def test_multi_base_unit_difference_counts_as_linked_sites() -> None:
    a = ["X"] * 10 + ["A"] + ["X"] * 19        # A differs from X at 4 linked bases
    members = _spans(a, 40, 1) + _spans(["X"] * 30, 40, 2)
    res = split_by_linked_sites(synth.allele(a), members, S, random.Random(1))
    assert res.basis == "linked_sites"
```

If the user prefers loci-level counting instead, merge sites closer than `UNIT` before `_linked` and invert the second test. Record the choice in the plan's decisions.

### Task 7

**C7.1** Test replacement. It only asks decidable fragments to be decided:

```python
def test_fragments_spanning_the_difference_are_assigned_correctly() -> None:
    refs = hybrid_references({"allele_1": synth.allele(A), "allele_2": synth.allele(B)}, synth.RD)
    tract_start, tract_end = (5 + 13) * 60, (5 + 13 + 21) * 60   # the 20/21-X tract in B
    wrong = decided = 0
    for truth, inner in (("allele_1", A), ("allele_2", B)):
        seq = synth.allele(inner)
        for i, start in enumerate(range(tract_start - 400, tract_start - 100, 20)):
            frag = synth.reads(seq[start:tract_end + 400], 1, err=0.02, seed=100 + i, flank_bp=0)[0].seq
            got = assign_read(rc(frag) if i % 2 else frag, refs, margin=3)
            decided += got.allele is not None
            wrong += got.allele not in (None, truth)
    assert wrong == 0 and decided == 30


def test_fragment_inside_shared_tract_is_undecided() -> None:
    refs = hybrid_references({"allele_1": synth.allele(A), "allele_2": synth.allele(B)}, synth.RD)
    frag = synth.allele(A)[(5 + 13) * 60:(5 + 13 + 15) * 60]
    assert assign_read(frag, refs, margin=3).allele is None
```

Check the tract offsets against `A`/`B` when implementing: pre-repeats 5, then 10 X, then `A A B`, then the X tract.

**C7.2** Off-target guard in `assign_read`: if `min(dist.values()) > max_error_rate * len(seq)`, return `Assignment("off_target", ...)`. `assign_reads` then gains an `"off_target"` key. Also assign when there is a single reference, subject to the same guard. Test: `"ACGT" * 400` gives `off_target`.

### Task 8

**C8.1** Type homopolymer events from the dictionary. `event_read_support` gains a `rd: RepeatDictionary` parameter. `rd.mutations[name]["changes"]` is the raw template (config.py:120-122).

```python
def homopolymer_event_run(mutation: dict, rd: RepeatDictionary, cons: str,
                          start: int, end: int) -> tuple[int, int, str] | None:
    """Consensus run changed by a templated single-base indel inside a run >= 4, else None."""
    changes = (rd.mutations.get(mutation.get("mutation_name") or "") or {}).get("changes", [])
    if len(changes) != 1 or changes[0].get("type") not in ("insert", "delete"):
        return None
    change = changes[0]
    if change["type"] == "insert" and len(change.get("sequence", "")) != 1:
        return None
    if change["type"] == "delete" and int(change["end"]) - int(change["start"]) != 0:
        return None
    pos = start + int(change["start"]) - 1
    base = change.get("sequence") or cons[max(pos - 1, start)]
    for s, e, b in _runs(cons[start:end], 4):
        if b == base and s + start <= pos <= e + start:
            return s + start, e + start, b
    return None
```

Only events for which this returns a run use the LLR path. Test: a substitution event in an X unit is scored with `kind == "competition"`, not `"homopolymer"`.

**C8.2** Status for homopolymer events (spec §5 thresholds plus strand consistency):

```python
def _hp_status(n: int, llr: float, alt_frac: float, strand_llr: dict[str, float],
               strand_n: dict[str, int], s: HybridSettings) -> str:
    if n < s.hp_min_reads:
        return "insufficient_depth"
    if llr < s.hp_llr_min or alt_frac < s.hp_min_alt_frac:
        return "not_supported"
    if any(strand_n[st] >= s.hp_min_strand_reads and strand_llr[st] < 0 for st in "+-"):
        return "discordant"
    return "supported"
```

`strand_llr[st]` is `homopolymer_llr` restricted to that strand's observations, and a strand with no reads contributes `0.0` (Review Focus #5). Tests: (a) all reads on `-` gives `supported` (the existing single-strand test plus a status assertion); (b) `+` strand LLR < 0 with ≥5 reads while the total is positive gives `discordant`.

**C8.3** Non-homopolymer events use competition, not the `% 60` rule. For each read, compare the piece `seq[t2q[start]:t2q[end]]` against the consensus window (template) and the parent unit `rd.repeats[mutation["closest_type"]]`. Count `alt` if `ed(piece, window) < ed(piece, parent)`, `ref` if greater, and `other` if equal. `supported` iff `n >= hp_min_reads`, `alt / n >= hp_min_alt_frac`, and `alt > ref`. Emit `kind: "competition"`. This is still provisional and tuned in benchmark plan Task 9.

### Task 9

**C9.1 Return type** (no side channels inside `alleles`):

```python
@dataclass
class HybridResult:
    alleles: dict[str, Any]                 # allele_1/allele_2 + ladder-compatible top-level keys only
    consensus_paths: dict[str, Path]
    block: dict[str, Any]                   # becomes summary["hybrid"]; never stored in alleles
    members: dict[str, list[tuple[str, str]]]  # allele -> (oriented seq, strand) for S10


def reconstruct_alleles(input_path: Path, output_dir: Path, rd: RepeatDictionary,
                        settings: RuntimeSettings) -> HybridResult: ...
```

Also write `hybrid_reads.json` and `hybrid_references.fa` as the plan says.

**C9.2 Pipeline branch.** Delete the Task 2 guard. In `execute_pipeline`, right after `configuration_record = write_run_configuration(...)` (pipeline.py:88-90) and **before** `igv_requested`/`check_tools` (pipeline.py:91-100):

```python
    if settings.run.engine == "hybrid":
        from muc_one_span.hybrid import reconstruct_alleles

        check_tools(["samtools"] if str(input_path).endswith(".bam") else [])
        hybrid = reconstruct_alleles(Path(input_path), out, rd, settings)
        (out / "alleles.json").write_text(json.dumps(hybrid.alleles, indent=2) + "\n")
        finish_run(out=out, input_path=input_path, rd=rd, settings=settings,
                   alleles_result=hybrid.alleles, consensus_paths=hybrid.consensus_paths,
                   vcf_paths={}, tool_versions={"poa_backend": hybrid.block["poa_backend"]},
                   configuration_record=configuration_record, report=report,
                   bam_path=None, fasta_path=out / "hybrid_references.fa",
                   annotate=partial(annotate_read_support, rd=rd, members=hybrid.members,
                                    settings=settings.hybrid),
                   extra_summary={"hybrid": hybrid.block})
        return
```

Do **not** call `annotate_selection_qc` on hybrid alleles, because it would overwrite `depth_status` with `not_assessed`.

**C9.3 `pipeline_tail.finish_run`** moves pipeline.py:152-228 verbatim. Its signature:

```python
def finish_run(*, out: Path, input_path: str, rd: RepeatDictionary, settings: RuntimeSettings,
               alleles_result: dict[str, Any], consensus_paths: dict[str, Path],
               vcf_paths: dict[str, Path], tool_versions: dict[str, str],
               configuration_record: dict[str, Any], report: bool, bam_path: Path | None,
               fasta_path: Path | None,
               annotate: Callable[[str, str, dict[str, Any]], dict[str, Any]] | None = None,
               extra_summary: dict[str, Any] | None = None) -> dict[str, Any]:
```

Import `classify_sequence`, `validate_mutations_against_vcf`, `parse_vcf_variants`, and `generate_report` **inside** the function body (keeps the patch targets in `test_pipeline_gates.py` and `test_wave1_cli.py`). `annotate(allele_key, sequence, result)` runs after classification and before `repeats.json` is written. `extra_summary` is merged into the top level of `summary`. With IGV on and `bam_path=None`, the report shows the notice and the preflight is skipped. The ladder call site passes `bam_path=bam, fasta_path=ref, annotate=None`.

**C9.4 Grouping/basis logic** (replaces the plan's `groups` loop):

```python
groups: list[tuple[list[SpanRead], str]] = []
selection = "resolved"
for peak in model.peaks:
    draft = draft_consensus(peak.members, h.n_poa, rng, backend)
    split = split_by_linked_sites(draft, peak.members, h, rng)
    if split.basis == "unconfirmed_single_site":
        selection = "unresolved_single_site"
    basis = split.basis if split.basis != "none" else ("length" if len(model.peaks) == 2 else "none")
    groups.extend((g, basis) for g in split.groups)
if len(groups) > 2:
    selection = "unresolved_max_alleles"   # never drop an allele silently
    groups = sorted(groups, key=lambda g: -len(g[0]))[:2]
if not groups:
    raise InsufficientEvidenceError("hybrid: no allele length peak passed the thresholds")
groups.sort(key=lambda g: statistics.median(m.length for m in g[0]))   # allele_1 = shorter
if any(r["reason"] in ("support_below_threshold", "max_alleles") for r in model.rejected):
    selection = "unresolved_rejected_peak" if selection == "resolved" else selection
if model.unassigned_fraction > h.max_unassigned_spanning_fraction and selection == "resolved":
    selection = "unresolved_unassigned_spanning"
```

`homozygous = len(groups) == 1 and selection == "resolved" and not residual`. Otherwise, set `homozygous = False` and `sequence_identity_status = "unresolved"` (P0 principle, calling.py:407-410). Then the duplicate `allele_2` alias (`candidate_duplicate_of`, `reconstruction_status: "not_separately_resolved"`) triggers the existing reconstruction reason (report.py:143-149).

**C9.5 `read_input`.** For FASTQ, stream records instead of calling `.read().splitlines()`. For BAM, stream `tools.run_tool_iter(["samtools", "fastq", "-F", "0x900", str(path)])` and parse it in groups of 4 lines.

**C9.6 Engine test additions** (to `test_heterozygous_dupc_sample_...`). With the fixed `synth.dupc()`:

```python
    assert "hybrid" not in alleles                               # C9.1: contract
    summary = json.loads((tmp_path / "summary.json").read_text())  # via execute_pipeline test below
```

Add `test_hybrid_pipeline_end_to_end(tmp_path)`. It calls `execute_pipeline(str(fq), str(out), None, "", 1, 10, 5.0, False, "ont", None, engine="hybrid", settings=DEFAULT_SETTINGS)` with `patch("muc_one_span.tools.check_tools")`, then asserts:

- `summary["hybrid"]["poa_backend"]` is set;
- the dupC mutation has `read_support["status"] == "supported"`, `vcf_support is False`, and `vcf_support_status == "not_applicable_read_consensus"`;
- `compute_clinical_decision(summary)["state"] == "PATHOGENIC"`;
- `load_observation(out).status` is not `"invalid_artifacts"`.

### Task 10 (reconciled with P0)

**C10.1 `clinical_gates.py`** (92 lines; room available):

```python
LOW_DEPTH_STATUSES = frozenset({"low", "insufficient"})
_DEPTH_BASIS_LABELS = {"primary_alignment_records": "primary alignments",
                       "spanning_reads": "spanning reads"}
_GENOTYPE_REASONS["residual_heterogeneity"] = (
    "residual read heterogeneity on this allele (possible unresolved mixture, chimera or mosaicism)"
)
```

In `allele_gate_reasons`, replace the depth block with:

```python
    if info.get("depth_status") in LOW_DEPTH_STATUSES:
        basis = info.get("depth_basis") or "primary_alignment_records"
        reasons.append(
            f"{label}: {info.get(basis)} {_DEPTH_BASIS_LABELS.get(basis, basis)}, below the "
            f"per-allele depth gate ({info.get('depth_threshold')})."
        )
```

This keeps the text "12 primary alignments" that test_clinical_decision.py:307,328 assert. Also show the secondary-mode fraction only when it is not None, and append `info.get("selection_detail")` when it is present. In `mutation_blockers`, when `read_support` is a dict whose status is not `supported`, append `f"read-level support {read_support.get('status')}"` instead of the VCF message.

**C10.2 `report.py`.** Line 98: `depth_assessed = any(a.get("depth_status") in ("adequate", *LOW_DEPTH_STATUSES) ...)`. Line 116: `if carrier.get("depth_status") in LOW_DEPTH_STATUSES:`. Import `LOW_DEPTH_STATUSES` from `clinical_gates`. **Do not** add a `summary["hybrid"]` branch. The hybrid sample-level blockers reach the decision through `selection_status` (C9.4) and `allele_genotype_status` (residual heterogeneity).

**C10.3 Tests** (append to `test_clinical_decision.py`; use `_gated_summary` and `state`):

```python
HYBRID_GATES = dict(RESOLVED_GATES, depth_basis="spanning_reads", spanning_reads=120,
                    length_status="consistent_with_consensus_contig",
                    allele_genotype_status="not_applicable_read_consensus", engine="hybrid")


def _hybrid_summary(**allele2: object) -> dict:
    s = _gated_summary()
    for k in ("allele_1", "allele_2"):
        s["alleles"][k].update(HYBRID_GATES)
    s["alleles"]["allele_2"].update(allele2)
    return s


def test_hybrid_resolved_is_negative() -> None:
    assert compute_clinical_decision(_hybrid_summary())["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"


@pytest.mark.parametrize("field,value", [
    ("depth_status", "insufficient"), ("selection_status", "unresolved_rejected_peak"),
    ("selection_status", "unresolved_unassigned_spanning"),
    ("selection_status", "unresolved_single_site"),
    ("allele_genotype_status", "residual_heterogeneity")])
def test_hybrid_gate_blocks_negative(field: str, value: str) -> None:
    assert compute_clinical_decision(_hybrid_summary(**{field: value}))["state"] == "INCONCLUSIVE"


def test_insufficient_depth_carrier_blocks_pathogenic() -> None:
    mut = dict(BASE, vcf_support=False, vcf_support_status="not_applicable_read_consensus",
               read_support={"status": "supported"})
    s = _hybrid_summary()
    s["classifications"]["allele_1"]["mutations"] = [mut]
    s["alleles"]["allele_1"]["depth_status"] = "insufficient"
    assert compute_clinical_decision(s)["state"] != "PATHOGENIC"


def test_hybrid_depth_message_names_spanning_reads() -> None:
    d = compute_clinical_decision(_hybrid_summary(depth_status="low", spanning_reads=12))
    assert any("12 spanning reads" in r for r in d["details"])
```

`test_clinical_decision.py` does not import `pytest` yet; add `import pytest`. `insufficient_depth` / `discordant` read-support statuses are covered through `mutation_blockers` in `test_clinical_gates.py`.

**C10.4 Evaluation.** In `models.Event`, add `read_support_status: str = field(default="unknown", compare=False)`. Change `supported` to `self.legacy_supported and self.support_status == "exact_sequence_concordance" or (self.frameshift and self.template_match and self.read_support_status == "supported")`, and leave `legacy_supported` unchanged (scoring.py:21,126 relies on it). In `artifacts.py:80-99`, validate that `read_support`, when present, is a dict with a string `status` (otherwise raise "malformed read support"), and pass `read_support_status`. `_statuses_allow_completion` (130-153) needs no change: `unresolved_single_site` and `read_consensus_low_depth` intentionally fall to ambiguous. The plan's evaluation test must use the real helper `write_valid_artifacts(tmp_path)` (test_evaluation_artifacts.py:47). It returns None and writes into `tmp_path`, so `root = tmp_path`. `load_observation(...)` returns an observation with `.status`, not `.state` (artifacts.py:259).

**C10.5 Harness forwarding.** `benchmarking.run_pipeline(..., engine: str | None = None)` appends `["--engine", engine]` only when set; the pinned argv tests at test_benchmark_tools.py:108-141 stay valid. `clinical_runner.run_case` reads `settings.get("engine")` and `settings.get("assay")` from the **hashed** settings dict (clinical_runner.py:125) and appends the options after `"--platform", "ont"` (:147-149) when present. `scripts/clinical_benchmark.py:126-131` adds `--engine/--assay` into that dict. `scripts/benchmark.py:30` gets `--engine`.

### Task 11

**C11.1** Use `DATA = Path(__file__).resolve().parents[1] / "data" / "generated" / "sample_close_51_58"`. The sample is untracked and absent in this worktree; generate it (`make generate-testdata`) or report the skip. Also skip when `samtools` is not on PATH.

**C11.2** Manifest loop: `for c in m["cases"].values(): print(c["case_id"], c["reads"], c["muconespan_platform"], sep="\t")`. Then use `--platform "$plat"` directly and drop the `hifi`/`ont` ternary.

**C11.3** PRJEB92208: add `"engine": "hybrid"` (and `"assay": "genomic"` for the WGS libraries) to the clinical harness settings, with a new output root. Record `poa_backend` and the extra's package versions in `hybrid-engine.json`.

---

## (c) Hybrid output contract against the P0 gates

What the engine must emit so that hybrid calls are gated the same way as ladder calls. "Task" is where it is produced.

| Field (where) | P0 reader | Hybrid value | Task |
| --- | --- | --- | --- |
| `alleles.allele_N.depth_status` | `report.py:98,116`; `allele_gate_reasons` | `adequate` (≥ `depth_adequate_spanning`), `low` (≥ `depth_low_spanning`), `insufficient`. Readers must accept `insufficient` (C10.1/C10.2). | 9 (+10) |
| `depth_basis`, `depth_threshold`, `spanning_reads` | `allele_gate_reasons` message | `"spanning_reads"`, `h.depth_adequate_spanning`, and the count. Spec's "or ≥40 assigned" is deferred; record the deviation. | 9 |
| `selection_status` (+ `secondary_mode_fraction: None`, `selection_detail`) | `allele_gate_reasons` (`startswith("unresolved")`) | `resolved`, `unresolved_rejected_peak`, `unresolved_unassigned_spanning`, `unresolved_single_site`, or `unresolved_max_alleles`. Sample-level conditions are stamped on **both** alleles (C9.4). | 9 |
| `length_status`, `reference_length` | `allele_gate_reasons` | `"consistent_with_consensus_contig"` with `reference_length = length`. The length is derived from the reported consensus by construction. `length_basis: "spanning_peak"`. | 9 |
| `allele_genotype_status`, `heterozygous_sites` | `allele_gate_reasons` via `_GENOTYPE_REASONS` | `"residual_heterogeneity"` if `residual_sites` is non-empty (gated, C10.1), else `"not_applicable_read_consensus"`. `heterozygous_sites = residual_sites`. `variant_filter: None`. | 9 (+10) |
| `phase_status`, `independent_haplotype_evidence` | `report.py:130-142` reconstruction reasons | `phased` / `no_informative_heterozygosity` / `unresolved_single_site`; `True` only for `length`/`linked_sites`. | 9 |
| `alleles.homozygous`, `sequence_identity_status`, `allele_2.candidate_duplicate_of`, `reconstruction_status` | `report.py:130-149`; `evaluation/artifacts.py:229-240` | `homozygous=True` only when the call is resolved, there is one group, and no residual sites. Otherwise `False` + `"unresolved"`, and the alias `not_separately_resolved` gets the INCONCLUSIVE reason. | 9 |
| No `hybrid` or `_members` key inside `alleles` | `artifacts.py:229` | Top-level `summary["hybrid"]` only. | 9 |
| `classifications.*.mutations[].template_match`, `mutation_name`, `frameshift`, `localization_status` | `mutation_blockers` | Unchanged; produced by `classify_sequence` (classify.py:385-397, 524-528). The synthetic dupC must use the dictionary template (C3.1). | 3/9 |
| `mutations[].vcf_support`, `vcf_support_status` | `mutation_supported`/`mutation_blockers` | `False`, `"not_applicable_read_consensus"`. | 9 (`annotate`) |
| `mutations[].read_support` | `mutation_supported` (`status == "supported"`) | `{kind, n, alt, ref, other, alt_frac, strand_alt_frac, llr?, strand_llr?, status}` with status ∈ {supported, insufficient_depth, discordant, not_supported, not_localized}. | 8/9 |
| `mutations[].support_status` (display only) | report banner "Evidence:" | `"read_level_" + status`. | 9 |
| `run_status`, `execution_status` | `_execution_context` | Same tail as the ladder (C9.3). | 9 |

**`read_support.status == "supported"` producer contract** (spec §3 S10, §5; C8):

- **Homopolymer event.** This is a dictionary template that is a single-base indel inside a run ≥4 (dupC is 7C→8C). Status is `supported` iff all of the following hold:
  - `n >= hp_min_reads` (20) reads assigned to the carrier allele cover the run;
  - the stutter-aware LLR `>= hp_llr_min` (10). The background is per sample and per strand, taken from the same-base runs of length n0 on the same consensus, excluding the event run; P(obs | n0+1) is the n0 profile shifted by +1 (hp_model.py:56-66, 130-140);
  - the fraction of reads with observed run length ≥ n0+1 is `>= hp_min_alt_frac` (0.30);
  - no strand with `>= hp_min_strand_reads` (5) reads has a negative strand LLR. Otherwise the status is `discordant`.

  A strand with zero reads never fails the event. With `n < 20`, the status is `insufficient_depth` and the event is blocked, never negative.
- **Other templated events.** Parent-vs-template edit-distance competition per read (C8.3). Status is `supported` iff `n >= 20`, `alt/n >= 0.30`, and `alt > ref`. This rule is provisional and tuned on the benchmark dev split.
- **Carrier gate still applies.** A `supported` event on an allele with `depth_status` in {low, insufficient} is blocked (C10.2), matching the ladder.

---

## (d) Dependencies

PyPI JSON, queried 2026-09-24:

| Package | Latest | Wheels | 3.10–3.13 | 3.14 | Notes |
| --- | --- | --- | --- | --- | --- |
| edlib (MIT) | 1.3.9.post1 (2024-09) | manylinux2014 x86_64, musllinux, macOS x86_64/universal2 | wheels | **no wheel**; the sdist **built and imported** on local CPython 3.14.4 (Linux x86_64, system gcc) | No linux aarch64 wheel. |
| pyspoa (MIT) | 0.3.2 (2025-12) | manylinux_2_28 x86_64 + aarch64 only | wheels | wheel; import and `poa(algorithm=1)` verified on 3.14 | **No macOS wheels** (the sdist needs cmake and a C++ compiler). |
| pyabpoa (MIT) | 1.5.7.1 (2026-09-17) | **sdist only** | build from source (unverified on 3.10–3.13) | built and verified on 3.14 locally | Needs a C compiler and zlib. The Docker builder (micromamba) has no compiler. Bioconda has binaries. |

Recommendations:

1. **Keep the edlib requirement identical in `bench` and `hybrid`.** Recommended: remove the `python_version < "3.14"` marker in **both** extras (coordinate with `feat/benchsim`). The sdist builds on 3.14 with the compiler that GitHub `ubuntu-24.04` runners have, and coverage is measured on 3.14 (test.yml:104-109). If the benchsim marker must stay, `hybrid` must use the same marker; then either move the coverage job to 3.13 or accept that hybrid tests are skipped (C3.3) on 3.14. Skipping would likely push coverage below 80% (about 89.6% → about 79–80%, estimated from the current `coverage.xml`: 5 406 lines at 91.9% and 1 944 branches at 83.1%).
2. The POA backend is an explicit setting (`hybrid.poa_backend`, default `pyabpoa`, because the prototype evidence uses pyabpoa), with no silent fallback. pyspoa is the alternative backend. A backend switch requires re-running the Task 11 regressions.
3. Drop `numpy` from the extra (C4.1).
4. Add `--extra hybrid` to both `UV_TEST` and `UV_QUALITY`, and add mypy overrides for `edlib`, `pyabpoa`, and `spoa`. Update the Dockerfile builder (compiler or bioconda pyabpoa) and the conda env. Run `make security-check` and `make build-check` (the wheel must not contain test data). `requires-python` stays `>=3.10`.
5. Prototype references exist (outside Git, not copied, and no in-house outputs were read):
   - `poa-prototype/proto.py` (26 kB) and `poa-prototype/hetsplit.py` (251 lines);
   - `homopolymer/hp_model.py` and `homopolymer/hp_extract.py`;
   - `inhouse/hybrid/assign_test.py` and `inhouse/hybrid/hybrid.py`.

   `hp_model.py` imports `scipy.stats.binom` (log10p only). Do **not** port that part; the LLR does not need scipy.

---

## (e) Recommended order

Pure-library tasks need no pipeline, CLI, or report changes and can run early or in parallel once 2a lands. They are marked ★.

1. **2a** Settings + extra + tooling (C2.1, C2.2, (d)). It unblocks everything.
2. **10a** ★ Gate reconciliation: C10.1, C10.2, C10.3. These are pure summary-dict tests. Doing it first means the engine is built against the final gates, and it fixes the latent `insufficient` gap for any producer.
3. **3** ★ Synth + align + spans (C3.x).
4. In parallel after 3:
   - **4** ★ lengths (C4.x);
   - **5** ★ POA + polish (C5.x);
   - **7** ★ assign (C7.x).
5. After 5:
   - **6** ★ phase (C6.1);
   - **8** ★ evidence (C8.x; needs `polish._runs`/`read_run_length`).
6. **2b** CLI move + `--engine/--assay` + pipeline kwargs + temporary guard (independent; any time before 9).
7. **9** Engine + `pipeline_tail` + branch (C9.x). Remove the guard.
8. **10b** Evaluation `Event.read_support_status`, loader validation, harness forwarding (C10.4, C10.5).
9. **11** Integration + regressions (C11.x; the outputs stay outside Git).
10. **12** Docs/CHANGELOG (the Task 12 row in (a)).

File-size headroom after the plan plus these corrections (estimates):

| File | Now → estimated |
| --- | --- |
| `cli.py` | 638 → ~535 |
| `cli_run.py` | new, ~135 |
| `settings.py` | 432 → ~510 |
| `pipeline.py` | 228 → ~175 |
| `pipeline_tail.py` | new, ~120 |
| `report.py` | 430 → ~432 |
| `clinical_gates.py` | 92 → ~115 |
| `evaluation/artifacts.py` | 263 → ~272 |
| `benchmarking.py` | 317 → ~320 |
| `clinical_runner.py` | 228 → ~232 |
| `hybrid/*.py` | each < 220 |
| `test_clinical_decision.py` | 330 → ~390 |
| `test_runtime_settings.py` | 334 → ~375 |

All stay under 649 lines. **Do not touch `clinical_scoring.py` (626) or `alleles.py` (641).** Neither needs changes for this plan.
