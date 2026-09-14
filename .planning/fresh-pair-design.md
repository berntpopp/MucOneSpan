# Fresh paired pipeline driver

Status: prepared only. No final data or scientific caller execution. Ignored
`fresh_generation/paired.py`, `paired_worker.py`, and `paired_check.py` implement
selection, subprocess execution, and deterministic selection checks. The actual
two-source smoke used a nonexistent input, exiting in Click argument validation
before tools/calling. No thresholds, calling defaults, or scientific settings
were tuned.

## Fixed selection and invocation order

The explicit generation manifest's `expected_cases`, not successful directories,
defines all 32 main inputs. Every case adds minority-H1, minority-H2, and partial
inputs, for **128 primary pairs / 256 primary caller invocations**. Failed or
missing generations remain inventory entries with not-attempted measurements
and truth directories; invalid/unavailable truth remains explicit at evaluation.
There is no directory discovery, survivor filtering, or scientific selection.

The four timing-repeat main cases are fixed as:

- `fresh_hifi_normal_60_61`
- `fresh_ont_normal_60_61`
- `fresh_hifi_dupc_120_140`
- `fresh_ont_dupc_120_140`

Each receives repeat 2 and repeat 3 in addition to its main invocation. This adds
eight pairs, making **136 total pairs / 272 total invocations**. Timing repeats
have separate output IDs and remain excluded from `primary_inventory.json`.
They cannot inflate primary accuracy denominators. Main and perturbed inputs
from a common simulation remain dependent observations.

Each pair runs sequentially with no concurrency. Even index is baseline first,
odd index candidate first. Primary order is generation-manifest case order with
main/minority-H1/minority-H2/partial inside each case; timing rounds 2 and 3 follow
in case order. There are no retries or overwrite/resume modes.

## Source and environment isolation

Baseline must be original commit `d8390b3c244ef8f3240af74b92db12b50dfc77d1` at
`/home/bernt-popp/development/MucOneSpan`, with no runtime source, pyproject or lock
diff. Existing docs-only user edits are allowed. Candidate is the production
worktree. Both callers execute in **separate subprocesses using the exact same
worktree uv Python**, same platform model, same external tool PATH and four
nominal threads. Each worker puts only the chosen source first on `sys.path` and
sets `PYTHONPATH` accordingly. It checks and records the actual imported paths of
CLI and all five scientific stage modules against that selected source tree.

The worker directly invokes the selected source's unmodified Click CLI and wraps
existing stage functions only for timing. This avoids the candidate maintained
benchmark runner's top-level imports of evaluation modules unavailable in the
baseline package. It does not alter baseline import search paths to borrow
candidate scientific modules. Simulator truth, source labels, and haplotype maps
are never supplied to the worker or caller.

Prepare records complete chosen source/resource hashes (excluding runtime
indexes/bytecode), interpreter and timeout executable, minimap2/samtools/bcftools/
Clair3-launcher hashes, model-file hashes, and paired driver hashes. Source
snapshots are rechecked before launch and each invocation, including additions
or removals. This is not a full shared-library/Python-environment dependency
closure; the locked Python environment and existing tool inventory remain
additional provenance. Each input hash is rechecked before both paired runs.

## Timing and failure evidence

Each run records CLI arguments, module import origins, interpreter, platform,
PATH, model, thread count, exit code, exception, CLI log, stage times and wall
time. Parent wall time separately includes worker process startup/import costs.
Workers record `resource.getrusage` before/after for self and waited children,
CPU deltas, and their separate Linux maxRSS values in KiB. Those peaks are not
summed and are not represented as concurrent whole-process-tree RSS. Resource
measurement after an externally killed worker can be unavailable, explicitly
retaining its incomplete measurement and parent failure record.

The bundled reference is reused; there is no new full ladder build. Reference
lookup is timed separately. Per-allele `samtools faidx` extraction/index commands
are timed and retained separately but are already included within calling time,
so they must not be added again. Implicit minimap2 indexing is included within
mapping time; it is not falsely presented as separately measured reference build.

All subprocesses use the existing `run_tool` abstraction with argument lists.
GNU timeout (default frozen 1800 seconds; TERM then KILL after ten seconds)
provides visible process failure and process-group termination. Timeout exits
124/137 yield timeout records even when stale running sidecars survive.
Inventories contain every expected run before execution and update after each
attempt. External failures produce nonzero driver exit and retain denominators.
Output roots and worker directories must not exist.

Each source output contains `measurements.json`, `inventory.json` (136 rows),
and `primary_inventory.json` (128 rows). Entries have absolute `truth_dir`, input,
result_dir, provenance and `run_record`; these work with `evaluate_inventory`
when the parallel measurement list is supplied, or with `scripts/evaluate.py`.
Truth is consumed only in this later offline evaluation. The run root retains
the copied generation manifest, coordinator freeze, paired settings and actual
ordered subprocess commands.

## Freeze and launch commands

Run from the production worktree after all source edits are frozen. `prepare`
itself invokes no caller and must use a fresh settings directory:

```bash
uv run --locked --all-extras python \
  tests/results/production_validation_20260914/fresh_generation/paired.py \
  --caller /home/bernt-popp/development/MucOneSpan/.worktrees/production-validation \
  prepare \
  --baseline /home/bernt-popp/development/MucOneSpan \
  --candidate /home/bernt-popp/development/MucOneSpan/.worktrees/production-validation \
  --tool-bin /home/bernt-popp/miniforge3/envs/env_clair3/bin \
  --prepared tests/results/production_validation_20260914/fresh_generation/paired_prepared_final \
  --timeout-seconds 1800
```

The command prints `paired_settings.json`'s immutable hash. Incorporate that hash
and the paired-script hashes into the coordinator freeze. After generation
finishes, fingerprint its completed generation manifest and run:

```bash
uv run --locked --all-extras python \
  tests/results/production_validation_20260914/fresh_generation/paired.py \
  --caller /home/bernt-popp/development/MucOneSpan/.worktrees/production-validation \
  run \
  --prepared tests/results/production_validation_20260914/fresh_generation/paired_prepared_final \
  --settings-sha256 PAIRED_SETTINGS_SHA256 \
  --freeze COORDINATOR_FREEZE_JSON \
  --freeze-sha256 COORDINATOR_FREEZE_SHA256 \
  --generation tests/results/production_validation_20260914/fresh_generation/data/generation_manifest.json \
  --generation-sha256 COMPLETED_GENERATION_MANIFEST_SHA256 \
  --output tests/results/production_validation_20260914/fresh_generation/paired_final
```

There is no implicit permission or automatic launch. The driver checks freeze
bytes but leaves the coordinator's freeze schema/approval policy to the root.

## Executed preparation checks

- `paired_check.py`: exact 128 primary / eight timing pairs; all absent generation
  entries retained; duplicate and incomplete generation manifests rejected.
- Ruff check and format check: all three paired scripts pass; each under 650 lines.
- Two real worker subprocesses with nonexistent development input: both exit 2,
  no external scientific tools or calling invoked. Module origins are baseline
  `/MucOneSpan/src` versus candidate `/MucOneSpan/.worktrees/production-validation/src`;
  both use `.worktrees/production-validation/.venv/bin/python3`. Full commands,
  logs, resources and measurements are under `fresh_generation/paired_smoke_badinput`.
- A 20ms timeout probe of a one-second Python sleep raises visible exit 124 through
  `run_tool`; exact probe retained alongside the smoke.
- `prepare` was exercised into `paired_prepared_development`, hash
  `d83adcda914e3050548df383378106eef603b02ec81badda6cf8e5cbf78f49cf`.
  This captures development source, **not final frozen settings**; final preparation
  must be repeated after remaining coordinator edits using the command above.

No complete pipeline performance, model invocation, scientific result, or
timeout of a full bioinformatics process tree is claimed from these smoke checks.
