# Legacy benchmark tool repair

## Scope and diagnosis

`scripts/benchmark.py` manually called pipeline internals, discarded the full
CLI artifacts, selected the first BAM, hard-coded four HiFi threads and returned
success for a missing data root. `scripts/batch_analyze.py` repeated that input
selection problem, invoked an unlocked nested `uv`, accepted failed samples in a
successful process, and called substring or wrong-event detections `TP_partial`.

The repair keeps the existing benchmark flags and batch positional arguments.
Both scripts use one in-process Click runner around the real `muconespan run`
command, with timed wrappers around the five existing stages. JSON inventory
entries can select `input`, `platform`, `model`, `threads`, `input_dir`,
`result_dir`, and `truth_dir`; supported discovery suffixes are BAM, FASTQ, FQ,
and their gzip forms. Discovery requires exactly one input. A declared per-sample
platform cannot conflict with the command-line platform.

## Output migration

`benchmark_results.json` remains an array and adds execution status, exact CLI
arguments, selected input, platform, model, threads, stage/wall timings, and
errors. Per-sample `measurement.json`, `cli.log`, and root `measurements.json`
make failed invocations auditable.

`batch_results.json` changes from the unversioned legacy array of
`TP`/`TP_partial`/`FN`/`FP`/`TN` rows to evaluation schema version 1. Its
`schema_migration` object names the removed legacy categories and points callers
to `samples`, `totals.metrics.event_{tp,fn,fp}`, and denominator-aware aggregate
ratios. This is deliberately a breaking data-schema correction: no legacy status
is emitted because it cannot represent wrong-site calls, positive-sample extras,
allele assignment ambiguity, or failed-sample denominators truthfully.

## Validation log

Focused tests were written and observed failing before production edits: the
first run reported 7 failures for the absent module and unchanged script entry
points. Two additional red runs proved that stale completed sidecars/nonzero
callers and simulator platform metadata were not yet handled.

- `uv run --locked --all-extras pytest tests/unit/test_benchmark_tools.py
  --no-cov -q`: 11 passed.
- Focused benchmark/evaluator/status regression selection: 60 passed on the
  final run.
- `make type-check file-size workflow-check`: passed; mypy checked 39 source
  files and the file-size gate checked 95 files.
- Both script `--help` invocations: exit 0.
- `make test-fast`: stopped at an unrelated concurrent
  `test_trim_exposes_failed_anchor_fallback` failure after 196 passes. The new
  benchmark tests had all passed in that run.
- `make quality`: focused lint passed, then the repository format gate reported
  pre-existing/concurrent formatting in `tests/unit/test_consensus.py` and
  `tests/unit/test_variant_support.py`, outside this task's ownership.

No external caller or expensive generated-data run was executed. The full CLI
path is exercised with its external tool check mocked at the boundary.
