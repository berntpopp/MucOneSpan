# Independent experiment-runner and benchmark handoff audit

2026-09-14. Read-only production audit of `experiments.py`,
`scripts/generate_testdata.py`, `scripts/benchmark.py`, `scripts/evaluate.py`,
benchmarking/evaluation adapters, and the maintained JSON design. Only this note
was authored by the auditor. Temporary fixtures and mocked simulator/caller tools
were used; no real generation, historical panel rerun, or scientific calling ran.

## Verdict

Two reproducible provenance bugs were reported promptly and fixed by the owning
agent. Independent reruns confirm both fixes. No remaining blocking finding was
identified in the bundled-dictionary generation → inventory → benchmark →
evaluation workflow exercised here. Custom-dictionary propagation remains a
scope limitation described below, not a validated end-to-end workflow.

Final inspected `experiments.py` SHA-256:
`6d531990e799fc8dbac776940fe57b22bf87eb6a1a380b47c0de00cb31327087`.

## F1: final simulator command could change the recorded model — fixed

The initial implementation verified config/model artifact hashes before each
command but not after the final `reads` command. A wrapper around the existing
unit-test `fake_simulator` wrote a different model after producing its final
read fixture. The case returned `completed`, while its recorded model SHA-256
was `9372c470eeadd5ecd9c3c74c2b3cb633f8e2f2fad799250a0f70d652b6b825e4`
and the final model hash was
`2e00766a86032c1526fd88d4e8289292e33c04ba395c576dea0db095d4179667`.

The owner added post-command verification. Repeating the same fixture now yields
`execution_failed` with `artifact changed during execution`, preserving the
failed inventory entry rather than publishing an input as completed.

## F2: recorded executable could differ from the executed one — fixed

The initial executable provenance used `shutil.which` on the ordinary PATH,
whereas `run_tool` removes Python virtualenv bin entries. A temporary PATH with
`project/.venv/bin/muconeup` before `external/bin/muconeup` caused the manifest to
record the former while the real tool abstraction would resolve the latter.

The owner now resolves the executable using the tool abstraction's cleaned PATH
and pins that absolute path. The independent repeated fixture verifies that the
manifest executable, version command, and every planned simulator command all
refer to the same external executable. Explicit relative executable paths remain
bound to the invocation directory before the config-directory cwd is applied.

This records executable/config/model hashes and successful version output. It
is not a complete Python package/environment lock; the manifest explicitly says
so. Transitive interpreter/module provenance should not be inferred from a
console-entry-script hash alone.

## Validation and legacy panel

The schema rejects duplicate JSON keys, nonfinite JSON values, nonobject roots,
unknown fields, invalid platform/config sections, unsafe names, duplicate cases,
boolean/nonpositive integer parameters where inappropriate, malformed diploid
lengths, and invalid/duplicate mutation targets. Platform-specific config paths
resolve relative to the design; explicit config overrides design entries, which
precede the environment fallback. Dry-run performs no tool execution or writes.

I parsed the baseline `d8390b3:scripts/generate_testdata.py` with Python AST and
compared every literal case against `examples/development-experiment.json`.
All **29 cases match exactly**, including names, order, 26 HiFi/three ONT platform
assignment, both lengths, every seed, mutation, target, and requested coverage.
The lower-template case remains 50; others remain 200.

Both simulator platform keys and explicit model-file handling were inspected
against the existing sibling MucOneUp source: HiFi maps to the simulator's
`pacbio` platform and `pacbio_params.model_file`; ONT maps to `ont` and
`ont_amplicon_params.model_file`. Temporary named `hifi.model`/`ont.model` files
were selected correctly and recorded separately. No model was installed or run.

## Actual CLI interoperability fixture

The audit invoked the actual `main` parsers for benchmark and evaluation, using
three explicit cases: successful HiFi, successful ONT, and an injected simulation
failure. The simulator fixture emits strict independently constructed diploid
truth plus FASTQ records; caller execution was intercepted deliberately.

- Generation statuses were `completed`, `completed`, `execution_failed`.
- All three rows were present in the generated inventory.
- A copy enriched with explicit per-platform Clair3 `model` paths was passed via
  benchmark's actual `--expected-samples` flag, with `--data-dir`, `--output-dir`,
  and `--threads 2`.
- The intercepted caller received exactly the successful sample's absolute FASTQ,
  correct `hifi`/`ont` platform, matching named Clair3 model, and two threads.
  It received no truth FASTA. The failed-generation sentinel path was not called.
- The two intercepted calls deliberately returned failure; benchmark retained
  statuses `execution_failed`, `execution_failed`, `not_attempted` and exited 1.
- Evaluation used the same explicit inventory with `--truth-root` and retained
  three rows: two `execution_failed` and one `invalid_truth`, exiting 1.
  Totals exposed three expected samples, two valid-truth samples, one invalid
  truth, and four missing truth alleles. Invalid truth cannot contribute unknown
  allele/event denominators; it remains explicitly visible at the sample level.

This verifies failure accounting and argument/input selection without claiming
successful real caller reconstruction. Usable-record counts preserve colliding
names and differ explicitly from requested template counts. Input discovery
rejects ambiguity; it does not accidentally choose a truth FASTA.

## Scope limitation: custom dictionaries

`design.repeat_dictionary` governs generation truth validation and is recorded
in the generation manifest. It is not forwarded through the generated inventory
to the maintained benchmark CLI, which does not accept a runtime config, or the
standalone evaluator, which loads the bundled dictionary. The bundled default
workflow is covered above. A custom dictionary requires deliberate coordinated
caller/evaluator configuration through the relevant APIs or future CLI support;
this audit does not claim that workflow is already automatic.

## Executed focused tests

```sh
uv run --locked --all-extras pytest -q tests/unit/test_experiments.py tests/unit/test_benchmark_tools.py tests/unit/test_evaluation_cli.py --no-cov
```

**61 passed, zero skipped**, after the owner's provenance fixes landed.
Independent temporary-fixture scripts additionally established exact 29-case
baseline parity, reproduced then verified both provenance corrections, and
asserted actual CLI argument selection and complete three-row failure accounting.
