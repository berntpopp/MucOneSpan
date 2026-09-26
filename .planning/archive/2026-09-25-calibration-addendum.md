# Addendum (2026-09-25): config-driven thresholds and calibration commands

Owner directive: everything is config-driven, with no hardcoded thresholds, magic
numbers or variables. Tunable thresholds get proper calibration commands.

This adds **Task 15** to `.planning/archive/2026-09-24-hybrid-engine-plan-v2.md` and a
constraint to every task (see `global-constraints.md` in the SDD workspace).

## Constraint (all tasks, both plans)

- Every tunable is a validated settings field with a documented default: a
  threshold, width, margin, multiplier, floor or divisor. It goes in
  `HybridSettings` or another `RuntimeSettings` section, or in benchsim
  config/profile JSON. Defaults preserve current behaviour.
- Domain constants come from bundled data. For example, the repeat-unit length
  comes from the repeat dictionary.
- Only structural literals (indexing, empty defaults) and named format
  constants (for example the Phred offset) remain.
- The final whole-branch review sweeps every module the branch touches for
  literals.

## Task 15: calibration commands

**Depends on:** the benchmark branch (`feat/benchsim`) merged to main and the
hybrid branch rebased on it, because calibration needs the benchsim truth,
`run` and `evaluate`. It runs after hybrid Tasks 1–14.

### 15a: `muconespan settings` (MucOneSpan CLI)

- `muconespan settings show [--section NAME] [--config FILE]` prints the
  effective settings as schema-1 JSON. The output loads with `--config`.
  Without `--config` it prints the defaults.
- `muconespan settings validate FILE` runs the strict loader and reports the
  first error. It exits non-zero on an invalid file.
- `examples/runtime-settings.json` is regenerated from `settings show`, and a
  test keeps it in sync.

### 15b: `benchsim calibrate` (scripts/benchsim.py)

```
benchsim calibrate --split {dev,val} --engine {ladder,hybrid} \
    --grid GRID.json [--config BASE.json] [--jobs N] [--name NAME]
```

- `GRID.json` maps dotted settings keys to value lists, for example
  `{"hybrid.smear_test_window_frac": [0.15, 0.25, 0.35]}`. It may also use
  `{"min", "max", "step"}`. Any `RuntimeSettings` field is allowed, so ladder
  thresholds calibrate the same way. Unknown keys and invalid values fail
  before any run. Every point is validated through the strict settings
  loader.
- The command refuses the sealed `test` split.
- Each grid point is written as an overlay config (base plus point,
  `schema_version` 1) named by its content hash. It runs through the existing
  `run` machinery with `--config`, then `evaluate`. Points are idempotent and
  resumable. Outputs go under
  `<out-root>/calibration/<split>/<name>/<point-hash>/`, outside Git.
- A `calibration.json` manifest records the grid, split, design seeds, base
  config hash, engine and tool versions, and per-point status.

### 15c: `benchsim calibrate-report`

```
benchsim calibrate-report --split {dev,val} --name NAME --objective OBJECTIVE.json
```

- `OBJECTIVE.json` declares the selection rule as data, with no defaults
  hardcoded in code. It has:
  - constraints such as `clinical_false_negative: {"max": 0}` or
    `smear_ambiguous_rate: {"max": 0.1}`;
  - a lexicographic `rank` list, for example
    `["-per_allele_exact", "inconclusive_rate"]`.
- Output is a ranked table (`calibration-report.md` and `.json`) with every
  point's metrics and cluster-bootstrap CIs, reusing the Task 8 statistics. It
  also writes `recommended-config.json`. That file is loadable by
  `muconespan run --config` and carries provenance (grid, objective, split,
  hashes).
- Confirmation step: `benchsim calibrate --split val --grid` with a
  single-point grid built from `recommended-config.json`. The report shows the
  dev→val shift. Shipped defaults change only in a separate, reviewed commit
  that cites the calibration report.

### 15d: fast stage calibration (length model)

- `benchsim calibrate --stage lengths` fits only the hybrid length model on
  each case's spanning reads. It scores peaks against the read truth: allele
  count, allele lengths, false/missed alleles and `smear_ambiguous` rate. No
  consensus or calling is run. This makes smear/peak threshold sweeps take
  seconds.

### Tests and docs

- Unit tests for each subcommand use synthetic manifests and mock engine runs:
  grid expansion, validation errors, refusal of `test`, resume, objective
  ranking and the recommended-config round trip.
- Docs: `docs/benchmark.md` gets a "Calibration" section, and
  `docs/guides/configuration.md` gets a `settings show/validate` section.

## Follow-up proposals (2026-09-25)

Each proposal becomes its own task after Task 15. Each ships only if the benchmark (standard set, dev/val) shows a gain:
- **P1:** per-read repeat typing as a second, independent structure source. Vote the unit letter per repeat position within each length peak, and cross-check the POA consensus.
- **P2:** per-repeat read support with an "any length change" count, and thresholds from a per-position background.
- **P3:** an explicit "unknown change at repeat N" state for events outside the dictionary. It is gate-relevant, so the call becomes INCONCLUSIVE with a located reason.
- **P4:** QC fields for the PCR length-bias estimate and the longest detectable allele.
- **P5:** per-allele re-calling. Align each allele's assigned reads to its own consensus. Then call residual base-level differences with a deep-learning caller: Clair3 by default, since it is already a dependency, or DeepVariant (BSD-3) as an optional backend. A homozygous call against the allele's own consensus corrects a consensus error. Heterozygous calls feed `residual_heterogeneity` as model-based evidence. Tools run through `tools.py`, and every threshold is a setting.
- **P6:** local haplotype re-vote for residual single-base consensus misses (a pileup split where one allele is counted partly as a substitution and partly as a neighbouring insertion). This is the remaining cause of the simpanel 77/80 sequence-exact result after Task 13b.

## Calibration requirements carried from Task 13/13b (controller rulings, 2026-09-25)

- `event_max_alternative_frac` (0.25) is relaxed only after the homopolymer background accounts for run length. Simulated HiFi shows about 31% deletion stutter at the true C8 run against about 8% at C7 runs. Calibrate on the MucOneUp 0.46.0 (v4) benchmark data.
- The minor read groups in generated negatives are simulator artefacts: systematic errors in lower-quality reads on one strand. Use them as input for `het_*` and `link_phi_min`.
- `benchmarks/clinical/prjeb92208/hybrid-engine.json` is refreshed at the end of Task 15, after any change to the defaults.
