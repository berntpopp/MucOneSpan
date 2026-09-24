# MucSim-Bench: realistic simulated benchmark

MucSim-Bench is the simulated benchmark used to compare MucOneSpan calling
engines (for example the current `ladder` engine against a candidate) on MUC1
VNTR structures with known truth. It designs stratified cases, simulates truth
haplotypes and reads with [MucOneUp](https://github.com/berntpopp/MucOneUp),
runs the engines, scores the results and applies a pre-registered decision
rule. Everything is driven by `scripts/benchsim.py`.

!!! warning "Simulation is not clinical validation"
    Results describe simulated reads only. The realism report measures how far
    the simulated reads are from public real data (the *sim-to-real gap*); it
    does not remove that gap. Do not quote MucSim-Bench numbers as clinical
    performance.

## Prerequisites

| Requirement | Used by | Notes |
| --- | --- | --- |
| MucOneUp **>= 0.45.0** (`muconeup`) | `generate` | Called as an external executable, never imported. The version is checked before any case is generated. |
| MucOneUp `config.json` | `generate` | `--muconeup-config` or `$MUCONEUP_CONFIG`. |
| MucOneUp read profiles | `generate` | `--muconeup-profiles <MucOneUp>/muc_one_up/data/read_profiles`. |
| pbsim3 (`pbsim`) and `ccs` | `generate` | MucOneUp calls them for read simulation, through the `tools` section of its `config.json`. |
| minimap2, samtools, bcftools, Clair3 + models | `run` | On `PATH`, as for `muconespan run`; ONT and HiFi models via `--model-ont` / `--model-hifi` or `$CLAIR3_MODEL_ONT` / `$CLAIR3_MODEL_HIFI`. |
| edlib (`bench` extra) | `realism` | `uv sync --extra bench`, or `make dev` (all extras). No edlib wheel exists for Python 3.14 yet; `realism` then fails with a clear error. |

## Profiles

| Profile | MucOneUp base profile | Status |
| --- | --- | --- |
| `ont_amplicon_r10` | `ont_r10_sup_amplicon_v1` | Calibrated against public PRJEB92208 amplicon aggregates. |
| `ont_genomic_targeted` | `ont_r10_genomic_v1` | Calibrated against public PRJEB92208 whole-genome aggregates. |
| `hifi_amplicon` | `hifi_amplicon_v1` | **Uncalibrated**: no public PacBio MUC1 data. Realism has no target section for it, so HiFi results show relative engine behaviour only. |

Each design layers artefact, error and PCR levels onto the base profile
(`smear`, `chimera`, `error`, `pcr`); the variant JSON and the SHA-256 of the
base profile are recorded for every case. The level values are settings (see
[Configuration](#configuration)).

## Configuration

Every tunable number of the benchmark is a validated setting in
`muc_one_span.benchsim.bench_config`. The defaults below apply unless
`--bench-config FILE` (given before the subcommand) names a JSON file with
`"schema_version": 1` and any subset of the sections. Unknown sections or fields,
duplicate keys and out-of-range values are rejected. The SHA-256 of the
effective settings is written to every `case.json` (`bench_config_sha256`),
`realism.json` and `report.json`.

```bash
python scripts/benchsim.py --bench-config my-bench.json design --split dev --n 30
```

| Setting | Default | Meaning |
| --- | --- | --- |
| `design.split_sizes` | dev 300, val 300, test 800, stress 100 | Cases per profile when `--n` is omitted |
| `design.normal_fraction` | 0.35 | Minimum share of normal cases per profile |
| `design.length_min` / `length_max` | 20 / 130 | Allele length range (repeat units) |
| `design.depths` | amplicon 5-2000, genomic 3-80 | Target spanning depth levels per profile |
| `design.delta_ranges` | `0_identical` ... `>20` (21-90) | Length difference classes |
| `design.compositions` | markov 0.75, real_derived 0.20, rare_units 0.05 | Structure source weights |
| `design.position_fraction` | 0.1 | Leading/trailing fraction for `first10` / `last10` |
| `design.pcr_levels`, `error_levels`, `chimera_levels` | all levels; chimera 0.01, 0.05 | Artefact factor levels |
| `design.smear_levels` | dev/val/test 0.05, 0.25; stress 0.5 | Smear rate levels per split |
| `design.offpeak_share_cap` | 0.5369 | Cap on max smear + max chimera for regular splits (real PRJEB92208 maximum of `span_off_gt1unit_frac`) |
| `design.offpeak_cap_exempt` | stress | Splits allowed above the cap |
| `amount.pcr_slope_per_unit` | calibrated 0.056, strong 0.112, none 0 | Minor-allele PCR share model |
| `amount.min_minor_share` | 0.05 | Floor on the minor-allele share when sizing amplicon templates |
| `amount.genomic_mc_draws` | 20000 | Monte-Carlo draws for genomic read counts |
| `amount.fragment_length_median` / `sigma` | 5000 / 0.5 | Fallback genomic fragment model (MucOneUp defaults) |
| `profiles.r10_pcr_alpha`, `strong_pcr_alpha_factor` | 9.27e-5, 2 | Strong PCR bias preset |
| `profiles.poor_error_scale`, `hifi_poor_accuracy_mean` | 1.5, 0.95 | `poor` error level |
| `structures.rare_fraction`, `rare_usage_max`, `stationary_steps` | 0.10, 0.01, 2000 | Rare-unit structures |
| `realism.*` | see the realism section | Metric definitions and tolerances |
| `report.alpha`, `ni_margin` | 0.05, 0.005 | Decision rule and interval level |
| `report.bootstrap_replicates`, `bootstrap_seed` | 2000, 0 | Cluster bootstrap |
| `run.threads` | 4 | Default `run --threads` |

Domain constants are not settings. The repeat-unit length and the conserved
unit IDs come from the bundled repeat dictionary. The conserved head (units
1-5) and tail (units 6-9) positions come from the bundled reference layout.

## Data layout outside Git

All datasets, reads, truth and results live outside the working tree. The
default `--out-root` is `<repository parent>/MucOneSpan-bench-data`; a path
inside the repository is refused.

```text
MucOneSpan-bench-data/
  designs_<split>.jsonl          # one design per line
  profiles/                      # generated profile variants
  <split>/
    manifest.jsonl               # one row per design, including failures
    realism.json, realism.md
    <design_id>/
      case.json                  # design, status, versions, realized depth, geometry
      truth/                     # MucOneUp haplotypes, structures, statistics
      reads/                     # FASTQ + read_truth.tsv.gz
  test/preregistration.jsonl     # append-only decision-rule ledger
  test/first_evaluation.json     # written once, when test truth is first read
  results/<split>/<engine>/      # engine output, inventory, evaluation.json
  results/<split>/report.json, report.md
```

## Splits, seeds and sealing

| Split | Cases per profile | Seeds |
| --- | --- | --- |
| `dev` | 300 | public salt; regenerable, used for development |
| `val` | 300 | public salt |
| `test` | 800 | **secret** salt file stored outside the working tree (`--salt-file`) |
| `stress` | 100 | public salt; hard corners (smear level above the real maximum) |

Every split and profile holds at least `design.normal_fraction` (35%) normal
cases. Regular splits use smear levels whose expected off-peak share (smear
plus chimera) stays within the real maximum; the `stress` split holds the
smear 0.5 level. Biological and read
seeds are derived from the salt and the design, so a split is reproducible
from its salt.

The `test` split is **sealed**: `evaluate`, `report` and `realism` refuse to
read `test` truth until the exact decision-rule text is pre-registered with
`preregister`. The first `evaluate` or `realism` of `test` writes
`test/first_evaluation.json`. The marker is written *before* scoring starts, so
a run that later fails still counts as unsealing; this is deliberately
conservative. After it exists, a new pre-registration is refused, so a rule
cannot be registered after the test truth has been seen. The audit (rule
SHA-256, registration and first-evaluation times) is copied into each
`evaluation.json` and `report.json`.

### Decision rule

Adopt a candidate engine over the baseline only if, for **every** profile, it
is superior on per-allele exact sequence (exact two-sided McNemar on paired
truth alleles, Holm-adjusted across the three primary metrics at alpha 0.05),
non-inferior on the false-positive `PATHOGENIC` rate among normal and benign
truths (Newcombe one-sided 95% upper bound of the difference below 0.5
percentage points), and calls no more pathogenic truths
`NO_PATHOGENIC_VARIANT_DETECTED` or `NO_CALL` than the baseline. The numbers
are `report.alpha` and `report.ni_margin`. The full text is `rule_text()` in
`muc_one_span.benchsim.report`, built from the report settings; its SHA-256 is
what `preregister` records, so changing a report setting needs a new
pre-registration.

## Subcommands

Examples assume `uv run --locked --all-extras` in front of
`python scripts/benchsim.py`, and `DATA=../MucOneSpan-bench-data`.

### `design`

Writes `designs_<split>.jsonl` with stratified factors (profile, event,
length difference class, depth, composition, PCR, smear, chimera, error).
Event targets never fall on the conserved head (units 1-5; unit 1 holds part
of the forward amplicon primer site) or the conserved tail (units 6-9, the last
four units). A drawn target that would land there (`first10` or `last10` of a
short allele) is clamped to the nearest allowed unit, and the design records
`target_clamped: true`. Clamping keeps the random draw sequence, so other
designs do not change. An allele too short to hold an event outside the head
and tail is `design_invalid`.

```bash
python scripts/benchsim.py design --split dev --n 30          # 30 per profile
python scripts/benchsim.py design --split test --salt-file ~/secrets/benchsim.salt
```

### `generate`

Simulates truth and reads for each design with MucOneUp and writes
`case.json`, `truth/`, `reads/` and the split `manifest.jsonl`. Cases are
written once; failed cases stay in the manifest with their status.

Amplicon template counts are sized so that the PCR-disadvantaged minor allele
reaches the design depth. With strong PCR bias and a large length difference
the minor-allele share can fall so low that the template count would explode.
The share is therefore floored at `amount.min_minor_share` (0.05). A floored
case records `amount_capped: true`, and every amplicon case records the floor
used (`min_minor_share`) in `case.json` and the manifest. The minor allele of
a floored case gets below-target depth, like real allelic dropout. The case
also records `target_clamped` from its design.

```bash
python scripts/benchsim.py generate --designs "$DATA/designs_dev.jsonl" --jobs 6 \
  --muconeup-config /path/to/MucOneUp/config.json \
  --muconeup-profiles /path/to/MucOneUp/muc_one_up/data/read_profiles
```

`--structure-pool` supplies local real-derived structures. Without it,
`real_derived` designs fall back to the Markov model and record
`composition_effective: "markov"`. `--flank-fasta` adds flanks for genomic
reads.

### `run`

Runs the engines over a split manifest. Cases that did not generate are kept as
`not_attempted`, so denominators never shrink. `run` is not resumable: it runs
every case again, so remove `results/<split>/<engine>/` before a clean rerun.
`--jobs` times `--threads` is the approximate core use.

```bash
python scripts/benchsim.py run --manifest "$DATA/dev/manifest.jsonl" --engines ladder \
  --model-ont /path/to/clair3/r1041_e82_400bps_sup_v500 \
  --model-hifi /path/to/clair3/models/hifi --threads 3 --jobs 4
```

### `preregister`

Appends the decision rule to `test/preregistration.jsonl` (only needed for
`test`).

```bash
python scripts/benchsim.py preregister
```

### `evaluate`

Scores each engine's results with `scripts/evaluate.py` into
`results/<split>/<engine>/evaluation.json`.

```bash
python scripts/benchsim.py evaluate --split dev --engines ladder
```

### `report`

Writes stratified tables and, when `--candidate` is given, the decision rule
to `results/<split>/report.json` and `report.md`: per-allele exact (metric 1),
case exact (metric 2), event recall and precision (metric 3), clinical
confusion per profile (metric 4), false positives and no-calls on normal
truths, and a failure atlas per design factor. Intervals are 95% cluster
bootstrap intervals over designs.

```bash
python scripts/benchsim.py report --split dev --baseline ladder
python scripts/benchsim.py report --split val --baseline ladder --candidate hybrid
```

### `realism`

Compares the simulated reads of a split with the public targets and writes
`<split>/realism.json` and `realism.md`.

```bash
python scripts/benchsim.py realism --split dev
```

## Realism report and the sim-to-real gap

Realism metrics are measured on the same scope as the targets: the VNTR
interval (motif-1 start to motif-9 end) of spanning reads for error rates,
C7 homopolymer accuracy and span offsets. Each metric is checked against the
public PRJEB92208 aggregates (`src/muc_one_span/benchsim/targets/prjeb92208_v1.json`)
with a tolerance from the `realism` settings:

| Setting | Default | Check |
| --- | --- | --- |
| `realism.error_rel_tol` | 0.20 | Error rates within 20% of the real median |
| `realism.c7_abs_tol` | 0.03 | C7 correct-length fraction, absolute |
| `realism.range_median_rel_tol` | 0.25 | Median of a range metric within 25% of the real median |
| `realism.jsd_max` | 0.1 | Span-offset histogram Jensen-Shannon distance |
| `realism.slope_abs_tol` | 0.01 | Allele-ratio slope per unit, absolute |

The metric definitions are settings as well:
- span-offset bins: `offset_bin_bp` 15, bins -12 to 4, size split at 55 units;
- C7: `c7_run_length` 7, `c7_extend_max` 3;
- peak widths: `on_peak_min_bp` 30, `on_peak_rel` 0.012, `off_peak_units` 1.5.

The bins must match the target file's `bin_lo_bp`, or `realism` fails with an
error.

A failed check is a known sim-to-real gap, not a reason to tune the caller. Report
it with the benchmark results. In-house genomic targets can be loaded
locally by path; they are never committed. `hifi_amplicon` has no target
section and is reported without a verdict.

## What may be committed

- Code, tests, design definitions and documentation.
- The public PRJEB92208 aggregate targets.
- Summary numbers from `dev` runs in `.planning/` (aggregates only).

Do not commit reads, truth, results, the `test` salt, patient reads,
in-house structures or aggregates, or any per-library output. Local inputs are
passed by path at run time and only their SHA-256 is recorded.
