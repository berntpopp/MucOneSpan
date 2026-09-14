# Evidence inventory and focused baseline — 2026-09-14

Scope: read-only inventory of the original checkout and prior ignored audit
artifacts, followed by five authorized fresh baseline executions against unchanged
production source. No simulator data were generated, no thresholds selected, and
no source files changed by this task.

## Verified provenance and retained inputs

All **155/155** entries in
`tests/results/deep_validation_20260914/provenance.json` match streaming SHA-256
recomputations (1 MiB blocks); no listed files are missing. All **77/77** input
hashes in the HiFi (44), ONT (3), heldout (6), and perturbation (24) measurements
also match. These are input/run counts, not 77 independent simulations.

The provenance manifest covers 141 of 309 files under `tests/data/generated/`:
47 simulation-statistics JSONs, 47 diploid truth FASTAs, and 47 structure text
files. Its remaining 14 hashes cover the repeat dictionary, ladder reference,
and six files in each Clair3 model directory. Input reads are fingerprinted in
measurements rather than that manifest.

It is **not a complete provenance inventory**. It omits 30 mutated-unit FASTAs,
47 read metadata TSVs, BAM indexes, simulator source/configuration, pbsim error
models, retained analysis code, and all 41 files under the previous heldout-data
root. Current supplementary fingerprints were recorded for the mutant-unit
FASTAs, heldout files, generation configs, and scripts. These establish their
current bytes, not historical immutability before this inventory.

All 44 original HiFi aligned BAMs, three original ONT FASTQs and six prior heldout
HiFi FASTQs are present. Original unaligned HiFi FASTQs, raw subreads, per-haplotype
CCS BAMs, pbsim per-read alignment truth, and molecule-to-source-haplotype mappings
are absent from the retained input directories. Original HiFi reads can be
extracted from their BAMs, but doing so cannot restore missing source assignment
truth. No real clinical sample cohort is available here.

All **53/53** full-data samples have duplicated read names: 2,950 duplicate records
beyond unique names across the set. BAM inspection counted primary records
(`samtools view -F 2304`); FASTQs were counted by records. Original BAM optional
tags are mapping tags only (`AS`, `NM`, `SA`, `cm`, `de`, `ms`, `nn`, `rl`, `s1`,
`s2`, `tp`), with no `HP` or explicit simulator source assignment. Renaming
records cannot recover the missing assignment. The metadata TSVs describe each
run and are not per-read truth tables.

## Available environment and reproducibility constraints

Current analysis environment: Python 3.12.9, Linux x86_64. Original pipeline
source commit is `d8390b3c244ef8f3240af74b92db12b50dfc77d1`; source/configuration
Git diff was empty for the fresh baseline. The unrelated documentation changes
in the original checkout were preserved.

| Purpose | Available installation | Verified version |
| --- | --- | --- |
| Analysis minimap2 | `~/miniforge3/envs/env_clair3/bin/minimap2` | 2.28-r1209 |
| Analysis samtools | same environment | 1.15.1 (HTSlib 1.17) |
| Analysis bcftools | same environment | 1.17 |
| Clair3 | same environment, `run_clair3.sh` | 1.0.10 |
| MucOneUp | `~/development/MucOneUp/.venv/bin/muconeup` | 0.44.5 |
| Simulator source | `~/development/MucOneUp` | clean commit `58f55a6e6040f4cad55d9dac5ac301b1ba5424b3` |
| pbsim3 | `~/miniforge3/envs/env_pacbio/bin/pbsim` | Conda metadata 3.0.5, build h9948957_2 |
| CCS | same PacBio environment | 6.4.0 |
| Generation minimap2 / samtools | same PacBio environment | 2.30-r1287 / 1.23 |

`pbsim --version` and `pbsim --help` reject those options; its Conda package
metadata supplies the version, and invoking it without arguments prints usage.
The help/version probe failure was retained, not represented as a successful
version response. A `pysam` inventory probe failed because the optional module is
not installed; the completed inventory uses existing samtools without adding a
dependency.

Both Clair3 model directories exist under
`~/miniforge3/envs/env_clair3/bin/models/{hifi,ont}` with matching recorded hashes.
Five pbsim model files are present under
`~/development/MucOneUp/reference/pbsim3/`; their resolved locations, sizes and
hashes are retained in `inventory_environment.json`. ERRHMM-SEQUEL and QSHMM-ONT-HQ
are accessible. Historical input metadata records MucOneUp 0.44.2; a clean current
0.44.5 checkout does not recreate that historical source/config state by itself.

Default shell PATH does not expose minimap2, bcftools, Clair3 or pbsim. Supply the
analysis environment prefix to subprocess PATH **after** uv environment setup,
so Clair3 retains its own Python interpreter. The worktree has its own `.venv`;
using an original-checkout Python installation without explicit import paths can
otherwise benchmark the wrong source.

## Focused baseline results

All five fresh pipelines completed. Evaluated per-sample rows exactly equal the
previous audit rows, including count errors, sequence distances, mutation calls,
extra calls, and haplotype assignment. This reproduction validates that the
retained problems still occur with the recorded original implementation.

| Sample | Truth / reported counts | Exact sequences | Wall seconds |
| --- | --- | --- | --- |
| HiFi bench5003 | 56/118 → 56/118 | 1/2 | 24.31 |
| HiFi dupC100/120 | 100/120 → 100/120 | 0/2 | 29.14 |
| HiFi same-length60/60 | 60/60 → 60/63 | 0/2 | 13.22 |
| HiFi normal60/80 | 60/80 → 60/80 | 1/2 | 8.77 |
| ONT normal60/80 | 60/80 → 60/81 | 1/2 | approximately 8.7 |

The four HiFi samples give 2/8 exact sequences and 7/8 exact counts, with 1 exact
mutation TP, 2 FN, and 1 extra VCF-supported mutation. ONT normal has no mutation
calls. These deliberately selected failures are development controls, not an
independent performance estimate. `bench5003` is 56/118 in stored truth despite
historical naming assumptions.

Fresh output root:
`tests/results/production_validation_20260914/baseline/`.
`commands.json` records exact subprocess argument lists, executable, PATH prefix,
source import path and run duration. Per-sample `measurement.json` records exact
CLI arguments, input hash, exit status and stage times; `cli.log` records output.
`{hifi,ont}_accuracy.json` contains full exact evaluations.
`inventory_baseline.json` fingerprints every baseline artifact and records the
complete equality comparison with previous evaluator rows.

## Commands for another isolated baseline

Use a fresh nonexistent result directory each time. Local evidence commands may
contain local installation paths; these are not production configuration defaults.
The existing ignored runner refuses existing sample output folders. It defaults
to four threads, uses the unmodified CLI, wraps stage functions only for timing,
and selects BAM before FASTQ. Explicit sample/platform selection is mandatory
because the input directory mixes platforms. It records a sample failure but
its own process may still exit zero, so check every `measurement.json` exit code.

From the original root, this form reproduces the baseline mechanics (substitute
an unused suffix for `baseline_repeat`):

```bash
uv run --locked --all-extras python - <<'PYRUN'
import os
import subprocess
import sys
from pathlib import Path
root = Path.cwd()
env = os.environ.copy()
env["PATH"] = str(Path.home() / "miniforge3/envs/env_clair3/bin") + os.pathsep + env["PATH"]
env["PYTHONPATH"] = str(root / "src")
command = [sys.executable, str(root / "tests/results/deep_validation_20260914/run_matrix.py"),
           "--data-dir", str(root / "tests/data/generated"),
           "--output", str(root / "tests/results/production_validation_20260914/baseline_repeat/hifi"),
           "--model", str(Path.home() / "miniforge3/envs/env_clair3/bin/models/hifi"),
           "--platform", "hifi"]
for name in ["sample_normal_60_80", "sample_dupc_100_120", "sample_bench_5003", "sample_homozygous_60_60"]:
    command.extend(["--sample", name])
subprocess.run(command, env=env, check=True)
PYRUN
```

For ONT use a distinct output subdirectory, model `models/ont`, platform `ont`,
and sample `sample_ont_normal_60_80`. For the production worktree change the
interpreter/current root/PYTHONPATH to that worktree and keep old input and runner
paths absolute. Verify the imported `muc_one_span.__file__` before execution.

Scoring command used for each platform:

```bash
uv run --locked --all-extras python results/classification-audit/evaluate_pipeline.py \
  tests/results/production_validation_20260914/baseline/hifi \
  --truth-root tests/data/generated \
  --output tests/results/production_validation_20260914/baseline/hifi_accuracy.json
```

## Fresh validation after the coordinator freezes settings

Do not rerun `mapping_audit/generate_heldout.py` directly: it hardcodes previous
seeds 9101–9106, previous output directories and simulator paths, and permits
existing directories rather than refusing them. Its six prior cases are already
observed development data. Likewise current ONT seeds 1001/1002/1006 are reused
biological designs, not new independent seeds.

Before generation, the coordinator should persist a frozen case manifest,
implementation/source hash, caller/consensus policy, scoring rules, seeds, model
fingerprints, and explicit no-call criteria. Use fresh never-generated seed IDs,
new output folders created with `exist_ok=False`, and an immutable copied config.
The accessible previous HiFi config is
`tests/results/deep_validation_20260914/heldout_data/generation_config.json`;
ONT generation config is `tests/results/deep_validation_20260914/ont_generation_config.json`.
Both specify explicit tools; human-reference alignment is omitted so FASTQ remains.

The verified simulator argument forms are:

```text
MUCONEUP --config CONFIG simulate --out-base SAMPLE --out-dir NEW_FOLDER
  --num-haplotypes 2 --fixed-lengths LENGTH1 --fixed-lengths LENGTH2
  --output-structure --seed FRESH_SEED
  [--mutation-name MUTATION --mutation-targets HAPLOTYPE,REPEAT]
MUCONEUP --config CONFIG reads amplicon NEW_FOLDER/SAMPLE.001.simulated.fa
  --out-dir NEW_FOLDER --out-base SAMPLE_reads --coverage TEMPLATE_COUNT
  --seed FRESH_SEED --platform pacbio
```

For ONT select `--platform ont` and freeze its error-model/accuracy parameters
explicitly. Confirm exact flags using the captured full `reads amplicon --help`.
Requested coverage is **initial template molecules**, not achieved depth; metadata
coverage currently remains 30 even when command requests 200. Derive achieved
read support from records and retained source assignment, not that metadata cell.

A complete new validation needs per-haplotype read truth captured before the
simulator merges haplotype reads. Source inspection identifies the PacBio boundary
at `amplicon_pipeline.py` stage 7 (`merge_bam_files(input_bams=hifi_bams)`), after
per-haplotype CCS generation; temporary `hifi_hap{i}_{j:04d}.bam` paths carry the
source identity. ONT similarly produces separate `ont_hap{i}` files before merging.
Capturing those outputs in an external wrapper can retain current diploid PCR
allocation without altering simulator source; record wrapper code/hash and prove
that capture preserves the sequence/quality multiset. Generate globally unique
record IDs with a sidecar linking them to original read ID and source haplotype.
Do not infer source assignment from ambiguous alignment to similar truth alleles.

This inventory has not generated these data or implemented that capture wrapper.
It establishes available inputs and commands; the frozen plan determines the
actual validation matrix and generation/capture implementation.

## Inventory artifacts and executed validation

Ignored root `tests/results/production_validation_20260914/` contains:

- `inventory_hashes.json`: 155 provenance and 77 input streaming hash checks;
  explicit coverage lists for omitted files.
- `inventory_environment.json`: exact executable/help/version responses and model
  file fingerprints, including failed pbsim version probe.
- `inventory_reads.json`: sample record counts, duplicate names, BAM tag inventory,
  supplementary current hashes and pbsim package metadata.
- `inventory_baseline.json`: unchanged source commit/diff, prior/current comparison,
  and fresh baseline artifact fingerprints.

No unit/integration suite was claimed for this read-only inventory; the five real
pipeline runs and two exact evaluator invocations are the executed behavior checks.
