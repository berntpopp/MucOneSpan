# Reproducible simulation experiments

Use `scripts/generate_testdata.py` to generate a declared HiFi or ONT experiment
with MucOneUp. The runner writes a complete inventory, validates generated truth,
and counts usable records. It does not run MucOneSpan or tune its scientific
settings. The commands below run from the MucOneSpan repository root.

## Requirements

Install MucOneUp with support for `reads amplicon --platform ont`; the maintained
workflow was exercised with version 0.44.5. PacBio simulation requires pbsim3 and
CCS; ONT uses single-pass pbsim3. Both require samtools. Optional simulator
alignment also requires minimap2 and a reference. Activate the appropriate tool
environment before invoking `uv run --locked --all-extras`.

Supply your own MucOneUp configuration with amplicon primers and explicit
`pacbio_params.model_file` and/or `ont_amplicon_params.model_file`. These are
sequencing simulation models, separate from the Clair3 models used for calling.
Use an absolute `--muconeup` executable path when it is installed in a virtualenv:
MucOneSpan's tool helper removes virtualenv directories from the external PATH.

```bash
uv run --locked --all-extras python scripts/generate_testdata.py \
  --design examples/ont-experiment.json \
  --config /path/to/muconeup/config.json \
  --muconeup /path/to/muconeup-environment/bin/muconeup \
  --output-dir tests/data/generated/ont-development --dry-run
```

Remove `--dry-run` to generate the declared cases. Dry runs print the planned
argument lists, create no files, and execute no tools; they do not verify that
models, executable or simulator configuration are usable. Real runs refuse an
existing experiment manifest, inventory, or any declared case directory. Choose
a new output directory to repeat an experiment, including after a failure.

## JSON experiment design

The default design is `examples/development-experiment.json`. It preserves the
historical 26 HiFi and three ONT names, seeds, mutation targets and requested
template counts. Edit a copy for new experiments. The small ONT example is a
**development smoke case**, not an accuracy study or an unexposed validation set.

`examples/future-validation-experiment.json` declares a larger 32-case panel
(16 HiFi and 16 ONT) with seeds 2026091401–2026091432. Those cases have not
been generated for the release evidence. Treat it as a future experiment design,
not a completed held-out validation result; preserve the protocol and document
any changes before generating or examining its outputs.

```json
{
  "schema_version": 1,
  "platform_configs": {
    "hifi": "configs/hifi-simulator.json",
    "ont": "configs/ont-simulator.json"
  },
  "cases": [
    {
      "sample": "ont_development_25_30",
      "platform": "ont",
      "seed": 19200902,
      "lengths": [25, 30],
      "mutation": "dupC",
      "targets": [[1, 10]],
      "requested_templates": 12
    }
  ]
}
```

Every case requires all seven fields. `sample` is unique and consists of letters,
digits, underscores, periods or hyphens, starting with a letter or digit.
`platform` is `hifi` or `ont`; the runner translates `hifi` to MucOneUp's `pacbio`.
`seed` is a nonnegative integer. `lengths` contains two positive total repeat
counts, including fixed units. `requested_templates` is a positive integer.

`mutation` names one mutation supported by the simulator and the truth
dictionary. `targets` contains distinct `[haplotype, repeat]` pairs: both indices
are 1-based, haplotype is 1 or 2, and the repeat must lie within its requested
length. Multiple positions may share the same mutation. For a normal case use
`"mutation": null` and `"targets": []`. Different mutation names within one case,
SNPs and arbitrary input structures are outside this runner's schema.

Optional `repeat_dictionary` selects a MucOneSpan-format dictionary used by the
strict truth adapter. It must describe the same repeat sequences, mutations and
flanks as the simulator. The bundled dictionary is used when absent. Generated
FASTA, ordered structures, mutation units and statistics must agree before a
case is marked completed; the requested lengths and mutation targets must also
match the validated truth.

Unknown fields, duplicate JSON keys, invalid types, nonfinite values and
out-of-range targets are errors. Booleans are not integers. Expand a seed or
coverage sweep into explicit case rows with unique names; there is no implicit
Cartesian expansion or automatic selection of favorable outputs.

Configuration precedence is:

1. Explicit CLI `--config`, applying to every platform.
2. The design's `platform_configs` entry for that case's platform.
3. The `MUCONEUP_CONFIG` environment variable.

Design paths resolve relative to the design JSON; CLI and environment paths
resolve relative to the invocation directory. Simulator commands execute from
the selected simulator configuration's directory, so its relative model and
reference paths have a stable base. There are no home-directory or sibling-repo
fallbacks. The runner does not infer a sequencing model from a sample name.

## Artifacts and limitations

`generation_manifest.json` records the design and its hash, exact command
argument lists and working directories, simulator version and executable hash,
configuration snapshots, explicit model and existing file-valued configuration
hashes, command results/errors, per-case generation wall time, output hashes and
actual usable record counts. Referenced artifacts are checked before and after each
simulation command; changes stop that case. The simulator executable is resolved
using the tool helper's cleaned PATH and pinned for version probing and generation. The manifest is updated as commands
run. Interrupted runs retain their last recorded state rather than claiming
completion. The hashes are not a lock of the simulator's complete Python and
external-tool environment; archive those environments separately for exact
reproduction.

`inventory.json` retains every declared case, its platform, truth directory and
explicit input. Failed cases retain a nonexistent input sentinel, so downstream
benchmarking records their failure instead of discovering leftover intermediates.
Failures return a nonzero CLI exit status. Correct remaining cases still run.
Outputs are never selected by truth agreement among alternatives: multiple
candidate read files are an error.

For amplicon simulation, MucOneUp's `--coverage` is the number of template
molecules before PCR allocation and CCS filtering. It is not achieved depth or
final read count. The manifest separately counts valid FASTQ records or non-secondary, non-supplementary
BAM records (including unmapped records with sequence and quality) with sequence and quality, retaining independent records with equal
QNAMEs. This is input usability, not evidence of full VNTR spanning or allele
coverage. Simulator output may be FASTQ when alignment is disabled.

MucOneUp 0.44.5 explicitly rejects `--track-read-source` for amplicon simulation,
even though help lists the generic option. The runner records source assignment
as unavailable and does not derive source haplotypes from QNAMEs. Template-level
source recovery needs a separately validated simulator extension. Simulator
statistics may omit their seed or retain stale per-repeat lengths after mutation;
the strict adapter preserves those warnings, while the generation manifest
records the actual requested seed and independently validates sequence truth.

## Calling and evaluation

For an ONT-only design, use the generated inventory directly:

```bash
uv run --locked --all-extras python scripts/benchmark.py \
  --data-dir tests/data/generated/ont-development \
  --expected-samples tests/data/generated/ont-development/inventory.json \
  --output-dir tests/results/ont-development \
  --clair3-model /path/to/clair3/ont-model --threads 4

uv run --locked --all-extras python scripts/evaluate.py \
  tests/results/ont-development \
  --truth-root tests/data/generated/ont-development \
  --expected-samples tests/data/generated/ont-development/inventory.json \
  --output tests/results/ont-development/evaluation.json
```

For mixed platforms, make an analysis inventory copy with a platform-appropriate
`model` on every row, and omit the shared `--clair3-model`; also unset
`CLAIR3_MODEL`. Preserve all rows, including failures, and archive the copy.
Repeat calling with separate result directories for comparisons.

Report HiFi and ONT separately. Include exact sequence/count/event recovery,
missing and extra alleles, positive extra calls, no-calls, unresolved negatives
and execution failures. Parameter exploration uses development data; it does not
provide independent validation. See [benchmarking](benchmarking.md) and
[limitations](../reference/limitations.md) for interpretation.


Custom `repeat_dictionary` validates generation truth only. The generated
inventory does not configure a caller dictionary, and the existing benchmark and
evaluate CLIs use the bundled dictionary. For a custom dictionary, invoke
`muconespan --config ... run` with a matching explicit reference and use the
Python evaluation API with that same dictionary; custom-dictionary end-to-end
benchmark CLI interoperability is not provided by this release.

The runner requires exactly one candidate read input per case. If your simulator
configuration retains both FASTQ and an aligned BAM, discovery fails explicitly;
disable optional simulator alignment for this workflow. A dry run does not test
which files the simulator will emit.
