# Migrating to 0.17.0

Version 0.17.0 makes the read-centric **hybrid engine the default** for every
input type (amplicon and genomic) and **deprecates the ladder engine**. This
page lists what changes for existing users and how to keep the previous
behaviour while it is still available.

## What changes

### `muconespan run` uses the hybrid engine by default

- `muconespan run` without `--engine` now runs `--engine hybrid`
  (`run.engine = "hybrid"`). The engine reconstructs each allele from the
  reads (motif anchoring, a length model, partial-order-alignment consensus,
  phase splitting, read assignment, polishing and per-event read-level
  support); see [Core Concepts](../getting-started/concepts.md#hybrid-engine).
- A FASTQ run needs no external tool: minimap2, Clair3 and bcftools are not
  called. A BAM input still needs `samtools` to extract the primary reads.
- `--clair3-model`, `--min-qual`, `--minimap2-preset`, `--platform`,
  `--min-coverage`, `--threads`, `--mapping-timeout` and `--reference` only
  apply to the ladder engine: the hybrid engine takes its thresholds from
  `hybrid.*` (for example `hybrid.depth_adequate_spanning` and
  `hybrid.min_peak_reads` instead of `--min-coverage`), runs single-threaded,
  needs no platform and builds its own references. A custom dictionary or
  reference layout therefore no longer needs `--reference` with the hybrid
  engine. When a hybrid run gets one of these options on the command line, or
  through a non-default value in a `--config` file, it prints
  `Warning: <option> is ignored by the hybrid engine; use --engine ladder (deprecated)`
  on stderr and records the option in `summary.json["ignored_options"]` and
  `run_configuration.json["ignored_options"]` (additive). A hybrid run records
  `resolved_minimap2_preset: null` and `model_selection: "not used (hybrid engine)"`
  in `run_configuration.json`. The run still completes.
- `--report-igv embedded|sidecar` is **rejected** by the hybrid engine (a
  hybrid run has no alignment tracks) before any output is written. The error
  names both remedies: `--engine ladder` (deprecated) or `--report-igv off`.
- `--assay {amplicon,genomic}` is recorded for provenance only
  (`summary["hybrid"]["assay"]`); it does not change any setting and is not
  auto-detected.

### The ladder engine is deprecated

- `--engine ladder` (or `run.engine: "ladder"` in a `--config` file) still
  works. The run prints one warning line on stderr:
  `Warning: The ladder engine (--engine ladder / run.engine = ladder) is deprecated ...`.
- `summary.json` gains an additive field, `deprecations`: an empty list for a
  hybrid run, and for a ladder run one record
  `{"setting": "run.engine", "value": "ladder", "status": "deprecated", "replacement": "hybrid", "message": "..."}`.
- The HTML report of a ladder run (`--report`) shows a "Deprecated" banner
  with the same message.
- **Removal policy.** The ladder engine and its ladder-only options stay
  available and tested throughout the 0.17 series. Removal happens no earlier
  than the next minor release, and it is announced in the changelog of the
  release before it. Until then ladder results keep their output schema.

### Output

- Hybrid runs add fields and never remove ladder fields: per-allele
  `spanning_reads`, `assigned_reads`, `depth_status`, `depth_basis`,
  `selection_status`, `selection_detail`, `split_basis`, `phase_status`,
  `residual_sites` and `consensus_concordance_fraction`; per-sample
  `summary["hybrid"]` (read categories, rejected length peaks, dimer and smear
  results, POA backend); per-event `read_support`. Files: `hybrid_reads.json`,
  `hybrid_references.fa`, `consensus_*.fa`.
- A hybrid run writes no `mapping.bam` and no per-allele VCF. Mutation
  support comes from `read_support` (`supported`, `insufficient_depth`,
  `discordant`, `not_supported`, `not_localized`) instead of VCF concordance.
- `summary.json["deprecations"]` and `summary.json["ignored_options"]` are new
  for both engines (empty lists when nothing applies).

### Clinical decisions can differ from the ladder

The two engines use different evidence, so the same sample can get a
different decision:

- The hybrid engine blocks a NEGATIVE result whenever it cannot resolve the
  sample: a gate-relevant rejected length peak, an unresolved heterozygous
  site or homopolymer run, residual read heterogeneity, or a per-allele
  spanning depth below `hybrid.depth_adequate_spanning`. Such samples are
  INCONCLUSIVE with a located reason.
- A PATHOGENIC result needs the event's own read support to be `supported`.
- Measured differences and limits are on the
  [limitations page](../reference/limitations.md#hybrid-engine).

### Ladder-visible changes since 0.16.1

Two ladder behaviours changed together with the hybrid work:

1. **Fail-closed depth rule.** When an allele carries a `depth_basis`, any
   `depth_status` other than `adequate` (including an unknown or missing
   value) now blocks a result. In a ladder run where one allele's
   `depth_status` is `not_assessed` while another allele's depth is assessed:
   - a result that would have been NEGATIVE is INCONCLUSIVE (0.16.1: NEGATIVE);
   - a frameshift on the `not_assessed` allele is no longer PATHOGENIC; it is
     reported as uncertain and the result is INCONCLUSIVE.
2. **Reason text.** The unresolved-selection reason no longer prints
   `secondary mode fraction None` when the fraction is missing; the fraction is
   printed only when it is known.

### Dependencies

- `edlib` and `pyabpoa` (the default POA backend) are **core dependencies**; a
  plain `pip install muc_one_span` installs them. `pyabpoa` is published as a
  source distribution only and needs a C compiler and zlib headers (`gcc`,
  `libc6-dev`, `zlib1g-dev` on Debian/Ubuntu).
- `pyspoa`, the alternative backend (`hybrid.poa_backend: "pyspoa"`), stays in
  the optional `hybrid` extra (`pip install 'muc_one_span[hybrid]'`). Selecting
  it without the extra fails with an `ImportError` that names the extra; there
  is no silent fallback. It has Linux wheels only; on macOS it builds with
  cmake and a C++ compiler. See the
  [installation guide](../getting-started/installation.md).

### Benchmark harnesses

- `scripts/benchmark.py --engine`, `scripts/clinical_benchmark.py run --engine`
  and `scripts/benchsim.py run|evaluate --engines` default to `hybrid`.
- `muc_one_span.benchmarking.run_pipeline` always passes `--engine` to the CLI,
  so a ladder benchmark stays a ladder benchmark. It passes `--clair3-model`,
  `--threads` and `--platform` only to the ladder engine (still recording them).
- `benchsim run` looks up a Clair3 model only for the ladder engine; a hybrid
  run needs none.
- `clinical_benchmark.py` hashes a ladder run without an `engine` key (the
  pre-0.17 hash) and passes `--engine ladder` to the worker explicitly. The
  worker gets `--platform`, `--clair3-model` and `--threads` only for a ladder
  run. `clinical_benchmark.py freeze --model` is optional: without it the frozen
  environment has no Clair3 model (only `samtools` is frozen) and can run the
  hybrid engine only; a ladder run on it stops with an error.
- `benchsim report --baseline` stays `ladder`: the decision rule compares a
  candidate against the ladder baseline.

## Keep the ladder engine for now

```bash
muconespan run \
  --engine ladder \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi
```

or in a runtime settings file:

```json
{"schema_version": 1, "run": {"engine": "ladder"}}
```

The ladder engine still needs minimap2, samtools, bcftools and Clair3.
