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
- `--clair3-model`, `--min-qual` and `--minimap2-preset` only apply to the
  ladder engine. A hybrid run accepts and ignores them.
- `--report-igv embedded|sidecar` is **rejected** by the hybrid engine (a
  hybrid run has no alignment tracks). A script that passes `--report-igv`
  must add `--engine ladder` or drop the option.
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
- `summary.json["deprecations"]` is new for both engines.

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

1. **Fail-closed depth rule.** A ladder run in which one allele's
   `depth_status` is `not_assessed` while another allele's depth is assessed is
   now INCONCLUSIVE (0.16.1: NEGATIVE).
2. **Reason text.** The unresolved-selection reason no longer prints
   `secondary mode fraction None` when the fraction is missing; the fraction is
   printed only when it is known.

### Dependencies

- `edlib`, `pyabpoa` and `pyspoa` are **core dependencies**; a plain
  `pip install muc_one_span` installs them. The `hybrid` extra is kept, so
  `pip install 'muc_one_span[hybrid]'` still works.
- `pyabpoa` is published as a source distribution only and needs a C compiler
  and zlib headers (`gcc`, `libc6-dev`, `zlib1g-dev` on Debian/Ubuntu).
  `pyspoa` has Linux wheels only; on macOS it builds with cmake and a C++
  compiler. See the [installation guide](../getting-started/installation.md).

### Benchmark harnesses

- `scripts/benchmark.py --engine`, `scripts/clinical_benchmark.py run --engine`
  and `scripts/benchsim.py run|evaluate --engines` default to `hybrid`.
- `muc_one_span.benchmarking.run_pipeline` always passes `--engine` to the CLI,
  so a ladder benchmark stays a ladder benchmark.
- `clinical_benchmark.py` hashes a ladder run without an `engine` key (the
  pre-0.17 hash) and passes `--engine ladder` to the worker explicitly.
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
