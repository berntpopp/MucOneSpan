# Benchmarking with MucOneUp

[MucOneUp](https://github.com/berntpopp/muconeup) generates simulated MUC1 VNTR
haplotypes and amplicon reads with known ground truth. Use it to evaluate allele
length detection and mutation classification separately from unit tests.

## Prerequisites

Run `make dev` for the locked MucOneSpan environment. The simulation workflow
also requires MucOneUp, its `config.json`, and pbsim3 in a suitable external tool
environment. Install minimap2 and samtools for alignment; full pipeline runs also
need bcftools, Clair3, and the appropriate model. The repository's
`conda/environment.yml` supplies the core alignment/consensus tools, but does not
install the simulator or Clair3.

Activate the tool environment before running commands. Point to your MucOneUp
configuration explicitly rather than relying on a developer's directory layout:

```bash
export MUCONEUP_CONFIG=/path/to/muconeup/config.json
make generate-testdata
```

The documented historical runs used MucOneUp 0.44.2. Record your exact simulator,
pbsim3, alignment, caller, and model versions when producing new results. The
current generator needs a MucOneUp version supporting the ONT amplicon option.

## Generated samples

`examples/development-experiment.json` defines sample names, seeds, mutation
targets and requested templates. The configurable runner
`scripts/generate_testdata.py` preserves its **26 HiFi samples and 3 ONT samples**
under `tests/data/generated/`.

| Group | Cases |
| --- | --- |
| Core mutations | dupC, dupA, insG, insCCCC, del18_31, normal |
| Length edge cases | Same-length 60/60, asymmetric 25/140, short 25/30, long 120/140 |
| Reproducibility | Additional dupC seeds and 40/50, 80/100, 100/120 lengths |
| Close alleles | 50/55, 50/57, 50/60 with mutation and normal controls |
| Long allele mutations | dupA and insG at 100/120, plus a normal control |
| Low coverage | dupC with 50 requested templates |
| ONT | dupC, dupA, and normal at 60/80 |

The default requests 200 template molecules; the low-coverage sample requests50.
These are pre-filter simulation requests, not retained reads or measured depth.
The manifest records actual usable reads and failed generation. See the
[simulation experiment guide](simulation-experiments.md) for custom JSON designs
and explicit per-platform simulator configurations. Ground truth
includes simulated haplotype FASTA and simulation statistics. Keep generated
reads and results out of Git; record parameters and summaries needed to reproduce
an experiment.

## Integration tests

```bash
make test-int
```

`make test-int` selects the `integration` marker. To request generated-data
end-to-end tests as well, run:

```bash
uv run --locked --all-extras pytest tests/integration -m "integration or e2e" --no-cov
```

Tests declare required tools and fixtures. Inspect skip reasons: missing tools,
models, or generated data can leave scientific behavior untested. pbsim3 is only
needed to generate reads; tests using existing reads do not need to rerun it.
The unit suite is independent and runs with `make ci-check`.

## Batch analysis

Set the model location appropriate for the platform and invoke the batch script
from the repository root:

```bash
export CLAIR3_MODEL=/path/to/clair3/models/hifi
uv run --locked --all-extras python scripts/batch_analyze.py \
  /path/to/hifi_samples tests/results/hifi --platform hifi
```

The generator writes HiFi and ONT samples together. Before a platform-specific
batch run, prepare an input directory containing only that platform's samples;
do not analyze ONT reads with the HiFi preset. For an ONT-only directory:

```bash
export CLAIR3_MODEL=/path/to/clair3/models/ont
uv run --locked --all-extras python scripts/batch_analyze.py \
  /path/to/ont_samples tests/results/ont --platform ont
```

The batch output includes `batch_results.json` and per-sample pipeline results.
Inspect error counts as well as TP/FN/FP/TN categories. A completed batch is not
necessarily an accurate or error-free batch. Avoid reusing stale outputs when
comparing implementations; use separate result directories.

For a single sample:

```bash
uv run --locked --all-extras muconespan run \
  --input /path/to/sample_reads.bam \
  --output-dir tests/results/single \
  --clair3-model "$CLAIR3_MODEL" \
  --platform hifi
```

## Interpretation and historical results

Assess repeat classification against known sequences, allele lengths against
ground truth, and full pipeline mutation calls separately. Record false positives,
false negatives, partial matches, errors, model/tool versions, seeds, and coverage.
Pay particular attention to same-length alleles, nearby lengths, long tandem
repeats, low coverage, and PCR bias.

The [historical benchmark report](https://github.com/berntpopp/MucOneSpan/blob/main/.planning/BENCHMARK_RESULTS.md)
records a 44-sample v0.3.0 run: 25/28 mutation detections and 15/16 correct normal
calls. That dataset includes additional experiments beyond the generator's core
catalog. Its observations include a same-length allele splitting failure and
missed single-base insertions in long repeats. These are version-specific
observations, not current validation results or clinical performance claims.

See [limitations](../reference/limitations.md),
[known mutations](../reference/mutations.md), and
[CLI reference](../reference/cli.md) for context. The longer
[historical testing notes](https://github.com/berntpopp/MucOneSpan/blob/main/.planning/TESTING_WITH_MUCONEUP.md)
retain past experiments; use the portable commands here for new work.

## Strict offline accuracy evaluation

The historical batch substring categories (`TP_partial`, or a template name found
anywhere in a sample) are not exact mutation accuracy. Use the versioned evaluator:

```bash
uv run --locked --all-extras python scripts/evaluate.py tests/results/run \
  --truth-root tests/data/generated --expected-samples samples.json \
  --output tests/results/run/evaluation.json
```

`samples.json` explicitly lists every expected sample. Objects can point to distinct
truth/output roots for perturbations and record exact input, platform, model and
execution metadata. Missing/invalid truth, execution failures and unattempted
samples remain visible. Nonzero status indicates evaluation/input/execution
failure; completed evaluation with poor biological accuracy is a valid report.

Compare total `length` to simulator total repeats; `contig_N` and
`canonical_repeats` exclude nine fixed units. Use actual truth FASTA/structures,
not filenames, stale per-repeat lengths or requested coverage metadata. Requested
simulator coverage may be template molecules before filtering. Measure retained
records and support instead of calling that value achieved sequencing depth.

Report exact full sequences, ordered structures and count pairs separately from
exact event annotation and broad mutation alarms. Keep positives' extra calls,
missing/extra alleles, ambiguity, no-calls and normal false positives. Report
platforms separately and preserve sample denominators. Related seeds and depth
subsets do not supply independent clinical performance estimates.

For classification timing, `scripts/benchmark_classification.py --help` documents
paired scalar/bit-vector runs with complete result equality. A local scorer gain
is not a whole-pipeline speedup. Pipeline comparisons require identical inputs,
model/reference hashes, thread budgets, stage/wall measurements and paired runs.

The current `scripts/batch_analyze.py` writes a versioned evaluation object in
`batch_results.json`, replacing the legacy array containing partial true-positive
labels. `scripts/validate_catalog.py` accepts both formats for display; accuracy
claims must use the strict schema. Exact annotation requires parent, name and
1-based total-repeat position on a one-to-one matched allele. All positive extras
remain false positives, including on mutation-positive samples.

Specificity uses confident TN/(TN+FP), retaining false positive alarms even when
sequence reconstruction is ambiguous. The confident-negative rate additionally
uses all normal controls as its denominator. Unresolved negatives are reported
separately and cannot improve either result by becoming no-calls. Candidate-output
call rate and resolved-reconstruction call rate answer different questions.

### Driver exit status

Batch and benchmark drivers return nonzero when any caller exits nonzero, including
a typed `insufficient_evidence` no-call. The offline evaluator can still complete
successfully and score that no-call as a scientific outcome. Use the evaluation
report's execution and no-call fields alongside the driver exit status.
