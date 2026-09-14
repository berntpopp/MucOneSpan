# Fresh paired analysis design — fixed before exposure

Scope: ignored offline `fresh_generation/analyze.py` and `analyze_check.py`, with
no production source or caller-setting changes and no new simulations/pipelines.
The analyzer will be tested on hand-built inputs before reading fresh results.

## Inputs and accounting

Use explicit per-arm primary inventories (128 expected entries), measurements
(136 expected entries), archived paired settings and generation manifest. The
four predetermined timing cases have three repeats; repeats2/3 are performance
observations only and never inflate scientific denominators. Validate unique IDs,
paired sample/metadata/input identities and exact expected measurement IDs. Missing
or failed outputs remain rows. Do not choose samples by output-directory discovery.

Use the frozen current evaluator for both arms, checking candidate source hashes
against archived paired settings. Preserve all sequence-optimal assignments and
conservative metric bounds. Write standard evaluator reports for baseline and
candidate so separate source-assignment tooling can consume them unchanged.

Enrich separate fresh read-source provenance from the archived generation run
records and hash-verified source sidecars. The legacy truth adapter's unavailable
label describes its own scope; it is not a conclusion about the new captured
source maps. Keep captured counts, usable records, requested templates, platform,
seed, generation failure and model/source hashes explicit.

## Scientific outputs

Main inputs and three dependent perturbation kinds are separate strata, separately
by HiFi/ONT. Any combined perturbation tables remain explicitly dependent, not
independent sample replicates. Report all assigned event annotations by exact
name/parent/repeat index/truth haplotype for every optimal sequence assignment,
including missing truth and extra predictions. Preserve supported/legacy/all-event
metrics as distinct endpoints from the official evaluator.

Report individual and independent sequence/structure recovery, count error >2,
missing and independently missing alleles, extras, normal TN/FP/unresolved,
call/no-call, per-sample gains/regressions and assignment ambiguity. Conditional
common-output and common-completed subsets are labeled selection-conditioned and
never replace complete-inventory results. Exact sequence-output matching is a
separate paired subset. Descriptive Wilson intervals apply to sample-level binary
endpoints only; no haplotype/event pooling as independent Bernoulli trials.

## Performance outputs

Keep all136 paired executions separate from the fixed four-case three-repeat
benchmark, and provide per-case repeat summaries. A speedup is baseline divided
by candidate, computed only for positive finite paired measurements from both
successfully completed runs; expose exclusions and execution outcomes. Summarize
paired ratios, paired differences and ratio of summed times separately. Do not
impute missing stages/resources as zero. Worker stage timing, subprocess wall,
self/children CPU and self/children peak RSS remain separate measures. Linux
waited-children maxRSS is not additive or concurrent process-tree peak memory.
Classification speedup cannot be described as whole-pipeline speedup.

## Pre-exposure checks

Hand-built cases must catch inventory mismatch/duplicates, missing measurements,
failed-call denominator retention, no favorable tie-breaking, exact event matching,
ratio orientation, zero/nonfinite metric exclusions, dependent-stratum selection,
repeat-only timing selection, normal unresolved denominators and Wilson boundaries.
End-to-end fixture analysis uses real strict truth/artifact adapters without
external bioinformatics tools. No actual final results are read in these checks.

## Implemented handoff and verification

The implementation is split by responsibility into ignored files under
`tests/results/production_validation_20260914/fresh_generation/`:
`analyze.py` (401 lines), `analyze_metrics.py` (274 lines), and
`analyze_check.py` (392 lines). No production modules were changed for this task.

Run only after the paired driver has finished, from the candidate worktree:

```bash
uv run --locked --all-extras python tests/results/production_validation_20260914/fresh_generation/analyze.py \
  --paired-root tests/results/production_validation_20260914/paired_final \
  --output tests/results/production_validation_20260914/fresh_analysis
```

The output directory must not already exist. The public CLI has no source-freeze
bypass and enforces 128 primary pairs and 136 execution pairs, including exactly
four predetermined cases with repeats2/3. A private function override exists only
for hand-built tests. Unknown measurement status, running workers, duplicated or
unplanned IDs, unbalanced final strata, inconsistent paired input hashes and
inconsistent repeat identity fail visibly. Not-attempted/generation-failed cases
remain inventory rows; their absent execution metadata is not imputed.

Outputs are standard official evaluator envelopes `evaluation_baseline.json` and
`evaluation_candidate.json`, plus `paired_analysis.json`. Original alternatives,
pairs, metrics, prediction events and run records are preserved. Per-sample
additions are inventory, all event-assignment interpretations, count errors>2,
matched-output count and separate fresh read-source provenance. The paired report
uses `primary_pairs=128`, `execution_pairs=136`, `timing_only_pairs=8`; it records
hashes for both analysis modules and all six archived input files. Separate
source-assignment tooling can consume the standard evaluation files unchanged.

Verification on 2026-09-14, before fresh result exposure:

- Ruff check and format check: passed for all three files.
- `uv run --locked --all-extras python .../analyze_check.py`: passed.
- CLI `--help`: passed.
- Hand-built end-to-end strict truth/artifact evaluation: passed, retaining all
  four primary rows and a failed candidate with stale output as missing/unresolved.
- Timing fixture: total-stage wall speedup4× versus classification8× remains
  separate; one failed pair excluded from ratios but included in wall burden.
- Exact 128/136 synthetic layout, duplicate/missing/unplanned entries, mismatched
  inputs/repeat indices, unknown status, normal unresolved denominators, Wilson
  boundaries, invalid truth retention, count error>2, sidecar tamper, modified
  frozen-source rejection, and both tied event assignments: passed.

The tests use temporary hand-built artifacts, no simulator, caller, external
bioinformatics tool, final seeds or actual fresh result data. They validate the
analysis contracts, not the future biological outcomes. Wilson intervals are
sample-level descriptive summaries and do not establish independent replicates.
Hash verification establishes archive integrity; detailed source-map semantics
and source-to-output compatibility belong to the separate assignment analysis.

Frozen handoff SHA256 values:

```text
8eb089edd83a3bed1752a8b88c39d3e8d719f73a83b5e4818f4244636d2acb5d  analyze.py
0f6f3d594919522420a5411dcf830fadb985f1098fb1594faf9e6b5634836f7b  analyze_metrics.py
ebdcbf7956b911d276f556d8460532419360f7d2fbe0eddfaf6b6ffec0328e10  analyze_check.py
```
