# MucOneUp validation and performance investigation

## Scope and status

Completed investigation of MucOneSpan 0.10.0 at commit d8390b3. Initial working tree
was clean. This is a testing/debugging study; changes to scientific algorithms
will be proposed from measured failures rather than bundled into the baseline.

## Evidence required

- Inventory simulator truth, input hashes, versions, model, seeds, coverage and
  actual read counts; identify mismatches between labels and generated data.
- Execute current unit/static checks and actual external-tool integration tests.
- Rerun all 44 existing HiFi datasets into a new ignored output directory.
- Measure sample-level detection separately from exact mutation name, haplotype,
  repeat index, false extra calls, count error, full repeat structure and sequence.
- Classify truth haplotypes directly to isolate classification from upstream loss.
- Investigate same-length, close-length, asymmetric, low-coverage and long VNTRs.
- Exercise ONT simulations where prerequisites permit; report unavailable or
  failed cases explicitly, including simulator/model mismatch.
- Profile stages and representative cases, then rank improvement experiments by
  observed cost/error contribution and sensitivity/specificity tradeoffs.
- Report denominators and uncertainty; simulations are a convenience challenge
  set, not independent clinical validation. Partial/wrong calls are not exact TPs.

## Workstreams

1. Parent: environment, fresh end-to-end run, integration checks, statistics,
   coverage perturbations, final synthesis.
2. truth_audit: truth conventions and benchmark measurement defects.
3. classification_audit: truth-only reconstruction, mutation accuracy, profiling.
4. mapping_audit: allele assignment, caller/consensus assumptions and bottlenecks.

## Initial verification

- `make dev`: passed, locked environment synchronized.
- `make ci-check`: passed; 239 unit tests, 86.11% branch-aware aggregate coverage.
- Existing generated data: 44 HiFi directories, 27 MiB; no ONT directories.
- Available external environments include Clair3 with HiFi and ONT checkpoints,
  minimap2, samtools, bcftools, and a separate pbsim3 environment.

## Completion audit

- Existing truth inventory and fingerprinting complete; documented simulator
  metadata/name defects and validated actual sequences.
- Current 44 HiFi + 3 new ONT + 6 held-out HiFi full pipelines complete; 24
  additional depth/name attempts include six explicit coverage failures.
- Exact truth-only classification: 106/106 haplotypes; dictionary context
  probes: 132/132. Full sequence, complete repeat structure, exact mutation
  identity/position/haplotype, extra calls, missing output and count errors
  are scored independently, with descriptive statistics and uncertainty.
- Unit/static checks passed (239 tests, 86.11% coverage); real tool tests
  passed (10); end-to-end tests gave six passes/two expected failures, no skips.
  Documentation build passed after the guide update.
- Ten cached raw-call-fixed QUAL/consensus conditions cover all 47 baseline
  HiFi/ONT samples; default-QUAL5 replay matches all complete baseline results.
- Three paired performance repetitions establish 8.07x faster difficult
  classification with exact bit-vector scoring and identical full output.
- Primary-only ablations, anchor-mode inference, held-out validation, read-name
  controls, depth replicates and a rejected unlocalized template screen are
  documented with their limitations.
- Independent agent review checked headline statistics, output denominators,
  no-call accounting and performance claims against artifacts.

The final deliverable is
[`2026-09-14-validation-report.md`](../2026-09-14-validation-report.md), with four
companion reports and ignored reproducibility artifacts. Recommended algorithm
changes remain proposals; implementing and validating a new production caller
was not required by this testing/debugging/brainstorming objective. No commit,
merge, deployment or release was performed.
