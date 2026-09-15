# Wave 2 handoff — paused at the user's request

The user requested issues, a **draft PR with the current code**, and a pause to
continue in another session. **Do not merge or release during this handoff.**

## Workspace

- Main remains at `c08000a0b97e9ee00e560a369ef320b4d06d3004` (v0.14.1).
- Draft branch: `feat/clinical-benchmark`; worktree: `.worktrees/wave2-clinical`.
- Separate branch/worktree: `fix/dupc-repeat-context` /
  `.worktrees/fix-dupc-repeat-context`. Its uncommitted correction, tests and docs
  were copied into the draft branch. Preserve and inspect before cleanup.
- External data root: the main checkout's sibling `MucOneSpan-wave2-data`, called
  `${DATA_ROOT}` below. Raw reads, full participant sequences, generated outputs
  and diagnostic scripts remain there, outside Git.
- Initial harness commit: `1deb5dcfd1dbeede0fbd3ca859cc9a6aab011d47`.
- No version bump, merge, annotated tag or release has occurred.

## Frozen baseline

All 20 ENA runs are inventoried: nine MUC1 amplicons, two WGS-labelled inputs,
nine ACAN exclusions, no PacBio. All 11 selected archives were validated against
ENA size/MD5; compressed total: 198,156,449 bytes.

Final `cohort-v2`: 11/11 executions completed, 6/11 callable (amplicon 5/9,
WGS 1/2), five ambiguous. Publication-reported control named `dupC` recovery is
0/4; supported recovery 0/4; callable-only recovery 0/1. Independently confirmed
diagnostic sensitivity and specificity are not estimable. The patient WGS
identity remains unresolved.

HG002 PCR: one of two literal alleles exact, lengths 3,900/900 bp. HG002 WGS:
both literal alleles exact, 3,900/4,638 bp. Both runs are ambiguous. Across the
two libraries: literal pair agreement 1/2; strict callable pair recovery 0/2.
Do not describe the strict gate failure as incorrect WGS DNA sequences.

Final scores equal engineering v1 scores; all 306 scientific artifact hashes
agree after documented normalization. Final cohort runtime: 847.38 seconds;
two requested threads, serial CPU execution, 1,800-second per-run budget.
Memory is a largest-process RSS lower bound. The original three-fixture HiFi/ONT
panel matched all 84 scientific hashes before the separately authorized fix.

## MP1 fix included in the draft

Issue #52: a correct 61-base B-repeat C-tract insertion was missing from exact
mutation templates. Fallback kept a 60-base window on an edit-distance tie and
shifted its terminal A into the following repeat.

The production correction changes only bundled `dupC.allowed_repeats` and its
citation: 12 exact 60-base contexts containing the same G + seven-C + non-C tract.
No Python algorithm, fuzzy matching, threshold or model changes. There are 26
new synthetic tests. Independent review found no blockers; an additional 1,156
normal ordered motif pairs produced zero dupC calls. This is not a clinical
specificity estimate.

Reclassification of all 22 frozen consensuses changes only MP1 allele 2's
mutation annotation list. The other 21 mutation lists are unchanged; full-field
equality was not claimed. A separate full MP1 pipeline rerun completes with
supported **B:dupC at repeat 17**. See `${DATA_ROOT}/mp1-fixed-end-to-end` and
`benchmarks/clinical/prjeb92208/mp1-classifier-fix.json`. Do not replace the
original 0/4 baseline with an inferred complete post-fix cohort rate.

## Issues and root causes

- #52: classifier correction included; closes only on a future merge.
- #53: distinct-length ONT candidates force GT1 despite unphased evidence.
  MP3 diagnostic GT2 replay yields dupC; choosing GT2 is not validated phasing.
- #54: MP2/MP4 long full-span reads enter the shorter major allele's BAM, while
  fragment-derived valleys become allele 2. Detector replay reproduces selection.
- #55: clinical displays exceed evidence, including an ambiguous in-frame HG002
  WGS change displayed as PATHOGENIC and ambiguous PCR displayed as negative.
- #56: complete the local VNTRPipeline comparison after compatibility repairs.
- #44 remains open: no PacBio, independently mapped diagnostic negative set, or
  published runtime measurements supporting the original requested comparison.

Details: `wave2-miss-debug.md`, `wave2-error-analysis.md`, and external
`debug-misses/{fix-replay,read-lineage,selection-replay}.json`.

## Comparator — incomplete and stopped

Source v1.0 (`f7f594e74ce9cde273deb269f06bb528f0586997`) was mounted read-only into
image `sha256:888583f8ef0b69af9c5b386c51c18a2f6a73b0d46dac06424ebd56b01845d5a8`.
The runtime was built in March 2026 and has newer dependencies: this is a
source-pinned comparison, not exact paper-environment reproduction. Clinical
model checkpoints match our selected model. Network was disabled, with two CPU
and 16-GiB limits; internal thread requests were unchanged. Shell/Perl files
required executable-bit adaptation; source bytes were preserved.

- `comparator-v1/ERR15277566`: packaging preflight failed with permission exit 126.
- `comparator-v1/cohort/ERR15277566`: original-header attempt exited zero without
  valid final FASTAs. Minimap2 `-y` copied ENA UUID comments into invalid SAM tags.
- `sources/header-comment-probe/result.json`: four-arm synthetic confirmation.
- `comparator-inputs`: normalized comments; original first-token IDs, bases,
  qualities and order preserved. Full original headers are retained locally.
- `comparator-v1/normalized/ERR15277566`: stopped at the user's request during
  downstream processing; partial outputs are not validated comparative results.
- Independent plotting blocker: missing R BiocManager and an unset CRAN mirror.
  Provision an explicitly pinned runtime dependency before claiming completion.

The comparator container has stopped and been removed. No other comparator
libraries have run. External drivers: `run_comparator.py` and
`prepare_comparator_headers.py`. They refuse existing attempt directories; use a
new output root for reruns. Read `wave2-comparator-audit.md` and the diagnostic
plan before resuming. Do not infer success from the wrapper's exit code alone.

## Validation

- `make ci-check`: 993 passed; 89.62% branch-aware coverage; quality gates passed.
- `make test-int`: 70 passed, eight deselected by the integration marker;
  external igv-reports 1.13.0 was selected.
- `make docs-check`: strict build passed.
- `make build-check`: distributions, isolated installation, CLI and resources passed.
- Earlier unchanged-dependency audit passed; earlier browser suite: 11 passed.
- Refresh GitHub Actions on the draft PR; do not assume they passed.
- Check log hashes: `benchmarks/clinical/prjeb92208/validation.json`.
  Local logs: `${DATA_ROOT}/draft-*.log`.

## Resume

1. Read the latest user instruction and draft PR; inspect CI, issues and worktrees.
2. Preserve all original data, cohorts and the current fix; do not resume release
   work without the user's requested continuation.
3. Complete the comparator and its artifact checks if still requested.
4. Design separate corrections for #53–55 from demonstrated causes. Diagnostic
   GT2 selection or lower thresholds are not validated scientific solutions.
5. Run appropriate post-fix cohort checks before claiming new performance.
6. After resumed review: version/release/integration gates, issue-scope audit,
   and task-owned worktree cleanup. Original scope is in `wave2-spec.md`;
   the latest pause/draft instruction takes precedence.
