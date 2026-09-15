
You are a principal computational genomics software engineer. Execute Wave 2 in
`${CHECKOUT}` using Superpowers. Resolve the verifiable
scope of issue #44:
https://github.com/berntpopp/MucOneSpan/issues/44

## Mission and authority

Build and actually run a reproducible benchmark of the available MUC1 data in
ENA PRJEB92208. Deliver an auditable baseline with honest event-detection,
sequence-reconstruction, execution and resource results before subsequent ONT
accuracy work. This prompt authorizes implementation and cohort execution,
review, a focused PR, and integration/release when the required checks pass.

Proceed through routine, reversible work without asking for reconfirmation.
Use parallel agents with disjoint file ownership where useful; the lead owns
integration and review. Make reasonable engineering choices and document them.
Ask only for genuinely missing information that blocks meaningful progress,
credentials, or costs beyond the available environment. Continue independent
work while waiting. Do not stop after writing a plan or running one pilot.

Follow the merge/release/cleanup instructions in section 8. Do not contact study
authors, upload patient reads/results, or post unrelated GitHub messages. Keep
raw human data and generated outputs outside Git. Commit only minimal accession,
truth-provenance and aggregate benchmark records needed for reproducibility.

## 1. Inspect, isolate and freeze the baseline

1. Read root/scoped `AGENTS.md`, `docs/development.md`, contributing guidance,
   the live issue body/comments, and existing evaluation/benchmark code and tests.
   Issue claims, proposed fixes and old line numbers are hypotheses.
2. Inspect status, branches, remotes and worktrees. Last verified main was release
   v0.14.1 at `c08000a0b97e9ee00e560a369ef320b4d06d3004`; refresh actual state.
   Use Superpowers using-git-worktrees to create an isolated branch/worktree from
   clean current main. Preserve any unrelated work and ignored scientific data.
   Do not restore rejected #49 experiments from the external archive.
3. Apply brainstorming/design, writing-plans, TDD, systematic-debugging,
   subagent-driven-development or executing-plans, code-review, and
   verification-before-completion as appropriate. Scope and execution are
   authorized here: record design decisions in `.planning/` and execute them.
4. Freeze caller commit, locked dependencies, tool/model versions and checksums,
   reference resources, settings, input manifest, preprocessing, and metric
   definitions before interpreting results. Capture a small existing scientific
   fixture panel to verify this benchmark work preserves caller semantics.
5. Inspect `scripts/evaluate.py`, `scripts/benchmark.py`, experiment runners and
   `src/muc_one_span/evaluation/`. The current evaluator expects strict MucOneUp
   simulation truth. Reuse execution/provenance utilities where suitable; add a
   distinct clinical event-truth adapter instead of manufacturing full simulated
   haplotype truth from a diagnostic label or weakening the existing contract.

## 2. Verify the dataset and truth before scoring

Primary sources:

- ENA: https://www.ebi.ac.uk/ena/browser/view/PRJEB92208
- Paper: https://doi.org/10.1038/s41598-025-30441-3
- Accessible paper PDF: https://www.nature.com/articles/s41598-025-30441-3.pdf
- Comparator: https://github.com/DHmeduni/VNTRPipeline

An ENA API inventory checked on 2026-09-15 contained 20 ONT runs, no PacBio
runs, and 11 runs labelled MUC1: nine AMPLICON and two WGS. Their ENA-generated
compressed FASTQ files totalled 198,156,449 bytes. This is a metadata snapshot,
not proof of input readiness. Refresh it and explain changes. The issue's claim
that this accession includes PacBio is unsupported by that snapshot.

Initial MUC1 accession checklist:

| Run | Alias | ENA strategy |
| --- | --- | --- |
| ERR15277562 | HG001_PCR_MUC1 | AMPLICON |
| ERR15277563 | HG002_PCR_MUC1 | AMPLICON |
| ERR15277564 | HG003_PCR_MUC1 | AMPLICON |
| ERR15277565 | HG004_PCR_MUC1 | AMPLICON |
| ERR15277566 | MP1_PCR_MUC1 | AMPLICON |
| ERR15277567 | MP2_PCR_MUC1 | AMPLICON |
| ERR15277568 | MP3_PCR_MUC1 | AMPLICON |
| ERR15277569 | MP4_PCR_MUC1 | AMPLICON |
| ERR15277570 | MP5_PCR_MUC1 | AMPLICON |
| ERR15277552 | HG002_WGS_MUC1 | WGS |
| ERR15277553 | MP1_WGS_MUC1 | WGS |

The paper describes four known positive controls and an additional diagnostic
case, MP5. HG002 has a referenced Q100 sequence benchmark. Repeated runs and
related participants require explicit accounting. Verify the exact orthogonal
confirmation for each scored claim in the paper/supplements and cited sources.
Do not assume Sanger confirmation, complete haplotype truth, or negative-control
status merely from an alias or the issue's wording.

Required inventory and truth ledger:

- Preserve study/sample/experiment/run accessions, aliases, biological identity
  where published, family/replicate relationships, platform/model, library
  strategy, basecalling provenance where available, URLs, sizes and checksums.
- Explicitly classify every study run as primary MUC1 amplicon, secondary MUC1
  WGS-derived input, or excluded non-MUC1, with the reason. Inspect whether the
  small WGS files already contain locus-selected reads; do not describe these as
  complete genomes or benchmark end-to-end genome processing without evidence.
- Separate confirmed event presence/absence, independently established full
  sequence, published comparator-only observations, and unknown truth. Every
  truth assertion must have a source location and stated resolution/limitations.
  MP5 remains exploratory unless independent confirmation is established.
- Record full-sequence truth only where an independently sourced versioned
  sequence supports it, with coordinates, orientation, repeat boundaries and
  checksums. Do not infer full VNTR sequence from a figure, repeat count or
  diagnosis. Resolve repeat-index conventions before assessing event positions.
- Preserve unknowns explicitly. A candidate caller/comparator output cannot
  serve as independent truth for itself. Healthy/reference labels alone do not
  prove absence of every target event. Do not invent a PacBio arm; document its
  absence and any resulting issue-scope limitation.

## 3. Divide implementation into reviewable units

Agree a small versioned manifest, truth-ledger and run-record interface first.
Record exact filenames and field contracts in the implementation plan.

- Worker A: ENA metadata inventory, checksum-verified downloader/preprocessor,
  and focused tests for those files.
- Worker B: event-truth validation, scoring/aggregation and focused tests; owns
  no execution or download code.
- Lead: baseline runner, existing-tool integration, measurements, cohort
  execution, results documentation, shared interfaces and combined validation.
- Use a separate reviewer after integration. If delegation is unavailable,
  execute these units sequentially with the same review gates.

Prefer focused helpers following current repository structure. Do not have
multiple workers edit a shared evaluator or report module concurrently. Preserve
legacy evaluation imports and schemas; avoid an unrelated framework rewrite.

## 4. Download, preprocess and run the complete eligible cohort

Implement a portable documented entry point that retrieves ENA metadata,
downloads selected inputs and checks expected byte counts and published MD5s;
record local SHA-256 as provenance. Use explicit configuration for data/result
roots, resource limits, tools and models. No committed machine-specific paths.

Downloads must be resumable or safely restartable, use temporary names until
validated, detect corruption and stale files, and never reuse mismatched content
as success. Test metadata failures, partial downloads, checksum mismatch,
duplicate accession/identity handling and successful offline replay with small
synthetic fixtures. Unit tests must not require ENA availability.

Inspect read headers/formats, primer and flank geometry, orientation, basecalling,
read lengths and evidence of prior filtering. Define deterministic preprocessing
from input properties before looking at diagnostic outcome. Preserve read IDs
and counts through every filter; record exclusions and exact commands. Do not
crop away informative bases, select reads by expected mutation, or convert noisy
ONT reads into an unsupported assay assumption to make a positive call appear.

Run a smoke case to validate wiring, then all eligible MUC1 amplicon inputs.
Also attempt all MUC1 WGS-derived inputs using an explicitly documented supported
input path, reporting this arm separately. If assay/format constraints prevent a
valid run, record the reason and retain it in the inventory; do not silently
force amplicon semantics or drop difficult cases. Leave unmet scope visible.

Capture per run: input/settings/resource hashes, caller/tool/model versions,
exact argument lists, start/end/status/exit code, logs, requested and observed
resources, wall time, peak memory with method/units, read counts, output hashes,
and actual analysis/evidence state. Bound subprocess time and concurrency using
existing tool abstractions. Handle interruption and prevent stale outputs from
being scored as a new successful run. Resume only after provenance validation.

Use available tools/models after verifying them. The previous session validated
minimap2, samtools, bcftools and Clair3 locally. It also found malformed embedded
VCF output with igv-reports 1.16.0 and verified 1.13.0. Recheck compatibility;
record report failures separately from valid scientific output, and preserve the
actual pipeline failure status. Do not weaken report validation to get green runs.

## 5. Score what the evidence supports

Predeclare detection endpoints and denominators before inspecting predictions:

- Confirmed-event recovery: exact event identity where truth supports it;
  separately assess position, allele assignment and full sequence when known.
- Negative event specificity only for independently established negative truth
  at the endpoint being evaluated. If no justified negative set exists, report
  specificity as not estimable with an explicit reason.
- Sequence/repeat-length concordance only for independently established sequence
  truth; normalize only documented coordinate, orientation and allele-order
  conventions. Preserve biological indels, ambiguity and novel motif differences.
- Execution success, callability and insufficient-evidence rates across the
  entire eligible inventory. Failed/no-call/unattempted cases must remain in
  accounting and must never become true negatives. Show confirmed positives
  missed because of failure/no-call in all-sample recovery, and label conditional
  callable-only metrics separately.
- Count unique biological samples separately from libraries/runs; show family
  relationships and technical replicates. Avoid claiming independent observations
  from PCR read counts or correlated individuals. Any confidence interval must
  state its assumptions and limitations for this small, related cohort.

Keep machine-readable per-run and aggregate results, numerator/denominator
counts, truth categories, exclusions and source references. A metric of zero is
valid; absent or invalid evidence must not be silently converted to zero or
success. Test matched events, wrong event despite positive alarm, missing truth,
unknown negatives, no-calls/failures, malformed/stale outputs, duplicates,
allele swaps and equivalent event representations where applicable.

Report discrepancies by stage: input/assay compatibility, read evidence/length,
calling, consensus/repeat classification, clinical display, and truth ambiguity.
This is a frozen baseline: do not tune thresholds, implement #20/#21/#47, adopt
POA, change clinical tiers, or enable rejected #49 settings. Preserve failed
cases as evidence for subsequent work. Correct harness bugs with regression
proof; version and rerun affected results if the harness changes.

## 6. Measure resources and document limitations

Measure analysis wall time and peak memory across the executed cohort, with
hardware, thread count, tool versions, input sizes and measurement scope. Separate
download/preprocessing from calling and report generation. State whether memory
measurement covers the process tree and how concurrent children are handled;
label an estimate or lower bound honestly. Use a consistent prespecified repeat
policy if reporting median performance rather than a single observed run.

Inspect the paper, supplements and versioned comparator documentation for
published runtime/memory measurements. Quote only supported figures and identify
hardware, input and stage differences. If a requested published measurement is
absent, say so; do not invent it or claim a matched speedup. A local comparator
run is optional additional evidence and must be pinned and fairly measured;
never present it as independent truth or silently expand this task into a
comparator rewrite. Document all unavailable comparisons explicitly.

Write `docs/guides/clinical-validation-results.md`, integrate it into the existing
docs navigation, and update the changelog. Include cohort inclusion/exclusion,
truth resolution, reproducible commands, settings/environment, per-run and
biological-sample tables, error analysis, resource results and limitations.
Separate historical simulation claims from these observations. A tiny benchmark
cannot establish universal clinical sensitivity, specificity or diagnostic safety.

## 7. Engineering, verification and review gates

- Python 3.10+; maximum 649 physical lines per authored code/config/template file.
- Keep CLI compatibility, scientific algorithms, default settings, output schemas
  and bundled resources stable except explicit additive benchmark interfaces.
- Keep subprocesses in the existing tool abstraction with argument lists and
  visible bounded failures. Do not add shell interpolation.
- Add meaningful deterministic unit tests and bounded integration tests for
  changed behavior. Use tiny synthetic fixtures rather than committing patient
  reads. Update `pyproject.toml` and `uv.lock` together if dependencies change.
- Run focused tests while developing, then `make ci-check`, `make test-int`,
  and `make docs-check`; run `make security-check` for dependency/security
  changes and `make build-check` for packaging/resource or release changes.
  Run browser checks if report rendering is changed. Report actual pass/fail/
  skip counts and missing prerequisites; a skipped behavior is not verified.
- Verify the existing representative fixture panel produces unchanged scientific
  results. Normalize only documented nonbiological metadata such as paths and
  timestamps. Do not require full 500-sample reruns merely for a separate harness,
  but rerun the appropriate existing evaluator regressions and any affected
  science panel if shared behavior changes.
- Obtain independent code/scientific review of identity mapping, truth circularity,
  denominators, failure handling, reproducibility and claims. Resolve substantive
  findings and review the final diff, generated evidence and Git status.

## 8. PR, merge, release and cleanup

Commit only intended code, tests, small curated provenance/results, docs and
changelog after required checks pass. Create a focused PR with the concrete
problem, actual cohort/results, validation and remaining limitations. Use
`Closes #44` only when its achievable scope and acceptance criteria are fulfilled
and unsupported original assumptions are explicitly reconciled. Otherwise use
`Related to #44`, state remaining work, and keep the issue open. Do not claim
PacBio validation or independent truth that the available data do not supply.

Watch all applicable required GitHub Actions on the latest PR revision. Fix
failures and rerun; do not merge while required checks are pending, failed or
unexplained skips. Once required checks pass and review findings are resolved,
merge using repository conventions. This prompt authorizes that merge.

Determine the appropriate patch or minor bump from the actual delivered public
behavior and repository conventions: documentation/benchmark corrections may
justify patch; an additive supported benchmark feature may justify minor.
Do not guess the next version from this prompt's historical baseline. Prepare
and validate the version/changelog/lock updates as required, integrate through
the repository workflow, and create an annotated version tag and GitHub release
on the verified release commit. Verify applicable main/tag/release Actions and
published artifacts. Release notes must describe supported findings and material
limitations without turning this cohort into a broad clinical-performance claim.

Verify the final GitHub state of #44 and report it. Do not silently close other
issues; list any already implemented but still-open issues as closure candidates.
Remove merged task branches and temporary worktrees only after verifying their
commits are integrated and work is preserved. Inspect other stale branches/files,
archive any valuable uncommitted work with a manifest/checksums before cleanup,
and delete only confirmed task-owned or demonstrably redundant leftovers.
Preserve unrelated work, ignored datasets and the existing external archive.

## Final delivery

Finish with:

1. Issue #44 acceptance table: each requested item, implementation, actual
   cohort/validation evidence, and remaining limitation.
2. Cohort counts, truth eligibility, detection/callability and sequence metrics
   with denominators, failures and unestimable endpoints explicit.
3. Reproduction commands; manifest/truth/results/documentation paths; data and
   output locations; baseline commit and environment provenance.
4. PR/merge, version/tag/release links, Actions status, actual issue state, and
   remaining branches/worktrees/uncommitted files or preserved archives.
5. Evidence-based recommendations for the next #20/#21/#46/#47 wave, with
   concrete observed failures and no premature choice of scientific fix.

Complete the authorized work instead of offering to continue. If external data,
truth or infrastructure prevents an acceptance item, preserve the useful result,
state the specific limitation, and do not label an incomplete benchmark complete.
