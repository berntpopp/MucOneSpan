# Wave 2 clinical benchmark implementation plan

> **For agentic workers:** Use superpowers:subagent-driven-development with disjoint ownership and independent review.

**Goal:** Execute and publish an auditable frozen PRJEB92208 benchmark, then integrate and release after validation.
**Architecture:** Add clinical inventory/preprocessing, event-truth/scoring, and supervised execution modules. Preserve simulation evaluation and scientific defaults.
**Tech Stack:** Python 3.10+, stdlib, existing tools.py, locked uv environment.
**Spec:** .planning/wave2-spec.md (full user requirements).

## Global constraints

Maximum 649 physical lines per authored code/config/template; no raw human data in Git; no scientific algorithm/default changes; retain all eligible runs in denominators; independent truth only. User authorizes design, execution, parallel work, PR, merge, release without reconfirmation.

## Interfaces (schema_version = 1)

JSON objects, UTF-8, finite numbers, explicit null for unknown values.
Inventory: {schema_version, study_accession, retrieved_at, source_url, runs: [...]}.
Each run preserves ENA metadata as `ena` mapping and has `run_accession`, `sample_accession`, `experiment_accession`, `alias`, `biological_sample`, `arm` (primary_amplicon, secondary_wgs, excluded), `inclusion_reason`, `platform`, `files` [{url, bytes, md5}]. Keep original ENA sample identifiers; biological identity is separately sourced by truth ledger. Replicated identities permitted and grouped; duplicate run accessions rejected.
Preparation record: {schema_version, run_accession, input_sha256, output_sha256, input_reads, output_reads, excluded_reads, read_lengths, headers, preprocessing, output_path, files: [{url, bytes, md5, sha256, path}]}. Identity-preserving validation/decompression only initially; no diagnostic-driven filters.
Truth ledger: {schema_version, sources: [{id, url, location, ...}], samples: [{biological_sample, run_accessions, relationships, event_truth, sequence_truth, limitations}]}.
`event_truth`: list of {event, status: present|absent|unknown, source_ids, resolution, ...}; event identity uses dictionary full mutation name where independently supportable (e.g. dupC); no inferred position/allele. Sequence truth null unless independently versioned sequences and boundary provenance obtained.
Run record: {schema_version, run_accession, biological_sample, arm, status, exit_code, started_at, ended_at, argv, provenance, measurements, outputs, analysis_state, report_state, events, error}. `status`: completed|insufficient_evidence|execution_failed|invalid_artifacts|unsupported_input|unattempted. `events` list {event, allele, repeat_index, ...}; absent valid evidence uses null, not empty list. Outputs include file hashes. Scorer only grants event credit for completed callable runs with validated artifact provenance. Separate recovery from failed-report scientific evidence descriptions.

## Task A — Inventory and preparation (Worker A)

Files: src/muc_one_span/clinical_data.py, tests/unit/test_clinical_data.py (split focused siblings if needed).
Interfaces: fetch_inventory(study: str) -> dict; validate_inventory(data: dict) -> dict; prepare_run(run: dict, data_root: Path) -> dict. Expose offline metadata conversion function and checksum-verified restartable downloads. Actual ENA fields determined from API; share schema adjustments with lead before changing consumers.
- [ ] Write failing deterministic tests: malformed metadata, duplicate accessions, repeated biological identity, truncated download, byte/MD5 mismatch, stale local file, successful offline replay, invalid FASTQ, preserved IDs/counts.
- [ ] Implement validated metadata and safe .part downloads with local SHA256 and bounded requests; preserve exclusions and classify every run.
- [ ] Implement deterministic format/read-length/header inspection and identity-preserving decompression.
- [ ] Run focused tests and Ruff; lead reviews interface/diff before integration.

## Task B — Truth validation and scoring (Worker B)

Files: src/muc_one_span/clinical_scoring.py, tests/unit/test_clinical_scoring.py.
Interfaces: validate_truth(data: dict, inventory: dict) -> dict; score_cohort(inventory: dict, truth: dict, records: list[dict]) -> dict.
- [ ] Write failing tests for exact event vs wrong positive alarm, absent/unknown truth, no-call/failure/unattempted, duplicate run records, stale/invalid artifacts, biological vs run counts, allele swaps and sequence equality if sequence truth exists.
- [ ] Implement strict separate clinical adapter and all-inventory denominators; callable-only endpoint separate; specificity null without confirmed negatives; sequence endpoint null without independent sequence truth.
- [ ] Require source-backed truth and explicit endpoint resolution; preserve truth ambiguity.
- [ ] Run focused tests/Ruff; lead reviews scientific contract.

## Task C — Freeze, supervise, execute (Lead)

Files: src/muc_one_span/clinical_runner.py (split provenance/observation helpers), scripts/clinical_benchmark.py, tests/unit/test_clinical_runner.py, tests/integration/test_clinical_benchmark.py, benchmarks/clinical/prjeb92208/*.json.
- [ ] Refresh issue, ENA, paper/supplements, comparator; independently verify truth provenance and source locations; record unsupported truth and PacBio scope explicitly.
- [ ] Freeze caller c08000a0b97e9ee00e560a369ef320b4d06d3004, lock/resources/config/model/tool hashes, hardware; save representative scientific fixture panel before harness implementation.
- [ ] Add failing runner tests for status, timeout/interruption, stale artifacts, resume provenance and real tiny supervised command integration.
- [ ] Implement CLI inventory/prepare/run/score, existing tool abstraction with bounded process supervision; unique output attempts and provenance validated resume; stage time and peak memory with stated scope.
- [ ] Inspect all input geometry and WGS origin; freeze preprocessing/endpoints before prediction review. Run smoke then all nine amplicons and attempt two WGS-derived inputs honestly.
- [ ] Validate outputs, retain all failures, aggregate minimal results and source-backed ledger. Do not commit raw reads/sequences from participants.

## Task D — Review and deliver (Lead + independent reviewer)

Files: docs/guides/clinical-validation-results.md, mkdocs.yml, CHANGELOG.md, version.py and uv.lock as release convention requires.
- [ ] Recheck fixture panel unchanged; run make ci-check, make test-int, make docs-check, make build-check and applicable security gate. Report all skips.
- [ ] Independent review of code, source truth, denominators, evidence and final claims; resolve substantive findings with regressions and rerun affected cohort records if harness changes.
- [ ] Commit intended artifacts; focused PR Related to #44 unless full achievable acceptance reconciled; inspect applicable latest-revision Actions; merge after passing.
- [ ] Select version from delivered supported feature; release workflow with annotated tag; verify main/tag/release Actions and artifacts.
- [ ] Audit actual issue status, related closure candidates, branches/worktrees; preserve unrelated archives/data; remove only integrated task-owned work.
- [ ] Final acceptance table, denominators, reproduction/provenance paths, PR/release/actions/issue state, next-wave evidence recommendations.

## Design decisions / alternatives

A standalone script would duplicate provenance and scoring; modifying the strict simulation truth schema would weaken its scientific contract. Additive package helpers with a portable script are selected. Single observed run per library, serial cohort execution, fixed threads and timeout; no median or confidence interval claims. Raw ENA FASTQ is validated/decompressed without sequence alteration unless input geometry proves a supported deterministic necessity. Unknown truth is not negative. No POA, threshold or clinical-tier tuning.

## Progress and validation

Initial main clean at c08000a; refreshed origin unchanged; isolated .worktrees/wave2-clinical. Task-local worktree exclusion in .git/info/exclude preserves tracked main. Baseline tests running; external data/results root will be explicit and outside Git.


## Current delivery steering (2026-09-16)

The user subsequently requested root-cause debugging, actual VNTRPipeline
comparison, an explicit MP1 classifier fix, issues, and a **draft PR with current
code**. Preserve the original frozen results and keep this delivery unmerged.
The separate fix worktree's exact-context dictionary correction is integrated
into the draft branch; end-to-end MP1 now recovers supported B:dupC at repeat17.
All22 frozen consensus mutation lists were replayed, with only MP1allele2 changed.
993combined unit tests pass (89.62%branch-aware). Issues52–56 track the fixed and
remaining defects. Full comparator execution and a complete post-fix cohort are
not yet validated. Release/version/merge/cleanup remain pending draft review.
