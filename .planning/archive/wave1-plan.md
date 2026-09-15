# Wave 1 design and implementation plan

> Execute with Superpowers debugging, TDD, subagent-driven development and independent review. User explicitly authorizes implementation and disjoint parallel ownership; no commits, push, merge, publishing, comments or issue closure.

**Goal:** Fix #41 mapping lifecycle and #39/#40 report evidence, status and percentage display.
**Base:** clean main 1f6c165; branch fix/wave1-mapping-reports in sibling MucOneSpan-wave1. Main's #49 dirty diff and untracked hashes are captured in /tmp/wave1-main-* for final preservation verification.
**Spec:** user Wave 1 attachment; live issue bodies read with gh, all have no comments.
**Architecture:** streaming external pipeline behind tools abstraction with isolated process groups and one monotonic finite budget; report boundary resolves allele tracks, status provenance and explicit percentage units. CLI and pipeline pass these interfaces without modifying science.
**Constraints:** Python 3.10+, files <=649 lines, no new dependencies unless necessary, no scientific threshold/default/preset changes, explicit arguments, offline assets and safe autoescaping.

## Interfaces and design decisions

- Mapping: map_reads gains keyword timeout=DEFAULT_SETTINGS.run.mapping_timeout; run.mapping_timeout defaults 3600 seconds and validates finite >0. CLI map/run --mapping-timeout overrides JSON. One pipeline budget includes bounded termination, reaping and drain joins, reserving cleanup time within it; index also bounded by remaining mapping budget. Keep SAM streaming and cleaned PATH.
- Reporting: generate_report gains vcf_paths dict[str,Path] or None and execution_status dict or None. A supplied plural mapping takes precedence over singular vcf_path (including empty mapping). Preserve --vcf and add repeatable --allele-vcf LABEL PATH; use explicit allele paths from summary when neither option supplied. Deduplicate shared resolved paths, label shared evidence without claiming phase.
- Loci: use actual summary.alleles[key].contig_name; selected missing reference contigs must produce explicit diagnostic, never silently switch to contig_1.
- Execution: explicit status > summary status; standalone reads only explicit --run-status or summary's sibling run_status.json. A malformed known status fails visibly. Missing status retains legacy decision behavior with unavailable execution status. Full pipeline passes analysis_completed only after successful analysis; sidecar stays running until callback/report succeeds. Pipeline updates summary after rendering; decorator updates sidecar terminal state; failed rendering keeps failure semantics. Known failure/interruption/insufficient/running blocks reassuring negatives, retains recorded mutations with execution warning.
- Percent: exact_match_pct is always 0–100, including 0.5 and 1. Invalid/missing values display unavailable; no magnitude-based inference. Unicode em dash for missing/null/empty nomenclature, autoescape stays on.

## Tasks and test matrix

- [x] A mapping: reproduce unbounded cleanup; mocked success, producer/consumer failure/hang, missing second executable, interruption, descendant survival, stderr pressure; bounded real process handshakes; HiFi/ONT/override commands; partial output cleanup.
- [x] B report: red reproductions of metadata/second VCF, dash escaping and 0.5/1 percent; single/two/shared VCF, absent inputs, missing contig, singular compatibility, embedded/sidecar sessions and visible browser loci/tracks; hostile HTML remains escaped; status conflicts/legacy and percentages 0/0.5/1/50/92/100/invalid.
- [x] Lead integration: failing tests for timeout JSON/CLI propagation and status precedence/lifecycle; split report CLI helper from near-limit cli.py; successful full callback, failed report, insufficient/failed/interrupted standalone provenance and relocated output.
- [x] Scientific panel: baseline archived clean source with real tools on HiFi normal, HiFi dupC and ONT dupC; compare consensus, repeats, mutation identity, genotype/phase and allele evidence after fix, normalizing only output paths/timestamps.
- [x] Independent review: hand reviewer actual diff/new files and original spec; resolve substantive findings and verify fixes.
- [x] Gates: make test-fast, make quality, make ci-check (>=80% branch-aware coverage), make test-int with tool PATH/model, make docs-check, uv run --locked --all-extras pytest tests/browser --no-cov, make build-check for templates. Exact counts and skips recorded.
- [x] Delivery: docs/changelog, final diff and main preservation verification, issue evidence table and PR text; archive plan only after all requirements audited.

## Evidence ledger

- Baseline make dev passed; make test-fast: 760 passed in 5.41s, /tmp/wave1-baseline-unit.log.
- Worker mapping owns tools/process/mapping and focused tests; worker reports owns report/IGV/templates and focused/browser tests; lead owns settings/CLI/pipeline/status/docs; baseline worker owns panel script/evidence.

## Completion audit

All required Linux work and gates completed; independent review approved. Final
results: 874 unit passes; 91.32% branch-aware coverage; 69 integration passes
(eight e2e-marker tests deselected, no skips); 11 browser passes; quality, docs
and distribution checks passed. Three final real pipelines preserved 84/84
scientific hashes after the last supervisor correction. Main #49 diff and nine
untracked hashes remain identical. Detailed issue evidence, exact commands, PR
text and explicit platform/upstream-tool limitations are in ../wave1-delivery.md.
