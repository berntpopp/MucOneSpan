# Wave 1 independent review

Reviewed against clean base `1f6c165` and the complete user attachment, using the
current tracked diff and new helper/test files. Review covers the mapping process
lifecycle, report provenance, VCF/locus plumbing, percent presentation, and the
scientific preservation evidence. No production files were edited by reviewer.

**Final verdict: approve specification compliance and code quality for the
reviewed Linux implementation, with igv-reports 1.13.0 verified externally.**
All substantive review findings were addressed. This verdict retains the
documented non-Linux reaping and igv-reports 1.16.0 limitations and is separate
from the lead's final combined validation results.

## Initial findings

| ID | Priority | Finding | Required resolution |
| --- | --- | --- | --- |
| R1 | High | `tool_pipeline._cleanup` only waits on direct `Popen` children. Orphaned grandchildren are killed but not reaped. Real integration `_is_running` explicitly accepts `/proc` state `Z`, so its passing assertion cannot establish the user's descendant reaping requirement. | Provide isolated descendant adoption/reaping on Linux and assert recorded descendant PIDs disappear. Avoid process-wide subreaper changes in the library caller or reaping unrelated children. State platform limits precisely. |
| R2 | High | The pipeline persists `summary.run_status=analysis_completed` before rendering. If rendering fails, its sidecar changes to `execution_failed` while summary remains successful for report decision purposes. Moving that summary without its sidecar permits a reassuring negative. | Persist final failure/interruption/insufficient provenance into existing summaries as well as the sidecar; regression must regenerate a standalone report from a relocated summary without the sidecar. Never mark callback completed before rendering succeeds. |
| R3 | High | `map_reads(input_path=out/'mapping.bam', output_dir=out)` now unlinks the source input in the initial artifact cleanup before FASTQ conversion. Previously BAM conversion happened first. This is a data-loss regression. | Reject source/output collisions before any mutation, or safely preserve conversion-first behavior. Cover the input's bytes surviving rejected invocation, including resolved path aliases. |
| R4 | Medium | `_cleanup` can raise before the handler builds diagnostics, hiding the original failure and both captured diagnostic tails behind a generic cleanup error. | Retain the primary failure and available tool diagnostics when cleanup is also incomplete. Add a deterministic cleanup-failure regression. |
| R5 | Low | `call --mapping-timeout` was added but its callback never uses that argument. | Remove the accidental unused flag/argument; retain the intended `map` and `run` configuration paths. |

The lead accepted R1–R5 and assigned fixes. R1–R5 and both supervisor follow-up
findings are resolved after reinspection and focused verification. No
substantive review blocker remains in the final reviewed diff.

## Reproductions and checks

- An isolated Python command using `unittest.mock.patch` around the external
  mapping boundary created `out/mapping.bam` containing original read bytes and
  invoked `map_reads` with that same input/output location. It printed
  `INPUT_EXISTS_WHEN_CONVERSION_STARTS False` and
  `INPUT_EXISTS_AFTER_FAILURE False`, directly reproducing R3.
- Static control-flow inspection establishes R4: `_cleanup(...)` precedes
  diagnostics construction with no exception handler. An injected cleanup
  failure exposed only `RuntimeError cleanup deadline expired`.
- Reviewer ran
  `uv run --locked --all-extras pytest tests/unit/test_report_wave1.py tests/unit/test_report_igv.py tests/unit/test_wave1_cli.py --no-cov -q`.
  At that moment the report worker had added three red edge regressions:
  **77 passed, 3 failed** (huge integer percentage and malformed list/dict status).
  These are ongoing TDD failures, not accepted completion evidence. Rerun after
  worker fixes is required.
- `git diff --check` was clean at initial inspection.
- After the lead and report edge fixes, the reviewer reran the focused command
  with `tests/unit/test_run_status.py` added: **89 passed, zero failures**.
- Reviewer independently exercised failure, interruption and insufficient
  evidence after an `analysis_completed` summary was written, copied each
  summary alone into a new directory, and rendered a standalone report. All
  three terminal statuses persisted and all three reports were inconclusive,
  with no sidecar present in the relocated directory. This verifies R2's
  portable-summary failure mode beyond merely asserting a stored JSON field.
- R5 reinspection confirmed only `map` and `run` retain the new timeout flag;
  `call` no longer exposes the unused option/argument.

## Findings that did not require changes

- Allele loci are selected from sorted `summary.alleles[*].contig_name`, with
  duplicate contigs removed. Requested missing reference contigs raise explicitly.
  Existing coordinate arithmetic remains unchanged.
- Plural VCF input overrides singular input, including an empty mapping. Paths
  are resolved, sorted by allele label, and deduplicated with a shared track name.
  Requested missing files raise; optional absent tracks remain optional.
- The standalone CLI reads status beside the input summary or from an explicit
  status argument, independently of the output HTML directory. Supplied status
  overrides summary status. Known incomplete execution preserves mutation rows
  while preventing a negative decision. Legacy behavior remains explicit.
- The template uses Unicode dashes and retains Jinja autoescaping. Report
  percentages preserve the producer's 0–100 units; text, width and ARIA agree for
  the required 0, 0.5, 1, 50, 92 and 100 cases.
- No threshold, mapping preset, genotype/phase selector, repeat identity or
  mutation classification changes appeared in the reviewed production diff.
  The recorded scientific panel covers three real full pipelines and compares
  84 artifacts, with no differences after documented normalization. This review
  inspected the runner/evidence scope; the lead owns the actual external gates.

## Review scope and limits

The lead owns the full final gate record, which must distinguish real
external/browser passes from skips. Reviewer checks are scoped and do not
substitute for the requested combined CI, external tool, browser, documentation,
packaging, and final scientific panel checks. Non-Linux orphan reaping remains
delegated to the OS with an explicit warning, preserving POSIX mapping
compatibility. Cleanup that cannot complete within its budget is reported
explicitly; uninterruptible kernel operations cannot be made reliably bounded
by user-space code.

## Supervisor follow-up review

The first isolated Linux supervisor implementation correctly confines the
subreaper flag and `waitpid(-1)` to its own process tree. It retains a POSIX
group-cleanup fallback elsewhere and makes the weaker orphan-reaping guarantee
visible. Two further issues were sent to the mapping worker and lead:

- Retrying `_collect` after interruption loses partial control-channel bytes,
  parsed result, EOF state and stderr tail. A deterministic reviewer reproduction
  supplied half of a result JSON message, injected `KeyboardInterrupt`, and
  supplied its remainder on retry: the observed result was `JSONDecodeError`.
  Keep collector state across retries.
- The caller's initial collector consumes the entire deadline before its
  emergency branch, leaving `wait(timeout=0)` unable to guarantee that even the
  supervisor is reaped. Reserve caller cleanup inside the one overall deadline,
  use earlier supervisor/collection finish cutoffs, and explicitly report any
  incomplete reaping. Do not add time beyond the original deadline. Protect
  bounded descendant cleanup from a retry signal interrupting the reaper itself.

Both follow-up findings were corrected and verified below.

### Reinspection and real process verification

The persistent `_Collector` now retains partial protocol bytes and result/EOF
state across retry. A regression based on the reviewer's exact split-message
reproduction passes. The supervisor confines orphan adoption/reaping to its own
process tree and suppresses a further SIGINT while reaping adopted descendants.
Input/output collision rejection precedes artifact deletion and preserves the
source bytes. Cleanup errors include the primary failure and both tool tails.

Reviewer ran
`uv run --locked --all-extras pytest tests/unit/test_tool_pipeline.py tests/unit/test_tool_supervisor.py tests/unit/test_mapping_lifecycle.py tests/integration/test_tool_pipeline.py --no-cov -q`:
**53 passed in 6.70 seconds**, including all ten real-process cases. The real
tests now reject zombies and require recorded tool/descendant PIDs to disappear.
The subsequently added cleanup-error diagnostic regression also passed.

An additional reviewer real-process check launched a descendant with its own
new session while retaining stderr, then failed its direct parent. The pipeline
returned the producer error after **1.981 seconds of a 2-second budget**, the
escaped-session descendant PID disappeared, and a separately launched unrelated
child remained alive. The reviewer cleaned up that unrelated child afterward.

One reserve issue was found during this reinspection: retry `_collect` could
still run until the final deadline, exhausting the initial caller reserve before
emergency kill/reap. An advancing-clock reproduction printed
`REAP_WAIT_REMAINING 0.0`. The worker was asked to retain a kill/reap reserve
after both collection cutoffs, with a clock-advancing regression.

The final implementation now separates supervisor completion, retry collection,
and emergency reaping cutoffs within the original deadline. The new regression
advances the clock through both collection timeouts and asserts positive final
kill/reap time. Reviewer independently ran final
`tests/unit/test_tool_supervisor.py tests/integration/test_tool_pipeline.py`
with locked all-extras pytest and no coverage: **25 passed in 6.62 seconds**
(15 unit, ten real process tests). Final `git diff --check` was clean.

## Final report guard review

Reinspected the guard added after real igv-reports 1.16.0 was observed joining
the first variant to `#CHROM`. The guard compares original versus decoded
embedded VCF headers, rejects malformed generated evidence, and removes the
malformed generated output. It does not rewrite source variants. The external
1.16.0 limitation and verified 1.13.0 version remain explicit in documentation.
Real-session tests now inspect the actual variant record and browser tests
query the live IGV feature coordinates, reference and alternate bases, in
addition to checking visible locus/track labels. No additional substantive
report finding resulted. Reviewer ran the new malformed-payload regression
plus `test_report_wave1.py`: **40 passed, zero failures**.
