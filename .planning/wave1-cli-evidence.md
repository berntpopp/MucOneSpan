# Wave 1 CLI and execution provenance evidence

## Reproductions and interface checks

- Baseline main RunSettings rejected mapping_timeout as an unknown constructor argument;
  the new setting validates finite positive numbers and rejects booleans.
- New CLI suite initially had 10 failing tests: map/run timeout values not propagated,
  sibling failure/insufficient/interrupted/running sidecars ignored in standalone report,
  plural VCF/status flags unavailable, full reporting callback missing explicit state.
  `/tmp/wave1-cli-red.log` records the actual red run.
- Relative VCF path reproduction: calling metadata retained an old working-directory-relative
  path while consensus_context supplied the absolute VCF. Standalone initially joined the
  wrong path to the summary folder. Red test `/tmp/wave1-vcfpath-red.log`; now prefer
  consensus_context.vcf_path, retaining summary-relative legacy behavior.
- Independent review found portable summaries retained analysis_completed after report
  rendering failed. RuntimeError and KeyboardInterrupt reproductions failed before the
  fix (`/tmp/wave1-status-review-red.log`). The decorator now synchronizes terminal
  provenance into existing readable summaries, including stale summaries on failed reruns.
  Invalid old JSON cannot mask original execution exceptions.
- Requested report ImportError now propagates as execution_failed with nonzero CLI outcome.
  Old test expecting graceful completed-without-report was updated to assert the requested
  artifact failure contract; mutation/scientific computation is unchanged.

## Implementation and decisions

CLI report loading moved to cli_report.py; public cli.report command/import retained.
map/run --mapping-timeout preserves JSON/CLI precedence. --vcf preserved and explicit
repeatable --allele-vcf LABEL PATH overrides singular; no option uses summary allele
VCF provenance. --run-status overrides summary sibling sidecar; report output destination
is never searched for status. Duplicate track labels and malformed sidecars fail visibly.

Pipeline rendering passes analysis_completed only after every analysis stage, with all
VCFs and actual allele metadata. Status sidecar stays running during render, then records
completed only after callback success. Terminal failure provenance is also persisted
inside summaries for standalone portability. Explicit status remains authoritative over
summary status; valid mutation evidence is retained with execution warnings.

## Verification

- `uv run --locked --all-extras pytest tests/unit/test_wave1_cli.py tests/unit/test_cli_settings.py --no-cov -q`: initial 39 passed.
- `uv run --locked --all-extras pytest tests/unit/test_wave1_cli.py tests/unit/test_cli_run.py tests/unit/test_runtime_settings.py --no-cov -q`: 112 passed at initial integration checkpoint.
- Final lead focused tests `tests/unit/test_run_status.py tests/unit/test_wave1_cli.py`: 32 passed.
- Configured mypy on cli.py, cli_report.py, pipeline.py, run_status.py and settings.py passed.
- Combined gates and final independent review are recorded in wave1-delivery.md.

## Real standalone CLI paths

With env_clair3 and vntyper executable directories before the default PATH,
`uv run --locked --all-extras muconespan report` generated both:

- `/tmp/wave1-standalone-embedded.html`: input is fixed real HiFi dupC summary,
  auto-discovers both absolute consensus VCF paths, loads detailed repeats, bundled
  reference and actual mapping BAM; embedded mode.
- `/tmp/wave1-standalone-sidecar.html`: same summary with explicit repeated
  `--allele-vcf allele_1 PATH --allele-vcf allele_2 PATH`; sidecar mode.

Both commands exited 0 and reported the artifact path. Full shell invocation/logs
are retained in the session and `/tmp/wave1-standalone-{embedded,sidecar}.log`.
Corresponding real session/browser tests validate loci, tracks and variant records.
