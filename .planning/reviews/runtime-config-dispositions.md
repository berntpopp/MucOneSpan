# Configuration audit dispositions before source freeze

The read-only internal audit reproduced six issues. This does not substitute for
the separately requested actual Claude Fable5.1 review.

| Finding | Fix and evidence |
| --- | --- |
| R1 prefix flank mismatch | Consensus anchors now use the same selected flank prefix as ladder construction. Real cached137 consensus FASTAs unchanged, and a synthetic flanking insertion now trims correctly. |
| R2 standalone dictionary ignored | Load configured or bundled dictionary before tools, forward it to consensus; additive --repeats-db works. This intentionally corrects prior standalone fixed-only trimming. |
| R3 CLI override validation bypass | Stage overrides validated through the same immutable settings contracts as JSON. threads0 and inverted ladder range regressions fail before fix and pass after. |
| R4 stale configuration on failed rerun | Invalidate prior run_configuration at invocation start, before file/layout validation. New run_status failure remains authoritative. |
| R5 tolerance coordinate off by anchor component | Preserve helper start-coordinate API, translate boundary expectation correctly, subtract actual sliced component lengths. Zero tolerance accepts expected boundary; oversized requested anchors use available components. |
| R6 missing actual phaser command metadata | Persist exact phase_command list before run_tool; null means not attempted. Mock invocation argv equivalence regression passes. |

Additional fixes: standalone configured defaults preserve programmatic default_map
when no file is supplied; ladder range alone does not require a custom full-run
reference; run provenance hashes the actual dictionary as well as reads/reference.
API defaults use DEFAULT_SETTINGS and one shared platform preset map. Per-tool
thread settings are not a strict process-tree CPU cap, as documented.

Root six-regression run:6failed/8passed before,14passed after; broader CLI tests
pass. Initial fullCI caught a mocked reference path without a file in the report
fallback test; write a valid FASTA fixture so hashing can execute, retaining the
same warning and exit assertions. Final make ci-check:626passed,93.59% coverage;
Ruff, format, mypy and649-line gate passed. make test-int:47passed,8deselected.
Strict docs build passed. Generated end-to-end/build checks were still finishing
at this note. Source edits remain uncommitted and all reserved seeds unused.

See consensus-anchor-config-fixes.md, read-phase-config-fixes.md and
api-defaults-evidence.md in the parent planning directory for focused evidence.
