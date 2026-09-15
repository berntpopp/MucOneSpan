# Wave 1 — Issue #41 mapping lifecycle evidence

## Scope and interface

Read the live issue body/comments, root/package/test instructions, developer guide,
contributing guide, complete mapping/tool implementations and relevant tests.
Issue #41 has no comments. The proposed `muc_one_span.process` API does not exist.

Worker A owns mapping/tool lifecycle helpers and focused tests; the lead owns
`RunSettings.mapping_timeout`, CLI/configuration/pipeline plumbing and docs.
Agreed interface: `map_reads(..., *, timeout=DEFAULT_SETTINGS.run.mapping_timeout)`;
the setting defaults to 3600 seconds. Existing positional mapping arguments,
HiFi `map-hifi`, ONT `lr:hq`, and explicit preset overrides are retained.

## Reproduced defects and test-first evidence

Before implementation, `uv run --locked --all-extras pytest
 tests/unit/test_mapping_lifecycle.py --no-cov -q` produced **7 failures**:

- Both RuntimeError and KeyboardInterrupt left stale/partial BAM and indexes.
- Four invalid timeout cases could not be rejected because timeout was absent.
- Mocked successful launches lacked both `start_new_session=True` and cleaned PATH.

The initial command with just artifact/timeout regressions had six failures;
adding the environment/session regression produced seven. Source inspection
confirmed unbounded `communicate`, stderr-thread `join`, and producer `wait`,
with only the direct producer killed when the consumer executable was missing.

Additional red-green cycles found and fixed boolean timeout acceptance, stopping
stderr drains before consuming final diagnostics, input/output collision data
loss in the first implementation, and lost partial supervisor JSON on interruption.
Independent review reproduced the initial implementation leaving orphan zombies;
the final Linux supervisor and strict PID-disappearance tests address that gap.

## Final design

- `tools.run_tool_pipeline` is the external-command boundary and applies existing
  virtualenv PATH cleanup. Streaming goes directly between pipe file descriptors;
  final command output goes to a supplied file or `/dev/null`.
- `tool_supervisor.py` starts a separate Python supervisor session on Linux and
  enables `PR_SET_CHILD_SUBREAPER` only there. It adopts and reaps orphaned
  descendants; it does not alter the application's subreaper state or wait for
  unrelated children. The caller receives group IDs and a structured result.
- `tool_pipeline.py` starts each command in an isolated session, monitors both
  return codes, and drains both stderr streams concurrently through stoppable
  nonblocking readers. It retains the final MiB per tool with replacement decoding.
- One absolute monotonic budget covers waits and joins. Up to 0.5 seconds (10%
  for short budgets) is reserved inside it for TERM/KILL cleanup. The caller
  reserves two small cutoffs (each at most 0.05 seconds, 1% for short budgets)
  for retry signaling and final supervisor reaping; retry cannot consume the
  final reap reserve. Interrupted protocol reads retain buffers and parsed state.
- Mapping conversion, streaming alignment/sort, and indexing consume the remaining
  overall mapping budget. Conversion within map_reads also streams to FASTQ.
- Failure/interruption removes output BAM and `.bam.bai`, `.bam.csi`, `.bai`, `.csi`
  indexes. Input/output collisions are rejected before unlinking the source.
- Primary errors, both tool stderr tails, and cleanup errors remain visible.

## Verification performed by Worker A

| Command/scope | Actual result |
| --- | --- |
| Mapping, lifecycle, worker, supervisor, tools units and `TestMapSubcommand` | 78 passed, 0 skipped (final complete focused command) |
| `tests/integration/test_tool_pipeline.py --no-cov -q` | 10 passed, 0 skipped; final rerun alongside 15 supervisor units: 25 passed in 6.60 seconds |
| Initial actual-tool `tests/integration/test_mapping.py --no-cov -q -rs` on default PATH | 1 passed, 2 skipped: minimap2 not installed on that PATH |
| Focused Ruff and formatting checks | Passed |
| Focused mypy over mapping, tools, tool_pipeline, tool_supervisor | Passed: 4 source files |
| Physical line counts | Every owned source/test file below 650; largest 270 lines at verification |

The ten real-process scenarios are success, heavy stderr, producer failure,
consumer failure, missing second executable, producer hang, consumer hang,
interruption, a descendant holding pipes after its parent exits, and a descendant
with closed inherited pipes after its parent fails. Child-created readiness
files provide handshakes, an outer eight-second watchdog bounds regressions, and
all recorded direct-tool/descendant PIDs must cease to exist (zombies fail).
The output check confirms the streamed payload is preserved exactly.

Unit tests cover command presets and argument lists, environment cleanup,
process isolation, early peer-failure detection, bounded escalation/waits/joins,
launch exceptions, partial BAM/index removal, shared stage budgets, invalid
budgets, stderr tails/errors, supervisor protocol splits, interruption retries,
subreaper adoption/reaping, primary-plus-cleanup diagnostics, and explicit
incomplete-cleanup errors.

The lead owns full-suite gates and baseline/fixed scientific-panel comparisons,
including actual minimap2/samtools/model execution in the external tool environment.

## Explicit limits

Linux provides the strong descendant-adoption/reaping guarantee. Other POSIX
platforms keep the compatible group-termination implementation and emit a warning
that orphan reaping belongs to the OS adopter. No non-Linux execution was tested.
An uninterruptible kernel task can defeat even SIGKILL within a finite budget;
cleanup incompleteness is reported rather than hidden behind an unbounded wait.
No sensitivity, classification, mutation, or clinical-tier policy changed.
