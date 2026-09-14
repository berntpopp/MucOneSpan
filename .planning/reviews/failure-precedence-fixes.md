# Evaluator failure precedence repair — 2026-09-14

## Reproduction

`load_observation` parsed `run_status.json` before consulting external execution
records. A truncated JSON sidecar, unsupported schema, unknown status or nonobject
JSON caused `invalid_artifacts`, even when an external record established exit1,
timeout or execution_failed. An unreadable sidecar had the same problem. Valid
running-sidecar failure precedence was already covered and retained.

Thirteen new regression cases were observed failing before implementation:
four invalid-sidecar types crossed with three known external failure records,
plus an unreadable sidecar. Existing successful summary/consensus files were
provided in the malformed-content cases to verify they cannot restore stale
predictions. Two nonfailure controls retained the expected invalid_artifacts
status.

## Change

Compute the existing external-failure predicate before reading the sidecar. In
the sidecar error handler only, that established failure now yields
`execution_failed`, the original external error and no predictions. An
`invalid_run_status_sidecar:<detail>` warning retains the damaged-sidecar evidence.
If no external failure is known, the sidecar error remains `invalid_artifacts`.

Valid sidecar handling is unchanged. In particular, valid typed
`insufficient_evidence` plus the CLI's exit1 remains insufficient_evidence, rather
than being flattened into a tool failure. The change does not relax artifact
schemas, rewrite scoring, or alter the priority of other valid sidecar statuses.

Only `src/muc_one_span/evaluation/artifacts.py` and
`tests/unit/test_evaluation_artifacts.py` were edited in this bounded repair.

## Verification

- Red: `uv run --locked --all-extras pytest tests/unit/test_evaluation_artifacts.py --no-cov -q`
  —13 failures /32 passes before the fix.
- Green: the same artifact suite —45 passed.
- Regression: `uv run --locked --all-extras pytest tests/unit/test_evaluation_artifacts.py tests/unit/test_evaluation_cli.py tests/unit/test_evaluation_scoring.py --no-cov -q`
  —66 passed in0.07s after final formatting.
- Ruff check and format check: pass for both edited files.
- Configured mypy: pass for artifacts.py.
- Physical lines: artifacts.py228; artifact tests270, both below649.

Source SHA256:
`7b95e31c24c64c1927c92fd0b9b2fb11fa6409beb10a362c49aba823700b47a0`.
Test SHA256:
`f96a2b19fb0bc2586de2ceb0709d798dcfb8627228526cff0861a6e328ab4a52`.

No external tools or heavy datasets were needed. Repository-wide CI remains
coordinated by root after concurrent changes settle.
