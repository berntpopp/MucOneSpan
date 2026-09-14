# Read-phasing provenance and configuration guide fixes

## Scope

- `read_phasing.phase_same_length_reads` now records the exact WhatsHap phase
  argument vector in `read_phasing.phase_command` before invoking the tool.
- Metadata returned without a phase invocation exposes `phase_command: null`, so
  consumers can distinguish an unattempted phase from an invoked command.
- The runtime configuration guide states that experimental read phasing is an
  explicit, default-disabled `calling.read_phase` opt-in, describes the actual
  reference layout and ladder fields, and distinguishes proven VCF absence from
  unavailable or ambiguous projection.
- The guide is linked from the MkDocs navigation.

## Test-first evidence

Before the implementation change, the focused regression failed with
`KeyError: 'phase_command'`:

```text
uv run --locked --all-extras pytest -q \
  tests/unit/test_read_phasing.py::test_collision_free_phase_retains_distinct_records_and_all_other_fields
1 failed
```

The same regression passed after the command was recorded.

## Final validation

```text
uv run --locked --all-extras pytest --no-cov -q tests/unit/test_read_phasing.py
31 passed in 0.04s

uv run --locked --all-extras ruff check \
  src/muc_one_span/read_phasing.py tests/unit/test_read_phasing.py
All checks passed!

uv run --locked --all-extras ruff format --check \
  src/muc_one_span/read_phasing.py tests/unit/test_read_phasing.py
2 files already formatted

uv run --locked --all-extras mypy src/muc_one_span/read_phasing.py
Success: no issues found in 1 source file

make docs-check
Documentation built in 0.29 seconds
```

The strict documentation build emitted only MkDocs Material's upstream MkDocs
2.0 notice and the repository's existing unnaved `includes/abbreviations.md`
informational message. No external phasing or sequencing jobs were run.
