# Test guidance

- Put deterministic Python tests in `unit/`; use `tmp_path` for outputs and mock
  external commands where the module under test looks them up.
- Put real tool tests in `integration/`, with `integration` or `e2e` markers and
  explicit availability checks for each external tool and required input/model.
  Missing optional fixtures should produce descriptive skips, not false passes.
- Assert observable scientific/CLI behavior, including errors and boundary cases.
  Avoid pass-only placeholders, hardcoded package versions, or assertions that
  merely repeat the implementation.
- Keep reusable fixtures in the narrowest appropriate `conftest.py`; preserve
  collection and parametrization when splitting files to remain below 650 lines.
- Run focused tests while developing, then `make ci-check`. Use `make test-int`
  for tool-dependent changes and report executed tests separately from skips.
- Generated MucOneUp datasets and outputs remain untracked. Generation and
  benchmark instructions live in `docs/guides/benchmarking.md`.
