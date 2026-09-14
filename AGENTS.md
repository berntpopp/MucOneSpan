# Repository guidance

This is the shared instruction source for Codex, Claude Code, and Gemini CLI.
Read any `AGENTS.md` inside directories you change; scoped instructions add to
this file. Keep vendor entry files as pointers rather than duplicate guidance.

## Start here

- Python 3.10+ package for MUC1 VNTR analysis of PacBio HiFi and ONT amplicons.
- Read [developer documentation](docs/development.md) for architecture and checks,
  and [contributing](.github/CONTRIBUTING.md) for the contribution workflow.
- Inspect the working tree and relevant code/tests before editing. Preserve
  unrelated changes. Keep the requested behavior and compatibility explicit.
- Put substantial implementation plans in `.planning/`; move completed plans
  into `.planning/archive/`. Record assumptions, decisions, and validation.

## Development contract

```bash
make dev          # Reproduce the locked environment, including all extras
make quality      # Ruff, formatting, configured mypy, file size
make test-fast    # Unit tests without external tools or coverage
make ci-check     # Quality + unit tests with >=80% branch-aware coverage
```

- Run `make ci-check` before committing or pushing code. Run the applicable
  `make test-int`, `make docs-check`, `make security-check`, and `make build-check`
  for changes to integration behavior, docs, dependencies/security, and packaging.
- Use the Makefile targets and `uv run --locked --all-extras` for ad hoc commands.
  Update `pyproject.toml` and `uv.lock` together when dependencies change.
- Authored code, configuration, and templates must contain **fewer than 650
  physical lines** (maximum 649, including blanks/comments). The file-size gate
  scans tracked and untracked nonignored files; generated/data files, lockfiles,
  and prose are outside the gate. Split by responsibility; preserve imports and
  behavior. Do not compress lines or add exclusions to evade the limit.
- Match Ruff configuration, annotate changed Python APIs, and document public
  behavior. Keep configured type checks passing without suppressing failures.
- Add focused tests for behavior changes. Unit tests must be deterministic and
  mock external tools; report integration skips and missing prerequisites.
- Preserve CLI flags, output schemas, bundled resources, and scientific semantics
  during refactors. Keep user-facing version values sourced from `version.py`.
- Keep subprocess execution in the existing tool abstraction, with argument
  lists and visible failures. Do not add shell interpolation or machine-specific
  paths. Never commit generated reads, tool indexes, results, or credentials.

## Completion

Review the final diff and report what changed, actual commands/results, and any
untested behavior. A skipped test is not validation of its behavior. Keep PRs
focused; update documentation and the changelog when relevant. Do not merge or
publish releases unless the user has authorized that action.
