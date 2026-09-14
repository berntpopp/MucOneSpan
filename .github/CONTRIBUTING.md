# Contributing to open-pacmuci

Open-pacmuci analyzes MUC1 VNTR amplicons from PacBio HiFi and ONT sequencing.
Shared repository conventions live in [AGENTS.md](../AGENTS.md); architecture,
scientific contracts, and detailed commands live in
[the developer guide](../docs/development.md).

## Set up

Install Python 3.10+ and [uv](https://docs.astral.sh/uv/), then:

```bash
git clone https://github.com/berntpopp/open-pacmuci.git
cd open-pacmuci
make dev
make hooks
make ci-check
```

`make dev` installs the locked environment with all extras. `make hooks` installs
pre-commit checks and a pre-push unit test gate using the same targets as CI.
When changing dependencies, update `pyproject.toml`, run `make lock`, review
`uv.lock`, and repeat `make dev`.

## Develop and validate

- Keep changes focused and preserve CLI/output compatibility during refactors.
- Use Ruff formatting, annotate changed Python APIs, and document public behavior.
- Keep authored code, configuration, and templates below 650 physical lines
  (maximum 649, including blanks and comments). Split by responsibility, retaining
  public imports where needed. `make file-size` enforces the scope.
- Add deterministic unit tests for changed behavior; mock external commands.
  Use integration tests for actual tool behavior and report skips honestly.
- Keep generated reads, tool indexes, benchmark results, and credentials untracked.
- Run `make ci-check` before committing or pushing. This runs lint, formatting,
  configured mypy, file-size and workflow checks, then unit tests with at least
  80% branch-aware coverage.

Run additional checks relevant to the change:

| Change | Check |
| --- | --- |
| Documentation or CLI documentation | `make docs-check` |
| Dependencies/security | `make security-check` |
| Packaging, bundled data, or templates | `make build-check` |
| Bioinformatics tool behavior | `make test-int` |
| All portable checks before broad maintenance/release work | `make check` |

Use `make test-fast` for rapid unit feedback. For external tools, create and
activate `conda/environment.yml`; Clair3 and its model need separate setup.
Simulation uses MucOneUp and pbsim3. See
[benchmarking](../docs/guides/benchmarking.md) for prerequisites and commands.
A green test command with skipped external tests does not validate those paths.

## Submit a pull request

Use a focused branch such as `feat/allele-detection`, `fix/vcf-coordinates`,
`docs/development`, or `chore/tooling`. Prefer Conventional Commit messages such
as `fix(mapping): preserve explicit alignment preset`.

Describe the concrete problem and resulting behavior, link related issues, and
include actual validation results and remaining limitations. Update docs for
workflow or user-facing changes and `CHANGELOG.md` under `[Unreleased]` when
relevant. Preserve unrelated work and review the final diff before submission.

Plans belong in `.planning/`, with completed plans in `.planning/archive/`.
Questions and bug reports belong in GitHub issues or discussions; search for
existing reports first.
