# Agent development modernization

## Objective

Maintain one accurate development contract for Claude Code, Codex, and Gemini
CLI; enforce small authored files and reproducible checks; prepare a tested draft
PR while preserving scientific behavior and public interfaces.

## Design decisions

- `AGENTS.md` is the shared source. Minimal `CLAUDE.md` and `GEMINI.md` import it.
  Package, test, and script instructions are scoped to their directories.
- Put architecture, scientific invariants, setup, and test details in maintained
  developer docs. Keep historical benchmark snapshots in `.planning/` and avoid
  presenting them as fresh evidence.
- Use `uv.lock` with all dependency groups and extras for local checks, hooks, and CI. Dependency
  changes must update the project manifest and lock together.
- Enforce fewer than 650 physical lines for authored code, configuration, and
  templates, including blank lines/comments. Scan tracked plus untracked
  nonignored files. Exclude generated/data artifacts, lockfiles, and prose.
- Split modules and tests along responsibilities without changing CLI behavior,
  import compatibility, outputs, scientific thresholds, or test coverage intent.
- Reuse Makefile targets across local development, hooks, and CI. Include scripts
  in lint/format/type checks. Never suppress failed commands.

## Implementation sequence

1. Audit current instructions, checks, file sizes, tests, and package resources.
2. Establish locked Make targets, file-size enforcement, hooks, and aligned CI.
3. Split oversized files and fix newly surfaced maintained-script quality issues.
4. Add shared/scoped guidance, minimal vendor pointers, developer and contribution
   docs, and a PR template recording actual validation and skips.
5. Run portable static/unit/docs/security/distribution checks and available real
   integration tests. Review the final diff and document missing prerequisites.
6. Commit the reviewed change and open the user-authorized draft PR.
7. Per the follow-up instruction, bump the minor version to 0.9.0 and merge the
   PR only after all applicable tests and GitHub Actions pass.
8. Per the performance follow-up, compare `../hum-clinical-reporting`, optimize
   dependency groups and container layers/cache, measure cold/incremental behavior,
   and require the optimized head to pass before merging all outstanding PRs.

## Validation contract

- `make dev`: reproduce the locked environment with all dependency groups and extras.
- `make ci-check`: lint, format, configured mypy, file-size/workflow checks, and
  unit tests with at least 80% branch-aware coverage.
- `make docs-check`: strict documentation build.
- `make security-check`: audit the locked dependency set.
- `make build-check`: build/validate distributions and bundled resources.
- `make test-int`: execute available tool integration tests; record skips with
  their tool/data/model prerequisites.
- Review hooks/workflows and compare final collection and coverage against the
  baseline. Record actual results in the final PR, not inferred success here.

## Audit findings addressed by the design

The old CLAUDE file duplicated policy and retained stale bwa-mem architecture,
old ladder defaults, local home paths, and an unsupported claim of strict mypy.
The published benchmark guide described only 10 samples despite 26 HiFi and 3
ONT definitions. The previous PR checklist described `make check` inconsistently.
The new guidance derives commands from the maintained Makefile and separates
unit checks, external integration evidence, and historical scientific results.
