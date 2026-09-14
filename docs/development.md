# Development

The repository uses one shared development contract in
[`AGENTS.md`](https://github.com/berntpopp/MucOneSpan/blob/main/AGENTS.md).
Claude Code and Gemini CLI import that file through their small vendor entry
files; Codex reads it directly. Directory-level `AGENTS.md` files cover package,
test, and script details. Change shared policy in its source rather than copying
it into each assistant's entry file.

## Environment and commands

Use Python 3.10 or newer and uv. From the repository root:

```bash
make dev
make ci-check
```

`make dev` synchronizes `uv.lock` with all dependency groups and extras, including development, report,
and documentation dependencies. Normal development uses the lock without
updating dependency resolution. Individual Make targets install only their required
group (`quality`, `test`, or `security`) or extra (`docs` or `report`); developer
tooling is managed with uv dependency groups rather than a published `dev` extra.
The local `dev` group includes all tool groups and pre-commit.

When intentionally changing dependencies, edit
`pyproject.toml`, run `make lock`, review the lock diff, and run `make dev` again.

| Command | Purpose |
| --- | --- |
| `make quality` | Ruff, formatting, configured mypy, file size, and workflow syntax |
| `make test-fast` | Unit tests without coverage or external tools |
| `make test-unit` | Unit tests with the 80% coverage gate |
| `make ci-check` | Quality and unit tests with at least 80% coverage |
| `make test-int` | Tests marked as tool-dependent integration tests |
| `make docs-check` | Strict documentation build |
| `make security-check` | Audit all locked extras against published Python advisories |
| `make build-check` | Distribution build and package validation |
| `make docker-test` | BuildKit runtime checks without image export/load |
| `make docker-build docker-smoke` | Build and test an image through Docker |
| `make format` | Apply Ruff formatting |
| `make lint-fix` | Apply supported Ruff fixes; review the resulting diff |

Ruff covers `src/`, `tests/`, and `scripts/`. Mypy covers the package and scripts
using the settings in `pyproject.toml`, including required function annotations.
`make quality` also validates GitHub Actions syntax and expressions.
Coverage includes branch measurement;
the 80% gate applies to the aggregate unit test result, not each individual file.
Do not weaken checks or swallow failures to make a change pass.

For a focused test:

```bash
uv run --locked --all-extras pytest tests/unit/test_config.py --no-cov
```

Install the repository hooks with `make hooks`. Hooks provide local feedback;
the checked-in Makefile targets and CI remain the common verification contract.

## Architecture

MucOneSpan reconstructs the published method described by Vrbacka et al. and supports HiFi and ONT
amplicons. The current implementation uses minimap2. The historical bwa-mem
description is not the implementation contract; see
[differences from the published method](getting-started/deviations.md).

| Component | Responsibility |
| --- | --- |
| `cli.py` | Click commands and pipeline orchestration |
| `config.py`, `data/repeats/` | Repeat dictionary, mutation definitions, packaged data |
| `ladder.py`, `data/reference/` | Synthetic reference contigs and bundled ladder |
| `mapping.py` | Read alignment, BAM processing, alignment statistics |
| `alleles.py` | Repeat-count clusters and allele detection |
| `calling.py`, `vcf.py` | Allele read extraction/remapping, Clair3 calls, VCF processing |
| `consensus.py` | Allele consensus sequences using bcftools |
| `classify.py`, `classify_types.py`, `repeat_alignment.py` | Repeat segmentation, nomenclature, mutation interpretation |
| `report.py`, `templates/` | Structured results and optional HTML report |
| `tools.py` | External command execution, environments, errors, tool versions |
| `scripts/` | Reference generation, simulation, benchmark and maintenance helpers |

Package resources live under `src/muc_one_span/data/` and
`src/muc_one_span/templates/`; there is no top-level runtime `data/` directory.
Resource changes need distribution checks, since an editable installation can
hide packaging omissions.

Authored code, configuration, and templates must stay **below 650 physical
lines**, including comments and blank lines. The gate checks tracked and
untracked nonignored files with these suffixes: `.py`, `.sh`, `.yml`, `.yaml`,
`.toml`, `.css`, `.js`, `.ts`, `.html`, and `.j2`, plus Makefiles and Dockerfiles.
Generated/data files, lockfiles, and prose are outside the gate. Split growing
modules into cohesive responsibilities before reaching 650. Retain public
imports where compatibility matters and retain meaningful tests when splitting
test files. Dense formatting and blanket exclusions defeat the purpose.

## Scientific contracts

- MUC1 repeat units are normally 60 bases. Ladder contig names encode canonical
  variable repeats; the pre- and after-repeat blocks are additional sequence.
  The generator defaults to `contig_1` through `contig_150` with 500-base flanks.
- The ladder combines flanks, pre-repeats 1–5, canonical repeats, and after-repeats
  6–9. Dictionary nomenclature and source provenance must remain traceable.
- Dictionary mutation coordinates use 1-based positions; Python sequence slices
  use 0-based positions. Check insertion anchors and inclusive deletion endpoints
  explicitly when changing coordinate conversion or VCF handling.
- Alleles may have equal or nearby repeat counts, and PCR bias can yield very
  unequal coverage. Keep expected behavior explicit for homozygous, close-length,
  asymmetric, low-coverage, and long-repeat cases.
- HiFi defaults to minimap2 `map-hifi`; ONT defaults to `lr:hq`. Platform selection
  also reaches Clair3. Preserve explicit `--minimap2-preset` overrides.
- Changes to thresholds, confidence scores, mutation naming, or allele assignment
  require behavior-specific evidence. Maintenance refactors should preserve them.

See [core concepts](getting-started/concepts.md),
[repeat nomenclature](reference/nomenclature.md), and
[limitations](reference/limitations.md) for domain background. Historical
benchmark observations are snapshots of a stated version and dataset, not
guarantees about a future change.

## Test layers and external tools

Unit tests run without bioinformatics executables or generated sequencing data.
Mock command calls at the boundary and assert returned results, command arguments,
errors, and coordinate behavior with small deterministic fixtures.

Integration tests require the tools exercised by each test: minimap2 and samtools
for alignment, bcftools for consensus/VCF work, and Clair3 with a suitable model
for variant calling. MucOneUp and pbsim3 are generation prerequisites; they are
not prerequisites for every integration test once reads exist. The conda tool
environment is separate from the Python environment created by uv.

```bash
conda env create -f conda/environment.yml
conda activate muconespan-tools
make test-int
```

The optional MucOneUp end-to-end tests use existing reads in
`tests/data/generated/`. With Clair3 and its dependencies on `PATH`, run:

```bash
# Preserve Clair3's Python interpreter on PATH while pytest uses the project env.
.venv/bin/python -m pytest tests/integration --no-cov
```

Clair3 tests discover its adjacent `models/hifi` directory or accept an explicit
`CLAIR3_MODEL` directory. The homozygous 60/60 and asymmetric 25/140 allele cases
have strict expected failures for the [documented limitations](reference/limitations.md).
They keep the original expected counts and tolerance; an unexpected pass requires
reviewing the limitation rather than silently retaining an obsolete expectation.

Inspect skip reasons. A successful invocation in an environment missing tools,
models, or generated data does not demonstrate that those behaviors work.
Report the commands, pass/fail/skip counts, and missing prerequisites in a PR.
See [benchmarking](guides/benchmarking.md) for portable generation and batch runs.

## CI and container efficiency

The workflow design was compared with the sibling `hum-clinical-reporting`
repository. Both benefit from locked dependency layers, small build contexts,
lock-keyed caches, explicit runner versions, and separate routine checks from
expensive tool checks. This pipeline retains Debian for its Conda/TensorFlow
bioinformatics stack; the reporting application's Alpine runtime is specific to
its own dependency set.

- Each Make target selects its required uv group or extra. A fresh Python 3.10
  unit environment installs 16 distributions, compared with 76 when all tooling
  and documentation dependencies were installed in every job.
- All five supported Python versions run the unit suite for runtime changes.
  Python 3.14 alone measures coverage; the other four avoid redundant coverage
  instrumentation and HTML/XML generation.
- PR file filters skip unrelated integration, package, dependency-audit, and
  container work. The always-running CI Gate rejects failures and unexpected
  skips. Branch pushes and manual runs execute all applicable checks.
- The solved integration Conda environment is cached by its specification.
  Python caches depend on `uv.lock`, so a version-only edit does not discard
  third-party packages. The editable project's metadata cache separately tracks
  `version.py`, keeping the installed version accurate after a bump.
- The Docker tools layer depends only on the pinned base and Conda inputs. The
  locked application and report dependencies live in a separate small virtualenv;
  external tool calls use the Conda interpreter. Source edits preserve the large
  tool layer. The smoke test verifies both interpreters, report rendering, and
  real Clair3 inference.
- PR container builds restore cache only and run the smoke test in a BuildKit
  test stage, without exporting/loading a Docker image. This stage adds no
  dependencies to the runtime it tests. Trusted main builds publish that same
  runtime ancestor and export final layers without duplicate builder environments.

Before optimization, the measured GitHub container build step took 424 seconds:
197 seconds were cache export and 119 seconds were image export/loading. The PR
BuildKit test target avoids both export tasks. Cold builds still need
to obtain the bioinformatics tools and models; their size is an inherent cost,
so compare cold and cached builds separately when measuring future changes.

## Changes and review

Keep implementation plans and decisions in `.planning/`, with completed plans
in `.planning/archive/`. Include the behavior being preserved, module boundaries,
relevant checks, and any remaining limitations. Plans are working documents;
developer instructions and published docs are the maintained reference.

Before review, inspect the diff for unrelated changes, generated artifacts,
hardcoded local paths, and compatibility breaks. Run `make ci-check` and any
additional checks relevant to the change. Update the changelog for user-visible
changes. PR descriptions should explain the concrete result and include actual
validation evidence, including skips and follow-up work.
