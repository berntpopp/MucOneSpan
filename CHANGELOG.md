# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- Docker build permission failure during unnecessary builder cleanup.
- Vulnerable Python dependencies in the runtime, test, and documentation lockfile.

### Changed
- Share concise agent instructions across Codex, Claude Code, and Gemini CLI.
- Enforce fewer than 650 physical lines in authored code, configuration, and templates.
- Use locked dependencies and matching local/CI lint, format, type, workflow,
  coverage, packaging, documentation, and dependency audit checks.
- Test Python 3.10–3.14 and build/test containers on pull requests before publishing.
- Replace placeholder integration tests with real tool and pipeline assertions;
  expose documented scientific limitations as strict expected failures.


## [0.8.0] - 2026-04-07

### Fixed
- ONT allele length off-by-one: derive canonical repeat count from AS-refined contig instead of noisy cluster center, with ±1 guard to avoid regression on long alleles (#18)

### Added
- CITATION.cff for machine-readable citation metadata
- ONT test data generation in `scripts/generate_testdata.py`
- BibTeX citation block in README

### Changed
- Upgrade codecov/codecov-action v5 to v6 (fixes Node.js 20 deprecation)
- Update project description to mention ONT support

## [0.7.0] - 2026-04-07

### Added
- `--platform` CLI option on `run`, `map`, and `call` subcommands to select sequencing platform (`hifi` or `ont`)
- `--minimap2-preset` CLI option on `run`, `map`, and `call` subcommands to override minimap2 alignment preset
- Auto-selection of minimap2 preset from platform (`hifi` -> `map-hifi`, `ont` -> `lr:hq`)
- ONT support: Clair3 `--platform=ont` and minimap2 `lr:hq` preset threaded through full pipeline
- `minimap2` added to `check_tools` for `call` subcommand (upfront dependency check)

## [0.6.0] - 2026-04-07

### Added
- Self-contained HTML report with modern UI/UX, hover tooltips, dark/light mode, print stylesheet
- Optional `[report]` extra for Jinja2 dependency (`pip install open-pacmuci[report]`)
- `--report` flag on `run` subcommand and standalone `report` subcommand
- Rich hover tooltips throughout report (repeat blocks, metrics, section headers)
- Parallel per-allele variant calling via ThreadPoolExecutor
- Streaming `run_tool_iter()` for memory-efficient SAM parsing
- Pipeline benchmark script (`scripts/benchmark.py`)

### Changed
- Streaming SAM parsing in `refine_peak_contig` and `_split_cluster_by_indel`
- Split thread budget across parallel alleles to avoid CPU oversubscription

## [0.5.0] - 2026-04-07

### Added
- TypedDict types for classification and allele results (documentation types)
- Dedicated `vcf.py` module for VCF parsing and filtering
- 18 new unit tests (mapping pipeline, CLI run, classify helpers, indel valley, VCF variants)
- Pre-push checklist and common pitfalls in CLAUDE.md

### Changed
- Decompose `classify_sequence()` into 3 focused helpers
- Raise CI coverage threshold from 70% to 80%
- Strict `make ci-check` (no swallowed mypy errors, coverage gate)

### Fixed
- Docker build: remove nonexistent `data/` COPY
- Dynamic version test (no more hardcoded version strings)
- Expand `__all__` in `calling.py` to include all public functions

## [0.4.0] - 2026-04-07

### Added
- Structured logging in all pipeline modules with --verbose/-v and --quiet/-q CLI flags
- Community health files: CONTRIBUTING.md, CHANGELOG.md, SECURITY.md, CODE_OF_CONDUCT.md, CODEOWNERS, issue/PR templates
- Multi-stage Docker build with micromamba slim base for smaller, faster, secure images
- BuildKit layer caching in Docker CI workflow
- Singularity/Apptainer definition for HPC environments
- GitHub Release automation on version tags
- Tool version recording in summary.json for reproducibility
- Conda version pinning (minimap2=2.28, samtools=1.21, bcftools=1.21, htslib=1.21)
- Validation script for repeat classification catalog
- Known limitations documentation page

### Fixed
- Replace assert statements with descriptive RuntimeError in classify.py and mapping.py
- Fix potential deadlock in mapping pipeline with concurrent stderr draining
- Use ERROR level for --quiet flag (previously WARNING, same as default)

## [0.3.0] - 2026-04-07

### Added
- Soft QUAL scoring with continuous confidence and boundary penalty

### Fixed
- Address review feedback on soft QUAL scoring
- Address Copilot review comments on soft QUAL scoring PR

## [0.2.0] - 2026-04-07

### Added
- MkDocs Material documentation site with GitHub Pages deploy
- Comprehensive README with deviations, usage, and validation

## [0.1.2] - 2026-04-07

### Added
- Indel-valley allele splitting for close allele pairs

## [0.1.1] - 2026-04-06

### Fixed
- Strip virtualenv PATH from subprocess calls for external tools
- Handle empty Clair3 VCFs and remove INFO/DP filter

## [0.1.0] - 2026-04-06

### Added
- Initial release of open-pacmuci pipeline
- Reference ladder generation (20-150 repeat units)
- Read mapping with minimap2
- Allele length detection with peak finding
- Variant calling with Clair3
- Consensus building with bcftools
- Repeat unit classification with Vrbacka nomenclature
- Full pipeline CLI with individual subcommands
- Docker image with conda-based tool stack
- Unit and integration test suite

## [0.0.0] - 2026-04-06

### Added
- Project scaffolding with uv, ruff, mypy, pytest, CI
- Initial pipeline implementation

[Unreleased]: https://github.com/berntpopp/open-pacmuci/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/berntpopp/open-pacmuci/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/berntpopp/open-pacmuci/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.0.0...v0.1.0
[0.0.0]: https://github.com/berntpopp/open-pacmuci/releases/tag/v0.0.0
