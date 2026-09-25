# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- MucSim-Bench, a realistic simulated benchmark (`scripts/benchsim.py`,
  `muc_one_span.benchsim`): stratified designs for three profiles (HiFi
  uncalibrated), generation with MucOneUp >= 0.45.0 kept outside Git, engine runs
  that keep failed cases in the denominator, cluster-bootstrap stratified reports,
  a failure atlas and a pre-registered decision rule for a sealed `test` split.
- Realism report against public PRJEB92208 aggregate targets; the optional
  `bench` extra provides edlib for it.
- Every benchmark tunable is a validated setting (`benchsim.bench_config`).
  `--bench-config` overrides it, and the settings hash is recorded per case
  and report. Regular splits cap smear plus chimera at the real off-peak
  maximum; smear 0.5 moved to `stress`.
- Amplicon template sizing floors the minor-allele share at
  `amount.min_minor_share`. Floored cases record `amount_capped`, and every
  amplicon case records the floor used (`min_minor_share`).
- Event targets are clamped off the conserved head (units 1-5) and tail
  (units 6-9) and recorded as `target_clamped`. Too-short alleles are
  `design_invalid`.
- MucSim-Bench guide (`docs/benchmark.md`).
- MucSim-Bench reason atlas. `evaluate` keeps each sample's clinical reason list
  (`clinical.reasons`) and evaluator `reconstruction_flags`; `report` tabulates
  normalised reasons by profile and design stratum and splits INCONCLUSIVE
  cases into resolvable and expected (split or depth below the caller's
  per-allele gate), with every threshold in the `atlas` bench-config section.
  Cases without a recorded depth are counted as `depth_unknown`.
- `generate` reuse is checked against `bench_generation_sha256` (design,
  amount, profiles and structures settings only), so report-, realism-, run- or
  atlas-only changes no longer invalidate generated cases.
- MucSim-Bench benchmark sets (`sets` bench-config section, `design --set`):
  the realistic `standard` headline set (depths 500-2000 ONT amplicon, 200-1000
  HiFi, 30-100 genomic; calibrated error; no or calibrated PCR bias; typical
  smear and chimera), a `clean` control and a `stress` set with the former
  harsh mix. Design IDs carry the set name, sets of a split share their
  haplotypes, `report` and `realism` write one section per set (headline
  first), and the decision rule (now v3) applies to the headline set only.
  Realism results are labelled indicative. The per-split smear levels and the
  `design` depth, PCR, error and chimera levels moved into the sets.
- `profiles.simulator_threads` bounds the threads of MucOneUp's simulator tools
  per generated case (written into every profile variant), so `generate --jobs`
  sets a predictable core budget.
- Benchmark sets also fix concatemer and off-target levels (`clean` has no
  molecule artefact at all); the generation hash covers only a case's own set
  levels and ignores `simulator_threads`; the reason atlas treats the `stress`
  set as expected (`atlas.expected_inconclusive_sets`).
  `bench_config_sha256` stays the full provenance hash.
- Owner-approved absolute targets in the decision rule (`targets` bench-config
  section, task 12e): a candidate must, on top of beating the baseline, clear
  a fixed floor or ceiling per bench set on `pathogenic_rate`,
  `inconclusive_rate` and `false_positive_rate` (`clean` >= 0.90 / <= 0.10 /
  <= 0, `standard` >= 0.80 / <= 0.20 / <= 0; `stress` has none), judged on the
  point estimate or a Clopper-Pearson CI bound (`targets.basis`). `report`
  evaluates every configured target for the candidate, pooled per set and per
  profile, and writes a pass/fail table to `report.md` and `report.json`
  (`decision.targets`); adoption now needs the relative rule and every named
  set's targets to pass. The decision rule is now v4; a changed threshold,
  comparator, basis or set membership needs a new pre-registration.
- Benchmark final-review fixes. `first_evaluation.json` records the rule that
  unsealed `test`, and only that rule is accepted afterwards. `generate` also
  refuses to reuse a case whose MucOneUp version, base read profile, MucOneUp
  config or `--flank-fasta` hash differs or is not recorded. A targeted set with
  no cases in the output root is reported as "not present" (it still blocks
  adoption), and `report.json` holds target tables for every engine, including
  reports without a decision. `report.json` also records harness and caller
  provenance (`run` writes `caller.json`). Manifests merge rows by `design_id`.
  A caller crash is `execution_failed`. A realism `IndexError` is recorded for
  its case. Engines scored on different cases give a clean `report` error. The
  `bench` extra now requires `edlib>=1.3.9` on every Python version, and on
  3.14 edlib builds from source.
- Task C1 owner ruling (2026-09-25): the relative false-positive
  non-inferiority margin (`report.ni_margin`) is dropped from the decision
  rule -- at the planned `test` size (280 normals per profile) its Newcombe
  upper bound could never clear a meaningful margin, even with 0 observed
  false positives in both engines. The false-positive `PATHOGENIC` rate is now
  judged solely by Part 2's absolute `false_positive_rate` targets, and is
  reported per profile with a Clopper-Pearson interval for information only.
  A bench config still naming `report.ni_margin` is rejected. The decision
  rule is now v5; a fresh `test` pre-registration is required.

## [0.16.1] - 2026-09-24

### Fixed

- A frameshift indel that Clair3's pileup stage calls (allele fraction >= 0.5, depth >= 10)
  but that is missing from the applied calls now blocks a negative result
  (`alleles[*].stage_concordance`). On simulated older ONT reads with the R9 model,
  NEGATIVE on pathogenic samples falls from 12/39 to 1/39; no other benchmark or
  PRJEB92208 decision changes (#72).
- Clinical decision thresholds are configurable and taken from the recorded run
  configuration: `clinical_decision.max_ambiguous_bases` (10) and
  `clinical_decision.legacy_min_total_reads` (30) (#73).
- Read-length splitting uses the configured reference layout and the dictionary repeat
  length instead of fixed values (#74).

### Added

- Settings `calling.stage_discordance_min_af`, `calling.stage_discordance_min_depth`, the
  `clinical_decision` section, 14 `allele_selection` heuristics formerly hardcoded, and
  `read_phasing.min_haplotype_reads`. All defaults equal the previous behaviour.
- Output field `alleles[*].stage_concordance`.

### Changed

- The experimental read-phased haplotag split's minimum-reads-per-haplotype floor is now
  controlled by `read_phasing.min_haplotype_reads` (default `5`) instead of being coupled
  to the VCF-filtering `min_dp` argument, which was never applied to filtering itself.
  Library callers who previously relied on passing a non-default `min_dp` to change the
  haplotag split must set `read_phasing.min_haplotype_reads` instead; default behaviour is
  unchanged (#74).
- When an allele's Clair3 partition has no `pileup.vcf.gz` (`stage_concordance.status =
  "not_assessed"`), a NEGATIVE report now carries an extra detail, `Quality caveat: Allele N:
  caller-stage concordance not assessed (...); Clair3 pileup-stage frameshift calls were not
  compared with the applied calls.`. The decision itself does not change (#72).

## [0.16.0] - 2026-09-24

### Fixed

- Distinct-length candidates no longer resolve a heterozygous genotype to the
  reference allele. Remaining heterozygous, phased or conflicting records use the
  IUPAC candidate with `allele_genotype_status` and no independent haplotype
  credit. This fixes the PRJEB92208 MP4 dupC false negative (#53). The same
  selection applies to each haplotype of experimental read-backed phase.
- Haploid genotypes on length-partitioned alleles use the allele-specific AD
  fraction instead of FORMAT/AF. The 0.5/0.2 cut-offs are now configurable (#65).
- VCF concordance mirrors `bcftools consensus -H I` for heterozygous indels and
  is judged per event. The new status `heterozygous_genotype_unresolved` is
  never support (#64).
- Clinical decisions require explicit support, a frameshift, exact template
  identity, resolved localization and adequate carrier-allele depth for
  PATHOGENIC. NEGATIVE additionally requires resolved allele selection,
  consistent length and adequate per-allele primary depth (#55, #63).
- `--min-qual`, `calling.haploid_majority` and `calling.haploid_min_qual` are
  honoured and the applied threshold is recorded. `consensus.haploid_*` are
  deprecated no-ops (#66).
- A requested IGV report is preflighted before mapping (#67). Cluster and
  remapped BAMs no longer share `allele_reads.bam` (#68).

### Added

- Settings `calling.haploid_alt_fraction`, `calling.haploid_ref_fraction`,
  `allele_selection.secondary_mode_min_fraction` and
  `allele_selection.min_allele_primary_records`.
- Additive output fields:
  - `alleles[*]`: `allele_genotype_status`, `heterozygous_sites`,
    `variant_filter`, `selection_status`, `secondary_mode_fraction`,
    `depth_status`, `depth_basis`, `depth_threshold`, `length_status`;
  - `vcf_projection.unresolved_genotype_edits`.

### Changed

- Untemplated (novel) frameshifts, and events on alleles below the per-allele
  depth gate, are reported INCONCLUSIVE with reasons instead of PATHOGENIC.
  Summaries without `vcf_support` no longer count as supported.
- Re-rendering summaries written before 0.16.0 that lack `template_match` or
  `mutation_name` yields INCONCLUSIVE instead of PATHOGENIC.
- Direct API callers of `filter_vcf(haploid_majority=True)` without
  `haploid_min_qual` no longer get an implied `min(min_qual, 4.0)`; `min_qual`
  is applied. Pass `haploid_min_qual=4.0` to keep the previous default threshold.
- Deprecated `consensus.haploid_*` settings are reported with a logged warning,
  which the CLI shows, instead of a hidden `DeprecationWarning`.

## [0.15.1] - 2026-09-23

### Changed

- Bump developer and CI dependencies:
  - `ruff` from 0.16.7 to 0.16.8 (#61).
  - `playwright` from 1.62.0 to 1.63.0 (#61).
  - `docker/setup-buildx-action` from 4.3.0 to 4.4.1 (#60).
  - `docker/build-push-action` from 7.3.0 to 7.4.0 (#60).
  - `codecov/codecov-action` from 7.0.0 to 7.1.1 (#60).

## [0.15.0] - 2026-09-16

### Added

- Reproducible ENA clinical-data benchmark with validated downloads, separate clinical
  truth scoring, bounded execution, provenance-checked replay and resource records (#44).
- PRJEB92208 ONT cohort baseline and independently sourced HG002 Q100 sequence
  comparison; reported positive controls, unknown negatives and unresolved sample
  identity remain explicit. No PacBio validation or general clinical-performance claim.
- Containerized VNTRPipeline v1.0 diagnostic comparator evaluation on PRJEB92208
  clinical inputs under pinned dependencies and offline BiocManager runtime (#56).

### Fixed

- Recognize exact seven-C-tract duplications across compatible MUC1 repeat
  backgrounds, preserving the parent repeat and 61-base boundary. This corrects
  B-repeat dupC being misreported as an A insertion in the following repeat (#52).
- Guard consensus genotype selection when phase status is unphased, disconnected,
  or conflicting, selecting IUPAC candidate consensus ("I") and clearing independent
  haplotype evidence instead of forcing alternate GT1 (#53).
- Normalize allele candidate valley splitting by contig reference length and filter
  out sub-biological degradation contigs (< 10 units), preventing long full-span
  alleles from being lost to short fragment clusters (#54).
- Calibrate clinical decision calling to require frameshift, unambiguous localization,
  and confirmed VCF support for PATHOGENIC calls, appropriately classifying in-frame
  insertions, unresolved duplicate candidates, and unphased heterozygous reconstructions
  as INCONCLUSIVE (#55).

## [0.14.1] - 2026-09-15

### Fixed

- Bound streaming minimap2/samtools mapping with a configurable finite timeout,
  isolated process groups, concurrent stderr draining and partial BAM/index cleanup (#41).
- Open report IGV sessions on actual allele contigs and include every allele VCF,
  deduplicating shared files. Preserve `--vcf`; add repeatable `--allele-vcf LABEL PATH` (#39).
  Reject malformed VCF payloads from affected external igv-reports versions
  instead of silently displaying empty evidence tracks (1.13.0 verified).
- Carry execution provenance into reports, prevent known unfinished/failed runs
  from appearing negative, and propagate requested report-generation failures.
  Preserve autoescaping with Unicode missing-value dashes and display
  `exact_match_pct` consistently as the producer's 0–100 percentage (#40).

## [0.14.0] - 2026-09-15

### Added

- Add `DurableLedger` in `src/muc_one_span/durable_ledger.py` with multi-attempt flock locking, broken symlink self-healing, atomic append with fsync, and process-group SIGTERM/SIGKILL cleanup in `src/muc_one_span/tools.py`.
- Add 3-state clinical decision banner (`PATHOGENIC`, `NO_PATHOGENIC_VARIANT_DETECTED`, `INCONCLUSIVE`) with distinct SVG geometric icons, WCAG 2.1 AA/AAA contrast, and allele multiplicity caveats in `src/muc_one_span/report.py` and templates.
- Add offline browser test suite (`tests/browser/test_report_browser.py`) with Playwright verifying zero CDN requests, zero console errors, keyboard navigability, and visual layout.
- Add comprehensive 500-dataset simulation benchmark evaluation across DEV (300), VAL (100), and TEST (100) splits.

### Fixed

- Eliminate bimodal midpoint bias in length candidate clustering (`src/muc_one_span/length_candidates.py` and `src/muc_one_span/alleles.py`) by replacing integer division with nearest-peak assignment and equidistant exclusion, preventing systematic $+1$ repeat shift for $\Delta \ge 2$ alleles.
- Conform HGVS nomenclature to 20.05 specification (`src/muc_one_span/nomenclature.py`), hiding unanchored transcript coordinates in reports and establishing `repeat_{idx}:c.{edit}` as the invariant clinical coordinate.
- Align `delinsAT` nomenclature to net +1 bp edit naming (`55delinsAT`) in `KNOWN_VARIANTS`.
- Support single-heterozygous-site haplotagging in `src/muc_one_span/read_phasing.py`.
- Derive sequence identity status from VCF concordance in `src/muc_one_span/calling.py` and eliminate shallow-copy aliasing on `allele_2`.

## [0.13.0] - 2026-09-15

### Added

- Add support for Oxford Nanopore Adaptive Sampling / genomic read simulation via
  NanoSim (`muconeup reads ont`) and PacBio HiFi genomic simulation via PBSIM3/CCS
  (`muconeup reads pacbio`), overcoming PCR length-dependent dropout in long alleles.
- Add 500-dataset stratified simulation catalog with 250 biological designs spanning
  7 clinical challenge categories across development (300 datasets), validation
  (100 datasets), and test (100 datasets) splits with dual cryptographic ledgers
  (`scripts/build_500_design.py`, `scripts/run_500_experiment.py`,
  `examples/experiment_500_designs.json`).
- Add official HGVS nomenclature engine (`src/muc_one_span/nomenclature.py`)
  conforming to HGVS recommendations (https://hgvs-nomenclature.org/stable/):
  supports 3'-most cDNA normalization (`NM_001204286.1`), 5'-roll ambiguity intervals,
  repeat-unit shorthand (`53C[7]>53C[8]`, `59dupC`), and clinical tiers (Tier A/B/C).
- Add offline, zero-CDN interactive IGV alignment browser into HTML reports via
  `--report-igv [embedded|sidecar|off]`, with vendored gzipped `igv.js` assets and
  SHA-256 integrity validation (`src/muc_one_span/report_igv.py`, `report_assets.py`).
- Add non-canonical/IUPAC nucleotide detection and validation warnings in consensus
  construction (`src/muc_one_span/consensus.py`).

### Fixed
- Eliminate IUPAC ambiguity code injection (`S`, `M`, `R`, `Y`) by forcing haploid
  consensus replay (`bcftools consensus -H 1`) on single-allele candidate calling
  in `src/muc_one_span/calling.py`, eliminating the primary root cause of sequence
  discordance on length-exact alleles.
- Fix negative coordinate indexing underflow in `igv_reports` by creating 1-based
  VNTR-spanning intervals on allele contigs.
- Fix standalone embedded IGV rendering in HTML reports by cleanly extracting
  container markup and initializing the offline variant track selector without
  unreplaced template placeholders.

## [0.12.0] - 2026-09-15

### Added

- Add pairwise read-dominance candidate length inference (`read_dominance.py`)
  with symmetric missing-score handling for candidates >= 6 repeats apart and
  primary alignment cluster filtering in allele detection (`alleles.py`).
- Add metric-aware ONT contig refinement (`metric="auto"`) minimizing mean indel
  length to eliminate systematic +1 repeat shift in long Oxford Nanopore contigs.
- Support multiallelic comma-separated `AF` fields in VCF parsing (`vcf.py`).
- Update reference ladder generation to use proximal left-flank sequence for
  exact anchor consensus trimming (`ladder.py`).
- Add 200 simulation dataset evaluation suite across PacBio HiFi and ONT amplicons
  with sealed and public ledgers (`scripts/run_200_experiment.py`, `scripts/evaluate_200_experiment.py`).

### Changed

- Complete diploid sequence reconstruction improved by +28.6 percentage points on
  PacBio HiFi (+20 net paired reconstructions) and +17.1 percentage points on ONT
  (+12 net paired reconstructions) across 140 development datasets, with 0 paired
  regressions on previously correct baselines.
- Maintained 100% normal control negative specificity (0 false alarms).

## [0.11.0] - 2026-09-14

### Added

- Add strict JSON simulation designs for reproducible MucOneUp HiFi/ONT experiments,
  preserving the 29 historical cases and providing a larger ungenerated example.
  Record requested templates separately from usable reads, exact commands, hashes,
  truth consistency and all failed cases; reject overwrites and changed artifacts.
- Restrict source-distribution contents to maintained files and verify both source
  and wheel archives exclude generated reads, results, indexes and private data.
- Document measured improvements, unsuccessful experiments and remaining failures.
  Fresh final scientific validation was not performed for this scoped release;
  general diploid reconstruction accuracy and two known length failures remain open.

- Add immutable validated runtime settings and global `--config` JSON input;
  explicit CLI options override file values. Record effective settings and
  configuration, read, reference and dictionary hashes in `run_configuration.json`.
- Centralize classification, confidence, length-selection, reference-layout,
  consensus and optional phasing tunables; reject invalid stage overrides.
- Correct boundary-anchor coordinates and sliced lengths, matching the actual
  ladder flank prefix; preserve cached development FASTAs and repair flanking
  indel trimming in regression fixtures.
- Standalone consensus now uses the configured or bundled repeat dictionary for
  anchor-aware trimming, matching the full pipeline; `--repeats-db` remains an
  additive override. Failed reruns invalidate stale configuration provenance.

- Add experimental opt-in WhatsHap phasing for unresolved same-length
  candidates, disabled by default after a development false-call regression.
- Version strict batch evaluation outputs and update catalog display compatibility.
- Report trimming provenance and exact VCF-projection failure reasons; prevent
  flank-only variants from establishing independent identical VNTR reconstructions.
- Add `length_selection_evidence` for excluded candidate-reference records and
  `phasing_selected_records_by_phase_set` for the experimental phaser's internal
  read selection, with explicit denominator limits.
- Distinguish literal missing outputs from `independent_missing_alleles` and
  `unproven_duplicate_alleles` in evaluation; unknown evidence states cannot count
  as completed negatives.
- Clear stale phase/VCF evidence on unresolved aliases during subcommand reruns.
- Label unavailable VCF projection separately from unsupported mutations and
  preserve dictionary-fit confidence when support cannot be assessed.

### Changed

- Exact bit-vector Levenshtein scoring preserves literal distances and traceback;
  classification benchmarks distinguish local from pipeline speedups.
- Signed net insertion/deletion length now determines downstream frameshift.
  Classification reports consumed spans, unresolved regions and heuristic-score
  semantics; uncertain sequence is not proof of complete reconstruction.
- VCF query and malformed-record failures propagate. Mutation support now requires
  exact sequence concordance through a replay-verified reference projection;
  position-only support callers receive an unavailable status. Confidence weights
  consequently change; they remain uncalibrated.
- Consensus explicitly selects sample/genotype alleles. Same-length calls preserve
  common phase-block evidence and expose unresolved phase; absence of heterozygous
  calls no longer establishes sequence homozygosity. Consensus context records the
  actual trimming coordinates, and final allele metadata is persisted.
- Allele fit counts are explicitly alignment records; primary record counts and
  actual selected-reference length are separate fields. Flat indel profiles no
  longer cause an arbitrary split.
- Strict offline evaluation includes failed/missing samples, tied one-to-one
  assignments, missing/extra alleles and exact event annotation rather than
  substring detection. Generated-data benchmarking metrics are versioned.


## [0.10.0] - 2026-09-14

### Changed
- Name the tool and repository MucOneSpan, following the companion MucOneUp conventions.
- Use `muconespan` for the CLI and `muc_one_span` for the Python distribution and imports
  (normalized distribution name: `muc-one-span`).
- Update package resources, containers, workflows, documentation, citations, and
  repository links consistently; preserve the scientific pipeline behavior.

## [0.9.0] - 2026-09-14

### Fixed
- Docker build permission failure during unnecessary builder cleanup.
- Vulnerable Python dependencies in the runtime, test, and documentation lockfile.
- Refresh editable package metadata when the dynamically loaded version changes.

### Changed
- Share concise agent instructions across Codex, Claude Code, and Gemini CLI.
- Enforce fewer than 650 physical lines in authored code, configuration, and templates.
- Use locked dependencies and matching local/CI lint, format, type, workflow,
  coverage, packaging, documentation, and dependency audit checks.
- Test Python 3.10–3.14 and build/test containers on pull requests before publishing.
- Reduce CI setup with targeted uv groups, one coverage run, cached tool environments,
  and checks selected by changed files; keep stable gates for every pull request.
- Preserve Docker tool layers across source changes, use locked application dependencies,
  and avoid uploading large build caches from pull requests.
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
- Optional `[report]` extra for Jinja2 dependency (`pip install muc_one_span[report]`)
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
- Initial release of MucOneSpan pipeline
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

[Unreleased]: https://github.com/berntpopp/MucOneSpan/compare/v0.16.1...HEAD
[0.16.1]: https://github.com/berntpopp/MucOneSpan/compare/v0.16.0...v0.16.1
[0.16.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.15.1...v0.16.0
[0.15.1]: https://github.com/berntpopp/MucOneSpan/compare/v0.15.0...v0.15.1
[0.15.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.14.1...v0.15.0
[0.14.1]: https://github.com/berntpopp/MucOneSpan/compare/v0.14.0...v0.14.1
[0.14.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.13.0...v0.14.0
[0.13.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.12.0...v0.13.0
[0.12.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.11.0...v0.12.0
[0.11.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/berntpopp/MucOneSpan/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/berntpopp/MucOneSpan/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.0.0...v0.1.0
[0.0.0]: https://github.com/berntpopp/MucOneSpan/releases/tag/v0.0.0
