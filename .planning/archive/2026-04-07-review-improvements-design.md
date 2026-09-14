# Review Improvements Design

**Date:** 2026-04-07
**Source:** `.planning/REVIEW_2026-04-07.md` (scored 8.5/10)
**Goals:** Maximize review score (target 9.5+) and achieve production/clinical readiness
**Approach:** Quick Wins First (Approach B) -- 4 phases, impact-first ordering

---

## Phase 1: Quick Wins -- Error Handling + Community Health

**Paths:** 2 (Production Error Handling) + 3 (Community & Project Health)
**Effort:** ~2 days
**Score impact:** Error Handling 8->9, Documentation 8->9, Git 8->9

### 1A. Error Handling

1. **Replace `assert` with guarded `raise`** (2 locations):
   - `classify.py:444` -- `assert best_result is not None` -> `RuntimeError` with position/remaining context
   - `mapping.py:134` -- `assert p1.stdout is not None` -> `RuntimeError` + kill process

2. **Add `logging` to all modules:**
   - `logger = logging.getLogger(__name__)` at top of each module
   - `alleles.py`: Log cluster detection, valley splitting, peak refinement
   - `calling.py`: Log Clair3 invocation, VCF filtering stats, genotype disambiguation
   - `classify.py`: Per-repeat classification at DEBUG, mutations at INFO
   - `mapping.py`: Pipeline start/completion, read counts
   - `tools.py`: Command execution at DEBUG
   - Keep `click.echo` for user-facing progress in `cli.py`

3. **Add `--verbose`/`-v` and `--quiet`/`-q` flags** to CLI main group:
   - Wire to `logging.basicConfig()` with level mapping
   - `-v` = INFO, `-vv` = DEBUG, `-q` = WARNING only, default = WARNING

4. **Fix potential deadlock in `mapping.py`:**
   - Read stderr before `p1.wait()` or use `communicate()` to avoid buffer fill deadlock

### 1B. Community Health

**All community files go in `.github/` to keep repo root clean, except `CHANGELOG.md` (convention).**

5. **`.github/CONTRIBUTING.md`** -- dev setup, code style, testing requirements, PR process, integration test data workflow

6. **`CHANGELOG.md`** (root) -- backfilled from git history, Keep a Changelog format, covering v0.1.0 through current

7. **`.github/SECURITY.md`** -- responsible disclosure policy, contact email, scope

8. **`.github/CODE_OF_CONDUCT.md`** -- Contributor Covenant v2.1

9. **`.github/ISSUE_TEMPLATE/bug_report.yml`** -- structured bug report with environment, steps to reproduce, expected/actual behavior

10. **`.github/ISSUE_TEMPLATE/feature_request.yml`** -- structured feature request

11. **`.github/pull_request_template.md`** -- checklist: tests added, `make check` passes, docs updated, CHANGELOG updated

12. **`.github/CODEOWNERS`** -- `* @berntpopp`

13. **Git tags for existing releases:**
    - `v0.1.0` on `68e83db`
    - `v0.1.1` on `2d5b1b6`
    - `v0.1.2` on `a937335`
    - `v0.2.0` on `3df1152`

---

## Phase 2: Reproducibility -- Packaging + Scientific Validation

**Paths:** 5 (Packaging, Distribution & Pinning) + 7 (Scientific Reproducibility)
**Effort:** ~2-3 days
**Score impact:** Packaging 8->9+, adds scientific trust dimension

### 2A. Docker Optimization

Current Dockerfile uses unpinned `condaforge/mambaforge:latest`, single-stage, no `.dockerignore`, runs as root. Redesign for small, fast, secure builds.

1. **Switch base to `mambaorg/micromamba:<pinned>-bookworm-slim`** -- ~80MB vs ~750MB+

2. **Multi-stage build:**
   - Stage 1 (builder): install conda env + pip install muc_one_span
   - Stage 2 (runtime): copy only resolved env, no build tools/source/cache

3. **Add `docker/.dockerignore`** -- exclude `.git/`, `.planning/`, `tests/`, `docs/`, `*.md`, `__pycache__/`, `.mypy_cache/`, `.ruff_cache/`

4. **Layer caching strategy:**
   - COPY `conda/environment.yml` + `pyproject.toml` first
   - Install deps
   - Then COPY source
   - Source changes don't bust dependency cache

5. **Non-root user** -- run as `$MAMBA_USER`

6. **OCI labels:**
   ```
   org.opencontainers.image.source
   org.opencontainers.image.description
   org.opencontainers.image.licenses
   ```

7. **Activate env via `ENV PATH`** instead of `mamba run` wrapping in ENTRYPOINT -- faster startup

8. **BuildKit caching in CI** -- add `cache-from`/`cache-to` with GitHub Actions cache in `docker.yml`

### 2B. Conda Pinning

9. **Pin exact versions in `conda/environment.yml`:**
   ```yaml
   - minimap2=2.28
   - samtools=1.21
   - bcftools=1.21
   - htslib=1.21
   ```

10. **Add `conda/environment-dev.yml`** -- loose ranges (`>=`) for development flexibility

11. **Generate and commit `conda/environment.lock.yml`** -- full resolved dependency tree via `conda env export --no-builds`

### 2C. Distribution (scoped -- no PyPI/conda for now)

12. **Add Singularity/Apptainer definition at `docker/muconespan.def`** -- bootstraps from the GHCR Docker image, no new top-level directory

13. **Add GitHub Release workflow `.github/workflows/release.yml`** -- triggered on `v*` tags, generates release notes via `softprops/action-gh-release`

14. **Record tool versions in pipeline output:**
    - Add `get_tool_versions()` to `tools.py`
    - Include tool versions dict in `summary.json` automatically
    - Captures minimap2, samtools, bcftools versions at runtime

### 2D. Scientific Reproducibility

15. **Add `scripts/validate_catalog.py`** -- runs classification on MucOneUp samples, produces structured validation report (expected vs observed structure strings, TSV/JSON)

16. **Add validation CI job** -- separate job in `test.yml`, runs after integration tests, produces validation artifact

17. **Add `docs/reference/limitations.md`** -- document PCR bias edge cases, Clair3 detection limits for long VNTRs, known failure modes

---

## Phase 3: Code Quality -- Test Coverage + Type Safety

**Paths:** 4 (Type Safety & Code Structure) + 1 (Test Coverage Hardening)
**Effort:** ~4-5 days
**Score impact:** Code Quality 9->10, Testing 8->9

### 3A. Type Safety (do first -- tests reference new types)

1. **Add `TypedDict` for classification results in `classify.py`:**
   - `RepeatClassification` -- type, match, confidence, closest_match, edit_distance, identity_pct, differences, classification, frameshift, mutation_name, parent_repeat, index
   - `SequenceClassification` -- structure, repeats, mutations_detected, cumulative_offset, allele_confidence, exact_match_pct

2. **Add `TypedDict` for allele results in `alleles.py`:**
   - `AlleleInfo` -- length, reads, canonical_repeats, contig_name, cluster_contigs
   - `AlleleResult` -- allele_1, allele_2, homozygous, same_length

3. **Decompose `classify_sequence()` (216 lines) into 3 helpers:**
   - `_forward_classify()` -- forward pass left-to-right with offset tracking
   - `_bidirectional_fallback()` -- backward anchoring when forward stalls
   - `_compute_summary()` -- confidence scores and summary statistics

4. **Extract VCF parsing from `calling.py` into `vcf.py`:**
   - Move: `parse_vcf_genotypes()`, `parse_vcf_variants()`, `filter_vcf()`
   - Keep in `calling.py`: `run_clair3()`, `call_variants_per_allele()`, `disambiguate_same_length_alleles()`

### 3B. Test Coverage Hardening

5. **`mapping.py` (50% -> 90%):** Mock `subprocess.Popen`, test FileNotFoundError paths for minimap2 and samtools, SIGPIPE handling, non-zero returncode paths

6. **`cli.py` (64% -> 85%):** CliRunner test for `run` subcommand with all pipeline stages mocked

7. **`classify.py` (67% -> 85%):** Test `_classify_backward` with known after-repeat anchors, test bidirectional fallback with sequences where forward pass stalls

8. **`alleles.py` (68% -> 85%):** Test `_split_cluster_by_indel` with mocked SAM, test `len(valleys) < 2` early return, test `abs(v2 - v1) < 3` rejection

9. **Raise CI threshold** from `--cov-fail-under=70` to `--cov-fail-under=80`

---

## Phase 4: Features -- HTML Report + Performance

**Paths:** 8 (HTML Report Output) + 6 (Performance & Scalability)
**Effort:** ~4-5 days
**Score impact:** CLI 9->10, Performance 9->10, significant adoption value

### 4A. HTML Report (optional dependency)

1. **Add `[report]` optional extra in `pyproject.toml`:**
   ```toml
   [project.optional-dependencies]
   report = ["jinja2>=3.1.0"]
   dev = [..., "types-jinja2>=2.11"]
   ```

2. **Add `src/muc_one_span/report.py`:**
   - `generate_report(summary, output_path, sample_name, tool_versions)` -> Path
   - Graceful `ImportError` with helpful message if Jinja2 not installed
   - Renders from `summary.json` data (primary), optionally enriched from full `repeats.json`

3. **Add `src/muc_one_span/templates/report.html.j2`:**
   Self-contained HTML report with these sections:

   **a) Header:**
   - Sample name, pipeline version, tool versions, generation timestamp

   **b) Allele Summary Panel:**
   - Side-by-side cards for allele 1 & 2
   - Repeat count (canonical + total), read depth, contig name
   - Homozygous/heterozygous badge
   - Read distribution bar (allele 1 vs allele 2 proportions)

   **c) Repeat Structure Visualization:**
   - Per-allele color-coded linear map of all repeat units
   - Each unit as a colored block: pre-repeats, canonical types (A-Z), variants, mutations, after-repeats
   - Tooltip on hover: type, confidence, edit distance, differences
   - Mutation positions highlighted with distinct marker (pattern/icon, not just color)
   - Visual legend for repeat type categories

   **d) Mutation Report Table:**
   - Per detected mutation: repeat index, parent type, mutation name
   - Frameshift flag, VCF support, VCF QUAL score, boundary flag
   - Template match vs computed distinction
   - Confidence with visual bar

   **e) Quality Metrics:**
   - `allele_confidence` as gauge/progress bar per allele
   - `exact_match_pct` as gauge/progress bar per allele
   - `cumulative_offset` indicator
   - Read coverage summary

   **f) Detailed Repeat Table (collapsible):**
   - Full `repeats[]` array: index, type, match method, confidence, closest match, edit distance, differences list

   **g) Reproducibility Footer:**
   - Pipeline version, tool versions, CLI arguments, timestamp

   **UI/UX standards:**
   - Modern, clean design: system fonts, generous whitespace, subtle borders
   - Responsive CSS grid (readable on laptop + printed A4)
   - WCAG AA color contrast minimum
   - Colorblind-safe: patterns/icons supplement colors
   - Collapsible sections for detailed data (scannable overview)
   - Print stylesheet (`@media print`) -- clean single-column, no interactive elements
   - Dark/light mode via `prefers-color-scheme`
   - Self-contained: zero external deps, inline CSS, <100KB target
   - Aesthetic inspiration: MultiQC / nf-core report conventions

4. **Wire into CLI:**
   - `--report` flag on `run` subcommand (opt-in, errors helpfully if Jinja2 missing)
   - Standalone `report` subcommand for regenerating from existing `summary.json`

5. **Add `tests/unit/test_report.py`:**
   - Test HTML generation from sample summary dict
   - Test self-containedness (no external URLs)
   - Test graceful error when Jinja2 not installed

### 4B. Performance

6. **Parallelize per-allele processing** in `calling.py`:
   - `ThreadPoolExecutor(max_workers=2)` for independent allele variant calling
   - Skip for homozygous or same_length (sequential disambiguation needed)

7. **Add streaming `run_tool_iter()` in `tools.py`:**
   - Yields stdout lines without buffering all in memory
   - Same PATH sanitization and error handling as `run_tool()`

8. **Use streaming in `alleles.py`:**
   - `refine_peak_contig` and `_split_cluster_by_indel` SAM parsing via `run_tool_iter()`

9. **Add `scripts/benchmark.py`:**
   - Times each pipeline stage on MucOneUp samples
   - Outputs per-stage timing table

---

## Implementation Order Summary

| Phase | Paths | Effort | Cumulative Score |
|-------|-------|--------|-----------------|
| 1. Quick Wins | 2 + 3 | ~2 days | 8.5 -> 8.7 |
| 2. Reproducibility | 5 + 7 | ~2-3 days | 8.7 -> 8.9 |
| 3. Code Quality | 4 + 1 | ~4-5 days | 8.9 -> 9.2 |
| 4. Features | 8 + 6 | ~4-5 days | 9.2 -> 9.5+ |

**Total estimated effort:** ~12-15 days
**Target score:** 9.5+ / 10

---

## Constraints & Decisions

- **No PyPI or conda releases** for now
- **Jinja2 is optional** -- `pip install muc_one_span[report]`
- **Report generation is opt-in** -- `--report` flag, not default
- **Clean repo root** -- community files in `.github/`, Singularity def in `docker/`, only `CHANGELOG.md` at root
- **Docker optimized** for small size (slim base, multi-stage), fast builds (layer caching, BuildKit), and security (non-root, pinned base)
- **All plans go in `.planning/`** per project convention
