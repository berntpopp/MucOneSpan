# Phase 1: Quick Wins -- Error Handling + Community Health

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix production safety gaps (assert statements, missing logging, potential deadlock) and add all community health files to prepare the project for external contributors and higher review scores.

**Architecture:** Two independent workstreams -- (A) error handling improvements touch `classify.py`, `mapping.py`, `cli.py`, `alleles.py`, `calling.py`, `tools.py`, `ladder.py`; (B) community files are all new files under `.github/` and root `CHANGELOG.md`. Both can be committed independently.

**Tech Stack:** Python stdlib `logging`, Click CLI options, Markdown, GitHub YAML templates

---

### Task 1: Replace assert in classify.py with proper error handling

**Files:**
- Modify: `src/muc_one_span/classify.py:444`
- Test: `tests/unit/test_classify.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_classify.py`:

```python
def test_classify_sequence_raises_on_no_match(mocker):
    """classify_sequence raises RuntimeError when no repeat match is found."""
    from muc_one_span.classify import classify_sequence
    from muc_one_span.config import RepeatDictionary

    # Create a minimal repeat dict with no matching sequences
    rd = RepeatDictionary(
        repeats={},
        mutations={},
        pre_repeat_ids=["1", "2", "3", "4", "5"],
        after_repeat_ids=["6", "7", "8", "9"],
        canonical_order=[],
    )

    # A sequence that can't match anything -- short enough to hit the fallback
    sequence = "A" * 60

    with pytest.raises(RuntimeError, match="Classification failed"):
        classify_sequence(sequence, rd)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/unit/test_classify.py::test_classify_sequence_raises_on_no_match -v --no-cov`
Expected: FAIL with `AssertionError` (the current assert fires instead of RuntimeError)

- [ ] **Step 3: Replace assert with guarded raise**

In `src/muc_one_span/classify.py`, replace line 444:

```python
# Before:
        assert best_result is not None

# After:
        if best_result is None:
            raise RuntimeError(
                f"Classification failed: no match found at position {pos} "
                f"(remaining: {remaining} bp, repeat index: {repeat_index})"
            )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/unit/test_classify.py::test_classify_sequence_raises_on_no_match -v --no-cov`
Expected: PASS

- [ ] **Step 5: Run full test suite to check for regressions**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All existing tests pass

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/classify.py tests/unit/test_classify.py
git commit -m "fix: replace assert with RuntimeError in classify_sequence

The assert statement would raise an unhelpful AssertionError in
production. Replace with a descriptive RuntimeError that includes
the position and remaining sequence length for debugging."
```

---

### Task 2: Replace assert in mapping.py and fix potential deadlock

**Files:**
- Modify: `src/muc_one_span/mapping.py:134-148`
- Test: `tests/unit/test_mapping.py`

- [ ] **Step 1: Write the failing test for assert replacement**

Add to `tests/unit/test_mapping.py`:

```python
def test_run_mapping_pipeline_stdout_none_raises(mocker):
    """_run_mapping_pipeline raises RuntimeError if p1.stdout is None."""
    from muc_one_span.mapping import _run_mapping_pipeline

    mock_p1 = mocker.MagicMock()
    mock_p1.stdout = None
    mock_p1.kill = mocker.MagicMock()
    mock_p1.wait = mocker.MagicMock()

    mocker.patch(
        "muc_one_span.mapping.subprocess.Popen",
        side_effect=[mock_p1],
    )

    with pytest.raises(RuntimeError, match="minimap2 process stdout was not captured"):
        _run_mapping_pipeline(
            input_path=Path("/tmp/test.fastq"),
            reference_path=Path("/tmp/ref.fa"),
            bam_path=Path("/tmp/out.bam"),
            threads=1,
        )

    mock_p1.kill.assert_called_once()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/unit/test_mapping.py::test_run_mapping_pipeline_stdout_none_raises -v --no-cov`
Expected: FAIL with `AssertionError` (the current assert fires instead of RuntimeError)

- [ ] **Step 3: Replace assert and fix stderr reading order**

In `src/muc_one_span/mapping.py`, replace lines 133-148:

```python
# Before (lines 133-148):
    # Allow p1 to receive SIGPIPE if p2 exits early
    assert p1.stdout is not None
    p1.stdout.close()

    _, p2_stderr = p2.communicate()
    p1.wait()

    if p1.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed with exit code {p1.returncode}.\n"
            f"stderr: {p1.stderr.read().decode() if p1.stderr else ''}"
        )
    if p2.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed with exit code {p2.returncode}.\nstderr: {p2_stderr.decode()}"
        )

# After:
    # Allow p1 to receive SIGPIPE if p2 exits early
    if p1.stdout is None:
        p1.kill()
        p1.wait()
        raise RuntimeError("minimap2 process stdout was not captured")
    p1.stdout.close()

    # Read stderr before wait() to avoid potential deadlock if stderr buffer fills
    p1_stderr = p1.stderr.read().decode() if p1.stderr else ""
    _, p2_stderr = p2.communicate()
    p1.wait()

    if p1.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed with exit code {p1.returncode}.\n"
            f"stderr: {p1_stderr}"
        )
    if p2.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed with exit code {p2.returncode}.\nstderr: {p2_stderr.decode()}"
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/unit/test_mapping.py::test_run_mapping_pipeline_stdout_none_raises -v --no-cov`
Expected: PASS

- [ ] **Step 5: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/mapping.py tests/unit/test_mapping.py
git commit -m "fix: replace assert with RuntimeError and fix potential deadlock in mapping

Replace assert p1.stdout is not None with a proper RuntimeError that
kills the process. Read p1.stderr before p1.wait() to prevent potential
deadlock if the stderr buffer fills."
```

---

### Task 3: Add logging to tools.py

**Files:**
- Modify: `src/muc_one_span/tools.py`
- Test: `tests/unit/test_tools.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_tools.py`:

```python
import logging


def test_run_tool_logs_command(mocker, caplog):
    """run_tool logs the command at DEBUG level."""
    mocker.patch(
        "muc_one_span.tools.subprocess.run",
        return_value=mocker.MagicMock(returncode=0, stdout="output"),
    )
    mocker.patch("muc_one_span.tools.shutil.which", return_value="/usr/bin/echo")

    with caplog.at_level(logging.DEBUG, logger="muc_one_span.tools"):
        from muc_one_span.tools import run_tool

        run_tool(["echo", "hello"])

    assert any("echo hello" in record.message for record in caplog.records)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/unit/test_tools.py::test_run_tool_logs_command -v --no-cov`
Expected: FAIL (no log records produced)

- [ ] **Step 3: Add logging to tools.py**

Add after the existing imports (after line 9) in `src/muc_one_span/tools.py`:

```python
import logging

logger = logging.getLogger(__name__)
```

Then in the `run_tool()` function, add a DEBUG log before the subprocess call:

```python
    logger.debug("Running: %s", " ".join(str(c) for c in cmd))
```

And in `check_tools()`, add an INFO log:

```python
    logger.info("All required tools found: %s", ", ".join(tools))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/unit/test_tools.py::test_run_tool_logs_command -v --no-cov`
Expected: PASS

- [ ] **Step 5: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/tools.py tests/unit/test_tools.py
git commit -m "feat: add logging to tools.py

Log command execution at DEBUG level and tool availability at INFO
level for improved observability when using the library API."
```

---

### Task 4: Add logging to remaining modules

**Files:**
- Modify: `src/muc_one_span/alleles.py`
- Modify: `src/muc_one_span/calling.py`
- Modify: `src/muc_one_span/classify.py`
- Modify: `src/muc_one_span/mapping.py`
- Modify: `src/muc_one_span/ladder.py`

Note: `consensus.py` already has logging. `config.py` and `version.py` are too simple to need it.

- [ ] **Step 1: Add logging to alleles.py**

Add after line 15 (after existing imports) in `src/muc_one_span/alleles.py`:

```python
import logging

logger = logging.getLogger(__name__)
```

Add log calls at key decision points:

In `detect_alleles()`, after peaks are found:
```python
    logger.info("Detected %d peak(s) from %d contigs", len(peaks), len(counts))
```

In `_split_cluster_by_indel()`, at the valley detection:
```python
    logger.debug("Indel valley splitting: found %d valleys", len(valleys))
```

In `refine_peak_contig()`, after refinement:
```python
    logger.debug("Refined peak contig: %s (AS=%.1f)", best_contig, best_score)
```

- [ ] **Step 2: Add logging to calling.py**

Add after line 5 (after existing imports) in `src/muc_one_span/calling.py`:

```python
import logging

logger = logging.getLogger(__name__)
```

Add log calls:

In `run_clair3()`, before invocation:
```python
    logger.info("Running Clair3 on %s", bam_path.name)
```

In `filter_vcf()`, after filtering:
```python
    logger.debug("VCF filtering: %d variants passed min_qual=%.1f", count, min_qual)
```

In `call_variants_per_allele()`, for disambiguation:
```python
    logger.info("Same-length alleles detected, using disambiguation")
```

- [ ] **Step 3: Add logging to classify.py**

Add after line 3 (after existing imports) in `src/muc_one_span/classify.py`:

```python
import logging

logger = logging.getLogger(__name__)
```

Add log calls:

In `classify_sequence()`, at start:
```python
    logger.info("Classifying sequence of %d bp", len(sequence))
```

Per-repeat at DEBUG:
```python
    logger.debug("Repeat %d: type=%s match=%s conf=%.3f", repeat_index, result["type"], result["match"], result.get("confidence", 0))
```

On mutation detection:
```python
    logger.info("Mutation detected at repeat %d: %s", repeat_index, result.get("mutation_name", "unknown"))
```

- [ ] **Step 4: Add logging to mapping.py**

Add after line 6 (after existing imports) in `src/muc_one_span/mapping.py`:

```python
import logging

logger = logging.getLogger(__name__)
```

Add log calls:

In `_run_mapping_pipeline()`:
```python
    logger.info("Starting minimap2 | samtools sort pipeline (%d threads)", threads)
    # ... after completion:
    logger.info("Mapping pipeline completed: %s", bam_path.name)
```

In `map_reads()`:
```python
    logger.info("Mapping reads from %s to %s", input_path.name, reference_path.name)
```

- [ ] **Step 5: Add logging to ladder.py**

Add after line 6 (after existing imports) in `src/muc_one_span/ladder.py`:

```python
import logging

logger = logging.getLogger(__name__)
```

In `generate_ladder()` (or equivalent function):
```python
    logger.info("Generating reference ladder with %d contigs", num_contigs)
```

- [ ] **Step 6: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass (logging doesn't change behavior)

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/alleles.py src/muc_one_span/calling.py src/muc_one_span/classify.py src/muc_one_span/mapping.py src/muc_one_span/ladder.py
git commit -m "feat: add structured logging to all pipeline modules

Each module now has a logger for diagnostic output at DEBUG/INFO
levels. This improves observability when using the library API and
when --verbose is enabled in the CLI (added in next commit)."
```

---

### Task 5: Add --verbose and --quiet CLI flags

**Files:**
- Modify: `src/muc_one_span/cli.py:14-18`
- Test: `tests/unit/test_cli.py`

- [ ] **Step 1: Write the failing test for --verbose**

Add to `tests/unit/test_cli.py`:

```python
import logging


def test_verbose_flag_sets_info_level(runner):
    """--verbose sets logging to INFO level."""
    result = runner.invoke(main, ["-v", "--version"])
    # Just verify the flag is accepted without error
    assert result.exit_code == 0


def test_double_verbose_sets_debug_level(runner):
    """--verbose --verbose sets logging to DEBUG level."""
    result = runner.invoke(main, ["-vv", "--version"])
    assert result.exit_code == 0


def test_quiet_flag_accepted(runner):
    """--quiet flag is accepted without error."""
    result = runner.invoke(main, ["-q", "--version"])
    assert result.exit_code == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_cli.py::test_verbose_flag_sets_info_level tests/unit/test_cli.py::test_double_verbose_sets_debug_level tests/unit/test_cli.py::test_quiet_flag_accepted -v --no-cov`
Expected: FAIL (no such option)

- [ ] **Step 3: Add --verbose and --quiet to the main group**

In `src/muc_one_span/cli.py`, replace lines 14-18:

```python
# Before:
@click.group()
@click.version_option(version=__version__, prog_name="muconespan")
def main() -> None:
    """muconespan: MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data."""

# After:
@click.group()
@click.version_option(version=__version__, prog_name="muconespan")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v for INFO, -vv for DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Suppress non-error output.")
@click.pass_context
def main(ctx: click.Context, verbose: int, quiet: bool) -> None:
    """muconespan: MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data."""
    import logging

    if quiet:
        level = logging.WARNING
    elif verbose >= 2:
        level = logging.DEBUG
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING

    logging.basicConfig(
        level=level,
        format="%(name)s %(levelname)s: %(message)s",
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_cli.py::test_verbose_flag_sets_info_level tests/unit/test_cli.py::test_double_verbose_sets_debug_level tests/unit/test_cli.py::test_quiet_flag_accepted -v --no-cov`
Expected: PASS

- [ ] **Step 5: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/cli.py tests/unit/test_cli.py
git commit -m "feat: add --verbose and --quiet CLI flags

-v enables INFO logging, -vv enables DEBUG, -q suppresses non-error
output. Wired to Python stdlib logging so all module loggers
respect the verbosity level."
```

---

### Task 6: Create CONTRIBUTING.md

**Files:**
- Create: `.github/CONTRIBUTING.md`

- [ ] **Step 1: Create the file**

```markdown
# Contributing to muconespan

Thank you for your interest in contributing to muconespan!

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/berntpopp/MucOneSpan.git
   cd MucOneSpan
   ```

2. Install with development dependencies:
   ```bash
   make dev
   ```

3. For integration tests, set up the conda environment with external tools:
   ```bash
   mamba env create -f conda/environment.yml
   ```

## Code Style

- **Linting:** [ruff](https://github.com/astral-sh/ruff) for both linting and formatting
- **Type checking:** [mypy](https://mypy-lang.org/) with strict configuration
- **Docstrings:** Google style with Args/Returns/Raises sections
- **Type hints:** Required on all function signatures

Run all quality checks:
```bash
make check
```

## Testing

### Unit tests (no external tools required)

```bash
make test-fast
```

### Integration tests (require minimap2, samtools, bcftools, Clair3)

```bash
export PATH="/path/to/conda/env/bin:$PATH"
uv run pytest tests/integration/ -v --no-cov
```

### Test requirements for PRs

- All new features must include unit tests
- All bug fixes must include a regression test
- Unit tests must pass: `make test-fast`
- Coverage must not decrease

## Pull Request Process

1. **Branch naming:** Use prefixes: `feat/`, `fix/`, `docs/`, `chore/`, `test/`
2. **Commit messages:** Use [Conventional Commits](https://www.conventionalcommits.org/): `feat:`, `fix:`, `docs:`, `chore:`, `test:`
3. **Before submitting:**
   - Run `make check` (lint, format, type check, tests)
   - Update `CHANGELOG.md` under the `[Unreleased]` section
   - Add or update documentation if applicable
4. **PR description:** Use the pull request template checklist

## Integration Test Data

Integration tests use simulated PacBio HiFi reads generated by [MucOneUp](https://github.com/berntpopp/muconeup). See `.planning/TESTING_WITH_MUCONEUP.md` for the data generation workflow.

## Questions?

Open an issue or start a discussion on the [GitHub repository](https://github.com/berntpopp/MucOneSpan).
```

- [ ] **Step 2: Commit**

```bash
git add .github/CONTRIBUTING.md
git commit -m "docs: add CONTRIBUTING.md with development guidelines"
```

---

### Task 7: Create CHANGELOG.md

**Files:**
- Create: `CHANGELOG.md` (root)

- [ ] **Step 1: Review git history for changelog entries**

Run: `git log --oneline --all` to confirm version tags and key changes.

- [ ] **Step 2: Create the file**

```markdown
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.3.0] - 2026-04-06

### Added
- Soft QUAL scoring with continuous confidence and boundary penalty

### Fixed
- Address review feedback on soft QUAL scoring

## [0.2.0] - 2025-12-XX

### Added
- MkDocs Material documentation site with GitHub Pages deploy
- Comprehensive README with deviations, usage, and validation

## [0.1.2] - 2025-XX-XX

### Added
- Indel-valley allele splitting for close allele pairs

## [0.1.1] - 2025-XX-XX

### Fixed
- Strip virtualenv PATH from subprocess calls for external tools
- Handle empty Clair3 VCFs and remove INFO/DP filter

## [0.1.0] - 2025-XX-XX

### Added
- Initial release of muconespan pipeline
- Reference ladder generation (20-150 repeat units)
- Read mapping with minimap2
- Allele length detection with peak finding
- Variant calling with Clair3
- Consensus building with bcftools
- Repeat unit classification with Vrbacka nomenclature
- Full pipeline CLI with individual subcommands
- Docker image with conda-based tool stack
- Unit and integration test suite

[Unreleased]: https://github.com/berntpopp/MucOneSpan/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/berntpopp/MucOneSpan/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/berntpopp/MucOneSpan/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/berntpopp/MucOneSpan/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/berntpopp/MucOneSpan/releases/tag/v0.1.0
```

**Note:** The `XX-XX` dates should be filled in from `git log --format="%ai" v0.1.0..v0.1.0` etc. The executor should look up exact dates from git tags.

- [ ] **Step 3: Commit**

```bash
git add CHANGELOG.md
git commit -m "docs: add CHANGELOG.md backfilled from git history"
```

---

### Task 8: Create SECURITY.md

**Files:**
- Create: `.github/SECURITY.md`

- [ ] **Step 1: Create the file**

```markdown
# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.3.x   | :white_check_mark: |
| < 0.3   | :x:                |

## Reporting a Vulnerability

If you discover a security vulnerability in muconespan, please report it responsibly.

**Do not open a public issue.**

Instead, please email: **bernt.popp@gmail.com**

Include:
- Description of the vulnerability
- Steps to reproduce
- Potential impact
- Suggested fix (if any)

You can expect an initial response within 72 hours.

## Scope

This policy covers:
- The muconespan Python package
- Docker images published to GHCR
- GitHub Actions workflows in this repository

This policy does **not** cover:
- External tools (minimap2, samtools, bcftools, Clair3)
- Third-party dependencies
```

**Note:** The executor should verify the contact email with the user or use the email from git config.

- [ ] **Step 2: Commit**

```bash
git add .github/SECURITY.md
git commit -m "docs: add SECURITY.md with responsible disclosure policy"
```

---

### Task 9: Create CODE_OF_CONDUCT.md

**Files:**
- Create: `.github/CODE_OF_CONDUCT.md`

- [ ] **Step 1: Create the file**

Use the [Contributor Covenant v2.1](https://www.contributor-covenant.org/version/2/1/code_of_conduct/) text. The full text is long -- use the standard template with the project contact email.

The file should start with:

```markdown
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone...
```

And include the full Contributor Covenant v2.1 text with the contact method set to the same email used in SECURITY.md.

- [ ] **Step 2: Commit**

```bash
git add .github/CODE_OF_CONDUCT.md
git commit -m "docs: add CODE_OF_CONDUCT.md (Contributor Covenant v2.1)"
```

---

### Task 10: Create issue templates and PR template

**Files:**
- Create: `.github/ISSUE_TEMPLATE/bug_report.yml`
- Create: `.github/ISSUE_TEMPLATE/feature_request.yml`
- Create: `.github/pull_request_template.md`
- Create: `.github/CODEOWNERS`

- [ ] **Step 1: Create bug report template**

`.github/ISSUE_TEMPLATE/bug_report.yml`:

```yaml
name: Bug Report
description: Report a bug in muconespan
title: "[Bug]: "
labels: ["bug"]
body:
  - type: textarea
    id: description
    attributes:
      label: Description
      description: A clear description of the bug.
    validations:
      required: true
  - type: textarea
    id: reproduce
    attributes:
      label: Steps to reproduce
      description: How to reproduce the behavior.
      placeholder: |
        1. Run `muconespan run --input ...`
        2. ...
    validations:
      required: true
  - type: textarea
    id: expected
    attributes:
      label: Expected behavior
      description: What you expected to happen.
    validations:
      required: true
  - type: textarea
    id: actual
    attributes:
      label: Actual behavior
      description: What actually happened. Include error messages or logs.
    validations:
      required: true
  - type: textarea
    id: environment
    attributes:
      label: Environment
      description: |
        - OS:
        - Python version:
        - muconespan version:
        - Installation method (pip/docker/conda):
      render: markdown
    validations:
      required: true
```

- [ ] **Step 2: Create feature request template**

`.github/ISSUE_TEMPLATE/feature_request.yml`:

```yaml
name: Feature Request
description: Suggest a new feature or improvement
title: "[Feature]: "
labels: ["enhancement"]
body:
  - type: textarea
    id: problem
    attributes:
      label: Problem
      description: What problem does this feature solve?
    validations:
      required: true
  - type: textarea
    id: solution
    attributes:
      label: Proposed solution
      description: How should this work?
    validations:
      required: true
  - type: textarea
    id: alternatives
    attributes:
      label: Alternatives considered
      description: Any alternative approaches you considered.
    validations:
      required: false
```

- [ ] **Step 3: Create PR template**

`.github/pull_request_template.md`:

```markdown
## Summary

<!-- Brief description of the changes -->

## Checklist

- [ ] `make check` passes (lint, format, type check, tests)
- [ ] Tests added for new functionality or bug fix
- [ ] `CHANGELOG.md` updated under `[Unreleased]`
- [ ] Documentation updated (if applicable)

## Related Issues

<!-- Link related issues: Fixes #123, Closes #456 -->
```

- [ ] **Step 4: Create CODEOWNERS**

`.github/CODEOWNERS`:

```
# Default owner for everything
* @berntpopp
```

- [ ] **Step 5: Commit**

```bash
git add .github/ISSUE_TEMPLATE/bug_report.yml .github/ISSUE_TEMPLATE/feature_request.yml .github/pull_request_template.md .github/CODEOWNERS
git commit -m "docs: add issue templates, PR template, and CODEOWNERS"
```

---

### Task 11: Run quality checks and final commit

- [ ] **Step 1: Run lint check**

Run: `uv run ruff check src/muc_one_span/`
Expected: No new errors

- [ ] **Step 2: Run type check**

Run: `uv run mypy src/muc_one_span/`
Expected: No new errors (the `click.pass_context` addition may need `ctx: click.Context` annotation -- already included in Step 3 of Task 5)

- [ ] **Step 3: Run full test suite with coverage**

Run: `uv run pytest tests/unit/ --cov=muc_one_span --cov-report=term-missing`
Expected: All pass, coverage >= 70% (unchanged or slightly improved)

- [ ] **Step 4: Fix any issues found and commit**

If lint/type/test issues are found, fix them and commit:

```bash
git add -u
git commit -m "fix: address lint and type check issues from Phase 1 changes"
```
