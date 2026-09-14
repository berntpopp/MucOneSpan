# Phase 4: Features -- HTML Report + Performance

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an optional self-contained HTML report rendered from pipeline JSON output with modern UI/UX, and improve performance with parallel allele processing and streaming tool output.

**Architecture:** Report feature uses Jinja2 as an optional dependency (`[report]` extra). Template is self-contained HTML with inline CSS. Performance improvements add `ThreadPoolExecutor` for parallel allele processing and a streaming variant of `run_tool`. The report subcommand and `--report` flag on `run` provide two entry points.

**Tech Stack:** Jinja2 (optional), HTML/CSS (inline, self-contained), Python `concurrent.futures`, subprocess streaming

---

### Task 1: Add optional Jinja2 dependency

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Add [report] optional extra**

In `pyproject.toml`, add a `report` optional dependency group after the existing `docs` group:

```toml
[project.optional-dependencies]
report = [
    "jinja2>=3.1.0",
]
dev = [
    "pytest>=8.0.0",
    "pytest-cov>=5.0.0",
    "pytest-mock>=3.14.0",
    "ruff==0.15.9",
    "mypy>=1.11.0",
    "types-PyYAML>=6.0",
    "types-jinja2>=2.11",
]
docs = [
    "mkdocs-material>=9.5",
    "mkdocstrings[python]>=0.25",
    "mkdocs-click>=0.8",
]
```

- [ ] **Step 2: Run uv sync to update lock file**

Run: `uv sync --extra dev --extra report`
Expected: Jinja2 installed

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "feat: add optional [report] dependency group for Jinja2

Install with: pip install muc_one_span[report]
Jinja2 is only required for HTML report generation."
```

---

### Task 2: Create report module with graceful import handling

**Files:**
- Create: `src/muc_one_span/report.py`
- Test: `tests/unit/test_report.py`

- [ ] **Step 1: Write failing tests**

Create `tests/unit/test_report.py`:

```python
"""Tests for HTML report generation."""

from __future__ import annotations

import pytest
from pathlib import Path


@pytest.fixture
def sample_summary() -> dict:
    """Minimal pipeline summary for testing."""
    return {
        "alleles": {
            "allele_1": {
                "length": 50,
                "reads": 554,
                "canonical_repeats": 41,
                "contig_name": "contig_41",
                "cluster_contigs": ["contig_41", "contig_42"],
            },
            "allele_2": {
                "length": 60,
                "reads": 432,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_51", "contig_52"],
            },
            "homozygous": False,
            "same_length": False,
        },
        "classifications": {
            "allele_1": {
                "structure": "1 2 3 4 5 X X X A B 6 7 8 9",
                "mutations": [
                    {
                        "repeat_index": 8,
                        "closest_type": "X",
                        "mutation_name": "dupC",
                        "template_match": True,
                        "frameshift": True,
                        "vcf_support": True,
                        "vcf_qual": 23.4,
                        "boundary": False,
                    }
                ],
            },
            "allele_2": {
                "structure": "1 2 3 4 5 X X X X X X 6 7 8 9",
                "mutations": [],
            },
        },
        "tool_versions": {
            "minimap2": "minimap2 2.28-r1209",
            "samtools": "samtools 1.21",
        },
        "pipeline_version": "0.3.0",
    }


@pytest.fixture
def sample_repeats() -> dict:
    """Detailed repeats data for testing."""
    return {
        "allele_1": {
            "structure": "1 2 3 4 5 X X X A B 6 7 8 9",
            "repeats": [
                {"type": "1", "match": "exact", "confidence": 1.0, "index": 1},
                {"type": "2", "match": "exact", "confidence": 1.0, "index": 2},
                {"type": "X", "match": "exact", "confidence": 1.0, "index": 6},
                {
                    "type": "unknown",
                    "match": "closest",
                    "closest_match": "A",
                    "edit_distance": 1,
                    "identity_pct": 98.3,
                    "confidence": 0.983,
                    "differences": [{"pos": 54, "ref": "G", "alt": "S", "type": "substitution"}],
                    "classification": "variant",
                    "index": 9,
                },
            ],
            "mutations_detected": [
                {
                    "repeat_index": 8,
                    "closest_type": "X",
                    "mutation_name": "dupC",
                    "template_match": True,
                    "frameshift": True,
                    "vcf_support": True,
                    "vcf_qual": 23.4,
                    "boundary": False,
                }
            ],
            "cumulative_offset": 1,
            "allele_confidence": 0.9971,
            "exact_match_pct": 92.9,
        },
    }


class TestGenerateReport:
    """Tests for generate_report function."""

    def test_creates_html_file(self, tmp_path, sample_summary):
        """Report generates a valid HTML file."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        result = generate_report(sample_summary, out, sample_name="test_sample")

        assert result == out
        assert out.exists()
        html = out.read_text()
        assert "<html" in html
        assert "test_sample" in html
        assert "muconespan" in html

    def test_report_self_contained(self, tmp_path, sample_summary):
        """Report must not reference external URLs."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        generate_report(sample_summary, out)
        html = out.read_text()

        # No external resource references
        assert "http://" not in html
        assert "https://" not in html

    def test_report_includes_allele_data(self, tmp_path, sample_summary):
        """Report includes allele detection data."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        generate_report(sample_summary, out, sample_name="test")
        html = out.read_text()

        assert "allele" in html.lower()
        assert "554" in html  # allele_1 reads
        assert "432" in html  # allele_2 reads

    def test_report_includes_mutations(self, tmp_path, sample_summary):
        """Report includes mutation data."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        generate_report(sample_summary, out)
        html = out.read_text()

        assert "dupC" in html
        assert "frameshift" in html.lower()

    def test_report_includes_tool_versions(self, tmp_path, sample_summary):
        """Report includes tool version metadata."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        generate_report(sample_summary, out)
        html = out.read_text()

        assert "minimap2" in html
        assert "2.28" in html

    def test_report_with_detailed_repeats(self, tmp_path, sample_summary, sample_repeats):
        """Report enriched with detailed repeat data."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        generate_report(sample_summary, out, detailed_repeats=sample_repeats)
        html = out.read_text()

        assert out.exists()
        assert "confidence" in html.lower()

    def test_report_creates_parent_dirs(self, tmp_path, sample_summary):
        """Report creates parent directories if they don't exist."""
        from muc_one_span.report import generate_report

        out = tmp_path / "nested" / "dir" / "report.html"
        generate_report(sample_summary, out)
        assert out.exists()

    def test_report_size_under_100kb(self, tmp_path, sample_summary):
        """Report file size should be under 100KB."""
        from muc_one_span.report import generate_report

        out = tmp_path / "report.html"
        generate_report(sample_summary, out)
        assert out.stat().st_size < 100_000


class TestReportMissingJinja2:
    """Tests for graceful handling when Jinja2 is not installed."""

    def test_import_error_message(self, tmp_path, sample_summary, monkeypatch):
        """Helpful error message when Jinja2 is missing."""
        import importlib
        import muc_one_span.report as report_module

        original_import = __builtins__.__import__ if hasattr(__builtins__, '__import__') else __import__

        def mock_import(name, *args, **kwargs):
            if name == "jinja2":
                raise ImportError("No module named 'jinja2'")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", mock_import)

        # Reload module to trigger import error path
        try:
            importlib.reload(report_module)
            out = tmp_path / "report.html"
            with pytest.raises(ImportError, match="report"):
                report_module.generate_report(sample_summary, out)
        finally:
            # Restore
            importlib.reload(report_module)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_report.py -v --no-cov`
Expected: FAIL (module doesn't exist)

- [ ] **Step 3: Create report.py**

Create `src/muc_one_span/report.py`:

```python
"""Self-contained HTML report generation.

Requires the ``jinja2`` package, which is an optional dependency.
Install with: ``pip install muc_one_span[report]``
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

try:
    from jinja2 import Environment, PackageLoader

    _HAS_JINJA2 = True
except ImportError:
    _HAS_JINJA2 = False

from muc_one_span.version import __version__


def generate_report(
    summary: dict,
    output_path: Path,
    sample_name: str = "unknown",
    tool_versions: dict[str, str] | None = None,
    detailed_repeats: dict | None = None,
) -> Path:
    """Generate a self-contained HTML report from pipeline results.

    Args:
        summary: Combined pipeline summary (from summary.json).
        output_path: Where to write the HTML file.
        sample_name: Sample identifier for the report header.
        tool_versions: Optional dict of tool name -> version string.
            If None, reads from summary["tool_versions"] if available.
        detailed_repeats: Optional detailed repeat data (from repeats.json)
            for the expanded repeat table.

    Returns:
        Path to the generated HTML report.

    Raises:
        ImportError: If Jinja2 is not installed.
    """
    if not _HAS_JINJA2:
        raise ImportError(
            "Jinja2 is required for report generation. "
            "Install with: pip install muc_one_span[report]"
        )

    env = Environment(
        loader=PackageLoader("muc_one_span", "templates"),
        autoescape=True,
    )
    template = env.get_template("report.html.j2")

    # Resolve tool versions from summary if not provided explicitly
    versions = tool_versions or summary.get("tool_versions", {})

    html = template.render(
        sample_name=sample_name,
        summary=summary,
        detailed_repeats=detailed_repeats,
        tool_versions=versions,
        pipeline_version=summary.get("pipeline_version", __version__),
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)
    return output_path
```

- [ ] **Step 4: Run tests to verify import/basic structure works**

Run: `uv run pytest tests/unit/test_report.py::TestGenerateReport::test_creates_html_file -v --no-cov`
Expected: FAIL (template not found yet -- that's expected, we create it next)

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/report.py tests/unit/test_report.py
git commit -m "feat: add report module with graceful Jinja2 import handling

generate_report() renders summary.json into self-contained HTML.
Raises helpful ImportError if Jinja2 not installed."
```

---

### Task 3: Create HTML report template

**Files:**
- Create: `src/muc_one_span/templates/report.html.j2`

- [ ] **Step 1: Create template directory**

Run: `mkdir -p src/muc_one_span/templates`

- [ ] **Step 2: Create the HTML template**

Create `src/muc_one_span/templates/report.html.j2`. This is a large file -- the executor should use the `frontend-design` skill or create a modern, self-contained HTML template with these requirements:

**Sections:**
1. **Header** -- sample name, pipeline version, generation date
2. **Allele Summary** -- side-by-side cards with repeat count, reads, contig, homozygous badge, read distribution bar
3. **Repeat Structure Visualization** -- per-allele color-coded linear map of repeat units (pre=blue, canonical=gray, variable=green, after=orange, mutation=red). Each unit is a small colored block. Tooltip with type + confidence.
4. **Mutation Report Table** -- repeat index, parent type, mutation name, frameshift badge, VCF support badge, QUAL score, boundary badge, confidence bar
5. **Quality Metrics** -- allele_confidence gauge, exact_match_pct gauge, cumulative_offset, read coverage
6. **Detailed Repeat Table** (collapsible via `<details>`) -- full repeat array with index, type, match, confidence, edit distance, differences
7. **Reproducibility Footer** -- pipeline version, tool versions, timestamp

**Design requirements:**
- Self-contained: zero external dependencies, all CSS inline in `<style>` block
- System fonts (`system-ui, -apple-system, sans-serif`)
- CSS custom properties for theming
- CSS Grid layout, responsive
- `prefers-color-scheme` for dark/light mode
- `@media print` stylesheet for clean printout
- WCAG AA color contrast
- Colorblind-safe: mutation blocks get a pattern overlay (diagonal stripes) in addition to color
- Max ~80KB total
- Semantic HTML with ARIA labels on interactive elements
- Tooltips via CSS `title` attributes (no JS required)
- Inspired by MultiQC / nf-core report aesthetics

**Template data available:**
- `{{ sample_name }}` -- sample ID
- `{{ pipeline_version }}` -- e.g. "0.3.0"
- `{{ generated_at }}` -- UTC timestamp
- `{{ tool_versions }}` -- dict of tool -> version string
- `{{ summary.alleles }}` -- allele detection result (allele_1, allele_2, homozygous, same_length)
- `{{ summary.classifications }}` -- per-allele structure + mutations
- `{{ detailed_repeats }}` -- optional detailed repeat data (repeats[], mutations_detected[], confidence scores)

**Minimal skeleton (the executor should expand this into the full template):**

```html
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>muconespan Report: {{ sample_name }}</title>
  <style>
    :root {
      --bg: #ffffff; --fg: #1a1a2e; --card-bg: #f8f9fa;
      --border: #dee2e6; --accent: #4361ee;
      --pre: #6baed6; --canonical: #969696; --variable: #74c476;
      --after: #fd8d3c; --mutation: #e6550d; --unknown: #bdbdbd;
      --success: #2dc653; --warning: #ff9f1c; --danger: #e71d36;
    }
    @media (prefers-color-scheme: dark) {
      :root {
        --bg: #1a1a2e; --fg: #e8e8e8; --card-bg: #16213e;
        --border: #3a3a5c;
      }
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: system-ui, -apple-system, sans-serif;
      background: var(--bg); color: var(--fg);
      max-width: 1000px; margin: 0 auto; padding: 2rem 1rem;
      line-height: 1.6;
    }
    /* ... complete CSS for all sections ... */
    .repeat-map { display: flex; flex-wrap: wrap; gap: 1px; margin: 1rem 0; }
    .repeat-unit {
      width: 14px; height: 40px; border-radius: 2px;
      position: relative; cursor: default;
    }
    .repeat-unit.mutation {
      background: var(--mutation);
      background-image: repeating-linear-gradient(
        45deg, transparent, transparent 3px,
        rgba(255,255,255,0.3) 3px, rgba(255,255,255,0.3) 6px
      );
    }
    .confidence-bar {
      height: 8px; border-radius: 4px; background: var(--border);
    }
    .confidence-bar-fill {
      height: 100%; border-radius: 4px; background: var(--accent);
    }
    table { width: 100%; border-collapse: collapse; margin: 1rem 0; }
    th, td { padding: 0.5rem; text-align: left; border-bottom: 1px solid var(--border); }
    .badge {
      display: inline-block; padding: 0.15rem 0.5rem;
      border-radius: 1rem; font-size: 0.75rem; font-weight: 600;
    }
    .badge-success { background: var(--success); color: white; }
    .badge-danger { background: var(--danger); color: white; }
    .badge-warning { background: var(--warning); color: white; }
    details { margin: 1rem 0; }
    summary { cursor: pointer; font-weight: 600; }
    @media print {
      body { max-width: 100%; padding: 0; }
      details { display: block; }
      details > summary { display: none; }
      .no-print { display: none; }
    }
  </style>
</head>
<body>
  <header>
    <h1>MUC1 VNTR Analysis Report</h1>
    <p><strong>{{ sample_name }}</strong> &mdash; {{ generated_at }}</p>
    <p>muconespan v{{ pipeline_version }}</p>
  </header>

  <!-- Allele Summary -->
  <!-- Repeat Structure -->
  <!-- Mutation Table -->
  <!-- Quality Metrics -->
  <!-- Detailed Repeat Table (collapsible) -->
  <!-- Reproducibility Footer -->

</body>
</html>
```

The executor should fully implement all sections following the design requirements above. Use the `frontend-design` skill for high-quality implementation.

- [ ] **Step 3: Run report tests**

Run: `uv run pytest tests/unit/test_report.py -v --no-cov`
Expected: All pass

- [ ] **Step 4: Commit**

```bash
git add src/muc_one_span/templates/report.html.j2
git commit -m "feat: add self-contained HTML report template

Modern, accessible report with allele summary, repeat structure
visualization, mutation table, quality metrics, and detailed repeat
data. Self-contained (inline CSS, no external deps), printable,
dark/light mode, colorblind-safe."
```

---

### Task 4: Wire report into CLI

**Files:**
- Modify: `src/muc_one_span/cli.py`

- [ ] **Step 1: Add --report flag to run subcommand**

In `src/muc_one_span/cli.py`, add option to the `run` command (after `--min-qual`):

```python
@click.option(
    "--report/--no-report",
    default=False,
    help="Generate HTML report (requires jinja2: pip install muc_one_span[report]).",
)
```

Add `report: bool` parameter to the `run()` function.

At the end of `run()`, after writing `summary.json` and before "Pipeline complete":

```python
    if report:
        try:
            from muc_one_span.report import generate_report

            report_path = out / "report.html"
            generate_report(
                summary,
                report_path,
                sample_name=Path(input_path).stem,
                detailed_repeats=all_results,
            )
            click.echo(f"Report: {report_path}")
        except ImportError as e:
            click.echo(f"Warning: {e}", err=True)
```

- [ ] **Step 2: Add standalone report subcommand**

Add after the `run` command:

```python
@main.command()
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Path to summary.json from a previous run.",
)
@click.option("--output", "-o", type=click.Path(), default="report.html", help="Output HTML path.")
@click.option("--sample-name", "-s", default=None, help="Sample name for report header.")
@click.option(
    "--repeats",
    type=click.Path(exists=True),
    default=None,
    help="Path to repeats.json for detailed repeat table.",
)
def report(input_path: str, output: str, sample_name: str | None, repeats: str | None) -> None:
    """Generate an HTML report from pipeline results."""
    import json

    try:
        from muc_one_span.report import generate_report
    except ImportError as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1) from e

    summary = json.loads(Path(input_path).read_text())
    name = sample_name or Path(input_path).parent.name

    detailed = None
    if repeats:
        detailed = json.loads(Path(repeats).read_text())

    out_path = generate_report(summary, Path(output), sample_name=name, detailed_repeats=detailed)
    click.echo(f"Report written to {out_path}")
```

- [ ] **Step 3: Write tests for CLI integration**

Add to `tests/unit/test_cli.py`:

```python
class TestReportSubcommand:
    """Tests for the report subcommand."""

    def test_report_from_summary_json(self, tmp_path):
        """report subcommand generates HTML from summary.json."""
        import json

        summary = {
            "alleles": {
                "allele_1": {"length": 50, "reads": 100, "canonical_repeats": 41,
                             "contig_name": "c41", "cluster_contigs": ["c41"]},
                "allele_2": {"length": 60, "reads": 80, "canonical_repeats": 51,
                             "contig_name": "c51", "cluster_contigs": ["c51"]},
                "homozygous": False, "same_length": False,
            },
            "classifications": {
                "allele_1": {"structure": "1 2 3 X 6 7 8 9", "mutations": []},
                "allele_2": {"structure": "1 2 3 X X 6 7 8 9", "mutations": []},
            },
        }
        summary_path = tmp_path / "summary.json"
        summary_path.write_text(json.dumps(summary))
        out = tmp_path / "report.html"

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["report", "--input", str(summary_path), "--output", str(out), "--sample-name", "test"],
        )

        assert result.exit_code == 0, result.output
        assert out.exists()
        assert "test" in out.read_text()
```

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/unit/test_cli.py::TestReportSubcommand -v --no-cov`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/cli.py tests/unit/test_cli.py
git commit -m "feat: add --report flag to run and standalone report subcommand

--report on run generates HTML alongside JSON output.
report subcommand regenerates from existing summary.json.
Both handle missing Jinja2 gracefully."
```

---

### Task 5: Add streaming run_tool_iter to tools.py

**Files:**
- Modify: `src/muc_one_span/tools.py`
- Test: `tests/unit/test_tools.py`

- [ ] **Step 1: Write failing test**

Add to `tests/unit/test_tools.py`:

```python
def test_run_tool_iter_yields_lines(mocker):
    """run_tool_iter yields stdout lines without buffering."""
    from unittest.mock import MagicMock
    from muc_one_span.tools import run_tool_iter

    mock_proc = MagicMock()
    mock_proc.stdout = iter(["line1\n", "line2\n", "line3\n"])
    mock_proc.wait.return_value = 0
    mock_proc.returncode = 0

    mocker.patch("muc_one_span.tools.subprocess.Popen", return_value=mock_proc)

    lines = list(run_tool_iter(["echo", "hello"]))
    assert lines == ["line1\n", "line2\n", "line3\n"]


def test_run_tool_iter_raises_on_failure(mocker):
    """run_tool_iter raises RuntimeError on non-zero exit."""
    from unittest.mock import MagicMock
    from muc_one_span.tools import run_tool_iter

    mock_proc = MagicMock()
    mock_proc.stdout = iter([])
    mock_proc.wait.return_value = 1
    mock_proc.returncode = 1

    mocker.patch("muc_one_span.tools.subprocess.Popen", return_value=mock_proc)

    with pytest.raises(RuntimeError, match="failed"):
        list(run_tool_iter(["failing_tool"]))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_tools.py::test_run_tool_iter_yields_lines tests/unit/test_tools.py::test_run_tool_iter_raises_on_failure -v --no-cov`
Expected: FAIL (function doesn't exist)

- [ ] **Step 3: Add run_tool_iter to tools.py**

Add to `src/muc_one_span/tools.py` after `run_tool()`:

```python
def run_tool_iter(
    cmd: list[str],
    cwd: str | None = None,
) -> "Generator[str, None, None]":
    """Run an external tool and yield stdout lines without buffering.

    Unlike :func:`run_tool`, this does not capture all stdout in memory.
    Use for commands that may produce large output (e.g., ``samtools view``).

    Args:
        cmd: Command and arguments.
        cwd: Working directory for the command.

    Yields:
        Lines of stdout as strings.

    Raises:
        FileNotFoundError: If the tool is not found on PATH.
        RuntimeError: If the command exits with non-zero status.
    """
    import typing

    env = os.environ.copy()
    env["PATH"] = _clean_path_for_externals(env.get("PATH", ""))

    logger.debug("Running (streaming): %s", " ".join(str(c) for c in cmd))

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        cwd=cwd,
    )

    assert proc.stdout is not None  # guaranteed by PIPE
    yield from proc.stdout

    proc.wait()
    if proc.returncode != 0:
        stderr = proc.stderr.read() if proc.stderr else ""
        raise RuntimeError(
            f"Command failed: {' '.join(str(c) for c in cmd)}\n"
            f"Exit code: {proc.returncode}\nstderr: {stderr}"
        )
```

Also add at the top of the file:

```python
from __future__ import annotations

from collections.abc import Generator
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_tools.py::test_run_tool_iter_yields_lines tests/unit/test_tools.py::test_run_tool_iter_raises_on_failure -v --no-cov`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/tools.py tests/unit/test_tools.py
git commit -m "feat: add streaming run_tool_iter for large output handling

Yields stdout lines without buffering all in memory. Uses same
PATH sanitization as run_tool. For use with samtools view on
large BAMs."
```

---

### Task 6: Parallelize per-allele variant calling

**Files:**
- Modify: `src/muc_one_span/calling.py`
- Test: `tests/unit/test_calling.py`

- [ ] **Step 1: Write failing test**

Add to `tests/unit/test_calling.py`:

```python
def test_call_variants_parallel_execution(mocker, tmp_path):
    """call_variants_per_allele processes alleles in parallel."""
    from muc_one_span.calling import call_variants_per_allele

    alleles = {
        "allele_1": {
            "length": 50,
            "reads": 100,
            "canonical_repeats": 41,
            "contig_name": "contig_41",
            "cluster_contigs": ["contig_41"],
        },
        "allele_2": {
            "length": 60,
            "reads": 80,
            "canonical_repeats": 51,
            "contig_name": "contig_51",
            "cluster_contigs": ["contig_51"],
        },
        "homozygous": False,
        "same_length": False,
    }

    vcf_a1 = tmp_path / "allele_1" / "variants.vcf.gz"
    vcf_a2 = tmp_path / "allele_2" / "variants.vcf.gz"
    vcf_a1.parent.mkdir(parents=True)
    vcf_a2.parent.mkdir(parents=True)
    vcf_a1.touch()
    vcf_a2.touch()

    mocker.patch("muc_one_span.calling.extract_allele_reads", return_value=tmp_path / "reads.bam")
    mocker.patch("muc_one_span.calling.run_clair3", return_value=tmp_path / "clair3.vcf.gz")
    mocker.patch("muc_one_span.calling.filter_vcf", side_effect=[vcf_a1, vcf_a2])

    result = call_variants_per_allele(
        bam_path=tmp_path / "mapped.bam",
        reference_path=tmp_path / "ref.fa",
        alleles=alleles,
        output_dir=tmp_path,
        threads=4,
    )

    assert "allele_1" in result
    assert "allele_2" in result
```

- [ ] **Step 2: Add parallel execution to call_variants_per_allele**

In `src/muc_one_span/calling.py`, modify `call_variants_per_allele` to use `ThreadPoolExecutor` when both alleles are independent (not same_length, not homozygous):

```python
from concurrent.futures import ThreadPoolExecutor

def call_variants_per_allele(...) -> dict[str, Path]:
    # ... existing docstring and same_length/homozygous checks ...

    if alleles.get("same_length") or alleles.get("homozygous"):
        # Sequential processing required for disambiguation
        return disambiguate_same_length_alleles(...)

    def _process_allele(allele_key: str) -> tuple[str, Path]:
        # ... existing per-allele logic (extract, clair3, filter) ...
        return allele_key, filtered_vcf

    allele_keys = ["allele_1", "allele_2"]

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(_process_allele, k): k for k in allele_keys}
        results = {}
        for future in futures:
            key, vcf_path = future.result()
            results[key] = vcf_path

    return results
```

- [ ] **Step 3: Run tests**

Run: `uv run pytest tests/unit/test_calling.py -v --no-cov`
Expected: All pass

- [ ] **Step 4: Commit**

```bash
git add src/muc_one_span/calling.py tests/unit/test_calling.py
git commit -m "feat: parallelize per-allele variant calling

Use ThreadPoolExecutor(max_workers=2) for independent allele
processing. Falls back to sequential for same_length/homozygous
disambiguation."
```

---

### Task 7: Use streaming in alleles.py

**Files:**
- Modify: `src/muc_one_span/alleles.py`

- [ ] **Step 1: Update refine_peak_contig to use run_tool_iter**

In `src/muc_one_span/alleles.py`, update the import:

```python
from muc_one_span.tools import run_tool, run_tool_iter
```

In `refine_peak_contig()`, replace the `run_tool` call for `samtools view` with `run_tool_iter` and process lines as they stream:

```python
# Before:
    sam_output = run_tool([
        "samtools", "view", str(bam_path), *cluster_contigs
    ])
    for line in sam_output.strip().splitlines():
        ...

# After:
    for line in run_tool_iter([
        "samtools", "view", str(bam_path), *cluster_contigs
    ]):
        line = line.strip()
        if not line:
            continue
        ...
```

Apply the same pattern to `_split_cluster_by_indel()`.

- [ ] **Step 2: Run existing tests**

Run: `uv run pytest tests/unit/test_alleles.py -v --no-cov`
Expected: All pass (tests mock run_tool, need to also mock run_tool_iter -- update mocks if needed)

- [ ] **Step 3: Commit**

```bash
git add src/muc_one_span/alleles.py
git commit -m "feat: use streaming SAM parsing in allele refinement

Replace run_tool with run_tool_iter for samtools view calls in
refine_peak_contig and _split_cluster_by_indel. Avoids buffering
full SAM output in memory for large BAMs."
```

---

### Task 8: Add benchmark script

**Files:**
- Create: `scripts/benchmark.py`

- [ ] **Step 1: Create benchmark script**

Create `scripts/benchmark.py`:

```python
#!/usr/bin/env python3
"""Benchmark pipeline stage timings on MucOneUp test samples.

Usage:
    export PATH="/path/to/conda/env/bin:$PATH"
    python scripts/benchmark.py [--data-dir tests/data/generated]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark pipeline stages")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tests/data/generated"),
        help="Directory containing generated test samples",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tests/results/benchmarks"),
        help="Output directory for benchmark results",
    )
    args = parser.parse_args()

    if not args.data_dir.exists():
        print(f"Error: {args.data_dir} not found. Generate test data first.")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

    from muc_one_span.config import load_repeat_dictionary
    from muc_one_span.tools import check_tools

    check_tools(["minimap2", "samtools", "bcftools", "run_clair3.sh"])
    rd = load_repeat_dictionary()

    samples = sorted(args.data_dir.glob("sample_*"))
    results = []

    for sample_dir in samples:
        bam_files = list(sample_dir.glob("*.bam"))
        if not bam_files:
            continue

        bam = bam_files[0]
        sample_name = sample_dir.name
        out = args.output_dir / sample_name
        out.mkdir(parents=True, exist_ok=True)

        timings: dict[str, float] = {}
        print(f"\n{sample_name}:")

        # Stage 1: Mapping
        from muc_one_span.ladder import generate_ladder
        from muc_one_span.mapping import get_idxstats, map_reads

        ref = out / "ladder.fa"
        t0 = time.perf_counter()
        generate_ladder(rd, ref)
        mapped = map_reads(bam, ref, out, threads=4)
        timings["mapping"] = time.perf_counter() - t0
        print(f"  mapping: {timings['mapping']:.2f}s")

        # Stage 2: Allele detection
        from muc_one_span.alleles import detect_alleles, parse_idxstats

        t0 = time.perf_counter()
        idxstats = get_idxstats(mapped)
        counts = parse_idxstats(idxstats)
        alleles_result = detect_alleles(counts, min_coverage=10, bam_path=mapped)
        timings["alleles"] = time.perf_counter() - t0
        print(f"  alleles: {timings['alleles']:.2f}s")

        # Stage 3: Variant calling
        from muc_one_span.calling import call_variants_per_allele

        t0 = time.perf_counter()
        vcf_paths = call_variants_per_allele(
            mapped, ref, alleles_result, out, threads=4
        )
        timings["calling"] = time.perf_counter() - t0
        print(f"  calling: {timings['calling']:.2f}s")

        # Stage 4: Consensus
        from muc_one_span.consensus import build_consensus_per_allele

        t0 = time.perf_counter()
        consensus = build_consensus_per_allele(ref, vcf_paths, alleles_result, out, repeat_dict=rd)
        timings["consensus"] = time.perf_counter() - t0
        print(f"  consensus: {timings['consensus']:.2f}s")

        # Stage 5: Classification
        from muc_one_span.classify import classify_sequence

        t0 = time.perf_counter()
        for allele_key, fa_path in consensus.items():
            fa_lines = fa_path.read_text().strip().splitlines()
            sequence = "".join(line for line in fa_lines if not line.startswith(">"))
            classify_sequence(sequence, rd)
        timings["classify"] = time.perf_counter() - t0
        print(f"  classify: {timings['classify']:.2f}s")

        timings["total"] = sum(timings.values())
        print(f"  TOTAL: {timings['total']:.2f}s")

        results.append({"sample": sample_name, "timings": timings})

    # Write summary
    summary_path = args.output_dir / "benchmark_results.json"
    summary_path.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nResults written to {summary_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add scripts/benchmark.py
git commit -m "feat: add pipeline benchmark script

Times each pipeline stage on MucOneUp samples. Outputs per-stage
timings to benchmark_results.json for performance tracking."
```

---

### Task 9: Run quality checks

- [ ] **Step 1: Run lint**

Run: `uv run ruff check src/muc_one_span/`
Expected: No errors

- [ ] **Step 2: Run type check**

Run: `uv run mypy src/muc_one_span/`
Expected: No errors

- [ ] **Step 3: Run full test suite with coverage**

Run: `uv run pytest tests/unit/ --cov=muc_one_span --cov-report=term-missing --cov-fail-under=80`
Expected: All pass, coverage >= 80%

- [ ] **Step 4: Fix any issues and commit**

```bash
git add -u
git commit -m "fix: address lint and type check issues from Phase 4 changes"
```
