# Phase 2: Reproducibility -- Packaging + Scientific Validation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Optimize Docker for small/fast/secure builds, pin all tool versions for reproducibility, add Singularity support, add release automation, record tool versions in output, and add scientific validation infrastructure.

**Architecture:** Docker rebuild (multi-stage, slim base), conda pinning, new workflows, tool version capture in `tools.py`, validation script and docs.

**Tech Stack:** Docker (multi-stage, micromamba), conda, GitHub Actions, Singularity/Apptainer, Python stdlib

---

### Task 1: Add .dockerignore

**Files:**
- Create: `docker/.dockerignore`

- [ ] **Step 1: Create .dockerignore**

Create `docker/.dockerignore`:

```
.git
.github
.planning
.mypy_cache
.ruff_cache
.pytest_cache
__pycache__
*.pyc
tests/
docs/
scripts/
*.md
!README.md
.pre-commit-config.yaml
.gitignore
mkdocs.yml
uv.lock
CLAUDE.md
CHANGELOG.md
```

**Note:** Since the Dockerfile is in `docker/` but the build context is the repo root (context: `.` in the workflow), this file should be at the repo root as `.dockerignore`. Check the docker.yml workflow `context` value.

Actually, looking at docker.yml line: `context: .` -- the build context is the repo root. So create `.dockerignore` at the repo root.

- [ ] **Step 2: Commit**

```bash
git add .dockerignore
git commit -m "chore: add .dockerignore to reduce Docker build context

Excludes .git, tests, docs, planning, and cache directories from
the Docker build context for faster builds."
```

---

### Task 2: Rewrite Dockerfile for slim, fast, secure builds

**Files:**
- Modify: `docker/Dockerfile`

- [ ] **Step 1: Read the current Dockerfile**

Current `docker/Dockerfile` (16 lines):
```dockerfile
FROM condaforge/mambaforge:latest AS base
COPY conda/environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -afy
RUN mamba install -n muconespan-tools -c bioconda clair3 && mamba clean -afy
COPY . /app
WORKDIR /app
RUN mamba run -n muconespan-tools pip install .
SHELL ["mamba", "run", "-n", "muconespan-tools", "/bin/bash", "-c"]
ENTRYPOINT ["mamba", "run", "-n", "muconespan-tools", "muconespan"]
CMD ["--help"]
```

- [ ] **Step 2: Rewrite with multi-stage build**

Replace `docker/Dockerfile` with:

```dockerfile
# ============================================================================
# muconespan -- Multi-stage Docker build
# Stage 1: Install conda env + pip package
# Stage 2: Copy resolved env into slim runtime image
# ============================================================================

# Pin base image for reproducibility
FROM mambaorg/micromamba:1.5.10-bookworm-slim AS builder

# Install conda environment first (layer cached unless environment.yml changes)
COPY --chown=$MAMBA_USER:$MAMBA_USER conda/environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm -f /tmp/environment.yml

# Install Clair3 into the same environment
RUN micromamba install -y -n muconespan-tools -c bioconda -c conda-forge clair3 && \
    micromamba clean --all --yes

# Install muconespan Python package
# Copy pyproject.toml first for layer caching (source changes don't bust dep cache)
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /app/
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ /app/src/
COPY --chown=$MAMBA_USER:$MAMBA_USER data/ /app/data/

WORKDIR /app
RUN micromamba run -n muconespan-tools pip install --no-cache-dir --no-deps . && \
    rm -rf /app

# ============================================================================
# Runtime stage -- slim image with only the resolved environment
# ============================================================================
FROM mambaorg/micromamba:1.5.10-bookworm-slim

# Copy the fully resolved conda environment from builder
COPY --from=builder /opt/conda/envs/muconespan-tools /opt/conda/envs/muconespan-tools

# Activate environment via PATH (faster than micromamba run wrapper)
ENV PATH="/opt/conda/envs/muconespan-tools/bin:${PATH}"

# Run as non-root user
USER $MAMBA_USER
WORKDIR /data

# Metadata
LABEL org.opencontainers.image.source="https://github.com/berntpopp/MucOneSpan" \
      org.opencontainers.image.description="muconespan: MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data" \
      org.opencontainers.image.licenses="MIT"

ENTRYPOINT ["muconespan"]
CMD ["--help"]
```

- [ ] **Step 3: Verify Dockerfile builds locally (optional -- requires Docker)**

Run: `docker build -f docker/Dockerfile -t muconespan:test .`

If Docker is not available locally, this will be verified in CI.

- [ ] **Step 4: Commit**

```bash
git add docker/Dockerfile
git commit -m "feat: rewrite Dockerfile for slim, fast, secure builds

- Switch from mambaforge:latest to micromamba:1.5.10-bookworm-slim
- Multi-stage build: builder installs deps, runtime copies only env
- Layer caching: pyproject.toml copied before source code
- Non-root user via \$MAMBA_USER
- ENV PATH activation instead of mamba run wrapper
- OCI metadata labels"
```

---

### Task 3: Update Docker CI workflow with BuildKit caching

**Files:**
- Modify: `.github/workflows/docker.yml`

- [ ] **Step 1: Update workflow with BuildKit cache**

Replace `.github/workflows/docker.yml` with:

```yaml
name: Docker

on:
  push:
    branches: [main]
    tags: ["v*"]

env:
  FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: true

jobs:
  build-and-push:
    name: Build & Push Docker Image
    runs-on: ubuntu-latest
    permissions:
      contents: read
      packages: write

    steps:
      - uses: actions/checkout@v6

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Log in to GHCR
        uses: docker/login-action@v4
        with:
          registry: ghcr.io
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v6
        with:
          images: ghcr.io/${{ github.repository }}
          tags: |
            type=ref,event=branch
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=sha

      - name: Build and push
        uses: docker/build-push-action@v7
        with:
          context: .
          file: docker/Dockerfile
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha
          cache-to: type=gha,mode=max
```

- [ ] **Step 2: Commit**

```bash
git add .github/workflows/docker.yml
git commit -m "feat: add BuildKit layer caching to Docker CI workflow

Uses GitHub Actions cache (type=gha) for faster rebuilds. Adds
docker/setup-buildx-action for BuildKit support."
```

---

### Task 4: Pin conda environment versions

**Files:**
- Modify: `conda/environment.yml`
- Create: `conda/environment-dev.yml`

- [ ] **Step 1: Pin exact versions in environment.yml**

Replace `conda/environment.yml`:

```yaml
name: muconespan-tools
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.10
  - minimap2=2.28
  - samtools=1.21
  - bcftools=1.21
  - htslib=1.21
```

- [ ] **Step 2: Create environment-dev.yml with loose ranges**

Create `conda/environment-dev.yml`:

```yaml
# Development environment with loose version ranges.
# Use this for local development where you want the latest compatible tools.
# For reproducible builds (CI, Docker, publications), use environment.yml.
name: muconespan-tools
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.10
  - minimap2>=2.28
  - samtools>=1.21
  - bcftools>=1.10
  - htslib
```

- [ ] **Step 3: Commit**

```bash
git add conda/environment.yml conda/environment-dev.yml
git commit -m "feat: pin conda tool versions for reproducibility

Pin minimap2=2.28, samtools=1.21, bcftools=1.21, htslib=1.21 in
environment.yml. Add environment-dev.yml with loose ranges for
development flexibility."
```

---

### Task 5: Add Singularity/Apptainer definition

**Files:**
- Create: `docker/muconespan.def`

- [ ] **Step 1: Create Singularity definition**

Create `docker/muconespan.def`:

```singularity
Bootstrap: docker
From: ghcr.io/berntpopp/muconespan:latest

%labels
    Author Bernt Popp
    Version 0.3.0
    Description muconespan: MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data

%help
    muconespan: MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data.

    Usage:
        singularity run muconespan.sif --help
        singularity run muconespan.sif run --input reads.bam --output-dir results/

    Build:
        singularity build muconespan.sif docker/muconespan.def

%runscript
    exec muconespan "$@"
```

- [ ] **Step 2: Commit**

```bash
git add docker/muconespan.def
git commit -m "feat: add Singularity/Apptainer definition for HPC environments

Bootstraps from the GHCR Docker image. Critical for HPC clusters
where Docker is typically unavailable."
```

---

### Task 6: Add GitHub Release workflow

**Files:**
- Create: `.github/workflows/release.yml`

- [ ] **Step 1: Create release workflow**

Create `.github/workflows/release.yml`:

```yaml
name: Release

on:
  push:
    tags: ["v*"]

env:
  FORCE_JAVASCRIPT_ACTIONS_TO_NODE24: true

permissions:
  contents: write

jobs:
  release:
    name: Create GitHub Release
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v6
        with:
          fetch-depth: 0

      - name: Create release
        uses: softprops/action-gh-release@v2
        with:
          generate_release_notes: true
```

- [ ] **Step 2: Commit**

```bash
git add .github/workflows/release.yml
git commit -m "feat: add GitHub Release automation on version tags

Triggers on v* tags and auto-generates release notes from commit
history using softprops/action-gh-release."
```

---

### Task 7: Record tool versions in pipeline output

**Files:**
- Modify: `src/muc_one_span/tools.py`
- Modify: `src/muc_one_span/cli.py`
- Test: `tests/unit/test_tools.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/unit/test_tools.py`:

```python
def test_get_tool_versions_returns_versions(mocker):
    """get_tool_versions captures version strings from external tools."""
    mocker.patch(
        "muc_one_span.tools.subprocess.run",
        return_value=mocker.MagicMock(returncode=0, stdout="minimap2 2.28-r1209\n"),
    )
    mocker.patch("muc_one_span.tools.shutil.which", return_value="/usr/bin/minimap2")

    from muc_one_span.tools import get_tool_versions

    versions = get_tool_versions(["minimap2"])
    assert "minimap2" in versions
    assert "2.28" in versions["minimap2"]


def test_get_tool_versions_handles_missing_tool(mocker):
    """get_tool_versions returns 'not found' for missing tools."""
    mocker.patch("muc_one_span.tools.shutil.which", return_value=None)

    from muc_one_span.tools import get_tool_versions

    versions = get_tool_versions(["nonexistent_tool"])
    assert versions["nonexistent_tool"] == "not found"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_tools.py::test_get_tool_versions_returns_versions tests/unit/test_tools.py::test_get_tool_versions_handles_missing_tool -v --no-cov`
Expected: FAIL (ImportError or AttributeError -- function doesn't exist)

- [ ] **Step 3: Add get_tool_versions to tools.py**

Add at the end of `src/muc_one_span/tools.py` (after `check_tools()`):

```python
def get_tool_versions(tools: list[str]) -> dict[str, str]:
    """Capture installed tool versions for reproducibility metadata.

    Attempts to run ``<tool> --version`` for each tool and captures the
    first line of output as the version string.

    Args:
        tools: List of tool names to query.

    Returns:
        Dictionary mapping tool name to version string, or "not found"
        if the tool is not available.
    """
    versions: dict[str, str] = {}
    for tool in tools:
        if not shutil.which(tool):
            versions[tool] = "not found"
            continue
        try:
            result = subprocess.run(
                [tool, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            first_line = result.stdout.strip().splitlines()[0] if result.stdout.strip() else ""
            if not first_line and result.stderr.strip():
                first_line = result.stderr.strip().splitlines()[0]
            versions[tool] = first_line or "unknown"
        except (subprocess.TimeoutExpired, FileNotFoundError, IndexError):
            versions[tool] = "unknown"
    return versions
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_tools.py::test_get_tool_versions_returns_versions tests/unit/test_tools.py::test_get_tool_versions_handles_missing_tool -v --no-cov`
Expected: PASS

- [ ] **Step 5: Wire tool versions into the run subcommand**

In `src/muc_one_span/cli.py`, in the `run()` function, after the `check_tools()` call (line 340), add:

```python
    from muc_one_span.tools import get_tool_versions

    tool_versions = get_tool_versions(["minimap2", "samtools", "bcftools", "run_clair3.sh"])
```

Then in the summary dict construction (around line 403), add tool_versions:

```python
    # Summary
    summary = {
        "alleles": alleles_result,
        "classifications": {
            k: {
                "structure": v["structure"],
                "mutations": v["mutations_detected"],
            }
            for k, v in all_results.items()
        },
        "tool_versions": tool_versions,
        "pipeline_version": __version__,
    }
```

- [ ] **Step 6: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/tools.py src/muc_one_span/cli.py tests/unit/test_tools.py
git commit -m "feat: record tool versions in pipeline output

Add get_tool_versions() to tools.py that captures version strings
from external tools. Include tool_versions and pipeline_version
in summary.json for reproducibility metadata."
```

---

### Task 8: Add validation script

**Files:**
- Create: `scripts/validate_catalog.py`

- [ ] **Step 1: Create the validation script**

Create `scripts/validate_catalog.py`:

```python
#!/usr/bin/env python3
"""Validate repeat classification against MucOneUp test samples.

Runs classification on all generated test samples and compares
observed structure strings against expected results. Produces a
TSV validation report.

Usage:
    python scripts/validate_catalog.py [--data-dir tests/data/generated]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate repeat classification catalog")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tests/data/generated"),
        help="Directory containing generated test samples",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("tests/results"),
        help="Directory containing batch_results.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV path (default: stdout)",
    )
    args = parser.parse_args()

    results_path = args.results_dir / "batch_results.json"
    if not results_path.exists():
        print(f"Error: {results_path} not found. Run batch_analyze.py first.", file=sys.stderr)
        sys.exit(1)

    results = json.loads(results_path.read_text())

    # Header
    header = [
        "sample",
        "expected_mutation",
        "status",
        "detected_mutations",
        "allele_1_len",
        "allele_2_len",
    ]

    rows = []
    for entry in results:
        detected = ", ".join(
            m.get("mutation_name", "unknown")
            for m in entry.get("details", {}).get("detected_mutations", [])
        )
        rows.append([
            entry["sample"],
            entry.get("expected_mutation", "none"),
            entry["status"],
            detected or "none",
            str(entry.get("details", {}).get("allele_1_structure_len", "")),
            str(entry.get("details", {}).get("allele_2_structure_len", "")),
        ])

    # Summary stats
    tp = sum(1 for r in results if r["status"] == "TP")
    tn = sum(1 for r in results if r["status"] == "TN")
    fp = sum(1 for r in results if r["status"] == "FP")
    fn = sum(1 for r in results if r["status"] == "FN")
    total = len(results)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

    out = args.output.open("w") if args.output else sys.stdout
    try:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
        out.write(f"\n# Summary: {total} samples, {tp} TP, {tn} TN, {fp} FP, {fn} FN\n")
        out.write(f"# Sensitivity: {sensitivity:.1%}, Specificity: {specificity:.1%}\n")
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add scripts/validate_catalog.py
git commit -m "feat: add validation script for repeat classification catalog

Reads batch_results.json and produces a TSV validation report with
sensitivity/specificity summary. Used for scientific validation."
```

---

### Task 9: Add known limitations documentation

**Files:**
- Create: `docs/reference/limitations.md`

- [ ] **Step 1: Create the limitations page**

Create `docs/reference/limitations.md`:

```markdown
# Known Limitations

## Allele Length Detection

### Homozygous Same-Length Alleles

When both alleles have the same repeat count (e.g., 60/60), the indel-valley splitting algorithm may incorrectly separate reads into two clusters. The pipeline reports `same_length: true` and uses a disambiguation strategy, but accuracy is reduced compared to heterozygous samples with distinct allele lengths.

### Extreme PCR Bias (Asymmetric Alleles)

For highly asymmetric allele pairs (e.g., 25/140), the shorter allele amplifies much more efficiently during PCR. The longer allele may have very low read coverage (<10x), making allele detection unreliable. The pipeline requires a minimum coverage threshold (default: 10 reads) per allele.

## Variant Calling

### Long VNTR Alleles (>100 repeats)

Clair3 has reduced sensitivity for detecting 1bp insertions within very long tandem repeat regions (>6kb). For alleles with >100 repeat units, single-nucleotide frameshifts (e.g., dupC) may not be detected by Clair3, leading to false negatives.

**Affected test samples:**
- `sample_dupc_100_120`: dupC mutation in 100-repeat allele not detected
- `sample_long_120_140`: mutations in 120+ repeat alleles not detected

### Boundary Repeat Mutations

Mutations detected in the last 3 repeat units of an allele receive a boundary penalty (0.5x confidence multiplier) because these positions are prone to alignment artifacts at the VNTR boundary.

## Classification

### Novel Repeat Types

Repeat sequences with >2 base substitutions from any known type are classified as `novel_repeat`. These may represent true biological variation or alignment/sequencing errors. The edit distance and identity percentage are reported for manual review.

### Bidirectional Classification Fallback

When the forward classification pass encounters a large indel (>30bp offset from the expected 60bp unit length), the algorithm falls back to bidirectional classification. This backward anchoring from the 3' end recovers downstream repeats but may produce different results than purely forward classification for the ambiguous region.
```

- [ ] **Step 2: Add to mkdocs.yml navigation**

Add the limitations page to the `nav:` section in `mkdocs.yml` under the Reference section.

- [ ] **Step 3: Commit**

```bash
git add docs/reference/limitations.md mkdocs.yml
git commit -m "docs: add known limitations page

Documents allele detection edge cases, Clair3 sensitivity limits
for long VNTRs, boundary penalties, and classification fallbacks."
```

---

### Task 10: Run quality checks

- [ ] **Step 1: Run lint check**

Run: `uv run ruff check src/muc_one_span/`
Expected: No new errors

- [ ] **Step 2: Run type check**

Run: `uv run mypy src/muc_one_span/`
Expected: No new errors

- [ ] **Step 3: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass

- [ ] **Step 4: Fix any issues and commit**

```bash
git add -u
git commit -m "fix: address lint and type check issues from Phase 2 changes"
```
