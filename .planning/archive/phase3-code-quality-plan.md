# Phase 3: Code Quality -- Type Safety + Test Coverage

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add TypedDict return types for classification and allele results, decompose the 216-line `classify_sequence()` into focused helpers, extract VCF parsing into its own module, and harden test coverage from 70% to 80%+.

**Architecture:** Type safety changes (Tasks 1-4) are done first since tests written in Tasks 5-8 reference the new types. VCF extraction (Task 4) creates a new `vcf.py` module and updates imports in `calling.py` and `cli.py`.

**Tech Stack:** Python `typing.TypedDict`, pytest, pytest-mock, Click CliRunner

---

### Task 1: Add TypedDict for classification results

**Files:**
- Modify: `src/muc_one_span/classify.py:16-18`

- [ ] **Step 1: Add TypedDict definitions**

In `src/muc_one_span/classify.py`, after the existing imports (line 18), add:

```python
from typing import TypedDict


class RepeatDifference(TypedDict):
    """A single difference between a repeat sequence and its closest reference."""

    pos: int
    ref: str
    alt: str
    type: str


class RepeatClassification(TypedDict, total=False):
    """Classification result for a single repeat unit."""

    type: str
    match: str  # "exact" or "closest"
    confidence: float
    closest_match: str
    edit_distance: int | float
    identity_pct: float
    differences: list[RepeatDifference]
    classification: str  # "mutation", "novel_repeat", or "variant"
    frameshift: bool
    mutation_name: str
    parent_repeat: str
    index: int


class MutationDetected(TypedDict, total=False):
    """A mutation detected during classification."""

    repeat_index: int
    closest_type: str
    mutation_name: str
    template_match: bool
    frameshift: bool
    differences: list[RepeatDifference]
    vcf_support: bool
    vcf_qual: float
    boundary: bool


class SequenceClassification(TypedDict):
    """Classification result for a full consensus sequence."""

    structure: str
    repeats: list[RepeatClassification]
    mutations_detected: list[MutationDetected]
    cumulative_offset: int
    allele_confidence: float
    exact_match_pct: float
```

- [ ] **Step 2: Update classify_repeat return type annotation**

Change `classify_repeat` signature (line 148):

```python
# Before:
def classify_repeat(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:

# After:
def classify_repeat(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> RepeatClassification:
```

- [ ] **Step 3: Update classify_sequence return type annotation**

Change `classify_sequence` signature (line 320):

```python
# Before:
def classify_sequence(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:

# After:
def classify_sequence(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> SequenceClassification:
```

- [ ] **Step 4: Run tests and type check**

Run: `uv run pytest tests/unit/test_classify.py --no-cov -q`
Expected: All pass (TypedDict is structurally compatible with dict)

Run: `uv run mypy src/muc_one_span/classify.py`
Expected: Pass or minor issues to fix

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/classify.py
git commit -m "feat: add TypedDict types for classification results

Add RepeatClassification, MutationDetected, SequenceClassification,
and RepeatDifference TypedDicts. Update classify_repeat and
classify_sequence return type annotations."
```

---

### Task 2: Add TypedDict for allele results

**Files:**
- Modify: `src/muc_one_span/alleles.py`

- [ ] **Step 1: Add TypedDict definitions**

In `src/muc_one_span/alleles.py`, after the existing imports (around line 15), add:

```python
from typing import TypedDict


class AlleleInfo(TypedDict):
    """Information about a single detected allele."""

    length: int
    reads: int
    canonical_repeats: int
    contig_name: str
    cluster_contigs: list[str]


class AlleleResult(TypedDict):
    """Result of allele detection for a sample."""

    allele_1: AlleleInfo
    allele_2: AlleleInfo
    homozygous: bool
    same_length: bool
```

- [ ] **Step 2: Update function return type annotations**

Change `_build_allele_info` (line 296):
```python
# Before:
def _build_allele_info(cluster: dict, best_contig: str | None = None) -> dict:
# After:
def _build_allele_info(cluster: dict, best_contig: str | None = None) -> AlleleInfo:
```

Change `detect_alleles` (line 316):
```python
# Before:
def detect_alleles(
    counts: dict[int, int],
    min_coverage: int = 10,
    bam_path: Path | None = None,
) -> dict:
# After:
def detect_alleles(
    counts: dict[int, int],
    min_coverage: int = 10,
    bam_path: Path | None = None,
) -> AlleleResult:
```

- [ ] **Step 3: Run tests and type check**

Run: `uv run pytest tests/unit/test_alleles.py --no-cov -q`
Expected: All pass

Run: `uv run mypy src/muc_one_span/alleles.py`
Expected: Pass or minor issues to fix

- [ ] **Step 4: Commit**

```bash
git add src/muc_one_span/alleles.py
git commit -m "feat: add TypedDict types for allele detection results

Add AlleleInfo and AlleleResult TypedDicts. Update _build_allele_info
and detect_alleles return type annotations."
```

---

### Task 3: Decompose classify_sequence into helpers

**Files:**
- Modify: `src/muc_one_span/classify.py:320-536`

This refactors the 216-line `classify_sequence()` into three focused helpers. The public API and return types remain identical.

- [ ] **Step 1: Extract _forward_classify helper**

Add before `classify_sequence()` (around line 318):

```python
def _forward_classify(
    sequence: str,
    repeat_dict: RepeatDictionary,
    unit_length: int,
    max_indel_probe: int,
) -> tuple[list[RepeatClassification], list[MutationDetected], list[str], int, int]:
    """Forward pass: classify repeats left-to-right with offset tracking.

    Args:
        sequence: Full consensus sequence.
        repeat_dict: The loaded repeat dictionary.
        unit_length: Expected repeat unit length (60bp).
        max_indel_probe: Maximum indel size to probe.

    Returns:
        Tuple of (repeats, mutations, labels, final_pos, cumulative_offset).
    """
    repeats: list[RepeatClassification] = []
    mutations: list[MutationDetected] = []
    labels: list[str] = []

    pos = 0
    repeat_index = 0
    cumulative_offset = 0

    while pos < len(sequence):
        repeat_index += 1

        remaining = len(sequence) - pos
        if remaining < unit_length // 2:
            break

        best_result: RepeatClassification | None = None
        best_dist: int | float = float("inf")
        best_window_size = unit_length

        # --- Phase 1: Check ALL probe sizes for exact match ---
        exact_found = False
        canonical_result: RepeatClassification | None = None

        for probe_size in _probe_sizes_generator(unit_length, max_indel_probe, remaining):
            window = sequence[pos : pos + probe_size]
            if repeat_dict.mutated_sequences and window in repeat_dict.mutated_sequences:
                parent, mname = repeat_dict.mutated_sequences[window]
                best_result = RepeatClassification(
                    type=f"{parent}:{mname}",
                    match="exact",
                    confidence=1.0,
                    mutation_name=mname,
                    parent_repeat=parent,
                )
                best_window_size = probe_size
                best_dist = 0
                exact_found = True
                break
            if window in repeat_dict.seq_to_id and canonical_result is None:
                canonical_result = RepeatClassification(
                    type=repeat_dict.seq_to_id[window],
                    match="exact",
                    confidence=1.0,
                )

        if not exact_found and canonical_result is not None:
            best_result = canonical_result
            best_window_size = unit_length
            best_dist = 0
            exact_found = True

        # --- Phase 2: Edit distance fallback ---
        if not exact_found:
            if remaining >= unit_length:
                window = sequence[pos : pos + unit_length]
                result = classify_repeat(window, repeat_dict)
                dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
                best_dist = dist
                best_result = result
                best_window_size = unit_length

            if best_dist > 0:
                for probe_size in range(
                    max(unit_length - max_indel_probe, unit_length // 2),
                    min(unit_length + max_indel_probe + 1, remaining + 1),
                ):
                    if probe_size == unit_length:
                        continue
                    window = sequence[pos : pos + probe_size]
                    result = classify_repeat(window, repeat_dict)
                    dist = 0 if result["match"] == "exact" else result.get("edit_distance", 999)
                    if dist < best_dist:
                        best_dist = dist
                        best_result = result
                        best_window_size = probe_size
                    if best_dist <= 1:
                        break

        if best_result is None:
            raise RuntimeError(
                f"Classification failed: no match found at position {pos} "
                f"(remaining: {remaining} bp, repeat index: {repeat_index})"
            )
        result = best_result
        advance = best_window_size

        if best_window_size != unit_length:
            net_indel = best_window_size - unit_length
            cumulative_offset += net_indel

        result["index"] = repeat_index

        if result["match"] == "exact":
            labels.append(result["type"])
            if result.get("mutation_name"):
                mutations.append(
                    MutationDetected(
                        repeat_index=repeat_index,
                        closest_type=result.get("parent_repeat", result["type"]),
                        mutation_name=result["mutation_name"],
                        template_match=True,
                        frameshift=True,
                    )
                )
        elif result.get("classification") == "mutation":
            labels.append(f"{result['closest_match']}m")
            mutations.append(
                MutationDetected(
                    repeat_index=repeat_index,
                    closest_type=result["closest_match"],
                    differences=result["differences"],
                    frameshift=result.get("frameshift", False),
                )
            )
        else:
            labels.append(f"?{result.get('closest_match', '?')}")

        repeats.append(result)
        pos += max(advance, 1)

    return repeats, mutations, labels, pos, cumulative_offset
```

- [ ] **Step 2: Extract _apply_bidirectional_fallback helper**

Add after `_forward_classify`:

```python
def _apply_bidirectional_fallback(
    sequence: str,
    repeat_dict: RepeatDictionary,
    repeats: list[RepeatClassification],
    mutations: list[MutationDetected],
    labels: list[str],
    forward_pos: int,
) -> tuple[list[RepeatClassification], list[MutationDetected], list[str]]:
    """Apply bidirectional fallback when forward pass leaves unconsumed sequence.

    Args:
        sequence: Full consensus sequence.
        repeat_dict: The loaded repeat dictionary.
        repeats: Repeats from forward pass (mutated in place).
        mutations: Mutations from forward pass (mutated in place).
        labels: Labels from forward pass (mutated in place).
        forward_pos: Position where forward pass stopped.

    Returns:
        Updated (repeats, mutations, labels) tuple.
    """
    unit_length = repeat_dict.repeat_length_bp

    if forward_pos >= len(sequence) - unit_length // 2:
        return repeats, mutations, labels

    backward = _classify_backward(sequence, repeat_dict, forward_pos)
    if not backward:
        return repeats, mutations, labels

    repeat_index = repeats[-1]["index"] if repeats else 0

    # Gap between forward and backward
    gap_start = forward_pos
    gap_end = backward[0][1]
    if gap_end > gap_start:
        gap_seq = sequence[gap_start:gap_end]
        gap_result = classify_repeat(gap_seq, repeat_dict)
        repeat_index += 1
        gap_result["index"] = repeat_index
        if gap_result.get("classification") == "mutation" or gap_result["match"] != "exact":
            labels.append(f"{gap_result.get('closest_match', '?')}m")
            mutations.append(
                MutationDetected(
                    repeat_index=repeat_index,
                    closest_type=gap_result.get("closest_match", "?"),
                    differences=gap_result.get("differences", []),
                    frameshift=gap_result.get("frameshift", False),
                )
            )
        else:
            labels.append(gap_result["type"])
        repeats.append(gap_result)

    # Append backward results
    for bwd_result, _bwd_start, _bwd_end in backward:
        repeat_index += 1
        bwd_result["index"] = repeat_index
        if bwd_result["match"] == "exact":
            labels.append(bwd_result["type"])
        else:
            labels.append(f"?{bwd_result.get('closest_match', '?')}")
        repeats.append(bwd_result)

    return repeats, mutations, labels
```

- [ ] **Step 3: Extract _compute_classification_summary helper**

Add after `_apply_bidirectional_fallback`:

```python
def _compute_classification_summary(
    repeats: list[RepeatClassification],
    mutations: list[MutationDetected],
    labels: list[str],
    cumulative_offset: int,
) -> SequenceClassification:
    """Compute confidence scores and build the final classification result.

    Args:
        repeats: Classified repeat units.
        mutations: Detected mutations.
        labels: Repeat type labels for structure string.
        cumulative_offset: Net indel offset from forward pass.

    Returns:
        Complete sequence classification result.
    """
    confidences = [r.get("confidence", 1.0) for r in repeats]
    exact_count = sum(1 for r in repeats if r.get("match") == "exact")
    allele_confidence = sum(confidences) / len(confidences) if confidences else 0.0
    exact_match_pct = (exact_count / len(repeats) * 100) if repeats else 0.0

    return SequenceClassification(
        structure=" ".join(labels),
        repeats=repeats,
        mutations_detected=mutations,
        cumulative_offset=cumulative_offset,
        allele_confidence=round(allele_confidence, 4),
        exact_match_pct=round(exact_match_pct, 1),
    )
```

- [ ] **Step 4: Rewrite classify_sequence to use helpers**

Replace the entire `classify_sequence` function body (lines 320-536) with:

```python
def classify_sequence(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> SequenceClassification:
    """Classify all repeat units in a consensus sequence.

    Uses offset-aware windowing: when a repeat contains an indel, the
    cumulative offset is tracked and subsequent window boundaries are
    shifted accordingly.  This corrects for frameshift propagation --
    a 1bp insertion at repeat 25 would otherwise misalign all downstream
    windows.

    Args:
        sequence: Full consensus sequence (flanking regions should be trimmed).
        repeat_dict: The loaded repeat dictionary.

    Returns:
        Classification with structure string, per-repeat details, and mutation report.
    """
    unit_length = repeat_dict.repeat_length_bp
    max_indel_probe = 30

    repeats, mutations, labels, pos, cumulative_offset = _forward_classify(
        sequence, repeat_dict, unit_length, max_indel_probe
    )

    repeats, mutations, labels = _apply_bidirectional_fallback(
        sequence, repeat_dict, repeats, mutations, labels, pos
    )

    return _compute_classification_summary(repeats, mutations, labels, cumulative_offset)
```

- [ ] **Step 5: Run tests to verify behavior is unchanged**

Run: `uv run pytest tests/unit/test_classify.py --no-cov -q`
Expected: All pass (refactor preserves exact behavior)

- [ ] **Step 6: Run full test suite and type check**

Run: `uv run pytest tests/unit/ --no-cov -q && uv run mypy src/muc_one_span/classify.py`
Expected: All pass

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/classify.py
git commit -m "refactor: decompose classify_sequence into focused helpers

Extract _forward_classify(), _apply_bidirectional_fallback(), and
_compute_classification_summary() from the 216-line classify_sequence().
No behavior change -- all existing tests pass."
```

---

### Task 4: Extract VCF parsing into vcf.py

**Files:**
- Create: `src/muc_one_span/vcf.py`
- Modify: `src/muc_one_span/calling.py`
- Modify: `src/muc_one_span/cli.py`
- Modify: `tests/unit/test_calling.py`
- Create: `tests/unit/test_vcf.py`

- [ ] **Step 1: Create vcf.py with extracted functions**

Create `src/muc_one_span/vcf.py`:

```python
"""VCF parsing and filtering utilities."""

from __future__ import annotations

from pathlib import Path

from muc_one_span.tools import run_tool


def filter_vcf(
    vcf_path: Path,
    reference_path: Path,
    output_dir: Path,
    min_qual: float = 0.0,
    min_dp: int = 0,
) -> Path:
    """Normalize and filter a VCF file with bcftools.

    Runs ``bcftools norm -f <reference>`` followed by
    ``bcftools view -f PASS`` (with optional quality filters) and indexes
    the result.

    Args:
        vcf_path: Path to input VCF (may be gzipped).
        reference_path: Path to reference FASTA for left-normalisation.
        output_dir: Directory for output files.
        min_qual: Minimum QUAL score to keep a variant (0 = no filter).
        min_dp: Minimum depth to keep a variant (0 = no filter).

    Returns:
        Path to the filtered, indexed VCF (``variants.vcf.gz``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    norm_vcf = output_dir / "normalized.vcf.gz"
    filtered = output_dir / "variants.vcf.gz"

    run_tool(
        [
            "bcftools",
            "norm",
            "-f",
            str(reference_path),
            "-o",
            str(norm_vcf),
            "-O",
            "z",
            str(vcf_path),
        ]
    )

    is_empty = not norm_vcf.exists() or norm_vcf.stat().st_size == 0

    view_cmd = [
        "bcftools",
        "view",
        "-f",
        "PASS",
    ]
    if not is_empty and min_qual > 0:
        view_cmd.extend(["-i", f"QUAL>={min_qual}"])

    view_cmd.extend(
        [
            "-o",
            str(filtered),
            "-O",
            "z",
            str(norm_vcf),
        ]
    )
    run_tool(view_cmd)
    run_tool(["bcftools", "index", str(filtered)])

    norm_vcf.unlink(missing_ok=True)

    return filtered


def parse_vcf_genotypes(vcf_path: Path) -> list[dict]:
    """Parse VCF to extract variant positions and genotypes.

    Returns list of dicts with keys: chrom, pos, ref, alt, genotype.
    """
    try:
        output = run_tool(
            [
                "bcftools",
                "query",
                "-f",
                "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n",
                str(vcf_path),
            ]
        )
    except RuntimeError:
        return []

    variants: list[dict] = []
    for line in output.strip().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        variants.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "ref": fields[2],
                "alt": fields[3],
                "genotype": fields[4],
            }
        )
    return variants


def parse_vcf_variants(vcf_path: Path) -> list[dict]:
    """Parse VCF to extract variant positions and quality scores."""
    try:
        output = run_tool(["bcftools", "query", "-f", "%POS\\t%QUAL\\n", str(vcf_path)])
    except RuntimeError:
        return []
    variants = []
    for line in output.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            try:
                variants.append({"pos": int(parts[0]), "qual": float(parts[1])})
            except ValueError:
                continue
    return variants
```

- [ ] **Step 2: Update calling.py to import from vcf.py**

In `src/muc_one_span/calling.py`:
- Remove `filter_vcf`, `parse_vcf_genotypes`, and `parse_vcf_variants` function definitions (lines 181-315)
- Add re-exports at the top for backwards compatibility:

```python
from muc_one_span.vcf import filter_vcf, parse_vcf_genotypes, parse_vcf_variants

__all__ = [
    "extract_allele_reads",
    "run_clair3",
    "filter_vcf",
    "parse_vcf_genotypes",
    "parse_vcf_variants",
    "disambiguate_same_length_alleles",
    "call_variants_per_allele",
]
```

And update internal uses of `filter_vcf` in `disambiguate_same_length_alleles` and `call_variants_per_allele` -- these already call `filter_vcf()` directly, which will now resolve via the import.

- [ ] **Step 3: Update cli.py import**

In `src/muc_one_span/cli.py` line 333, update:

```python
# Before:
from muc_one_span.calling import call_variants_per_allele, parse_vcf_variants

# After:
from muc_one_span.calling import call_variants_per_allele
from muc_one_span.vcf import parse_vcf_variants
```

- [ ] **Step 4: Create test_vcf.py by moving relevant tests**

Create `tests/unit/test_vcf.py` with tests for the three extracted functions. Copy the relevant test classes from `test_calling.py`:
- `TestFilterVcf`
- `TestFilterVcfEmptyVcf`
- `TestFilterVcfQuality`
- `TestParseVcfGenotypes`

Update the patch targets from `muc_one_span.calling.run_tool` to `muc_one_span.vcf.run_tool`.

- [ ] **Step 5: Update test_calling.py**

Remove the test classes that were moved to `test_vcf.py`. Keep tests for:
- `TestExtractAlleleReads`
- `TestRunClair3`
- `TestCallVariantsPerAllele`
- `TestExtractAndRemapReads`
- `TestDisambiguateSameLengthAlleles`

- [ ] **Step 6: Run full test suite**

Run: `uv run pytest tests/unit/ --no-cov -q`
Expected: All pass

Run: `uv run mypy src/muc_one_span/`
Expected: Pass

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/vcf.py src/muc_one_span/calling.py src/muc_one_span/cli.py tests/unit/test_vcf.py tests/unit/test_calling.py
git commit -m "refactor: extract VCF parsing into dedicated vcf.py module

Move filter_vcf, parse_vcf_genotypes, parse_vcf_variants from
calling.py to vcf.py. Re-export from calling.py for compatibility.
Tests split correspondingly."
```

---

### Task 5: Add mapping.py tests (50% -> 90%)

**Files:**
- Modify: `tests/unit/test_mapping.py`

- [ ] **Step 1: Write tests for _run_mapping_pipeline**

Add to `tests/unit/test_mapping.py`:

```python
from unittest.mock import MagicMock, patch


class TestRunMappingPipeline:
    """Tests for _run_mapping_pipeline subprocess handling."""

    def test_minimap2_not_found(self):
        """Raises FileNotFoundError when minimap2 is missing."""
        from muc_one_span.mapping import _run_mapping_pipeline

        with patch("muc_one_span.mapping.subprocess.Popen", side_effect=FileNotFoundError):
            with pytest.raises(FileNotFoundError, match="minimap2"):
                _run_mapping_pipeline(
                    input_path=Path("/tmp/test.fastq"),
                    reference_path=Path("/tmp/ref.fa"),
                    bam_path=Path("/tmp/out.bam"),
                    threads=1,
                )

    def test_samtools_not_found(self):
        """Raises FileNotFoundError when samtools is missing."""
        from muc_one_span.mapping import _run_mapping_pipeline

        mock_p1 = MagicMock()
        mock_p1.stdout = MagicMock()

        def popen_side_effect(*args, **kwargs):
            if "minimap2" in args[0]:
                return mock_p1
            raise FileNotFoundError

        with patch("muc_one_span.mapping.subprocess.Popen", side_effect=popen_side_effect):
            with pytest.raises(FileNotFoundError, match="samtools"):
                _run_mapping_pipeline(
                    input_path=Path("/tmp/test.fastq"),
                    reference_path=Path("/tmp/ref.fa"),
                    bam_path=Path("/tmp/out.bam"),
                    threads=1,
                )

    def test_minimap2_nonzero_exit(self):
        """Raises RuntimeError when minimap2 exits with non-zero code."""
        from muc_one_span.mapping import _run_mapping_pipeline

        mock_p1 = MagicMock()
        mock_p1.stdout = MagicMock()
        mock_p1.stderr = MagicMock()
        mock_p1.stderr.read.return_value = b"error output"
        mock_p1.returncode = 1
        mock_p1.wait.return_value = 1

        mock_p2 = MagicMock()
        mock_p2.communicate.return_value = (b"", b"")
        mock_p2.returncode = 0

        with patch(
            "muc_one_span.mapping.subprocess.Popen",
            side_effect=[mock_p1, mock_p2],
        ):
            with pytest.raises(RuntimeError, match="minimap2 failed"):
                _run_mapping_pipeline(
                    input_path=Path("/tmp/test.fastq"),
                    reference_path=Path("/tmp/ref.fa"),
                    bam_path=Path("/tmp/out.bam"),
                    threads=1,
                )

    def test_samtools_nonzero_exit(self):
        """Raises RuntimeError when samtools exits with non-zero code."""
        from muc_one_span.mapping import _run_mapping_pipeline

        mock_p1 = MagicMock()
        mock_p1.stdout = MagicMock()
        mock_p1.stderr = MagicMock()
        mock_p1.stderr.read.return_value = b""
        mock_p1.returncode = 0
        mock_p1.wait.return_value = 0

        mock_p2 = MagicMock()
        mock_p2.communicate.return_value = (b"", b"samtools error")
        mock_p2.returncode = 1

        with patch(
            "muc_one_span.mapping.subprocess.Popen",
            side_effect=[mock_p1, mock_p2],
        ):
            with pytest.raises(RuntimeError, match="samtools sort failed"):
                _run_mapping_pipeline(
                    input_path=Path("/tmp/test.fastq"),
                    reference_path=Path("/tmp/ref.fa"),
                    bam_path=Path("/tmp/out.bam"),
                    threads=1,
                )

    def test_successful_pipeline(self, tmp_path):
        """Successful pipeline completes without errors."""
        from muc_one_span.mapping import _run_mapping_pipeline

        mock_p1 = MagicMock()
        mock_p1.stdout = MagicMock()
        mock_p1.stderr = MagicMock()
        mock_p1.stderr.read.return_value = b""
        mock_p1.returncode = 0
        mock_p1.wait.return_value = 0

        mock_p2 = MagicMock()
        mock_p2.communicate.return_value = (b"", b"")
        mock_p2.returncode = 0

        with patch(
            "muc_one_span.mapping.subprocess.Popen",
            side_effect=[mock_p1, mock_p2],
        ):
            _run_mapping_pipeline(
                input_path=tmp_path / "test.fastq",
                reference_path=tmp_path / "ref.fa",
                bam_path=tmp_path / "out.bam",
                threads=4,
            )

        mock_p1.stdout.close.assert_called_once()
```

- [ ] **Step 2: Run new tests**

Run: `uv run pytest tests/unit/test_mapping.py::TestRunMappingPipeline -v --no-cov`
Expected: All pass

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_mapping.py
git commit -m "test: add comprehensive tests for _run_mapping_pipeline

Test minimap2/samtools FileNotFoundError, non-zero exit codes,
SIGPIPE handling, and successful completion. Targets 90% coverage
for mapping.py."
```

---

### Task 6: Add cli.py run subcommand test (64% -> 85%)

**Files:**
- Modify: `tests/unit/test_cli.py`

- [ ] **Step 1: Write run subcommand test**

Add to `tests/unit/test_cli.py`:

```python
class TestRunSubcommand:
    """Tests for the full pipeline run subcommand."""

    def test_run_full_pipeline(self, tmp_path):
        """run subcommand orchestrates all pipeline stages."""
        from unittest.mock import patch

        runner = CliRunner()
        input_file = tmp_path / "reads.fastq"
        input_file.write_text("@read1\nACGT\n+\nIIII\n")
        output_dir = tmp_path / "results"

        mock_alleles = {
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

        mock_classification = {
            "structure": "1 2 3 4 5 X X X 6 7 8 9",
            "repeats": [{"type": "X", "match": "exact", "confidence": 1.0, "index": 1}],
            "mutations_detected": [],
            "cumulative_offset": 0,
            "allele_confidence": 1.0,
            "exact_match_pct": 100.0,
        }

        consensus_fa = tmp_path / "consensus_a1.fa"
        consensus_fa.write_text(">allele_1\nACGT\n")

        with (
            patch("muc_one_span.cli.check_tools"),
            patch("muc_one_span.cli.get_tool_versions", return_value={}),
            patch("muc_one_span.cli.map_reads", return_value=tmp_path / "mapped.bam"),
            patch("muc_one_span.cli.get_idxstats", return_value="contig_41\t100\t50\t0\n"),
            patch("muc_one_span.cli.parse_idxstats", return_value={41: 100, 51: 80}),
            patch("muc_one_span.cli.detect_alleles", return_value=mock_alleles),
            patch(
                "muc_one_span.cli.call_variants_per_allele",
                return_value={"allele_1": tmp_path / "a1.vcf.gz"},
            ),
            patch(
                "muc_one_span.cli.build_consensus_per_allele",
                return_value={"allele_1": consensus_fa},
            ),
            patch("muc_one_span.cli.classify_sequence", return_value=mock_classification),
            patch("muc_one_span.cli.validate_mutations_against_vcf", return_value=mock_classification),
            patch("muc_one_span.cli.parse_vcf_variants", return_value=[]),
            patch("muc_one_span.cli._bundled_reference", return_value=tmp_path / "ref.fa"),
        ):
            result = runner.invoke(
                main,
                ["run", "--input", str(input_file), "--output-dir", str(output_dir)],
            )

        assert result.exit_code == 0, result.output
        assert "Pipeline complete" in result.output
```

**Note:** The exact patch targets depend on where the `run()` function imports from. Since `run()` uses lazy imports inside the function body, patches need to target the `muc_one_span.cli` namespace. The executor should verify and adjust patch targets based on the actual import structure at execution time.

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/unit/test_cli.py::TestRunSubcommand -v --no-cov`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_cli.py
git commit -m "test: add run subcommand integration test

Tests full pipeline orchestration with all stages mocked. Verifies
output files are written and pipeline completes successfully."
```

---

### Task 7: Add classify.py tests for backward classification (67% -> 85%)

**Files:**
- Modify: `tests/unit/test_classify.py`

- [ ] **Step 1: Write tests for _forward_classify and _apply_bidirectional_fallback**

Add to `tests/unit/test_classify.py`:

```python
class TestForwardClassify:
    """Tests for the _forward_classify helper."""

    def test_single_exact_repeat(self, repeat_dict):
        """Forward pass classifies a single known repeat."""
        from muc_one_span.classify import _forward_classify

        # Get a known repeat sequence
        first_id = list(repeat_dict.repeats.keys())[0]
        seq = repeat_dict.repeats[first_id]

        repeats, mutations, labels, pos, offset = _forward_classify(
            seq, repeat_dict, repeat_dict.repeat_length_bp, 30
        )

        assert len(repeats) == 1
        assert repeats[0]["match"] == "exact"
        assert repeats[0]["type"] == first_id
        assert len(mutations) == 0
        assert offset == 0

    def test_two_exact_repeats(self, repeat_dict):
        """Forward pass classifies two concatenated known repeats."""
        from muc_one_span.classify import _forward_classify

        ids = list(repeat_dict.repeats.keys())[:2]
        seq = repeat_dict.repeats[ids[0]] + repeat_dict.repeats[ids[1]]

        repeats, mutations, labels, pos, offset = _forward_classify(
            seq, repeat_dict, repeat_dict.repeat_length_bp, 30
        )

        assert len(repeats) == 2
        assert labels[0] == ids[0]
        assert labels[1] == ids[1]


class TestApplyBidirectionalFallback:
    """Tests for _apply_bidirectional_fallback."""

    def test_no_fallback_when_fully_consumed(self, repeat_dict):
        """No fallback when forward pass consumed all sequence."""
        from muc_one_span.classify import _apply_bidirectional_fallback

        seq = repeat_dict.repeats[list(repeat_dict.repeats.keys())[0]]
        repeats = [{"type": "X", "match": "exact", "confidence": 1.0, "index": 1}]
        mutations = []
        labels = ["X"]

        result_repeats, result_mutations, result_labels = _apply_bidirectional_fallback(
            seq, repeat_dict, repeats, mutations, labels, len(seq)
        )

        assert len(result_repeats) == 1  # no additional repeats added


class TestComputeClassificationSummary:
    """Tests for _compute_classification_summary."""

    def test_empty_repeats(self):
        """Summary handles empty repeat list."""
        from muc_one_span.classify import _compute_classification_summary

        result = _compute_classification_summary([], [], [], 0)
        assert result["structure"] == ""
        assert result["allele_confidence"] == 0.0
        assert result["exact_match_pct"] == 0.0

    def test_all_exact_matches(self):
        """Summary computes 100% confidence for all exact matches."""
        from muc_one_span.classify import _compute_classification_summary

        repeats = [
            {"type": "X", "match": "exact", "confidence": 1.0, "index": 1},
            {"type": "A", "match": "exact", "confidence": 1.0, "index": 2},
        ]
        labels = ["X", "A"]

        result = _compute_classification_summary(repeats, [], labels, 0)
        assert result["allele_confidence"] == 1.0
        assert result["exact_match_pct"] == 100.0
        assert result["structure"] == "X A"
```

- [ ] **Step 2: Run new tests**

Run: `uv run pytest tests/unit/test_classify.py::TestForwardClassify tests/unit/test_classify.py::TestApplyBidirectionalFallback tests/unit/test_classify.py::TestComputeClassificationSummary -v --no-cov`
Expected: All pass

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_classify.py
git commit -m "test: add tests for classify_sequence decomposed helpers

Test _forward_classify, _apply_bidirectional_fallback, and
_compute_classification_summary independently."
```

---

### Task 8: Add alleles.py tests for indel valley splitting (68% -> 85%)

**Files:**
- Modify: `tests/unit/test_alleles.py`

- [ ] **Step 1: Write tests for _split_cluster_by_indel**

Add to `tests/unit/test_alleles.py`:

```python
class TestSplitClusterByIndel:
    """Tests for _split_cluster_by_indel valley splitting."""

    def test_returns_none_when_no_bam(self):
        """Returns None when BAM has no reads for cluster contigs."""
        from muc_one_span.alleles import _split_cluster_by_indel

        cluster = {
            "center": 50,
            "total_reads": 100,
            "contigs": [(50, 100)],
        }

        with patch("muc_one_span.alleles.run_tool", return_value=""):
            result = _split_cluster_by_indel(Path("/tmp/test.bam"), cluster)

        assert result is None

    def test_returns_none_when_fewer_than_two_valleys(self):
        """Returns None when reads don't form two distinct indel clusters."""
        from muc_one_span.alleles import _split_cluster_by_indel

        cluster = {
            "center": 50,
            "total_reads": 100,
            "contigs": [(50, 100)],
        }

        # All reads have similar indel lengths -- no valley
        sam_lines = "\n".join(
            _make_sam_line(f"read{i}", "contig_50", cigar="3600M", as_score=3000)
            for i in range(20)
        )

        with patch("muc_one_span.alleles.run_tool", return_value=sam_lines):
            result = _split_cluster_by_indel(Path("/tmp/test.bam"), cluster)

        # With all-0 indel lengths, there's only one valley -> None
        assert result is None

    def test_splits_when_two_clear_valleys(self):
        """Returns two sub-clusters when reads form distinct indel groups."""
        from muc_one_span.alleles import _split_cluster_by_indel

        cluster = {
            "center": 50,
            "total_reads": 40,
            "contigs": [(49, 10), (50, 20), (51, 10)],
        }

        # Group A: reads with 0bp indels (matching contig)
        # Group B: reads with ~180bp indels (3-repeat mismatch)
        sam_group_a = "\n".join(
            _make_sam_line(f"readA{i}", "contig_50", cigar="3600M", as_score=3500)
            for i in range(20)
        )
        sam_group_b = "\n".join(
            _make_sam_line(f"readB{i}", "contig_50", cigar="3420M180I", as_score=2500)
            for i in range(20)
        )

        with patch("muc_one_span.alleles.run_tool", return_value=f"{sam_group_a}\n{sam_group_b}"):
            result = _split_cluster_by_indel(Path("/tmp/test.bam"), cluster)

        if result is not None:
            assert len(result) == 2
            assert all("center" in sub for sub in result)
            assert all("total_reads" in sub for sub in result)
```

- [ ] **Step 2: Run new tests**

Run: `uv run pytest tests/unit/test_alleles.py::TestSplitClusterByIndel -v --no-cov`
Expected: All pass (may need adjustment based on exact _split_cluster_by_indel behavior)

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_alleles.py
git commit -m "test: add tests for _split_cluster_by_indel valley splitting

Test no-reads, single-valley, and two-valley scenarios for the
indel-based allele disambiguation."
```

---

### Task 9: Raise CI coverage threshold

**Files:**
- Modify: `.github/workflows/test.yml:68`

- [ ] **Step 1: Verify current coverage meets new threshold**

Run: `uv run pytest tests/unit/ --cov=muc_one_span --cov-report=term-missing`
Expected: Coverage >= 80%

- [ ] **Step 2: Update threshold**

In `.github/workflows/test.yml`, line 68:

```yaml
# Before:
        run: pytest tests/unit/ --cov=muc_one_span --cov-report=xml --cov-report=term-missing --cov-fail-under=70

# After:
        run: pytest tests/unit/ --cov=muc_one_span --cov-report=xml --cov-report=term-missing --cov-fail-under=80
```

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/test.yml
git commit -m "chore: raise CI coverage threshold from 70% to 80%

Coverage hardening in Phase 3 brings actual coverage above 80%."
```

---

### Task 10: Run quality checks

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
git commit -m "fix: address lint and type check issues from Phase 3 changes"
```
