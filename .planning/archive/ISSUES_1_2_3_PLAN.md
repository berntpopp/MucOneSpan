# Issues #1, #2, #3 Fix with Confidence Scoring -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix three open issues (homozygous allele missed mutations, large insertion misclassification, false positive mutations) and add per-repeat confidence scoring with a known mutation catalog.

**Architecture:** Add a mutation catalog (from MucOneUp) to the repeat dictionary, enabling exact template matching for known mutations. Add confidence scoring to classification output. Strengthen VCF filtering, add anchor-based flanking trim, and add VCF-backed mutation validation to eliminate false positives. Disambiguate same-length alleles using Clair3 heterozygous genotype calls. Add bidirectional classification fallback for novel large mutations.

**Tech Stack:** Python 3.10+, Click CLI, bcftools, samtools, Clair3, pytest

---

## File Structure

| File | Responsibility | Tasks |
|------|---------------|-------|
| `src/muc_one_span/data/repeats/repeats.json` | Repeat + mutation definitions | 1 |
| `src/muc_one_span/config.py` | Load repeat dictionary + pre-compute mutated sequences | 1 |
| `src/muc_one_span/classify.py` | Classification with confidence, templates, bidirectional | 2, 3, 5, 8 |
| `src/muc_one_span/calling.py` | VCF quality filters + same-length disambiguation | 4, 7 |
| `src/muc_one_span/consensus.py` | Anchor-based flanking trim | 6 |
| `src/muc_one_span/alleles.py` | Same-length allele detection | 7 |
| `src/muc_one_span/cli.py` | Wire new features into CLI | 7 |
| `tests/unit/test_config.py` | Test mutation loading | 1 |
| `tests/unit/test_classify.py` | Test confidence, templates, bidirectional | 2, 3, 5, 8 |
| `tests/unit/test_calling.py` | Test VCF filters | 4 |
| `tests/unit/test_consensus.py` | Test anchor trim | 6 |
| `tests/unit/test_alleles.py` | Test same-length semantics | 7 |

---

### Task 1: Add mutation catalog to repeat dictionary

**Files:**
- Modify: `src/muc_one_span/data/repeats/repeats.json`
- Modify: `src/muc_one_span/config.py`
- Modify: `tests/unit/test_config.py`

- [ ] **Step 1: Write test for mutation loading**

Add to `tests/unit/test_config.py`:

```python
class TestMutationCatalog:
    """Tests for mutation catalog loading and sequence pre-computation."""

    def test_mutations_loaded(self):
        """Repeat dictionary loads mutation definitions."""
        rd = load_repeat_dictionary()
        assert hasattr(rd, "mutations")
        assert "dupC" in rd.mutations

    def test_mutation_has_required_fields(self):
        """Each mutation has allowed_repeats and changes."""
        rd = load_repeat_dictionary()
        dupc = rd.mutations["dupC"]
        assert "allowed_repeats" in dupc
        assert "changes" in dupc
        assert "X" in dupc["allowed_repeats"]

    def test_mutated_sequences_precomputed(self):
        """Pre-computed mutated sequences map sequence -> (repeat_id, mutation_name)."""
        rd = load_repeat_dictionary()
        assert hasattr(rd, "mutated_sequences")
        assert len(rd.mutated_sequences) > 0
        # dupC on X should produce a 61bp sequence
        found_dupc = any(
            name == "dupC" and repeat == "X"
            for seq, (repeat, name) in rd.mutated_sequences.items()
        )
        assert found_dupc

    def test_dupc_sequence_is_61bp(self):
        """dupC on X inserts C at position 59, producing 61bp."""
        rd = load_repeat_dictionary()
        for seq, (repeat_id, mut_name) in rd.mutated_sequences.items():
            if repeat_id == "X" and mut_name == "dupC":
                assert len(seq) == 61
                # The inserted C is at position 59 (0-indexed)
                x_seq = rd.repeats["X"]
                assert seq == x_seq[:60] + "C"
                break
        else:
            pytest.fail("dupC on X not found in mutated_sequences")

    def test_only_real_mutations_loaded(self):
        """Synthetic mutations (type=synthetic) are excluded."""
        rd = load_repeat_dictionary()
        assert "bigDel" not in rd.mutations
        assert "snpA_testing" not in rd.mutations
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_config.py::TestMutationCatalog -v --no-cov`
Expected: FAIL -- `mutations` attribute does not exist on RepeatDictionary

- [ ] **Step 3: Add mutations section to repeats.json**

Add the 13 real mutations from MucOneUp `config.json` to `src/muc_one_span/data/repeats/repeats.json`. Add a `"mutations"` key at the top level, after `"canonical_repeat"`. Only include mutations with `"type": "real"` (exclude synthetic test mutations).

```json
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X"],
      "changes": [{"type": "insert", "start": 60, "sequence": "C"}],
      "citation": "Kirby et al. 2013, PMID:23396133"
    },
    "dupA": {
      "allowed_repeats": ["X"],
      "changes": [{"type": "insert", "start": 60, "sequence": "A"}],
      "citation": "Olinger et al. 2020, PMID:32647000"
    },
    "insG": {
      "allowed_repeats": ["X"],
      "changes": [{"type": "insert", "start": 59, "sequence": "G"}],
      "citation": "Olinger et al. 2020, PMID:32647000"
    },
    "insCCCC": {
      "allowed_repeats": ["X"],
      "changes": [{"type": "insert", "start": 60, "sequence": "CCCC"}],
      "citation": "Vrbacka et al. 2025"
    },
    "insC_pos23": {
      "allowed_repeats": ["A", "E"],
      "changes": [{"type": "insert", "start": 23, "sequence": "C"}],
      "citation": "Vrbacka et al. 2025"
    },
    "insG_pos58": {
      "allowed_repeats": ["B", "X"],
      "changes": [{"type": "insert", "start": 58, "sequence": "G"}],
      "citation": "Vrbacka et al. 2025"
    },
    "insG_pos54": {
      "allowed_repeats": ["B", "J"],
      "changes": [{"type": "insert", "start": 54, "sequence": "G"}],
      "citation": "Vrbacka et al. 2025"
    },
    "insA_pos54": {
      "allowed_repeats": ["A", "H"],
      "changes": [{"type": "insert", "start": 54, "sequence": "A"}],
      "citation": "Vrbacka et al. 2025"
    },
    "del18_31": {
      "allowed_repeats": ["X"],
      "changes": [{"type": "delete", "start": 18, "end": 31}],
      "citation": "Vrbacka et al. 2025"
    },
    "ins16bp": {
      "allowed_repeats": ["C"],
      "changes": [{"type": "insert", "start": 57, "sequence": "GGGCTCCACCGCCCCC"}],
      "citation": "Vrbacka et al. 2025"
    },
    "ins25bp": {
      "allowed_repeats": ["A", "B", "J", "K", "N", "S", "X"],
      "changes": [{"type": "insert", "start": 31, "sequence": "CAGGCCGGCCCCGGGCTCCGGACAC"}],
      "citation": "Saei et al. 2023, PMID:37456840"
    },
    "delinsAT": {
      "allowed_repeats": ["X"],
      "changes": [{"type": "delete_insert", "start": 54, "end": 56, "sequence": "AT"}],
      "citation": "Olinger et al. 2020, PMID:32647000"
    },
    "delGCCCA": {
      "allowed_repeats": ["5", "5C", "6", "6p", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "V", "W", "X"],
      "changes": [{"type": "delete", "start": 1, "end": 5}],
      "citation": "Saei et al. 2023, PMID:37456840"
    }
  }
```

- [ ] **Step 4: Update RepeatDictionary dataclass and loader in config.py**

Add `mutations` and `mutated_sequences` fields to `RepeatDictionary`. Add `_apply_mutation()` helper and `_precompute_mutated_sequences()`. Update `load_repeat_dictionary()`.

```python
@dataclass
class RepeatDictionary:
    repeats: dict[str, str]
    repeat_length_bp: int
    pre_repeat_ids: list[str]
    after_repeat_ids: list[str]
    canonical_repeat: str
    flanking_left: str
    flanking_right: str
    vntr_region: str
    source: str
    mutations: dict[str, dict]
    mutated_sequences: dict[str, tuple[str, str]]


def _apply_mutation(sequence: str, changes: list[dict]) -> str:
    """Apply mutation changes to a repeat sequence.

    Uses 1-based indexing (matching MucOneUp/Vrbacka conventions).
    """
    result = list(sequence)
    # Apply changes in reverse order to preserve indices
    for change in sorted(changes, key=lambda c: c["start"], reverse=True):
        start = change["start"] - 1  # convert to 0-based
        if change["type"] == "insert":
            result.insert(start + 1, change["sequence"])
        elif change["type"] == "delete":
            end = change["end"] - 1  # 0-based inclusive
            del result[start:end]
        elif change["type"] == "delete_insert":
            end = change["end"] - 1
            result[start:end] = list(change["sequence"])
    return "".join(result)


def _precompute_mutated_sequences(
    repeats: dict[str, str],
    mutations: dict[str, dict],
) -> dict[str, tuple[str, str]]:
    """Pre-compute mutated sequences for all (repeat, mutation) combinations.

    Returns:
        Dict mapping mutated DNA sequence -> (parent_repeat_id, mutation_name).
    """
    result: dict[str, tuple[str, str]] = {}
    for mut_name, mut_def in mutations.items():
        for repeat_id in mut_def["allowed_repeats"]:
            if repeat_id not in repeats:
                continue
            seq = repeats[repeat_id]
            mutated = _apply_mutation(seq, mut_def["changes"])
            result[mutated] = (repeat_id, mut_name)
    return result
```

Update `load_repeat_dictionary()` to load mutations and pre-compute:

```python
def load_repeat_dictionary(path: Path | None = None) -> RepeatDictionary:
    # ... existing loading code ...
    mutations_raw = data.get("mutations", {})
    # Filter to real mutations only
    mutations = {k: v for k, v in mutations_raw.items() if v.get("type", "real") != "synthetic"}
    # Note: our JSON doesn't have a "type" field -- all entries are real since
    # we only copied real mutations from MucOneUp. Keep filter for safety.

    mutated_seqs = _precompute_mutated_sequences(data["repeats"], mutations)

    flanking = data["flanking_hg38"]
    return RepeatDictionary(
        repeats=data["repeats"],
        repeat_length_bp=data["repeat_length_bp"],
        pre_repeat_ids=data["pre_repeat_ids"],
        after_repeat_ids=data["after_repeat_ids"],
        canonical_repeat=data["canonical_repeat"],
        flanking_left=flanking["left"],
        flanking_right=flanking["right"],
        vntr_region=flanking["vntr_region"],
        source=data["source"],
        mutations=mutations,
        mutated_sequences=mutated_seqs,
    )
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_config.py -v --no-cov`
Expected: All PASS

- [ ] **Step 6: Run full test suite to check no regressions**

Run: `uv run pytest tests/unit/ -v --no-cov`
Expected: All PASS (existing tests may need minor updates if RepeatDictionary constructor changed)

- [ ] **Step 7: Commit**

```bash
git add src/muc_one_span/data/repeats/repeats.json src/muc_one_span/config.py tests/unit/test_config.py
git commit -m "feat: add known MUC1 mutation catalog with pre-computed templates"
```

---

### Task 2: Add confidence scoring to classify.py

**Files:**
- Modify: `src/muc_one_span/classify.py`
- Modify: `tests/unit/test_classify.py`

- [ ] **Step 1: Write tests for confidence scoring**

Add to `tests/unit/test_classify.py`:

```python
class TestConfidenceScoring:
    """Tests for per-repeat confidence scores."""

    def test_exact_match_confidence_is_1(self, repeat_dict):
        """Exact match to known repeat type has confidence 1.0."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_repeat(x_seq, repeat_dict)
        assert result["confidence"] == 1.0

    def test_close_match_confidence_uses_identity(self, repeat_dict):
        """Near-match (ed=1) confidence equals identity_pct / 100."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]  # 61bp, ed=1
        result = classify_repeat(mutated, repeat_dict)
        assert 0.9 < result["confidence"] < 1.0
        assert result["confidence"] == result["identity_pct"] / 100

    def test_distant_match_has_lower_confidence(self, repeat_dict):
        """High edit distance produces lower confidence."""
        # Random 60bp sequence
        seq = "ACGT" * 15
        result = classify_repeat(seq, repeat_dict)
        assert result["confidence"] < 0.9


class TestSequenceConfidenceSummary:
    """Tests for per-allele confidence summary in classify_sequence."""

    def test_all_exact_gives_confidence_1(self, repeat_dict):
        """All exact matches produce allele_confidence=1.0."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_sequence(x_seq * 3, repeat_dict)
        assert result["allele_confidence"] == 1.0
        assert result["exact_match_pct"] == 100.0

    def test_mixed_match_reduces_confidence(self, repeat_dict):
        """One mutation among exact matches reduces allele_confidence below 1.0."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]
        result = classify_sequence(x_seq + mutated + x_seq, repeat_dict)
        assert result["allele_confidence"] < 1.0
        assert result["exact_match_pct"] < 100.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_classify.py::TestConfidenceScoring -v --no-cov`
Expected: FAIL -- `confidence` key not in result

- [ ] **Step 3: Add confidence to classify_repeat()**

In `src/muc_one_span/classify.py`, modify `classify_repeat()`:

For exact matches, add `"confidence": 1.0` to the return dict.

For non-exact matches, compute confidence as `identity_pct / 100`:

```python
    # At the end of classify_repeat(), before return:
    result["confidence"] = identity_pct / 100
```

For exact matches:
```python
    if sequence in seq_to_id:
        return {"type": seq_to_id[sequence], "match": "exact", "confidence": 1.0}
```

- [ ] **Step 4: Add summary metrics to classify_sequence()**

At the end of `classify_sequence()`, before the return statement, compute:

```python
    confidences = [r.get("confidence", 1.0) for r in repeats]
    exact_count = sum(1 for r in repeats if r.get("match") == "exact")
    allele_confidence = sum(confidences) / len(confidences) if confidences else 0.0
    exact_match_pct = (exact_count / len(repeats) * 100) if repeats else 0.0

    return {
        "structure": " ".join(labels),
        "repeats": repeats,
        "mutations_detected": mutations,
        "cumulative_offset": cumulative_offset,
        "allele_confidence": round(allele_confidence, 4),
        "exact_match_pct": round(exact_match_pct, 1),
    }
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_classify.py -v --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/classify.py tests/unit/test_classify.py
git commit -m "feat: add per-repeat confidence scoring and allele summary metrics"
```

---

### Task 3: Exact-match-first probing with mutation templates

**Files:**
- Modify: `src/muc_one_span/classify.py`
- Modify: `tests/unit/test_classify.py`

- [ ] **Step 1: Write tests for mutation template matching**

Add to `tests/unit/test_classify.py`:

```python
class TestMutationTemplateMatching:
    """Tests for exact matching against known mutation templates."""

    def test_dupc_exact_template_match(self, repeat_dict):
        """dupC on X matches the pre-computed 61bp template exactly."""
        x_seq = repeat_dict.repeats["X"]
        dupc_seq = x_seq[:60] + "C"  # 61bp: insert C at end
        result = classify_repeat(dupc_seq, repeat_dict)
        assert result["match"] == "exact"
        assert result["type"] == "X:dupC"
        assert result["confidence"] == 1.0

    def test_ins16bp_exact_template_match(self, repeat_dict):
        """16bp insertion on C matches the 76bp template exactly."""
        c_seq = repeat_dict.repeats["C"]
        ins_seq = c_seq[:57] + "GGGCTCCACCGCCCCC" + c_seq[57:]  # 76bp
        result = classify_repeat(ins_seq, repeat_dict)
        assert result["match"] == "exact"
        assert result["type"] == "C:ins16bp"
        assert result["confidence"] == 1.0

    def test_del18_31_exact_template_match(self, repeat_dict):
        """14bp deletion on X matches the 46bp template exactly."""
        x_seq = repeat_dict.repeats["X"]
        del_seq = x_seq[:17] + x_seq[30:]  # 46bp: delete pos 18-31 (1-based)
        result = classify_repeat(del_seq, repeat_dict)
        assert result["match"] == "exact"
        assert result["type"] == "X:del18_31"

    def test_sequence_with_dupc_classifies_correctly(self, repeat_dict):
        """classify_sequence finds dupC at correct window size."""
        x_seq = repeat_dict.repeats["X"]
        dupc_seq = x_seq[:60] + "C"  # 61bp
        full = x_seq + dupc_seq + x_seq
        result = classify_sequence(full, repeat_dict)
        assert len(result["repeats"]) == 3
        assert result["repeats"][1]["type"] == "X:dupC"
        assert result["repeats"][1]["match"] == "exact"

    def test_sequence_with_16bp_ins_classifies_correctly(self, repeat_dict):
        """classify_sequence handles 76bp mutated repeat with template match."""
        x_seq = repeat_dict.repeats["X"]
        c_seq = repeat_dict.repeats["C"]
        ins_seq = c_seq[:57] + "GGGCTCCACCGCCCCC" + c_seq[57:]
        full = x_seq + ins_seq + x_seq
        result = classify_sequence(full, repeat_dict)
        assert len(result["repeats"]) == 3
        assert result["repeats"][1]["type"] == "C:ins16bp"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_classify.py::TestMutationTemplateMatching -v --no-cov`
Expected: FAIL -- mutation templates not checked

- [ ] **Step 3: Update classify_repeat() to check mutation templates**

Modify `classify_repeat()` in `src/muc_one_span/classify.py`:

```python
def classify_repeat(
    sequence: str,
    repeat_dict: RepeatDictionary,
) -> dict:
    # O(1) exact-match lookup via reverse map (sequence -> ID)
    seq_to_id = {seq: rid for rid, seq in repeat_dict.repeats.items()}
    if sequence in seq_to_id:
        return {"type": seq_to_id[sequence], "match": "exact", "confidence": 1.0}

    # Check mutation templates (variable-length exact matches)
    if hasattr(repeat_dict, "mutated_sequences") and sequence in repeat_dict.mutated_sequences:
        parent_repeat, mut_name = repeat_dict.mutated_sequences[sequence]
        return {
            "type": f"{parent_repeat}:{mut_name}",
            "match": "exact",
            "confidence": 1.0,
            "mutation_name": mut_name,
            "parent_repeat": parent_repeat,
        }

    # ... rest of existing edit distance logic (unchanged) ...
```

- [ ] **Step 4: Update classify_sequence() with exact-match-first probing**

In the probing loop of `classify_sequence()`, add a first pass that checks ALL probe sizes for exact matches (against both standard repeats and mutation templates) before falling back to edit distance:

```python
        # --- Phase 1: Check ALL probe sizes for exact match ---
        exact_found = False
        for probe_size in _probe_sizes_generator(unit_length, max_indel_probe, remaining):
            window = sequence[pos : pos + probe_size]
            # Check standard repeats
            if window in seq_to_id:
                best_result = {"type": seq_to_id[window], "match": "exact", "confidence": 1.0}
                best_window_size = probe_size
                best_dist = 0
                exact_found = True
                break
            # Check mutation templates
            if hasattr(repeat_dict, "mutated_sequences") and window in repeat_dict.mutated_sequences:
                parent, mname = repeat_dict.mutated_sequences[window]
                best_result = {
                    "type": f"{parent}:{mname}",
                    "match": "exact",
                    "confidence": 1.0,
                    "mutation_name": mname,
                    "parent_repeat": parent,
                }
                best_window_size = probe_size
                best_dist = 0
                exact_found = True
                break

        # --- Phase 2: Edit distance fallback (only if no exact match) ---
        if not exact_found:
            # ... existing probing logic unchanged ...
```

Add the helper:

```python
def _probe_sizes_generator(
    unit_length: int, max_indel_probe: int, remaining: int
) -> list[int]:
    """Generate probe sizes: canonical first, then small-to-large."""
    sizes = [min(unit_length, remaining)]
    for ps in range(
        max(unit_length - max_indel_probe, unit_length // 2),
        min(unit_length + max_indel_probe + 1, remaining + 1),
    ):
        if ps != unit_length:
            sizes.append(ps)
    return sizes
```

Build `seq_to_id` once at the top of `classify_sequence()`:

```python
    seq_to_id = {seq: rid for rid, seq in repeat_dict.repeats.items()}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_classify.py -v --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/classify.py tests/unit/test_classify.py
git commit -m "feat: exact-match-first probing with mutation templates (fixes #2 for known mutations)"
```

---

### Task 4: Strengthen VCF quality filters

**Files:**
- Modify: `src/muc_one_span/calling.py`
- Modify: `tests/unit/test_calling.py`

- [ ] **Step 1: Write test for VCF quality filtering**

Add to `tests/unit/test_calling.py`:

```python
class TestFilterVcfQuality:
    """Tests for VCF quality filter parameters."""

    @patch("muc_one_span.calling.run_tool", return_value="")
    def test_filter_vcf_includes_quality_expression(self, mock_run_tool, tmp_path):
        """filter_vcf passes QUAL and DP filter to bcftools view."""
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        filter_vcf(vcf, ref, tmp_path, min_qual=15.0, min_dp=5)

        # Find the bcftools view call
        view_calls = [
            c[0][0] for c in mock_run_tool.call_args_list
            if c[0][0][:2] == ["bcftools", "view"]
        ]
        assert len(view_calls) == 1
        view_cmd = view_calls[0]
        assert "-i" in view_cmd
        i_idx = view_cmd.index("-i")
        expr = view_cmd[i_idx + 1]
        assert "QUAL" in expr
        assert "DP" in expr

    @patch("muc_one_span.calling.run_tool", return_value="")
    def test_filter_vcf_default_params(self, mock_run_tool, tmp_path):
        """filter_vcf works with default parameters (backward compatible)."""
        vcf = tmp_path / "input.vcf.gz"
        vcf.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        # Should not raise with no extra args
        filter_vcf(vcf, ref, tmp_path)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_calling.py::TestFilterVcfQuality -v --no-cov`
Expected: FAIL -- `filter_vcf()` does not accept `min_qual`/`min_dp`

- [ ] **Step 3: Add quality parameters to filter_vcf()**

Modify `filter_vcf()` in `src/muc_one_span/calling.py`:

```python
def filter_vcf(
    vcf_path: Path,
    reference_path: Path,
    output_dir: Path,
    min_qual: float = 0.0,
    min_dp: int = 0,
) -> Path:
```

Change the bcftools view command to include quality filtering:

```python
    # Build filter expression
    view_cmd = [
        "bcftools",
        "view",
        "-f",
        "PASS",
    ]
    # Add quality filters if specified
    filters = []
    if min_qual > 0:
        filters.append(f"QUAL>={min_qual}")
    if min_dp > 0:
        filters.append(f"INFO/DP>={min_dp}")
    if filters:
        view_cmd.extend(["-i", " && ".join(filters)])

    view_cmd.extend([
        "-o",
        str(filtered),
        "-O",
        "z",
        str(norm_vcf),
    ])
    run_tool(view_cmd)
```

- [ ] **Step 4: Update call_variants_per_allele to pass quality params**

Add `min_qual` and `min_dp` parameters to `call_variants_per_allele()` and pass them through to `filter_vcf()`. Use defaults of `min_qual=15.0` and `min_dp=5`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_calling.py -v --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/calling.py tests/unit/test_calling.py
git commit -m "feat: add VCF quality filters (QUAL, DP) to reduce false positives (fixes #3)"
```

---

### Task 5: VCF-backed mutation validation

**Files:**
- Modify: `src/muc_one_span/classify.py`
- Modify: `tests/unit/test_classify.py`

- [ ] **Step 1: Write tests for VCF validation**

Add to `tests/unit/test_classify.py`:

```python
from muc_one_span.classify import validate_mutations_against_vcf


class TestVcfMutationValidation:
    """Tests for VCF-backed mutation validation."""

    def test_confirmed_mutation_keeps_high_confidence(self, repeat_dict):
        """Mutation with VCF support keeps confidence unchanged."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]
        result = classify_sequence(x_seq + mutated + x_seq, repeat_dict)

        # Simulate VCF with a variant at repeat 2 position
        vcf_variants = [{"pos": 560, "qual": 25.0}]  # flank(500) + 60bp

        validated = validate_mutations_against_vcf(
            result, vcf_variants=vcf_variants, flank_length=500, unit_length=60
        )
        mut = validated["mutations_detected"][0]
        assert mut.get("vcf_support") is True

    def test_unsupported_mutation_gets_low_confidence(self, repeat_dict):
        """Mutation without VCF support gets reduced confidence."""
        x_seq = repeat_dict.repeats["X"]
        mutated = x_seq[:59] + "C" + x_seq[59:]
        result = classify_sequence(x_seq + mutated + x_seq, repeat_dict)

        # No VCF variants at all
        validated = validate_mutations_against_vcf(
            result, vcf_variants=[], flank_length=500, unit_length=60
        )
        mut = validated["mutations_detected"][0]
        assert mut.get("vcf_support") is False

    def test_no_mutations_unchanged(self, repeat_dict):
        """Sequence with no mutations passes through unchanged."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_sequence(x_seq * 3, repeat_dict)

        validated = validate_mutations_against_vcf(
            result, vcf_variants=[], flank_length=500, unit_length=60
        )
        assert validated["mutations_detected"] == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_classify.py::TestVcfMutationValidation -v --no-cov`
Expected: FAIL -- `validate_mutations_against_vcf` not found

- [ ] **Step 3: Implement validate_mutations_against_vcf()**

Add to `src/muc_one_span/classify.py`:

```python
def validate_mutations_against_vcf(
    classification_result: dict,
    vcf_variants: list[dict] | None = None,
    flank_length: int = 500,
    unit_length: int = 60,
) -> dict:
    """Cross-validate detected mutations against VCF variant positions.

    For each mutation, check if a VCF variant overlaps the repeat's
    genomic position. Adds ``vcf_support`` flag and adjusts confidence.

    Args:
        classification_result: Output from :func:`classify_sequence`.
        vcf_variants: List of dicts with ``pos`` (int) and ``qual`` (float).
            If None, VCF validation is skipped (standalone classify mode).
        flank_length: Flanking bp on each side of the contig.
        unit_length: Expected repeat unit length (60bp).

    Returns:
        Updated classification result with VCF validation annotations.
    """
    result = classification_result.copy()
    result["mutations_detected"] = [m.copy() for m in result.get("mutations_detected", [])]

    if vcf_variants is None:
        return result

    for mutation in result["mutations_detected"]:
        repeat_idx = mutation["repeat_index"]
        # Map repeat index to contig coordinates
        repeat_start = flank_length + (repeat_idx - 1) * unit_length
        repeat_end = repeat_start + unit_length + 30  # allow for indels

        # Check if any VCF variant overlaps this repeat
        supporting = [
            v for v in vcf_variants
            if repeat_start <= v["pos"] <= repeat_end
        ]
        mutation["vcf_support"] = len(supporting) > 0
        mutation["vcf_qual"] = max((v["qual"] for v in supporting), default=0.0)

        # Adjust confidence in the corresponding repeat
        repeat_result = result["repeats"][repeat_idx - 1]
        base_confidence = repeat_result.get("confidence", 1.0)
        if supporting:
            best_qual = mutation["vcf_qual"]
            vcf_score = 1.0 if best_qual >= 20 else 0.7
        else:
            vcf_score = 0.3
        repeat_result["confidence"] = round(base_confidence * vcf_score, 4)

    # Recompute allele_confidence
    confidences = [r.get("confidence", 1.0) for r in result["repeats"]]
    result["allele_confidence"] = round(
        sum(confidences) / len(confidences), 4
    ) if confidences else 0.0

    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_classify.py -v --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/classify.py tests/unit/test_classify.py
git commit -m "feat: add VCF-backed mutation validation with confidence adjustment"
```

---

### Task 6: Anchor-based flanking trim

**Files:**
- Modify: `src/muc_one_span/consensus.py`
- Modify: `tests/unit/test_consensus.py`

- [ ] **Step 1: Write tests for anchor-based trim**

Add to `tests/unit/test_consensus.py`:

```python
from muc_one_span.config import load_repeat_dictionary


class TestAnchorBasedTrim:
    """Tests for anchor-based flanking trim."""

    def test_anchor_trim_handles_indel_in_flank(self, tmp_path):
        """Anchor trim finds correct boundary despite indel in flanking."""
        rd = load_repeat_dictionary()
        # Build a sequence: left_flank (with 1bp insertion) + pre-repeat 1 + X + after-repeat 9 + right_flank
        left_flank = rd.flanking_left[-500:]
        right_flank = rd.flanking_right[:500]
        vntr = rd.repeats["1"] + rd.repeats["X"] + rd.repeats["9"]
        # Insert a fake base in the left flank (simulating Clair3 false positive)
        corrupted_left = left_flank[:250] + "A" + left_flank[250:]
        sequence = corrupted_left + vntr + right_flank  # 501 + vntr + 500

        fasta = tmp_path / "full.fa"
        fasta.write_text(f">contig\n{sequence}\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 500, output, repeat_dict=rd)

        trimmed = output.read_text().strip().splitlines()
        trimmed_seq = "".join(ln for ln in trimmed if not ln.startswith(">"))
        # Should contain the VNTR, not be offset by the flank indel
        assert rd.repeats["1"][:30] in trimmed_seq

    def test_anchor_trim_fallback_without_repeat_dict(self, tmp_path):
        """Without repeat_dict, falls back to fixed-position trim."""
        sequence = "A" * 100 + "VNTR_SEQ" + "A" * 100
        fasta = tmp_path / "full.fa"
        fasta.write_text(f">contig\n{sequence}\n")
        output = tmp_path / "trimmed.fa"

        trim_flanking(fasta, 100, output, repeat_dict=None)

        trimmed = output.read_text().strip().splitlines()
        trimmed_seq = "".join(ln for ln in trimmed if not ln.startswith(">"))
        assert trimmed_seq == "VNTR_SEQ"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_consensus.py::TestAnchorBasedTrim -v --no-cov`
Expected: FAIL -- `trim_flanking()` does not accept `repeat_dict`

- [ ] **Step 3: Add _find_anchor() helper and update trim_flanking()**

Add to `src/muc_one_span/consensus.py`:

```python
from muc_one_span.config import RepeatDictionary


def _find_anchor(
    sequence: str,
    anchor: str,
    expected_pos: int,
    tolerance: int = 50,
) -> int | None:
    """Find anchor sequence near expected position.

    Searches for an exact substring match within ``expected_pos +/- tolerance``.
    Returns the position where the anchor ENDS (i.e., the start of the VNTR).

    Returns None if not found.
    """
    search_start = max(0, expected_pos - tolerance)
    search_end = min(len(sequence), expected_pos + tolerance + len(anchor))
    region = sequence[search_start:search_end]
    idx = region.find(anchor)
    if idx >= 0:
        return search_start + idx + len(anchor)
    return None
```

Update `trim_flanking()` signature and body:

```python
def trim_flanking(
    consensus_fasta: Path,
    flank_length: int,
    output_path: Path,
    repeat_dict: RepeatDictionary | None = None,
) -> Path:
```

Add anchor logic before the existing fixed-position trim:

```python
    # Try anchor-based trimming if repeat_dict provided
    left_trim = flank_length
    right_trim = len(sequence) - flank_length if flank_length > 0 else len(sequence)

    if repeat_dict is not None and flank_length > 0:
        # Left anchor: last 20bp of left flank + first 20bp of pre-repeat "1"
        left_anchor = repeat_dict.flanking_left[-20:] + repeat_dict.repeats["1"][:20]
        anchor_pos = _find_anchor(sequence, left_anchor, flank_length)
        if anchor_pos is not None:
            # anchor_pos points to end of anchor (= start of repeat "1" + 20bp)
            # We want to start at the beginning of repeat "1"
            left_trim = anchor_pos - 20  # back up to start of repeat "1"

        # Right anchor: last 20bp of after-repeat "9" + first 20bp of right flank
        right_anchor = repeat_dict.repeats["9"][-20:] + repeat_dict.flanking_right[:20]
        right_anchor_pos = _find_anchor(
            sequence, right_anchor, len(sequence) - flank_length
        )
        if right_anchor_pos is not None:
            right_trim = right_anchor_pos - 20  # end after repeat "9"

    vntr = sequence[left_trim:right_trim]
```

- [ ] **Step 4: Update build_consensus_per_allele() to pass repeat_dict**

Add `repeat_dict` parameter to `build_consensus_per_allele()` and pass it through to `trim_flanking()`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_consensus.py -v --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add src/muc_one_span/consensus.py tests/unit/test_consensus.py
git commit -m "feat: anchor-based flanking trim resilient to consensus indels (fixes #3)"
```

---

### Task 7: Same-length allele disambiguation (Issue #1)

**Files:**
- Modify: `src/muc_one_span/alleles.py`
- Modify: `src/muc_one_span/calling.py`
- Modify: `src/muc_one_span/consensus.py`
- Modify: `src/muc_one_span/cli.py`
- Modify: `tests/unit/test_alleles.py`

- [ ] **Step 1: Write tests for same_length flag**

Add to `tests/unit/test_alleles.py`:

```python
class TestSameLengthDetection:
    """Tests for same_length allele detection."""

    def test_single_cluster_sets_same_length(self):
        """Single cluster sets same_length=True."""
        counts = parse_idxstats(IDXSTATS_HOMOZYGOUS)
        result = detect_alleles(counts, min_coverage=10)
        assert result["same_length"] is True

    def test_same_length_keeps_allele2_data(self):
        """same_length allele_2 has full data (not reads=None)."""
        counts = parse_idxstats(IDXSTATS_HOMOZYGOUS)
        result = detect_alleles(counts, min_coverage=10)
        assert result["allele_2"]["contig_name"] is not None

    def test_different_lengths_no_same_length(self):
        """Different-length alleles have same_length=False."""
        counts = parse_idxstats(IDXSTATS_TWO_PEAKS)
        result = detect_alleles(counts, min_coverage=10)
        assert result.get("same_length", False) is False
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_alleles.py::TestSameLengthDetection -v --no-cov`
Expected: FAIL -- `same_length` key not in result

- [ ] **Step 3: Update detect_alleles() in alleles.py**

Change the two homozygous code paths to use `same_length=True` and `homozygous=False`:

For single cluster (line ~284):
```python
    if len(clusters) < 2:
        return {
            "allele_1": allele_1,
            "allele_2": {**allele_1},
            "homozygous": False,
            "same_length": True,
        }
```

For two clusters with same length (line ~293):
```python
    if allele_1["length"] == allele_2["length"]:
        allele_1["reads"] += allele_2["reads"]
        return {
            "allele_1": allele_1,
            "allele_2": {**allele_1},
            "homozygous": False,
            "same_length": True,
        }
```

For different lengths, add `same_length=False`:
```python
    return {
        "allele_1": allele_1,
        "allele_2": allele_2,
        "homozygous": False,
        "same_length": False,
    }
```

- [ ] **Step 4: Add parse_vcf_genotypes() to calling.py**

```python
def parse_vcf_genotypes(vcf_path: Path) -> list[dict]:
    """Parse VCF to extract variant positions and genotypes.

    Returns list of dicts with keys: chrom, pos, ref, alt, genotype.
    """
    try:
        output = run_tool([
            "bcftools", "query", "-f",
            "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n",
            str(vcf_path),
        ])
    except RuntimeError:
        return []

    variants = []
    for line in output.strip().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        variants.append({
            "chrom": fields[0],
            "pos": int(fields[1]),
            "ref": fields[2],
            "alt": fields[3],
            "genotype": fields[4],
        })
    return variants
```

- [ ] **Step 5: Add disambiguate_same_length_alleles() to calling.py**

```python
def disambiguate_same_length_alleles(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = "",
    threads: int = 4,
    min_qual: float = 15.0,
    min_dp: int = 5,
) -> dict[str, Path]:
    """Disambiguate same-length alleles using Clair3 genotype calls.

    Runs Clair3 on all reads, checks for heterozygous (0/1) variants.
    If found, creates two VCFs: one with only hom-alt variants (WT allele),
    one with all variants forced to hom-alt (mutant allele).

    Returns dict mapping allele key to filtered VCF path.
    Updates alleles["homozygous"] based on findings.
    """
    allele_info = alleles["allele_1"]
    contig_name = allele_info["contig_name"]
    cluster_contigs = allele_info["cluster_contigs"]
    merged_dir = output_dir / "merged"

    # Remap all cluster reads to peak contig
    merged_bam = _extract_and_remap_reads(
        bam_path, cluster_contigs, contig_name,
        reference_path, merged_dir, threads,
    )

    # Run Clair3
    contig_ref = merged_dir / f"{contig_name}.fa"
    clair3_dir = merged_dir / "clair3"
    raw_vcf = run_clair3(merged_bam, contig_ref, clair3_dir,
                          model_path=clair3_model, threads=threads)
    filtered_vcf = filter_vcf(raw_vcf, contig_ref, merged_dir,
                               min_qual=min_qual, min_dp=min_dp)

    # Check genotypes
    variants = parse_vcf_genotypes(filtered_vcf)
    het_variants = [v for v in variants if "/" in v["genotype"]
                    and v["genotype"] not in ("0/0", "1/1", "./.", ".|.")]

    results: dict[str, Path] = {}

    if not het_variants:
        # Truly homozygous -- no het variants found
        alleles["homozygous"] = True
        results["allele_1"] = filtered_vcf
        return results

    # Compound heterozygous -- split into two VCFs
    alleles["homozygous"] = False

    # allele_1 (WT): exclude het variants, keep only hom-alt
    wt_vcf = output_dir / "allele_1" / "variants.vcf.gz"
    wt_vcf.parent.mkdir(parents=True, exist_ok=True)
    run_tool([
        "bcftools", "view",
        "-i", "GT=\"1/1\" || GT=\"1|1\"",
        "-o", str(wt_vcf), "-O", "z",
        str(filtered_vcf),
    ])
    run_tool(["bcftools", "index", str(wt_vcf)])
    results["allele_1"] = wt_vcf

    # allele_2 (mutant): include all variants, force het to hom-alt
    # Use bcftools +setGT to convert 0/1 -> 1/1, then output
    mut_vcf_raw = output_dir / "allele_2" / "variants_raw.vcf.gz"
    mut_vcf = output_dir / "allele_2" / "variants.vcf.gz"
    mut_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Copy all variants
    run_tool([
        "bcftools", "view",
        "-o", str(mut_vcf), "-O", "z",
        str(filtered_vcf),
    ])
    run_tool(["bcftools", "index", str(mut_vcf)])
    results["allele_2"] = mut_vcf

    return results
```

- [ ] **Step 6: Update call_variants_per_allele() to handle same_length**

Add at the beginning of `call_variants_per_allele()`:

```python
    if alleles.get("same_length"):
        return disambiguate_same_length_alleles(
            bam_path, reference_path, alleles, output_dir,
            clair3_model, threads,
        )
```

- [ ] **Step 7: Update calling.py skip logic for homozygous**

Change the existing skip logic to check both `homozygous` and `same_length`:

```python
        if alleles.get("homozygous") and not alleles.get("same_length") and allele_key == "allele_2":
            continue
```

- [ ] **Step 8: Run tests**

Run: `uv run pytest tests/unit/ -v --no-cov`
Expected: All PASS (some existing tests may need updating for new `same_length` field)

- [ ] **Step 9: Commit**

```bash
git add src/muc_one_span/alleles.py src/muc_one_span/calling.py tests/unit/test_alleles.py
git commit -m "feat: disambiguate same-length alleles using Clair3 het genotypes (fixes #1)"
```

---

### Task 8: Bidirectional classification fallback

**Files:**
- Modify: `src/muc_one_span/classify.py`
- Modify: `tests/unit/test_classify.py`

- [ ] **Step 1: Write tests for bidirectional classification**

Add to `tests/unit/test_classify.py`:

```python
class TestBidirectionalClassification:
    """Tests for bidirectional classification fallback."""

    def test_novel_large_insertion_classified_via_bidirectional(self, repeat_dict):
        """Novel large insertion (not in templates) uses bidirectional fallback."""
        x_seq = repeat_dict.repeats["X"]
        pre1 = repeat_dict.repeats["1"]
        after9 = repeat_dict.repeats["9"]
        # Insert 20bp of random sequence into X (novel, not in templates)
        novel_mutated = x_seq[:30] + "A" * 20 + x_seq[30:]  # 80bp
        full = pre1 + x_seq + novel_mutated + x_seq + after9
        result = classify_sequence(full, repeat_dict)
        # Should classify pre1, 2x X, after9 correctly, with the mutated one in between
        labels = result["structure"].split()
        assert labels[0] == "1"
        assert labels[-1] == "9"
        # The novel mutation should still be detected (even if not exact)
        assert len(result["mutations_detected"]) >= 1 or any(
            "m" in label for label in labels
        )

    def test_forward_only_when_no_large_mutation(self, repeat_dict):
        """No bidirectional needed when all repeats classify well."""
        x_seq = repeat_dict.repeats["X"]
        result = classify_sequence(x_seq * 5, repeat_dict)
        assert all(r["match"] == "exact" for r in result["repeats"])
        assert result["allele_confidence"] == 1.0
```

- [ ] **Step 2: Run tests to verify current behavior**

Run: `uv run pytest tests/unit/test_classify.py::TestBidirectionalClassification -v --no-cov`
Expected: May pass or fail depending on current fallback behavior

- [ ] **Step 3: Implement backward pass in classify_sequence()**

Add `_classify_backward()` helper to `classify.py`:

```python
def _classify_backward(
    sequence: str,
    repeat_dict: RepeatDictionary,
    stop_pos: int,
) -> list[tuple[dict, int, int]]:
    """Classify repeats from 3' end backward, anchored on after-repeats.

    Returns list of (result, start_pos, end_pos) tuples in forward order.
    Stops when reaching stop_pos or when confidence drops.
    """
    unit_length = repeat_dict.repeat_length_bp
    after_ids = list(reversed(repeat_dict.after_repeat_ids))
    seq_to_id = {seq: rid for rid, seq in repeat_dict.repeats.items()}

    results: list[tuple[dict, int, int]] = []
    pos = len(sequence)

    for expected_id in after_ids:
        if pos - unit_length < stop_pos:
            break
        window = sequence[pos - unit_length : pos]
        if window in seq_to_id and seq_to_id[window] == expected_id:
            result = {"type": expected_id, "match": "exact", "confidence": 1.0}
            results.append((result, pos - unit_length, pos))
            pos -= unit_length
        else:
            break

    # Continue backward through canonical region
    max_indel_probe = 30
    while pos - unit_length // 2 > stop_pos:
        remaining_back = pos - stop_pos
        if remaining_back < unit_length // 2:
            break

        # Try exact match at canonical size first
        best_result = None
        best_size = unit_length
        best_dist = float("inf")

        for probe_size in _probe_sizes_generator(unit_length, max_indel_probe, remaining_back):
            start = pos - probe_size
            if start < stop_pos:
                continue
            window = sequence[start:pos]
            if window in seq_to_id:
                best_result = {"type": seq_to_id[window], "match": "exact", "confidence": 1.0}
                best_size = probe_size
                best_dist = 0
                break
            if hasattr(repeat_dict, "mutated_sequences") and window in repeat_dict.mutated_sequences:
                parent, mname = repeat_dict.mutated_sequences[window]
                best_result = {"type": f"{parent}:{mname}", "match": "exact", "confidence": 1.0}
                best_size = probe_size
                best_dist = 0
                break

        if best_dist > 0:
            # Edit distance fallback
            window = sequence[max(stop_pos, pos - unit_length) : pos]
            best_result = classify_repeat(window, repeat_dict)
            best_size = len(window)
            best_dist = best_result.get("edit_distance", 999)

        if best_result is None or best_dist > 3:
            break

        results.append((best_result, pos - best_size, pos))
        pos -= best_size

    results.reverse()
    return results
```

Modify `classify_sequence()` to use bidirectional fallback when a repeat has high edit distance. After the existing forward pass loop, if the sequence is not fully consumed and a high-ed repeat was found:

```python
    # After the forward pass loop ends:
    if pos < len(sequence) - unit_length // 2:
        # Forward pass stopped -- try backward pass
        backward = _classify_backward(sequence, repeat_dict, pos)
        if backward:
            # Gap between forward and backward = mutated region
            gap_start = pos
            gap_end = backward[0][1]
            if gap_end > gap_start:
                gap_seq = sequence[gap_start:gap_end]
                gap_result = classify_repeat(gap_seq, repeat_dict)
                gap_result["index"] = repeat_index
                repeat_index += 1
                # Label as mutation
                if gap_result.get("classification") == "mutation" or gap_result["match"] != "exact":
                    labels.append(f"{gap_result.get('closest_match', '?')}m")
                    mutations.append({
                        "repeat_index": gap_result["index"],
                        "closest_type": gap_result.get("closest_match", "?"),
                        "differences": gap_result.get("differences", []),
                        "frameshift": gap_result.get("frameshift", False),
                    })
                else:
                    labels.append(gap_result["type"])
                repeats.append(gap_result)

            # Append backward results
            for bwd_result, bwd_start, bwd_end in backward:
                repeat_index += 1
                bwd_result["index"] = repeat_index
                if bwd_result["match"] == "exact":
                    labels.append(bwd_result["type"])
                else:
                    labels.append(f"?{bwd_result.get('closest_match', '?')}")
                repeats.append(bwd_result)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/unit/test_classify.py -v --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/classify.py tests/unit/test_classify.py
git commit -m "feat: bidirectional classification fallback for novel large mutations (fixes #2)"
```

---

### Task 9: Wire features into CLI

**Files:**
- Modify: `src/muc_one_span/cli.py`

- [ ] **Step 1: Update `run` command to pass repeat_dict to consensus**

In `cli.py`, update the consensus step in the `run` command:

```python
    # Step 4: Build consensus
    click.echo("Step 4/5: Building consensus...")
    consensus_paths = build_consensus_per_allele(
        ref, vcf_paths, alleles_result, out, repeat_dict=rd
    )
```

- [ ] **Step 2: Update `run` command to pass VCF info to validation**

After classification, add VCF validation when VCF paths are available:

```python
    # Step 5: Classify repeats
    click.echo("Step 5/5: Classifying repeats...")
    all_results: dict[str, dict] = {}
    for allele_key, fa_path in consensus_paths.items():
        fa_lines = fa_path.read_text().strip().splitlines()
        sequence = "".join(line for line in fa_lines if not line.startswith(">"))
        result = classify_sequence(sequence, rd)

        # VCF-backed validation if VCF available
        if allele_key in vcf_paths:
            from muc_one_span.classify import validate_mutations_against_vcf
            vcf_variants = _parse_vcf_for_validation(vcf_paths[allele_key])
            result = validate_mutations_against_vcf(
                result, vcf_variants=vcf_variants,
            )

        all_results[allele_key] = result
        click.echo(f"  {allele_key}: {result['structure']}")
        if result.get("allele_confidence") is not None:
            click.echo(f"    confidence: {result['allele_confidence']:.2f}")
```

Add helper:

```python
def _parse_vcf_for_validation(vcf_path: Path) -> list[dict]:
    """Parse VCF variants for mutation validation."""
    from muc_one_span.tools import run_tool
    try:
        output = run_tool([
            "bcftools", "query", "-f",
            "%POS\\t%QUAL\\n",
            str(vcf_path),
        ])
    except (RuntimeError, FileNotFoundError):
        return []
    variants = []
    for line in output.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            variants.append({"pos": int(parts[0]), "qual": float(parts[1])})
    return variants
```

- [ ] **Step 3: Update `consensus` subcommand to pass repeat_dict**

```python
    fastas = build_consensus_per_allele(
        Path(reference),
        vcf_paths,
        alleles_data,
        out,
        repeat_dict=rd if repeats_db else None,
    )
```

Add `--repeats-db` option to the consensus subcommand if not already present.

- [ ] **Step 4: Run full test suite**

Run: `uv run pytest tests/unit/ -v --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/cli.py
git commit -m "feat: wire mutation catalog, confidence scoring, and VCF validation into CLI"
```

---

### Task 10: Quality checks and final cleanup

- [ ] **Step 1: Run linter**

Run: `uv run ruff check src/ tests/`
Fix any issues.

- [ ] **Step 2: Run formatter**

Run: `uv run ruff format src/ tests/`

- [ ] **Step 3: Run type checker**

Run: `uv run mypy src/muc_one_span/`
Fix any type errors.

- [ ] **Step 4: Run full test suite with coverage**

Run: `uv run pytest tests/unit/ --cov=muc_one_span --cov-report=term-missing`
Expected: All PASS, coverage >= 70%

- [ ] **Step 5: Commit any fixes**

```bash
git add -u
git commit -m "chore: fix lint, format, and type check issues"
```
