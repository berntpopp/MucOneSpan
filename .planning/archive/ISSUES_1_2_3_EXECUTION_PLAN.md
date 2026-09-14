# Execution Plan: Issues #1, #2, #3 Fix with Confidence Scoring

## Context

Three open GitHub issues need fixing in the MucOneSpan pipeline:
- **Issue #1**: Homozygous (same-length) alleles miss mutations because Clair3 is never run on them
- **Issue #2**: Large insertions (e.g., 16bp, 25bp) misclassified because only 60bp windows are probed
- **Issue #3**: False positive mutations from low-quality VCF calls and flanking region indels

The plan in `.planning/ISSUES_1_2_3_PLAN.md` adds a mutation catalog, confidence scoring, template matching, VCF quality filters, VCF validation, anchor-based trim, same-length disambiguation, and bidirectional classification.

## Branch Setup

1. Create branch `fix/issues-1-2-3-confidence-scoring` from `main`
2. Set up `.worktrees/` directory (add to `.gitignore` first)
3. Create worktree and install dev dependencies with `make dev`
4. Verify baseline tests pass with `make test-fast`

## Critical Bug Fix in Plan

**`_apply_mutation()` has an off-by-one error in the plan's delete/delete_insert handling.**

The plan converts 1-based inclusive `end` to 0-based inclusive (`end - 1`), then uses it in Python slice notation which expects exclusive end. This deletes one fewer element than intended.

**Fix:** Use `end = change["end"]` (1-based inclusive end == 0-based exclusive end) instead of `end = change["end"] - 1`:

```python
# WRONG (plan):
end = change["end"] - 1  # 0-based inclusive
del result[start:end]     # misses last element

# CORRECT:
end = change["end"]       # 1-based inclusive → 0-based exclusive
del result[start:end]     # correct range
```

Similarly fix `delete_insert` case. And fix the test `test_del18_31_exact_template_match`: use `x_seq[:17] + x_seq[31:]` (not `x_seq[30:]`) for correct 46bp result.

## Execution Order (Sequential, 10 Tasks)

Tasks are executed in plan order. Each task follows TDD: write tests, verify they fail, implement, verify they pass.

### Task 1: Add mutation catalog to repeat dictionary
- **Files:** `repeats.json`, `config.py`, `test_config.py`
- Add 13 known MUC1 mutations to `repeats.json`
- Add `mutations` and `mutated_sequences` fields to `RepeatDictionary`
- Add `_apply_mutation()` (with fixed indexing) and `_precompute_mutated_sequences()`
- Update `load_repeat_dictionary()` to load and pre-compute
- **Commit:** `feat: add known MUC1 mutation catalog with pre-computed templates`

### Task 2: Add confidence scoring to classify.py
- **Files:** `classify.py`, `test_classify.py`
- Add `"confidence": 1.0` for exact matches in `classify_repeat()`
- Add `"confidence": identity_pct / 100` for non-exact matches
- Add `allele_confidence` and `exact_match_pct` summary to `classify_sequence()`
- **Commit:** `feat: add per-repeat confidence scoring and allele summary metrics`

### Task 3: Exact-match-first probing with mutation templates
- **Files:** `classify.py`, `test_classify.py`
- Add `_probe_sizes_generator()` helper
- Update `classify_repeat()` to check `mutated_sequences` dict
- Update `classify_sequence()` with Phase 1 (exact match across all probe sizes) before Phase 2 (edit distance fallback)
- Build `seq_to_id` reverse map once at top of `classify_sequence()`
- **Commit:** `feat: exact-match-first probing with mutation templates (fixes #2 for known mutations)`

### Task 4: Strengthen VCF quality filters
- **Files:** `calling.py`, `test_calling.py`
- Add `min_qual` and `min_dp` parameters to `filter_vcf()`
- Add `bcftools view -i "QUAL>=X && INFO/DP>=Y"` filter expression
- Update `call_variants_per_allele()` to pass quality params (defaults: `min_qual=15.0`, `min_dp=5`)
- **Commit:** `feat: add VCF quality filters (QUAL, DP) to reduce false positives (fixes #3)`

### Task 5: VCF-backed mutation validation
- **Files:** `classify.py`, `test_classify.py`
- Add `validate_mutations_against_vcf()` function
- Cross-validates detected mutations against VCF variant positions
- Adds `vcf_support` flag and adjusts confidence based on VCF evidence
- **Commit:** `feat: add VCF-backed mutation validation with confidence adjustment`

### Task 6: Anchor-based flanking trim
- **Files:** `consensus.py`, `test_consensus.py`
- Add `_find_anchor()` helper for substring search near expected position
- Add `repeat_dict` parameter to `trim_flanking()` for anchor-based boundary detection
- Fallback to fixed-position trim without `repeat_dict`
- Update `build_consensus_per_allele()` to pass `repeat_dict`
- **Commit:** `feat: anchor-based flanking trim resilient to consensus indels (fixes #3)`

### Task 7: Same-length allele disambiguation (Issue #1)
- **Files:** `alleles.py`, `calling.py`, `consensus.py`, `cli.py`, `test_alleles.py`
- Change `detect_alleles()` to return `same_length=True, homozygous=False` for single-cluster/same-length cases
- Add `parse_vcf_genotypes()` to calling.py
- Add `disambiguate_same_length_alleles()` using Clair3 het genotype detection
- Update `call_variants_per_allele()` to route same_length cases
- **Commit:** `feat: disambiguate same-length alleles using Clair3 het genotypes (fixes #1)`

### Task 8: Bidirectional classification fallback
- **Files:** `classify.py`, `test_classify.py`
- Add `_classify_backward()` helper for 3'-anchored backward pass
- After forward pass stalls, try backward pass anchored on after-repeats
- Gap between forward and backward = mutated region
- **Commit:** `feat: bidirectional classification fallback for novel large mutations (fixes #2)`

### Task 9: Wire features into CLI
- **Files:** `cli.py`
- Pass `repeat_dict` to `build_consensus_per_allele()`
- Add VCF-backed validation after classification
- Add `_parse_vcf_for_validation()` helper
- Display confidence scores in output
- **Commit:** `feat: wire mutation catalog, confidence scoring, and VCF validation into CLI`

### Task 10: Quality checks and final cleanup
- Run `ruff check`, `ruff format`, `mypy`
- Run full test suite with coverage (target >= 70%)
- Fix any issues
- **Commit:** `chore: fix lint, format, and type check issues`

## Verification

After all tasks complete:
1. `make lint` -- no errors
2. `make format` -- no changes
3. `make type-check` -- no errors
4. `make test` -- all pass, coverage >= 70%
5. Review git log for clean commit history

## Key Files Reference

| File | Current Lines | Changes |
|------|--------------|---------|
| `src/muc_one_span/data/repeats/repeats.json` | 62 | Add `mutations` section (~70 lines) |
| `src/muc_one_span/config.py` | ~110 | Add 2 fields, 2 functions, update loader |
| `src/muc_one_span/classify.py` | ~342 | Add confidence, templates, VCF validation, bidirectional |
| `src/muc_one_span/calling.py` | ~314 | Add quality params, genotype parsing, disambiguation |
| `src/muc_one_span/consensus.py` | ~158 | Add anchor-based trim, repeat_dict param |
| `src/muc_one_span/alleles.py` | ~305 | Add same_length flag |
| `src/muc_one_span/cli.py` | ~386 | Wire new features |
