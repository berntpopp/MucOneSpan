# Milestone 2 low-coverage length evidence accounting

Date: 2026-09-14. Scope: Claude milestone-2 finding F2 only.

## Reproduction

Three focused tests were added before implementation. All three failed with
`KeyError: 'length_selection_evidence'`, while the four pre-existing allele
evidence tests passed. The cases require:

- a contig with `0 < idxstats alignment records < min_coverage` to remain
  visible after the existing coverage filter;
- primary records on two excluded contigs to be counted by one targeted
  `samtools view` query, including two separate primary records with the same
  QNAME and excluding secondary and unmapped records;
- a third passing cluster to remain visible even though the existing caller
  selects only the first two; and
- an empty excluded set to report zero records without requiring a BAM.

## Bounded diagnostic behavior

Each returned allele now contains a nested `length_selection_evidence` object.
It records the configured minimum coverage, identifies the unit as
`alignment_records_not_molecules`, keeps `molecule_count` null, lists every
sub-threshold nonzero contig and its idxstats alignment-record count, and lists
passing clusters beyond the two selected by the existing algorithm.

When an existing BAM is supplied, one narrow `samtools view BAM CONTIG...`
invocation covers all excluded contigs. Each mapped primary SAM record is counted
individually; QNAMEs are not deduplicated. Secondary, supplementary and unmapped
records are excluded from that primary count. Without a BAM, primary counts are
null when excluded contigs exist. An empty excluded set has zero alignment and
primary records.

The evidence is nested under `allele_1` and `allele_2`; no new top-level mapping
is introduced. Existing length values, cluster order, first-two selection,
single-candidate duplication, same-length status, and phase inputs are unchanged.
The new metadata is not read by a cutoff, support gate, phase rule, or caller
decision. It does not establish molecule counts, joint read assignment, or a
one-versus-two genotype. Experimental read phasing remains disabled by default.

## Validation

- Red: `uv run --locked --all-extras pytest tests/unit/test_allele_evidence.py
  -q --no-cov` — **3 failed, 4 passed**, all failures on the absent nested field.
- Focused after implementation: same command — **7 passed**.
- Complete allele unit scope after preserving the nonexistent-BAM mock boundary:
  `uv run --locked --all-extras pytest tests/unit/test_alleles.py
  tests/unit/test_allele_evidence.py -q --no-cov` — **39 passed**.
- Focused Ruff and formatting checks passed; configured mypy reported no issues
  in `alleles.py`; `git diff --check` passed. The source and focused test files
  are 574 and 134 lines. Repository-wide final gates remain coordinator work.

Exact regressions:

- `tests/unit/test_allele_evidence.py::test_length_selection_evidence_retains_dropped_and_unselected_candidates`
- `tests/unit/test_allele_evidence.py::test_excluded_primary_records_count_colliding_qnames_in_one_query`
- `tests/unit/test_allele_evidence.py::test_empty_excluded_evidence_has_zero_counts_without_bam`
