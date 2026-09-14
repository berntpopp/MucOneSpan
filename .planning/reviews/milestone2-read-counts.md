# Milestone 2: selected read-count metadata

Date: 2026-09-14. Bounded follow-up to Claude review F1/F8, restricted to
`src/muc_one_span/read_phasing.py` and `tests/unit/test_read_phasing.py`.

The review correctly identified that aggregate read-list size hides a 10:1 or
12:1 phase split. The coordinator authorized additive selected-read metadata,
with no support cutoff and no acceptance, default, or phasing algorithm changes.
The existing experimental read-phasing option remains disabled by default.

## Implemented metadata contract

The existing `phasing_assigned_records` total is preserved. After existing
read-list validation, metadata additionally reports:

```json
{
  "primary_records": 30,
  "phasing_assigned_records": 11,
  "phasing_selected_records_by_phase_set": {
    "1": {"0": 10, "1": 1}
  },
  "phasing_read_count_scope": "internal_selection_not_all_primary_records"
}
```

Outer keys are phase-set strings; inner keys retain WhatsHap's zero-based
haplotype labels. Each represented phase set includes both labels, with zero
when one has no selected rows. Counts are accumulated only after the existing
identity, sample, source, haplotype, position, and variant-count checks succeed.
Malformed assignments still fail. Empty valid unresolved lists yield `{}`;
attempts without a validated list retain `null` for the new count map.

These counts describe the tool's output read list after internal selection and
possible downsampling. They are not all-original-read support, original molecule
counts, or evidence that an unrepresented biological allele is absent. A true
30:2 source mixture can become 14:1 internally; no minimum support rule may be
calibrated against these numbers as though they were full-data assignments.
Different phase sets remain separate because their local haplotype labels need
not represent the same biological haplotypes.

The parser now returns the existing validated total plus the grouped counts.
The tool command, candidate VCF identity check, common-phase-set acceptance
policy, source-map handling, and consensus selection remain unchanged. This
addresses observability in F1/F8, not the broader scientific phase-acceptance
concern. No improvement in allele recovery or false-positive rate is claimed.

## Validation

Tests were changed before production code. The first focused run produced
**4 failed, 26 passed**, all four failures being the missing count-map key.
After implementation, all 30 read-phasing unit cases passed. New checks cover
30 primary records versus 11 selected rows, both 10:1 and 11:0 selected splits,
an empty unresolved list, and disconnected phase sets with independent counts.
The imbalanced fixtures explicitly retain the existing acceptance behavior;
they do not assert a new support threshold.

Commands and results:

- `uv run --locked --all-extras pytest tests/unit/test_read_phasing.py --no-cov -q`:
  30 passed.
- `uv run --locked --all-extras pytest tests/unit/test_read_phasing.py tests/unit/test_phasing.py tests/unit/test_calling.py --no-cov -q`:
  65 passed, including the calling/default-path checks.
- `uv run --locked --all-extras ruff check src/muc_one_span/read_phasing.py tests/unit/test_read_phasing.py`:
  passed.
- `uv run --locked --all-extras ruff format --check src/muc_one_span/read_phasing.py tests/unit/test_read_phasing.py`:
  both files already formatted.
- `uv run --locked --all-extras mypy src/muc_one_span/read_phasing.py`:
  passed.
- Initial integration invocation without the configured tool PATH: 17 skipped
  because bcftools was absent from that PATH; this was not counted as validation.
- `PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH uv run --locked --all-extras pytest tests/integration/test_read_phasing.py --no-cov -q -rs`:
  **17 passed** using the existing temporary real-tool fixtures (3.61 seconds).

The coordinator runs repository-wide final checks. No commit, final validation
data generation, final caller run, experiment-driver change, or default change
was performed in this follow-up.
