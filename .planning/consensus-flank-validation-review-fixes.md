# Runtime review F1/F5: known flank extents and stale phase metadata

2026-09-14. Follow-up to the actual external review stored in
`.planning/reviews/runtime-configuration-response.json`. These changes address
accepted nondefault configuration and metadata correctness; they do not promote
the rejected length-model prototypes or change scientific calling thresholds.

## F1: reject inconsistent reference/trim extents

The review reproducer has four-base left/right dictionary flanks but requests
six bases. The ladder previously truncated to four, while the trimmer still used
six-base boundary coordinates. At zero tolerance the 12-base VNTR
`GTTGACGTCACA` silently became eight bases, `TGACGTCA`.

`ConsensusSettings.validate_flanks(left, right, *, flank_length=None)` now
requires the effective nonnegative integer extent to fit both known dictionary
flank strings. An explicit argument overrides the setting. A zero extent allows
empty strings. With a supplied dictionary an empty flank is known to have length
zero; a positive requested extent therefore fails rather than silently shrinking.

The shared method runs in `build_contig`, in `generate_ladder_fasta` before output
creation, in `trim_flanking` before reading/writing FASTA, and in
`build_consensus_per_allele` before output-directory creation or external tools.
The coordinator wires the same method into pipeline preflight. Without a supplied
dictionary, the trimmer retains the explicit positional-flank contract because
actual reference flank lengths are unavailable.

The selected contract is explicit rejection, not implicit rescaling to different
left/right lengths. The error identifies the requested extent and both available
lengths. Valid explicit flank overrides retain precedence over oversized defaults.

## F5: clear stale rebuilt-candidate phase status

After rebuilding a single output with `independent_haplotype_evidence=False`,
`vntr_phase_status` is now explicitly `unresolved`. A stale
`distinct_genotype_candidates` value can no longer survive from a previous run.
The existing identical-sequence status and positive independent-evidence branch
are preserved. FASTA sequence generation is unchanged by this metadata fix.

## Executed evidence

Before implementation, **six failures and 23 passing controls** reproduced:
three left/right extent combinations, silent VNTR truncation, external tools
being reached before rejecting the invalid extent, and stale phase metadata.
After implementation:

- **147 focused unit tests passed**, zero skips.
- **20 real bcftools/samtools integration tests passed**, zero skips.
- Focused Ruff and mypy passed. Initial Ruff findings concerned test regex/raw
  strings and nested context-manager formatting; these were fixed and rechecked.
- All **137 available cached full consensuses**, selected from the explicit
  77-input development inventory, produce FASTA bytes and trim coordinates
  identical to their cached outputs. The 17 absent allele slots remain absent
  and are not counted as parity successes.
- The complete default ladder remains byte-identical with SHA-256
  `2349a88c8435790aff3deb1c446901888ba9f862d7d4e9f00ed59bd8e8272e34`.

```sh
uv run --locked --all-extras pytest -q tests/unit/test_reference_settings.py tests/unit/test_consensus_settings.py tests/unit/test_consensus.py tests/unit/test_ladder.py tests/unit/test_runtime_settings.py --no-cov
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH uv run --locked --all-extras pytest -q tests/integration/test_phase_consensus.py tests/integration/test_variant_concordance.py --no-cov -rs
uv run --locked --all-extras ruff check src/muc_one_span/settings.py src/muc_one_span/consensus.py src/muc_one_span/ladder.py tests/unit/test_reference_settings.py tests/unit/test_consensus_settings.py
uv run --locked --all-extras mypy src/muc_one_span/settings.py src/muc_one_span/consensus.py src/muc_one_span/ladder.py
```

Cached parity used the same explicit inventory and temporary-output procedure
recorded in `consensus-anchor-config-fixes.md`. No stored results were overwritten,
no final seed was generated, no new scientific pipeline was called, and no commit
or tag was created by this subtask. The earlier anchor correction's documented
behavior difference on synthetic flank-insertion cases remains intentional;
these F1/F5 changes add no further sequence difference for valid configurations.
