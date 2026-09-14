# Classification runtime-settings evidence

Date: 2026-09-14. Scope: classification and heuristic confidence sections of
the schema-1 runtime configuration. No final reserved seed was read or generated.

## Interface and precedence

`classify_repeat` and `classify_sequence` accept keyword-only
`ClassificationSettings | None`; `validate_mutations_against_vcf` and the
compatibility-exported `_qual_to_confidence` accept `ConfidenceSettings | None`.
Omission uses the corresponding immutable object in `DEFAULT_SETTINGS`.
Forward/backward fitting receives the already resolved classification object, so
validated configuration is not rebuilt inside sequence-search loops.

The configured classification fields control their pre-existing decisions:

- `max_indel_probe`: variable-length probe interval;
- `max_fit_edit_distance`: strict-fit acceptance and backward recovery limit;
- `novel_repeat_edit_distance`: variant versus novel-repeat label;
- `minimum_unit_fraction`: minimum processable/probed repeat fraction;
- `early_stop_edit_distance`: alternate-window search stopping point; and
- `strict_segmentation`: default segmentation policy.

The explicit `strict_segmentation` argument now uses a `None` sentinel and
overrides the settings object when supplied. The private forward helper retains
its existing explicit `max_indel_probe` argument as an override for compatibility.

Confidence settings control the low/high QUAL breakpoints, below/low/high
weights, replay-verified-absence weight, terminal-repeat interval, and boundary
weight. Existing explicit `boundary_repeats` and `boundary_penalty` arguments use
`None` sentinels and override the settings object. N1 semantics are unchanged:
only replay status `absent` receives `absent_weight`; unavailable or ambiguous
projection applies no VCF-evidence multiplier, while the independent configured
boundary multiplier can still apply.

No codon, coordinate, edit-distance recurrence, tie, genotype, sequence, event,
or support-status rule became configurable. No mutable runtime global was added.

## Test-first evidence

The new behavior suite was first run before `settings.py` was available and
failed at collection. After the immutable class definitions became available,
the pre-configuration source snapshot was isolated with those definitions and
the same seven tests produced **7 failures**, each on the absent `settings`
keyword. This distinguishes missing integration from missing shared types.

`tests/unit/test_classification_settings.py` now observes every override:

- novelty threshold changes only the label at fixed edit distance;
- strict fit threshold stops a one-edit repeat, and explicit
  `strict_segmentation=False` overrides a strict configured default;
- minimum unit fraction changes whether a 30-base terminal unit is processed;
- zero indel probe width prevents exact variable-length template selection;
- early-stop distance changes which alternate probe wins;
- all QUAL breakpoints/weights change the continuous mapping; and
- absence/boundary weights change confidence, with both legacy boundary
  arguments tested as explicit overrides.

Focused final command:

```text
uv run --locked --all-extras pytest tests/unit/test_classify.py \
  tests/unit/test_classification_evidence.py \
  tests/unit/test_classification_settings.py \
  tests/unit/test_repeat_alignment.py tests/unit/test_variant_support.py \
  -q --no-cov
82 passed
```

## Isolated default parity

The immutable pre-configuration source is
`tests/results/production_validation_20260914/configuration/before/src`; its
sibling `manifest.json` records the review-fixed source/resource hashes. Separate
Python processes put either that snapshot or the current worktree first on
`PYTHONPATH`. They classified these six existing cached, trimmed consensuses:

- HiFi difficult `sample_bench_5003/allele_2`;
- HiFi `sample_dupc_40_50/allele_1`;
- HiFi same-length `sample_homozygous_60_60/allele_1`;
- HiFi long `sample_long_120_140/allele_2`;
- HiFi normal `sample_normal_60_80/allele_1`; and
- ONT `sample_ont_dupa_60_80/allele_1`.

The comparison also includes every exact dictionary repeat, a one-substitution
case for every repeat, every mutation template, default QUAL weights at
0/4.9/5/12.5/20/100, and exact/absent/unavailable/localization-ambiguous support
confidence cases including the unchanged boundary penalty.

The first comparison found only a serialized type mismatch: the settings
dataclass default `weight_high=1` returned integer `1`, while the compatibility
helper historically returned float `1.0`. Returning a float preserves the API.
After that repair, `before.json` and `after.json` are byte identical with SHA-256:

```text
be5beba006dc6588b279e03b41066e2dc5b5410e0147e045d7a87fd486e58154
```

Ignored artifacts are retained under
`tests/results/production_validation_20260914/configuration/classification_default_parity/`.
This is cached development equivalence, not fresh validation or scientific
validation of nondefault choices.

## Static checks and limits

Focused Ruff, formatting, configured mypy, `git diff --check`, and file-size
checks pass. `classify.py` remains below the 649-line maximum, so no private helper
was moved and its existing import paths remain intact. Full configuration,
integration, and repository gates remain coordinator work.
