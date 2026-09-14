# Development ablation of minimum coverage

Recorded 2026-09-14. Decision: retain the default `--min-coverage 10`.
Reducing it to 5 recovered one additional exact allele across 154 truth alleles,
with no improvement in complete diploid reconstruction, while normal-sample
annotation false alarms increased from 4 to 7 out of 30 normal inputs. This is
completed development evidence, not a fresh or blinded validation result.

## Evidence and scope

All paths below are relative to ignored
`tests/results/production_validation_20260914/phasing_debug/`:

- `length_thresholds.py` and `length_thresholds.json`: cached mapping-based allele
  detection at coverage 10, 5, and 3; 231 rows for 77 existing inputs.
- `coverage_pipeline.py`: complete pipeline reruns of all 77 inputs at coverage 5.
- `coverage5/inventory.json` and `coverage5/measurements.json`: all expected inputs
  and measured runs, including group, platform, truth location, and exact CLI args.
- `compare_coverage.py` and `coverage5/comparison.json`: paired evaluator results
  for default-10 artifacts in `../development/` and coverage-5 artifacts.

The actual comparison file is inside `coverage5/`, not directly in
`phasing_debug/`. Its four groups are original HiFi (44), original ONT (3), an
existing challenge group named `heldout` (6 HiFi), and existing perturbations
(24 HiFi). Thus the platform totals are **74 HiFi and 3 ONT**, not two equally
sized platform cohorts. The six historical `heldout` cases were inspected during
development and do not constitute unseen final validation here.

The coverage-5 driver used Click's `CliRunner` in process, four threads, the
platform-specific Clair3 model, and explicit `--min-coverage 5`. Baseline CLI args
omit the option and use default 10. The driver ran all inventory entries in order,
not only previous successes. All 77 coverage-5 runs exited 0; their recorded wall
times sum to 474.381 seconds. This is descriptive runtime only: this driver was
not the separately frozen, alternating-order paired timing protocol and does not
provide isolated process RSS/CPU or a controlled performance comparison.

No runs, source edits, driver edits, or data generation were performed while
writing this note. Verification read saved JSON, checked inventory/run identity
sets and paired sample order, and recomputed per-sample metric sums.

## Full pipeline results

Each cell containing an arrow means coverage **10 → 5**. Exact alleles use the
saved evidence-aware sequence metric. Calls include ambiguous reconstructions;
they are not equivalent to complete recovery or a confidently negative result.

| Development group | Inputs / truth alleles | Exact alleles | Exact diploid samples | Calls | Exact repeat counts | Annotation event TP / FP | Supported event TP / FP |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Original HiFi | 44 / 88 | 39 → 39 | 10 → 10 | 44 → 44 | 82 → 83 | 24/14 → 24/14 | 19/0 → 20/0 |
| Original ONT | 3 / 6 | 3 → 3 | 0 → 0 | 3 → 3 | 3 → 3 | 2/0 → 2/0 | 2/0 → 2/0 |
| Historical challenge, HiFi | 6 / 12 | 1 → 1 | 0 → 0 | 6 → 6 | 5 → 5 | 2/0 → 2/0 | 0/0 → 0/0 |
| Perturbations, HiFi | 24 / 48 | 7 → 8 | 0 → 0 | 18 → 24 | 27 → 39 | 7/11 → 10/24 | 7/0 → 7/0 |
| All | 77 / 154 | 50 → 51 | 10 → 10 | 71 → 77 | 117 → 130 | 35/25 → 38/38 | 28/0 → 29/0 |

The truth-event denominator remains 47: 28 original HiFi, 2 ONT, 3 challenge,
and 14 perturbation events. Supported events require the evaluator's exact
sequence-concordant VCF support policy; annotation events include calls without
that support. Neither event count is the number of extra reconstructed alleles.
The inventory has zero cardinality `extra_alleles` at both settings, which does
not establish that every returned allele sequence is correct.

| Platform, including its challenge inputs | Exact alleles | Complete diploid recovery | Calls | Resolved reconstruction statuses |
| --- | ---: | ---: | ---: | ---: |
| HiFi | 47/148 → 48/148 | 10/74 → 10/74 | 68/74 → 74/74 | 14/74 → 13/74 |
| ONT | 3/6 → 3/6 | 0/3 → 0/3 | 3/3 → 3/3 | 0/3 → 0/3 |

Overall exact-allele recovery is 50/154 (32.47%) → 51/154 (33.12%). Complete
diploid recovery remains 10/77 (12.99%). Call rate rises from 71/77 (92.21%) to
77/77 (100%), while six insufficient-evidence results become outputs. Missing
truth alleles decrease from 17 to 2, but most newly returned alleles are not exact.
The truth denominator remains 154 throughout.

The perturbation-only conditional exact-sequence fraction actually falls from
7/32 returned alleles (21.88%) to 8/47 (17.02%), while all-truth exact recovery
increases from 7/48 (14.58%) to 8/48 (16.67%). Both must be shown with call rate;
neither conditional accuracy nor more calls alone determines the preferred default.

## Normal-sample false alarms and unresolved negatives

These counts are sample-level annotation alarms on mutation-negative truth,
recomputed from the saved per-sample `normal_false_positive` fields. The normal
denominator includes every normal input, including unresolved and no-call results.

| Group | Normal inputs | Annotation false alarms | Confident true negatives | Unresolved negatives | Supported false alarms |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original HiFi | 16 | 1 → 1 | 5 → 5 | 10 → 10 | 0 → 0 |
| Original ONT | 1 | 0 → 0 | 0 → 0 | 1 → 1 | 0 → 0 |
| Historical challenge, HiFi | 3 | 0 → 0 | 0 → 0 | 3 → 3 | 0 → 0 |
| Perturbations, HiFi | 10 | 3 → 6 | 1 → 0 | 6 → 4 | 0 → 0 |
| All | 30 | 4 → 7 | 6 → 5 | 20 → 18 | 0 → 0 |

For HiFi alone, this is 4/29 → 7/29 normal inputs with annotation alarms;
ONT remains 0/1, but its one normal input is unresolved at both settings.
Across platforms, annotation alarms rise from 13.33% to 23.33% of the fixed
normal-input denominator. Supported false alarms remain zero, but most normal
inputs are unresolved. It would be incorrect to report 100% normal specificity
from those zero supported alarms: saved supported-specificity assessable
denominators are only 6 at coverage 10 and 5 at coverage 5.

The original HiFi alarm is `sample_gap4_40_44` at both settings. The existing
perturbation alarms at both settings are
`sample_normal_60_80__n40_seed1701`, `sample_normal_60_80__n40_seed1703`, and
`sample_normal_60_80__unique_names`. Coverage 5 additionally alarms on
`sample_normal_60_80__n20_seed1701`, `sample_normal_60_80__n20_seed1702`, and
`sample_normal_60_80__n40_seed1702`. The last changes from a confident negative
to an ambiguous reconstruction with two false annotation events.

## What changed in individual development cases

- `sample_asymmetric_25_140`, and its unique-name perturbation, change candidate
  lengths from 25/28 to 25/140. The long allele length is recovered, but the long
  sequence remains inexact; both samples retain one exact allele and no exact
  diploid reconstruction. The unique-name case also gains a false annotation event.
- `sample_dupa_100_120` gains one supported true mutation event without gaining
  an exact allele or complete diploid reconstruction.
- The three `sample_dupc_60_80__n10_seed1701/1702/1703` inputs previously fail
  the coverage gate. At coverage 5 they emit candidate lengths 60/63 rather than
  the truth 60/80. Seed1702 recovers one exact allele and one supported event;
  the other two each add one true and one false annotation event without exact
  allele recovery.
- `sample_dupc_60_80__n20_seed1703` loses a previously supported true event at 5.
  Thus the perturbation supported-event total staying at 7/14 hides one gain and
  one loss. This is not a uniform improvement.
- The three normal n10 inputs also change from insufficient evidence to ambiguous
  outputs, with no exact allele gain. None becomes a confident normal negative.

## Cached detection at 3 is a separate, incomplete experiment

The cached screen calls `detect_alleles` on existing mapping BAMs for all 77 inputs
at 10, 5, and 3. It does not run consensus, mutation classification, or evaluation
at 3. The 231 saved rows therefore support candidate-length comparisons only.

Across the three settings, 15 inputs change their length result or error; the
other 62 are unchanged. At 10, six n10 perturbations have no contig meeting the
threshold; their maximum observed counts are 6–8 reads. At 5 all six return
candidates. Lowering 5 further to 3 changes candidate lengths in five additional
comparisons: dupC n10 seeds1702/1703, dupC n20 seed1703, normal n10 seed1701,
and normal n20 seed1703 become 60/80. That does not demonstrate recovery of the
corresponding full sequences or absence of false mutation calls.

This parameter controls coverage for length-candidate detection. It is not a
validated minimum support floor for an independently reconstructed phased allele.
In particular, any future phase-support threshold must count original eligible
reads, not a downsampled polishing subset that can turn true 30:2 support into 14:1.

## Decision and limits

Retain default 10 for the freeze. Coverage 5 improves some candidate lengths and
call rate, but provides only one extra exact allele, no complete-pair gain, three
additional normal-sample alarms, and a supported-event regression hidden by an
offsetting gain. Its zero supported false-event count does not erase the larger
annotation burden or resolve the normal-negative denominator.

These existing development inputs and perturbations are correlated and were
already inspected. There are only three ONT inputs and no ONT depth-perturbation
cohort. No repeated timing experiment, independent final-seed validation, or full
pipeline coverage-3 experiment is claimed. Saved comparisons were generated
against development artifacts, not a pair of immutable source freezes; incidental
pipeline nondeterminism or development changes cannot be excluded by this driver.
The results justify declining a default change, not a universal optimality claim
for coverage 10.

## Evidence fingerprints

SHA256 values checked when preparing this note:

| File relative to `phasing_debug/` | SHA256 |
| --- | --- |
| `coverage_pipeline.py` | `aab88a99e1b8bfa2628495279d01611d39d53062a82da5838334c71a42262440` |
| `compare_coverage.py` | `091664ee4e1b48bc2ed3ac0abde454fc08f6435e914156424b67d7c91450893f` |
| `length_thresholds.py` | `0431cb5f2fa9851e6fb7e8bc15532d9ee1fb53e39d3fa9a26752a62b8337aee6` |
| `length_thresholds.json` | `4ff3e062045c8297dbffd705ca63a25b1411a7f6ebc744d2a9a93d345d41d0bb` |
| `coverage5/inventory.json` | `66182019aba7b89c9d0d6317be714fe7e25585f7f1c2eff8d94f3d0b007364ac` |
| `coverage5/measurements.json` | `d91adae633a57753201ec503be1f387b7fe50578fceabe8ad0b6b0989a1f5be3` |
| `coverage5/comparison.json` | `40a0246a5eebe6014fbdfd129750f6cb7bfd3c68a8e5cb63d9802328c702e44e` |
