# MucSim-Bench dev pilot (Task 12)

Run on 2026-09-25; the file name keeps the plan date (2026-09-24). Branch `feat/benchsim`.
Data and results live outside Git in `../MucOneSpan-bench-data` (regenerable).
This file holds only aggregate numbers: no reads, truth sequences, per-library
values or patient data.

> **Superseded settings (Task 12 fix round 1).** After this pilot, all
> benchmark tunables moved to `benchsim.bench_config`. Two defaults changed:
> - Event targets are now clamped to units 6..L-4, from the bundled reference
>   layout (head 1-5, tail 6-9), instead of 5..L-5.
> - Smear 0.5 moved to the `stress` split; dev/val/test now use 0.05 and 0.25.
>
> The numbers below come from the earlier settings. The pilot is regenerated
> in the next round.

## Engine under test

The pilot runs the **`ladder` engine of this branch**, which is the
MucOneSpan v0.15.1 caller. The branch predates the Phase 0 clinical-safety fixes.
These numbers are therefore the **pre-P0 baseline**, not the current caller.
A v0.16.0 comparison runs after P0 merges (ledger ruling).

## Setup

| Item | Value |
| --- | --- |
| Split | `dev`, `--n 30` per profile, public salt (default) = 90 designs |
| MucOneUp | 0.45.0 (local checkout `config.json`, built-in read profiles) |
| Read simulators | pbsim3 and ccs via MucOneUp `tools` (conda env) |
| Caller tools | minimap2 2.28-r1209, samtools 1.15.1, bcftools 1.17, Clair3 v1.0.10 |
| Clair3 models | ONT `r1041_e82_400bps_sup_v500`; HiFi Clair3 bundled `hifi` model |
| `--structure-pool` | not supplied: `real_derived` designs fell back to Markov (recorded per case) |
| `--flank-fasta` | not supplied |
| Parallelism | generate `--jobs 6`; run `--jobs 4 --threads 3` |

`CLAIR3_MODEL_ONT` / `CLAIR3_MODEL_HIFI` were not set, so both models were passed by flag.

## Commands and runtimes (wall clock)

| Step | Command | Wall | Result |
| --- | --- | --- | --- |
| design | `benchsim design --split dev --n 30` | < 1 s | 90 designs |
| generate | `benchsim generate --designs designs_dev.jsonl --jobs 6 ...` | 7 min 14 s (118 CPU-min) | 90 ok, 0 failed |
| realism | `benchsim realism --split dev` | 2 min 56 s | 0 case failures; checks below |
| run | `benchsim run --engines ladder --threads 3 --jobs 4 ...` | 20 min 57 s (133 CPU-min) | 84 completed, 6 insufficient_evidence |
| evaluate | `benchsim evaluate --split dev --engines ladder` | 4 s | exit 0 |
| report | `benchsim report --split dev --baseline ladder` | < 1 s | tables only (no candidate) |

Per-case caller wall time (s): HiFi amplicon median 22 (max 360); ONT amplicon
median 12 (max 959; the depth-2000 cases dominate); ONT genomic median 25 (max 66).

## First attempt and the bug it exposed

The first generation, before commit `efc6227`, gave 83 ok, 5 `design_invalid` and
2 `generation_failed` cases out of 90. All 7 failures had an event target on a
conserved unit:

- Both `generation_failed` cases (HiFi amplicon) targeted **unit 1**. The event rewrote
  the unit that holds part of the forward amplicon primer site. The amplicon
  primer search then found no forward site (`primer sites not unique ([], [...])`).
- All 5 `design_invalid` cases (1 HiFi, 2 ONT amplicon, 2 ONT genomic) targeted the
  **last unit** of the allele (`first10` or `last10` of short alleles). Truth
  validation then failed with "reconstructed sequence/flanks do not match FASTA".

The fix (`efc6227`, tested first) clamps event targets to units 5..L-5. It
matches the rare-unit rule, which never touches units 1-4 or the last 5.
Clamping keeps the random draw sequence, so only the targets of 22 of 90
designs changed. The pre-fix data are archived in
`../MucOneSpan-bench-data/_prefix_2026-09-25/` (outside Git). All numbers
below come from the post-fix regeneration.

A second small fix (`2400d25`): `report` without `--candidate` printed
"Decision: NOT ADOPTED". It now says that no decision was evaluated.

## Generation outcome

- 90/90 `ok`; 0 `design_invalid`; 0 `generation_failed`.
- `amount_capped` (minor-allele share floored at 0.05 for template sizing, `c47dc49`): HiFi 4, ONT amplicon 3, genomic 0.

### Realized versus target depth (ledger ruling: within 30% for most cases)

Realized spanning depth per allele compared with the design depth:

| Profile | Alleles within +/-30% | Cases, both alleles within | Minor allele within | Median ratio (q10-q90) |
| --- | --- | --- | --- | --- |
| hifi_amplicon | 47/60 (78%) | 21/30 | 24/30 | 0.94 (0.75-2.45) |
| ont_amplicon_r10 | 50/60 (83%) | 21/30 | 28/30 | 1.02 (0.90-2.80) |
| ont_genomic_targeted | 36/60 (60%) | 14/30 | 21/30 | 1.07 (0.67-1.50) |
| all | 133/180 (74%) | 56/90 (62%) | 73/90 (81%) | |

Verdict: **most alleles (74%) and most minor alleles (81%) fall within 30%.**
The ruling's criterion is met at the allele level. It is not met for "both alleles" in genomic cases. The deviations follow the design:

- **Amplicon over-target on the major allele.** Template counts are sized so
  the minor (PCR-disadvantaged) allele reaches the target. The major allele
  then overshoots, up to 2.8x at q90.
- **Genomic depth 3 and 6.** 0/10 and 2/10 alleles fall within 30%. With counts
  of 2-8 reads, one read is already 17-33%. Poisson sampling makes a 30%
  window unattainable at these depths. From depth 10 up, 34/40 genomic alleles
  are within 30%.

Realized depth is stored in every `case.json`, so analyses can stratify on it.

## Realism report (dev, 30 cases per profile)

| Profile | Checks passed | Failing checks |
| --- | --- | --- |
| ont_amplicon_r10 | 10/21 | insertion + and all; deletion - and all; error - and all; span_off_gt1unit_frac, span_between_alleles_frac, span_below_short_frac; span-offset JS lt55u and ge55u |
| ont_genomic_targeted | 9/15 | mismatch -; insertion - and all; deletion + and -; error all |
| hifi_amplicon | no verdict | no public target section (uncalibrated profile) |

Aggregate values (sim versus public PRJEB92208 target):

- ONT amplicon error rate, all strands: 0.0263 versus 0.0214 (+23%; the tolerance is 20%).
  Minus strand 0.0217 versus 0.0172. Deletions on the minus strand: 0.0080 versus 0.0049.
  Mismatch rates pass. C7 homopolymer accuracy passes on both strands (0.520/0.876 versus 0.517/0.888).
- ONT amplicon span-offset histogram: Jensen-Shannon distance 0.210 (< 55 units)
  and 0.271 (>= 55 units), against a limit of 0.1. This is the **known gap**
  (Task 9 probes: 0.24-0.45).
- ONT amplicon range checks: the medians sit close to the real medians
  (off-by->1-unit 0.280 versus 0.270; below-short 0.258 versus 0.236).
  The checks fail because the pooled designs span a wider range than real
  libraries. Minimum 0.0 comes from low-smear designs, maximum 0.686 from
  `smear 0.5` designs, against a real range of 0.109-0.537. The
  between-alleles fraction is 0 in most designs (real median 0.023).
  The stratified artefact levels therefore exceed the observed real range at
  both ends. Part of this is intentional stress, but the `smear 0.5` level lies
  above every real library.
- ONT amplicon allelic-ratio slope: -0.061 versus -0.056 per unit (pass).
- ONT genomic error rate, all strands: 0.0050 versus 0.0040 (+25%). Insertions
  0.0016 versus 0.0011; deletions pass in pooled form but not per strand. This is
  the **known gap**: genomic error rates sit above the 2-library WGS target.
  Spanning fraction passes (median 0.316 versus 0.420).

Per the benchmark guide, these are sim-to-real gaps to report with results,
not caller-tuning targets. Upstream calibration (MucOneUp profiles) remains a
follow-up.

## Caller results (ladder = v0.15.1, pre-P0)

Metrics are pooled, with 95% cluster-bootstrap intervals over designs:

| Profile | Per-allele exact | Case exact | Critical FN | FP on normal (0/11 each) | Inconclusive |
| --- | --- | --- | --- | --- | --- |
| hifi_amplicon | 7/60 = 0.117 [0.033, 0.217] | 2/30 | 1/30 | 0/11 | 25/30 |
| ont_amplicon_r10 | 4/60 = 0.067 [0.017, 0.133] | 0/30 | 2/30 | 0/11 | 25/30 |
| ont_genomic_targeted | 7/60 = 0.117 [0.033, 0.217] | 2/30 | 3/30 | 0/11 | 23/30 |
| all | 18/180 = 0.100 [0.050, 0.156] | 4/90 | 6/90 | 0/33 | 73/90 |

- Event recall 7/57 (0.12); supported-event TP 2, supported-event FP 1.
  dupC recall is 0/2 in every profile.
- Clinical confusion (all profiles): pathogenic truths gave 3 PATHOGENIC,
  48 INCONCLUSIVE and 6 NO_CALL calls. Normal truths gave 8
  NO_PATHOGENIC_VARIANT_DETECTED and 25 INCONCLUSIVE calls. There were no
  false PATHOGENIC calls.
- The 6 critical false negatives are exactly the 6 `insufficient_evidence` runs
  ("No contig has >= 10 mapped reads"). All 6 have design depth 3 or 5,
  which is below the caller's 10-read minimum per contig: HiFi 1, ONT
  amplicon 2, ONT genomic 3.
- Predicted allele count: 62/90 cases called 2 alleles, 22 called 1 and 6 called 0.
  Missing alleles total 34/180. Of the 22 single-allele calls, 17 are
  `0_identical`, `0_different` or delta-1 designs, where alleles of the same
  or nearly the same length are merged.
- Exactness concentrates in `>20` length differences (per-allele 0.25-0.63) and
  is 0 in every delta class <= 5 for the amplicon profiles.
- Inspection of the depth >= 150 amplicon cases shows the main failure. The
  ladder picks smear or chimera products as alleles: for example, called 100
  versus truth 42 repeats, or 22 versus 109. It also reports many ambiguous
  units (`?X`) under the `poor` error level.
  This is consistent with the caller weaknesses the hybrid engine and P0
  address. It is not evidence of a scoring fault: truth repeat counts and
  structures in `evaluation.json` match the design.
- Earlier MucOneUp-only validation reported 22-30% HiFi diploid exact
  (`docs/guides/validation-results.md`). That validation had no artefact
  layers. The drop to 7% here shows the effect of realistic artefacts.

The per-factor failure atlas is in `report.md` (not committed). Its main
signals: case failure is 100% for smear 0.05 and 0.5, and 87% for 0.25; the
atlas for PCR, chimera and error level is flat (93-98%). The atlas cannot
separate factors at n=90.

## Concerns and follow-ups

1. Stratified artefact levels (`smear 0.5`, and the low extreme) fall outside the
   real PRJEB92208 range. Decide whether `smear 0.5` stays as a stress level
   or moves to the `stress` split.
2. The amplicon span-offset shape (JS 0.21-0.27) and the genomic and amplicon
   indel rates (+20-60% on some strands) remain uncalibrated. This work belongs
   upstream in MucOneUp.
3. The HiFi profile is uncalibrated, so its numbers show relative engine
   behaviour only.
4. `real_derived` ran as Markov (no structure pool): `composition_effective` is
   Markov in 84 cases and rare_units in 6.
5. Benign stratum is empty until issue #69.
6. The 6 critical false negatives are all minimum-coverage aborts at design
   depth 3-5. The P0 and hybrid comparison should report them separately from
   miscalls.

## v0.16.0 ladder on dev v2 (Task 12b)

Run on 2026-09-24 (system clock) on `feat/benchsim` after merging `origin/main` (v0.16.0,
merge commit `4dd3ccc`). The caller is `muconespan 0.16.0`, which has the P0
clinical-safety gates. Data are in `../MucOneSpan-bench-data/v2` (outside Git).
The older `dev` and `_prefix_2026-09-25` directories were not touched.

### Setup and commands

The settings are the current bench-config defaults, SHA-256 `032ba4f1…4a7d94`.
The same hash is recorded in every `case.json` and in `report.json`. Tools,
models and MucOneUp are as in the first pilot.

| Step | Command (out-root `v2`) | Wall | Result |
| --- | --- | --- | --- |
| design | `design --split dev --n 30` | < 1 s | 90 designs |
| generate | `generate --jobs 4` | 9 min 01 s (112 CPU-min) | 90 ok, 0 failed |
| realism | `realism --split dev` | 1 min 41 s | 0 case failures |
| run | `run --engines ladder --threads 3 --jobs 4` | 14 min 35 s (112 CPU-min) | 84 completed, 6 insufficient_evidence |
| evaluate | `evaluate --split dev --engines ladder` | 3 s | exit 0 |
| report | `report --split dev --baseline ladder` | < 1 s | tables and reason atlas |

Median per-case caller wall time is 11 s for HiFi amplicon (max 273 s),
12 s for ONT amplicon (max 652 s) and 8 s for ONT genomic (max 17 s).

Generation: 23/90 designs had a clamped event target. `amount_capped` applied
to 4 HiFi and 3 ONT amplicon cases. `real_derived` again fell back to Markov
(84 Markov, 6 rare-unit cases). 138/180 alleles (77%) reached within 30% of the
target depth.

Realism is unchanged in kind:
- ONT amplicon passes 10/21 checks. The span-offset JS distance is 0.209/0.247,
  against a limit of 0.1.
- ONT genomic passes 10/15 checks.
- The amplicon range checks still fail, although smear 0.5 is gone. The
  known indel-rate and span-shape gaps remain.

### Side by side with the v0.15.1 baseline (indicative, not paired)

The designs changed between the two runs, so the same design IDs are not the
same cases:
- event targets moved to units 6..L-4;
- dev smear levels are now {0.05, 0.25};
- the settings hash differs.

The comparison is therefore indicative only. No paired test applies.

| Metric (all profiles) | v0.15.1 (dev, first pilot) | v0.16.0 (dev v2) |
| --- | --- | --- |
| Per-allele exact | 18/180 = 0.100 [0.050, 0.156] | 22/180 = 0.122 [0.067, 0.178] |
| Case exact | 4/90 | 5/90 |
| Event recall | 7/57 | 12/57 |
| FP PATHOGENIC on normal | 0/33 | 0/33 |
| Critical FN | 6/90 (all insufficient_evidence) | 6/90 (all insufficient_evidence, design depth 3-5) |
| INCONCLUSIVE | 73/90 | **81/90** |
| NEGATIVE on normal truths | 8/33 | 1/33 |
| PATHOGENIC on pathogenic truths | 3/57 | 2/57 |
| Cases with 2 / 1 / 0 predicted alleles | 62 / 22 / 6 | 56 / 28 / 6 |

Per profile (v0.16.0): per-allele exact is 0.117 for HiFi, 0.133 for ONT
amplicon and 0.117 for ONT genomic. INCONCLUSIVE is 28/30, 28/30 and 25/30.

INCONCLUSIVE rose from 73 to 81. The v0.16.0 gates explain this. A negative
now needs:
- resolved allele selection;
- a consistent length;
- adequate per-allele primary depth;
- no heterozygous call left inside a length partition.

Normal-truth NEGATIVE calls dropped from 8 to 1. There are still no false
PATHOGENIC calls, and the critical false negatives are still only the 6
minimum-coverage aborts. The gates work as designed; the reconstruction did
not improve.

### INCONCLUSIVE reason atlas (v0.16.0)

Atlas settings are the defaults:
- decisions: INCONCLUSIVE;
- expected when the split is `stress` or the lowest realized allele depth is
  below 30 (the caller's per-allele gate).

Every case counts once per reason.

| Profile | INCONCLUSIVE | Expected (depth) | Resolvable | Depth unknown |
| --- | --- | --- | --- | --- |
| hifi_amplicon | 28/30 | 14 | 14 | 0 |
| ont_amplicon_r10 | 28/30 | 12 | 16 | 0 |
| ont_genomic_targeted | 25/30 | 17 | 8 | 0 |
| all | 81/90 | 43 | 38 | 0 |

Of the 38 resolvable cases, 22 have a pathogenic truth and 16 a normal truth.
Every case has a recorded realized depth, so none is depth-unknown.

Top causes (cases out of 81). Keys are verbatim from `report.json`
(fix round 1 re-render; the counts did not change):

| # | Reason key | Cases | Expected / resolvable / depth unknown |
| --- | --- | --- | --- |
| 1 | `evaluator: ambiguous_reconstruction` | 75 | 39 / 36 / 0 |
| 1 | `gate: reconstruction incomplete; independent biological haplotype evidence not established` | 75 | 39 / 36 / 0 |
| 3 | `evaluator: iupac_bases` | 71 | 35 / 36 / 0 |
| 4 | `gate: high number of ambiguous consensus bases (#) detected` | 60 | 27 / 33 / 0 |
| 5 | `gate: # primary alignments, below the per-allele depth gate (#)` | 50 | 33 / 17 / 0 |
| 6 | `gate: observed sequence variant (<variant>) is inconclusive: localization ambiguous` | 49 | 24 / 25 / 0 |
| 7 | `gate: observed sequence variant (<variant>) is inconclusive: no explicit sequence-level support (localization_ambiguous)` | 47 | 23 / 24 / 0 |
| 8 | `gate: heterozygous call left within the length-partitioned allele; consensus uses unresolved (iupac) selection` | 46 | 19 / 27 / 0 |
| 9 | `gate: observed sequence variant (<variant>) is inconclusive: event identity not established (no exact dictionary template)` | 35 | 18 / 17 / 0 |
| 10 | `gate: reported length # differs from the consensus contig length #` | 31 | 9 / 22 / 0 |
| 11 | `gate: observed sequence variant (<variant>) is inconclusive: carrying allele is below the per-allele depth gate` | 30 | 23 / 7 / 0 |
| 12 | `evaluator: missing_allele` | 28 | 20 / 8 / 0 |
| 12 | `evaluator: unresolved_allele_alias` | 28 | 20 / 8 / 0 |
| 14 | `gate: allele selection unresolved (unresolved_secondary_mode; secondary mode fraction #)` | 24 | 7 / 17 / 0 |
| 15 | `gate: allele selection unresolved (unresolved_unselected_clusters; secondary mode fraction #)` | 13 | 6 / 7 / 0 |

Readings:

- **Reconstruction dominates.** 36 of the 38 resolvable cases have an IUPAC
  consensus and are `ambiguous_reconstruction`. The other 2 fail only the
  per-allele depth gate.
- **The per-allele depth gate hits resolvable cases with high realized depth.**
  "Resolvable" uses the simulator's per-allele spanning depth (truth), not the
  caller's primary-record count.
  17 resolvable cases have fewer than 30 primary alignments on an allele,
  even though their lowest realized spanning depth is 30-1984. The ladder
  loses reads during allele partitioning or selection; the sample itself has
  enough reads.
- **Missing alleles are mostly expected-depth cases.** 20 of the 28 are
  expected. 18 of the 28 are `0_identical` or `0_different` designs, and most
  of the rest are low-depth designs.
- The expected half (43, all through the depth condition) can only become
  definitive with more reads. Gates are not relaxed.
- The v0.15.1 baseline shows the same evaluator picture: 73 INCONCLUSIVE, all
  `ambiguous_reconstruction`, 71 with IUPAC bases, 22 missing an allele.
  I reproduced it offline from the archived result directories. Its caller
  gate reasons cannot be recovered from the v0.16.0 code, so they are not
  compared.

Full reason × profile × stratum tables (depth, smear, chimera, delta class,
event position) are in `v2/results/dev/report.md` (not committed). At n=90 they
do not separate factors beyond the depth effect.

### Which causes the hybrid engine is designed to address

These mappings follow `2026-09-23-hybrid-engine-spec.md`.

| Cause (atlas) | Hybrid mechanism |
| --- | --- |
| IUPAC consensus, `ambiguous_reconstruction`, ambiguous consensus bases, heterozygous call left in the length partition | S3/S7 POA consensus plus polishing gives an ACGT consensus with no Clair3 IUPAC selection. S4 splits equal-length or Δ1 alleles only on ≥ 2 linked sites. |
| Independent haplotype evidence not established | `independent_haplotype_evidence` is true for length or linked-site splits; homozygous cases get `no_informative_heterozygosity`. |
| Per-allele depth gate in resolvable, high-depth cases | S6 all-read assignment by edit-distance competition; adequate means ≥ 30 spanning or ≥ 40 assigned reads. |
| Missing allele / unresolved alias, allele selection unresolved, reported length ≠ contig length | S2 smear-aware length model with recorded `rejected_peaks`. The length comes from the consensus itself. |
| Variant localization ambiguous / no explicit support / event identity | S10 per-event read-level support (`read_support.status`), classified on an ACGT consensus. |
| Depth below the gate (expected cases) | Not addressed by design. These stay INCONCLUSIVE (`depth_status` low or insufficient). |

## v0.16.0 on standard / clean / stress (Task 12c)

The owner asked for a realistic headline benchmark rather than one built from
edge cases. The benchmark now has named **sets** (see `docs/benchmark.md`,
"Benchmark sets"):

- `standard` (headline): ONT amplicon at 500/1000/2000 reads, HiFi at
  200/500/1000, genomic at 30/60/100; calibrated error; PCR bias none or
  calibrated; amplicon smear 0.24 and chimera 0.023 (the calibrated MucOneUp
  profile rates, about the real median).
- `clean` (control): the same haplotypes at the top standard depth, with no PCR
  bias, smear or chimera.
- `stress`: the harsh mix. The dev v2 data are reused as stress (v2 ≈ stress:
  the same factor mix, without the smear 0.5 level).

The caller is ladder v0.16.0 throughout. Aggregates only.

### Setup

- Out-roots: `../MucOneSpan-bench-data/v3/{standard,clean}`, both fresh, with
  dev n = 30 per profile (90 cases per set). Stress is `v2`, re-reported with
  the set-aware report (legacy designs map to `sets.legacy` = stress; the
  numbers are identical to Task 12b).
- Generation ran at `11755cc`; settings SHA-256 `0ef978d8…`, generation hash
  `8fed16a7…`. Realism, run, evaluate and report ran at `35849d3`, whose
  report settings hash `7543ea45…` differs only by `profiles.simulator_threads`,
  added after generation. v3 generation therefore used MucOneUp's
  config-default 8 simulator threads per case.
- Both sets generated in parallel with `--jobs 6` each. The later steps ran one
  set at a time, as the controller required after the machine load reached
  about 160: `run --threads 3 --jobs 4`.
- `standard` generation took 45 min 05 s. One HiFi case (1000 reads, alleles
  62/122 units, calibrated PCR) needed about 27,000 templates because of the
  minor-share floor, and spent about 40 min in `ccs`. `clean` generation took
  33 min 42 s. Depth was not lowered. All 180 cases were `ok`. Realism took
  4 min 20 s (standard) and 6 min 23 s (clean), with 0 failures. `run` took
  26 min 44 s (standard) and 31 min 13 s (clean), with 90/90 completed in
  each set.

### Results per set (ladder v0.16.0, dev, 90 cases per set)

| Set | Per-allele exact | Case exact | Event recall | FP on normals | INCONCLUSIVE | expected / resolvable / depth unknown |
| --- | --- | --- | --- | --- | --- | --- |
| **standard** | 12/180 = 0.067 [0.033, 0.106] | 0/90 | 11/57 | 0/33 | 86/90 | 5 / 81 / 0 |
| clean | 23/180 = 0.128 [0.067, 0.194] | 10/90 | 15/57 | 0/33 | 81/90 | 0 / 81 / 0 |
| stress (v2) | 22/180 = 0.122 [0.067, 0.178] | 5/90 | 12/57 | 0/33 | 81/90 | 43 / 38 / 0 |

Intervals are 95% cluster-bootstrap intervals from `report.json`.

| Set / profile | Per-allele exact | Event recall | PATHOGENIC on pathogenic truths | INCONCLUSIVE | NEGATIVE on normals |
| --- | --- | --- | --- | --- | --- |
| standard ONT amplicon | 4/60 | 3/19 | 1/19 | 29/30 | 0/11 |
| standard HiFi amplicon | 7/60 | 5/19 | 2/19 | 27/30 | 1/11 |
| standard ONT genomic | 1/60 | 3/19 | 0/19 | 30/30 | 0/11 |
| clean ONT amplicon | 14/60 | 6/19 | 3/19 | 26/30 | 1/11 |
| clean HiFi amplicon | 9/60 | 6/19 | 1/19 | 26/30 | 3/11 |
| clean ONT genomic | 0/60 | 3/19 | 1/19 | 29/30 | 0/11 |
| stress ONT amplicon | 8/60 | 2/19 | 0/19 | 28/30 | 1/11 |
| stress HiFi amplicon | 7/60 | 7/19 | 0/19 | 28/30 | 0/11 |
| stress ONT genomic | 7/60 | 3/19 | 2/19 | 25/30 | 0/11 |

There were no critical false negatives in standard or clean (stress: 6, all
`insufficient_evidence` at design depth 3-5), and no false PATHOGENIC calls in
any set.

**Standard vs clean, paired.** Both sets simulate the same haplotypes, so their
truth alleles pair. Per-allele exact, as discordant pairs (standard-only /
clean-only):
- ONT amplicon 1 / 11 (exact McNemar p = 0.006);
- HiFi 2 / 4;
- genomic 1 / 0;
- all profiles 4 / 15 (p = 0.019).

The typical-level artefacts (smear, chimera, PCR bias) and the lower standard
depths together cost the ladder most of its ONT amplicon exact alleles. Depth
and artefacts are confounded in this contrast. Two-thirds of standard cases
are below the clean depth, and genomic, which has no artefacts, barely
changes.

### Top INCONCLUSIVE causes (reason atlas, cases)

| Reason key | standard (86) | clean (81) | stress (81) |
| --- | --- | --- | --- |
| `evaluator: ambiguous_reconstruction` | 84 | 76 | 75 |
| `evaluator: iupac_bases` | 82 | 74 | 71 |
| `gate: heterozygous call left within the length-partitioned allele; consensus uses unresolved (iupac) selection` | 77 | 57 | 46 |
| `gate: reported length # differs from the consensus contig length #` | 76 | 34 | 31 |
| `gate: high number of ambiguous consensus bases (#) detected` | 69 | 62 | 60 |
| `gate: allele selection unresolved (unresolved_secondary_mode; secondary mode fraction #)` | 62 | 33 | 24 |
| `gate: observed sequence variant (<variant>) is inconclusive: localization ambiguous` | 57 | 42 | 49 |
| `gate: # primary alignments, below the per-allele depth gate (#)` | 29 | 27 | 50 |

At realistic depth, almost every INCONCLUSIVE is **resolvable**: 81/86 in
standard and 81/81 in clean, against 38/81 in stress. Sample depth no longer
explains them. The per-allele depth-gate reason still appears (29 standard and
27 clean cases), but it comes from called alleles with few primary records. One
example is a 7-read allele at 96 units in a 73/73-unit case with about 1000
reads per allele. It is not a lack of sample reads. What explains the
INCONCLUSIVE calls is the reconstruction:
IUPAC consensus, secondary length modes and length-selection disagreement. The
secondary-mode and length-disagreement gates roughly double from clean to
standard. Smear and chimera products at the real median rate form extra
clusters that the ladder cannot set aside:
- in the standard ONT amplicon set, 24/30 cases have a selection-unresolved
  reason and 26/30 a length disagreement;
- in the ONT genomic profile, 27/30 standard cases are selection-unresolved.
  That profile has no molecule artefacts, so the cause there is partial-read
  and depth behaviour, not PCR artefacts.

Per-allele exact does not rise with depth in any set: stress ONT amplicon at
>= 150 reads is 1/20 exact alleles, and standard ONT amplicon at 1000 reads is
0/20.

### Standard ONT amplicon vs real PRJEB92208 (v0.16.0)

| | Real PRJEB92208 amplicons | standard ONT amplicon | clean ONT amplicon |
| --- | --- | --- | --- |
| PATHOGENIC on pathogenic truths | 3/5 = 0.60 [0.15, 0.95] | 1/19 = 0.05 [0.00, 0.26] | 3/19 = 0.16 [0.03, 0.40] |
| INCONCLUSIVE | 6/11 = 0.55 [0.23, 0.83] | 29/30 = 0.97 [0.83, 1.00] | 26/30 = 0.87 [0.69, 0.96] |
| FP PATHOGENIC on normals | 0 | 0/11 | 0/11 |

(Clopper-Pearson 95% intervals.)

**The simulation is not in the same range as the real data.** It is now
harsher than real. The INCONCLUSIVE intervals of the real data and the standard
set do not overlap, and the PATHOGENIC rate is lower (the intervals overlap
only at the edges). The v2/stress data showed the same gap: 0/19 PATHOGENIC and
28/30 INCONCLUSIVE.

Known differences, which remain unresolved:
1. The event mix: 13 dictionary events with about 1-2 dupC per profile,
   against mostly the classic dupC in real MPs.
2. Depth: 500-2000 reads, against 10k-48k in real libraries.
3. The simulator's span-shape and strand-error gaps (below). Simulated spans
   are too narrow and shifted short, and the >= +53 bp tail from
   concatemers/chimeras is 3× the real share, a likely source of the spurious
   longer clusters seen in selection-unresolved cases.

Nothing was tuned to close the gap. PRJEB92208 is a check, not the target.

### Realism (indicative)

| Check | standard | clean | stress (v2) |
| --- | --- | --- | --- |
| ONT amplicon checks passed | 14/21 | 12/21 | 10/21 |
| ONT amplicon total error vs target 0.02138 | 0.02231 (+4%) | 0.02232 (+4%) | 0.02653 (+24%) |
| ONT amplicon deletion - strand vs 0.00485 | 0.00743 (+53%) | 0.00743 | 0.00783 |
| ONT amplicon mismatch + strand vs 0.00811 | 0.00570 (-30%) | 0.00570 | 0.00731 |
| Span-offset JS < 55 u / >= 55 u (limit 0.1) | 0.232 / 0.270 | 0.462 / 0.450 | 0.209 / 0.247 |
| Off by > 1 unit, median (real 0.2695) | 0.2686 | 0.0244 | 0.2054 |
| ONT genomic checks passed (2 HG002 WGS libraries; not targeted enrichment) | 11/15 | 10/15 | 10/15 |

**Correction to the premise.** The "calibrated error about 25% above target"
figure came from v2, where half of the cases use `error: poor` (×1.5). With
calibrated error only, the total ONT amplicon error is +4% of target and passes
the check. The remaining error gap is strand-specific: minus-strand deletions
are +53% and plus-strand mismatches -30%. The span-offset JS distance
(0.23-0.27) is the real outstanding realism gap. The standard smear setting
reproduces the real median off-peak share (0.2686 vs 0.2695). The upstream
draft is `.planning/2026-09-25-muconeup-realism-issue.md`.

### Fix round 1: artefact-free clean ONT amplicon (`v3/clean2`)

The v3 `clean` set was **not artefact-free**. The profile variant overrode only
smear and chimera, so the ONT amplicon cases kept the base profile's
concatemer rate (0.024) and off-target fraction (0.3). Their off-by->1-unit
median of 0.0244 was the concatemer rate.

Sets now fix `concatemer_levels` and `offtarget_levels` as well. `clean` sets
both to 0; `standard` and `stress` keep the profile values. I regenerated only
the 30 clean ONT amplicon cases into a fresh out-root, `v3/clean2`, at commit
`dfd36b9`. The design IDs, haplotypes and read seeds are the same as in v3
clean.

| Step | Wall | Result |
| --- | --- | --- |
| generate `--jobs 6` (2 simulator threads each) | 2 min 06 s | 30/30 ok |
| realism | 3 min 42 s | 0 failures; off-by->1-unit median 0.0, off-target 0.0 |
| run `--threads 3 --jobs 4` | 26 min 33 s | 30/30 completed |

| ONT amplicon | Per-allele exact | Case exact | Event recall | PATHOGENIC on pathogenic truths | INCONCLUSIVE | NEGATIVE on normals | FP |
| --- | --- | --- | --- | --- | --- | --- | --- |
| standard | 4/60 | 0/30 | 3/19 | 1/19 | 29/30 | 0/11 | 0/11 |
| v3 clean (concatemer and off-target left) | 14/60 | 7/30 | 6/19 | 3/19 | 26/30 | 1/11 | 0/11 |
| **clean2 (artefact-free)** | **21/60 = 0.35 [0.20, 0.50]** | 9/30 | 6/19 | **4/19 = 0.21 [0.06, 0.46]** | **21/30 = 0.70 [0.51, 0.85]** | 5/11 | 0/11 |

- The per-allele interval is the cluster bootstrap from `report.json`. The
  decision intervals are Clopper-Pearson.
- The clean2 atlas has 21/21 resolvable cases.
- Paired per-allele exact, as discordant pairs:
  - v3 clean-only / clean2-only: 0 / 7 (McNemar p = 0.016);
  - standard-only / clean2-only: 0 / 17 (p < 0.001).
- Selection-unresolved reasons fall from 24/30 (standard) and 15/30 (v3 clean)
  to 3/30. Length disagreements fall from 26 and 3 to 1.
- The remaining clean2 INCONCLUSIVE calls are IUPAC reconstructions (20/21).
  18 of them lose an allele (`missing_allele` / `unresolved_allele_alias`).

**Corrected conclusions.**
- Concatemers and off-target reads alone, at the base profile rates, cost the
  ladder a third of its clean exact alleles. Together with typical smear and
  chimera they explain most of the standard-vs-clean drop in the ONT amplicon
  set. This supports the concatemer-tail explanation for the spurious longer
  clusters.
- Against real PRJEB92208, the artefact-free clean2 ONT amplicon set overlaps
  the real intervals: PATHOGENIC 4/19 [0.06, 0.46] vs 3/5 [0.15, 0.95], and
  INCONCLUSIVE 21/30 [0.51, 0.85] vs 6/11 [0.23, 0.83].
- The **standard** headline set is still out of range (29/30 INCONCLUSIVE,
  [0.83, 1.00]). Real data contain the same artefact classes, so the likely
  reason is that simulated artefacts are more disruptive to the ladder than
  real ones. The span-shape and concatemer-tail gaps are in the upstream draft.
  The event mix and 10-50× lower depth also remain.
- The earlier "v3 clean" rows above describe a control that still had
  concatemers. Use clean2 for the ONT amplicon clean control. The HiFi and
  genomic clean rows are unaffected: HiFi has no base concatemer or off-target
  rate, and genomic gets no molecule step.
- The reason atlas now treats `stress`-set cases as expected
  (`atlas.expected_inconclusive_sets`). Re-reported, stress is 81 expected / 0
  resolvable; the depth condition alone would have given 43. Standard (5 / 81)
  and clean (0 / 81) are unchanged.
