# Proposed execution order: P0 fixes, benchmark, hybrid engine (2026-09-24)

Inputs:
- `2026-09-23-deep-review-roadmap.md` (P0–P4)
- `2026-09-23-hybrid-engine-plan.md` (12 tasks)
- `2026-09-24-realistic-benchmark-plan.md` (12 tasks)
- MucOneUp v0.45.0 (released 2026-09-24)

Status: proposal; needs owner approval. Merging to MucOneSpan `main` needs explicit approval each time.

## Principles

1. **Clinical safety first, on the engine people run today.** The ladder engine stays
   the default until the sealed benchmark says otherwise, so its known false negative
   must not wait for the new engine.
2. **Build the measuring instrument before tuning the new engine.** The prototype's
   failures (44-bp Δ window, 5-read peak, polishing at 43% support) were found on
   unrealistic simulations. Its thresholds should be chosen on the realistic `dev`/`val`
   splits: smear 8–52%, strand-specific stutter, chimeras, off-target reads.
3. **One sealed test, evaluated once**, after the decision rule is pre-registered.

## Phase 0: P0 clinical safety on v0.15.x (days; release v0.16.0)

Branch `fix/p0-clinical-safety`, issues first, TDD, one commit per issue:

| Order | Item | Where | Test |
| --- | --- | --- | --- |
| 0.1 | **D1 MP4 false negative.** `single_heterozygous_unordered` is not in the unphased set, so the consensus uses GT allele 1 (REF for `0/1`) and `independent_haplotype_evidence=True`. Decide per event from the allele-specific AD fraction: the reads were already split by length, so ALT ≥ 0.5 means the allele carries it. A fraction in the ambiguous band (0.2–0.5) gives an unphased consensus plus INCONCLUSIVE, never REF + NEGATIVE. Fix the test that asserts the bug. | `calling.py:528-539`; `tests/unit/test_distinct_calling.py:157-160` | MP4-shaped synthetic VCF (`G>GC`, AD 192,275) → dupC called on allele 2 |
| 0.2 | **D5** haploid rule: use the AD fraction, not FORMAT/AF. | `vcf.py:113-130` | AF 0.448 vs AD 0.59 fixture |
| 0.3 | **D4** evaluate VCF support per event, not globally; het indels applied by `-H I` are not "unresolvable". | `variant_support.py:42-51` | two-event fixture |
| 0.4 | **Hybrid-plan Task 1** (explicit-support clinical decision): a mutation without support is not "supported". This is the same safety theme and a prerequisite for both engines, so it moves here from the hybrid plan. | `report.py:85-100` | as in hybrid plan Task 1 |
| 0.5 | **D8 (#55)** NEGATIVE only when allele selection is resolved (no super-cluster) and spanning depth is adequate; PATHOGENIC only on event identity + frameshift + support. | `report.py:92-221` | super-cluster fixture → INCONCLUSIVE |
| 0.6 | D3, D7, D9, D10 small correctness fixes. | per roadmap | per issue |

Validation: `make ci-check`, `make test-int`, then a rerun of the PRJEB92208 cohort (local outputs).
The target is MP4 dupC detected, with no new positives in HG001–HG004/MP5. Release after
owner approval of the merge.

## Phase 1: Benchmark instrument (benchmark plan Tasks 1–9, then 10–11 with the ladder engine)

- Tasks 1–9 are independent of any engine.
- Run Task 10 with `--engines ladder` only; it adds the `engine=` forwarding
  (a no-op for `ladder`) that hybrid Task 10 will reuse.
- Deliverables:
  - a `dev` split (900 cases) with a **realism report** against the public PRJEB92208 targets;
  - a first **ladder baseline** on `dev`, comparing v0.15.1 with the Phase 0 build. This is the first quantitative evidence for P0 on realistic data (MP4-like cases, the INCONCLUSIVE rate).
- Realism failures are fixed **upstream in MucOneUp** (profile recalibration → 0.45.x), not
  in MucOneSpan.

## Phase 2: Hybrid engine (hybrid plan Tasks 2–9), developed on `dev`, thresholds on `val`

- Tasks 2–6 (settings, spans, lengths, POA/polish, phase split) are pure library stages
  with synthetic tests. They touch files disjoint from Phase 1, so **they may run in
  parallel with Phase 1** once Phase 0 is merged.
- Tasks 7–9 (assignment, evidence/LLR, engine orchestration) run after Phase 1, so the
  review-focus cases can be checked on realistic `dev` cases as well as synthetic ones.
  The review-focus cases are 40–50% smear, a 44-bp Δ, a 5-read long allele, and one strand only.
- Threshold selection (peak separation, callable depth gate, homopolymer LLR cutoffs) uses
  `val`. Every look is logged in the ledger.
- Hybrid Task 10 (clinical gating on hybrid evidence, evaluation acceptance) merges with
  benchmark Task 10's `engine` forwarding.

## Phase 3: Sealed evaluation and decision

1. `benchsim preregister` records the decision rule (benchmark spec §6).
2. Generate `test` (2,400 cases) with the secret salt outside the repository.
3. Run `ladder` (Phase 0 build) and `hybrid`, then `evaluate` and `report`: McNemar + Holm, and
   non-inferiority on the normal false-positive rate. No increase in false NO_PATHOGENIC is allowed.
4. Hybrid Task 11 real-data regressions: the frozen simpanel/held-out sets, PRJEB92208
   (public), and the in-house samples (local only; LB-level outputs never in Git).
5. Hybrid Task 12 docs. Flip the default to `--engine hybrid` **only if** the rule passes
   for every profile. HiFi conclusions are reported separately (uncalibrated profile).

## Dependency sketch

```
Phase 0 (P0 fixes + hybrid T1) ──┬──> Phase 1 (bench T1–T9) ──> bench T10/T11 (ladder) ──┐
                                 └──> hybrid T2–T6 ───────────> hybrid T7–T9 ─> T10 ─────┴─> Phase 3
```

## Not in scope yet

- Read-assignment accuracy scoring waits until the hybrid engine writes per-read assignments.
- MucOneUp follow-ups #125–#129 are not blockers for the benchmark:
  - #125 affects `pacbio_params.coverage` only when `read_simulation.coverage` is unset; benchsim always passes `--coverage`.
  - #127 affects PacBio WGS manifests, which benchsim does not use.
