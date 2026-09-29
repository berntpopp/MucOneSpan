# Session handoff (2026-09-24)

## A. MucOneSpan deep review — DONE (planning only, nothing committed)
- Roadmap: `.planning/2026-09-23-deep-review-roadmap.md` (defects D1–D10, priorities P0–P4).
- Specs/plans: `.planning/2026-09-23-hybrid-engine-{spec,plan}.md`,
  `.planning/2026-09-23-realistic-benchmark-spec.md` (genomic profile = decided: extend MucOneUp).
- Evidence (outside Git, patient data local only): `../MucOneSpan-review-20260923/`
  (poa-prototype, homopolymer, simpanel, heldout, realprofile/targets.json, inhouse/).
- Key numbers: POA prototype 80/80 dev, 76/80 held-out vs 45–49/80 ladder; dupC 4/4 PRJEB92208;
  MP4 false-NEGATIVE bug (`calling.py:528-539`); one in-house sample (ID kept local) dupC motif 33 (82-unit allele, 7 reads).
- Open: benchmark PLAN not yet written (was blocked on MucOneUp work; now unblocked).

## B. MucOneUp realistic read simulation — IN PROGRESS (almost releasable)
- Repo `../MucOneUp`, branch `feat/realistic-read-simulation` (pushed to origin), 24 commits, version bumped to 0.45.0.
- Issues filed/addressed: #97, #100–#111 (all referenced with Closes in PR body).
- `uv run make ci-check`: 1549 passed, 11 skipped; `mkdocs build --strict` passes.
- PR body drafted: `../MucOneSpan-review-20260923/muconeup-pr-body.md`
- An independent code review of `main...HEAD` was running when the session ended; its result is lost — rerun it.

### Remaining steps (user authorized: document in PR/commits, bump, tag, release)
1. Rerun independent review of `git diff main...HEAD`; fix real findings (TDD, one commit each).
2. `uv run make ci-check` + `uv run --extra docs mkdocs build --strict`.
3. Open PR (`gh pr create` with the body), wait for GitHub CI green.
4. Merge PR, tag `v0.45.0` on main, `gh release create v0.45.0` with changelog notes
   (docs/about/changelog.md [0.45.0]); verify Zenodo/PyPI workflows if they trigger.
5. Then back to MucOneSpan: write the benchmark plan using MucOneUp 0.45.0
   (`--read-profile ont_r10_sup_amplicon_v1`, `reads ont --simulator pbsim3-fragments`, truth TSV).

### Key design decisions (MucOneUp)
- Everything opt-in; legacy output byte-identical (verified md5 with real pbsim3).
- pbsim3 ONT models floor at ~3.5–4% error in templ mode -> added empirical error channel
  (Sequencer strategy: PbsimRun | EmpiricalSequencer); ONT built-in profiles use it.
- PCR preset madritsch2025_r10 alpha = 9.27e-5 (issue #104 comment corrects 7.7e-5).
- Profiles hold aggregate public (PRJEB92208) stats only; in-house aggregates NOT shipped (needs owner approval).
- Known gaps: e2e gate wants `pbsim3` binary (conda installs `pbsim`); HiFi profile uncalibrated;
  fragment reads clip at 10 kb flanks (use --flank-fasta); #88 (VNTR downsample architecture) out of scope.

## C. Update (2026-09-24, later session): MucOneUp v0.45.0 RELEASED
- Three independent reviews (2 Claude + Codex) of `main...HEAD`; 13 verified findings filed
  (#112–#124) and fixed test-first, one commit each; follow-ups filed #125–#129 (not blockers).
  Notable: pbsim3 rejects `--accuracy-sd` (#115); off-target reads crossed artificial junctions (#112);
  NanoSim `--no-align` still aligned (#116); config beat profile for fragment lengths/PCR (#114).
- Legacy byte-identity re-verified with real pbsim3/ccs vs main: ONT md5 301d7171…, HiFi 1f500a8b….
- ci-check 1590 passed / 11 skipped; mkdocs --strict ok; PR #130 CI green → merged (eabb721),
  tag v0.45.0, GitHub release published. Docker (release+main), docs, tests all green; Zenodo v0.45.0 DOI 10.5281/zenodo.22927441.
- MucOneSpan planning written (not committed, awaiting review):
  `.planning/2026-09-24-realistic-benchmark-plan.md`, `.planning/2026-09-24-execution-order.md`.
- Next: owner reviews both; then Phase 0 (P0 clinical fixes incl. D1 MP4 at calling.py:528-539).

## D. Execution state (2026-09-24, "finish end to end")
Worktrees (each has a SDD ledger in .superpowers/sdd/<plan>/progress.md — trust ledgers + git log after compaction):
- ../MucOneSpan-p0  branch fix/p0-clinical-safety (v0.16.0): ALL tasks done + final review + fix wave; PR #71 open, CI green, MERGEABLE.
  WAITING ON OWNER: (1) accept MP5 PATHOGENIC→INCONCLUSIVE; (2) approve merge + tag v0.16.0 + release. Issues filed: #63-#68, #69 (benign dictionary), #70 (AD min depth).
- ../MucOneSpan-bench branch feat/benchsim (from 8a96697 planning commit): Tasks 1-10 complete; Task 11 implemented (5b5390a), in review; ruling: metric 1 per-allele pairs (fix round pending). Task 12 (docs + dev pilot) remaining; then final whole-branch review.
  Known realism gaps for Task 12/upstream MucOneUp: amplicon span-offset JS 0.24-0.45 (>0.1), genomic error rates above WGS target.
- ../MucOneSpan-hybrid branch feat/hybrid-engine (stacked on fix/p0-clinical-safety 2675e9b): pre-flight of hybrid plan vs post-P0 code running → .planning/2026-09-24-hybrid-plan-preflight.md; then SDD over hybrid Tasks 2-12 (Task 1 done in P0).
- MucOneUp v0.45.0 released; Zenodo DOI 10.5281/zenodo.22927441.

## E. Resume point (after session loss, 2026-09-25 morning)
- Machine did NOT crash (uptime 22 h, no OOM); the Claude Code process exited and killed background agents.
- P0: DONE 2026-09-24 — owner accepted MP5 → INCONCLUSIVE; PR #71 merged (f0b98ac), tag v0.16.0, GitHub release published. Hybrid branch base 2675e9b is now in main.
- Bench (../MucOneSpan-bench, feat/benchsim @495dc25): Task 12 in progress. Uncommitted docs edits (CHANGELOG.md, docs/development.md, mkdocs.yml, docs/benchmark.md) in worktree. Dev pilot data in ../MucOneSpan-bench-data/dev: 89/90 cases (82 ok, 5 design_invalid, 2 generation_failed), no manifest yet → rerun `benchsim generate` (idempotent), then realism, run --engines ladder, evaluate, report, pilot report .planning/2026-09-24-benchsim-dev-pilot.md. Then final whole-branch review + fix wave (include ledger minors and the cross-plan edlib marker ruling), push, PR.
- Hybrid (../MucOneSpan-hybrid, feat/hybrid-engine @59ef562): executing .planning/2026-09-24-hybrid-engine-plan-v2.md via SDD; Task 2 partially done (uncommitted Makefile/pyproject/settings/uv.lock/tests/unit/test_hybrid_settings.py) — resume Task 2, then Tasks 3–14.
- Keep parallel load moderate (≤ ~12 jobs total) during pilot runs.

## F. Resume point (2026-09-25, midday) — READ THIS FIRST
Trust the SDD ledgers + `git log` over any memory. Ledgers:
- Bench: ../MucOneSpan-bench/.superpowers/sdd/2026-09-24-realistic-benchmark-plan/progress.md
- Hybrid: ../MucOneSpan-hybrid/.superpowers/sdd/2026-09-24-hybrid-engine-plan-v2/progress.md
Owner directives (binding, in both global-constraints.md): fully config-driven, NO magic numbers; calibration commands (hybrid Task 15, .planning/2026-09-25-calibration-addendum.md in the hybrid worktree, untracked — commit with Task 15). INCONCLUSIVE 73/90 in pilot must improve: cause is ladder IUPAC reconstruction (71/73), not gates — do not relax gates.
- P0/v0.16.0: DONE (merged f0b98ac, tagged, released).
- Bench (feat/benchsim): Tasks 12, 12b, 12c complete (sets standard/clean/stress; clean2 = artefact-free ONT amplicon overlaps real PRJEB92208). MucOneUp realism issue filed: berntpopp/MucOneUp#131. Task 12d running: 2×2 {old caller, v0.16.0} × {legacy MucOneUp reads, new clean reads} → .planning/2026-09-25-ladder-2x2-diagnosis.md (old-caller worktrees ../MucOneSpan-oldcaller-<tag>). Then final whole-branch review + ONE fix wave → ask owner before push/PR. OPEN QUESTION to owner: v0.16.1 ladder fixes (read loss at allele split; hardcoded ambiguous_bases>10, LEGACY_MIN_TOTAL_READS; concatemer sensitivity; 77→142 length error).
- Hybrid (feat/hybrid-engine): Tasks 2-11 complete (end-to-end `run --engine hybrid` works; smoke case exact). Task 12 (evaluation acceptance) running; then 13 (integration + real PRJEB92208/in-house regression), 14 (docs), 15 (calibration). Preliminary hybrid benchmark running → ../MucOneSpan-bench-data/hybrid-preliminary.md. MucOneUp issues filed: #131 (realism), #132 (C8/A/T homopolymers error-free → dupC sensitivity overstated).
- Keep load ≤ ~12 jobs.

## G. Owner decisions (2026-09-25 afternoon)
- MucOneUp #132 (homopolymer stutter coverage): authorized PR + merge + bump + tag + release → released as **0.46.0** (breaking profile-format change ⇒ minor bump). Worktree ../MucOneUp-hp. #131 follows as separate PR.
- v0.16.1 ladder fixes AND hybrid focus: plan being written → .planning/2026-09-25-v0.16.1-ladder-plan.md; issues first, TDD, PR, release.
- Decision-rule absolute targets approved: clean ≥0.90 PATHOGENIC / ≤0.10 INCONCLUSIVE / 0 FP; standard ≥0.80 / ≤0.20 / 0 FP → bench Task 12e.
- 2×2 diagnosis done (bench .planning/2026-09-25-ladder-2x2-diagnosis.md): no regression; hard designs + stricter metrics. Safety: R9-era ONT reads → 12/39 pathogenic NEGATIVE in v0.16.0.
- Never mention the external caller repo/paper that was evaluated (memory rule); ideas P1–P5 in hybrid calibration addendum.
- Old-caller worktrees ../MucOneSpan-oldcaller-{v0.8.0,v0.14.0,v0.15.1} can be removed when no longer needed.

## H. RESUME HERE (2026-09-25, snapshot for a new session)
Read order: this section → the three SDD ledgers → `git log` in each worktree. Trust ledgers + git over memory.
Ledgers: ../MucOneSpan-bench/.superpowers/sdd/2026-09-24-realistic-benchmark-plan/progress.md ·
../MucOneSpan-hybrid/.superpowers/sdd/2026-09-24-hybrid-engine-plan-v2/progress.md
Binding owner rules (also in memory): fully config-driven, no magic numbers; never mention the evaluated external caller repo/paper; issues first → TDD → documented PR → bump/tag/release; ask before push/merge/release unless authorized (authorized so far: MucOneUp 0.46.0 release only).

### 1. MucOneUp 0.46.0 (#132) — RELEASED 2026-09-25
- PR #133 merged (86334cf), tag v0.46.0, GitHub release published; local ../MucOneUp on 86334cf, `muconeup --version` = 0.46.0 (editable install).
- Post-release workflows all green (Docker tag+main, Test & Quality, Docs, Dependency Graph); no PyPI workflow exists. Zenodo DOI appears automatically.
- Next MucOneUp: #131 (per-strand error rates, span-length shape) as separate PR; small fix: Makefile `ci-check` uses bare `pytest` (needs `uv run`). Worktree ../MucOneUp-hp can be removed (`git worktree remove`).

### 2. Hybrid engine (feat/hybrid-engine, ../MucOneSpan-hybrid @ c23f36a)
- Tasks 2–12 complete and reviewed. Fix F1 committed (6a40be9 infix_hit None-start crash; 3080dc9 confidence = read-support evidence; c23f36a polish-corrects test) — report `.superpowers/sdd/.../fix-f1-report.md`; **needs review** (scoped: 6d8e012..c23f36a), then rerun ont_genomic_targeted cases.
- Next: Task 13 (integration test + real PRJEB92208 + in-house regression; carry Task 11 ⚠️: evidence gets only assigned+oriented spanning reads; check R9-era ONT dupC false-negative risk), Task 14 docs (include unresolved_group_size, het_af_min floor, not_assessed rule, detection limits), Task 15 calibration (.planning/2026-09-25-calibration-addendum.md, untracked, commit with Task 15; P1–P5 follow-ups). Final whole-branch review: deferred minors in ledger (Anchors "1"/"9" hardcode, test literals, etc.).
- Preliminary benchmark (pre-F1, ef59e6f): ../MucOneSpan-bench-data/hybrid-preliminary.md — standard HiFi allele exact 88% vs ladder 12%, ONT 90% vs 7%, PATHOGENIC 18/19 & 16/19, FP 0; dupC numbers inflated until MucOneUp 0.46.0 data is regenerated.

### 3. Benchmark (feat/benchsim, ../MucOneSpan-bench @ 8dcb63e)
- Tasks 1–12d complete. Task 12e implemented (9fb79d5, absolute targets; INCONCLUSIVE rate over all cases) — **needs review** (8dcb63e..9fb79d5).
- Then: regenerate benchmark data with MucOneUp 0.46.0 (C8 stutter realistic), rerun ladder + hybrid; final whole-branch review + ONE fix wave; ask owner before push/PR.

### 4. v0.16.1 ladder (main) — approved; PLAN READY
- Plan: .planning/2026-09-25-v0.16.1-ladder-plan.md (Tasks 0–8). Scope: Fix A (safety: block NEGATIVE when Clair3 pileup saw a frameshift the applied calls dropped — ms_ont_sub 12/39 → 1/39 NEGATIVE on pathogenic, 0 changes elsewhere by replay), Fix B (clinical_decision settings for ambiguous_bases>10 and LEGACY_MIN_TOTAL_READS), Fix C (ladder allele-selection literals → settings; fixes hardcoded 9 fixed repeats). Deferred to hybrid: mixed-partition split, span-length assignment, positive read support, genomic length doubling (docs only).
- Next: Task 0 (worktree ../MucOneSpan-v0161 from f0b98ac, baseline snapshot, **file issues A/B/C** — filing is part of the approved flow but the plan marks it "owner approval required": confirm with owner), then SDD Tasks 1–8; PRJEB92208 no-regression gate (Task 7); PR; release only after owner OK.
- Owner questions in plan: Fix C in the patch or v0.16.2; accept residual R9 FN (pair_5178); min_depth 10 vs 30.

### Housekeeping
- Old-caller worktrees ../MucOneSpan-oldcaller-{v0.8.0,v0.14.0,v0.15.1}: remove with `git worktree remove` when no longer needed.
- External evaluation files live only in ~/development/external/ (never copy into repos).
- Keep machine load ≤ ~12 threads (other projects run on this host).

### Owner decisions (2026-09-25, v0.16.1 plan questions)
- Fix C ships in v0.16.1 (Tasks 4–5 kept).
- Residual R9 FN (pair_5178) accepted; document in limitations.md (NEGATIVE from R9-model ONT runs not validated); hybrid S10 covers it later.
- `calling.stage_discordance_min_depth` default = 10.
- Filing issues A/B/C and commenting on #58/#54: AUTHORIZED. PR push, merge, tag, release: still ask.
- (2026-09-25 evening) OWNER DECISION: hybrid becomes the DEFAULT engine (supersedes "ladder stays default until the benchmark decision rule passes"). Ladder stays available via `--engine ladder`. Implemented as hybrid Task 16 after Task 15 calibration; ships as v0.17.0 (v0.16.1 stays a ladder-only patch).
- (2026-09-25) OWNER AUTHORIZED v0.16.1: push fix/v0.16.1-ladder, open PR (Fixes #72/#73/#74), merge after CI green, tag v0.16.1, GitHub release. Stop if CI fails.

### v0.16.1 RELEASED (2026-09-25): PR #75, merge f19d63a, tag v0.16.1. Rulings made during execution (from SDD ledger):
- Ruling: accept R1–R17 in preflight-rulings.md (budgets raised within the 649 cap; tests strengthened; allele_2 criterion scoped to own partition) — the plan's budgets/tests were estimates that can't be met or prove nothing — cost if wrong: slightly larger files than planned.
- Ruling: calling.py budget raised from ≤600 to ≤610 (hard cap 649) so Task 5's ~2 lines fit without an unrelated extraction — cost if wrong: none (well under cap).
- Ruling (amends R1, which was self-contradictory: alleles.py had no __all__): drop __all__; re-export moved names with the redundant-alias form (`from .ladder_clusters import X as X`), which ruff treats as an explicit re-export — star-import surface identical to 64c43ce — cost if wrong: none.
- Ruling: final fix wave = I1 + M1 + M2 + M3 (cheap, tested); M4 left (history rewrite not worth it; PR body links #73); T5 minors (min 0, forwarding test, early tolerance) → follow-up issue after release — cost if wrong: none.
- Preflight rulings R1–R17 (budgets, tests, allele_2 scope, caveat for missing pileup): see ledger copy below.
  # Pre-flight scan rulings (v0.16.1) — binding on implementers; carried per task
  Scan by opus Plan agent vs main @ f0b98ac (full table summarized in progress.md).
  
  - R1 (T4/T5 alleles.py budget): ~109 lines move, not ~130. Budgets: after T4 ≤550, after T5 ≤575. `alleles.py` re-exports moved names via a complete `__all__` (no F401, star-imports unchanged).
  - R2 (settings.py budget): new fields may use a table-driven (name, minimum) validation loop — no line compression; budget ≤560.
  - R3 (test_alleles.py 520→~620): Task 5 tests go into new tests/unit/test_allele_selection_settings.py (min_haplotype_reads test into test_haplotag_calling.py).
  - R4 (T2 vs T6/T7): stage_concordance is required on each allele that has its own Clair3 partition (unphased same-length alias shares allele_1's). Drop the dead alias pop-tuple edit (calling.py:431-441).
  - R5 (T5 tests prove nothing): use three valleys (e.g. 6, 12, 20 with 6 lowest) so default vs non-default pick differ; use a unit length other than 60.
  - R6 (T5 vs Issue C #74): one parametrized test per new split setting proving a non-default value is honoured (read_length_split_min_reads, bin_bp, min_delta_bp, max_delta_bp, offset_bp, unit_tolerance_bp).
  - R7 (T1 cleanup test): mocked bcftools norm must create the -o file so "normalized file is gone" can fail; add a test for the final_vcf_unavailable branch.
  - R8 (T1 work_dir): default work_dir = pileup_vcf.parent.
  - R9 (replay.py): save every row (not only counts/changed rows) so Task 6 can check run-by-run before == baseline. Task 0 baseline must be re-saved with rows before Task 6 (re-run --baseline-only with main venv).
  - R10 (T6 claim): replay only re-scores summaries; it does not prove Task 5 neutral — reword; Task 5 Step 3 (allele replay) is the proof.
  - R11 (v2 stress): n = 84 (6 pre-existing runs without summary.json); state why.
  - R12 (T3 resolver): extract a shared section builder from load_settings (missing keys → defaults; unknown keys → ValueError); do not duplicate the unknown-key logic.
  - R13 (T5 wiring): detect_alleles passes reference_layout to _split_cluster_by_indel. Neighbour window `2 * bin` factor: structural constant, name it and say so in the docstring.
  - R14 (T8 CHANGELOG): fix link references — add [0.16.1] and [0.16.0], [Unreleased] compares from v0.16.1.
  - R15 (T3): clinical_gates.LEGACY_MIN_TOTAL_READS loses its consumer → remove it (no aliases for dead constants) unless a test/import needs it; then keep as documented alias of the settings default.
  - R16 (T2 safety): missing pileup.vcf.gz → status not_assessed must not pass silently: log a warning and add a quality-caveat reason (does not block NEGATIVE — same as other not_assessed semantics).
  - R17 (T1 normalisation wording): final VCF multi-ALT is split not trimmed → may over-flag (safe); correct open-risk 2 wording in docs/PR.
  - report.py budget ≤460 (was ≤450).
- Follow-up (needs owner OK to file): T5 minors — dominance_zero_primary_extra_reads min 0; forwarding test covers only last dominance call; early unit_tolerance validation; M4 commit 64c43ce says (#B).

## I. RESUME HERE (2026-09-25 night — session stopped on API spend limit; weekly limit resets 2026-09-28)
- v0.16.1 RELEASED (PR #75, f19d63a). Plan + ledger archived in .planning/archive/.
- OWNER DECISIONS this session: hybrid becomes DEFAULT (hybrid Task 16 after calibration → v0.17.0); bench FP criterion = absolute 0-FP (drop relative margin); INCONCLUSIVE-rate denominator = all cases (owner confirmed).
- Hybrid (../MucOneSpan-hybrid @ d8057ec): Tasks 2–14, 13b, 15a–d done+reviewed; benchsim merged locally (a1bd0b0). Task 15e (SAFETY: equal-length single-event het → NEGATIVE, case dev-clean-ont_amplicon_r10-0005) STOPPED mid-way with UNCOMMITTED changes — see hybrid ledger last line. Then Task 15-run (calibration on v4; run-length-aware stutter before relaxing event_max_alternative_frac), refresh prjeb92208/hybrid-engine.json, Task 16 default switch, final review.
- v4 hybrid baseline (bea5c13): standard P 84.2% / I 27.8% / FP 0; clean 86.0/23.3/0; clean2 78.9/10.0/0; dupC recall 100%; ladder ~7–26% P. Details ../MucOneSpan-bench-data/v4/README-v4.md.
- Bench (../MucOneSpan-bench @ d2e3a06): final review done, fix wave committed — needs scoped re-review, then C1 implementation, re-merge into hybrid, owner OK for push/PR.
- Follow-ups needing owner OK to file: v0.16.1 T5 minors (see archive ledger); MucOneUp #131; remove old worktrees (oldcaller ×3, MucOneUp-hp, MucOneSpan-v0161).
- OWNER ANSWERS (2026-09-25, Q&A): INCONCLUSIVE denominator = all cases (confirmed). Simpanel 77/80 does NOT block the default switch (document; P6 later). Hybrid + bench when finished: push + open PRs (bench first, then hybrid), then ASK before merge/v0.17.0 release. v0.16.1 follow-up issue filed → #76. MucOneUp #131 after hybrid calibration (issue → TDD → PR → ask before release). Old worktrees removed (oldcaller ×3, MucOneUp-hp, MucOneSpan-v0161).
- OWNER ANSWERS (2026-09-25, fill-in list):
  1. Hybrid is the default for ALL input types (amplicon and genomic).
  2. Ladder engine: DEPRECATE (deprecation warning on --engine ladder; docs mark legacy; removal in a later release).
  3. Version v0.17.0.
  4. Sealed test split: run BOTH before the v0.17.0 release (formal gate, pre-registered) and after the release (re-check with the released build, same pre-registered rule hash).
  5. Calibration rank: hard constraints 0 FP and 0 NEGATIVE on pathogenic; then lowest INCONCLUSIVE rate, then most alleles sequence-exact.
  6. Calibrated default changes may land without asking (dev gain, val confirmed, commit cites report; visible in PR).
  7. event_max_alternative_frac relaxed only after the run-length-aware stutter model.
  8. Equal-length single-event het: PATHOGENIC when supported; otherwise INCONCLUSIVE with a located reason.
  9. ../MucOneSpan-p0 removed.
  10. No cost reduction: keep current model selection and parallelism.
  11. No extra v0.17.0 scope; P1–P6 after.
- Bench PR #77 open (feat/benchsim @ 8c889a0, CI green) — awaiting owner merge. After merge, hybrid should merge main (not feat/benchsim) to pick up bench fixes.
- (2026-09-25) Owner: 'please merge' → PR #77 (feat/benchsim) merged into main.
- OWNER (2026-09-26): fix the INCONCLUSIVE cause (ONT amplicon dimer peaks at 2× length + inter-allele smear) BEFORE v0.17.0 → hybrid Task 15h.
- (2026-09-26) Hybrid PR #78 open (feat/hybrid-engine @ bb6c986), CI 13/13 green; awaiting owner merge/release. OWNER: sealed test split BEFORE release at reduced size (~150 cases/profile) — running.
- (2026-09-27) SEALED TEST test150: NOT ADOPTED (2 hybrid NEG-on-pathogenic safety gaps: 0105 flags not blocking NEGATIVE; 0138 real allele relabelled smear by 15h; HiFi INCONCLUSIVE per-profile over target). Ladder v0.16.1 had 4 FP dupC on equal-length ONT normals. README: ../MucOneSpan-bench-data/test150/README-test150.md. Next: hybrid Task 15i safety fixes → new sealed split. PR #78 must not be merged as-is.
- (2026-09-27) OWNER: ladder FP → issue #79 filed; fix via v0.17.0 hybrid default, no ladder patch.

## J. RESUME HERE (2026-09-27 — clean session switch after hybrid Task 15i)
Read this section, then the hybrid SDD ledger: ../MucOneSpan-hybrid/.superpowers/sdd/2026-09-24-hybrid-engine-plan-v2/progress.md (trust ledger + git over memory). Plan now archived at ../MucOneSpan-hybrid/.planning/archive/2026-09-24-hybrid-engine-plan-v2.md (plan-path file in the SDD dir points there).

State:
- Released: MucOneSpan v0.16.1 (PR #75). Merged: bench harness PR #77. Issues: #76 (v0.16.1 follow-ups), #79 (ladder false dupC on equal-length ONT normals; fix = hybrid default in v0.17.0, no ladder patch).
- Hybrid PR #78 (feat/hybrid-engine) OPEN, NOT merged. Local HEAD b24a677 = PR head bb6c986 + Task 15i (not pushed yet).
- Sealed test split test150 (150/profile, rule 9249bf48) was run at bb6c986 → NOT ADOPTED (README ../MucOneSpan-bench-data/test150/README-test150.md). That split is now UNSEALED for this rule (non-blind).
- Task 15i (b24a677, reviewed/approved) closed both hybrid NEGATIVE-on-pathogenic paths (0105 → INCONCLUSIVE located; 0138 → PATHOGENIC). Dev/val/panels/PRJEB unchanged, FP 0.
- Remaining blocker: HiFi INCONCLUSIVE per profile (test: clean .193 vs ≤.10; standard .213 vs ≤.20).

Next (in order, SDD with reviews; owner decisions already given):
1. Hybrid Task 15j — brief: ../MucOneSpan-hybrid/.superpowers/sdd/2026-09-24-hybrid-engine-plan-v2/task-15j-brief.md (read-quality-aware site detection; dev-only calibration, val confirm; fold in 15i minors). Opus implementer + opus review.
2. Push feat/hybrid-engine to update PR #78 (push authorized; merge NOT).
3. NEW sealed test split with a NEW secret salt (owner: reduced size 150/profile; same rule 9249bf48 unless changed): new root e.g. ../MucOneSpan-bench-data/test150b/, preregister first, generate (use --jobs up to 4–6 with simulator_threads 2 within ≤12 threads; do NOT run a second generate in parallel on the same set — it collided last time), run ladder + hybrid, evaluate, report. benchsim run has no resume and is slow for ladder (~6 h / 900 cases at 2×3 threads) — consider --jobs 4 --threads 3.
4. If ADOPTED: ask owner to merge PR #78 + tag/release v0.17.0 (owner must approve merge/release); then the post-release sealed re-check (owner wants both before and after).
5. After v0.17.0: MucOneUp #131; P1–P6 follow-ups; #76; deferred minors → one follow-up issue.

Binding owner rules: fully config-driven, no magic numbers; never mention the evaluated external caller repo/paper; issues first → TDD → documented PR → bump/tag/release; ask before merge/release; ≤ ~12 threads; no cost reduction on models; hybrid default for ALL inputs; ladder deprecated; v0.17.0; calibration objective 0 FP + 0 NEG-on-pathogenic, then lowest INCONCLUSIVE, then most sequence-exact; calibrated defaults may land without asking (dev gain + val confirm + cited report).

## K. RESUME HERE (2026-09-29 morning — handover to another agent)
- Hybrid branch feat/hybrid-engine (../MucOneSpan-hybrid): CODE FROZEN at 1d2c416 (engine). Pushed head 32f732e (CI all green). Local commits 4b98acb..1d2c416 not pushed. Uncommitted docs refresh (CHANGELOG, docs/benchmark.md, configuration.md, limitations.md, prjeb92208/hybrid-engine.json) from final numbers run calibfinal/ffinal5 (done 02:18).
- test150b (../MucOneSpan-bench-data/test150b): SEALED. Ladder run done (be5c3c5). Rule v7 preregistered (sha c77d513f…, config bench-config-test150b-v7.json). Hybrid run started 07:20 as systemd unit test150b-hybrid from ../MucOneSpan-gen @ 1d2c416 (run_hybrid.sh).
- Ledger: ../MucOneSpan-hybrid/.superpowers/sdd/2026-09-24-hybrid-engine-plan-v2/progress.md (all rulings/owner decisions).
- (2026-09-29) SEALED TEST test150b: ADOPTED under rule v7 (sha c77d513f...). Hybrid vs ladder: 0 FP (0/318), 0 NEG-on-pathogenic (0/582), pooled standard I 0.178 (<=0.20), clean I 0.122 (<=0.15), pathogenic standard 0.883 (>=0.80), clean 0.928 (>=0.90), allele-exact superior on all 3 profiles (p<1e-74). PR #78 pushed (3919878), CI 17/17 all green. Owner approved: PR #78 MERGED into main (aafb2a5). Tagging and publishing v0.17.0.
