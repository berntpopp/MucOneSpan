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
  tag v0.45.0, GitHub release published. Docker/docs workflows and Zenodo: see final session report.
- MucOneSpan planning written (not committed, awaiting review):
  `.planning/2026-09-24-realistic-benchmark-plan.md`, `.planning/2026-09-24-execution-order.md`.
- Next: owner reviews both; then Phase 0 (P0 clinical fixes incl. D1 MP4 at calling.py:528-539).
