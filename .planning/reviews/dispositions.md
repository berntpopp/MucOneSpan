# Review finding dispositions

Actual reviewer: installed Claude Code2.1.268, explicitly requested claude-fable-5-1;
modelUsage confirms requested model. Spec/plan original text retained separately.
Reviews were read completely; findings are checked rather than treated as facts.

| Finding | Disposition and evidence |
|---|---|
| Empty review JSON | Rejected as final failure: files were read during running processes; both final responses nonempty success with Fable modelUsage. Save process exit status too. |
| Stop rule loses downstream mutations | Confirmed: 106 cached consensuses, five named events lost; ignored classification/stop-rule-regressions.json. Strict mode rejected as default. Candidate default retains scanning plus ambiguity metadata; regression fixture added. |
| Unowned coordinate projection | Resolved interface: consensus_context actual reference/full/trim/sample/haplotype; shared with exact replay/reversion support. Unit flank-indel and wrong-trim tests pass. Real normalization tests pending. |
| Haplotype selectors error on unphased GT | Rejected: real bcftools1.17 integration GT0/1 insertion returns REF under-H1 and ALT under-H2, both exit0. No claim unphased multilocus selection establishes biological phase. |
| No read-backed phase producer | Confirmed limitation: PS/GT evidence path implemented and tested; ordinary Clair3 lacks phase. Read-backed phase remains unmet; honest mixed candidate is a safety change. |
| Duplicate observations inflate exact pair | Addressed in evaluator source metadata/independent genotype selection tests; duplicated single sequence cannot earn exact diploid success. |
| Missing status producer | Coordinator implementing schema1 run_status.json; loader tested against stale output and insufficient evidence. |
| Raw normalization mismatches | Addressed by whole-consensus replay and event reversion; no raw coordinate equality. Full real-tool mutation matrix pending. |
| Frameshift/score semantic migration | Signed-net and exact-support corrections documented; old confidence tests now supply true replay-verified indel rather than unrelated nearby positions. |
| Scorer equivalence/independence | Same-process scalar/bitvector swap; scalar oracle independent recurrence, baseline provenance and near-exact/>5kb controls requested. Exact sequence metric uses literal string equality. |
| File-size and ownership | Split classification summaries; preserve public helper imports; exclusive worker paths. Coordinator owns CLI/legacy scripts/integration assertion migration. |
| Ladder distal flank defect | Confirmed source: build_contig uses left[:500], trimming uses left[-20:]. No reference change silently introduced; known limitation explicitly retained, reference ablation still unmet. |
| Uncovered reference sequence | Valid concern; explicit unverified reference confidence implemented. Full coverage masking/reconstruction objective remains pending/unmet unless validated. |
| Diploid prior caused low QUAL | Causality unproven by supplied evidence; useful experimental hypothesis, no production prior change or causal claim. |
| Flat valley acceptance feasibility | Targeted flat regression passes; full development mapping comparison pending. No automatic promotion of pooled/anchor inference. |
| Generalization/assignment truth | Documented 53/53 QNAME collisions/no source tags, seed variation limitations. Fresh source capture and unseen designs planned before freeze. |

## Milestone 1 — received and reproduced

Actual Claude Fable 5.1 review completed successfully (modelUsage verified), using
an immutable tracked/untracked snapshot and manifest. No finding is closed by
review approval alone.

1. **Accepted HIGH evaluator alias mismatch.** Six failing regression cases first;
   allow only explicitly declared unresolved aliases, retain observed predictions
   and missing-allele denominator. Unknown key mismatches remain invalid. Details
   in milestone1-evaluator-fixes.md; 46 related tests pass.
2. **Partially accepted HIGH consensus measurement gap; causal claim rejected by
   reproduction.** Actual bcftools1.17 0/1 insertion produces ALT both without -H
   and with explicit -s SAMPLE -H I. Added real-tool regression; full cached-panel
   equivalence experiment pending. Heterozygous-indel projection is deliberately
   unavailable rather than guessed; its reason now appears per classification.
3. **Accepted MED flank-only evidence.** Two identical trimmed VNTRs no longer
   receive independent haplotype credit. Both sequences remain visible; real
   bcftools flank-only heterozygosity fixture passes.
4. **Accepted MED stale running failure.** External nonzero/timeout dominates
   incomplete sidecar, producing execution_failed. Running without completion is
   explicit invalid/incomplete, never success. Tested before and after fix.
5. **Accepted MED projection diagnostics.** Add per-allele vcf_projection status
   and reason: missing context, genotype ambiguity, REF/contig/trim/replay mismatch,
   overlapping records or unsupported alleles. Three observed failing regression
   cases now pass. The exact support rule remains conservative and unchanged.
6. **Accepted LOW provenance/report wording.** Left/right trim methods expose
   exact anchor versus fixed fallback or short untrimmed sequence. Report labels
   evidence state and alignment records. Strict valley limitation retained:
   genuine tied minima may be missed; all 53 exposed count pairs unchanged.

Coordinator additional finding: specificity previously omitted ambiguous normal
false positives from its denominator. A failing regression reproduced this.
Specificity now uses TN/(TN+FP), retaining positive alarms even when reconstruction
is ambiguous; unresolved negative controls remain separate and are not TN. Also
report confident negatives/all normal controls and resolved reconstruction call
rate. This fixes evaluator accounting, not caller accuracy.

Root focused checks: 31 consensus/projection units, 15 real bcftools phase tests,
14 evaluator scoring units, 10 report tests and 23 calling units passed. First
make test-int: 30 passed, 8 generated-data e2e deselected; new read-phasing tests
were still under construction and require a later complete run.

Milestone1 finding2 cached reproduction completed: all94 original HiFi/ONT allele
VCF/reference pairs are byte-identical under old default and new explicit IUPAC
commands, including17 heterozygous 0/1 indels across9 HiFi allele VCFs. Actual
bcftools1.17 commands/hashes saved. The predicted allele selection regression is
rejected for this matrix, not generalized to other versions or missing genotypes.
See milestone1-consensus-equivalence.md.

Pre-freeze complete gates after integration: make ci-check406 passed,91.87%
branch-aware coverage; make test-int47 passed,8 e2e deselected; explicit generated
selection6 passed,2 unchanged strict expected failures,47 integration deselected;
make docs-check passed. Earlier transient global failures arose from concurrent
incomplete files and newly introduced external hook mock boundaries, subsequently
resolved without changing expected scientific results. Final checks will run
again if further source changes are necessary.

## Independent internal integration audit (separate from Claude)

The bounded reviewer used actual reproductions, recorded in
independent-integration-audit.md. Accepted: prevent unproven duplicate copies from
inflating independent sequence recovery (retain literal equality explicitly),
exclude extra alleles from resolved call rate, require exact concordance status
for strict supported endpoints (retain legacy comparator), and reject a shared
explicit model for mixed-platform benchmark inventories. Tests first reproduced
these cases. Full input/reference existence checking now occurs within the run
status decorator, so invalid-file reruns replace stale completed status while
preserving Click exit2. General command-line syntax errors still require the
runner's explicit nonzero invocation record; no claim that unentered callbacks
can record execution. Experimental phaser read-list validation is being hardened.

No production parameter has changed from the phase-debugging ablations. The
singleton-floor proposal is not accepted on a downsampled read list: an actual
true30:2read fixture becomes14:1selected despite exact phase. Missing-allele and
no-call accounting must remain visible for any future support-floor policy.
