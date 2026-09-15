# Checkpoint 3 Implementation Review Synthesis: Candidate Implementation & Development Benchmark

**Date:** 2026-09-15  
**Reviewers:**
- Claude Code CLI (`claude-fable-5-1`, high-effort plan mode)
- Codex CLI (`gpt-6-astra`, read-only sandbox mode)
**Coordinating Agent:** Antigravity Senior Bioinformatician / Software Architect

---

## 1. Executive Summary & Verdict

Both independent external reviewers delivered a unanimous verdict:
**VERDICT: REVISE AND RESUBMIT**

While both reviewers confirmed substantial, regression-free improvements on the PacBio HiFi platform (HiFi net paired exact-sample reconstruction: +14, exceeding the $\ge +8$ gate; exact sequence reconstruction net: +21; annotation FDR reduced from 44.9% to 25.6%), both reviewers identified critical blockers that preclude opening the protected 60-dataset validation split:

1. **Benchmark Input Defect (Duplicate Read Names in Simulated FASTQs):**
   Every development (and simulated) FASTQ contains duplicate read names (`@S/1/ccs` ... `@S/N/ccs` in HiFi; `@S_1` ... `@S_N` in ONT) generated independently for each haplotype template during amplicon simulation. These colliding identifiers cause score merging across distinct physical molecules in read-dominance scoring.
2. **Read-Dominance Missing-Score Inference Flaw:**
   Treating a missing alignment on $c_1$ as decisive dominance for $c_2$ ($d_2 += 1$) failed on `ctrl_ident_80_80_hifi`, accepting a spurious 75-repeat candidate because minimap2 suppressed secondary alignments on $c_1$.
3. **ONT Reconstruction Gate Missed:**
   ONT net exact reconstructions achieved only +1 (gate requires $\ge +5$). The primary failure mode is a systematic $+1$ repeat unit overcall on alleles $\ge 80$ repeats (e.g. calling 81 instead of 80, or 101 instead of 100) due to ONT indel noise shifting the weighted cluster center.
4. **Haploid VCF Re-Genotyping Issues:**
   `vcf.py` silently clamped user QUAL to `min(min_qual, 4.0)` rather than exposing a configurable parameter, crashed on multiallelic `AF` strings (`"0.5,0.4375"`), and contributed to four lost true-positive mutation events (`mut_eq_inscccc_60_60_h1_hifi`, `mut_asym_ins16bp_20_90_h2_ont`, `mut_gap1_insg_50_51_h1_ont`, `mut_gap3_dupc_60_63_h1_ont`).
5. **Config & Architecture Contract:**
   Heuristic thresholds (dominance margins, read floors, AF thresholds) were hardcoded constants instead of typed immutable settings. Production `alleles.py` did not invoke the functions in `length_candidates.py`.

---

## 2. Reviewer Scorecard & Gate Audit

| Frozen Gate (Spec §3) | Required | Baseline (140 Dev) | Candidate (140 Dev) | Claude Finding | Codex Finding | Adjudication |
|---|---|---|---|---|---|---|
| HiFi Paired Exact Reconstruction | $\ge +8$ | 4/70 (5.7%) | 18/70 (25.7%) | +14 (14W, 0L) — **Pass** | +14 (14W, 0L) — **Pass** | **PASS** (+14 net) |
| ONT Paired Exact Reconstruction | $\ge +5$ | 9/70 (12.9%) | 10/70 (14.3%) | +1 (1W, 0L) — **Fail** | +1 (1W, 0L) — **Fail** | **FAIL** (+1 net, needs +4 more) |
| HiFi Exact Recovery Gain | $\ge +12$ pp | 5.7% | 25.7% | +20.0 pp — **Pass** | +20.0 pp — **Pass** | **PASS** (+20.0 pp) |
| Exact Count Pairs Wins > Losses | Wins > Losses | 42/140 | 48/140 | 6W, 0L — **Pass** | 6W, 0L — **Pass** | **PASS** (6W, 0L, 0 regressions) |
| Mutation Alarm Sensitivity | $\ge$ Baseline | 72/96 (75.0%) | 69/96 (71.9%) | 72→69 — **Fail** | 72→69 — **Fail** | **FAIL** (Lost 4 TP events) |
| Exact Annotation Accuracy | $\ge +5$ pp | 53/96 (55.2%) | 53/96 (55.2%) | 0 pp — **Fail** | 0 pp — **Fail** | **FAIL** (Flat: 4 lost, 4 gained) |
| Control False Alarms | $\le 1$ per plat. | HiFi: 2, ONT: 0 | HiFi: 2, ONT: 0 | HiFi 2, ONT 0 — **Fail** | HiFi 2, ONT 0 — **Fail** | **FAIL** (`ctrl_ident_80_80_hifi` spurious split) |
| Zero New Control Alarms | 0 new | — | 0 new | Pass (9→2 alarms) | Pass (9→2 alarms) | **PASS** |
| Confident Negative Control Rate | $\ge 65\%$ | 11.4% (5/44) | 25.0% (11/44) | 25.0% — **Fail** | 22.7% / 27.3% — **Fail** | **FAIL** (Too many unresolved negatives) |
| File Size Gate | $< 650$ lines | Max 649 | Max 595 | Verified | Verified | **PASS** |
| Software & CI Tests | $\ge 80\%$ br. cov | Pass | 92.3% br. cov | Verified static | Verified static | **PASS** |

---

## 3. Actionable Remediations & Technical Plan

To resolve the reviewers' objections and meet all acceptance criteria before opening the protected split:

### Action 1: Benchmark Read-Identifier Sanitization
- Modify the FASTQ generation / preparation in `scripts/run_200_experiment.py` (and input preparation) so every read record receives a globally unique, deterministic, neutral identifier (e.g. `@{sample_id}_{haplotype}_{read_idx}` or `@{qname}_r{idx}`).
- Verify sequence and Phred quality score identity before/after renaming.
- Regenerate development and final FASTQs and update manifests.

### Action 2: Robust Read-Dominance & Secondary Suppression Handling
- In `src/muc_one_span/read_dominance.py`:
  - A read observed *only* on $c_2$ cannot be assumed dominant over $c_1$ without verification. If $c_1$ is the primary contig, missing $c_1$ in the ladder BAM may be due to minimap2 secondary alignment suppression (`-N 5` default limit).
  - Treat missing scores on either contig as *incomparable / non-decisive* by default, unless verified by explicit pairwise alignment or alignment score bounds.
  - Require that any candidate $c_2$ have genuine, decisive, pairwise-validated read support ($AS(c_2) - AS(c_1) \ge \Delta$) with physical reads that were mapped against both candidates.
  - Add `ctrl_ident_80_80_hifi` rejection regression test to verify spurious 75-repeat allele is 100% rejected.

### Action 3: ONT Length Noise & Systematic $+1$ Bias Resolution
- Analyze ONT read alignment distributions in `alleles.py`:
  - Why does ONT systematically call 81 instead of 80 on long alleles?
  - Examine the clustering and peak selection: In long ONT reads ($\ge 80$ VNTRs, $\approx 5$ kb), homopolymer deletions and indel errors cause minimap2 alignments to truncate or under/over-estimate repeat unit counts.
  - Evaluate discrete mode / integer rounding vs weighted average on ONT clusters, or evaluate peak contig assignment using modal repeat count rather than noisy weighted means.

### Action 4: Haploid VCF Genotyping Refinements
- In `src/muc_one_span/vcf.py`:
  - Handle multiallelic `AF` strings safely by parsing comma-separated lists and evaluating the primary alternative allele.
  - Remove the silent `min(min_qual, 4.0)` override. Expose `haploid_min_qual: float = 4.0` in `ConsensusSettings` as a typed setting.
  - Investigate and recover the 4 lost mutation events (`mut_eq_inscccc_60_60_h1_hifi`, `mut_asym_ins16bp_20_90_h2_ont`, `mut_gap1_insg_50_51_h1_ont`, `mut_gap3_dupc_60_63_h1_ont`).

### Action 5: Architecture & Settings Compliance
- Move all dominance thresholds (`score_margin_hifi`, `score_margin_ont`, `min_dominant_reads_hifi`, `min_dominant_reads_ont`, `min_dominance_ratio`) into `AlleleSettings` in `src/muc_one_span/settings.py`.
- Wire or clean up `src/muc_one_span/length_candidates.py` to ensure all exported functions are either integrated into production or cleanly scoped.
- Ensure runtime benchmarking is formally measured and reported within the $\le +15\%$ gate.
