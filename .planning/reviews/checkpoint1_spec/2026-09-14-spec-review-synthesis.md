# Review Synthesis: 200-Dataset Experimental Specification & Acceptance Criteria

**Date:** 2026-09-14  
**Checkpoint:** 1 — Experimental Specification & Acceptance-Criteria Freeze  
**Reviewers:** 
- Claude Code (`claude-fable-5-1`, session `eb3dd2a4-8994-4c2e-855a-42962f65bc96`)
- Codex CLI (`gpt-6-astra`, thread `01a0a195-1b13-7783-b7bd-77d97cf3cf7b`)
**Consensus Verdict:** **REVISE AND RESUBMIT** (Clear, unanimous agreement on blocking defects and concrete remedies).

---

## 1. Executive Synthesis & Core Findings

Both external models conducted deep, independent, uninfluenced reviews of `.planning/2026-09-14-200-dataset-experimental-spec.md`. Their findings converged strongly on the primary structural risks, with complementary insights on statistical scoring and software contracts.

### 1.1 Unanimous Critical Dispositions & Action Items

| Issue | Claude Finding | Codex Finding | Risk / Consequence | Agreed Resolution / Remedy |
|---|---|---|---|---|
| **Absolute Gates vs Historical Baseline** | C1 (CRITICAL) | F01 (CRITICAL) | The proposed absolute gates (≥92% counts, ≥37.7% diploid sequence) were calibrated to the historical 44-sample HiFi set. The new 70-design development panel intentionally over-indexes on difficult cases (identical 60/60, gap 1–3, extreme asymmetry 25/140, low coverage). v0.11.0 baseline will be substantially lower on this panel. | **Re-anchor all gates to paired empirical deltas on the new panel:** Execute frozen v0.11.0 on all 140 dev datasets as Step 1. Require statistically defensible net paired improvement (discordant pairs win > loss) + strict zero-regression safety rules. Separate HiFi (N=70) and ONT (N=70) denominators. |
| **Identical Haplotype Contract Conflict** | C3 (CRITICAL), M3 (MAJOR) | F02 (CRITICAL) | Native `muconeup simulate --fixed-lengths 60 60` draws haplotypes independently from a Markov model, producing non-identical sequences. If identical sequences are generated, the codebase's strict evaluator (`scoring.py`) refuses duplicate credit without independent physical evidence, capping count accuracy at 90% and confident negatives at 68.2%. | (1) Use explicit authored structure files for identical sequence designs. (2) Update evaluation contracts to distinguish unique sequence discovery, homozygous single-sequence call (multiplicity unresolved), and independent two-haplotype resolution. For identical designs, the target is correct single sequence without false split, not false two-allele assertion. |
| **Sealed Final-Validation Isolation & Leakage** | C2 (CRITICAL) | F08, F09 (MAJOR) | Since `muconeup` is deterministic given seed and config, publishing `bio_seed` or generation commands in a public run ledger leaks the final ground truth to development agents. Sequential token assignment (`val_001..010`) also leaks normal controls. | **Dual-ledger architecture:** (a) Public execution ledger contains only randomized neutral tokens, platform, input FASTQ hash, run status, timings. (b) Sealed truth ledger contains seeds, truth sequences, mutation labels, and token permutations; sealed ledger SHA-256 is committed publicly before generation. Read records are shuffled so concat order does not leak haplotype. |
| **Seed Derivation & Collision Risk** | C4 (CRITICAL) | F03 (MAJOR) | `reads amplicon` derives PBSIM3 seeds as `S+1` (H1) and `S+2` (H2), and CCS seeds as `S+100*hap+bam`. Consecutive top-level seeds across designs share internal PBSIM seeds. | Establish disjoint seed domains with stride ≥ 1,000 for `bio_seed` (1,000,000+), `hifi_seed` (2,000,000+), and `ont_seed` (3,000,000+). Record derived component seeds in the sealed ledger. |
| **Boundary Repeat Mutations** | M7 (MAJOR) | F06 (MAJOR) | MUC1 fixed units (pre 1–5, after 6–9) do not permit most mutations due to parent unit restrictions (e.g. `insC_pos23` requires A/E; `ins16bp` requires C). Non-strict MucOneUp silently alters the parent repeat before mutating. | Define boundary mutations as either (a) terminal canonical repeats (total repeat indices 6–8 or L-6..L-4), or (b) explicit `delGCCCA` on fixed repeat 5/6. Strictly validate and record the exact parent repeat. |
| **Normal Specificity & Safety Gate** | M5 (MAJOR) | F11 (MAJOR) | "Maximum 1 false alarm" across 10 final controls is 90% specificity, which violates an absolute ≥95% gate. "Zero tolerance" vs "maximum 1" was ambiguous. | Re-specify: (1) Zero new false-positive mutation alarms relative to baseline on paired controls. (2) Absolute cap: ≤1 false alarm on dev HiFi (22 controls, ≥95.5%), 0 on final HiFi (10 controls, 100%). Report confident negative rate and unresolved controls separately. |
| **Minority Support & PCR Bias Tiers** | M2 (MAJOR) | F07 (MAJOR) | Default PCR bias at 200 requested templates yields ~4 usable 140-repeat reads in 25/140. Medium/low tiers yield 0–1 reads. Conditioning gates on "≥5 spanning reads" lets candidate failures escape the denominator. | Define minority strata by expected minority molecule tiers (low: 1–3, mid: 5–8, high: 10–15). Keep all intended datasets in unconditional denominators; report coverage-conditioned analysis as secondary. |
| **Inner Tuning / Selection Split** | M11 (MAJOR) | — | Factorial screening of parameter matrices across all 140 dev datasets risks overfitting. | Use a fixed 50/20 tuning/selection split within dev designs (100 tuning datasets, 40 selection datasets across platforms) or 5-fold cross-validation. |

---

## 2. Synthesis of Architectural Recommendations (Families A–F)

Both models provided rich, complementary guidance for the algorithmic investigations:

### Family A — Reference Architecture
- **Proximal Flank Defect:** Verified that the current ladder uses `flanking_left[:500]` (distal), whereas amplicon reads terminate adjacent to the VNTR (`flanking_left[-500:]`), forcing artificial soft-clipping and preventing biological left-flank anchor alignment.
- **Action:** Implement `flanking_left[-500:]`. Conduct a 2×2 factorial comparison: (distal vs proximal flank) × (ladder vs X-only discovery).
- **Coordinate Integrity:** Update trimming anchors, VCF projection, and coordinate mapping to ensure exact position preservation under the proximal reference.

### Family B — Iterative Allele-Specific References
- **Mechanism:** Hard EM with per-read realignment log-likelihood ratios, dominance margin $\delta$, minimum dominant read floor $m$, bounded iterations, and rollback if oscillation or minority collapse occurs.
- **Cross-Fitting:** Partition reads by deterministic hash (A/B) to assess stability and prevent error self-reinforcement.
- **Discriminating Test:** Challenge with 60/80 (10% minority) and 60/60 identical (verify convergence to a single reference without invented divergence).

### Family C — Length Inference & Molecule Assignment (Highest Priority)
- **Read-Dominance Scoring:** Realign reads against candidate references and evaluate per-read alignment score (AS) margins. In 60/60, 103/103 reads preferred 60 over 63; dominance scoring directly suppresses the false 63 candidate.
- **Minority Allele Retention:** Apply a length-implausibility test: four 8.4 kb reads cannot arise from a 25-repeat template via sequencing indels; retain them as a candidate even if below a global coverage threshold, with reconstruction flagged as low-support rather than discarded.
- **Error-Tolerant Flank Anchors:** Use edit-distance $k$-bounded anchor matching for spanning read discovery.

### Family D — Consensus and Phasing
- **Per-Haploid Calling:** After read assignment, each allele pool is haploid. Run Clair3 in haploid mode (`--haploid_precise` / `--haploid_sensitive`) to eliminate spurious heterozygous calls that trigger false splits.
- **Connected Component Phasing:** Require a single connected phase set across heterozygous variants before asserting phase; otherwise report per-block haplotypes or unresolved status. Never emit unphased IUPAC bases.
- **All-Read Assignment:** Use `whatshap haplotag` across all aligned reads rather than an internal downsampled subset.

### Family E — Technology Parameters
- **Minimap2:** Increase `-r` and `-z` so large indels do not fracture alignment chains; test `--secondary=yes -N 20 -p 0.9`.
- **Clair3:** Fix model compatibility; screen candidate depth, QUAL breakpoints, and amplicon-specific flags (`--var_pct_full`, `--var_pct_homo`).
- **Inner Selection:** Validate parameter choices on the 50/20 tuning/selection dev split to avoid simulator-specific overfitting.

### Family F — Performance Engineering
- **Whole-Pipeline Profiling:** Profile baseline v0.11.0 first; Clair3 dominates after exact bit-vector classification.
- **Concrete Optimizations:**
  1. Restrict Clair3 to the VNTR interval using a targeted BED file.
  2. Cache ladder `.mmi` minimap2 indexes keyed by reference SHA-256 and minimap2 version.
  3. Stream SAM/FASTQ conversions and compute contig statistics in a single pass without redundant disk writes.
  4. Ensure strict scientific output equivalence on benchmarked samples.

---

## 3. Discriminating Experimental Plan & Order of Operations

Based on the synthesized review recommendations, the implementation will proceed in this strict, non-circular order:

1. **Amend Specification & Publish Allocations:**
   - Update `.planning/2026-09-14-200-dataset-experimental-spec.md` with all C1–C4, M1–M13, F01–F12 resolutions.
   - Author the explicit 100-design allocation table with seeded stratified split (70 dev, 30 final).
   - Implement the dual-ledger protocol and sealed custodian commitment.
2. **Execute Separate Pilot Validation:**
   - Run a 4-dataset pilot (2 designs × 2 platforms) including one authored structure, one extreme gap, and one low-coverage case to verify all generation, capture, and truth-verification mechanics.
3. **Generate 140 Development Datasets & Seal 60 Final Datasets:**
   - Generate biological truth once per design; simulate HiFi and ONT reads with disjoint seeds.
   - Seal the 60 final validation datasets; verify custodian hash.
4. **Establish Empirical Baseline on Development Panel:**
   - Execute frozen v0.11.0 on all 140 development datasets. Record exact sample-level failure atlas.
5. **Causal Stage Diagnosis & Oracle Ceilings:**
   - Test truth read assignment and truth reference oracles in dev to measure theoretical upper bound.
6. **Implement & Test Algorithmic Improvements (Families A, C, D, B, E, F):**
   - Family A: Proximal flank reference ablation.
   - Family C: Read-dominance candidate selection to fix 60/60 and 25/140.
   - Family D: Haploid Clair3 calling per assigned allele + connected-component phasing.
   - Family B: Bounded iterative refinement with minority protection.
   - Family E: Parameter optimization on the inner tuning split.
   - Family F: Targeted performance optimizations.
7. **Freeze Candidate & Verify Development Gates:**
   - Run complete candidate pipeline on all 140 dev datasets; verify paired improvement against v0.11.0 baseline.
8. **Unseal & Evaluate Protected Final-Validation Panel:**
   - Execute baseline and candidate on the 60 protected final datasets without tuning.
9. **Final Reviews & Documentation:**
   - Claude Fable 5.1 and Codex GPT-6 Astra reviews on final diff and report.
