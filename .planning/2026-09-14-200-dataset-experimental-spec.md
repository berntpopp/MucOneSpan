# MucOneSpan 200-Dataset Experimental Specification & Acceptance Criteria (Frozen Revision)

**Date:** 2026-09-14  
**Status:** Frozen Revision 2 (Post External Review Synthesis)  
**Baseline Comparator:** MucOneSpan v0.11.0 (`9535f7ee02033cda3da9b22c5b73752e5473b016`)  
**Simulator:** MucOneUp 0.44.5 (`c60618ed6d4786591ad7f9283a72b7d1a5504625`)  
**Authoritative Scope:** 100 distinct diploid biological designs × 2 platforms (PacBio HiFi and ONT) = 200 sequencing datasets (140 development, 60 protected final-validation).

---

## 1. Experimental Invariants & Design Principles

### 1.1 Core Biological Invariants
1. **100 Distinct Biological Designs:** Each design declares a diploid biological truth.
2. **Single Truth Generation:** Biological truth is simulated exactly once per design via `muconeup simulate`.
3. **Paired Platform Sequencing:** The identical diploid FASTA is piped to PacBio HiFi (`reads amplicon --platform pacbio`) and ONT (`reads amplicon --platform ont`).
4. **Truth Identity Verification:** Bitwise SHA-256 checksums verify that HiFi and ONT simulations originate from identical biological truth sequences.
5. **Deterministic PCR Allocation:** Simulator PCR allocation runs in deterministic mode (`pcr_bias.stochastic: false`), ensuring identical initial template allocations per allele between platforms for a given design.
6. **Disjoint Seed Domains:** To prevent internal PBSIM3/CCS seed collisions (`S+1`, `S+2`, `S+100*hap`), top-level seeds use non-overlapping ranges with stride ≥ 1,000:
   - Biological truth generation: `bio_seed = 1_000_000 + design_id * 10`
   - HiFi read simulation: `hifi_seed = 2_000_000 + design_id * 10`
   - ONT read simulation: `ont_seed = 3_000_000 + design_id * 10`

### 1.2 Split Structure & Protection Protocol
- **Development Split:** 70 biological designs × 2 platforms = 140 sequencing datasets (IDs `dev_001` through `dev_070`).
- **Protected Final-Validation Split:** 30 biological designs × 2 platforms = 60 sequencing datasets (IDs `val_001` through `val_030`).
- **Inner Development Split:** Development designs are partitioned into a 50-design tuning set (100 datasets) and a 20-design selection set (40 datasets) for hyperparameter evaluation.
- **Dual-Ledger Architecture:**
  - *Public Execution Ledger (`ledger_public.jsonl`):* Contains only randomized neutral sample tokens, platform, input FASTQ SHA-256, run status, wall/CPU timings, and peak memory.
  - *Sealed Truth Ledger (`ledger_sealed.jsonl`):* Contains design IDs, biological seeds, sequencing seeds, true lengths, mutation identities, truth FASTA hashes, and the token permutation mapping.
  - *Cryptographic Commitment:* The SHA-256 of the sealed ledger is published before any dataset generation.
  - *Caller Blinding:* The caller and development agents access only caller-facing read files with randomized neutral tokens. Final truth is sealed until candidate freeze.

---

## 2. 100-Design Stratification Matrix

The 100 biological designs are strictly allocated into disjoint strata:

| Category | Dev Count (70) | Final Count (30) | Total (100) | Detailed Design Specification |
|---|---|---|---|---|
| **Mutation-Negative Controls (Equal Identical)** | 7 | 3 | 10 | Authored identical repeat chains: 30/30, 40/40, 50/50, 60/60, 70/70, 80/80, 100/100. Evaluates false-splitting and over-clustering. Target: single sequence negative without false split. |
| **Mutation-Negative Controls (Equal Variant)** | 5 | 2 | 7 | Equal repeat counts with distinct Markov-drawn repeat compositions (e.g. 50/50, 60/60, 80/80). Evaluates composition discrimination at equal length. |
| **Mutation-Negative Controls (Heterozygous Gaps)** | 10 | 5 | 15 | Normal controls across gaps: Gap 1 (60/61, 40/41), Gap 2 (60/62, 50/52), Gap 3 (60/63), typical (40/50, 60/80), asymmetric (25/100). |
| **Total Negative Controls** | **22** | **10** | **32** | ≥30 controls constraint satisfied (32 total). |
| **Equal Lengths with Mutation** | 7 | 3 | 10 | Same repeat count (e.g. 60/60) with mutation on H1 or H2 (`dupC`, `dupA`, `insG`, `insCCCC`, etc.). |
| **Gap 1 Alleles (Mutation Positive)** | 4 | 2 | 6 | Repeat count difference of 1 (e.g. 50/51, 60/61, 75/76) with mutation. |
| **Gap 2 Alleles (Mutation Positive)** | 4 | 2 | 6 | Repeat count difference of 2 (e.g. 40/42, 60/62, 80/82) with mutation. |
| **Gap 3 Alleles (Mutation Positive)** | 4 | 2 | 6 | Repeat count difference of 3 (e.g. 45/48, 60/63, 70/73) with mutation. |
| **Wider Heterozygous Gaps (Mutation Positive)** | 15 | 6 | 21 | Typical gaps: 30/50, 40/60, 50/70, 60/80, 70/90, 80/100, 100/120. |
| **Extreme Asymmetry (Minority Tiers)** | 10 | 5 | 15 | 25/140, 25/100, 30/120, 20/90, 35/130. Stratified across expected minority molecule tiers (Tier 1: 1–3, Tier 2: 5–8, Tier 3: 10–15). |
| **Length Extremes** | 4 | 0 | 4 | Short (20–25 repeats) and very long (>120 repeats: 120/140, 130/150). |
| **Total Mutation Positive** | **48** | **20** | **68** | Covers all 13 supported mutations on both H1 and H2. |
| **Total Biological Designs** | **70** | **30** | **100** | Exactly 100 designs. |

### 2.1 Mutation Coverage & Boundary Rules
- **All 13 Supported Types Covered:** `dupC` (10), `dupA` (6), `insG` (6), `insCCCC` (6), `insC_pos23` (5), `insG_pos58` (5), `insG_pos54` (5), `insA_pos54` (5), `del18_31` (5), `ins16bp` (5), `ins25bp` (5), `delinsAT` (5), `delGCCCA` (5).
- **Haplotype Balance:** Mutated on H1 (34 designs), Mutated on H2 (34 designs).
- **Boundary Repeat Rule:** 16 designs test mutations in boundary regions:
  - 12 designs place mutations in terminal variable canonical repeats (total repeat indices 6–8 or L-6..L-4).
  - 4 designs test `delGCCCA` on fixed repeat 5/6 (biologically validated constant-region boundary).
- **Strict Parent Unit Enforcement:** All mutation targets are validated to match the parent unit requirements of the mutation definition (e.g. `insC_pos23` on A/E, `ins16bp` on C). Permissive silent parent substitution is forbidden.

### 2.2 Coverage Tiers
- **High Tier (200 requested templates):** ~90–160 usable reads. Standard evaluation of diploid resolution and sensitivity.
- **Moderate Tier (60 requested templates):** ~25–45 usable reads. Tests robustness at typical clinical amplicon depths.
- **Low Tier (20 requested templates):** ~8–15 usable reads. Explicitly evaluates the `insufficient_evidence` / no-call contract and minority dropout safety.

---

## 3. Predeclared Numerical Acceptance Criteria

All acceptance criteria are evaluated strictly on paired runs of the candidate algorithm versus frozen baseline v0.11.0 on the **newly generated datasets**, with HiFi (N=70 dev, N=30 final) and ONT (N=70 dev, N=30 final) evaluated and reported separately.

### 3.1 Primary Scientific Gates (Development Panel: N=70 HiFi, N=70 ONT)
1. **Paired Diploid Sequence Reconstruction (Primary Scientific Endpoint):**
   - Candidate must achieve a statistically significant net paired improvement over v0.11.0 baseline:
     $$\text{Discordant Wins} - \text{Discordant Losses} \ge 8 \text{ on HiFi}$$
     $$\text{Discordant Wins} - \text{Discordant Losses} \ge 5 \text{ on ONT}$$
   - HiFi exact diploid sequence recovery must increase by at least +12.0 percentage points over the new-panel v0.11.0 baseline.
2. **Exact Allele Count Pairs:**
   - Candidate must achieve a net positive paired delta: $\text{Wins} > \text{Losses}$.
   - Zero regression on previously exact distinct-length pairs in the baseline.
3. **Mutation Detection & Annotation:**
   - Overall mutation sensitivity must meet or exceed baseline sensitivity on both platforms.
   - Exact annotation accuracy (matching parent, name, 1-based total-repeat position, and assigned haplotype) must increase by at least +5 percentage points over baseline.
   - Extra mutation calls on positive samples (Annotation False Discovery Rate) must decrease relative to baseline.
4. **Normal Control Specificity (Strict Safety Gate):**
   - **Zero new false-positive mutation alarms** on normal controls relative to baseline.
   - Absolute false alarm cap: $\le 1$ false alarm across the 22 dev HiFi controls ($\ge 95.5\%$ specificity); $\le 1$ across 22 dev ONT controls.
   - Confident-negative rate $\ge 65\%$ on controls, with unresolved negatives reported transparently.

### 3.2 Challenging Strata Gates
1. **Homozygous Same-Length (Identical 60/60):**
   - Correctly report single-sequence negative with multiplicity unresolved (no false split into 59/60 or 60/63) in $\ge 85\%$ of identical cases.
2. **Extreme Asymmetry (25/140):**
   - Detect both lengths in $\ge 60\%$ of cases with $\ge 5$ usable minority reads.
   - In low-support minority cases (<3 reads), emit explicit low-coverage/unresolved status without manufacturing a false second allele.
3. **Gap 1 Resolution (60/61):**
   - Achieve correct separate length detection without extra alleles in $\ge 80\%$ of HiFi Gap 1 cases having $\ge 20$ usable reads per allele.

### 3.3 Protected Final-Validation Confirmation Rules (N=30 HiFi, N=30 ONT)
- Final validation is a confirmatory execution under sealed conditions.
- **Pass Rule:** The candidate point estimate must exceed baseline for diploid sequence recovery and count pairs on both platforms.
- **Safety Rule:** Zero new false-positive alarms on the 10 final normal controls (100% specificity).
- At most 2 net regressions across the 30 designs.

### 3.4 Engineering & Performance Gates
1. **Runtime Efficiency:**
   - Whole-pipeline wall time must not increase by more than 15% under matched thread budgets.
   - Alternating paired runs (N ≥ 5) measuring median wall time, user CPU, system CPU, and maximum process-tree RSS.
2. **Resource Bound:**
   - Peak RSS $\le 8$ GB per pipeline execution.
3. **Software Contracts:**
   - `make ci-check` passing with 100% unit tests passed and $\ge 80\%$ branch-aware coverage.
   - Author line-length gate: strictly $\le 649$ physical lines per authored file.
   - All parameters configuration-driven via typed immutable settings.

---

## 4. Competing Algorithmic Hypotheses (Families A–F)

1. **Family A (Reference Architecture):**
   - *Hypothesis:* The current ladder uses distal left flank `flanking_left[:500]`, creating an artificial mismatch with amplicon reads that terminate adjacent to the VNTR (`flanking_left[-500:]`). Replacing with proximal flank eliminates artificial soft-clipping and establishes biological trimming anchors.
   - *Test:* 2×2 comparison (distal vs proximal flank) × (ladder vs X-only discovery reference).
2. **Family B (Iterative Allele-Specific References):**
   - *Hypothesis:* Two-pass refinement with competitive realignment of all eligible reads resolves internal repeat mismatch without self-reinforcing errors if governed by a dominance margin $\delta$ and minimum dominant read floor $m$.
   - *Test:* 60/80 with 10% minority (verify minority retention); 60/60 identical (verify single reference convergence).
3. **Family C (Length Inference & Molecule Assignment):**
   - *Hypothesis:* Replacing crude indel-valley slicing with per-read realignment score (AS) dominance resolves 60/60 false splitting (reads strictly prefer 60 over 63) and preserves coherent minority reads in 25/140.
   - *Test:* Realignment dominance test on 60/60 and 25/140 cached BAMs.
4. **Family D (Consensus and Phasing):**
   - *Hypothesis:* Running Clair3 in haploid mode (`--haploid_precise` / `--haploid_sensitive`) after read assignment eliminates spurious heterozygous calls. Enforcing single connected phase set across heterozygous variants prevents false phasing.
   - *Test:* Phased synthetic cis/trans mixtures; verify 0 false events induced on collapsed clusters.
5. **Family E (Technology Parameters):**
   - *Hypothesis:* Tuning minimap2 chaining/gap penalties (`-r`, `-z`, `--secondary=yes -N 20 -p 0.9`) and Clair3 amplicon thresholds on the 50/20 inner dev split improves indel sensitivity without overfitting.
   - *Test:* Factorial screening on inner tuning split; evaluate on inner selection split.
6. **Family F (Performance Engineering):**
   - *Hypothesis:* Restricting Clair3 to the VNTR interval via BED, caching `.mmi` ladder index, and streaming SAM/FASTQ processing reduces pipeline runtime by 20–35% while preserving 100% scientific output equivalence.
   - *Test:* Paired alternating timings with output hash equivalence verification.
