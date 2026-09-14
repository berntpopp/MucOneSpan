# Independent Review Request: Experimental Specification & Acceptance Criteria Freeze

**Target Artifact:** `.planning/2026-09-14-200-dataset-experimental-spec.md`  
**Repository State:** MucOneSpan branch `improve/simulation-200` (at commit `9535f7ee02033cda3da9b22c5b73752e5473b016`, tag `v0.11.0`)  
**Context:** 
We are designing, generating, evaluating, and improving MucOneSpan using 200 new MucOneUp sequencing datasets (100 distinct biological designs × 2 platforms: PacBio HiFi and ONT).
The 200 datasets are pre-split by biological design into 140 development datasets (70 designs × 2 platforms) and 60 protected final-validation datasets (30 designs × 2 platforms).

Please perform an exceptionally rigorous bioinformatics, statistical, and software engineering review of the specification.

## Document Under Review: `.planning/2026-09-14-200-dataset-experimental-spec.md`

Please inspect the file at `.planning/2026-09-14-200-dataset-experimental-spec.md`.

## Specific Review Questions

1. **Dataset Definition & Randomness:**
   - Does generating each biological truth once and feeding the exact same diploid FASTA into HiFi and ONT read simulation (with distinct sequencing error seeds) properly satisfy biological pairing while preserving sequencing error independence?
   - Are the seed policies (`bio_seed`, `hifi_seed`, `ont_seed`) sufficient to prevent artifactual cross-contamination?

2. **Stratification & Coverage Matrix:**
   - Are the 32 mutation-negative controls (22 dev, 10 final), 13 supported mutation types across both haplotypes, boundary-repeat mutations, equal lengths (identical and variant), repeat gaps (1, 2, 3, wider, asymmetric), and coverage tiers well-specified?
   - Are there hidden biases or gaps in the stratification matrix?

3. **Protection of Final Validation:**
   - Is the protocol for blinding the caller and development agents with neutral identifiers and sealing final truth/reads until post-candidate freeze statistically robust?

4. **Predeclared Numerical Acceptance Criteria:**
   - Are the primary scientific endpoints (diploid sequence recovery, exact count pairs, mutation detection/annotation, normal specificity safety gate) defensible and ambitious yet realistic?
   - Are the challenging strata gates (60/60 identical, 25/140 extreme asymmetry, 60/61 gap 1) properly formulated?

5. **Competing Algorithmic Families (A–F) Brainstorming:**
   - What specific failure modes, mathematical formulation issues, or architectural risks exist for:
     - Family A (Reference architecture: proximal flank vs distal flank, X-only vs ladder)
     - Family B (Iterative allele-specific reference refinement)
     - Family C (Length inference & molecule assignment: replacing indel-valleys with span-evidence + read-dominance scoring)
     - Family D (Consensus and phasing: WhatsHap linkage, phase-set connectivity, IUPAC resolution)
     - Family E (Platform-specific parameter matrices for minimap2, Clair3, WhatsHap)
     - Family F (Pipeline performance optimization)
   - What discriminating experiments should be prioritized?

## Required Review Response Format

Please structure your review as:
1. **Summary Verdict:** (APPROVED / APPROVED WITH AMENDMENTS / REVISE AND RESUBMIT)
2. **Detailed Findings:** Each finding formatted as:
   - **Finding ID & Severity:** [CRITICAL | MAJOR | MINOR | ADVISORY]
   - **Location:** (Section / Line)
   - **Failure Scenario / Risk:** What breaks, why, and under what conditions.
   - **Suggested Validation / Remedy:** Concrete recommendation.
3. **Independent Brainstorming & Architectural Recommendations:** Concrete technical proposals for Families A–F.
