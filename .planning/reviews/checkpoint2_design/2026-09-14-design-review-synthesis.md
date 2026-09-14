# Checkpoint 2 Architectural Design Review Synthesis

**Date:** 2026-09-14  
**Subject:** Architectural Design Specification for MucOneSpan Reconstruction  
**Reviewers:**
1. Claude Code CLI (`claude-fable-5-1`, 667.2s, exit code 0)
2. Codex CLI (`gpt-6-astra`, 461.5s, exit code 0)
**Empirical Baseline:** 140 development datasets evaluated under frozen v0.11.0 with active Clair3 v1.0.10 neural network models.

---

## 1. Executive Summary & Consensus

Both independent external reviewers conducted exhaustive, code-level and empirical reviews of the proposed architecture. The consensus across Claude Code, Codex CLI, and our empirical baseline investigation is clear:
1. **Empirical Baseline Established:** The initial run revealed an environment defect (PATH ordering picking up `env_pacbio`'s Python 3.9 without TensorFlow for Clair3) and a metric key mismatch in the diagnosis script. Both were immediately repaired. The verified v0.11.0 baseline across all 140 development datasets achieves **14/140 exact reconstructions** (4/70 HiFi, 10/70 ONT), with length inference being the primary failure mode (82 datasets: 43 HiFi, 39 ONT) followed by variant calling/annotation (44 datasets: 23 HiFi, 21 ONT).
2. **Family A (Flank Correction):** The switch from distal `flanking_left[:flank_length]` to proximal `flanking_left[-flank_length:]` is mathematically and biologically correct for the Wenzel 2018 amplicons. However, because amplicon reads span only ~6 bp of left flank and ~24 bp of right flank, the flank cancels across contigs during initial ladder mapping. Its critical value is in **exact anchor consensus trimming** and preventing silent fallback to fixed coordinates. The reference ladder, consensus trimming, and configuration guard must be updated in lockstep.
3. **Family C (Read-Dominance Length Inference):** Naive secondary alignment margin testing is mathematically circular when conditioned on argmax secondary contigs. Instead, candidate discovery must be separated from pairwise rescoring:
   - Identify candidate lengths from spanning reads.
   - Perform pairwise rescoring of all spanning reads against candidate pairs using an affine-gap margin $\delta = 0.5 \cdot g$ (~42-43).
   - Require a minimum support count of dominant spanning reads ($|D_2| \ge 3$ HiFi, $\ge 4$ ONT) and a low ambiguous fraction. If the margin or support floor is not met, report a single candidate with multiplicity unresolved rather than forcing a false split.
4. **Family D (Consensus & Phasing):** Clair3's `--haploid_precise` discards heterozygous calls and emits GT `1`, which causes `phasing.py:phase_evidence` to fail as `non_diploid`. Partitioning reads first by dominance, validating partition purity, and applying genotype and depth checks preserves sensitivity and avoids IUPAC masking while preventing spurious split alleles.
5. **Modularity & 649-Line Enforcement:** To prevent bloat in `alleles.py` (currently 597 lines) and `calling.py` (447 lines), new functionality is partitioned into clean, dedicated modules:
   - `src/muc_one_span/length_candidates.py`
   - `src/muc_one_span/read_dominance.py`
   - `src/muc_one_span/read_partition.py`

---

## 2. Reviewer Findings & Dispositions

| Area | Reviewer Observation (Claude / Codex) | Evaluation & Disposition | Action |
|---|---|---|---|
| **Baseline Environment** | Both noted Clair3 pileup logs failed with `ModuleNotFoundError: No module named 'tensorflow'` due to PATH order, and diagnosis read incorrect `eval_row` keys. | **Valid & Critical.** Fixed environment PATH prepending `env_clair3/bin` and corrected metric schema to read `eval_row["metrics"]`. Full rerun across 140 dev sets completed. | Baseline re-established and locked. |
| **Family A: Proximal Flank** | Reads carry only ~6 bp left / ~24 bp right flank. Distal vs proximal cancels in initial mapping, but is vital for exact anchor trimming. Trimming anchor in `consensus.py` must match. | **Valid.** The change is biologically necessary and required for exact anchor trimming. Update `ladder.py` and `consensus.py` together; update bundled reference ladder. | Implement proximal flank in `ladder.py` and `consensus.py`. |
| **Family C: Dominance Test** | Testing margins on reads already assigned to secondary contigs is circular. Mean margin $\Delta AS$ can be skewed by a single outlier. Minimap2 `-N 10 -p 0.8` leaves scores missing. | **Valid.** Replace with pairwise rescoring of spanning reads against candidate pairs. Count distinct dominant reads ($S(r, c_2) - S(r, c_1) > \delta$) rather than average margin. | Implement pairwise candidate rescoring in `read_dominance.py`. |
| **Family C: Margin Threshold $\delta$** | Universal raw AS margin is invalid across platforms and lengths. $\delta$ must derive from the affine gap minimum $g$: $\delta = 0.5 \cdot g$ (43 for HiFi, 42 for ONT). | **Valid & Principled.** Derives directly from minimap2 scoring presets (`map-hifi` and `lr:hq`). | Use $\delta = 43$ (HiFi) and $\delta = 42$ (ONT) derived from gap cost. |
| **Family D: Haploid Calling** | Clair3 `--haploid_precise` outputs GT `1` (which breaks `phase_evidence`), drops real variants below 100% AF, while `sensitive` keeps 0/1. | **Valid.** Do not rely blindly on `--haploid_precise`. Use read partitioning, filter variants by allele depth/fraction and purity, and format VCF records for consensus replay. | Implement partition-aware variant calling and consensus. |
| **Family E: Filtering & QUAL** | Raising min QUAL to 8 or 6 drops true mutations (e.g. dupC at QUAL 3.22). Gain precision from partition purity instead. | **Valid.** Keep min QUAL $\le 5.0$. Eliminate false positives through read partition purity and spanning filters. | Maintain `min_qual = 5.0`. |
| **Architecture & Line Gate** | `alleles.py` is at 597 lines; adding rescoring directly will breach the 649-line gate. | **Valid.** Extract candidate discovery and read dominance into dedicated modules. | Create `length_candidates.py` and `read_dominance.py`. |

---

## 3. Implementation Plan for Candidate Pipeline

### Phase 1: Core Modular Components
1. **`src/muc_one_span/length_candidates.py`:**
   - Candidate length discovery from spanning read spans and primary/secondary contig clusters.
   - Handles equal-length plateau detection and gap 1–3 candidate generation.
2. **`src/muc_one_span/read_dominance.py`:**
   - Pairwise alignment scoring of spanning reads against candidate contigs ($c_1, c_2$).
   - Computes $D_1 = \{r : S_1 - S_2 > \delta\}$, $D_2 = \{r : S_2 - S_1 > \delta\}$, $A = \{r : |S_1 - S_2| \le \delta\}$.
   - Applies decision rule: if $|D_2| \ge m_{\min}$ and $|D_2| / (|D_1| + |D_2|) \ge \theta_{\text{dom}}$, confirm two alleles; otherwise retain single candidate with multiplicity unresolved.
3. **`src/muc_one_span/ladder.py` & `consensus.py`:**
   - Use proximal left flank `flanking_left[-flank_length:]`.
   - Update left trimming anchor to proximal sequence in lockstep.
   - Assert `exact_anchor` trimming method and validate coordinates.

### Phase 2: Verification and Quality
- Unit tests for all new modules with deterministic mock fixtures.
- Verify line count $\le 649$ physical lines for every authored file.
- Run `make ci-check` ($\ge 80\%$ branch coverage).

### Phase 3: Empirical Candidate Evaluation on 140 Dev Datasets
- Run candidate evaluation pipeline.
- Verify paired criteria against baseline:
  - $\text{Wins} - \text{Losses} \ge +8$ (HiFi)
  - $\text{Wins} - \text{Losses} \ge +5$ (ONT)
  - Zero baseline regressions on previously correct samples.
  - Zero false positive mutation alarms on normal controls.
