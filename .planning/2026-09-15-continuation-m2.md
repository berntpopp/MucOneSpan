# Continuation Record: MucOneSpan Scientific Reliability & Clinical Reporting
**Timestamp:** 2026-09-15T11:22:00+02:00  
**Active Revision (git sha):** `948eb3a` (working tree dirty with M1 implementations)  
**Branch:** `feat/scientific-reliability-reporting`  
**Current Milestone:** Milestone M2 (Scientific Deconvolution & Baseline Failure Atlas)

---

## 1. Dirty-State & File Summary
Modified tracked files:
- `scripts/run_500_experiment.py`: Bounded concurrency, staging isolation, atomic rename, SHA256 copy verification, test split unseal guard.
- `scripts/run_eval_pipeline.py`: Split filtering, frozen denominator accounting, worker isolation, resume progress logging (`eval_progress.jsonl`).
- `src/muc_one_span/evaluation/artifacts.py`: Added `interrupted` status support, exit code signal translation.
- `src/muc_one_span/run_status.py`: Added `KeyboardInterrupt` handling emitting `interrupted` status.
- `src/muc_one_span/tools.py`: Added process group registry, timeouts, and `killpg(SIGKILL)` cancellation.
- `tests/unit/test_run_status.py`: Added unit tests for interrupted status handling.
- `tests/unit/test_tools.py`: Added unit tests for timeouts and process termination.

New untracked files:
- `src/muc_one_span/inventory.py`: Frozen sample inventory record, deterministic token generation `hash(seed, design_name)`, integrity validation.
- `src/muc_one_span/durable_ledger.py`: Dual-ledger transactional manager with `.run.lock` and `.ledger.lock`.
- `tests/unit/test_inventory.py`: Unit tests for inventory creation and validation.
- `tests/unit/test_durable_ledger.py`: Unit tests for locking, atomic writing, and fail-closed loading.
- `.planning/reviews/m0-*`, `.planning/reviews/m1-*`: Adversarial review packets and disposition matrices for M0 and M1.

Line Count & Coverage:
- Authored files <= 649 lines: 100% compliant across all 135 files.
- Branch-aware unit test coverage: **90.90%** (750 passed; gate is >=80%).

---

## 2. Configuration & Inventory Hashes
- **Designs Specification:** `examples/experiment_500_designs.json` (250 designs: 150 dev, 50 val, 50 test; 7 categories).
- **Primary Simulator Configuration:** `tests/results/production_validation_20260914/fresh_generation/prepared_v2/generation_config.json` (SHA256: `445f65166c34025d73e004f2b84dbf29bfcbfb853606ce8431ee179362db4714`).
- **Primary Evaluated Modes:** `amplicon_hifi` (250 datasets) + `genomic_ont` (250 datasets) = 500 primary datasets.
- **Supplemental Mode:** `genomic_pacbio` (tracked outside the 500 primary denominator).

---

## 3. Milestone Completion Status
| Milestone | Description | Status | Verification / Artifacts |
|---|---|---|---|
| **M0** | Scope baseline, host audit, and design packet | **COMPLETE** | `.planning/2026-09-15-m0-design-spec.md`, `.planning/reviews/m0-review.md`, `.planning/reviews/m0-dispositions.md` |
| **M1** | Concurrency, durable ledgers, process registry, blinded isolation | **COMPLETE** | `.planning/reviews/m1-review.md`, `.planning/reviews/m1-dispositions.md`, `tests/data/experiment_500/ledger_{sealed,public}.jsonl` |
| **M2** | Scientific deconvolution & failure atlas | **IN PROGRESS** | Initial 30-dataset failure atlas produced (`tests/data/experiment_500/eval_benchmark_30/failure_atlas.md`); root causes isolated. |
| **M3** | Clinical report redesign & Playwright suite | **PENDING** | Designs prepared; awaiting M2 deconvolution completion. |
| **M4** | Frozen validation & sealed holdout unsealing | **PENDING** | Awaiting M2/M3 completion. |
| **M5** | Final review, documentation & delivery | **PENDING** | Awaiting M4 completion. |

---

## 4. Key Scientific Findings & Failure Atlas Root Causes (DEV Benchmark)
From initial DEV benchmark evaluation (`eval_benchmark_30`):
1. **Clair3 Homopolymer Indel Blind Spot in Mixed 50/50 Clusters:**
   In equal-length alleles (e.g. C4 `50_50`, `60_60`), reads from both haplotypes mix into a single cluster. The heterozygous `dupC` insertion (AF ~ 0.50) is predicted as `RefCall` (`0/0`, QUAL ~21) by Clair3's HiFi neural network pileup model, suppressing the mutation. Conversely, in distinct-length alleles (e.g. C7 `60_80`), length partitioning separates the alleles, Clair3 sees AF 1.0, and calls `dupC` with 100% exact sequence concordance.
2. **Unphased Genotype Collapse:**
   Clair3 outputs unphased `0/1` by default because WhatsHap output phasing was omitted (`--enable_phasing`), causing `phasing.py` to trigger fallback to a single unphased IUPAC candidate (`haplotypes: ["I"]`), leaving `allele_2` unresolved.
3. **C3 Near-Equal Length Peak Collapse:**
   Alleles differing by only 1 repeat (40/41, 50/51, 60/61, 70/71) merge into a single slightly larger length peak during kernel density / histogram peak calling.
4. **C2 Asymmetric ONT Loss:**
   Extreme length alleles (e.g. 140 repeats) suffer coverage dropout and indel smear in ONT NanoSim, failing length detection.

---

## 5. Next Immediate Actions
1. Implement Option Beta scientific fixes:
   - Introduce read-based k-mer variant separation / WhatsHap phasing before consensus calling for equal-length and near-equal-length amplicons.
   - Add `--enable_phasing` and tuned indel parameters to Clair3 invocation in `src/muc_one_span/calling.py`.
   - Implement adaptive bandwidth clustering for ONT read length distributions in `src/muc_one_span/length_candidates.py`.
   - Implement HGVS 20.05 transcript coordinate validator and fallback in `src/muc_one_span/nomenclature.py`.
2. Generate full DEV split (300 datasets: 150 designs × 2 primary platforms) via `scripts/run_500_experiment.py --split dev`.
3. Run evaluation pipeline and generate before/after DEV ablation metrics report.
4. Dispatch Milestone M2 review packet to Claude Fable 5.1 Red Team.
