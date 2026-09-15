# Milestone M1: Adversarial Review Dispositions
**Review Source:** Claude Fable 5.1 Red Team (`claude-fable-5-1` via `claude -p`)  
**Date:** 2026-09-15  
**Context:** Review of Milestone M1 Systems Infrastructure, Ledger Durability, and Concurrency  

---

| Finding ID | Severity | Area | Status | Disposition & Action Plan |
|---|---|---|---|---|
| **1** | High | Systems | **ACCEPTED** | **Process Group Registry & Cancellation:** Maintain a thread-safe registry of active process group IDs. On interrupt (SIGINT/SIGTERM) or unhandled error, immediately signal and terminate all registered child process groups (`os.killpg`), cancel pending futures, and shutdown executors with `wait=False`. |
| **2** | High | Systems | **ACCEPTED** | **Run Lock & Read-Modify-Write Safety:** (a) Acquire non-blocking exclusive run lock (`.run.lock` via `fcntl.flock`) at batch startup to reject concurrent overlapping runs. (b) Inside `commit_entry`, re-read existing entries under file lock before merging and atomically rewriting. |
| **3** | High | Systems | **ACCEPTED** | **Atomic Staging & Corruption Defense:** Simulate all reads and truth files into an isolated temporary directory/file, verify integrity (non-empty, valid record structure, quickcheck), and atomically rename into target location only on exit code 0. Eliminate in-place header rewrites. |
| **4** | High | UX/Genomics | **ACCEPTED** | **Blinded Path Integrity & Leakage Prevention:** Ensure public ledger references the exact blinded path (`blinded_reads/{token}/{platform}/{token}_{platform}.fastq`) and assert that no design names, categories, or length numbers leak into public artifacts. |
| **5** | High | Systems | **ACCEPTED** | **Verified Atomic Blinded Copies:** Copy raw reads to blinded path via temp file + fsync + atomic rename; verify SHA256 matches raw read SHA256 before committing. |
| **6** | High | Genomics | **ACCEPTED** | **Runtime Split Integrity & Test Split Guard:** Call `validate_inventory_integrity` at startup in both generator and evaluator. Fail closed on split conflicts, missing designs, or test split access without explicit `--unseal-test` flag. |
| **7** | Med | Systems | **ACCEPTED** | **Reconciliation & Fail-Closed Loading:** Fail closed if existing ledger lines are malformed. Automatically reconstruct public ledger from sealed ledger on startup if out of sync. |
| **8** | Med | Systems/Genomics | **ACCEPTED** | **Interrupted Status in Artifacts:** Update `src/muc_one_span/evaluation/artifacts.py` to accept `interrupted` as a recognized run status. Map negative exit codes (killed by signal) to `execution_failed` with the signal name. |
| **9** | Med | Systems | **ACCEPTED** | **Runner Resume Replay & Atomic Reports:** Replay `eval_progress.jsonl` on resume; verify matching read hashes before skipping. Fsync report files before directory fsync. Use standard tool abstraction. |
| **10** | Med | Genomics | **ACCEPTED** | **Platform Alias Denominator Accounting:** Map platform aliases (`hifi`, `ont`) to canonical platform modes in denominator expectations. Fail fatally if `--inventory` path is missing. |
| **11** | Med | Systems/Genomics | **ACCEPTED** | **Provenance & Biological Truth Validation:** Populate `git_sha` and `config_sha256` in all ledger entries. Validate simulated truth FASTA (2 haplotypes, matching lengths) before running read simulation. |
| **12** | Med | Systems | **ACCEPTED** | **Tool Timeouts & Samtools Memory Limits:** Add default timeouts to tool calls; set `-m 500M` on `samtools sort`. Set OpenMP, BLAS, and MKL thread caps before any simulation or evaluation stage. |
| **13** | Med | UX | **ACCEPTED** | **Order-Invariant Token Hashing:** Compute sample tokens deterministically from `hash(seed, design_name)` so token assignment is invariant to design file ordering. |
