Act as the Fable 5.1 Red Team. Review the supplied systems implementation for parallel scheduler, ledger durability, and resume behavior. Return evidence, counterexamples, and required tests. Do not modify files or execute commands.

# Review Packet: Milestone M1 — Parallel Scheduler, Ledger Durability, and Resume Behavior

Review the systems implementation described and code in tree:
1. `src/muc_one_span/tools.py`: process group creation (`start_new_session=True`), timeouts, `killpg(os.getpgid(proc.pid), signal.SIGKILL)`, pipe cleanup.
2. `src/muc_one_span/run_status.py`: distinct `interrupted` status on `KeyboardInterrupt`.
3. `src/muc_one_span/inventory.py`: frozen inventory of 500 primary samples (dev=300, val=100, test=100) + supplementary panel, split-leakage prevention, deterministic tokens.
4. `src/muc_one_span/durable_ledger.py`: authoritative ledger, file lock (`.ledger.lock`), temporary file flush + `os.fsync`, atomic `replace`, directory fsync (`fsync_dir`), verified artifact checks (`compute_sha256`), dual-ledger atomic export.
5. `scripts/run_500_experiment.py`: single truth generation per design, parallel read simulations with worker thread caps (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`), serialized authoritative ledger commits.
6. `scripts/run_eval_pipeline.py`: manifest split support (`dev`, `val`, `test`), denominator accounting against frozen inventory, worker exception isolation (preventing batch aborts), incremental progress tracking with `eval_progress.jsonl`.
7. Empirical verification: pilot generation across all three modes (`amplicon_hifi`, `genomic_ont`, `genomic_pacbio`) completed successfully.

Critique for systems failure modes: race conditions, deadlocks, orphaned child processes, partial writes, filesystem corruption, resume divergence, memory oversubscription.
