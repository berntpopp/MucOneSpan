Static review complete. I read all six implementation files, their unit tests, the M0 spec and dispositions, and the evaluation consumer that reads the run-status sidecar. Nothing was executed, so every counterexample below is reasoned from the code, not reproduced. The verdict: the durability primitives are individually sound, but the scheduler, the resume path, and the blinding export have defects that will seal corrupt or mislabeled data into the authoritative ledger, and the pilot run exercises none of the components under review.

## Findings, ranked

**1. Interrupt and first-error abort do not stop the generation batch, and forced kills orphan the simulators. High.**
- Evidence: `scripts/run_500_experiment.py:490-505` re-raises on the first worker failure inside the executor context. Leaving the context calls shutdown with wait only, so every queued task still runs to completion, but the commit loop has already exited. `src/muc_one_span/tools.py:91` puts each simulator in its own session, so terminal SIGINT never reaches it. The kill path at `tools.py:102-107` only fires when the exception is raised in the thread blocked on `communicate`, and KeyboardInterrupt is delivered to the main thread only.
- Counterexample: Run dev with 300 datasets, press Ctrl+C after 10 commits. The process prints nothing further but keeps simulating the remaining 290 for hours, and none of those results is committed. The user then sends a second Ctrl+C or SIGKILL. Up to 8 simulators keep running under init, still writing into the read directories. A restart then finds their half-written FASTQs and seals them via finding 3.
- Required test: Integration test with a fake simulator that sleeps and writes slowly. Send SIGINT to the script, assert it exits within a bounded time, assert every child process group is gone, and assert the ledger contains exactly the entries committed before the signal.

**2. The ledger lock serializes writes but not read-modify-write, and there is no run-level lock. High.**
- Evidence: `src/muc_one_span/durable_ledger.py:73` loads the ledger once in the constructor without the lock. `commit_entry` at lines 141-146 mutates the in-memory map and rewrites the whole file from memory. Nothing prevents two generation processes from sharing an output directory.
- Counterexample: Terminal A runs `--split dev --platforms amplicon_hifi`, terminal B runs `--split dev --platforms genomic_ont` at the same time. Each commit overwrites the other process's rows. The final sealed ledger holds only the rows of whichever process committed last plus whatever it loaded at startup. The eval runner then counts the missing datasets as failed. Two overlapping selections would also drive two simulators into the same output files.
- Required test: Two subprocesses each commit disjoint entries to the same directory concurrently. Assert the union survives. Add a test that a second process on the same output directory refuses to start while a run lock is held.

**3. Resume accepts any non-empty file as a finished dataset, including partial output from a killed simulator or from the in-place header rewrite. High.**
- Evidence: `scripts/run_500_experiment.py:189-203` returns the first non-empty glob match with no exit-code marker, sentinel, or hash check. The truth check at line 117 is size-only. `_sanitize_fastq_headers` at lines 156-172 rewrites the FASTQ in place with a plain write and skips work whenever the first header already carries the suffix. `is_verified_complete` at `durable_ledger.py:137-139` then verifies the recorded hash, so a sealed corrupt file verifies as complete forever.
- Counterexamples: SIGKILL during simulation leaves a FASTQ truncated at a record boundary. Restart picks it up, counts its records, hashes it, seals it. Truncated mid-record instead: the record counter raises, the batch aborts, and every subsequent resume hits the same file and aborts again, a permanent wedge that needs manual cleanup. Kill during the header rewrite: the first header is already rewritten, so the truncated file is committed unchanged. The same size-only logic applies to a truncated truth FASTA, which then feeds every platform for that design.
- Required test: Crash-resume matrix using a fake simulator with an injectable kill point. For each of mid-simulation, mid-rewrite, and mid-truth, assert the resume path either regenerates or fails closed, never commits. Design fix to test against: simulate into a temporary directory and rename on success, or write a completion marker after exit code zero.

**4. The public ledger points at a path that does not exist and leaks the design name into the blinded artifact. High.**
- Evidence: `durable_ledger.py:170` builds the public path from the raw read file name. The worker at `run_500_experiment.py:301` names the blinded copy from token and platform. The raw name is derived from the design name, for example `c1_homo_30_30_dev`, which encodes category, both allele lengths, and split.
- Counterexample: A blinded evaluator opens the public ledger, sees `blinded_reads/sample_0183/amplicon_hifi/c1_homo_30_30_dev_amplicon_hifi.fastq`, learns the truth from the filename, and finds the file missing anyway. The test-split seal is broken by construction.
- Required test: After commit, for every public row assert the referenced file exists and that neither design name, category, nor lengths appear anywhere in the public ledger text.

**5. The blinded copy is neither atomic nor verified, and a truncated copy is kept forever. High.**
- Evidence: `run_500_experiment.py:302-303` skips the copy whenever the destination exists, and the copy itself is a direct write.
- Counterexample: Kill during the copy. The destination exists but is short. Every later run keeps it. The sealed ledger hash refers to the raw file, so nothing ever notices the blinded panel is corrupt.
- Required test: Pre-create a truncated blinded file, run the worker, assert the file is replaced and its hash equals the raw file's hash.

**6. Split-leakage prevention is never invoked at runtime, and nothing seals the test split. High.**
- Evidence: `validate_inventory_integrity` at `src/muc_one_span/inventory.py:120` has no caller outside `tests/unit/test_inventory.py`. Both scripts call `build_expected_inventory` and proceed. `--split test` is accepted by both scripts with no guard, contrary to disposition G6.
- Counterexample: Hand-edit the designs file to move one design from test to dev. Generation proceeds. The old ledger rows still say test, the new inventory says dev, and the eval runner's denominator marks the design as missing from dev while silently retaining the stale test row.
- Required test: Script-level test that a designs file with a duplicated name, 249 entries, or a per-design split conflict exits non-zero before any subprocess is launched. Add a test that the ledger's split for each design matches the inventory's split at eval time.

**7. Dual-ledger export is not atomic across the pair, and malformed rows are silently dropped on the next rewrite. Medium.**
- Evidence: Two separate replaces at `durable_ledger.py:178-179`. Commits run in the main thread, which is exactly where SIGINT lands. `_load_existing` at lines 116-117 warns and skips unparseable rows, and the next commit rewrites the file without them.
- Counterexample: Ctrl+C between the two replaces leaves the sealed ledger one row ahead of the public one, and nothing reconciles them at startup. A single corrupted row disappears permanently after the next commit.
- Required test: Inject a failure between the replaces and assert the next constructor regenerates the public ledger from the sealed one. Assert that a malformed row causes construction to fail rather than proceed.

**8. The new interrupted status is rejected by the evaluation consumer, and signal kills leave the sidecar at running. Medium.**
- Evidence: `src/muc_one_span/evaluation/artifacts.py:175-179` accepts only four statuses, so an `interrupted` sidecar is treated as an invalid sidecar with a warning. `run_status.py:42` only catches KeyboardInterrupt, so SIGTERM and OOM SIGKILL leave `running`. `run_eval_pipeline.py:135-148` then falls through to the mapping check and attributes the death to a pipeline stage.
- Counterexample: Memory pressure OOM-kills one child. Its sidecar says running, its exit code is negative, and the failure atlas records `initial_mapping` as the first failing stage.
- Required test: Sidecar with `interrupted` loads without the invalid-sidecar warning. A negative exit code with a running sidecar must diagnose as execution with the signal number.

**9. The evaluation runner bypasses the tool abstraction and has the same non-cancellable executor; its progress file is write-only. Medium.**
- Evidence: `run_eval_pipeline.py:287-289` uses `subprocess.run` with no session, no timeout, and no group kill, against the repository rule. Lines 454-464 repeat the shutdown pattern from finding 1. `eval_progress.jsonl` is appended at lines 462-471, never read, never truncated across runs, and never fsynced. The final reports at lines 515-560 are plain writes followed by a directory fsync, which does not make their contents durable.
- Counterexample: Two runs into the same output directory produce duplicate progress rows. A crash mid-batch loses all in-memory results because nothing replays the progress file, and the skip-existing path re-evaluates stale summaries without checking they came from the current reads hash.
- Required test: Replay test proving a restarted run reconstructs its results from the progress file and re-runs only missing samples. Duplicate-row test across two invocations.

**10. Denominator accounting is silently disabled for the platform aliases and when the inventory path is wrong. Medium.**
- Evidence: `run_eval_pipeline.py:477-480` compares the alias `hifi` or `ont` against inventory mode names, so every expected record is skipped. Line 415 makes the whole denominator conditional on the inventory file existing. `--limit` at line 438 truncates the ledger but not the expectations.
- Counterexample: `--split test --platform ont` with an empty ledger reports zero evaluated and zero missing. A typo in `--inventory` yields the same silent result.
- Required test: Alias run with an empty ledger must yield 50 missing rows. A missing inventory path must be fatal.

**11. Provenance and truth validation are missing from the generation path. Medium.**
- Evidence: `git_sha` and `config_sha256` on `LedgerEntry` are never populated at `run_500_experiment.py:305-321`. The default config search at lines 30-41 reaches outside the repository. Truth generation at lines 107-153 checks only file existence, whereas `experiments.py:388-401` validates haplotype count, lengths, and mutation events.
- Counterexample: A sibling checkout of the simulator config changes. Reads regenerate differently, the ledger cannot say which config produced them, and a mutation-name mismatch produces a wrong truth that seals silently.
- Required test: Committed entries must carry non-null provenance. A truth FASTA with the wrong haplotype count or lengths must fail before read simulation.

**12. No caller passes a timeout, and the memory budget from disposition S3 is not enforced. Medium.**
- Evidence: All 40 call sites of the tool runner pass no timeout. `mapping.py:121-128` runs `samtools sort` without a memory cap. The eval runner allows 8 parallel samples with 4 threads each. The OMP and BLAS caps do not govern NanoSim's multiprocessing pool or Clair3's TensorFlow thread pool. `run_500_experiment.py:289-290` sets the caps globally from worker threads, omits the MKL cap, and Phase 1 truth generation runs before any cap is set.
- Counterexample: One hung simulator holds a worker forever, and executor shutdown never returns. Eight Clair3 runs plus uncapped sort buffers can exceed the 45 GiB available.
- Required evidence: Measured peak RSS for one sample per mode and for the full parallel configuration. Timeout test with a fake tool that sleeps past the limit.

**13. Portability and tokenization deviate from the accepted dispositions. Medium.**
- Evidence: Machine paths remain at `run_eval_pipeline.py:36`, `run_eval_pipeline.py:278`, and `run_500_experiment.py:60` despite S4. Tokens at `inventory.py:52-54` are positional over file order, so reordering the designs file silently re-tokenizes everything already distributed. Tokens are per design rather than per dataset as the spec states, so a blinded evaluator can pair the two modes of one design.
- Required test: Token for a design must be a function of seed and design name only, invariant to file order.

**Low items.** The timeout unit test mocks the group kill and never verifies that a grandchild dies. The logging unit test mocks `subprocess.run`, which the runner no longer uses. The header-rewrite idempotency check depends on the first read name lacking `_r`, so a contig name containing that string disables sanitization. `--workers` is clamped to 8 without warning. Startup hashes every artifact serially.

## Empirical verification claim

The pilot at `run_500_experiment.py:324-401` constructs no ledger, uses no executor, never resumes, and never invokes the evaluation runner. It therefore provides no evidence for the scheduler, durability, or resume claims in packet items 4 through 6. The eval runner's pilot split reads the 500-dataset ledger, which has no pilot rows, so it evaluates nothing. Unverified and needing evidence from an actual run: which files each simulator mode leaves in the output directory, since the first-match glob may select an intermediate BAM or unaligned-read file; whether `samtools` is on the generation script's cleaned PATH when a BAM is produced; and the per-mode file sizes that the in-memory header rewrite must hold.

## Required test harness

Beyond the per-finding unit tests above, the milestone needs one integration harness under `tests/integration/` using the existing marker conventions and a fake simulator script, covering: SIGINT and SIGKILL delivery at each phase with assertions on child process groups and ledger state; two-process concurrency on one output directory; a crash-resume matrix with injectable kill points; a property check that the public ledger references only existing, non-leaking paths; and one minimal full run of two designs through generation, resume, and evaluation with the denominator asserted at both stages.
