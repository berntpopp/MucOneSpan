Review complete. Findings below are grounded in the current tree on `feat/scientific-reliability-reporting`; each is marked Confirmed (verified in code or artifacts) or Speculative (plausible, needs measurement). Method note: everything was read-only. One survey subagent ran a few read-only directory listings and greps before its shell access was revoked; no files were modified or created. The empty `.planning/reviews/m0-review.md` was left untouched.

**Summary of findings**

| ID | Sev | Area | Claim |
|---|---|---|---|
| G1 | High | Genomics/UX | BENIGN as specified is either unreachable or unsafe: the pipeline cannot establish two alleles from one observed length, and amplicon dropout of long alleles is documented |
| G2 | High | Genomics | Emitted HGVS c. coordinates are repeat-unit offsets, not transcript coordinates; several nomenclature paths are inconsistent or empty |
| G3 | High | Genomics | Post-hoc haploid genotype override silently deletes variants with AF below 0.5 and lowers the QUAL gate |
| S1 | High | Systems | No timeouts, no process-group kill, Ctrl-C durably recorded as execution failure |
| S2 | High | Systems | Evaluation runner has no durable state; one worker exception discards all results; resume trusts a bare file-exists check |
| U1 | High | UX | Report hides the multiplicity caveat and renders a fabricated 50/50 allele balance for single-cluster results |
| U2 | High | UX | Proposed CSP is incompatible with the eval-based IGV loader and cannot be a header for file:// reports |
| G4 | Med | Genomics | "Spanning reads" and "usable depth" are not implemented; counts include secondary alignments |
| G5 | Med | Genomics | Inventory and metric design: C6 contains no SNVs, genomic C7 collapses to 10x, no confidence intervals, test split underpowered per category |
| G6 | Med | Genomics | Truth sealing is nominal only; test split not selectable, thresholds already tuned on structurally identical designs |
| S3 | Med | Systems | Resource budget undercounts nested allele workers, samtools sort memory, and Clair3 TensorFlow threads |
| S4 | Med | Systems | Hardcoded operator paths in the runner and persisted into artifacts |
| S5 | Med | Systems | Atomic-commit plan lacks directory fsync, locking, and cross-file transaction semantics |
| U3 | Med | UX | Raw third-party HTML injection and BAM read names reach the DOM; JS-context escaping is by accident |
| U4 | Med | UX | Print loses collapsed evidence, all tooltip explanations, and dark-mode colours survive to paper |
| U5 | Med | UX | Current palette and interactions fail the stated AAA target; full AAA conformance is not achievable for this content |
| G7 | Low | Genomics | bcftools consensus stderr discarded, skipped overlapping variants invisible |
| S6 | Low | Systems | Provenance gaps: no git SHA, Clair3 model unhashed, simulator config path missing |
| U6 | Low | UX | Acceptance tests as written can pass with an empty viewer or a broken wheel |

## Genomics

**G1. High, Confirmed. The BENIGN state contradicts the pipeline's own evidence model and the documented dropout behaviour.**
- Evidence: single-cluster and equal-length results are hard-coded `homozygous: False`, `allele_multiplicity_status: "unresolved"`, with `allele_2` a copy of `allele_1` at `src/muc_one_span/alleles.py:604-634`; `src/muc_one_span/calling.py:278-281` reasserts this. `docs/reference/limitations.md:15` states the longer allele of 25/140 pairs yields 0 to 1 usable reads at typical template depths. `docs/reference/limitations.md:21-25` records Clair3 missing dupC in alleles over 100 repeats, with named expected-failure samples. Per-allele "usable depth" does not exist: counts are alignment records including secondaries, and `molecule_count` is null by design.
- Counterexample: amplicon_hifi, lengths 25/140, dupC on the 140 allele. One length is observed with hundreds of records, reconstruction of that allele is complete and clean. Under the spec's rule the header is either BENIGN (false reassurance) or, if "full diploid reconstruction" is enforced strictly, every C1 homozygous design becomes INCONCLUSIVE and the negative call is never usable. The spec does not say which.
- Required tests: run the existing strict expected-failure fixtures for 25/140 and `sample_dupc_100_120` through the decision table and assert the header is never BENIGN. Pre-register the expected verdict per category for dev, including the target INCONCLUSIVE rate for C1 and C7. Add a coverage-imbalance panel with the mutant allele at 20 percent, 10 percent, and 5 percent of reads and report the false-negative rate as a function of imbalance.
- Remedy: replace BENIGN with "No pathogenic variant detected within tested scope" and state the scope. Gate the negative state on two independently reconstructed alleles, per-allele primary spanning depth, allele length inside the validated sensitivity envelope per platform, and zero reference-filled spans. Single observed length in amplicon mode must carry a "possible allele dropout" caveat and cannot be a confident negative.

**G2. High, Confirmed. HGVS output is syntactically valid and numerically wrong.**
- Evidence: `format_hgvs_cdna` at `src/muc_one_span/nomenclature.py:228-243` prefixes the transcript onto a repeat-unit position, and the unit test at `tests/unit/test_nomenclature.py:121` asserts the string `NM_001204286.1:c.59dupC`. Position 59 of that transcript is in the signal-peptide coding region, not the VNTR, so a validator would accept the string as a different variant. The design's own example, c.530dupC, is a third coordinate for the same event. `g.` is never computed and the genomic constant at line 24 is unused. Novel indels reach `enrich_mutation_record` without a name (`src/muc_one_span/classify.py:400-410`) and exit with `hgvs_cdna` equal to `NM_001204286.1:c.` and `event_type` "substitution" (`nomenclature.py:330, 382-402`). All parents are named against the X unit regardless of `closest_type` (`nomenclature.py:405-410`). The `delinsAT` change is applied as delete 1, insert 2 in `src/muc_one_span/config.py:71-76`, named as `54_56delinsAT` implying delete 3, and documented as "delete 2, insert 2" in `docs/reference/mutations.md:18`. VCF indels are left-aligned by `bcftools norm` (`src/muc_one_span/vcf.py:50-63`) while names are 3'-rolled; nothing reconciles the two. The changelog claims conformance to the HGVS "stable" site while the mission requires archived version 20.05.
- Counterexample: dupC in repeat 5 and dupC in repeat 95 produce byte-identical HGVS strings.
- Required tests: a property test that HGVS output differs across repeat indices or is explicitly marked unresolved. A golden test for every catalogued mutation that applies the emitted description to an offline copy of the transcript sequence and checks the local context is the expected repeat. One authoritative `delinsAT` definition with a single test used by config, nomenclature, and docs. A test for the no-name path asserting no empty description is ever emitted.
- Remedy: make repeat-relative notation the primary clinical field, emit `c.` only when derived from the transcript sequence with a validator, and print "transcript coordinate not resolved" otherwise. Record the HGVS version string in the report.

**G3. High, Confirmed. The haploid majority override can delete a real pathogenic call.**
- Evidence: `src/muc_one_span/vcf.py:103-141` rewrites GT from AF and sets `0/0` when the maximum AF is below 0.5; lines 83-87 lower the effective QUAL threshold to 4.0 when haploid mode is on; `src/muc_one_span/calling.py:431-440` hard-codes `haploid_majority=True` for the per-allele path. Clair3 is never run with its haploid flags. Clusters are formed by gap of 5 units or valley separation of 3 units (`src/muc_one_span/settings.py:83-85`), so near-equal alleles are one cluster.
- Counterexample: C3 design 45/46 with dupC on one allele, genomic_ont. Both alleles land in one cluster; dupC AF is about 0.5 before ONT homopolymer deletion error at the 7C tract lowers it further; GT becomes `0/0`; consensus is wild type; classification is negative with a clean-looking reconstruction.
- Required tests: unit test with a VCF row AF 0.45, QUAL 30, asserting the variant is retained and flagged as mixed rather than erased. Integration test over all C3 designs with mutations on both platforms asserting the outcome is INCONCLUSIVE or correct, never a confident negative. Report the count of overridden genotypes per allele in `summary.json`.
- Remedy: treat AF between roughly 0.2 and 0.8 on a supposedly haploid cluster as contamination evidence that forces INCONCLUSIVE. Do not lower QUAL silently. Evaluate Clair3's native haploid modes on dev instead of string-rewriting genotypes.

**G4. Medium, Confirmed. Depth and "spanning" are not what the documents say.**
- Evidence: idxstats counts include secondary alignments (`src/muc_one_span/alleles.py:62`); read-dominance excludes only supplementary records (`src/muc_one_span/read_dominance.py:60`); `min_dominance_ratio` is 0.01 (`settings.py:92`), so 3 or 4 dominant records at one percent establish a second candidate. A tracked artifact shows 799 records against 132 primary. No flank-anchoring check exists anywhere.
- Required tests: synthetic BAM with N primary spanning reads plus 5N secondaries, assert the reported depth equals N. Synthetic case with 3 dominant reads at one percent, assert no "two selected candidates" without a low-support warning.
- Remedy: define usable depth as primary, non-supplementary, both-flank-anchored records and use it in the decision table; re-derive the dominance ratio from dev evidence.

**G5. Medium, Confirmed with one Speculative element. Inventory and metric design will mislead.**
- Evidence: C6 simulates only `del18_31`, `ins16bp`, `ins25bp` (`scripts/build_500_design.py:364-402`) and the truth adapter rejects SNPs (`src/muc_one_span/evaluation/truth.py:157-158`), so "SNV" in the category name is false. Genomic modes pass `max(10, templates // 5)` as coverage (`scripts/run_500_experiment.py:236, 253`), so every C7 genomic design gets 10x or 12x and standard designs get 40x. No confidence intervals exist (`src/muc_one_span/evaluation/scoring.py:41-46`). Failed normals leave the specificity denominator (`scoring.py:196-201, 345-348`) while failed positives stay in the sensitivity denominator, so execution failures inflate one and depress the other. The test split has 50 designs, roughly 7 per category. Speculative: if MucOneUp genomic coverage is total diploid depth, per-allele depth at 40x is Binomial and about half of alleles fall below a 20x per-allele gate.
- Required tests: verify MucOneUp coverage semantics on one design and record it. Add Wilson intervals to every ratio. A power statement per category before any test-split run. Add a dev-only SNV panel or rename C6.
- Remedy: pre-register acceptance as intervals, report INCONCLUSIVE rate jointly with sensitivity so a system cannot buy sensitivity with no-calls, and stratify by achieved primary depth rather than requested templates.

**G6. Medium, Confirmed as a process gap. Sealing is a filename.**
- Evidence: the runner reads the sealed ledger with truth mapping in the same process that runs inference (`scripts/run_eval_pipeline.py:39, 332`), the split flag accepts only `dev`, `final`, `pilot` (line 56) so `val` and `test` are unreachable and `final` selects zero rows and exits 0. Read-dominance margins and support floors were calibrated on the earlier 200-design dev set, which used the same length pairs the 500-design categories reuse.
- Required tests: a script that counts test designs whose (lengths, mutation, target) signature also appears in dev and fails above a stated threshold. A runner test asserting `--split test` refuses to run without an unseal record containing the frozen config hash and timestamp.
- Remedy: keep test truth in a separate directory the runner cannot read by default, write an unseal event into the ledger, and freeze the candidate config hash before unsealing.

**G7. Low, Confirmed.** `bcftools consensus` stderr is discarded on exit 0 (`src/muc_one_span/consensus.py:48-64`, `src/muc_one_span/tools.py:86-92`), so "overlaps with another variant, skipping" produces a valid-looking FASTA. Test with two overlapping indels in one VCF and assert a warning or INCONCLUSIVE. Parse the applied-variant count and compare with the record count.

**Architecture note.** Adopting Alpha for M1 while shipping the decision header in M3 means the header will be issued in exactly the regimes Alpha is documented to mishandle, since M2 validates only on dev. Route near-equal, severe-asymmetry, and ONT-smear regimes to INCONCLUSIVE by policy until M2 evidence exists, and test that policy. For Beta, the gate predicate is itself a classifier and needs its own sensitivity and specificity report; a missed gate is a silent fallback to Alpha's failure mode.

## Systems

**S1. High, Confirmed. Hangs and interrupts are not survivable.**
- Evidence: `run_tool` has no timeout and no `start_new_session` (`src/muc_one_span/tools.py:75-82`); `run_tool_iter` waits on the child in `finally` without closing stdout or killing it (`tools.py:142-147`), so a consumer that stops early deadlocks on a full pipe. `record_run_status` catches `BaseException` and writes `execution_failed` on Ctrl-C (`src/muc_one_span/run_status.py:42-44`), which the evaluator then treats as a scientific failure. A hung Clair3 stalls the runner forever, and only the direct child would ever be killed.
- Required tests: inject a sleeping fake tool and assert a per-stage timeout fires and the process group is gone. Send SIGINT mid-sample and assert the status is "interrupted" and resume reprocesses it.
- Remedy: per-stage timeouts, new session plus `killpg`, a distinct interrupted status, and generator cleanup that terminates the child.

**S2. High, Confirmed. The runner loses everything on one error and trusts stale artifacts.**
- Evidence: results are written once at the end with plain `write_text` (`scripts/run_eval_pipeline.py:380-385`); `load_truth` is uncaught at line 302 and any ledger `KeyError` propagates through `fut.result()` at line 364, aborting after all subprocess work is done; `--skip-existing` checks only that `summary.json` exists (line 209) and returns 0, so a truncated file from a killed run or a file produced under `--src-dir` baseline code is scored as the candidate.
- Required tests: crash injection at sample k, assert k minus 1 durable rows and a resume that skips exactly those. Truncate a `summary.json` and assert `--skip-existing` reruns or marks invalid artifacts. Run baseline with one source dir, then candidate with another and `--skip-existing`, assert refusal.
- Remedy: append-only JSONL with per-row fsync, each row carrying pipeline version, config hash, input hash, and an artifact hash; resume keyed on that tuple; worker exceptions captured as rows, never raised into the collector.

**S3. Medium, Confirmed with Speculative magnitude. The 8 by 4 budget is not what runs.**
- Evidence: each sample spawns a nested two-worker pool (`src/muc_one_span/calling.py:450`), remapping uses the full thread count per allele (`calling.py:403-411`), `samtools sort` runs without `-m` (default 768 MiB per thread), and no OMP or TensorFlow thread variable is set anywhere in code. The design's plan to set `OMP_NUM_THREADS=1` on all children would also throttle Clair3 inference. The mission prompt assumes 128 GB RAM; the machine has 45 GiB available.
- Counterexample: 8 samples, 2 alleles each, 4 sort threads each at 768 MiB is about 48 GiB of sort buffers alone, before Clair3 model memory.
- Required tests: pilot at full budget with RSS and load sampling, assert peak RSS under 80 percent of available and no OOM kills; check Clair3 exit codes and joblog retries.
- Remedy: a semaphore around the Clair3 stage, `samtools sort -m` set explicitly, per-stage thread caps that account for the nested pool, and TensorFlow intra-op and inter-op thread settings equal to the per-allele budget.

**S4. Medium, Confirmed. Portability and path leakage.**
- Evidence: developer home paths in `scripts/run_eval_pipeline.py:28-31, 241`; default ledger points at a nonexistent, gitignored directory (lines 39, 45). Absolute paths are persisted into `summary.json`, `run_configuration.json`, and consensus context (`src/muc_one_span/cli_settings.py:112-124`, `src/muc_one_span/consensus.py:285-292`), visible in the tracked artifact under `tests/data/experiment_500/eval_benchmark_30`.
- Required tests: a container run with no such prefix; a test that no artifact meant for sharing or any HTML contains `/home/` or `/Users/`.
- Remedy: config-driven tool roots, run-relative paths plus hashes in artifacts.

**S5. Medium, Confirmed. The atomic-commit design is incomplete.**
- Evidence: the only atomic writer lacks fsync (`src/muc_one_span/experiments.py:178-181`); no `fcntl` or lock file exists anywhere; the plan replaces sealed and public ledgers independently, which the mission itself notes is not a transaction.
- Required tests: kill -9 between the two renames and assert the consistency check detects divergence; start a second runner concurrently and assert it fails fast with a lock error.
- Remedy: one authoritative JSONL, exports regenerated from it, `flock` on the ledger, and directory fsync after rename.

**S6. Low, Confirmed.** No git SHA is recorded at runtime; the Clair3 model directory is never hashed and is logged as "external default unverified" when unset; `failure_atlas` has no provenance; the simulator config default falls back to a sibling repository path (`scripts/run_500_experiment.py:27-38`) and the primary path does not exist; the NanoSim model name appears only in spec prose. Remedy: git SHA plus dirty flag, model hash, simulator versions, and config hash in every ledger row.

## UX and reporting

**U1. High, Confirmed. The report already hides the one caveat that matters.**
- Evidence: the template never reads `candidate_duplicate_of`, `allele_multiplicity_status`, `observed_length_candidates`, `reconstruction_status`, `vcf_projection`, or `ambiguous_bases`; single-cluster results render two identical allele cards, and `total_reads` sums the duplicated allele (`src/muc_one_span/templates/report.html.j2:33`) so the balance bar shows a fabricated 50/50. The only pathogenicity wording is a hover tooltip on the Tier A badge (lines 207-213). Tier assignment uses a synthetic `support_reads` proxy, not read counts (`src/muc_one_span/nomenclature.py:333`), and a "Tier A" label invites confusion with AMP/ASCO/CAP tiers.
- Counterexample: an allele-dropout sample renders as a balanced biallelic normal.
- Required tests: golden render of a single-cluster summary asserting visible text "one length observed; second allele not established" in the header region and no 50/50 bar. Golden render of a dropout case asserting the header is not a negative state.
- Remedy: decision table driven by explicit evidence fields and covered by unit tests; a QC block in the header with per-allele primary spanning depth, multiplicity status, reference-filled spans, ambiguous bases, and projection status; rename tiers to internal evidence levels.

**U2. High, Confirmed. The offline security model as written cannot work.**
- Evidence: embedded mode decompresses igv.js and runs it through indirect `eval` (`report.html.j2:462-471`), which requires `'unsafe-eval'`; the spec's CSP omits it, so IGV would silently degrade to the "library unavailable" message and a console-error test could still pass with an empty viewer. A report opened from disk has no HTTP headers, so only a meta CSP applies and it cannot carry frame-ancestors or sandbox. No CSP exists today. `'self'` on file:// does not cover the data and blob loads igv.js performs.
- Required tests: Playwright opens the packaged report via file:// and via a loopback server, counts `securitypolicyviolation` events and asserts zero, asserts the IGV browser has at least one loaded track with rendered reads, asserts zero attempted requests using route-abort counters plus context offline mode, and asserts no page errors or unhandled rejections.
- Remedy: either inline igv.js as a plain script with a hash in `script-src`, or keep the compressed loader and document `'unsafe-eval'` as an accepted risk; enumerate `connect-src`, `img-src`, and `style-src` explicitly; avoid `'unsafe-inline'` for scripts.

**U3. Medium, Confirmed. The injection surface is narrower than feared but real.**
- Evidence: Jinja autoescape is on (`src/muc_one_span/report.py:86-90`) and the JSON literal escaper handles `<` and, via ASCII-only output, U+2028 (`src/muc_one_span/report_assets.py:148-160`). But `igv_content` is injected with `| safe` after positional slicing of third-party HTML (`report_assets.py:163-184`), igv.js renders read names and tags in popups, and the base64 payload is HTML-escaped inside a JS string context. The CSS include is autoescaped, so any future `<` or `&` in CSS silently corrupts it.
- Counterexample: a FASTQ whose read name is an image tag with an error handler, opened in a report whose CSP allows inline scripts, executes on popup.
- Required tests: adversarial fixture with hostile read names and a hostile sample file name, render, click a read in Playwright, assert no dialog and no injected global; assert the extracted fragment is a balanced DOM subtree.
- Remedy: sanitize read names before IGV export or verify igv.js popup escaping for the pinned version, hash-based `script-src`, and parser-based extraction with a balance assertion.

**U4. Medium, Confirmed. Printing drops evidence.**
- Evidence: the print block only avoids breaks inside open details (`src/muc_one_span/templates/report.css:317-325`); the per-repeat table is a closed `<details>` (`report.html.j2:309`) and print CSS cannot open it; tooltips stay at opacity zero, so every explanation including the Tier A citation vanishes; dark-mode variables are not reset; no `@page` rule and no per-page sample identifier; the IGV canvas prints blank.
- Required tests: Playwright `page.pdf()` under print emulation and dark colour-scheme emulation; parse the PDF and assert the per-repeat table text, the citation text, and the sample ID on every page; sample background pixels and assert white.
- Remedy: render evidence expanded by default and collapse only under screen media, force the light palette in print, add page margins and running headers, and replace the canvas with a print-only text alternative.

**U5. Medium, Confirmed. AAA is neither met nor fully attainable.**
- Evidence: muted text on card background is about 6.5:1 and the warning badge about 5.6:1 (`report.css:7, 110-119`), and the palette comment says "WCAG AA"; tooltips are hover-only and non-focusable, failing keyboard and content-on-hover criteria; no focus styles, no main landmark, no skip link; `tabular-nums` is absent. Beyond contrast, AAA includes reading-level and visual-presentation criteria a clinical report cannot meet, so a blanket "AAA target" is a misleading claim.
- Required tests: axe-core in Playwright with AAA tags, computed-style contrast across light, dark, and print, and a keyboard-only traversal that reaches every tooltip's content.
- Remedy: scope the claim to enhanced contrast plus AA elsewhere, make tooltips focusable and toggleable, add a text alternative for IGV, and record which criteria were checked by automation, by inspection, and manually.

**U6. Low, Confirmed. Acceptance criteria can pass on broken output.**
- Evidence: the "self-contained" test only runs with IGV off (`tests/unit/test_report.py:85-92`); the wheel check does not require the igv asset or the IGV template (`scripts/check_distribution.py:48-53`); `create_report` is not in the preflight tool list (`src/muc_one_span/pipeline.py:86`) so a missing binary fails after all science stages complete; the 80 percent coverage gate sits 11 points below the 91.34 percent baseline.
- Remedy: extend the self-contained test to the embedded report, require the asset in the wheel, preflight `create_report` when IGV is requested, count console errors together with page errors and CSP violations, confine horizontal scroll to the IGV region and assert no document-level overflow at 320 px, and set the coverage gate at baseline minus one point.

Two corrections to the packet itself: the environment section lists a 128 GB assumption in the mission prompt but the spec correctly measured 45 GiB available, so budgets must use the measured value; and the category counts in the design file do match the spec's 25/45/40/50/35/35/20 split, contrary to an intermediate survey result I discarded after checking design 167.
