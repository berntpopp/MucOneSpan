# Wave 2 frozen baseline decisions

## Before clinical prediction interpretation

Caller baseline: clean main c08000a0b97e9ee00e560a369ef320b4d06d3004 (v0.14.1). The separate harness uses the identical tracked caller code/resources; each invocation verifies source/resource/lock/tool/model hashes. Generated `.fai`/other runtime indexes are excluded from source identity because they are not tracked input resources. Initial preflight caught an index mistakenly included in source hashing; no clinical process had executed, and freeze was corrected before the smoke invocation.

The clean-main panel completed all three normal HiFi/dupC HiFi/dupC ONT fixtures and captured84 scientific artifact hashes. Fixed two threads. Panel tool models follow its established fixture settings; clinical cohort model differs as detailed below. Repeat panel after implementation with the same fixture models/settings and compare normalization of only paths/timestamps.

## Cohort and preprocessing

All20 ENA runs retained in inventory: nine primary MUC1 amplicons, two secondary MUC1 WGS-labelled FASTQs, nine excluded ACAN. All11 selected files downloaded with byte count and MD5 validation; local SHA256 retained. Total compressed size198156449 bytes, unchanged from supplied snapshot. No PacBio runs. Decompress only, preserve all IDs/counts/bases/qualities and both orientations. No length/quality/expected-mutation filtering or cropping. Short reads remain in denominators. See source and input audits for empirical primer/anchor geometry.

WGS arm is exploratory use of the existing ONT FASTQ path on the deposited reads, with no conversion to an asserted PCR assay. Broad whole-genome throughput is not measured. One patient WGS identity is unresolved between ENA and publication; execution remains eligible and truth identity is separate.

## Endpoints

- Execution completed: CLI exit0 plus validated current artifacts and completed sidecar. Actual nonzero failures retained, including report failures.
- Callable: execution completed AND existing strict observation parser says completed. Ambiguous reconstruction is separately visible, never a negative.
- Exact named-event recovery: literal full mutation identity (dupC) among callable runs, regardless positional assignment. Supported named-event recovery separately requires existing exact sequence/VCF support boolean.
- MP1–MP4 are publication-reported known dupC controls. Score their recovery under `reported_known_control`; no paper-linked assay establishes an independent-confirmation sensitivity denominator. Do not borrow Wenzel family truth without participant mapping.
- Negative specificity is not estimable: no endpoint-specific independently confirmed negative set is established. HG aliases do not create negatives.
- MP5 event and published motif positions/allele assignments are comparator-only, not independent truth.
- HG002 independent Q100v1.1 sequence benchmark: exact motif1-to-motif9 intervals on coding-forward strand, both full anchors retained. Preserve the maternal18bp biological insertion. Allow only documented anchor interval/orientation and allele order. Report literal sequence comparison and independent recovery separately; allele assignment uses every minimum-edit-distance match, with length errors and exact equality. Do not infer repeat count by length/60 for noncanonical indels.
- All-sample denominators include failures/no-calls/unattempted. Conditional callable denominator explicitly separate. Biological identities and relatedness reported separately; no binomial confidence interval or independent PCR-read interpretation.

## Execution/resource settings

Serial clinical cohort; two requested threads;1800-second total supervised lifecycle budget per library. One observed execution per run in the final cohort; smoke excluded from performance aggregates. No threshold/POA/clinical tier/read-phase/strict-segmentation changes. Report enabled, existing default IGV off; report phase timed separately. Download/preprocessing timing separate.

CPU execution with CUDA_VISIBLE_DEVICES empty. Installed minimap2 2.28-r1209, samtools1.15.1, bcftools1.17, Clair31.0.10. Study describes Dorado0.7.1 SUP and SQK-LSK114; exact original model unavailable. Selected Rerio r1041_e82_400bps_sup_v500 before clinical predictions from Dorado0.7.1 contemporary SUPv5.0.0 default documentation. Archived source pins/checksums and synthetic checkpoint-load proof accompany model provenance. No model comparison on diagnostic outcomes. The installed generic ont alias was older R9.4 HAC and is not selected for this cohort.

Resource measurement: isolated worker wall time and per-stage timing using existing benchmark wrappers; RSS=max(RUSAGE_SELF,RUSAGE_CHILDREN) on Linux inKiB. This is the largest individual-process high-water value, a lower bound on concurrent process-tree peak; no claim of exact process-tree peak. Single observations, not medians. Caller process tree is supervised by existing Linux subreaper, with bounded termination and descendant cleanup. Full exact commands, timestamps, CLI codes, sidecars, logs and hashes retained outside Git.

## Review state

Independent data review found and reproduced collision/restart/atomic-output/schema bugs; fixed with synthetic regressions and independently re-reviewed. Runner review found identity binding, malformed-worker, timeout-precedence and per-run continuation issues; corrections and regressions in progress before full cohort. Scoring review found unnamed de novo annotations must remain valid evidence; named endpoint matching excludes only unnamed identities. No caller scientific code was modified to address these harness bugs.
