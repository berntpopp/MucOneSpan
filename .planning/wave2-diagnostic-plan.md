# Wave 2 diagnostic extension

User steering (2026-09-16): debug why positive controls are missed and compare
with the study's VNTRPipeline. Preserve completed frozen cohort-v2 and original
release objective; complete diagnostic evidence before finalizing publication.

## Questions and ownership

1. Scoring investigator: prove MP1 insertion representation and MP3 genotype
   selection causes by tracing final artifacts and minimal isolated replay;
   inspect MP2/MP4 candidate losses. Own `wave2-miss-debug.md` only.
2. Source investigator: paper/source-version/defaults and algorithm comparison;
   distinguish demonstrated code behavior from hypotheses. Own
   `wave2-comparator-audit.md` only.
3. Lead: run comparator, bind inputs/source/image/settings, validate scientific
   outputs and summarize the comparison. Separate reviewer audits conclusions.

## Prespecified comparator execution

Use source tag v1.0 (f7f594e74ce9cde273deb269f06bb528f0586997), closest
publication-era tagged source; paper provides no exact commit/image. Mount its
unaltered scripts read-only into the installed image pinned by digest
`sha256:888583f8ef0b69af9c5b386c51c18a2f6a73b0d46dac06424ebd56b01845d5a8`.
The runtime image was built March 2026 and contains v2 scripts before overlay;
its Canu build is master r10516 fb43b3c. This is a source-pinned diagnostic
comparison in available pinned dependencies, not exact paper-environment replay.
Registry image tags v1.0/v2.0 are unavailable. Record actual image/tool hashes.

Start with MP1 to check wiring, then process remaining eligible inputs with the
same published source defaults: `-v MUC1 -r t2t`, `-p pcr` for nine amplicons,
`-p wgs` for two WGS-labelled inputs. `DELETE_TMP=N` preserves evidence. No manual
allele lengths, frequency changes, forced phasing or outcome-based read filters.
Source-internal preprocessing is part of the comparator, starting from the same
validated FASTQs. Container limit2CPUs/16GiB; hardcoded internal32-thread requests
remain unchanged and Canu autodetects resources. No matched speedup claim.
Per-case lifecycle budget1800seconds, serial execution; retain failed/partial cases.
Mount only the selected reads read-only and task-owned output read-write. Network
is disabled during clinical execution. Existing MucOneSpan tool supervision wraps
container execution, and explicit container cleanup handles daemon-owned processes.

## Comparison and interpretation

Record actual final FASTA availability, independently normalized HG002 sequence
concordance, candidate lengths, algorithm branches, allele read counts and LoF/QC
outputs. Broad LoF alarms are not automatically exact named dupC recovery.
Comparator sequences/annotations are comparison evidence, never independent
patient truth. Validate output artifacts separately from wrapper exit status.
Separate original frozen clinical metrics from diagnostic replay results. No
scientific production edits or threshold tuning are included in this investigation.
