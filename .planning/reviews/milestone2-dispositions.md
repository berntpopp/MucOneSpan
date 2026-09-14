# Actual Claude Fable 5.1 milestone 2 dispositions

The fresh CLI session completed successfully; `milestone2-response.json` records
`is_error=false` and the requested `claude-fable-5-1` reviewer in modelUsage.
The CLI also recorded a small auxiliary Haiku operation; Fable produced the
review. It reviewed immutable `reviews/milestone2b`, including untracked code.
The source has since received the fixes below, all before reserved seed use.

| Finding | Reproduction and disposition | Validation / remaining limit |
| --- | --- | --- |
| F1 high, read-support acceptance | Accepted scientific limitation. Existing 12:1 microvariation split and extra event are reproduced. Read phase remains library-only, opt-in and default disabled. Add per-phase-set selected-read counts and explicit internal-selection denominator. Reject the suggested universal minimum on current evidence: true 30:2 becomes selected14:1, and all45 independent sequence fixtures are correct before applying floors. | `read-phasing-cached-validation.md`, `milestone2-read-counts.md`. Joint all-original-read assignment and an evidence-based one-versus-two genotype model remain unmet; metadata is not a solution to that scientific requirement. Single-site GT selection uses the caller's genotype, not a new claim of read-level validation. |
| F2 high, discarded length evidence | Accepted diagnostic gap and underlying known inference limitation. Preserve excluded candidate-reference alignment counts and, when available, mapped-primary counts; retain distinct colliding QNAME records. No threshold or phase gate is justified by record counts alone. | `milestone2-length-accounting.md`; full77-input coverage5 trial improves count pairs but adds three normal false-alarm samples. Default10 retained. Low-coverage minority reconstruction remains unmet. |
| F3 medium, IUPAC indel projection | Reproduced and retained as explicit conservative unavailability, not a claim of absent variant support. Whole-replay proof rejects ambiguous genotype selection; the known consensus can still contain ALT. No supported-recall improvement is claimed from this restriction. Per-allele reason metadata is implemented. | Actual bcftools1.17 cached94-pair equivalence and heterozygous-indel fixture pass. Resolving each ambiguous selection jointly without exponential enumeration or inferring genotype from the desired event needs a separate tested projection design. Supported and raw event sensitivity remain separate endpoints, with this limitation reported. |
| F4 medium, duplicate recovery credit | The source snapshot predates the independent-accuracy fix; primary sequence/structure accuracy already credits one unproven identical observation. Reproduced remaining visibility gap and added `independent_missing_alleles`, `unproven_duplicate_alleles` and warnings. Preserve `missing_alleles` as literal output cardinality for compatibility. | `test_identical_consensus_copies_require_independent_genotype_evidence`; two literal identical candidates give literal2/2, independent1/2, independent-missing1, sample-exact0. |
| F5 medium, unknown status fail-open | Reproduced unknown/null/future status strings becoming completed. Replace deny lists with allow lists on both allele and classification metadata; absent legacy fields remain compatible. | `milestone2-f5-fixes.md`; predictions retained, status ambiguous, no normal TN. |
| F6 medium, attribution | Accepted. All three additional exact sequences are the majority short allele after removing IUPAC ambiguity; each paired minority candidate is twenty repeats shorter than its assigned truth. Neither phasing nor a support floor reconstructs the omitted long allele. | Existing sample table reports signed0/-20 errors and full alternatives retain count errors. Final report includes count-error strata and missing/recoverable counts; no new algorithmic recovery claim. |
| F7 low/medium, legacy unresolved alias | Reproduced `KeyError` for missing alias pointer in legacy-shaped two-allele input. Same-length calling now explicitly writes `candidate_duplicate_of=allele_1` when allele2 is not separately resolved. | `test_empty_variants_leave_identity_unresolved`, all24 calling units pass. |
| F8 low, fixture asymmetry | Accepted metadata-test gap; add 10:1 selected records with more original primary records and multiple phase sets. Actual45-fixture asymmetric experiment already covers true10:1/30:2 rather than only10:10. | `milestone2-read-counts.md`; fixture correctness does not validate erroneous input GT or omitted original alleles. |
| F9 low, docs/minimum depth | Correct stale development text to describe the experimental default-disabled helper. The claimed CLI `--min-dp` option does not exist in current `cli.py`; the Python API parameter is already explicitly documented as not applied in `vcf.py`. | Source inspection and CLI help; no new flag, depth behavior or implicit quality threshold. |
| F10 low, absolute random paths | Retained for exact local artifact provenance and collision-safe attempt isolation. Semantic benchmark comparison separates paths from sequences/events, records hashes, and never claims byte-identical complete pipeline directories. | Reproducibility noise is disclosed; stripping paths would impair traceability without improving scientific output. |

## Corrected review inputs

The request gave an incorrect shortcut path for extra-event debugging. The actual
existing directory is `tests/results/production_validation_20260914/read_phasing/
production_cached/extra_event_debug/`; its `report.json`, `primary.sam`, source map,
and six VCF/read-list outputs are preserved. The singleton is original-input
record14 and matches the short-read group by exact orientation-normalized sequence
and quality. Those joins prove original-record identity, not unavailable simulator
haplotype identity. The 14short/5long/1partial anchor result is diagnostic span
evidence with explicit uncertainty, not an independently known molecular label.

Old on-disk evaluation reports must not provide final specificity numbers. The
frozen evaluator will re-score explicit baseline/candidate inventories, including
every failed sample; fresh runs will carry current projection diagnostics.
Earlier development outputs that lack reason fields will be labeled legacy
reason-unavailable rather than silently rewritten as newly executed results.

Fresh scoped review of the fixes and final diff/scientific report remain required.

Coordinator follow-up while the scoped review was running: a retained valid
insufficient-evidence sidecar also overrode a newer timeout or exit2 record.
Two failing fixtures reproduced this. Nonzero codes other than the documented
coverage exit1, and explicit timeout status, now retain execution_failed and the
new invocation error. Valid typed insufficient-evidence with exit1 remains
compatible. All90 artifact/scoring/CLI evaluator tests pass. This later source
change is outside the milestone2-fixes snapshot and must be included explicitly
in the final review and final source freeze; no reserved seed has been used.
