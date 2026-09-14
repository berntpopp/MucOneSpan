Use the explicitly selected claude-fable-5-1 as an independent adversarial reviewer.
Fresh session: inspect only bounded fixes to your preceding milestone review,
with at most 2500 words of actionable final findings. Do not re-audit unrelated
modules or propose arbitrary scientific thresholds. Read-only tools only.

Immutable snapshot:
/home/bernt-popp/development/MucOneSpan/tests/results/production_validation_20260914/reviews/milestone2-fixes
manifest.json includes untracked files and tracked.patch is against d8390b3.
Read .planning/reviews/milestone2-dispositions.md and the F5, read-counts and
length-accounting notes, then exact source/test changes in evaluation/artifacts.py,
evaluation/scoring.py, alleles.py, read_phasing.py, calling.py and docs/development.md.

F1 optional read phase remains experimental Python-only read_phase=True, default
false, noCLIflag. Adds counts perPS/hap with explicit internally-selected-read
denominator; tests30primary/11selected10:1. No universal floor: actual true30:2
becomes selected14:1 and known-truth45 original genotype pairs all exact. Minimum
floors remove real minority alleles; joint original-read model remains unmet.
F2 low-coverage exclusions now have per-allele nested diagnostics and one narrow
samtools primary-record query, retaining QNAME collisions; passing groups beyond
top2 also explicit. No threshold/phase gate or scientific accuracy claim.
F4 primary independent seq/structure metric already fixed before this review;
now explicit independent_missing_alleles/unproven_duplicate_alleles+warning. Literal
missing_alleles retains cardinality, so no schema confusion with independent missing.
F5 allow-lists on both allele and classification phase/reconstruction fields keep
absent legacy fields compatible but explicit unknown/null/future values ambiguous,
never TN; predictions retained. Eighteen new failing regressions reproduced.
F7 legacy unresolved allele2 now explicitly points candidate_duplicate_of=allele_1;
existing legacy-shaped empty-variant regression failed then passes.
F9 docs corrected; --min-dp does not exist in CLI, API no-op already documented.

F3 conservative support limitation deliberately retained: a selected IUPAC
heterozygous indel prevents full replay proof and yields projection_unavailable
with explicit reason, not negative experimental evidence. No supported sensitivity
gain claimed. Biological whole-allele identity is not inferred from annotation.
Question whether this is documented unmet scope versus a new false contract; do
not implement ALT selection heuristics without proving multiallelic/version behavior.

Correct raw debug path (previous request had a wrong shortcut):
/home/bernt-popp/development/MucOneSpan/tests/results/production_validation_20260914/read_phasing/production_cached/extra_event_debug/report.json
It includes exact sequence/quality source joins and all six parameter interventions.
The singleton is original record14, from short-read group, not the long/partial
record. These are source-record matches, not unavailable original simulator labels.

All reserved final seeds remain UNUSED. Fresh generation/paired128input evaluation
will follow source/scoring freeze; no final success is claimed. Last fullCI before
these final review fixes441passed92.16%, latest focused allele39, readphase/calling65,
artifact/scoring88 pass; actualphaseintegration17pass. Root will rerun whole gates.
For actionable bugs give severity, exactlocation, reproducer/evidence and validation.
Distinguish fixed claims, accepted unmet scientific methods, and genuine introduced
bugs. Fresh final diff/scientific report review remains a separate required session.
