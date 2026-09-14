# Asymmetric 25/140 failure: bounded causal investigation

2026-09-14, exposed cached development inputs only. Investigation stopped after
the user changed scope to finish the current work and plan future experiments.
No production change, candidate promotion, new simulator seed, caller run or
length-model expansion was performed. The existing expected failure remains.

Artifacts are ignored `tests/results/production_validation_20260914/length_failure_fix/asymmetric_causal.py`
and `asymmetric_causal.json`. The script reads cached BAMs through the existing
tool abstraction, using samtools view only. Identity matches use exact sequence
and quality in a canonical orientation, never QNAME. Original biological source
truth remains unavailable. All174 input sequence/quality identities in this case
are distinct, and all174 cached mapping primaries matched one original record.

The old pipeline calls25/28. Candidate25 has734 alignment records but168 remapped
primary records:166 with a span rounding to25, one rounding to24, and one without
a complete anchor pair. Candidate28 has278 alignment records but exactly one
remapped primary, original zero-based record173, QNAME S/136/ccs, length1531bp,
exact anchor span1501bp (approximately25 repeats). Its remapping against contig19
contains a178bp deletion. This is short-read evidence assigned to a false longer
candidate, not evidence for the140-repeat allele.

Five long original records (zero-based0,2,3,4,5) are entirely absent from both
calling groups. Four have exact spans8418,8400,8423,8416bp at every tested cached
anchor allowance0/1/2, each rounding to140. The fifth read is8439bp long and lacks
a complete accepted anchor pair. Primary mappings fall on contigs129–131;
contigs128–134 receive at most five total alignments each and are excluded by
minimum coverage10 before the false short component is selected. Removing
secondary alignments alone is not a repair: the correct reference fit can itself
be secondary, and the long component still has fewer than ten observations.

The sequence/quality-identical unique-QNAME control repeats25/28 and the same
source-record partition (minor count changes734→733 and278→279 are alignment
observations). Renaming does not repair this failure.

A credible future correction needs joint truth-independent candidate selection
and original-record assignment: preserve the coherent rare long-length component
while avoiding false modes from repeat-length noise, then establish whether its
coverage permits reconstruction. A hard minimum floor, generic ±1 suppression,
outlier deletion or lower global coverage threshold is not justified here;
previous experiments either lose the minority or create false positives, and
adjacent60/61 alleles remain a counterexample to fixed suppression. Partial
records must retain explicit uncertainty. No model was promoted, no count/event
accuracy improvement is claimed, and no expected-failure assertion was removed.
