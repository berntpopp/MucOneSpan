# Fresh milestone fix-review dispositions

Actual reviewer: Claude Code CLI with explicitly requested `claude-fable-5-1`;
source snapshot and complete response are retained in milestone2-fixes artifacts.
No final scientific acceptance is inferred from review.

| Finding | Disposition and evidence |
| --- | --- |
| N1 unavailable projection shown as negative support | Reproduced, fixed; [dedicated evidence](milestone2-n1-fixes.md). Only verified absent support applies the absent weight; unavailable/ambiguous retain base dictionary fit and explicit report badges. |
| N2 stale phase evidence on unresolved alias | Reproduced with prepopulated second-allele evidence. Clear VCF/GT/source/context/read-phase fields, set nonindependent current evidence, skip alias in standalone consensus unconditionally. `test_calling.py` and `test_cli.py` regressions pass. |
| N3 zero selected reads does not mean observed haplotype | Documented in limitations. Internally selected read list can expose genotype-complement haplotype with zero reads; no universal floor introduced. Joint all-original-read model remains unmet, optional phasing default disabled. |
| N4 min_dp no-op described as active filter | Corrected calling docstring to explicit compatibility no-op; no nonexistent CLI option claimed. |
| N5 excluded primary count null/aggregate interpretation | Limitations document null as not queried (including absent BAM), zero as observed none; total mixes tails and separate candidates, per-contig detail required. |
| N6 status allow-list regression coverage | Added explicit insufficient_evidence and not_separately_resolved cases to parameterized artifact tests. |
| N7 changelog and prose | Added named diagnostics and missing-independent-allele fields, repaired limitations wording. |

Subsequent configuration changes expose optional phasing through explicit JSON
configuration as well as the library; the historical reviewer statement of
library-only access describes its immutable snapshot. Fresh configuration review
and final scientific/diff review remain required. Timeout/abnormal exits also
override stale insufficient-evidence sidecars, covered by focused artifact tests.
