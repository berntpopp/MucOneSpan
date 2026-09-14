Use the explicitly selected claude-fable-5-1 as an independent adversarial reviewer.
Fresh read-only session. Review the bounded runtime configuration extension and
recent N1/N2 fixes; do not spend the review re-auditing unrelated old workstreams.
Return actionable findings with severity, exact source location, failure scenario,
evidence/reproducer and suggested validation, at most2500words.

Immutable snapshot:
/home/bernt-popp/development/MucOneSpan/tests/results/production_validation_20260914/reviews/runtime-configuration
manifest.json SHA256 dc80d83a3b06b9652bd38ec8bd28153bce21a681ff5c519df5e1f0ce04e5b200
219 files include untracked implementation/tests/docs; tracked.patch against
base d8390b3c244ef8f3240af74b92db12b50dfc77d1. Read source ONLY from snapshot;
raw development evidence can be accessed in sibling PV directories.

Requirements/spec .planning/2026-09-14-runtime-configuration-spec.md and user guide
docs/guides/configuration.md; internal audit plus dispositions in .planning/reviews.
No machine-specific production defaults, sample/truth/seed branches or mutable
global runtime settings. Strict JSON unknown/type/range/duplicate/nonfinite checks,
relative path resolution, explicitCLI > config >defaults. Fixed scientific
coordinate/scoring/genotype contracts are not user-adjustable knobs.

Focus source settings.py, cli_settings.py, cli.py, pipeline.py; config integration
in classify/classification_summary/alleles/ladder/consensus/calling/read_phasing;
evaluation fixed_repeat_count and generic count diagnostic. Check accepted fields
actually affect their claimed stages, explicit overrides remain compatible,
reference layout uses chosen9IDs not category12IDs, stale provenance cannot leak,
and pipeline extraction preserves stage execution/failure semantics.

Internal audit fixed standalone dictionary ignoring, invalidCLIstage overrides,
stale configuration after failed rerun, anchor start-vs-boundary coordinates,
short/oversized sliced anchors, actual phaser argv provenance. Standalone consensus
now deliberately uses bundled/configured dictionary for anchor-aware trimming
like fullrun (priorfixed-only behavior), documented migration. Reference leftflank
still originalprefix; changing the biological reference is separate unmet work.
Actual cached137fullconsensuses across77developmentinputs:137FASTAidentical,
only left_trim_method now exact_anchor;17missingalleles explicitlycounted. Synthetic
500bpflank+1bp insertioncorrected181→180VNTRbp, so universalparityis NOTclaimed.
Scoring/classification completeoutputparity and length/generatorparity in notes.

Previous actualFable N1 fixed: unavailable/localizationambiguoussupport retains
base dictionary-fitconfidence plus independentboundarypenalty; only verifiedabsent
applies absentweight. Report distinguishes statuses. N2 fixed: unresolvedalias clears
staleVCF/GT/source/context/readphase and standalone skipsaliasunconditionally.
See milestone2-fixes-dispositions.md and milestone2-n1-fixes.md.
Experimentalphasing now configorPythonopt-in, defaultfalse due faileddevelopment
falsecallgate. No readthreshold invented. Minimumcoverage10/QUAL5 remain frozen
planneddefaults. Threads per-tool, NOTstrictwholeprocesstreecap.

Actual current gates: make ci-check626passed93.59%coverage; lint/format/mypy/filegate
passed; make test-int47passed8deselected; make docs-checkpassed. Generated-dataE2E
and buildcheck stillfinishing, finalfreshscientific/performanceevaluationpending.
FirstfullCI caught a fixturemockingmissingreference; fixedfixturewritesactualFASTA
for newhashing, didnotweaken reportwarningassertions. Internalred6tests→green14CLI.
All reservedfinalseeds remain UNUSED; no tuning on finaldata possible yet.
Do not mistake passingtests/review for a scientificmethod promotion. Report real
bugs that need fixing before freeze; retain unmetmethods explicitly instead of
inventing thresholds. Finaldiffandscientificreportreview will be separate.
