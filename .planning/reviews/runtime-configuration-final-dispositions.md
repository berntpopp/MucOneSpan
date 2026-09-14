# Actual Claude configuration review dispositions

Reviewer: installed Claude Code2.1.268, explicit claude-fable-5-1, fresh read-only
session against snapshot manifestdc80d83a3b06b9652bd38ec8bd28153bce21a681ff5c519df5e1f0ce04e5b200.
Response retained locally as runtime-configuration-response.json; modelUsage
confirms claude-fable-5-1. Tiny auxiliary Haiku operation was not a substituted
reviewer. Six internal fixes and previousN1/N2 verified by the reviewer.

| Finding | Disposition and actual evidence |
| --- | --- |
| F1 requested flank exceeds available dictionary flank | Reproduced then rejected with shared validate_flanks before ladder/consensus tools/output. Explicit0 remainsvalid, explicitoverrideswin. Sixredtests,147focusedunit+20realtooltests green;137cachedconsensusbytes/coordinates anddefaultladderhashunchanged. See consensus-flank-validation-review-fixes.md. |
| F2 invalid run values leak bareValueError | Three actualCliRunner regressions (threads0,coverage0,QUALnan) reproducedexit1; nowBadParameter/exit2 withinrecordedrun, no map invocation. |
| F3 recorded auto preset becomes explicit onreuse | ActualroundtriptestreproducedHiFi map-hifi persistingafterplatformONT override; preserve nullautoselectioninsettings and record resolved_minimap2_preset separately. CLIrelativeexplicitmodel/referencepaths normalizedforrecording. |
| F4 global config failure leaves prior attempt files | Documented andtestedchosenpolicy: malformedglobalconfig exits2 before runcallback, outputuntouched. Callerexitrecordrequiredby evaluation; previousfilescannotprove successfornewinvocation. |
| F5 stale vntr_phase_status | Addedstalealiasfieldtoregression (failedbeforefix), clearitduringdisambiguation; consensuswithoutindependentphasewritesunresolved. No scoring/sequencechange. |
| F6 custom layout backward hint order | Acceptedlimitationdocumented: classifier backwardhintusesdictionarycategoryorder; generalmatcherfollows. Layoutcontrolsselectedfixedcounts/constructedsequence/anchors, arbitrarylayoutreconstructionnotvalidated. |

Reviewer suggestion of exact reference-length validation is deferred: arbitrary
valid repeat dictionaries can contain variable-length units; a naive count×unit
formula would falsely reject them. Currentexplicit-referencecompatibilityremains
usercontractanddocumentedlimitation; no claimofautomaticallyverifiedreferenceidentity.

User subsequently requested finishing with currentvalidatedwork and publishing
0.11.0, plusfutureMucOneUpHiFi/ONTexperimentdesign. Reservedfreshpanelisnotexecuted;
release doesnotclaim freshscientificvalidation. FinalactualClaude diff/reportreview
will include these fixes, newexperimentworkflow and packagingguard. No review
approval replaces actualchecks.
