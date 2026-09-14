# Release scope, 0.11.0

The user explicitly changed the task on 2026-09-14: finish with the current work,
document experiments, add a tag/release, and support designing additional MucOneUp
simulations including ONT. This authorizes release preparation, versioning,
commit/tag publication and the repository release workflow. It does not authorize
merging the isolated branch into main; no merge is planned. Earlier prohibition
on release is superseded by this explicit instruction.

Close the ambitious original workflow as superseded with unmet endpoints, not as
scientifically complete. The 32 reserved simulation seeds and128-input/272-run final
panel were prepared but never executed. They remain available for a future frozen
experiment; no fresh validation or broad accuracy improvement is claimed for this
release. Preserve exact-scoring and tested correctness/configuration improvements.
The two strict expected length failures remain; do not remove or relax their tests.

Before publication: reproduce and resolve actual Claude configuration findings,
finish a strict configurable simulation-design runner and documentation, execute
small development-only HiFi/ONT smoke cases, run repository gates and final actual
Claude Fable 5.1 review against the complete source/diff/report including untracked
files. Document unresolved scientific methods, evidence denominators, configuration
semantics and all actual checks. Archive only closed/superseded plans with their
unchecked requirements explicitly retained. Keep simulation/results/model/index
artifacts untracked. Release notes must disclose the limits above.

Use version 0.11.0 for additive configuration/evaluation interfaces and documented
scientific output corrections. Version remains sourced from version.py; verify the locked dynamic project metadata and built wheel/resources
(the root lock entry does not embed a version, so uv.lock needs no change). Publish annotated
v0.11.0 on the isolated reviewed commit without moving main. The existing tag
workflow creates a GitHub release and may run the repository's Docker publication;
retain clear curated release notes and verify the resulting tag/release revision.
