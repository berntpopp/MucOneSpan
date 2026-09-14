# Maintained MucOneUp experiment workflow

Scope: user requested finish current work and reproducible future HiFi/ONT experiments.
No additional length-model expansion, no production inference promotion, no reserved
fresh seeds. Implemented bounded generation-only runner; calling/evaluation reuse
existing inventoried scripts.

Local evidence: `/home/bernt-popp/development/MucOneUp/.venv/bin/muconeup --version`
returned 0.44.5. `simulate --help` and `reads amplicon --help` plus local
`cli/commands/{simulate,reads}.py` and `read_simulator/{amplicon,ont_amplicon}_pipeline.py`
confirm repeated `--fixed-lengths`, repeated `--mutation-targets`, explicit seed,
`reads amplicon --platform pacbio|ont --coverage N`, and optional aligned BAM or
unaligned FASTQ output. Amplicon coverage is total pre-filter template molecules.
The generic `--track-read-source` help option is explicitly rejected by amplicon.
No invented source assignments or achieved-depth claims.

Implemented strict JSON schema1, per-case sample/platform/seed/diploid lengths/
one mutation with explicit targets/requested_templates, per-platform simulator
config paths and optional compatible truth dictionary. CLI config > platform
config > MUCONEUP_CONFIG. Relative design paths bind to design; simulator runs
with cwd equal to selected config directory. Absolute executable recommended.
No machine defaults, automatic grid or automatic caller launch. Original29cases
and their seeds/template counts moved verbatim to development JSON.

Generation writes expected inventory and manifest before tools. Every failed case
remains; unknown/malformed configs, command errors, invalid/mismatched truth,
ambiguous/missing input and empty records fail explicitly. Existing case results
are refused. Exact command vectors/results/cwd, version probe including failure,
config snapshots, file/model/input/output hashes, positive retained records and
strict truth warnings are archived. File-valued configuration artifacts are
rechecked before each tool stage. Full Python/external environment is not locked;
read source truth remains unavailable. Future maintained extension could add an
explicit environment manifest or separately validated source capture; neither is
silently claimed by this implementation.

TDD: initial absent-module red; 30 initial green; actual strict mutated truth
fixtures caught Event.name API mismatch (2red→34green); executable-relative and
failed-version provenance tests red; root review malformed section and changed
model tests red; final44unit tests pass before smoke. Ruff and configured mypy
pass for experiments.py/generate_testdata.py. Documentation and final broader
checks are owned/coordinated by root; no release approval claimed here.

Authorized real smoke, new development seed19200902 (already exposed after this
run): two cases each25/30 dupC haplotype1 repeat10,12requested templates. ToolPATH
prefixed env_pacbio, absolute local MucOneUp executable, existing normalized
`PV/fresh_generation/prepared_v2/generation_config.json`, outputs
`PV/experiment_runner_smoke_19200902` and design
`PV/experiment_runner_smoke_design_19200902.json` (PV means
`tests/results/production_validation_20260914`). Both completed with strict truth
validation: HiFi10usable FASTQ records, ONT12. Case wall2.47s/0.36s is generation
only. Seven artifacts hashed percase. FASTQ reflects no configured alignment
reference. Adapter warnings missing_seed/stale_repeat_lengths retained; requested
seed and command explicit in manifest. No pipeline accuracy or calibration claim.

Maintained guide: docs/guides/simulation-experiments.md. JSON examples are portable
and require an explicit simulator configuration. Root updates existing benchmark
navigation/text; no automatic migration of external result directories.

Root-added future-validation-experiment.json declares32cases (16perplatform),
reserved unused seeds2026091401–32. Guide explicitly says not generated release
evidence and no completed held-out result. Runner unit coverage91% branch-aware
(44tests); targeted run aggregatecoverage23% is irrelevant to globalCI (root owns
fullCI). Documentation final checks pending root.

## Final internal review dispositions

Two additional reviewer reproductions failed before fixes: a model modified by
the final reads command was incorrectly accepted as completed; uncleaned PATH
lookup could hash a virtualenv executable while run_tool executed another tool.
Frozen configuration/file/model hashes are now verified after every command as
well as before. Final-stage mutation leaves execution_failed and the failed-input
sentinel. Executable lookup uses existing _clean_path_for_externals on the same
PATH as run_tool and pins the absolute resolved program for every command and
version/hash metadata. Explicit relative executable stays bound to invocation
directory. Both adversarial fixtures pass after fixes. Final focused validation:
46experiment +15truth +12benchmark tests (73total), Ruff/format and configured
mypy pass. No simulator rerun required: authorized smoke used explicit absolute
executable and unchanged artifacts; guard changes do not alter generated reads.
