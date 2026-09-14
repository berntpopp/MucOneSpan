# Runtime configuration documentation and examples

Owned documentation changes:

- `docs/guides/configuration.md`: schema, global CLI usage, explicit CLI > file >
  central-default precedence, relative paths, all parameter defaults/ranges,
  immutable Python settings, custom reference consistency, provenance and
  scientific invariants.
- `examples/runtime-settings.json`: complete canonical defaults generated with
  `json.dumps(settings_as_dict(DEFAULT_SETTINGS), indent=2)`, including selected
  reference IDs and ladder range. No manually copied default values in this file.

No production source, tests or MkDocs navigation were edited for this task.
The coordinator owns navigation and final CLI integration verification.

## Examples and verification

The guide's minimal schema-one configuration loads as DEFAULT_SETTINGS. Its
Python example loads a file and constructs/derives immutable settings with eight
threads. Executed both examples in a temporary directory; all assertions passed.
Loaded the repository's full JSON example with load_settings and compared both
the typed result and serialized object to DEFAULT_SETTINGS; both matched.

The guide's shell commands specify the global option before the command:

```bash
muconespan --config settings.json run --input reads.bam --output-dir results --threads 8
muconespan --config settings.json ladder --output matching_ladder.fa
muconespan --config settings.json run --input reads.bam --reference matching_ladder.fa --output-dir results
```

These shell examples were not used to launch a pipeline or generate a ladder
for this documentation task. The coordinator is integrating and testing those
entry points independently. At initial inspection, global Click configuration
and the provenance helper existed, while pipeline forwarding/custom-reference
checks were still being edited. The guide documents the agreed completed
contract; the coordinator must verify the following before final freeze:

1. Explicit command options override file values for run and stage commands;
   `--no-report` overrides a true configured report value.
2. Global configured ladder min/max, dictionary/layout and flank values are
   forwarded; an explicitly provided ladder argument still wins.
3. Custom dictionary/layout/flanks in run require a matching explicit reference;
   supplied reference contents cannot be inferred from the settings alone.
4. Effective settings and input/configuration/reference hashes are written to
   `run_configuration.json` before external tools; early invalid paths/settings
   have visible failure status and cannot retain stale success artifacts.
5. `calling.read_phase` is a false-default experimental config opt-in, without a
   dedicated CLI flag. When WhatsHap is invoked, command/version/overrides remain
   recorded and None parameters do not add overrides.

The guide explicitly says confidence weights are heuristics, not calibrated
probabilities. It does not claim that arbitrary nondefault scientific settings
are validated. Evaluation acceptance endpoints, coordinate/genotype conventions,
unit-cost edit-distance semantics, codon size, supported-event identity and
failure handling remain fixed contracts rather than tuning parameters.

## Documentation build

`make docs-check` passed on 2026-09-14 (`mkdocs build --strict`, exit0). MkDocs
reported the new guide was not yet in navigation, as expected because the
coordinator owns that edit. The installed Material package also printed its
upstream MkDocs2 notice; no build warning or failure was reported. No final fresh
validation data, seeds, simulator, caller or external bioinformatics tools were
used.

File sizes: guide218 lines; generated example68 lines (below the649-line cap).

```text
648bed972503ec59b0c7a30510199e4fac41278d825b343abed383cb44f839b5  docs/guides/configuration.md
323eafb647cc77b0458df2061d4cd598bf9f1f42fac93c9c87b59e0d0959cbf5  examples/runtime-settings.json
```
