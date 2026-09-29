# CLI Reference

Complete command-line interface documentation for MucOneSpan.

---

::: mkdocs-click
    :module: muc_one_span.cli
    :command: main
    :prog_name: muconespan
    :depth: 2

---

## Examples

See [Quick Start](../getting-started/quickstart.md) for complete usage examples.

### Hybrid engine (default)

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --assay amplicon
```

`run` stays registered as `muconespan run` (it moved internally to its own
module, `cli_run.py`, without changing the command or its import path).
The hybrid engine is the default since 0.17.0 and does not accept
`--report-igv`; see the
[configuration guide](../guides/configuration.md#hybrid-engine)
for every `hybrid.*` setting and dependency note, and
[Known Limitations](../reference/limitations.md#hybrid-engine)
for measured detection limits.

### Ladder engine (deprecated)

```bash
muconespan run \
  --engine ladder \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi
```

`--engine ladder` prints a deprecation warning and records it in
`summary.json["deprecations"]`. It is removed no earlier than the next minor
release; see the [migration guide](../guides/migration.md).
