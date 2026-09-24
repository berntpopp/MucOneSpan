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

### Hybrid engine (experimental)

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --engine hybrid \
  --assay amplicon
```

`run` stays registered as `muconespan run` (it moved internally to its own
module, `cli_run.py`, without changing the command or its import path).
`--engine hybrid` needs the optional `hybrid` extra and does not accept
`--report-igv`; see the
[configuration guide](../guides/configuration.md#hybrid-engine-experimental)
for every `hybrid.*` setting and dependency note, and
[Known Limitations](../reference/limitations.md#hybrid-engine-experimental)
for measured detection limits.
