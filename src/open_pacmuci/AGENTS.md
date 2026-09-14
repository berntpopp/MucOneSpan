# Package guidance

- Keep CLI parsing/orchestration distinct from scientific computation and file
  formats. Use existing module boundaries described in `docs/development.md`.
- Maintain Python 3.10 compatibility and annotate changed APIs. Preserve public
  import paths when splitting modules, including existing test patch targets.
- Use `tools.py` for external commands; retain its environment cleanup and error
  handling. Avoid capturing unbounded sequencing data in new subprocess wrappers.
- Resolve packaged repeat definitions, references, and templates through package
  resources. Test packaging when changing resources; editable installs can mask
  missing wheel files.
- Treat coordinate bases, repeat counts, allele assignment, and mutation naming
  as scientific contracts. Verify changes against focused fixtures with known
  sequences and expected coordinates. Do not silently tune thresholds in refactors.
- Check both HiFi and ONT paths when changing mapping/calling defaults, and retain
  explicit preset overrides. Test report changes with the optional report extra.
- Run `make ci-check`; add integration/build checks where the changed behavior
  depends on bioinformatics tools or packaging.
