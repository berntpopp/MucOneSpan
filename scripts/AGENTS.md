# Script guidance

- Scripts are maintained Python: Ruff, formatting, mypy, and the maximum 649-line
  rule apply here too. Annotate entry points and use a `__main__` guard.
- Expose configuration through arguments or documented environment variables;
  do not assume a developer's home directory or conda installation path.
- Reuse package logic where practical. Preserve command arguments, file formats,
  sample seeds, and benchmark interpretation during maintenance refactors.
- Exit nonzero on command or validation failure. Avoid `shell=True`, command
  interpolation, and silently treating failed samples as successful results.
- Write large generated data/results to ignored directories. Document external
  prerequisites and expected output locations in the benchmarking guide.
- Run `make ci-check`, plus focused script tests/help or representative fixture
  runs as appropriate. Do not regenerate expensive datasets for unrelated edits.
