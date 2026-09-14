# Packaging leak fix evidence

## Problem reproduced

The pre-fix local `muc_one_span-0.10.0.tar.gz` was **349 MB** with **35,826**
archive entries. The new validator rejected it with **35,703 allowlist
violations**. Its contents included 35,136 `tests/results` entries, 260
`.planning` entries, 91 review entries, generated sequencing/tool outputs,
indexes, coverage/site artifacts, and machine-specific paths. Of the planning
entries, 243 were untracked. The 104 KB wheel did not contain planning or result
artifacts.

The test-first focused run failed during collection before implementation because
`validate_sdist` and `validate_wheel` did not exist. Synthetic tar/zip regressions
then covered generated results, BAM/FASTQ data including compressed FASTQ,
reference indexes, bytecode, private planning content, parent traversal, tar
symbolic/hard links, and zip symbolic links.

## Implemented boundary

`pyproject.toml` now gives Hatch an explicit sdist allowlist for maintained
package source, root test fixtures plus unit/integration tests, scripts, docs,
examples, and release/build metadata (`README`, `LICENSE`, `CHANGELOG`,
`CITATION.cff`, `Makefile`, `mkdocs.yml`, `pyproject.toml`, `uv.lock`, and optional
`.python-version`). Generated reads/alignments, indexes, and bytecode are excluded.
`.planning`, generated test data/results, local site/coverage output, repository
workflows, and container/development material are outside the allowlist.

`scripts/check_distribution.py` now requires exactly one wheel and one sdist,
checks safe archive paths and member types, rejects links and generated/private
content, requires packaged runtime resources, and retains the isolated wheel
installation, version, CLI, and report smoke tests. Hatch force-includes the
maintained root `.gitignore` as VCS metadata; this single backend-added file is
explicitly permitted by the validator.

## Validation

```text
uv run --locked --all-extras pytest --no-cov -q \
  tests/unit/test_distribution_contents.py
9 passed in 0.02s

uv run --locked --all-extras ruff check \
  scripts/check_distribution.py tests/unit/test_distribution_contents.py
All checks passed!

uv run --locked --all-extras ruff format --check \
  scripts/check_distribution.py tests/unit/test_distribution_contents.py
2 files already formatted

uv run --locked --all-extras mypy scripts/check_distribution.py
Success: no issues found in 1 source file

make build-check
Successfully built dist/muc_one_span-0.11.0.tar.gz
Successfully built dist/muc_one_span-0.11.0-py3-none-any.whl
muconespan, version 0.11.0
Wheel/sdist contents, wheel installation, CLI and bundled resources passed
```

The final inspected artifacts at that source state were:

| Artifact | Bytes | Entries | Forbidden entries | SHA-256 |
| --- | ---: | ---: | ---: | --- |
| `muc_one_span-0.11.0.tar.gz` | 315,862 | 130 | 0 | `b89a339a053646253dfc4837d0325e2b21dbbd97cb7243f15751a674b6a7f453` |
| `muc_one_span-0.11.0-py3-none-any.whl` | 111,480 | 45 | 0 | `66f2a4f6e972ced3938d0a4b353d1de3ea6d9d2c01661c7404e227189b478331` |

The sdist therefore changed from 349 MB/35,826 entries with 35,703 violations to
315,862 bytes/130 entries with none. These hashes are evidence for the packaging
fix run, not final release hashes: source changed afterward for the experiments
and documentation workstreams. The release coordinator must run `make
build-check` again after the final source freeze and inspect the rebuilt archive
counts/hashes before tagging or uploading artifacts.

The historical `.planning/release-readiness-audit.md` correctly recorded the
then-current blockers. Its version blocker is now resolved by the coordinated
`0.11.0` update, and its unsafe-sdist blocker is resolved by this allowlist and
validator. Retain that audit as historical evidence; final release gates,
staging/commit review, clean artifact rebuild, tag verification, and workflow
monitoring remain coordinator responsibilities.

## Release rebuild after the source freeze

`make build-check` passed again after all production and public documentation edits.
The wheel and sdist content validator and isolated wheel installation/CLI/resource
checks passed. Final review can still require a new rebuild if source changes.

| Artifact | Bytes | Archive entries | SHA-256 |
| --- | ---: | ---: | --- |
| `muc_one_span-0.11.0-py3-none-any.whl` | 111,701 | 45 | `1fd722f1b61b8e0fbbca6a5afd6317c6fda62cdae2d1d4d811ea0ad538a7a9ad` |
| `muc_one_span-0.11.0.tar.gz` | 320,780 | 132 | `f11c3a41d5ef7163f7063c9fab40ced413e5e6011c96ccf1c9f69e3559ab9fd3` |
