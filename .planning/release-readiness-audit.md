# Release readiness audit

Date: 2026-09-14. This is a read-only audit of the current release process and
working tree. It does not authorize or record a commit, push, tag, image publish,
or GitHub release.

## Current state and blockers

1. **A new version is required before tagging.** The package reports `0.10.0`,
   and signed annotated tag `v0.10.0` already exists at
   `78fa9cfcebedd268d7888677ba0a804ad2f638f3`. The current branch HEAD is
   `d8390b3c244ef8f3240af74b92db12b50dfc77d1`, with all production-validation
   changes still in the working tree. Tagging HEAD now would omit the work, and
   pushing another `v0.10.0` is not possible without rewriting an existing
   release tag.
2. **The tag workflow does not run or require the test workflow.**
   `.github/workflows/release.yml` runs for `v*`, checks only that the tag text is
   `v` plus `src/muc_one_span/version.py`, and creates a GitHub release with
   generated notes. `.github/workflows/test.yml` runs on selected branch pushes,
   pull requests, and manual dispatch, not tag pushes. The release job can finish
   independently of the CI Gate.
3. **A tag also publishes a container.** `.github/workflows/docker.yml` runs its
   BuildKit test stage and publishes GHCR tags derived from the semantic version,
   minor version, and commit SHA. It does not publish `latest`. The GitHub release
   and container workflows are independent, so the release can exist even if the
   container job later fails.
4. **The current local sdist must not be published.** `dist/muc_one_span-0.10.0.tar.gz`
   is 349 MB with 35,826 entries. It includes 35,136 `tests/results` entries, 260
   `.planning` entries, 91 review entries, local documentation/coverage outputs,
   and machine-specific paths. Of its planning entries, 243 are not tracked in
   Git. This happens because no Hatch sdist include/exclude policy is configured;
   the current `check_distribution.py` validates the wheel only.
5. **The wheel is bounded, but the source archive policy is incomplete.** The
   current 104 KB wheel has 44 package/dist-info entries and no `.planning` or
   `tests/results` content. It does include `muc_one_span/AGENTS.md` because that
   file lives inside the configured package directory. The Docker context
   allowlist excludes planning, results, and `AGENTS.md` files.

The release workflow does not build, upload, or publish the local wheel/sdist.
It creates a GitHub release only; GitHub supplies source archives from tracked tag
contents. Therefore ignored results and untracked review responses will not enter
the tag archive unless they are deliberately added to Git. The 17 already-tracked
planning files will remain in GitHub source archives. All `.planning/reviews`
files are currently untracked; do not stage them unless public release is an
explicit decision.

## Version decision and required edits

The change set adds a versioned runtime configuration API and substantial
backwards-compatible behavior/reporting changes. Under the repository's stated
Semantic Versioning policy, `0.11.0` is the clearest next version. If the
maintainer chooses a different version, use the same value consistently in:

- `src/muc_one_span/version.py` (`__version__` and dynamic package version);
- `CHANGELOG.md`: move current Unreleased content under `[0.11.0] - 2026-09-14`,
  retain an empty `[Unreleased]`, and add/update comparison links;
- `CITATION.cff`: `version` and `date-released`;
- `README.md` installation tag;
- `docs/getting-started/installation.md` installation tag and version example;
- `docs/about/citation.md` software citation version;
- `docker/muconespan.def` base GHCR tag and `%labels` version.

`pyproject.toml` is dynamic and has no literal package version. The editable root
entry in `uv.lock` also has no version, so a version-only edit does not require a
new dependency resolution. A locked check/build should still verify that wheel
metadata and CLI report the new version.

## Source-distribution decision

Before any wheel/sdist publication, add and review an explicit Hatch sdist policy
or build from an isolated clean checkout and inspect the archive. A sensible
policy should exclude at least `.planning`, `tests/results`, top-level `results`,
`site`, `htmlcov`, `coverage.xml`, and private tool/review artifacts. Building in
a clean worktree removes ignored/untracked data but still includes committed
planning unless the sdist policy excludes it.

The current GitHub release workflow does not need a distribution artifact. Do
not attach the existing local 0.10.0 artifacts to the new release. If PyPI or
GitHub wheel/sdist publication is intended, that is a separate missing release
capability: no trusted-publishing or artifact-upload step exists today.

## Release prerequisites

- Decide which planning/evidence files are public deliverables. Stage files by an
  explicit allowlist; leave private review prompts/responses, stderr captures,
  model-access records, ignored results, and machine-local artifacts untracked.
- Review the complete staged diff and confirm the release commit contains every
  new production module, test, guide, example, and changelog entry. The audit saw
  25 modified tracked paths and 80 untracked paths before this note was added.
- Make the version edits above and check that `git diff --check` remains clean.
- Run `make check` on the exact release commit. This covers quality, the unit
  coverage gate, strict docs, dependency audit, and package build. Run
  `make test-int` and the applicable cached generated-data checks separately;
  `make check` does not run generated end-to-end data.
- Inspect the final wheel file list and a clean-build sdist file list. A passing
  `make build-check` alone does not constrain sdist contents.
- Merge/push the release commit and wait for the branch `CI Gate` to succeed
  before pushing the tag. Record the exact commit SHA used for the tag.
- Confirm GHCR write permissions and monitor both `Release` and `Docker` tag
  workflows. Verify the GitHub release points to the intended commit and that the
  exact version container tag exists.

The coordinator separately reported 626 passing unit tests at 93.59% coverage,
47 passing integration tests, six generated end-to-end passes with two unchanged
expected failures, and a passing build before this audit. Those are useful
pre-release results but should be tied to the final release commit after version
and staging changes.

## Suggested commands after the release commit is prepared

```bash
# Verify version surfaces and the exact staged tree.
uv run --locked --all-extras python -c \
  'from muc_one_span import __version__; print(__version__)'
git diff --cached --check
git status --short

# Run gates on the exact release commit.
make check
make test-int

# Build in a clean worktree and inspect both artifacts.
git worktree add --detach /tmp/muconespan-release-check HEAD
cd /tmp/muconespan-release-check
make build-check
unzip -Z1 dist/*.whl
tar -tzf dist/*.tar.gz

# After the release commit is pushed and its branch CI Gate is green.
git tag -s v0.11.0 -m "MucOneSpan 0.11.0"
git show --no-patch --show-signature v0.11.0
git push origin v0.11.0
gh run list --branch v0.11.0
gh release view v0.11.0
```

Remove the temporary worktree after artifact inspection. Do not use a dirty local
build as the release artifact source.
