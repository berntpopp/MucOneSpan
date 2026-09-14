# MucOneSpan Rename Implementation Plan

**Goal:** Rename the repository and complete current codebase to MucOneSpan.
**Architecture:** Match MucOneUp naming: MucOneSpan display/repository, muc_one_span Python module/project metadata, muconespan CLI and container name. Retain scientific behavior and paper attribution by authors/title.
**Tech stack:** Python, uv, Click, MkDocs, GitHub Actions, Docker/BuildKit.
**Spec:** User-approved name and complete rename, following sibling conventions.

## Constraints

- No previous branding in current tracked content or filenames; preserve Git history.
- Python 3.10–3.14, existing scientific semantics and fewer than 650 authored physical lines.
- Publish the new identity as version 0.10.0 after passing checks.
- Preserve the optimized CI and Docker build design.

## Tasks

- [x] Rename Python source directory, imports, entry points, metadata, resources, tests, scripts and configuration; refresh lock metadata without unrelated upgrades.
- [x] Rewrite display names, repository/Pages links, installation examples, citation metadata, archived plans and container filenames; preserve original paper attribution.
- [x] Validate with make ci-check, docs-check, build-check, security-check, real-tool integration/full suite and Docker smoke.
- [x] Review diffs and scan every tracked path and text file for previous identity.
- [ ] Rename GitHub repository, update description/homepage/remotes, create PR and verify all Actions.
- [ ] Merge after checks pass, verify main publication, tag/release 0.10.0 and verify tag Actions and new container references.
- [ ] Synchronize and rename the local checkout directory; archive this completed plan.

## Validation evidence

- Fresh renamed-checkout environment: 239 unit tests passed, 86.11% branch-aware coverage, all lint/format/type/workflow/file-size checks pass.
- Full suite with external bioinformatics tools and MucOneUp fixtures: 255 passed, 2 documented strict expected failures.
- Isolated wheel installation verifies distribution metadata, console entry point, version, packaged data and report rendering.
- Strict documentation build and complete dependency audit pass.
- BuildKit and host Docker smoke pass, including exact Clair3 SNP calling.
- Independent review confirms normalized scientific ASTs and bundled reference data unchanged.
- All 119 current project files scanned: zero previous-name matches in content or filenames.
- Repository and local checkout renamed; remote points to the canonical MucOneSpan URL.
