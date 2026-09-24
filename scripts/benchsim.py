#!/usr/bin/env python3
"""MucSim-Bench: design and generate simulated MUC1 benchmark cases with MucOneUp.

Examples:
    python scripts/benchsim.py design --split dev --n 30
    python scripts/benchsim.py design --split test --salt-file ~/secrets/mucsim_salt.txt
    python scripts/benchsim.py generate --designs ../MucOneSpan-bench-data/designs_dev.jsonl \\
        --muconeup-config "$MUCONEUP_CONFIG" --muconeup-profiles <MucOneUp>/muc_one_up/data/read_profiles

Generated designs, truth and reads go to ``--out-root`` (default
``<repository parent>/MucOneSpan-bench-data``); paths inside the repository
or any of its worktrees are refused. ``test`` designs need a secret salt
file that resolves outside the repository.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import ModuleType
from typing import Any

from muc_one_span.benchsim.design import Design, build_split
from muc_one_span.benchsim.generate import GenerateContext, generate_case, write_manifest
from muc_one_span.benchsim.muconeup import BUILTIN_PROFILE, require_muconeup
from muc_one_span.benchsim.profiles import builtin_profile_dir, write_variant
from muc_one_span.benchsim.realism import aggregate as realism_aggregate
from muc_one_span.benchsim.realism import case_metrics
from muc_one_span.benchsim.realism_targets import compare as realism_compare
from muc_one_span.benchsim.realism_targets import load_targets
from muc_one_span.benchsim.report import (
    RULE_TEXT,
    decide,
    normalize_rows,
    preregister,
    render_markdown,
    render_tables,
    require_preregistered,
    stratified_table,
)
from muc_one_span.benchsim.run_cases import run_split
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.tools import run_tool

DATA_DIR_NAME = "MucOneSpan-bench-data"
DEFAULT_N = {"dev": 300, "val": 300, "test": 800, "stress": 100}
DEFAULT_SALT = "mucsim-bench-v1"
HERE = Path(__file__).resolve().parent
TABLE_STRATA = (("profile",), ("profile", "delta_class"), ("profile", "depth"))
TABLE_METRICS = ("allele_exact", "critical_false_negative", "inconclusive", "no_call")


def _load_evaluate() -> ModuleType:
    """``scripts/evaluate.py`` (a script, not a package module)."""
    spec = spec_from_file_location("benchsim_evaluate", HERE / "evaluate.py")
    assert spec and spec.loader
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


evaluate_run = _load_evaluate().run


def _git_path(flag: str) -> Path | None:
    try:
        out = run_tool(["git", "rev-parse", flag], cwd=str(HERE)).strip()
    except (FileNotFoundError, RuntimeError):
        return None
    return (HERE / out).resolve() if out else None


def repo_root() -> Path:
    """Top level of the Git checkout that contains this script."""
    return _git_path("--show-toplevel") or HERE.parent


def protected_roots() -> list[Path]:
    """This checkout and the main checkout owning its Git directory (covers all worktrees)."""
    roots = [repo_root().resolve()]
    common = _git_path("--git-common-dir")
    if common is not None:
        roots.append(common.parent)
    return roots


def ensure_outside(path: Path, what: str) -> Path:
    """Resolve ``path``; exit if it lies inside the repository or any of its worktrees."""
    resolved = Path(path).expanduser().resolve()
    for root in protected_roots():
        if resolved.is_relative_to(root):
            raise SystemExit(f"{what} must resolve outside the repository ({root}): {resolved}")
    return resolved


def out_root(args: argparse.Namespace) -> Path:
    """``--out-root`` or ``<repo parent>/MucOneSpan-bench-data``, never inside Git."""
    default = repo_root().parent / DATA_DIR_NAME
    return ensure_outside(args.out_root or default, "--out-root")


def read_salt(args: argparse.Namespace) -> str:
    """Salt for seed derivation; ``test`` requires a salt file outside the repository."""
    if args.salt_file is None:
        if args.split == "test":
            raise SystemExit("--salt-file is required for the test split")
        return str(args.salt)
    path = ensure_outside(args.salt_file, "--salt-file")
    salt = path.read_text().strip()
    if not salt:
        raise SystemExit(f"empty salt file: {path}")
    return salt


def cmd_design(args: argparse.Namespace) -> int:
    """Write ``designs_<split>.jsonl``."""
    salt = read_salt(args)
    mutations = (
        args.mutations.split(",") if args.mutations else sorted(load_repeat_dictionary().mutations)
    )
    n = args.n if args.n is not None else DEFAULT_N[args.split]
    designs = build_split(args.split, n, salt, mutations)
    if args.output is not None:
        out = ensure_outside(args.output, "--output")
    else:
        out = out_root(args) / f"designs_{args.split}.jsonl"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(json.dumps(d.to_dict(), sort_keys=True) + "\n" for d in designs))
    print(f"wrote {len(designs)} designs to {out}")
    return 0


def _run_all(designs: list[Design], ctx: GenerateContext, jobs: int) -> list[dict[str, Any]]:
    if jobs <= 1:
        return [generate_case(d, ctx) for d in designs]
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(generate_case, designs, [ctx] * len(designs)))


def cmd_generate(args: argparse.Namespace) -> int:
    """Generate every design, write per-split manifests, exit 1 on generation failures."""
    root = out_root(args)
    config = args.muconeup_config or os.environ.get("MUCONEUP_CONFIG")
    if not config:
        raise SystemExit("--muconeup-config (or MUCONEUP_CONFIG) is required")
    lines = Path(args.designs).read_text().splitlines()
    designs = [Design.from_dict(json.loads(line)) for line in lines if line.strip()]
    version = require_muconeup(args.muconeup)
    profile_dir = builtin_profile_dir(args.muconeup_profiles)
    for design in designs:  # write shared variant files once, before workers start
        base = profile_dir / f"{BUILTIN_PROFILE[design.profile]}.json"
        write_variant(base, design, root / "profiles")
    ctx = GenerateContext(
        args.muconeup,
        Path(config).resolve(),
        root,
        profile_dir,
        args.flank_fasta.resolve() if args.flank_fasta else None,
        args.structure_pool.resolve() if args.structure_pool else None,
        version,
    )
    cases = _run_all(designs, ctx, args.jobs)
    split_of = {d.design_id: d.split for d in designs}
    for split in sorted(set(split_of.values())):
        rows = [c for c in cases if split_of[c["design_id"]] == split]
        print(f"manifest: {write_manifest(root / split, rows)}")
    counts = Counter(c["status"] for c in cases)
    print("status: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 1 if counts.get("generation_failed") else 0


def _model_lookup(models: dict[str, str | None], platform: str) -> str:
    """Picklable ``model_for`` (a bound closure is not, and breaks ``--jobs`` > 1)."""
    model = models.get(platform)
    if not model:
        raise SystemExit(
            f"no --model-{platform} (or $CLAIR3_MODEL_{platform.upper()}) for platform {platform!r}"
        )
    return model


def cmd_run(args: argparse.Namespace) -> int:
    """Run every engine over a split's manifest, keeping every case's denominator."""
    manifest = Path(args.manifest)
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")
    engines = [e.strip() for e in args.engines.split(",") if e.strip()]
    if not engines:
        raise SystemExit("--engines must name at least one engine")
    split = manifest.resolve().parent.name
    default_results = manifest.resolve().parent.parent / "results" / split
    results_root = ensure_outside(args.results_root or default_results, "--results-root")
    models: dict[str, str | None] = {
        "ont": args.model_ont or os.environ.get("CLAIR3_MODEL_ONT"),
        "hifi": args.model_hifi or os.environ.get("CLAIR3_MODEL_HIFI"),
    }
    records = run_split(
        manifest, engines, results_root, partial(_model_lookup, models), args.threads, args.jobs
    )
    counts = Counter((r["engine"], r["status"]) for r in records)
    for engine in engines:
        line = ", ".join(f"{k}={v}" for (e, k), v in sorted(counts.items()) if e == engine)
        print(f"{engine}: {line}")
    print(f"results: {results_root}")
    return 0


def _prereg_path(root: Path) -> Path:
    return root / "test" / "preregistration.jsonl"


def _guard_sealed(split: str, root: Path) -> None:
    """``test`` truth stays sealed until the decision rule is pre-registered."""
    if split != "test":
        return
    try:
        require_preregistered(_prereg_path(root), RULE_TEXT)
    except PermissionError as exc:
        raise SystemExit(f"test split is sealed: {exc}") from exc


def _results_root(args: argparse.Namespace, root: Path) -> Path:
    return ensure_outside(args.results_root or root / "results" / args.split, "--results-root")


def _engines(text: str) -> list[str]:
    engines = [e.strip() for e in text.split(",") if e.strip()]
    if not engines:
        raise SystemExit("--engines must name at least one engine")
    return engines


def _write(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def cmd_preregister(args: argparse.Namespace) -> int:
    """Append the decision rule to ``<out-root>/test/preregistration.jsonl``."""
    path = _prereg_path(out_root(args))
    print(f"pre-registered rule sha256 {preregister(RULE_TEXT, path)} in {path}")
    return 0


def cmd_evaluate(args: argparse.Namespace) -> int:
    """Score each engine's results with ``scripts/evaluate.py`` into ``evaluation.json``."""
    root = out_root(args)
    _guard_sealed(args.split, root)  # before any truth is read
    results, split_dir = _results_root(args, root), root / args.split
    code = 0
    for engine in _engines(args.engines):
        engine_dir = results / engine
        ns = argparse.Namespace(
            result_root=engine_dir,
            truth_root=split_dir,
            expected_samples=engine_dir / f"inventory_{args.split}.json",
        )
        report, rc = evaluate_run(ns)
        _write(engine_dir / "evaluation.json", report)
        print(f"{engine}: {engine_dir / 'evaluation.json'} (exit {rc})")
        code = max(code, rc)
    return code


def _cases(split_dir: Path) -> dict[str, dict[str, Any]]:
    manifest = split_dir / "manifest.jsonl"
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")
    rows = [json.loads(line) for line in manifest.read_text().splitlines() if line.strip()]
    return {row["design_id"]: row for row in rows}


def _tables(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    tables = {}
    for by in TABLE_STRATA:
        for metric in TABLE_METRICS:
            tables[f"{metric} by {' x '.join(by)}"] = stratified_table(rows, by, metric)
        normals = [r for r in rows if r["normal"] or r["benign"]]
        tables[f"false_positive (normal+benign) by {' x '.join(by)}"] = stratified_table(
            normals, by, "false_positive"
        )
    return tables


def cmd_report(args: argparse.Namespace) -> int:
    """Stratified tables per engine and the decision rule: ``report.json`` + ``report.md``."""
    root = out_root(args)
    _guard_sealed(args.split, root)
    results, cases = _results_root(args, root), _cases(root / args.split)
    engines = [args.baseline] + ([args.candidate] if args.candidate else [])
    rows: dict[str, list[dict[str, Any]]] = {}
    for engine in engines:
        path = results / engine / "evaluation.json"
        if not path.is_file():
            raise SystemExit(f"missing {path}; run `benchsim evaluate` first")
        rows[engine] = normalize_rows(json.loads(path.read_text()), cases)
    decision = decide(rows, args.baseline, args.candidate) if args.candidate else None
    tables = {engine: _tables(r) for engine, r in rows.items()}
    report = {
        "split": args.split,
        "decision": decision,
        "tables": tables,
        "engines": {engine: {"rows": r} for engine, r in rows.items()},
    }
    _write(results / "report.json", report)
    parts = [render_markdown(decision or {"profiles": {}, "adopt": False})]
    parts += [f"## Engine `{e}`\n\n{render_tables(t)}" for e, t in tables.items()]
    (results / "report.md").write_text("\n".join(parts))
    print(f"report: {results / 'report.json'}")
    return 0


def cmd_realism(args: argparse.Namespace) -> int:
    """Task 9 realism metrics over a split, aggregated and compared per profile."""
    root = out_root(args)
    _guard_sealed(args.split, root)
    split_dir = root / args.split
    per_profile: dict[str, list[dict[str, Any]]] = {}
    failures = []
    for design_id, case in sorted(_cases(split_dir).items()):
        if case.get("status") != "ok":
            failures.append({"design_id": design_id, "error": f"status {case.get('status')}"})
            continue
        try:
            metrics = case_metrics(split_dir / design_id, args.muconeup_config, args.flank_fasta)
        except (OSError, ValueError, KeyError) as exc:
            failures.append({"design_id": design_id, "error": f"{type(exc).__name__}: {exc}"})
            continue
        per_profile.setdefault(case["profile"], []).append(metrics)
    targets = load_targets()
    profiles: dict[str, Any] = {}
    for profile, cases in sorted(per_profile.items()):
        agg = realism_aggregate(cases)
        try:
            checks = realism_compare(agg, targets, profile)
        except KeyError:
            checks = None  # no public target section for this profile
        profiles[profile] = {"n_cases": len(cases), "aggregate": agg, "compare": checks}
    _write(
        split_dir / "realism.json",
        {"split": args.split, "profiles": profiles, "failures": failures},
    )
    lines = [f"# Realism: {args.split}", ""]
    for profile, res in profiles.items():
        lines += [f"## {profile} ({res['n_cases']} cases)", ""]
        if res["compare"] is None:
            lines += ["No public target section.", ""]
            continue
        lines += ["| check | pass |", "|---|---|"]
        lines += [f"| {k} | {v['pass']} |" for k, v in res["compare"].items()] + [""]
    lines += [f"Failures: {len(failures)}"]
    (split_dir / "realism.md").write_text("\n".join(lines) + "\n")
    print(f"realism: {split_dir / 'realism.json'} ({len(failures)} failures)")
    return 0


def parser() -> argparse.ArgumentParser:
    """Command-line interface."""
    result = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = result.add_subparsers(dest="action", required=True)
    design = commands.add_parser("design", help="write designs_<split>.jsonl")
    design.add_argument("--split", choices=sorted(DEFAULT_N), required=True)
    design.add_argument("--n", type=int, help="cases per profile (default: split size)")
    design.add_argument("--salt", default=DEFAULT_SALT, help="public salt for non-test splits")
    design.add_argument("--salt-file", type=Path, help="secret salt file (required for test)")
    design.add_argument("--mutations", help="comma-separated event names (default: dictionary)")
    design.add_argument("--out-root", type=Path, help=f"default: <repo parent>/{DATA_DIR_NAME}")
    design.add_argument("--output", type=Path, help="default: <out-root>/designs_<split>.jsonl")
    design.set_defaults(func=cmd_design)
    gen = commands.add_parser("generate", help="simulate truth and reads for designs")
    gen.add_argument("--designs", type=Path, required=True)
    gen.add_argument("--jobs", type=int, default=1, help="parallel cases (process pool)")
    gen.add_argument("--muconeup", default="muconeup", help="MucOneUp executable")
    gen.add_argument("--muconeup-config", type=Path, help="default: $MUCONEUP_CONFIG")
    gen.add_argument("--muconeup-profiles", type=Path, help="MucOneUp data/read_profiles dir")
    gen.add_argument("--flank-fasta", type=Path, help="extra left/right flanks (genomic)")
    gen.add_argument("--structure-pool", type=Path, help="local real-derived structures")
    gen.add_argument("--out-root", type=Path, help=f"default: <repo parent>/{DATA_DIR_NAME}")
    gen.set_defaults(func=cmd_generate)
    run_cmd = commands.add_parser("run", help="run caller engines over a split's manifest")
    run_cmd.add_argument(
        "--manifest", type=Path, required=True, help="<out-root>/<split>/manifest.jsonl"
    )
    run_cmd.add_argument("--engines", default="ladder", help="comma-separated (default: ladder)")
    run_cmd.add_argument("--results-root", type=Path, help="default: <out-root>/results/<split>")
    run_cmd.add_argument("--model-ont", help="caller model for ont (default: $CLAIR3_MODEL_ONT)")
    run_cmd.add_argument("--model-hifi", help="caller model for hifi (default: $CLAIR3_MODEL_HIFI)")
    run_cmd.add_argument("--threads", type=int, default=4)
    run_cmd.add_argument(
        "--jobs", type=int, default=1, help="parallel (engine, case) pairs (process pool)"
    )
    run_cmd.set_defaults(func=cmd_run)
    prereg = commands.add_parser("preregister", help="pre-register the decision rule for test")
    prereg.add_argument("--out-root", type=Path, help=f"default: <repo parent>/{DATA_DIR_NAME}")
    prereg.set_defaults(func=cmd_preregister)
    for name, func, text in (
        ("evaluate", cmd_evaluate, "score engine results into <engine>/evaluation.json"),
        ("report", cmd_report, "stratified tables and decision rule: report.json + report.md"),
        ("realism", cmd_realism, "realism metrics vs targets: realism.json + realism.md"),
    ):
        sp = commands.add_parser(name, help=text)
        sp.add_argument("--split", choices=sorted(DEFAULT_N), required=True)
        sp.add_argument("--out-root", type=Path, help=f"default: <repo parent>/{DATA_DIR_NAME}")
        sp.set_defaults(func=func)
        if name == "evaluate":
            sp.add_argument("--engines", default="ladder", help="comma-separated engines")
        if name == "report":
            sp.add_argument("--baseline", default="ladder")
            sp.add_argument("--candidate", help="engine to decide on (omit: tables only)")
        if name in ("evaluate", "report"):
            sp.add_argument("--results-root", type=Path, help="default: <out-root>/results/<split>")
        if name == "realism":
            sp.add_argument("--muconeup-config", type=Path, help="only for cases without geometry")
            sp.add_argument("--flank-fasta", type=Path, help="flanks used at generation (genomic)")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point."""
    args = parser().parse_args(argv)
    code: int = args.func(args)
    return code


if __name__ == "__main__":
    sys.exit(main())
