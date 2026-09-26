#!/usr/bin/env python3
"""MucSim-Bench: design and generate simulated MUC1 benchmark cases with MucOneUp.

Examples:
    python scripts/benchsim.py design --split dev --n 30 --set standard
    python scripts/benchsim.py design --split test --salt-file ~/secrets/mucsim_salt.txt
    python scripts/benchsim.py generate --designs ../MucOneSpan-bench-data/designs_dev_standard.jsonl \\
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
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import ModuleType
from typing import Any

from muc_one_span.benchsim.atlas import build_atlas, render_atlas
from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, load_bench_config
from muc_one_span.benchsim.calibration import (
    CALIBRATION_SPLIT,
    CalibrationRequest,
    EvaluateFn,
    check_split,
    run_calibration,
)
from muc_one_span.benchsim.calibration_cli import add_calibration_parsers
from muc_one_span.benchsim.calibration_grid import LENGTHS_STAGE
from muc_one_span.benchsim.calibration_lengths import evaluate_lengths_stage, run_lengths_stage
from muc_one_span.benchsim.calibration_report import build_report
from muc_one_span.benchsim.design import Design, build_split
from muc_one_span.benchsim.generate import (
    GenerateContext,
    StaleCaseError,
    generate_case,
    write_manifest,
)
from muc_one_span.benchsim.muconeup import BUILTIN_PROFILE, require_muconeup
from muc_one_span.benchsim.profiles import builtin_profile_dir, write_variant
from muc_one_span.benchsim.provenance import (
    caller_record,
    git_commit,
    report_provenance,
    write_caller,
)
from muc_one_span.benchsim.realism import aggregate as realism_aggregate
from muc_one_span.benchsim.realism import case_metrics
from muc_one_span.benchsim.realism_targets import INDICATIVE_NOTE, load_targets, target_note
from muc_one_span.benchsim.realism_targets import compare as realism_compare
from muc_one_span.benchsim.report import (
    decide,
    first_evaluation,
    mark_first_evaluation,
    normalize_rows,
    preregister,
    render_markdown,
    require_preregistered,
    rule_text,
)
from muc_one_span.benchsim.report_tables import build_tables, render_engine_tables
from muc_one_span.benchsim.run_cases import run_split
from muc_one_span.benchsim.targets import render_targets, targets_by_set
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.settings import DEFAULT_SETTINGS
from muc_one_span.tools import run_tool

DATA_DIR_NAME = "MucOneSpan-bench-data"
SPLITS = sorted(DEFAULT_BENCH_CONFIG.design.split_sizes)
DEFAULT_SALT = "mucsim-bench-v1"
HERE = Path(__file__).resolve().parent
DEFAULT_ENGINE = DEFAULT_SETTINGS.run.engine


def _load_evaluate() -> ModuleType:
    """``scripts/evaluate.py`` (a script, not a package module)."""
    spec = spec_from_file_location("benchsim_evaluate", HERE / "evaluate.py")
    assert spec and spec.loader
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ``scripts/evaluate.py::run``; loaded lazily by ``cmd_evaluate`` (tests patch it).
evaluate_run: Callable[[argparse.Namespace], tuple[dict[str, Any], int]] | None = None


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
    """Write ``designs_<split>_<set>.jsonl``."""
    salt = read_salt(args)
    mutations = (
        args.mutations.split(",") if args.mutations else sorted(load_repeat_dictionary().mutations)
    )
    sizes = args.bench.design.split_sizes
    if args.split not in sizes:
        raise SystemExit(f"no size configured for split {args.split!r}")
    bench_set = args.set or args.bench.sets.default
    if bench_set not in args.bench.sets.definitions:
        known = ", ".join(args.bench.sets.definitions)
        raise SystemExit(f"unknown set {bench_set!r} (configured: {known})")
    n = args.n if args.n is not None else sizes[args.split]
    try:
        designs = build_split(args.split, n, salt, mutations, args.bench, bench_set)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if args.output is not None:
        out = ensure_outside(args.output, "--output")
    else:
        out = out_root(args) / f"designs_{args.split}_{bench_set}.jsonl"
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
        write_variant(base, design, root / "profiles", args.bench.profiles)
    ctx = GenerateContext(
        args.muconeup,
        Path(config).resolve(),
        root,
        profile_dir,
        args.flank_fasta.resolve() if args.flank_fasta else None,
        args.structure_pool.resolve() if args.structure_pool else None,
        version,
        args.bench,
    )
    try:
        cases = _run_all(designs, ctx, args.jobs)
    except StaleCaseError as exc:
        raise SystemExit(str(exc)) from exc
    split_of = {d.design_id: d.split for d in designs}
    for split in sorted(set(split_of.values())):
        rows = [c for c in cases if split_of[c["design_id"]] == split]
        print(f"manifest: {write_manifest(root / split, rows)}")
    counts = Counter(c["status"] for c in cases)
    print("status: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 1 if counts.get("generation_failed") else 0


def _models(args: argparse.Namespace) -> dict[str, str | None]:
    return {
        "ont": args.model_ont or os.environ.get("CLAIR3_MODEL_ONT"),
        "hifi": args.model_hifi or os.environ.get("CLAIR3_MODEL_HIFI"),
    }


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
    threads = args.threads if args.threads is not None else args.bench.run.threads
    extra = {"config": args.config.resolve()} if args.config else {}
    model_for = partial(_model_lookup, _models(args))
    records = run_split(manifest, engines, results_root, model_for, threads, args.jobs, **extra)
    commit = git_commit(HERE, run_tool)
    for engine in engines:  # which caller build produced these results
        write_caller(results_root / engine, caller_record(engine, commit))
    counts = Counter((r["engine"], r["status"]) for r in records)
    for engine in engines:
        line = ", ".join(f"{k}={v}" for (e, k), v in sorted(counts.items()) if e == engine)
        print(f"{engine}: {line}")
    print(f"results: {results_root}")
    return 0


def _rule(args: argparse.Namespace) -> str:
    """The decision rule under the effective settings (headline set and targets included)."""
    return rule_text(args.bench.report, args.bench.sets.headline, args.bench.targets)


def _prereg_path(root: Path) -> Path:
    return root / "test" / "preregistration.jsonl"


def _guard_sealed(split: str, root: Path, rule: str) -> dict[str, Any] | None:
    """``test`` truth stays sealed until the decision rule is pre-registered.

    Returns the matched pre-registration entry for ``test`` (``None`` otherwise).
    """
    if split != "test":
        return None
    try:
        return require_preregistered(_prereg_path(root), rule)
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
    try:
        digest = preregister(_rule(args), path)
    except PermissionError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"pre-registered rule sha256 {digest} in {path}")
    return 0


def cmd_evaluate(args: argparse.Namespace) -> int:
    """Score each engine's results with ``scripts/evaluate.py`` into ``evaluation.json``."""
    root = out_root(args)
    audit = _guard_sealed(args.split, root, _rule(args))  # before truth is read
    if audit is not None:  # marked before scoring starts (conservative: a failed run counts)
        audit["test_first_evaluated_at"] = mark_first_evaluation(
            _prereg_path(root), audit["sha256"]
        )
    results, split_dir = _results_root(args, root), root / args.split
    run = evaluate_run or _load_evaluate().run
    code = 0
    for engine in _engines(args.engines):
        engine_dir = results / engine
        ns = argparse.Namespace(
            result_root=engine_dir,
            truth_root=split_dir,
            expected_samples=engine_dir / f"inventory_{args.split}.json",
        )
        report, rc = run(ns)
        if audit is not None:
            report["preregistration"] = audit
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


def _set_sections(
    args: argparse.Namespace, rows: dict[str, list[dict[str, Any]]]
) -> dict[str, dict[str, Any]]:
    """Tables and reason atlas per benchmark set, headline set first."""
    sets = args.bench.sets
    present = {r["bench_set"] for engine_rows in rows.values() for r in engine_rows}
    out: dict[str, dict[str, Any]] = {}
    for name in sets.order(present):
        sub = {
            e: [r for r in engine_rows if r["bench_set"] == name] for e, engine_rows in rows.items()
        }
        definition = sets.definitions.get(name)
        out[name] = {
            "headline": name == sets.headline,
            "description": definition.description if definition else None,
            "n_cases": len(sub[args.baseline]),
            "tables": {e: build_tables(r, args.bench.report) for e, r in sub.items()},
            "atlas": {e: build_atlas(r, args.split, args.bench.atlas) for e, r in sub.items()},
        }
    return out


def cmd_report(args: argparse.Namespace) -> int:
    """Tables and atlas per set and engine, rule on the headline set: report.json + report.md."""
    root = out_root(args)
    audit = _guard_sealed(args.split, root, _rule(args))
    if audit is not None:
        audit["test_first_evaluated_at"] = first_evaluation(_prereg_path(root))
    results, cases = _results_root(args, root), _cases(root / args.split)
    engines = [args.baseline] + ([args.candidate] if args.candidate else [])
    rows: dict[str, list[dict[str, Any]]] = {}
    for engine in engines:
        path = results / engine / "evaluation.json"
        if not path.is_file():
            raise SystemExit(f"missing {path}; run `benchsim evaluate` first")
        try:
            rows[engine] = normalize_rows(
                json.loads(path.read_text()), cases, args.bench.sets.legacy
            )
        except (KeyError, ValueError) as exc:
            raise SystemExit(f"cannot normalize {path}: {exc}") from exc
    headline = args.bench.sets.headline
    decision, header = None, {}
    if args.candidate and any(r["bench_set"] == headline for r in rows[args.baseline]):
        try:
            decision = decide(rows, args.baseline, args.candidate, args.bench)
        except ValueError as exc:
            raise SystemExit(
                f"`{args.baseline}` and `{args.candidate}` were scored on different cases "
                f"({exc}); evaluate both engines on the same manifest"
            ) from exc
        header = decision
    elif args.candidate:
        header = {
            "note": f"Decision: not evaluated (no `{headline}` cases; the rule applies "
            f"to the `{headline}` set only)."
        }
    sets = _set_sections(args, rows)
    alpha = args.bench.report.alpha
    targets = {e: targets_by_set(r, args.bench.targets, alpha) for e, r in rows.items()}
    report = {
        "split": args.split,
        "preregistration": audit,
        "provenance": report_provenance(
            git_commit(HERE, run_tool), {e: results / e for e in engines}, cases
        ),
        "bench_config_sha256": args.bench.sha256(),
        "decision": decision,
        "headline_set": headline,
        "set_order": list(sets),
        "sets": sets,
        "targets": targets,  # every engine, every targeted set; absent sets not present
        "engines": {engine: {"rows": r} for engine, r in rows.items()},
    }
    _write(results / "report.json", report)
    parts = [render_markdown(header)]
    if decision is None:
        parts += [
            f"## Absolute targets: engine `{e}` (descriptive; no decision)\n\n{render_targets(t)}"
            for e, t in targets.items()
            if render_targets(t)
        ]
    for name, section in sets.items():
        label = f"Set `{name}`" + (" (headline)" if section["headline"] else "")
        about = f"{section['description']}\n\n" if section["description"] else ""
        parts += [
            f"## {label}: engine `{e}` ({section['n_cases']} cases)\n\n{about}"
            f"{render_atlas(section['atlas'][e], args.bench.atlas)}"
            f"{render_engine_tables(t, args.bench.report)}"
            for e, t in section["tables"].items()
        ]
    (results / "report.md").write_text("\n".join(parts))
    print(f"report: {results / 'report.json'}")
    return 0


def cmd_realism(args: argparse.Namespace) -> int:
    """Task 9 realism metrics over a split, aggregated and compared per profile."""
    root = out_root(args)
    audit = _guard_sealed(args.split, root, _rule(args))
    if audit is not None:  # realism reads test truth: unseals under this rule
        mark_first_evaluation(_prereg_path(root), audit["sha256"])
    split_dir = root / args.split
    per_set: dict[str, dict[str, list[dict[str, Any]]]] = {}
    failures = []
    for design_id, case in sorted(_cases(split_dir).items()):
        if case.get("status") != "ok":
            failures.append({"design_id": design_id, "error": f"status {case.get('status')}"})
            continue
        try:
            metrics = case_metrics(
                split_dir / design_id, args.muconeup_config, args.flank_fasta, args.bench.realism
            )
        except (OSError, ValueError, KeyError, IndexError) as exc:
            failures.append({"design_id": design_id, "error": f"{type(exc).__name__}: {exc}"})
            continue
        name = (case.get("design") or {}).get("bench_set") or args.bench.sets.legacy
        per_set.setdefault(name, {}).setdefault(case["profile"], []).append(metrics)
    targets = load_targets()
    sets: dict[str, Any] = {}
    for name in args.bench.sets.order(per_set):
        profiles: dict[str, Any] = {}
        for profile, cases in sorted(per_set[name].items()):
            agg = realism_aggregate(cases, args.bench.realism)
            try:
                checks = realism_compare(agg, targets, profile, args.bench.realism)
            except KeyError:
                checks = None  # no public target section for this profile
            except ValueError as exc:
                raise SystemExit(f"realism settings do not match the targets: {exc}") from exc
            note = target_note(targets, profile)
            profiles[profile] = {
                "n_cases": len(cases),
                "note": note,
                "aggregate": agg,
                "compare": checks,
            }
        sets[name] = {"profiles": profiles}
    _write(
        split_dir / "realism.json",
        {
            "split": args.split,
            "bench_config_sha256": args.bench.sha256(),
            "indicative": INDICATIVE_NOTE,
            "sets": sets,
            "failures": failures,
        },
    )
    lines = [f"# Realism: {args.split}", "", INDICATIVE_NOTE, ""]
    for name, section in sets.items():
        lines += [f"## Set `{name}`", ""]
        for profile, res in section["profiles"].items():
            lines += [f"### {profile} ({res['n_cases']} cases)", "", f"Scope: {res['note']}.", ""]
            if res["compare"] is None:
                lines += ["No public target section.", ""]
                continue
            lines += ["| check | pass |", "|---|---|"]
            lines += [f"| {k} | {v['pass']} |" for k, v in res["compare"].items()] + [""]
    lines += [f"Failures: {len(failures)}"]
    (split_dir / "realism.md").write_text("\n".join(lines) + "\n")
    print(f"realism: {split_dir / 'realism.json'} ({len(failures)} failures)")
    return 0


def cmd_calibrate(args: argparse.Namespace) -> int:
    """Run and score every point of a settings grid (resumable; test refused).

    ``--stage lengths`` (Task 15d) fits only the hybrid length model on each case's
    spanning reads instead of running the full pipeline, so it needs ``--engine
    hybrid`` (the length model is a hybrid-engine concept; refused otherwise, before
    any point is run).
    """
    check_split(args.split)
    if args.stage == LENGTHS_STAGE and args.engine != "hybrid":
        raise SystemExit("--stage lengths requires --engine hybrid")
    evaluate_fn: EvaluateFn
    if args.stage == LENGTHS_STAGE:
        evaluate_fn = evaluate_lengths_stage

        def run_point(manifest: Path, results: Path, config: Path) -> Any:
            return run_lengths_stage(manifest, results, config, engine=args.engine)
    else:
        threads = args.threads if args.threads is not None else args.bench.run.threads
        model_for = partial(_model_lookup, _models(args))
        evaluate_fn = evaluate_run or _load_evaluate().run

        def run_point(manifest: Path, results: Path, config: Path) -> Any:
            return run_split(
                manifest, [args.engine], results, model_for, threads, args.jobs, config=config
            )

    request = CalibrationRequest(
        args.split,
        args.name or args.grid.stem,
        args.engine,
        args.grid,
        args.config,
        out_root(args),
        stage=args.stage,
    )
    return run_calibration(request, run_point, evaluate_fn)


def cmd_calibrate_report(args: argparse.Namespace) -> int:
    """Rank a calibration under OBJECTIVE.json and write the recommended config."""
    check_split(args.split)
    root = out_root(args)
    shift = (args.shift_from, _cases(root / CALIBRATION_SPLIT)) if args.shift_from else None
    cases = _cases(root / args.split)
    return build_report(root, args.split, args.name, args.objective, cases, args.bench, shift)


def parser() -> argparse.ArgumentParser:
    """Command-line interface."""
    result = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    result.add_argument(
        "--bench-config",
        type=Path,
        help="JSON overriding MucSim-Bench settings (default: built-in, docs/benchmark.md)",
    )
    commands = result.add_subparsers(dest="action", required=True)
    design = commands.add_parser("design", help="write designs_<split>.jsonl")
    design.add_argument("--split", choices=SPLITS, required=True)
    design.add_argument("--n", type=int, help="cases per profile (default: split size)")
    design.add_argument("--set", help="benchmark set (default: bench config sets.default)")
    design.add_argument("--salt", default=DEFAULT_SALT, help="public salt for non-test splits")
    design.add_argument("--salt-file", type=Path, help="secret salt file (required for test)")
    design.add_argument("--mutations", help="comma-separated event names (default: dictionary)")
    design.add_argument("--out-root", type=Path, help=f"default: <repo parent>/{DATA_DIR_NAME}")
    design.add_argument(
        "--output", type=Path, help="default: <out-root>/designs_<split>_<set>.jsonl"
    )
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
    run_cmd.add_argument(
        "--engines", default=DEFAULT_ENGINE, help=f"comma-separated (default: {DEFAULT_ENGINE})"
    )
    run_cmd.add_argument("--results-root", type=Path, help="default: <out-root>/results/<split>")
    run_cmd.add_argument("--model-ont", help="caller model for ont (default: $CLAIR3_MODEL_ONT)")
    run_cmd.add_argument("--model-hifi", help="caller model for hifi (default: $CLAIR3_MODEL_HIFI)")
    run_cmd.add_argument("--threads", type=int, help="default: bench config run.threads")
    run_cmd.add_argument("--config", type=Path, help="runtime settings for every run (--config)")
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
        sp.add_argument("--split", choices=SPLITS, required=True)
        sp.add_argument("--out-root", type=Path, help=f"default: <repo parent>/{DATA_DIR_NAME}")
        sp.set_defaults(func=func)
        if name == "evaluate":
            sp.add_argument("--engines", default=DEFAULT_ENGINE, help="comma-separated engines")
        if name == "report":
            sp.add_argument("--baseline", default="ladder")
            sp.add_argument("--candidate", help="engine to decide on (omit: tables only)")
        if name in ("evaluate", "report"):
            sp.add_argument("--results-root", type=Path, help="default: <out-root>/results/<split>")
        if name == "realism":
            sp.add_argument("--muconeup-config", type=Path, help="only for cases without geometry")
            sp.add_argument("--flank-fasta", type=Path, help="flanks used at generation (genomic)")
    out_help = f"default: <repo parent>/{DATA_DIR_NAME}"
    add_calibration_parsers(commands, SPLITS, out_help, cmd_calibrate, cmd_calibrate_report)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point."""
    args = parser().parse_args(argv)
    try:
        args.bench = load_bench_config(args.bench_config)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"--bench-config: {exc}") from exc
    code: int = args.func(args)
    return code


if __name__ == "__main__":
    sys.exit(main())
