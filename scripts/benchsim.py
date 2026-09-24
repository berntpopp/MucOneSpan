#!/usr/bin/env python3
"""MucSim-Bench: design and generate simulated MUC1 benchmark cases with MucOneUp.

Examples:
    python scripts/benchsim.py design --split dev --n 30
    python scripts/benchsim.py design --split test --salt-file ~/secrets/mucsim_salt.txt
    python scripts/benchsim.py generate --designs ../MucOneSpan-bench-data/designs_dev.jsonl \\
        --muconeup-config "$MUCONEUP_CONFIG" --muconeup-profiles <MucOneUp>/muc_one_up/data/read_profiles

Generated designs, truth and reads go to ``--out-root`` (default
``../MucOneSpan-bench-data``), outside Git. ``test`` designs need a secret salt
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
from pathlib import Path
from typing import Any

from muc_one_span.benchsim.design import Design, build_split
from muc_one_span.benchsim.generate import GenerateContext, generate_case, write_manifest
from muc_one_span.benchsim.muconeup import BUILTIN_PROFILE, require_muconeup
from muc_one_span.benchsim.profiles import builtin_profile_dir, write_variant
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.tools import run_tool

DEFAULT_OUT = Path("../MucOneSpan-bench-data")
DEFAULT_N = {"dev": 300, "val": 300, "test": 800, "stress": 100}
DEFAULT_SALT = "mucsim-bench-v1"


def repo_root() -> Path:
    """Top level of the Git checkout that contains this script."""
    here = Path(__file__).resolve().parent
    try:
        return Path(run_tool(["git", "rev-parse", "--show-toplevel"], cwd=str(here)).strip())
    except (FileNotFoundError, RuntimeError):
        return here.parent


def read_salt(args: argparse.Namespace) -> str:
    """Salt for seed derivation; ``test`` requires a salt file outside the repository."""
    if args.salt_file is None:
        if args.split == "test":
            raise SystemExit("--salt-file is required for the test split")
        return str(args.salt)
    path = Path(args.salt_file).expanduser().resolve()
    if path.is_relative_to(repo_root().resolve()):
        raise SystemExit(f"--salt-file must resolve outside the repository: {path}")
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
    out = args.output or args.out_root / f"designs_{args.split}.jsonl"
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
    config = args.muconeup_config or os.environ.get("MUCONEUP_CONFIG")
    if not config:
        raise SystemExit("--muconeup-config (or MUCONEUP_CONFIG) is required")
    lines = Path(args.designs).read_text().splitlines()
    designs = [Design.from_dict(json.loads(line)) for line in lines if line.strip()]
    version = require_muconeup(args.muconeup)
    profile_dir = builtin_profile_dir(args.muconeup_profiles)
    out_root = Path(args.out_root).resolve()
    for design in designs:  # write shared variant files once, before workers start
        base = profile_dir / f"{BUILTIN_PROFILE[design.profile]}.json"
        write_variant(base, design, out_root / "profiles")
    ctx = GenerateContext(
        args.muconeup,
        Path(config).resolve(),
        out_root,
        profile_dir,
        args.flank_fasta.resolve() if args.flank_fasta else None,
        args.structure_pool.resolve() if args.structure_pool else None,
        version,
    )
    cases = _run_all(designs, ctx, args.jobs)
    split_of = {d.design_id: d.split for d in designs}
    for split in sorted(set(split_of.values())):
        rows = [c for c in cases if split_of[c["design_id"]] == split]
        print(f"manifest: {write_manifest(out_root / split, rows)}")
    counts = Counter(c["status"] for c in cases)
    print("status: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 1 if counts.get("generation_failed") else 0


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
    design.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
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
    gen.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    gen.set_defaults(func=cmd_generate)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point."""
    args = parser().parse_args(argv)
    code: int = args.func(args)
    return code


if __name__ == "__main__":
    sys.exit(main())
