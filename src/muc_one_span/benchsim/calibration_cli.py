"""Argument parsers of ``benchsim calibrate`` and ``benchsim calibrate-report``.

The command functions live in ``scripts/benchsim.py`` (they need its
out-root guard, model lookup and evaluator loader); this module only declares
the arguments, so the script stays within the file-size limit.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
from pathlib import Path

from muc_one_span.benchsim.calibration_grid import DEFAULT_STAGE, STAGES

Command = Callable[[argparse.Namespace], int]


def add_calibration_parsers(
    commands: argparse._SubParsersAction[argparse.ArgumentParser],
    splits: Sequence[str],
    out_root_help: str,
    calibrate: Command,
    calibrate_report: Command,
) -> None:
    """Register both subcommands on the benchsim subparser group."""
    cal = commands.add_parser(
        "calibrate", help="run and score a settings grid: calibration/<split>/<name>/"
    )
    cal.add_argument("--split", choices=splits, required=True, help="dev or val (test is sealed)")
    cal.add_argument("--engine", required=True, help="engine the overlays select (run.engine)")
    cal.add_argument("--grid", type=Path, required=True, help="GRID.json: settings key -> values")
    cal.add_argument(
        "--stage",
        choices=STAGES,
        default=DEFAULT_STAGE,
        help="'lengths': fit only the hybrid length model (fast; needs --engine hybrid)",
    )
    cal.add_argument("--config", type=Path, help="base runtime settings (default: built-in)")
    cal.add_argument("--name", help="calibration name (default: the grid file stem)")
    cal.add_argument("--out-root", type=Path, help=out_root_help)
    cal.add_argument("--model-ont", help="caller model for ont (default: $CLAIR3_MODEL_ONT)")
    cal.add_argument("--model-hifi", help="caller model for hifi (default: $CLAIR3_MODEL_HIFI)")
    cal.add_argument("--threads", type=int, help="default: bench config run.threads")
    cal.add_argument("--jobs", type=int, default=1, help="parallel cases per point (process pool)")
    cal.set_defaults(func=calibrate)
    rep = commands.add_parser(
        "calibrate-report", help="rank a calibration under OBJECTIVE.json; recommend a config"
    )
    rep.add_argument("--split", choices=splits, required=True, help="dev or val (test is sealed)")
    rep.add_argument("--name", required=True, help="calibration name")
    rep.add_argument("--objective", type=Path, required=True, help="OBJECTIVE.json")
    rep.add_argument("--shift-from", help="dev calibration name to compare with (--split val)")
    rep.add_argument("--out-root", type=Path, help=out_root_help)
    rep.set_defaults(func=calibrate_report)
