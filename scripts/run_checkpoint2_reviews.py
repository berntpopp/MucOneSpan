#!/usr/bin/env python3
"""Run independent external model reviews for Checkpoint 2 (Architecture Design)."""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REVIEWS_DIR = ROOT / ".planning" / "reviews" / "checkpoint2_design"
REVIEWS_DIR.mkdir(parents=True, exist_ok=True)

REQUEST_FILE = REVIEWS_DIR / "review-request-design.md"
DESIGN_FILE = ROOT / ".planning" / "2026-09-14-architectural-design.md"
PROMPT_FILE = REVIEWS_DIR / "prompt.txt"


def prepare_prompt() -> None:
    request_text = REQUEST_FILE.read_text()
    design_text = DESIGN_FILE.read_text()
    prompt = f"""You are acting as an exceptionally rigorous principal bioinformatician, software architect, and algorithm engineer.

Please perform an independent, in-depth, highly critical review of the following architectural design proposal for MucOneSpan VNTR reconstruction.

=== REVIEW REQUEST & QUESTIONS ===
{request_text}

=== ARCHITECTURAL DESIGN PROPOSAL ===
{design_text}

Please provide your rigorous critique, answering each specific question, highlighting any subtle failure modes, edge cases, or mathematical risks, and recommending concrete adjustments.
"""
    PROMPT_FILE.write_text(prompt)


def run_claude() -> None:
    print("[Reviewer 1/2] Invoking Claude Code (model: claude-fable-5-1)...")
    out_json = REVIEWS_DIR / "claude-response.json"
    err_log = REVIEWS_DIR / "claude-stderr.log"

    cmd = [
        "claude",
        "-p",
        "--model",
        "claude-fable-5-1",
        "--effort",
        "high",
        "--permission-mode",
        "plan",
        "--tools",
        "Read,Glob,Grep",
        "--output-format",
        "json",
    ]
    t0 = time.monotonic()
    with PROMPT_FILE.open("r") as fin, out_json.open("w") as fout, err_log.open("w") as ferr:
        res = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=ferr, cwd=str(ROOT))
    elapsed = time.monotonic() - t0
    print(f"  Claude finished in {elapsed:.1f}s (Exit code: {res.returncode})")


def run_codex() -> None:
    print("[Reviewer 2/2] Invoking Codex CLI (model: gpt-6-astra)...")
    out_jsonl = REVIEWS_DIR / "codex-events.jsonl"
    err_log = REVIEWS_DIR / "codex-stderr.log"

    cmd = [
        "codex",
        "exec",
        "--model",
        "gpt-6-astra",
        "--sandbox",
        "read-only",
        "--json",
        "-",
    ]
    t0 = time.monotonic()
    with PROMPT_FILE.open("r") as fin, out_jsonl.open("w") as fout, err_log.open("w") as ferr:
        res = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=ferr, cwd=str(ROOT))
    elapsed = time.monotonic() - t0
    print(f"  Codex finished in {elapsed:.1f}s (Exit code: {res.returncode})")


def main() -> int:
    prepare_prompt()
    # Run independently
    run_claude()
    run_codex()
    print("=== All Checkpoint 2 reviews completed! ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
