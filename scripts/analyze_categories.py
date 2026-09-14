#!/usr/bin/env python3
import json
from pathlib import Path

cand_path = Path("tests/results/candidate_dev_140/evaluation_report.json")
cand = json.loads(cand_path.read_text()) if cand_path.exists() else {"samples": []}

categories: dict[str, dict[str, int]] = {}
for s in cand["samples"]:
    name = s["sample"]
    metrics = s["metrics"]
    seq_exact = metrics["sequence_exact"]["min"]
    all_seq_exact = metrics["all_sequences_exact"]["min"]
    count_exact = metrics["count_exact"]["min"]
    all_count_exact = metrics["all_counts_exact"]["min"]

    parts = name.split("_")
    cat = parts[1]  # asym, eq, eqvar, gap1, gap2, gap3, bnd, typ, etc.
    if cat not in categories:
        categories[cat] = {
            "total": 0,
            "all_seq_exact": 0,
            "seq_exact": 0,
            "all_count_exact": 0,
            "count_exact": 0,
        }
    categories[cat]["total"] += 1
    categories[cat]["all_seq_exact"] += all_seq_exact
    categories[cat]["seq_exact"] += seq_exact
    categories[cat]["all_count_exact"] += all_count_exact
    categories[cat]["count_exact"] += count_exact

header = (
    f"{'Category':<15} | {'Total':<5} | {'AllSeq':<6} | {'SeqExact':<8} | "
    f"{'AllCount':<8} | {'CountExact':<10}"
)
print(header)
print("-" * len(header))
for cat, stats in sorted(categories.items()):
    print(
        f"{cat:<15} | {stats['total']:5d} | {stats['all_seq_exact']:6d} | "
        f"{stats['seq_exact']:8d} | {stats['all_count_exact']:8d} | {stats['count_exact']:10d}"
    )
