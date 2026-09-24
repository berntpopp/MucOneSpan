"""MucOneUp >= 0.45 per-read truth manifests ({base}_read_truth.tsv.gz).

Each manifest row records where one simulated read came from: which haplotype
and template molecule it was drawn from, its kind (full/fragment/smear/
chimera/offtarget/...), strand, source coordinates, and any homopolymer
sequencing-error edits applied. `realized_depth` turns those rows into the
spanning depth actually achieved per haplotype, distinct from the nominal
`--coverage` requested at simulation time.
"""

from __future__ import annotations

import gzip
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

COLUMNS = (
    "read_id",
    "hap",
    "molecule",
    "kind",
    "strand",
    "src_start",
    "src_end",
    "n_hp_edits",
    "hp_edits",
    "detail",
)


@dataclass(frozen=True)
class ReadTruth:
    """One row of a MucOneUp per-read truth manifest."""

    read_id: str
    hap: int
    molecule: int
    kind: str
    strand: str
    src_start: int
    src_end: int
    hp_edits: str
    detail: str


def load_read_truth(path: Path) -> list[ReadTruth]:
    """Parse a gzip-compressed per-read truth manifest.

    Args:
        path: Path to a `{base}_read_truth.tsv.gz` manifest.

    Returns:
        One `ReadTruth` per data row, in file order.

    Raises:
        ValueError: If the header does not match the MucOneUp >= 0.45 columns,
                    or a `read_id` appears more than once.
    """
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if tuple(header) != COLUMNS:
            raise ValueError(f"{path}: unexpected read truth header {header}")
        rows: list[ReadTruth] = []
        seen: set[str] = set()
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] in seen:
                raise ValueError(f"{path}: duplicate read id {f[0]}")
            seen.add(f[0])
            rows.append(
                ReadTruth(
                    f[0],
                    int(f[1]),
                    int(f[2]),
                    f[3],
                    f[4],
                    int(f[5]),
                    int(f[6]),
                    f[8],
                    f[9],
                )
            )
    return rows


def realized_depth(
    rows: Sequence[ReadTruth], span: dict[int, tuple[int, int]] | None
) -> dict[int, int]:
    """Count spanning reads actually realized per haplotype.

    Args:
        rows: Per-read truth rows from `load_read_truth`.
        span: `None` for amplicon reads, where depth is the count of
              `kind == "full"` reads per haplotype. Otherwise a mapping of
              haplotype to the `[start, end)` interval a fragment must cover
              to count as spanning.

    Returns:
        Mapping of haplotype id to realized spanning depth.
    """
    haps = sorted({r.hap for r in rows} | set(span or {}))
    depth = dict.fromkeys(haps, 0)
    for r in rows:
        if span is None:
            depth[r.hap] += r.kind == "full"
        elif r.kind == "fragment" and r.hap in span:
            lo, hi = span[r.hap]
            depth[r.hap] += r.src_start <= lo and r.src_end >= hi
    return depth


def composition(rows: Sequence[ReadTruth]) -> dict[str, Any]:
    """Summarize read kind/strand/haplotype composition of a truth manifest.

    Args:
        rows: Per-read truth rows from `load_read_truth`.

    Returns:
        Dict with `n_reads`, `kind_frac`, `strand_frac`, and `reads_per_hap`.
    """
    n = len(rows) or 1
    kinds = Counter(r.kind for r in rows)
    strands = Counter(r.strand for r in rows)
    haps = Counter(r.hap for r in rows)
    return {
        "n_reads": len(rows),
        "kind_frac": {k: v / n for k, v in sorted(kinds.items())},
        "strand_frac": {k: v / n for k, v in sorted(strands.items())},
        "reads_per_hap": dict(sorted(haps.items())),
    }
