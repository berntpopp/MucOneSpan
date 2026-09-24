"""S4 site table: per-read alleles at consensus sites, candidate sites and site events.

Ported from the prototype ``hetsplit.py`` (``site_table``, ``candidate_sites``) with two
spec-S4 additions: a candidate must be strand-consistent, and features that touch on
the consensus (adjacent columns, a homopolymer run and its neighbouring column, an
insertion slot and its neighbouring column) are merged into one event, so one sequence
change never counts as two linked sites.

Every tunable is a validated ``HybridSettings`` field taken from the ``settings``
argument; nothing is defaulted from ``DEFAULT_SETTINGS``.
"""

from __future__ import annotations

import statistics
from collections import Counter
from typing import Any

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.settings import HybridSettings

AF_DECIMALS = 3  # reporting precision of a site's allele fraction (format, not a tunable)
GAP = "-"  # deletion symbol in projected columns (``align.Columns``)
# Event grid (structural): column p sits at GRID_COLUMN * p and the insertion slot
# before column p half a column earlier, so features touch when <= one column apart.
GRID_COLUMN = 2
Site = tuple[str, int]
Meta = dict[Site, tuple[str, int]]


def features(
    cons: str, reads: list[str], settings: HybridSettings
) -> tuple[list[dict[Site, Any]], Meta]:
    """Per-read alleles at every column, insertion slot and homopolymer-run site."""
    runs = _runs(cons, settings.phase_run_min_len)
    in_run = {i for s, e, _ in runs for i in range(s, e)}
    feats = []
    for read in reads:
        proj = global_columns(read, cons)
        f: dict[Site, Any] = {}
        for pos, base in enumerate(proj.cols):
            if pos not in in_run:
                f[("col", pos)] = base
                if pos - 1 not in in_run:
                    f[("ins", pos)] = proj.ins.get(pos, "")
        for s, e, b in runs:
            f[("run", s)] = read_run_length(read, proj.t2q, s, e, b)
        feats.append(f)
    return feats, {("run", s): (b, e - s) for s, e, b in runs}


def _run_background(
    counts: dict[Site, Counter[Any]], meta: Meta, window: int
) -> dict[Site, dict[int, float]]:
    """Per run site, the noise fraction of reads reporting each nearby run length x.

    The background is the median, over the *other* runs of the same base and length, of
    the fraction of reads showing the same length shift (leave-one-out: a heterozygous
    run must not be its own background, which it would be if its length is unique in
    the consensus). With no such peer, the other runs of the same base are used; with
    none at all the background is 0 and only het_af_min applies.
    """
    shifts: dict[Site, dict[int, float]] = {}
    for site, (_base, length) in meta.items():
        c = counts.get(site, Counter())
        tot = sum(c.values()) or 1
        shifts[site] = {
            d: c.get(length + d, 0) / tot
            for d in range(-window, window + 1)
            if d and length + d >= 0
        }
    out: dict[Site, dict[int, float]] = {}
    for site, (base, length) in meta.items():
        others = [s for s, (b, _n) in meta.items() if s != site and b == base]
        peers = [s for s in others if meta[s][1] == length] or others
        out[site] = {}
        for d in shifts[site]:
            vals = [shifts[p][d] for p in peers if d in shifts[p]]
            out[site][length + d] = statistics.median(vals) if vals else 0.0
    return out


def _run_minor(
    c: Counter[Any], major: Any, noise: dict[Any, float], settings: HybridSettings
) -> tuple[Any, int] | None:
    """Best-scoring run length among all that clear max(het_af_min, mult * background)."""
    tot = sum(c.values())
    best: tuple[float, Any, int] | None = None
    for allele, n in c.items():
        if allele == major or n < settings.phase_min_minor_reads:
            continue
        bg = settings.phase_run_bg_multiplier * noise.get(allele, 0.0)
        if n / tot >= max(settings.het_af_min, bg) and (best is None or n / tot - bg > best[0]):
            best = (n / tot - bg, allele, n)
    return None if best is None else (best[1], best[2])


def _column_minor(c: Counter[Any], major: Any, settings: HybridSettings) -> tuple[Any, int] | None:
    """The top non-major allele if it clears het_af_min (gap alleles need more)."""
    rest = [(a, n) for a, n in c.most_common() if a != major]
    if not rest or rest[0][1] < settings.phase_min_minor_reads:
        return None
    minor, n_minor = rest[0]
    af, af_min = n_minor / sum(c.values()), settings.het_af_min
    if af < af_min or (GAP in (major, minor) and af < settings.phase_gap_af_factor * af_min):
        return None
    return minor, n_minor


def _strand_consistent(
    feats: list[dict[Site, Any]],
    strands: list[str],
    site: Site,
    minor: Any,
    settings: HybridSettings,
) -> bool:
    """Minor AF >= het_af_min on every strand with >= hp_min_strand_reads reads."""
    tally: dict[str, list[int]] = {}
    for f, strand in zip(feats, strands, strict=True):
        if site in f:
            t = tally.setdefault(strand, [0, 0])
            t[0] += 1
            t[1] += f[site] == minor
    return all(
        n_minor >= settings.het_af_min * tot
        for tot, n_minor in tally.values()
        if tot >= settings.hp_min_strand_reads
    )


def candidates(
    feats: list[dict[Site, Any]], strands: list[str], meta: Meta, settings: HybridSettings
) -> list[dict[str, Any]]:
    """Strand-consistent sites whose minor allele clears the run/column noise rules."""
    counts: dict[Site, Counter[Any]] = {}
    for f in feats:
        for site, allele in f.items():
            counts.setdefault(site, Counter())[allele] += 1
    bg = _run_background(counts, meta, settings.phase_run_bg_window)
    out = []
    for site, c in counts.items():
        major, _ = c.most_common(1)[0]
        if site[0] == "run":
            pick = _run_minor(c, major, bg[site], settings)
        else:
            pick = _column_minor(c, major, settings)
        if pick is None or not _strand_consistent(feats, strands, site, pick[0], settings):
            continue
        tot = sum(c.values())
        out.append(
            {
                "site": site,
                "major": major,
                "minor": pick[0],
                "af": round(pick[1] / tot, AF_DECIMALS),
                "n": tot,
            }
        )
    return out


def top_site(sites: list[dict[str, Any]]) -> dict[str, Any]:
    """The most balanced site (largest af * (1 - af)); the reference for orientation."""
    return max(sites, key=lambda s: s["af"] * (1 - s["af"]))


def _grid_span(site: Site, meta: Meta) -> tuple[int, int]:
    kind, pos = site
    if kind == "run":
        return GRID_COLUMN * pos, GRID_COLUMN * (pos + meta[site][1] - 1)
    grid = GRID_COLUMN * pos - (GRID_COLUMN // 2 if kind == "ins" else 0)
    return grid, grid


def events(sites: list[dict[str, Any]], meta: Meta) -> list[list[int]]:
    """Indices of ``sites`` grouped into events of touching features, in consensus order."""
    order = sorted(range(len(sites)), key=lambda i: _grid_span(sites[i]["site"], meta))
    out: list[list[int]] = []
    hi = 0
    for i in order:
        lo, end = _grid_span(sites[i]["site"], meta)
        if out and lo - hi <= GRID_COLUMN:
            out[-1].append(i)
            hi = max(hi, end)
        else:
            out.append([i])
            hi = end
    return out
