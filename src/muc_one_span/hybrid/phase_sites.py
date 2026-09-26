"""S4 site table: per-read alleles at consensus sites, candidate sites and site events.

Ported from the prototype ``hetsplit.py`` (``site_table``, ``candidate_sites``) with two
spec-S4 additions: a candidate must show no strand bias (a strand-bias test, not a
per-strand AF floor, so unbiased imbalanced heterozygotes pass), and features that touch on
the consensus (adjacent columns, a homopolymer run and its neighbouring column, an
insertion slot and its neighbouring column) are merged into one event, so one sequence
change never counts as two linked sites.

Every tunable is a validated ``HybridSettings`` field taken from the ``settings``
argument; nothing is defaulted from ``DEFAULT_SETTINGS``.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from typing import Any

from muc_one_span.hybrid.align import global_columns
from muc_one_span.hybrid.polish import _runs, read_run_length
from muc_one_span.hybrid.run_strand import run_mixture
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


def site_counts(feats: list[dict[Site, Any]]) -> dict[Site, Counter[Any]]:
    """Per site, how many reads show each allele."""
    counts: dict[Site, Counter[Any]] = {}
    for f in feats:
        for site, allele in f.items():
            counts.setdefault(site, Counter())[allele] += 1
    return counts


def modal_meta(counts: dict[Site, Counter[Any]], meta: Meta) -> Meta:
    """Run sites keyed to their most frequently observed length, not the draft's.

    A POA draft can carry a run one base short or long (heavy stutter); taking such a
    run's draft length as its true length would make its reads look like stutter and
    poison the background and stutter profiles of every run it is a peer of. The
    reads' modal length is the run's best length estimate. Runs no read observes keep
    their draft length.
    """
    return {
        site: (base, counts[site].most_common(1)[0][0] if counts.get(site) else length)
        for site, (base, length) in meta.items()
    }


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
    c: Counter[Any],
    major: Any,
    noise: dict[Any, float],
    settings: HybridSettings,
    multiplier: float,
) -> tuple[Any, int] | None:
    """Best-scoring run length among all that clear max(het_af_min, multiplier * background).

    ``multiplier`` is ``phase_run_bg_multiplier`` for a candidate site and
    ``phase_run_safety_multiplier`` for the NEGATIVE-blocking tier (``run_excess_sites``).
    """
    tot = sum(c.values())
    best: tuple[float, Any, int] | None = None
    for allele, n in c.items():
        if allele == major or n < settings.phase_min_minor_reads:
            continue
        bg = multiplier * noise.get(allele, 0.0)
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


def _depletion_p(n_minor: int, n_strand: int, minor_total: int, total: int) -> float:
    """One-sided Fisher exact p: P(<= n_minor minor reads among n_strand strand reads).

    Hypergeometric lower tail of the allele x strand 2x2 table, conditional on its
    margins (``total`` reads of which ``minor_total`` carry the minor allele).
    """
    lo = max(0, n_strand - (total - minor_total))
    tail = sum(
        math.comb(minor_total, i) * math.comb(total - minor_total, n_strand - i)
        for i in range(lo, n_minor + 1)
    )
    return tail / math.comb(total, n_strand)


def _absent_on_a_strand(
    feats: list[dict[Site, Any]], strands: list[str], site: Site, minor: Any, s: HybridSettings
) -> tuple[bool, dict[str, list[int]]]:
    """(minor absent from a strand with >= hp_min_strand_reads reads, strand tally)."""
    tally: dict[str, list[int]] = {}
    for f, strand in zip(feats, strands, strict=True):
        if site in f:
            t = tally.setdefault(strand, [0, 0])
            t[0] += 1
            t[1] += f[site] == minor
    absent = any(m == 0 and n >= s.hp_min_strand_reads for n, m in tally.values())
    return absent, tally


def _column_consistent(
    feats: list[dict[Site, Any]], strands: list[str], site: Site, minor: Any, s: HybridSettings
) -> bool:
    """False when a column/insertion site's minor allele shows strand bias.

    Biased means: the minor is absent from a strand with >= hp_min_strand_reads reads,
    or a one-sided Fisher exact test finds it depleted on a strand at
    phase_strand_bias_alpha. An unbiased imbalanced heterozygote passes.
    """
    absent, tally = _absent_on_a_strand(feats, strands, site, minor, s)
    total = sum(t[0] for t in tally.values())
    minor_total = sum(t[1] for t in tally.values())
    return not absent and all(
        _depletion_p(m, n, minor_total, total) >= s.phase_strand_bias_alpha
        for n, m in tally.values()
    )


def _run_consistent(
    feats: list[dict[Site, Any]],
    strands: list[str],
    meta: Meta,
    site: Site,
    pick: tuple[int, int],
    s: HybridSettings,
) -> bool:
    """False when a run site's minor length is stutter or strand-biased (``run_strand``).

    ``pick`` is (major, minor) length. The site is kept only when the minor is observed
    on every strand with >= hp_min_strand_reads reads, its stutter-deconvolved mixture
    weight reaches het_af_min, and the stutter-aware strand test is not significant at
    phase_strand_bias_alpha (strand-asymmetric stutter is not strand bias). A run
    whose base forms no other run in the consensus has no stutter model and is tested
    like a column.
    """
    base = meta[site][0]
    if not any(b == base for p, (b, _n) in meta.items() if p != site):
        # No other run of this base: its stutter cannot be modelled (background 0 in
        # _run_background), so the column rules apply to the observed lengths.
        return _column_consistent(feats, strands, site, pick[1], s)
    if _absent_on_a_strand(feats, strands, site, pick[1], s)[0]:
        return False
    weight, p = run_mixture(feats, strands, meta, site, pick, s)
    return weight >= s.het_af_min and p >= s.phase_strand_bias_alpha


def candidates(
    feats: list[dict[Site, Any]], strands: list[str], meta: Meta, settings: HybridSettings
) -> list[dict[str, Any]]:
    """Strand-consistent sites whose minor allele clears the run/column noise rules."""
    counts = site_counts(feats)
    modal = modal_meta(counts, meta)
    bg = _run_background(counts, modal, settings.phase_run_bg_window)
    out = []
    for site, c in counts.items():
        major, _ = c.most_common(1)[0]
        if site[0] == "run":
            pick = _run_minor(c, major, bg[site], settings, settings.phase_run_bg_multiplier)
            keep = pick is not None and _run_consistent(
                feats, strands, modal, site, (major, pick[0]), settings
            )
        else:
            pick = _column_minor(c, major, settings)
            keep = pick is not None and _column_consistent(feats, strands, site, pick[0], settings)
        if pick is None or not keep:
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


def run_excess_sites(
    feats: list[dict[Site, Any]], meta: Meta, settings: HybridSettings
) -> list[dict[str, Any]]:
    """Run sites whose minor length clears the lower safety floor (Task 15g).

    The floor is ``max(het_af_min, phase_run_safety_multiplier * background)`` with the
    same leave-one-out peer background as a candidate site, but with none of the
    candidate's other tests (split floor, stutter-deconvolved weight, strand bias).
    A site here is never split on or scored as an event; the engine uses it only to
    keep an unsplit equal-length peak from a negative call, because a heterozygous run
    below the split floor (heavy stutter at a true longer run) would otherwise merge
    both alleles into one consensus. Largest excess over the floor first.
    """
    counts = site_counts(feats)
    modal = modal_meta(counts, meta)
    bg = _run_background(counts, modal, settings.phase_run_bg_window)
    mult = settings.phase_run_safety_multiplier
    out = []
    for site, c in counts.items():
        if site[0] != "run":
            continue
        major, _ = c.most_common(1)[0]
        pick = _run_minor(c, major, bg[site], settings, mult)
        if pick is None:
            continue
        tot = sum(c.values())
        excess = pick[1] / tot - mult * bg[site].get(pick[0], 0.0)
        out.append(
            (
                excess,
                {
                    "site": site,
                    "major": major,
                    "minor": pick[0],
                    "af": round(pick[1] / tot, AF_DECIMALS),
                    "n": tot,
                },
            )
        )
    out.sort(key=lambda item: -item[0])
    return [site for _excess, site in out]


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
