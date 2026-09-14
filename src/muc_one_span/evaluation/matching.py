"""Maximum-cardinality, minimum global sequence-distance allele assignment."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations, permutations

from muc_one_span.repeat_alignment import edit_distance

from .models import Assignment, PredictedAllele, TruthHaplotype


def match_alleles(
    truth: Sequence[TruthHaplotype], predictions: Sequence[PredictedAllele]
) -> Assignment:
    """Return every optimal assignment; event labels never affect the optimum.

    Intended for diploid/small-ploidy evaluation. Missing/extra alleles are
    represented by indices absent from each alternative, including no-call cases.
    """
    if not truth:
        raise ValueError("at least one truth haplotype is required")
    matrix = tuple(tuple(edit_distance(t.sequence, p.sequence) for p in predictions) for t in truth)
    n = min(len(truth), len(predictions))
    best: int | None = None
    alternatives = []
    for ts in combinations(range(len(truth)), n):
        for ps in permutations(range(len(predictions)), n):
            pairs = tuple(zip(ts, ps, strict=True))
            cost = sum(matrix[ti][pi] for ti, pi in pairs)
            if best is None or cost < best:
                best, alternatives = cost, [pairs]
            elif cost == best:
                alternatives.append(pairs)
    return Assignment(best if best is not None else 0, tuple(alternatives), matrix)
