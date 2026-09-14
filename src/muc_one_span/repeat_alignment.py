"""Sequence alignment and indel accounting for repeat classification."""

from __future__ import annotations


def edit_distance(s1: str, s2: str) -> int:
    """Compute Levenshtein edit distance between two sequences.

    Args:
        s1: First sequence.
        s2: Second sequence.

    Returns:
        Minimum number of single-character edits (insert, delete, substitute).
    """
    return _edit_distance_bitvector(s1, s2)


def _edit_distance_bitvector(s1: str, s2: str) -> int:
    """Compute exact Levenshtein distance with Myers' bit-vector recurrence."""
    if s1 == s2:
        return 0
    if not s1:
        return len(s2)
    if not s2:
        return len(s1)

    # Fewer pattern positions mean smaller integer bit operations. Distance is
    # symmetric, so the shorter input can always be used as the pattern.
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    pattern_masks: dict[str, int] = {}
    for position, symbol in enumerate(s1):
        pattern_masks[symbol] = pattern_masks.get(symbol, 0) | (1 << position)

    positive = ~0
    negative = 0
    distance = len(s1)
    final_bit = 1 << (len(s1) - 1)

    for symbol in s2:
        matches = pattern_masks.get(symbol, 0)
        vertical = matches | negative
        horizontal = (((matches & positive) + positive) ^ positive) | matches
        positive_horizontal = negative | ~(horizontal | positive)
        negative_horizontal = positive & horizontal

        if positive_horizontal & final_bit:
            distance += 1
        elif negative_horizontal & final_bit:
            distance -= 1

        positive_horizontal = (positive_horizontal << 1) | 1
        negative_horizontal <<= 1
        positive = negative_horizontal | ~(vertical | positive_horizontal)
        negative = positive_horizontal & vertical

    return distance


def characterize_differences(ref: str, query: str) -> list[dict]:
    """Characterize specific differences between reference and query sequences.

    Uses Needleman-Wunsch style traceback to identify individual
    substitutions, insertions, and deletions with positions.

    Args:
        ref: Reference sequence.
        query: Query sequence.

    Returns:
        List of difference dicts with keys: pos, ref, alt, type.
    """
    if ref == query:
        return []

    m, n = len(ref), len(query)

    # Build full DP matrix for traceback
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if ref[i - 1] == query[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1])

    # Traceback
    diffs: list[dict] = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and ref[i - 1] == query[j - 1]:
            i -= 1
            j -= 1
        elif i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + 1:
            # Substitution
            diffs.append(
                {
                    "pos": i,  # 1-based position in reference
                    "ref": ref[i - 1],
                    "alt": query[j - 1],
                    "type": "substitution",
                }
            )
            i -= 1
            j -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + 1:
            # Insertion in query
            diffs.append(
                {
                    "pos": i + 1,  # 1-based reference position before insertion
                    "ref": "",
                    "alt": query[j - 1],
                    "type": "insertion",
                }
            )
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + 1:
            # Deletion from reference
            diffs.append(
                {
                    "pos": i,
                    "ref": ref[i - 1],
                    "alt": "",
                    "type": "deletion",
                }
            )
            i -= 1
        else:
            break

    diffs.reverse()
    return diffs


def _compute_net_indel(diffs: list[dict]) -> int:
    """Compute net indel offset from a list of differences.

    Insertions add bases (positive offset), deletions remove bases
    (negative offset).  The net offset tells us how much the sequence
    shifted relative to the reference frame.

    Returns:
        Net indel length: positive = sequence is longer, negative = shorter.
    """
    offset = 0
    for d in diffs:
        if d["type"] == "insertion":
            offset += len(d["alt"])
        elif d["type"] == "deletion":
            offset -= len(d["ref"])
    return offset
