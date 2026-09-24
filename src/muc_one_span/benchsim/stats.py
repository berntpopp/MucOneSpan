"""Statistics for the realistic benchmark: exact CIs, McNemar, Holm, non-inferiority.

scipy is not a project dependency, so the regularized incomplete beta function
(needed for the exact Clopper-Pearson interval) is implemented directly with
`math.lgamma` plus a continued-fraction expansion (Numerical Recipes
`betacf`/`betai`), inverted by bisection. The standard-normal quantile (needed
for the Newcombe hybrid-score non-inferiority bound) uses Acklam's rational
approximation to the inverse normal CDF. Both are pure `math`, no `random`
beyond `cluster_bootstrap`'s explicit `random.Random(seed)`.

Reference values (Clopper-Pearson bounds, exact McNemar p-value) were checked
against R 4.5.2 `binom.test` before the tests were written; see
task-8-report.md for the commands and output.
"""

from __future__ import annotations

import math
import random
import statistics
from collections.abc import Callable, Sequence
from typing import Any

_MAX_ITER = 200
_EPS = 3e-16
_FPMIN = 1e-300
_BISECT_ITER = 80


def _betacf(a: float, b: float, x: float) -> float:
    """Continued-fraction expansion for the incomplete beta function (Numerical Recipes)."""
    qab, qap, qam = a + b, a + 1, a - 1
    c = 1.0
    d = 1.0 - qab * x / qap
    d = _FPMIN if abs(d) < _FPMIN else d
    d = 1.0 / d
    h = d
    for m in range(1, _MAX_ITER + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = _FPMIN if abs(d) < _FPMIN else d
        c = 1.0 + aa / c
        c = _FPMIN if abs(c) < _FPMIN else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = _FPMIN if abs(d) < _FPMIN else d
        c = 1.0 + aa / c
        c = _FPMIN if abs(c) < _FPMIN else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < _EPS:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta function I_x(a, b) (the Beta(a, b) CDF at x)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _beta_quantile(a: float, b: float, target: float) -> float:
    """Bisect for p such that `_betai(a, b, p) == target` (a Beta(a, b) quantile)."""
    lo, hi = 0.0, 1.0
    for _ in range(_BISECT_ITER):
        mid = (lo + hi) / 2
        if _betai(a, b, mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact (Clopper-Pearson) two-sided binomial confidence interval for k/n."""
    lower = 0.0 if k == 0 else _beta_quantile(k, n - k + 1, alpha / 2)
    upper = 1.0 if k == n else _beta_quantile(k + 1, n - k, 1 - alpha / 2)
    return lower, upper


_LOG_HALF = math.log(0.5)
_MCNEMAR_TOL = 1e-9  # absolute log-probability slack so symmetric pmf[i]/pmf[n-i]
# ties (equal up to lgamma's floating-point rounding) are never split apart.


def _log_binom_pmf(n: int, k: int) -> float:
    """log P(X = k) for X ~ Binomial(n, 0.5), via math.lgamma (no big-int conversion)."""
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + n * _LOG_HALF


def mcnemar_exact(b: int, c: int) -> float:
    """Exact two-sided McNemar test: binomial(n=b+c, p=0.5) test on the discordant pairs.

    Two-sided p-value is the sum of all Binomial(n, 0.5) point probabilities
    that are <= the observed count's probability (matches R `binom.test`'s
    two-sided method), capped at 1.0. Computed in log-space (`math.lgamma`,
    as `_betai` already does) rather than via `math.comb(n, i) * 0.5**n`:
    for n above ~1030, `math.comb`'s exact integer result is too large to
    convert to a Python float at all (`OverflowError`), which the log-space
    binomial pmf never hits (every log-pmf value is <= 0, so `math.exp` of
    it never overflows).
    """
    n = b + c
    if n == 0:
        return 1.0
    log_pmf = [_log_binom_pmf(n, i) for i in range(n + 1)]
    observed = log_pmf[b]
    total = sum(math.exp(lp) for lp in log_pmf if lp <= observed + _MCNEMAR_TOL)
    return min(total, 1.0)


def holm(pvalues: dict[str, float]) -> dict[str, float]:
    """Holm-Bonferroni step-down adjustment: monotone, capped at 1.0."""
    n = len(pvalues)
    order = sorted(pvalues, key=lambda name: pvalues[name])
    adjusted: dict[str, float] = {}
    running_max = 0.0
    for i, name in enumerate(order):
        running_max = max(running_max, (n - i) * pvalues[name])
        adjusted[name] = min(running_max, 1.0)
    return adjusted


def _norm_ppf(p: float) -> float:
    """Inverse standard-normal CDF (probit), Acklam's rational approximation (~1e-9 error)."""
    a = (
        -3.969683028665376e01,
        2.209460984245205e02,
        -2.759285104469687e02,
        1.383577518672690e02,
        -3.066479806614716e01,
        2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01,
        1.615858368580409e02,
        -1.556989798598866e02,
        6.680131188771972e01,
        -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e00,
        -2.549732539343734e00,
        4.374664141464968e00,
        2.938163982698783e00,
    )
    d = (
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e00,
        3.754408661907416e00,
    )
    p_low, p_high = 0.02425, 1 - 0.02425
    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1
        )
    if p <= p_high:
        q = p - 0.5
        r = q * q
        num = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        den = ((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1
        return num / den
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
        (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1
    )


def _wilson_bounds(x: int, n: int, z: float) -> tuple[float, float]:
    """Wilson score interval bounds for x/n at critical value z."""
    phat = x / n
    denom = 1 + z * z / n
    center = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    return center - half, center + half


def noninferior(
    fp_new: int, n_new: int, fp_ref: int, n_ref: int, margin: float = 0.005, alpha: float = 0.05
) -> dict[str, Any]:
    """Newcombe hybrid-score non-inferiority test on p_new - p_ref (one-sided upper bound).

    Uses one-sided Wilson score bounds at z = `_norm_ppf(1 - alpha)` (z=1.6449
    for alpha=0.05) combined per Newcombe's method 10; noninferior iff the
    upper bound of the difference is strictly below `margin`.

    Raises `ValueError` if `n_new` or `n_ref` is 0 (undefined proportion),
    rather than letting the division raise `ZeroDivisionError`.
    """
    if n_new == 0 or n_ref == 0:
        raise ValueError("noninferior requires n_new > 0 and n_ref > 0")
    z = _norm_ppf(1 - alpha)
    p_new, p_ref = fp_new / n_new, fp_ref / n_ref
    diff = p_new - p_ref
    _, u_new = _wilson_bounds(fp_new, n_new, z)
    l_ref, _ = _wilson_bounds(fp_ref, n_ref, z)
    upper = diff + math.sqrt((u_new - p_new) ** 2 + (p_ref - l_ref) ** 2)
    return {"diff": diff, "upper": upper, "noninferior": upper < margin}


def _percentile(sorted_values: list[float], pct: float) -> float:
    """Linear-interpolation percentile (R type-7 / numpy default) of pre-sorted values."""
    idx = (len(sorted_values) - 1) * pct / 100
    lo_i, hi_i = math.floor(idx), math.ceil(idx)
    if lo_i == hi_i:
        return sorted_values[int(idx)]
    frac = idx - lo_i
    return sorted_values[lo_i] + (sorted_values[hi_i] - sorted_values[lo_i]) * frac


def cluster_bootstrap(
    rows: Sequence[dict[str, Any]],
    key: str,
    value: Callable[[dict[str, Any]], float],
    n: int = 2000,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Cluster bootstrap: point estimate plus (2.5th, 97.5th) percentile CI.

    Point estimate is the plain mean of `value(row)` over all `rows`. Each of
    the `n` bootstrap replicates resamples cluster keys with replacement
    (`len(groups)` draws from `random.Random(seed)`) and recomputes the
    statistic over every row of the resampled clusters — each selected
    cluster's own rows are used intact, not themselves resampled, since the
    cluster (not the row) is the exchangeable unit. Clusters with identical
    composition (e.g. every cluster contributing the same multiset of
    values) therefore give a degenerate CI (`lo == hi == point`): resampling
    which cluster is picked cannot change the pooled composition. That is
    correct cluster-bootstrap behaviour, not a bug — see
    `test_cluster_bootstrap_identical_clusters_is_degenerate`.

    Raises `ValueError` if `rows` is empty (no point estimate is defined).
    """
    if not rows:
        raise ValueError("cluster_bootstrap requires at least one row")
    groups: dict[Any, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(row[key], []).append(row)
    keys = list(groups.keys())
    point = statistics.fmean(value(row) for row in rows)
    rng = random.Random(seed)
    estimates = []
    for _ in range(n):
        values = [
            value(row)
            for chosen_key in rng.choices(keys, k=len(keys))
            for row in groups[chosen_key]
        ]
        estimates.append(statistics.fmean(values))
    estimates.sort()
    return point, _percentile(estimates, 2.5), _percentile(estimates, 97.5)
