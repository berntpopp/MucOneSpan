"""Exact CIs, McNemar, Holm, non-inferiority, cluster bootstrap.

Reference values verified against R (R 4.5.2) `binom.test` before writing
these assertions; see task-8-report.md for the verification commands and
output. All four brief values matched R to well within the stated
tolerances, so no test values were changed from the brief.
"""

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG
from muc_one_span.benchsim.stats import (
    clopper_pearson,
    cluster_bootstrap,
    holm,
    mcnemar_exact,
    noninferior,
)

R = DEFAULT_BENCH_CONFIG.report
A, M = R.alpha, R.ni_margin


def test_clopper_pearson_reference() -> None:
    lo, hi = clopper_pearson(59, 59, alpha=A)
    assert lo == pytest.approx(0.9394, abs=1e-4) and hi == 1.0
    assert clopper_pearson(0, 280, alpha=A)[1] == pytest.approx(0.01308, abs=1e-4)
    lo, hi = clopper_pearson(7, 20, alpha=A)
    assert (lo, hi) == pytest.approx((0.1539, 0.5922), abs=1e-4)


def test_mcnemar_exact() -> None:
    assert mcnemar_exact(0, 0) == 1.0
    assert mcnemar_exact(1, 9) == pytest.approx(0.02148, abs=1e-5)


def test_mcnemar_exact_at_test_split_scale_no_overflow() -> None:
    # n = b + c = 2400 (the plan's `test` split size). The old
    # math.comb(n, i) * 0.5**n implementation raised OverflowError above
    # n ~= 1030 (comb(n, n//2) too large to convert to float); this must
    # work at benchmark scale. Reference: R binom.test(1150, 2400)$p.value
    # and binom.test(1200, 2400)$p.value.
    assert mcnemar_exact(1150, 1250) == pytest.approx(0.04327505, abs=1e-6)
    assert mcnemar_exact(1200, 1200) == 1.0


def test_holm_monotone() -> None:
    adj = holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert adj == pytest.approx({"a": 0.03, "b": 0.06, "c": 0.06})


def test_holm_caps_at_one() -> None:
    # Raw Holm-adjusted p-values (2 * 0.8 = 1.6, 1 * 0.9 = 0.9) both exceed
    # or approach 1.0 once the running max is applied; both must be capped.
    adj = holm({"a": 0.9, "b": 0.8})
    assert adj == {"a": 1.0, "b": 1.0}


def test_noninferiority() -> None:
    assert (
        noninferior(0, 280, 0, 280, margin=M, alpha=A)["noninferior"] is False
    )  # CI too wide at n=280
    assert noninferior(0, 2000, 0, 2000, margin=M, alpha=A)["noninferior"] is True


def test_noninferiority_rejects_zero_n() -> None:
    with pytest.raises(ValueError, match="n_new"):
        noninferior(0, 0, 0, 280, margin=M, alpha=A)
    with pytest.raises(ValueError, match="n_new"):
        noninferior(0, 280, 0, 0, margin=M, alpha=A)


def test_noninferiority_alpha_widens_or_narrows_bound() -> None:
    # Same data, stricter (smaller) alpha uses a larger z, so a wider upper bound.
    strict = noninferior(2, 2000, 0, 2000, margin=M, alpha=A / 50)
    default = noninferior(2, 2000, 0, 2000, margin=M, alpha=A)
    assert strict["upper"] > default["upper"]
    assert strict["noninferior"] is False
    assert default["noninferior"] is True


def test_cluster_bootstrap_heterogeneous_clusters() -> None:
    # 20 clusters, each internally homogeneous (all-0 or all-1); which clusters
    # get resampled varies the pooled mean, so the CI is non-degenerate.
    rows = [{"g": g, "ok": 1 if g % 2 else 0} for g in range(20) for _ in range(2)]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1, alpha=A)
    assert mean == 0.5 and lo < 0.5 < hi


def test_cluster_bootstrap_identical_clusters_is_degenerate() -> None:
    # Every cluster has the same composition (one 0, one 1), so no matter
    # which clusters a cluster-only bootstrap draws, the pooled mean is
    # always exactly 0.5 -- a point-mass CI is the correct behaviour here.
    rows = [{"g": i // 2, "ok": i % 2} for i in range(40)]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1, alpha=A)
    assert mean == 0.5
    assert lo == hi == 0.5


def test_cluster_bootstrap_single_replicate() -> None:
    rows = [{"g": 0, "ok": 1}, {"g": 1, "ok": 0}]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=1, seed=0, alpha=A)
    assert mean == 0.5
    assert lo == hi  # a single replicate has no percentile spread


def test_cluster_bootstrap_rejects_empty_rows() -> None:
    with pytest.raises(ValueError, match="at least one row"):
        cluster_bootstrap([], "g", lambda r: r["ok"], n=1, seed=0, alpha=A)


def test_cluster_bootstrap_interval_follows_alpha() -> None:
    rows = [{"g": g, "ok": 1 if g % 2 else 0} for g in range(20) for _ in range(2)]
    _, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1, alpha=A)
    _, lo_w, hi_w = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1, alpha=A / 5)
    assert lo_w <= lo and hi <= hi_w and (lo_w, hi_w) != (lo, hi)
