"""Exact CIs, McNemar, Holm, non-inferiority, cluster bootstrap.

Reference values verified against R (R 4.5.2) `binom.test` before writing
these assertions; see task-8-report.md for the verification commands and
output. All four brief values matched R to well within the stated
tolerances, so no test values were changed from the brief.
"""

import pytest

from muc_one_span.benchsim.stats import (
    clopper_pearson,
    cluster_bootstrap,
    holm,
    mcnemar_exact,
    noninferior,
)


def test_clopper_pearson_reference() -> None:
    lo, hi = clopper_pearson(59, 59)
    assert lo == pytest.approx(0.9394, abs=1e-4) and hi == 1.0
    assert clopper_pearson(0, 280)[1] == pytest.approx(0.01308, abs=1e-4)
    lo, hi = clopper_pearson(7, 20)
    assert (lo, hi) == pytest.approx((0.1539, 0.5922), abs=1e-4)


def test_mcnemar_exact() -> None:
    assert mcnemar_exact(0, 0) == 1.0
    assert mcnemar_exact(1, 9) == pytest.approx(0.02148, abs=1e-5)


def test_holm_monotone() -> None:
    adj = holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert adj == pytest.approx({"a": 0.03, "b": 0.06, "c": 0.06})


def test_noninferiority() -> None:
    assert noninferior(0, 280, 0, 280)["noninferior"] is False  # CI too wide at n=280
    assert noninferior(0, 2000, 0, 2000)["noninferior"] is True


def test_noninferiority_alpha_widens_or_narrows_bound() -> None:
    # Same data, stricter (smaller) alpha uses a larger z, so a wider upper bound.
    strict = noninferior(2, 2000, 0, 2000, alpha=0.001)
    default = noninferior(2, 2000, 0, 2000, alpha=0.05)
    assert strict["upper"] > default["upper"]
    assert strict["noninferior"] is False
    assert default["noninferior"] is True


def test_cluster_bootstrap_heterogeneous_clusters() -> None:
    # 20 clusters, each internally homogeneous (all-0 or all-1); which clusters
    # get resampled varies the pooled mean, so the CI is non-degenerate.
    rows = [{"g": g, "ok": 1 if g % 2 else 0} for g in range(20) for _ in range(2)]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1)
    assert mean == 0.5 and lo < 0.5 < hi


def test_cluster_bootstrap_identical_clusters_is_degenerate() -> None:
    # Every cluster has the same composition (one 0, one 1), so no matter
    # which clusters a cluster-only bootstrap draws, the pooled mean is
    # always exactly 0.5 -- a point-mass CI is the correct behaviour here.
    rows = [{"g": i // 2, "ok": i % 2} for i in range(40)]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=500, seed=1)
    assert mean == 0.5
    assert lo == hi == 0.5


def test_cluster_bootstrap_single_replicate() -> None:
    rows = [{"g": 0, "ok": 1}, {"g": 1, "ok": 0}]
    mean, lo, hi = cluster_bootstrap(rows, "g", lambda r: r["ok"], n=1, seed=0)
    assert mean == 0.5
    assert lo == hi  # a single replicate has no percentile spread
