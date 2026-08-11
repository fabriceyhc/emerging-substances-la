"""
Tests for the Kelsall-Diggle relative-risk surface.

The estimator has three places it can be quietly wrong: the pseudo-count
(which must pull toward *no difference*, not toward zero risk), the
low-density mask (without which the log-ratio is unbounded and the colour
scale is set by empty desert), and the tolerance p-values. All three are
checked here against constructions where the answer is known.

The grid/geometry class is not tested — it reads the county boundary file and
is exercised end to end by `relrisk surface`.
"""

from __future__ import annotations

import numpy as np
import pytest

from emerging import relrisk as rr


def _grids(shape=(60, 60), rate=0.2, p=0.25, seed=0):
    """Uniform deaths, a constant fraction of which are cases."""
    rng = np.random.default_rng(seed)
    allg = rng.poisson(rate * 40, size=shape).astype(float)
    cases = rng.binomial(allg.astype(int), p).astype(float)
    return cases, allg


def test_log_rr_is_zero_when_the_fraction_is_constant() -> None:
    """The estimand is a ratio against the county rate, so a uniform field
    must come out flat at 0 regardless of how dense it is."""
    allg = np.full((40, 40), 10.0)
    cases = allg * 0.3
    r = rr.log_rr(cases, allg, sigma_cells=3.0)
    assert np.nanmax(np.abs(r)) < 1e-9


def test_pseudo_count_is_centred_on_no_difference() -> None:
    """A cell at exactly the global rate must read 0 at any local density.

    This is the property that separates a sane pseudo-count from one that
    drags sparse areas toward zero risk.
    """
    for density in (1.0, 10.0, 500.0):
        allg = np.full((30, 30), density)
        cases = allg * 0.05
        r = rr.log_rr(cases, allg, sigma_cells=2.0)
        assert np.nanmax(np.abs(r)) < 1e-9, f"density {density}"


def test_zero_cases_gives_a_finite_negative_value() -> None:
    """Not -27 from an epsilon floor, which is what the first version did."""
    allg = np.full((30, 30), 20.0)
    cases = np.zeros_like(allg)
    cases[15, 15] = 1.0          # keep the global rate non-zero
    r = rr.log_rr(cases, allg, sigma_cells=2.0)
    far = r[5, 5]
    assert np.isfinite(far) and far < 0
    assert far > -6, "an empty area should not dominate the colour scale"


def test_low_density_areas_are_masked() -> None:
    allg = np.zeros((40, 40))
    allg[20:24, 20:24] = 50.0    # all the deaths in one block
    cases = allg * 0.2
    r = rr.log_rr(cases, allg, sigma_cells=1.5, min_local=5.0)
    assert np.isnan(r[0, 0]), "empty corner must be masked, not estimated"
    assert np.isfinite(r[22, 22])


def test_planted_hotspot_is_recovered() -> None:
    cases, allg = _grids(p=0.10, seed=3)
    cases[28:33, 28:33] += allg[28:33, 28:33] * 0.6      # local excess
    r = rr.log_rr(cases, allg, sigma_cells=2.5)
    peak = np.unravel_index(np.nanargmax(r), r.shape)
    assert 26 <= peak[0] <= 34 and 26 <= peak[1] <= 34
    assert np.nanmax(r) > 0.5


def test_lcv_prefers_the_smoothness_actually_present() -> None:
    """A broad, gentle gradient must select a larger bandwidth than a tight
    hotspot does."""
    rng = np.random.default_rng(11)
    n = 1200
    xy = rng.uniform(0, 20_000, size=(n, 2))

    broad = 0.05 + 0.30 * (xy[:, 0] / 20_000)            # gradient over 20 km
    tight = np.where(np.linalg.norm(xy - 10_000, axis=1) < 1500, 0.6, 0.05)
    h_broad, _ = rr.lcv_bandwidth_km(xy, rng.random(n) < broad,
                                     candidates=(0.5, 1.0, 2.0, 4.0, 8.0))
    h_tight, _ = rr.lcv_bandwidth_km(xy, rng.random(n) < tight,
                                     candidates=(0.5, 1.0, 2.0, 4.0, 8.0))
    assert h_broad > h_tight, (
        f"CV chose {h_broad} km for a 20 km gradient and {h_tight} km for a "
        f"1.5 km hotspot; it is not responding to the scale of the signal")


def test_tolerance_is_calibrated_pointwise() -> None:
    """Under random labelling with no real excess, about 5% of cells should
    fall below p = 0.05. These are pointwise tests, so that is the target."""
    rng = np.random.default_rng(5)
    shape = (40, 40)
    allg = rng.poisson(6, size=shape).astype(float)
    n_all = int(allg.sum())
    idx = np.repeat(np.arange(shape[0] * shape[1]), allg.astype(int).ravel())

    class _G:
        def __init__(self):
            self.shape = shape
        def bin(self, i):
            # Deliberately integer: log_rr must coerce, because
            # gaussian_filter otherwise returns a truncated integer surface.
            return np.bincount(i, minlength=shape[0] * shape[1]).reshape(shape)

    g = _G()
    n_case = n_all // 5
    pick = rng.choice(n_all, size=n_case, replace=False)
    obs = rr.log_rr(g.bin(idx[pick]).astype(float), allg, 2.0)
    p = rr.tolerance(idx[pick], idx, g, 2.0, 199, rng, obs)
    ok = np.isfinite(obs)
    rate = float((p[ok] <= 0.05).mean())
    assert 0.01 < rate < 0.15, (
        f"{rate:.3f} of cells significant at a nominal 0.05 under the null")
