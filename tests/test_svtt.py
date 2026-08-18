"""
Tests for the spatial-variation-in-temporal-trends scan.

Two things have to be right and neither is visible in the output. The trend
MLE is a hand-rolled Newton solve, so it is checked against the score equation
it claims to solve and against data with a known planted slope. The p-values
come from a Monte Carlo maximum, so they are checked for calibration on data
generated from H0 — one common trend, arbitrary spatial variation in level,
which is the null that matters here and the one a level-based scan would fail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging.analysis import svtt


def _t(T: int) -> np.ndarray:
    t = np.arange(T, dtype=float)
    return t - t.mean()


def test_fit_beta_satisfies_its_own_score_equation() -> None:
    """The MLE must equate population-weighted and case-weighted mean time."""
    rng = np.random.default_rng(0)
    T = 6
    t = _t(T)
    P = rng.uniform(5_000, 50_000, size=(40, T))
    true_b = rng.uniform(-0.3, 0.3, size=40)
    lam = P * np.exp(true_b[:, None] * t[None, :]) * 1e-3
    c = rng.poisson(lam).astype(float)

    C, TC = c.sum(1), (c * t).sum(1)
    b = svtt._fit_beta(C, TC, P, t)

    w = P * np.exp(b[:, None] * t[None, :])
    lhs = (w * t).sum(1) / w.sum(1)
    rhs = TC / C
    assert np.allclose(lhs, rhs, atol=1e-8), "Newton did not reach the root"


def test_fit_beta_recovers_a_planted_trend() -> None:
    rng = np.random.default_rng(1)
    T = 8
    t = _t(T)
    P = np.full((1, T), 200_000.0)
    b_true = 0.25
    lam = P * np.exp(b_true * t[None, :]) * 2e-4
    c = rng.poisson(lam).astype(float)
    b = svtt._fit_beta(c.sum(1), (c * t).sum(1), P, t)[0]
    assert abs(b - b_true) < 0.05


def test_llr_is_zero_when_the_trends_agree() -> None:
    """Same slope inside and outside, wildly different levels."""
    T = 6
    t = _t(T)
    P_in = np.full((1, T), 50_000.0)
    p_all = np.full(T, 1_000_000.0)
    b = 0.2
    # Inside is 10x the rate of outside, but on the identical slope.
    c_in = 50 * np.exp(b * t)
    c_all = c_in + 200 * np.exp(b * t)
    llr, b_i, b_o = svtt.svtt_llr(c_in[None, :], P_in, c_all, p_all, t)
    assert b_i[0] == pytest.approx(b_o[0], abs=1e-6)
    assert llr[0] == pytest.approx(0.0, abs=1e-6)


def test_llr_fires_when_only_the_slope_differs() -> None:
    T = 6
    t = _t(T)
    P_in = np.full((1, T), 200_000.0)
    p_all = np.full(T, 2_000_000.0)
    c_in = 100 * np.exp(0.35 * t)
    c_all = c_in + 900 * np.exp(-0.10 * t)
    llr, b_i, b_o = svtt.svtt_llr(c_in[None, :], P_in, c_all, p_all, t)
    assert b_i[0] > b_o[0]
    assert llr[0] > 20


def test_llr_is_one_sided() -> None:
    """A circle falling faster than the county is not a signal."""
    T = 6
    t = _t(T)
    P_in = np.full((1, T), 200_000.0)
    p_all = np.full(T, 2_000_000.0)
    c_in = 100 * np.exp(-0.35 * t)
    c_all = c_in + 900 * np.exp(0.10 * t)
    llr, b_i, b_o = svtt.svtt_llr(c_in[None, :], P_in, c_all, p_all, t)
    assert b_i[0] < b_o[0]
    assert llr[0] == 0.0


def _toy_circles(n_zip: int = 16, seed: int = 0):
    rng = np.random.default_rng(seed)
    side = int(np.sqrt(n_zip))
    cen = pd.DataFrame(
        {"x": [(z % side) * 5000.0 for z in range(n_zip)],
         "y": [(z // side) * 5000.0 for z in range(n_zip)]},
        index=[f"z{z:02d}" for z in range(n_zip)])
    pop = rng.uniform(20_000, 60_000, size=n_zip)
    return cen, pop, svtt.Circles(cen, pop, max_fraction=0.25)


def test_circles_respect_the_population_cap() -> None:
    cen, pop, circ = _toy_circles()
    biggest = np.cumsum(pop[circ.order], axis=1)[
        np.arange(len(pop)), circ.k_len - 1]
    assert (biggest <= 0.25 * pop.sum() + 1e-9).all()


def test_planted_cluster_is_recovered() -> None:
    """A zip whose rate climbs while the rest falls must be found."""
    cen, pop, circ = _toy_circles(seed=2)
    T = 5
    t = _t(T)
    P = np.repeat(pop[:, None], T, axis=1)
    rng = np.random.default_rng(4)
    lam = P * 6e-4 * np.exp(-0.15 * t)[None, :]
    lam[0] = P[0] * 6e-4 * np.exp(0.45 * t)      # one zip going the other way
    lam[0] *= 6                                   # enough cases to be scoreable
    cases = rng.poisson(lam).astype(float)

    llr, b_i, b_o = svtt.svtt_llr(circ.accumulate(cases), circ.accumulate(P),
                                  cases.sum(0), P.sum(0), t)
    llr = np.where(circ.valid.ravel(), llr, 0.0)
    ci, ki = circ.labels()
    best = int(llr.argmax())
    assert llr[best] > 0
    assert "z00" in {circ.zips[m] for m in circ.order[ci[best], :ki[best] + 1]}


def test_null_calibration() -> None:
    """One common trend, spatial variation in level only -> nominal rejection.

    The level variation is deliberately extreme (a 5x spread across zips).
    A level-based scan would light up on this; SVTT must not, because nothing
    about the *slope* varies.
    """
    cen, pop, circ = _toy_circles(seed=7)
    T, alpha, n_rep, n_perm = 5, 0.20, 60, 99
    t = _t(T)
    P = np.repeat(pop[:, None], T, axis=1)
    p_all = P.sum(0)
    valid = circ.valid.ravel()
    p_cyl = circ.accumulate(P)

    b_g = -0.12
    level = np.linspace(1.0, 5.0, len(pop))          # the confound
    q = (P * level[:, None] * np.exp(b_g * t)[None, :])
    q = (q / q.sum()).ravel()
    total = 3000
    rng = np.random.default_rng(21)

    rejects = 0
    for r in range(n_rep):
        cases = rng.multinomial(total, q).reshape(P.shape).astype(float)
        llr, _, _ = svtt.svtt_llr(circ.accumulate(cases), p_cyl,
                                  cases.sum(0), p_all, t)
        obs = np.where(valid, llr, 0.0).max()
        sim_rng = np.random.default_rng(5000 + r)
        null = np.empty(n_perm)
        for b in range(n_perm):
            sim = sim_rng.multinomial(total, q).reshape(P.shape).astype(float)
            s, _, _ = svtt.svtt_llr(circ.accumulate(sim), p_cyl,
                                    sim.sum(0), p_all, t)
            null[b] = np.where(valid, s, 0.0).max()
        rejects += ((1 + (null >= obs).sum()) / (n_perm + 1)) <= alpha

    rate = rejects / n_rep
    # Binomial sd at p=0.2, n=60 is 0.052; this is roughly +/-3.5 sd.
    assert 0.02 < rate < 0.38, (
        f"rejection rate {rate:.3f} at nominal {alpha}: SVTT is miscalibrated, "
        f"or it is responding to variation in level rather than in trend")
