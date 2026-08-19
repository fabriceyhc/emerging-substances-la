"""
Tests for the spatial-concentration score (`docs/GPS_V2_DESIGN.md` #3).

The load-bearing claim in `core/concentration.py` is the closed-form null
moments: the whole point of `_moments` is that it replaces a permutation test,
so if the algebra is wrong there is nothing left to catch it -- the sweep runs
fine and produces plausible, systematically wrong z-scores. Every other test
here is secondary to `test_moments_match_simulation`, which checks the
derivation against the random-labelling null it claims to describe by actually
drawing from it.

The second failure this pins down is quieter: `sweep_concentration` feeds an
as-of sweep, so a window that reached forward even one quarter would let a
substance's future geography raise its past score. That is the same class of
leak `_ever_independent` had (`trends._position_credit`), and
`test_sweep_is_trailing_only` is the same style of check -- extend the data
and assert the earlier answers do not move.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from emerging.core import concentration as cn


def _ring(n: int, seed: int = 0) -> np.ndarray:
    """n points on a circle: every point has the same neighbour structure, so
    the graph is regular and hand-checkable."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.c_[np.cos(t), np.sin(t)] * 1000.0


# --------------------------------------------------------------------------
# The null moments -- the derivation this module rests on
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n_cases", [5, 10, 30, 100])
def test_moments_match_simulation(n_cases: int) -> None:
    """Exact E[S] and Var[S] against 4,000 draws from the random-labelling null.

    This is the test that makes the closed form usable in place of a
    permutation run. Tolerances are Monte Carlo error at 4,000 replicates
    (~2% on the mean, ~4% on the sd), not slack for an approximation -- the
    formula is exact, and a real algebra error moves these by far more.
    """
    rng = np.random.default_rng(11)
    xy = rng.normal(scale=5000.0, size=(400, 2))
    h, a2, p2 = cn.neighbour_graph(xy, frac=0.05)
    n = len(xy)

    sim = np.empty(4000)
    for i in range(len(sim)):
        u = np.zeros(n)
        u[rng.choice(n, n_cases, replace=False)] = 1.0
        sim[i] = u @ (h @ u)

    mean, var = cn._moments(a2, p2, np.array([float(n_cases)]), float(n))
    assert mean[0] == pytest.approx(sim.mean(), rel=0.05, abs=0.05)
    assert np.sqrt(var[0]) == pytest.approx(sim.std(), rel=0.08)


def test_moments_degenerate_when_every_death_is_a_case() -> None:
    """C == N means the case set *is* the population: S is deterministic."""
    xy = _ring(50)
    h, a2, p2 = cn.neighbour_graph(xy, frac=0.1)
    mean, var = cn._moments(a2, p2, np.array([50.0]), 50.0)
    assert mean[0] == pytest.approx(a2)
    assert var[0] == pytest.approx(0.0, abs=1e-6)


def test_neighbour_graph_is_symmetric_with_empty_diagonal() -> None:
    """`_moments` derives `a22` assuming a symmetric, zero-diagonal kernel;
    an asymmetric one would need a separate mutual-neighbour term."""
    h, a2, p2 = cn.neighbour_graph(_ring(40), frac=0.1)
    assert (abs(h - h.T) > 0).nnz == 0
    assert h.diagonal().sum() == 0
    assert set(np.unique(h.data)) <= {1.0}
    # A regular graph: a2 = N*deg and p2 = N*deg*(deg-1).
    deg = np.asarray(h.sum(axis=1)).ravel()
    assert a2 == pytest.approx(deg.sum())
    assert p2 == pytest.approx((deg * (deg - 1)).sum())


# --------------------------------------------------------------------------
# Withholding: the two regimes where the score must not be read
# --------------------------------------------------------------------------

def _frame(xy: np.ndarray, assign: dict[str, list[int]]) -> tuple:
    bg = pd.DataFrame({"CaseNumber": [f"c{i}" for i in range(len(xy))],
                       "x": xy[:, 0], "y": xy[:, 1]})
    rows = [(f"c{i}", s) for s, idx in assign.items() for i in idx]
    return bg, pd.DataFrame(rows, columns=["CaseNumber", "substance"])


def test_low_count_substances_are_withheld_not_scored() -> None:
    """Below `MIN_CASES` the statistic is too discrete for a z-score: at C = 3
    it can only be 0, 2, 4 or 6, so one coincidence is four sds."""
    rng = np.random.default_rng(3)
    xy = rng.normal(scale=5000.0, size=(300, 2))
    bg, cases = _frame(xy, {"rare": [0, 1, 2], "common": list(range(50))})
    out = cn.window_scores(bg, cases)

    assert out.loc["rare", "withheld"] == "n<min_cases"
    assert out.loc["rare", "z"] == 0.0
    assert out.loc["common", "withheld"] == ""


def test_dominant_substance_is_withheld() -> None:
    """Past `MAX_FRACTION` the random-labelling null draw becomes the case set
    itself, so the contrast degenerates -- the same guard, and the same
    constant, `core/spatial.py` applies to methamphetamine and fentanyl."""
    rng = np.random.default_rng(4)
    xy = rng.normal(scale=5000.0, size=(200, 2))
    bg, cases = _frame(xy, {"everywhere": list(range(150))})
    out = cn.window_scores(bg, cases)

    assert out.loc["everywhere", "withheld"] == "case_fraction>max"
    assert out.loc["everywhere", "z"] == 0.0


def test_clustered_substance_scores_above_dispersed_one() -> None:
    """The statistic's actual job, on data where the answer is known."""
    rng = np.random.default_rng(5)
    far = rng.normal(scale=20000.0, size=(300, 2))
    tight = rng.normal(loc=(60000, 0), scale=200.0, size=(30, 2))
    xy = np.vstack([far, tight])
    bg, cases = _frame(xy, {
        "clustered": list(range(300, 330)),
        "dispersed": list(range(0, 300, 10)),
    })
    out = cn.window_scores(bg, cases)

    assert out.loc["clustered", "z"] > 5
    assert out.loc["dispersed", "z"] < out.loc["clustered", "z"]


# --------------------------------------------------------------------------
# The squash, and as-of safety
# --------------------------------------------------------------------------

def test_concentration_weight_is_bounded_and_monotone() -> None:
    z = np.array([-10.0, -1.0, 0.0, 1.96, 5.0, 50.0])
    w = cn.concentration_weight(z)
    assert (w >= 0).all() and (w <= 1).all()
    assert (np.diff(w) >= 0).all()
    # A dispersed substance is left alone, never penalised.
    assert w[0] == 0.0 and w[2] == 0.0
    assert w[3] == pytest.approx(0.95, abs=0.01)


def test_sweep_is_trailing_only() -> None:
    """Appending later quarters must not change any earlier as-of score.

    The leak this guards against is the one `_ever_independent` actually had:
    a cross-window aggregate that lets the future raise the past. Here the
    window is trailing by construction, and this asserts it stays that way.
    """
    rng = np.random.default_rng(6)
    qs = pd.date_range("2020-01-01", periods=8, freq="QS")
    n = 120
    bg = pd.DataFrame({
        "CaseNumber": [f"c{i}" for i in range(n)],
        "quarter": rng.choice(qs, n),
        "x": rng.normal(scale=5000.0, size=n),
        "y": rng.normal(scale=5000.0, size=n),
    })
    cases = pd.DataFrame({
        "CaseNumber": bg["CaseNumber"], "quarter": bg["quarter"],
        "substance": rng.choice(["a", "b"], n),
    })

    full = cn.sweep_concentration(bg, cases, list(qs))
    early_q = list(qs[:6])
    early = cn.sweep_concentration(
        bg[bg["quarter"].isin(early_q)], cases[cases["quarter"].isin(early_q)],
        early_q)

    key = ["substance", "as_of"]
    merged = early.merge(full, on=key, suffixes=("_early", "_full"))
    assert not merged.empty
    assert np.allclose(merged["z_early"], merged["z_full"])
