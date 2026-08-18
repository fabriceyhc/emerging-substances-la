"""
Tests for the space-time scan.

Same priority as `test_treescan.py`: the output is p-values, so calibration is
the property that matters and a miscalibrated scan still returns a plausible
table. `test_null_calibration` is why this file exists. The rest pin the two
likelihood ratios against hand computation and guard the one optimisation that
could silently diverge — `accumulate_batch`, which exists only for speed and
must agree with the straightforward `accumulate` it replaced.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging.analysis import spacetime as st


def _toy(n_zip: int = 16, n_q: int = 12, per_cell: int = 4, seed: int = 0):
    """A small square grid of zips with a uniform scatter of deaths."""
    rng = np.random.default_rng(seed)
    side = int(np.sqrt(n_zip))
    rows = []
    quarters = list(pd.date_range("2023-01-01", periods=n_q, freq="QS"))
    for z in range(n_zip):
        cx, cy = (z % side) * 5000.0, (z // side) * 5000.0
        for q in quarters:
            for _ in range(per_cell):
                rows.append({"zip": f"z{z:02d}", "quarter": q,
                             "x": cx + rng.normal(0, 200),
                             "y": cy + rng.normal(0, 200)})
    bg = pd.DataFrame(rows)
    geo = st.Geometry(bg, quarters, max_fraction=0.25, max_window=4)
    bg["_zi"] = bg["zip"].map({z: i for i, z in enumerate(geo.zips)})
    bg["_qi"] = bg["quarter"].map({q: i for i, q in enumerate(quarters)})
    return bg.reset_index(drop=True), geo, quarters


def test_bernoulli_llr_matches_hand_computation() -> None:
    # 40 cases of 100 deaths inside; 60 of 900 outside; 100 of 1000 overall.
    c, n, C, N = np.array([40.0]), np.array([100.0]), 100.0, 1000.0
    expected = (40 * np.log(0.4) + 60 * np.log(0.6)
                + 60 * np.log(60 / 900) + 840 * np.log(840 / 900)
                - (100 * np.log(0.1) + 900 * np.log(0.9)))
    assert st._bernoulli_llr(c, n, C, N) == pytest.approx(expected)


def test_poisson_llr_matches_hand_computation() -> None:
    c, e, C = np.array([20.0]), np.array([5.0]), 100.0
    expected = 20 * np.log(20 / 5) + 80 * np.log(80 / 95)
    assert st._poisson_llr(c, e, C) == pytest.approx(expected)


@pytest.mark.parametrize("fn,args", [
    (st._bernoulli_llr, (np.array([5.0]), np.array([100.0]), 100.0, 1000.0)),
    (st._poisson_llr, (np.array([2.0]), np.array([5.0]), 100.0)),
])
def test_llr_is_zero_for_a_deficit(fn, args) -> None:
    """Both scans are high-rate only; a cold spot is not a signal."""
    assert fn(*args) == pytest.approx([0.0])


def test_trailing_windows_are_cumulative_from_the_end() -> None:
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert st._trailing(x, 3) == pytest.approx([4.0, 7.0, 9.0])


def test_batched_accumulate_matches_the_loop() -> None:
    """`accumulate_batch` is a speed optimisation and nothing else."""
    _, geo, quarters = _toy()
    rng = np.random.default_rng(3)
    grids = rng.poisson(1.0, size=(5, len(geo.zips), len(quarters))).astype(float)
    one = np.stack([geo.accumulate(g) for g in grids])
    assert np.allclose(one, geo.accumulate_batch(grids))


def test_circles_respect_the_size_cap() -> None:
    bg, geo, _ = _toy()
    biggest = geo.n_cyl[np.arange(len(geo.zips)), geo.k_len - 1, -1]
    assert (biggest <= st.MAX_SPATIAL_FRACTION * geo.N + 1e-9).all()
    assert geo.valid.sum() == int(geo.k_len.sum())


def test_planted_cluster_is_recovered() -> None:
    """A real excess in a known place and window must be found there."""
    bg, geo, quarters = _toy(per_cell=6)
    grid = np.zeros((len(geo.zips), len(quarters)))
    # Baseline everywhere, then a heavy excess in zip 0 over the last 2q.
    grid[:, :] = 1.0
    grid[0, -2:] = 30.0
    C = float(grid.sum())
    res = st.scan_substance(geo, grid, C, 199, np.random.default_rng(1),
                            bg, bg.index.to_numpy(), model="interaction")
    assert geo.zips[0] in res["zips"].split("|")
    assert res["window_q"] == 2
    assert res["p_value"] <= 0.05


@pytest.mark.parametrize("model", ["bernoulli", "interaction"])
def test_null_calibration(model: str) -> None:
    """On data drawn from the null, rejection must sit near the nominal rate.

    The two models have different nulls — random labelling against the
    background for `bernoulli`, permutation of case times for `interaction` —
    so each has to be checked on data generated its own way.
    """
    bg, geo, quarters = _toy(per_cell=6)
    all_idx = bg.index.to_numpy()
    cell = bg["_zi"].to_numpy() * len(quarters) + bg["_qi"].to_numpy()
    n_cells = len(geo.zips) * len(quarters)
    rng = np.random.default_rng(17)
    C, alpha, n_rep = 120, 0.20, 60

    rejects = 0
    for r in range(n_rep):
        if model == "bernoulli":
            pick = rng.choice(all_idx, size=C, replace=False)
            grid = np.bincount(cell[pick], minlength=n_cells).reshape(
                len(geo.zips), len(quarters)).astype(float)
        else:
            zi = rng.integers(0, len(geo.zips), size=C)
            qi = rng.integers(0, len(quarters), size=C)
            grid = np.zeros((len(geo.zips), len(quarters)))
            np.add.at(grid, (zi, qi), 1.0)
        res = st.scan_substance(geo, grid, float(grid.sum()), 99,
                                np.random.default_rng(900 + r), bg, all_idx,
                                model=model)
        rejects += res["p_value"] <= alpha

    rate = rejects / n_rep
    # Binomial sd at p=0.2, n=60 is 0.052; this band is about +/-3.5 sd.
    assert 0.02 < rate < 0.38, (
        f"[{model}] rejection rate {rate:.3f} at nominal {alpha}: the scan is "
        f"miscalibrated and every p-value it reports is wrong")
