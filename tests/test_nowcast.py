"""
Tests for the reporting-delay nowcast.

The dangerous failure here is not a crash, it is a plausible correction
computed from a curve that does not apply. Two guards therefore have to work:
the settled-era drift check, which rejects a vintage pair that differs by case
definition rather than by accrual, and `transfer_check`, which rejects a curve
that the target extract does not obey. Both are tested on synthetic data where
the right answer is known by construction.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging import nowcast as nc


def _write(tmp_path, name, counts: dict[str, int]):
    """A minimal LACME-shaped CSV: one DeathDate row per death."""
    rows = []
    for month, n in counts.items():
        rows += [f"{month}-15"] * n
    p = tmp_path / name
    pd.DataFrame({"DeathDate": rows}).to_csv(p, index=False)
    return p


# Spans the settled era (<= 2022-12) so the drift check has something to read.
TRUE = {f"{y}-{m:02d}": 100 for y in (2021, 2022, 2023) for m in range(1, 13)}
TRUE |= {f"2024-{m:02d}": 100 for m in range(1, 9)}
# Completeness by months elapsed from 2024-08.
PI = {0: 0.30, 1: 0.60, 2: 0.80, 3: 0.90, 4: 0.95}


def _vintages(tmp_path):
    early = {}
    asof = pd.Period("2024-08", "M")
    for m, n in TRUE.items():
        d = (asof - pd.Period(m, "M")).n
        early[m] = int(round(n * PI.get(d, 1.0)))
    return _write(tmp_path, "early.csv", early), _write(tmp_path, "late.csv", TRUE)


def test_pava_is_non_decreasing_and_leaves_monotone_input_alone() -> None:
    y = np.array([0.3, 0.7, 0.65, 0.9, 0.88, 1.0])
    out = nc._pava(y, np.ones_like(y))
    assert np.all(np.diff(out) >= -1e-12)
    already = np.array([0.2, 0.4, 0.6, 0.9])
    assert np.allclose(nc._pava(already, np.ones(4)), already)


def test_pava_preserves_the_weighted_total() -> None:
    y = np.array([0.9, 0.1, 0.5])
    w = np.array([2.0, 3.0, 1.0])
    assert (nc._pava(y, w) * w).sum() == pytest.approx((y * w).sum())


def test_delay_curve_recovers_a_planted_curve(tmp_path) -> None:
    early, late = _vintages(tmp_path)
    curve, meta = nc.delay_curve(early, late, max_delay=8)
    assert meta["settled_ok"], "a clean pair must pass the drift check"
    got = dict(zip(curve["delay_m"], curve["completeness"]))
    for d, p in PI.items():
        assert got[d] == pytest.approx(p, abs=0.02), f"delay {d}"
    assert got[6] == pytest.approx(1.0, abs=0.02)


def test_settled_drift_rejects_a_definitional_change(tmp_path) -> None:
    """A pair whose settled era moves is not measuring delay."""
    early, late = _vintages(tmp_path)
    inflated = {m: int(n * 1.05) for m, n in TRUE.items()}
    late2 = _write(tmp_path, "late2.csv", inflated)
    _, meta = nc.delay_curve(early, late2, max_delay=8)
    assert not meta["settled_ok"]
    assert meta["settled_drift"] < -0.01


def test_delay_curve_refuses_a_pair_with_no_settled_era(tmp_path) -> None:
    recent = {f"2024-{m:02d}": 100 for m in range(1, 9)}
    a = _write(tmp_path, "a.csv", recent)
    b = _write(tmp_path, "b.csv", recent)
    with pytest.raises(ValueError, match="settled-era drift check"):
        nc.delay_curve(a, b, max_delay=8)


def test_apply_curve_inverts_the_thinning(tmp_path) -> None:
    early, late = _vintages(tmp_path)
    curve, _ = nc.delay_curve(early, late, max_delay=8)
    obs, _ = nc._monthly(early)
    est = nc.apply_curve(obs, pd.Period("2024-08", "M"), curve)
    recent = est[est["delay_m"] <= 4]
    assert np.allclose(recent["estimated"], 100, atol=6)


def test_transfer_check_accepts_an_extract_that_obeys_the_curve(tmp_path) -> None:
    early, late = _vintages(tmp_path)
    curve, _ = nc.delay_curve(early, late, max_delay=8)
    obs, _ = nc._monthly(early)          # thinned by exactly this curve
    chk = nc.transfer_check(obs, pd.Period("2024-08", "M"), curve,
                            baseline_lo=6, baseline_hi=30)
    near = chk[chk["delay_m"] <= 4]
    assert ((near["ratio"] > 0.65) & (near["ratio"] < 1.35)).all()


def test_transfer_check_rejects_a_cliff(tmp_path) -> None:
    """The real case: full strength, then a pull boundary.

    This is the guard that stopped a wrong 23% upward correction to 2025Q4.
    """
    early, late = _vintages(tmp_path)
    curve, _ = nc.delay_curve(early, late, max_delay=8)
    cliff = dict(TRUE)
    cliff["2024-07"] = 25          # would be 60 under the curve
    cliff["2024-08"] = 5           # would be 30
    obs, _ = nc._monthly(_write(tmp_path, "cliff.csv", cliff))
    chk = nc.transfer_check(obs, pd.Period("2024-08", "M"), curve,
                            baseline_lo=6, baseline_hi=30)
    near = chk[chk["delay_m"] <= 4]
    bad = ((near["ratio"] > 1.35) | (near["ratio"] < 0.65)).sum()
    assert bad >= 2, "a cliff must be flagged as not obeying the curve"
