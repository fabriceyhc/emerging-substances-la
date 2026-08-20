"""
Tests for the model-the-curve family (`emerging/analysis/aberration.py`),
the derivative-threshold alternative to EB05 surveyed in
`docs/LITERATURE_REVIEW.md` sec 1.5.

The failure modes worth guarding here are the ones the benchmark writeup
(`docs/findings/benchmark.md`) actually found by running the real data, not
generic crash-proofing:

- `ears_sweep`'s zero-baseline floor is the whole reason a genuinely novel
  substance doesn't divide by zero -- if the floor stopped applying, a
  substance silent in its baseline window would either crash or read `inf`.
- `nb_trend_sweep`'s `own` reading exists specifically to separate "share
  rose" from "count rose" -- the Methamphetamine finding (share/derivative
  reads positive purely from a shrinking county denominator, while the raw
  count is flat or falling) is exactly what `nb_trend_own_z` must be able to
  tell apart from `nb_trend_z`, or the dual gate built on top of it is inert.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging.analysis import aberration
from emerging.analysis.aberration import (_confirmed_flags, _emergent_flags,
                                          _poisson_trend_slope, ears_sweep,
                                          nb_trend_sweep)
from emerging.analysis.trends import EMERGENT_THRESHOLD

RECENT, BASELINE = 2, 4
TOTAL = RECENT + BASELINE


def _counts_denom(counts_by_q: dict[str, list[int]], denom: list[int]):
    """Build a (counts, denom) pair with exactly `TOTAL` quarters, so both
    sweeps produce exactly one as-of row per substance."""
    qs = pd.date_range("2024-01-01", periods=TOTAL, freq="QS")
    rows = []
    for sub, vals in counts_by_q.items():
        assert len(vals) == TOTAL
        for q, n in zip(qs, vals):
            if n:
                rows.append({"substance": sub, "quarter": q, "n": n})
    counts = pd.DataFrame(rows, columns=["substance", "quarter", "n"])
    denom_s = pd.Series(denom, index=qs, name="N")
    return counts, denom_s


def test_ears_flags_a_step_change_against_a_flat_baseline() -> None:
    """A substance flat through the baseline and elevated in the recent
    window must score a large positive z; a substance flat throughout scores
    near zero."""
    counts, denom = _counts_denom(
        {"riser": [5, 5, 5, 5, 30, 30], "flat": [5, 5, 5, 5, 5, 5]},
        denom=[100] * TOTAL,
    )
    sw = ears_sweep(counts, denom, recent_quarters=RECENT, baseline_quarters=BASELINE)
    z = sw.set_index("substance")["ears_z"]
    assert z["riser"] > 3.0
    assert abs(z["flat"]) < 1e-6


def test_ears_floors_the_zero_baseline_variance() -> None:
    """A substance entirely silent in the baseline window has sample SD == 0
    by construction -- the Poisson noise floor must keep the score finite
    (never inf/nan) rather than dividing by zero, since 'silent baseline,
    active recent window' is exactly what a genuinely novel substance looks
    like."""
    counts, denom = _counts_denom(
        {"novel": [0, 0, 0, 0, 4, 4]}, denom=[100] * TOTAL,
    )
    sw = ears_sweep(counts, denom, recent_quarters=RECENT, baseline_quarters=BASELINE)
    z = sw.set_index("substance")["ears_z"]
    assert np.isfinite(z["novel"])
    assert z["novel"] > 0


def test_poisson_trend_slope_sign_and_significance() -> None:
    """A clearly increasing count series must fit a positive slope with a
    Wald z well past the 1.96 line; a flat series must not."""
    offset = np.full(8, 100.0)
    rising = np.array([2, 3, 4, 6, 9, 14, 20, 30], dtype=float)
    flat = np.full(8, 10.0)

    b_rise, se_rise, _ = _poisson_trend_slope(rising, offset)
    b_flat, se_flat, _ = _poisson_trend_slope(flat, offset)

    assert b_rise > 0
    assert b_rise / se_rise > 1.96
    assert abs(b_flat) < 1e-6
    assert abs(b_flat / se_flat) < 1.96


def test_nb_trend_separates_share_drift_from_a_real_count_change() -> None:
    """The Methamphetamine finding in docs/findings/benchmark.md: a
    substance's own count can be perfectly flat while its *share* rises
    purely because the county denominator shrinks. `nb_trend_z` (share) must
    read that as a rise; `nb_trend_own_z` (raw count, no offset) must not --
    otherwise the dual gate built from the two is not testing anything a
    single reading didn't already test."""
    counts, denom = _counts_denom(
        {"flat_count": [20, 20, 20, 20, 20, 20]},
        denom=[400, 400, 400, 400, 200, 200],  # county total halves in the recent window
    )
    sw = nb_trend_sweep(counts, denom, recent_quarters=RECENT, baseline_quarters=BASELINE)
    row = sw.set_index("substance").loc["flat_count"]
    assert row["nb_trend_z"] > 1.96          # share reads as rising
    assert abs(row["nb_trend_own_z"]) < 1.96  # raw count reads as flat


def test_nb_trend_own_and_share_agree_on_a_genuine_rise() -> None:
    """When the county denominator is flat, a real rise in a substance's own
    count must show up as a significant slope on *both* readings -- the dual
    gate should not suppress a substance that is genuinely rising by both
    measures, only one that clears share alone."""
    counts, denom = _counts_denom(
        {"riser": [5, 5, 5, 5, 40, 40]}, denom=[500] * TOTAL,
    )
    sw = nb_trend_sweep(counts, denom, recent_quarters=RECENT, baseline_quarters=BASELINE)
    row = sw.set_index("substance").loc["riser"]
    assert row["nb_trend_z"] > 1.96
    assert row["nb_trend_own_z"] > 1.96


def test_confirmed_flags_requires_truly_consecutive_quarters_not_just_rows() -> None:
    """The bug this guards: rolling over a substance's own *filtered* rows
    (n_recent==0 quarters are missing entirely from a sweep) would treat two
    alarming quarters either side of a real gap as consecutive. Two
    substances with identical hot/cold quarter *labels* but a real missing
    quarter in one of them must not both confirm."""
    qs = pd.date_range("2024-01-01", periods=5, freq="QS")
    sw = pd.DataFrame({
        "substance": ["contig"] * 3 + ["gapped"] * 2,
        "as_of": [qs[0], qs[1], qs[2], qs[0], qs[2]],  # "gapped" skips qs[1]
        "score": [5.0, 5.0, 5.0, 5.0, 5.0],
    })
    flags = _confirmed_flags(sw, "score", threshold=1.96, min_consecutive=3)
    assert flags[sw["substance"] == "contig"].iloc[-1]       # 3 truly consecutive
    assert not flags[sw["substance"] == "gapped"].any()      # only 2, with a gap


def test_confirmed_flags_needs_the_full_run_length() -> None:
    qs = pd.date_range("2024-01-01", periods=4, freq="QS")
    sw = pd.DataFrame({
        "substance": ["x"] * 4, "as_of": qs, "score": [5.0, 5.0, 0.1, 5.0],
    })
    flags = _confirmed_flags(sw, "score", threshold=1.96, min_consecutive=3)
    # Only 2 consecutive hot quarters at a time anywhere in this series.
    assert not flags.any()


def test_confirmed_flags_marks_every_quarter_from_the_kth_onward_not_just_the_kth() -> None:
    """A sustained rise should stay confirmed for as long as it holds, not
    flicker back to unconfirmed -- `score_model`'s episode/lag logic reads
    contiguous alarm spans, so a single confirmed quarter surrounded by
    unconfirmed ones on either side of an otherwise-continuing hot streak
    would misread as multiple short episodes."""
    qs = pd.date_range("2024-01-01", periods=5, freq="QS")
    sw = pd.DataFrame({
        "substance": ["x"] * 5, "as_of": qs, "score": [5.0] * 5,
    })
    flags = _confirmed_flags(sw, "score", threshold=1.96, min_consecutive=3)
    assert list(flags) == [False, False, True, True, True]


def _episode_sweep(rows) -> pd.DataFrame:
    """A minimal (substance, as_of, score, n_baseline) sweep -- enough to
    drive `_emergent_flags` without going through the full `nb_trend_sweep`
    pipeline."""
    return pd.DataFrame(rows, columns=["substance", "as_of", "score", "n_baseline"])


def test_emergent_flags_separates_a_rare_riser_from_an_already_large_one() -> None:
    """The Methamphetamine-vs-Xylazine distinction this function exists for:
    a substance whose baseline was already far above EMERGENT_THRESHOLD stays
    unflagged even while its slope is significant, and a substance whose
    baseline was near zero gets flagged, at the identical alarm z."""
    qs = pd.date_range("2024-01-01", periods=3, freq="QS")
    baseline_established = 4 * EMERGENT_THRESHOLD  # avg/quarter well above the line
    baseline_rare = 0.0
    sw = _episode_sweep([
        ("established", qs[0], 3.0, baseline_established * 4),
        ("established", qs[1], 3.0, baseline_established * 4),
        ("rare", qs[0], 3.0, baseline_rare * 4),
        ("rare", qs[1], 3.0, baseline_rare * 4),
    ])
    flags = _emergent_flags(sw, "score", threshold=1.96, baseline_quarters=4)
    assert not flags[sw["substance"] == "established"].any()
    assert flags[sw["substance"] == "rare"].all()


def test_emergent_flags_holds_the_onset_baseline_for_the_whole_episode() -> None:
    """Ported from `trends.alarm_history`'s own anchoring rule: a substance's
    baseline is read once, at the episode's first alarm quarter, and held for
    the rest of the episode -- not recomputed fresh every quarter, which
    would let a real, sustained rise "grow into" a higher rolling baseline
    and flip from emergent to not partway through its own episode."""
    qs = pd.date_range("2024-01-01", periods=3, freq="QS")
    sw = _episode_sweep([
        # Onset baseline is near zero (emergent); n_baseline creeps up in
        # later rows as if the substance's own growth were feeding back into
        # a naively-refreshed baseline window -- must not matter here.
        ("riser", qs[0], 3.0, 0.0),
        ("riser", qs[1], 3.0, 20 * EMERGENT_THRESHOLD),
        ("riser", qs[2], 3.0, 40 * EMERGENT_THRESHOLD),
    ])
    flags = _emergent_flags(sw, "score", threshold=1.96, baseline_quarters=4)
    assert flags.all()


def test_emergent_flags_false_outside_any_alarm_episode() -> None:
    """`emergent` is an annotation on top of being in alarm, never a
    standalone signal -- a quarter that never breaches the threshold must
    read `False` regardless of how low its baseline is."""
    qs = pd.date_range("2024-01-01", periods=2, freq="QS")
    sw = _episode_sweep([("quiet", qs[0], 0.1, 0.0), ("quiet", qs[1], 0.1, 0.0)])
    flags = _emergent_flags(sw, "score", threshold=1.96, baseline_quarters=4)
    assert not flags.any()


def test_nb_trend_sweep_emergent_column_matches_direct_calculation() -> None:
    """End-to-end: a substance absent through the baseline and alarming in
    the recent window must come back `emergent=True`; a substance whose
    baseline was already far above `EMERGENT_THRESHOLD` must not, even
    though both clear `nb_trend_z`'s own alarm line."""
    counts, denom = _counts_denom(
        {
            "rare_riser": [0, 0, 0, 0, 30, 30],
            # ~25/quarter baseline, an order of magnitude above the line.
            "established": [100, 100, 100, 100, 400, 400],
        },
        denom=[1000] * TOTAL,
    )
    sw = nb_trend_sweep(counts, denom, recent_quarters=RECENT, baseline_quarters=BASELINE)
    assert "n_baseline" in sw.columns and "emergent" in sw.columns
    row = sw.set_index("substance")
    assert row.loc["rare_riser", "nb_trend_z"] > 1.96
    assert bool(row.loc["rare_riser", "emergent"])
    assert row.loc["established", "nb_trend_z"] > 1.96
    assert not bool(row.loc["established", "emergent"])


def test_weighted_recompute_is_used_instead_of_the_passed_counts(monkeypatch) -> None:
    """`weighted=True` must recompute `load_quarterly` fresh per as-of quarter
    (the cross-case-aggregate discipline `sweep_eb05` also follows), not
    silently fall back to slicing the unweighted `counts` table passed in --
    checked by monkeypatching `load_quarterly` to a value distinguishable
    from the fixture's own counts, so the test fails if the weighted branch
    is ever skipped."""
    counts, denom = _counts_denom({"x": [1] * TOTAL}, denom=[10] * TOTAL)
    qs = sorted(denom.index)
    calls = []

    def fake_load_quarterly(mentions_path, cutoff, weighted, quiet, role_discount):
        calls.append((cutoff, weighted, role_discount))
        weighted_counts = pd.DataFrame(
            {"substance": ["x"] * TOTAL, "quarter": qs, "n": [99] * TOTAL})
        return weighted_counts, None

    monkeypatch.setattr(aberration, "load_quarterly", fake_load_quarterly)
    sw = ears_sweep(counts, denom, recent_quarters=RECENT, baseline_quarters=BASELINE,
                    weighted=True, role_discount=True)
    assert len(calls) == 1
    assert calls[0][1] is True and calls[0][2] is True
    # n_recent must reflect the weighted stub's count (99*2), not the
    # unweighted fixture's (1*2) -- proof the weighted table was actually used.
    assert int(sw.set_index("substance").loc["x", "n_recent"]) == 99 * RECENT


def test_sweeps_raise_when_the_window_does_not_fit() -> None:
    counts, denom = _counts_denom({"x": [1] * TOTAL}, denom=[10] * TOTAL)
    with pytest.raises(ValueError):
        ears_sweep(counts, denom.iloc[:-1], recent_quarters=RECENT,
                  baseline_quarters=BASELINE)
    with pytest.raises(ValueError):
        nb_trend_sweep(counts, denom.iloc[:-1], recent_quarters=RECENT,
                       baseline_quarters=BASELINE)
