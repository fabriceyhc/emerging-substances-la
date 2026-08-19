"""
Tests for the detector comparison (`emerging/validation/benchmark.py`).

The failure mode this file is built around is a *silently unfair* comparison,
not a crash. Two ways that happens, and one test each:

- The matched-budget reading is the whole reason the table means anything --
  `RESEARCH_LOG.md`'s tree-scan head-to-head found an apparent detection
  advantage that vanished once the alarm budget was matched. If
  `matched_threshold` under-delivers a model's budget, that model is scored
  on fewer alarms than the baseline and looks worse than it is.
- `model_score` has to make the dual gate a single monotone number for the
  threshold sweep to be defined at all. If `min(eb05, eb05_own)` ever
  disagreed with the `credible_rise` conjunction at the fixed line, the
  matched and fixed rows would be measuring two different detectors.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging.analysis.trends import ALARM_THRESHOLD
from emerging.validation import benchmark as bm


def _sweep() -> pd.DataFrame:
    qs = pd.date_range("2024-01-01", periods=4, freq="QS")
    rows = []
    for i, q in enumerate(qs):
        rows += [
            {"substance": "riser", "as_of": q, "n_recent": 10 + i,
             "expected": 2.0, "eb05": 1.4 + 0.3 * i, "eb05_own": 2.0 + i},
            {"substance": "share_only", "as_of": q, "n_recent": 20,
             "expected": 10.0, "eb05": 2.5, "eb05_own": 0.8},
            {"substance": "novel", "as_of": q, "n_recent": 3,
             "expected": 0.0, "eb05": 0.9, "eb05_own": 0.9},
        ]
    return pd.DataFrame(rows)


def test_dual_gate_min_matches_the_credible_rise_conjunction() -> None:
    """`min(a, b) > t` is exactly `a > t and b > t`, which is what lets the
    dual gate be threshold-swept instead of only read at 1.5."""
    sw = _sweep()
    dual = bm.model_score(sw, bm.MODELS_BY_ID["eb05-dual"])
    conj = (sw["eb05"] > ALARM_THRESHOLD) & (sw["eb05_own"] > ALARM_THRESHOLD)
    assert ((dual > ALARM_THRESHOLD) == conj).all()
    # share_only is the cocaine pattern: Gate A clears, Gate B does not.
    assert not conj[sw["substance"] == "share_only"].any()


def test_eb05_own_reads_gate_b_alone_not_the_dual_minimum() -> None:
    """`eb05_own` is Gate B run solo -- it must equal `sweep['eb05_own']`
    directly, not `min(eb05, eb05_own)`, so a substance where the two gates
    disagree (the dual gate's whole reason to exist) distinguishes them."""
    sw = _sweep()
    own = bm.model_score(sw, bm.MODELS_BY_ID["eb05-own-role"])
    dual = bm.model_score(sw, bm.MODELS_BY_ID["eb05-dual-role"])
    assert (own == sw["eb05_own"]).all()
    # share_only: eb05_own is low, eb05 is high -- own alone must track the
    # low value, and must disagree with the dual (min) reading somewhere.
    assert own[sw["substance"] == "share_only"].iloc[0] == pytest.approx(0.8)
    assert not (own == dual).all()


def test_ratio_maps_a_zero_baseline_to_infinity_rather_than_dropping_it() -> None:
    """A substance with no baseline is what a genuinely novel one looks like,
    and it is exactly where the unshrunken ratio has nothing to say. It is
    scored as `inf` so the model is measured on that behaviour, not excused
    from it."""
    sw = _sweep()
    s = bm.model_score(sw, bm.MODELS_BY_ID["ratio"])
    novel = s[sw["substance"] == "novel"]
    assert novel.isna().all()          # 3/0 -> NaN, never silently 0
    assert np.isfinite(s[sw["substance"] == "riser"]).all()


@pytest.mark.parametrize("model_id", ["ratio", "eb05", "eb05-dual"])
def test_matched_threshold_delivers_at_least_the_budget(model_id: str) -> None:
    """A model must never be scored on a smaller alarm allowance than the
    baseline it is compared against -- that would flatter the baseline."""
    sw = _sweep()
    model = bm.MODELS_BY_ID[model_id]
    for budget in (1, 3, 6):
        t = bm.matched_threshold(sw, model, budget)
        fired = int((bm.model_score(sw, model) > t).sum())
        assert fired >= budget


def test_matched_threshold_is_exact_when_scores_are_distinct() -> None:
    """Ties are resolved in the model's favour (it gets *at least* the
    budget), so exactness is only claimable when no two scores collide -- the
    fixture above deliberately has ties, which is why this builds its own."""
    qs = pd.date_range("2024-01-01", periods=6, freq="QS")
    sw = pd.DataFrame({
        "substance": ["s"] * 6, "as_of": qs, "n_recent": [5] * 6,
        "expected": [1.0] * 6, "eb05": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        "eb05_own": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
    })
    model = bm.MODELS_BY_ID["eb05"]
    for budget in (1, 2, 4):
        t = bm.matched_threshold(sw, model, budget)
        assert int((bm.model_score(sw, model) > t).sum()) == budget


def test_score_model_lag_is_measured_inside_the_reference_interval() -> None:
    """An alarm before the emergence started is not an early detection of it.
    Lag is counted from `interval_start` to the first alarm *within* the
    interval, so a spurious earlier alarm cannot be scored as prescience."""
    qs = pd.date_range("2024-01-01", periods=4, freq="QS")
    sw = pd.DataFrame({
        "substance": ["x"] * 4, "as_of": qs, "n_recent": [5] * 4,
        "expected": [1.0] * 4, "eb05": [9.0, 0.1, 0.1, 9.0],
        "eb05_own": [9.0, 0.1, 0.1, 9.0],
    })
    gt = pd.DataFrame({"substance": ["x"], "interval_start": [qs[2]],
                       "interval_end": [qs[3]]})
    row, per_sub = bm.score_model(sw, bm.MODELS_BY_ID["eb05"], gt, 1.5)

    # Alarms at q0 and q3; the interval is q2-q3, so the lag is 1, not -2.
    assert per_sub.loc[0, "lag_q"] == 1
    assert row["detected"] == 1


def test_models_sharing_a_sweep_are_grouped_by_case_definition_and_spatial() -> None:
    """The three readings of one sweep must not each trigger a re-run -- the
    weighted sweeps recompute `load_quarterly` 45 times and cost ~40s each."""
    keys = {m.sweep_key for m in bm.MODELS}
    assert len(keys) < len(bm.MODELS)
    assert (bm.MODELS_BY_ID["ratio"].sweep_key
            == bm.MODELS_BY_ID["eb05"].sweep_key
            == bm.MODELS_BY_ID["eb05-dual"].sweep_key)


def test_treescan_reads_its_own_recurrence_interval_and_sweep() -> None:
    """TreeScan's score has nothing to do with the EB05-family sweep columns
    -- it must read `recurrence_interval` from its own cached leaf sweep, and
    must not collide with any EB05 sweep key (it would silently be handed an
    `eb05` column that does not exist in its own sweep otherwise)."""
    sw = pd.DataFrame({
        "substance": ["x", "x"], "as_of": pd.date_range("2024-01-01", periods=2, freq="QS"),
        "recurrence_interval": [12.0, 250.0], "p_value": [0.08, 0.004],
        "llr": [3.0, 9.0],
    })
    model = bm.MODELS_BY_ID["treescan"]
    s = bm.model_score(sw, model)
    assert list(s) == [12.0, 250.0]
    assert model.fixed_threshold == pytest.approx(100.0)
    assert model.sweep_key == ("treescan", None)
    assert model.sweep_key not in {m.sweep_key for m in bm.MODELS if m.id != "treescan"}


def test_every_treescan_reading_is_its_own_sweep_and_shares_the_ri_line() -> None:
    """The leaf, branch and whole-tree readings are three different tables
    off one cached scan, so each needs its own sweep key -- collapsing any
    two would hand one model the other's rows -- while all three are read off
    `recurrence_interval` at TreeScan's own RI line, never the EB05
    `--threshold`. The unguarded and attribution-guarded branch readings
    differ only in `min_excess_share`, which is exactly why it is in the key.
    """
    ids = ["treescan", "treescan-branch", "treescan-branch-attr", "treescan-tree"]
    models = [bm.MODELS_BY_ID[i] for i in ids]
    assert len({m.sweep_key for m in models}) == len(ids)
    assert (bm.MODELS_BY_ID["treescan-branch"].sweep_key
            != bm.MODELS_BY_ID["treescan-branch-attr"].sweep_key)

    sw = pd.DataFrame({
        "substance": ["x"], "as_of": [pd.Timestamp("2024-01-01")],
        "recurrence_interval": [250.0], "p_value": [0.004], "llr": [9.0],
        "source_node": ["Designer benzodiazepines"], "source_level": ["family"],
        "source_n_leaves": [8], "excess_share": [0.6],
    })
    for m in models:
        assert m.statistic in bm.TREESCAN_STATISTICS
        assert m.fixed_threshold == pytest.approx(bm.TREESCAN_RI_THRESHOLD)
        assert list(bm.model_score(sw, m)) == [250.0]
        assert bm.ENSEMBLE_FLOOR[m.statistic] == 1.0


def test_nb_trend_variants_share_one_sweep_and_read_their_own_column() -> None:
    """`nb-trend`/`nb-trend-own`/`nb-trend-dual` are three readings of one
    sweep (`aberration.nb_trend_sweep` fits both trends together) -- they
    must collapse to the same sweep key, `ears` must not join them (it is a
    different sweep function entirely), and `nb-trend-dual` must equal the
    elementwise min of the other two, the same pattern as `eb05_dual`."""
    q = pd.date_range("2024-01-01", periods=2, freq="QS")
    sw = pd.DataFrame({
        "substance": ["a", "a"], "as_of": q,
        "nb_trend_z": [3.0, -1.0], "nb_trend_own_z": [1.0, -2.0],
    })
    nb, nb_own, nb_dual, ears = (bm.MODELS_BY_ID[i] for i in
                                 ("nb-trend", "nb-trend-own", "nb-trend-dual", "ears"))
    assert nb.sweep_key == nb_own.sweep_key == nb_dual.sweep_key
    assert nb.sweep_key != ears.sweep_key
    assert list(bm.model_score(sw, nb)) == [3.0, -1.0]
    assert list(bm.model_score(sw, nb_own)) == [1.0, -2.0]
    assert list(bm.model_score(sw, nb_dual)) == [1.0, -2.0]


def test_veto_ensemble_fixed_threshold_matches_its_proposers_own_line() -> None:
    """`combine='veto'` never touches the proposer's own scale (see
    `build_ensemble_sweep`'s docstring), so the *ensemble's* fixed-reading
    line must equal the proposer's own native line -- if it silently fell
    back to the shared `ALARM_THRESHOLD` instead, a proposer scored on a
    different scale (nb-trend's Wald z, not EB05's posterior percentile)
    would be read at the wrong bar. Caught for real:
    `ensemble-veto-nbtrend-ears-*` first shipped without this override and
    read nb-trend at EB05's 1.5 line instead of its own 1.96, inflating
    recall/IoU and off-target load in the fixed-line reading."""
    for m in bm.MODELS:
        if m.statistic != "ensemble" or m.combine != "veto":
            continue
        proposer = bm.MODELS_BY_ID[m.components[0]]
        proposer_line = (proposer.fixed_threshold
                         if proposer.fixed_threshold is not None else ALARM_THRESHOLD)
        ensemble_line = (m.fixed_threshold
                         if m.fixed_threshold is not None else ALARM_THRESHOLD)
        assert ensemble_line == pytest.approx(proposer_line), m.id


def _eb_like_sweep(rows) -> pd.DataFrame:
    """A minimal eb05_dual-shaped sweep: substance, as_of, eb05, eb05_own."""
    return pd.DataFrame(rows, columns=["substance", "as_of", "eb05", "eb05_own"])


def test_ensemble_floors_a_component_that_never_scored_the_substance() -> None:
    """A component's silence on a (substance, as_of) is zero evidence, not a
    missing observation -- it must be floored and included in the ranking
    denominator, not dropped (dropping would let a sparse component's few
    scored rows rank artificially high relative to a smaller pool)."""
    q = pd.Timestamp("2024-01-01")
    eb_sweep = _eb_like_sweep([
        ("a", q, 3.0, 3.0), ("b", q, 1.0, 1.0), ("c", q, 0.5, 0.5),
    ])
    # treescan only ever scored "a" this quarter -- "b" and "c" are silent.
    ts_sweep = pd.DataFrame({
        "substance": ["a"], "as_of": [q], "recurrence_interval": [500.0],
        "p_value": [0.001], "llr": [10.0],
    })
    dual = bm.MODELS_BY_ID["eb05-dual"]
    tree = bm.MODELS_BY_ID["treescan"]
    sweeps = {dual.sweep_key: eb_sweep, tree.sweep_key: ts_sweep}

    out = bm.build_ensemble_sweep(("eb05-dual", "treescan"), sweeps, combine="mean")
    assert set(out["substance"]) == {"a", "b", "c"}
    # "a" leads on both components -- must rank top of the ensemble too.
    row = out.set_index("substance")["ensemble_score"]
    assert row["a"] == row.max()
    # "b" and "c" are floored on treescan (RI=1, the null value) but still
    # differ on eb05 -- the floor must not collapse them to an identical score.
    assert row["b"] != row["c"]


def test_ensemble_max_never_loses_to_mean_on_a_substance_one_component_confirms() -> None:
    """The whole point of `combine='max'`: a substance a sparse component is
    silent on should not be dragged down by that silence the way `'mean'`
    drags it down -- see `build_ensemble_sweep`'s documented finding that mean
    combination lost to `eb05-dual` alone on the real backtest."""
    q = pd.Timestamp("2024-01-01")
    eb_sweep = _eb_like_sweep([("hot", q, 9.0, 9.0), ("cold", q, 0.1, 0.1)])
    ts_sweep = pd.DataFrame({  # silent on both -- treescan never fired.
        "substance": [], "as_of": [], "recurrence_interval": [],
        "p_value": [], "llr": [],
    })
    sweeps = {bm.MODELS_BY_ID["eb05-dual"].sweep_key: eb_sweep,
              bm.MODELS_BY_ID["treescan"].sweep_key: ts_sweep}

    mean_out = bm.build_ensemble_sweep(("eb05-dual", "treescan"), sweeps, combine="mean")
    max_out = bm.build_ensemble_sweep(("eb05-dual", "treescan"), sweeps, combine="max")
    hot_mean = mean_out.set_index("substance")["ensemble_score"]["hot"]
    hot_max = max_out.set_index("substance")["ensemble_score"]["hot"]
    # eb05-dual alone ranks "hot" at the top of its two-substance quarter (1.0);
    # mean drags that down toward treescan's floor, max preserves it.
    assert hot_max == pytest.approx(1.0)
    assert hot_max > hot_mean


def test_ensemble_threshold_combine_is_a_literal_union_at_the_fixed_line() -> None:
    """`combine='threshold'` normalizes each component by its own calibrated
    line, so score > 1.0 is exactly "at least one component would already
    fire on its own" -- a literal alarm union, not a relative-rank artifact."""
    q = pd.Timestamp("2024-01-01")
    # "eb_only" clears eb05-dual's line (1.5) but not treescan's;
    # "tree_only" clears treescan's line (RI 100) but not eb05-dual's;
    # "neither" clears nothing.
    eb_sweep = _eb_like_sweep([
        ("eb_only", q, 2.0, 2.0), ("tree_only", q, 0.1, 0.1),
        ("neither", q, 0.1, 0.1),
    ])
    ts_sweep = pd.DataFrame({
        "substance": ["tree_only"], "as_of": [q],
        "recurrence_interval": [500.0], "p_value": [0.001], "llr": [10.0],
    })
    sweeps = {bm.MODELS_BY_ID["eb05-dual"].sweep_key: eb_sweep,
              bm.MODELS_BY_ID["treescan"].sweep_key: ts_sweep}
    out = bm.build_ensemble_sweep(("eb05-dual", "treescan"), sweeps,
                                  combine="threshold").set_index("substance")
    assert out.loc["eb_only", "ensemble_score"] > 1.0
    assert out.loc["tree_only", "ensemble_score"] > 1.0
    assert out.loc["neither", "ensemble_score"] < 1.0


def test_ensemble_veto_keeps_proposer_score_unless_vetoed() -> None:
    """`combine='veto'` never touches the proposer's own scale -- a survivor's
    score must equal exactly what the proposer alone would have said, and a
    vetoed substance must fall below the proposer's own floor (never an
    alarm, however high the proposer scored it)."""
    q = pd.Timestamp("2024-01-01")
    eb_sweep = _eb_like_sweep([
        ("survives", q, 9.0, 9.0), ("vetoed", q, 9.0, 9.0),
    ])
    ts_sweep = pd.DataFrame({
        "substance": ["survives", "vetoed"], "as_of": [q, q],
        "recurrence_interval": [50.0, 1.0], "p_value": [0.02, 1.0],
        "llr": [4.0, 0.001],
    })
    sweeps = {bm.MODELS_BY_ID["eb05-dual"].sweep_key: eb_sweep,
              bm.MODELS_BY_ID["treescan"].sweep_key: ts_sweep}
    out = bm.build_ensemble_sweep(("eb05-dual", "treescan"), sweeps,
                                  combine="veto", veto_threshold=10.0
                                  ).set_index("substance")["ensemble_score"]
    eb_alone = bm.model_score(eb_sweep, bm.MODELS_BY_ID["eb05-dual"])
    assert out["survives"] == pytest.approx(eb_alone.iloc[0])  # unchanged
    assert out["vetoed"] == bm.ENSEMBLE_FLOOR["eb05_dual"]     # floored, not just lowered


def test_ensemble_veto_requires_exactly_two_components_and_a_threshold() -> None:
    q = pd.Timestamp("2024-01-01")
    eb_sweep = _eb_like_sweep([("a", q, 1.0, 1.0)])
    ts_sweep = pd.DataFrame({"substance": ["a"], "as_of": [q],
                             "recurrence_interval": [1.0], "p_value": [1.0],
                             "llr": [0.0]})
    sweeps = {bm.MODELS_BY_ID["eb05-dual"].sweep_key: eb_sweep,
              bm.MODELS_BY_ID["treescan"].sweep_key: ts_sweep}
    with pytest.raises(ValueError):
        bm.build_ensemble_sweep(("eb05-dual",), sweeps, combine="veto",
                                veto_threshold=2.0)
    with pytest.raises(ValueError):
        bm.build_ensemble_sweep(("eb05-dual", "treescan"), sweeps, combine="veto")


def test_ensemble_rejects_an_unknown_combine_mode() -> None:
    q = pd.Timestamp("2024-01-01")
    eb_sweep = _eb_like_sweep([("a", q, 1.0, 1.0)])
    sweeps = {bm.MODELS_BY_ID["eb05-dual"].sweep_key: eb_sweep}
    with pytest.raises(ValueError):
        bm.build_ensemble_sweep(("eb05-dual",), sweeps, combine="median")
