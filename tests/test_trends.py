"""
Tests for the GPS v2 prototypes (`docs/GPS_V2_DESIGN.md`): the position-
weighted case definition and its dispersion correction, and the dual gate.

The dangerous failure for the credit scheme was not a crash -- it ran fine
and produced numbers -- it was a systematically inflated false-alarm rate
from an overdispersed weighted sum being fed to a Poisson-based fit as if it
were a plain count. `test__credit_dispersion_is_one_for_binary_credit` pins
down the property the fix (and the {0,1} fallback) actually rests on: for any
0/1-valued credit, E[c^2] == E[c] identically, so phi == 1 no matter how the
keep/discard split varies case to case -- not merely close to 1 empirically.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging.analysis import trends as tr


# --------------------------------------------------------------------------
# Position credit (component 1)
# --------------------------------------------------------------------------

def test_credit_from_order_solo_named_trailing() -> None:
    order = {
        "solo": ["Fentanyl"],
        "lead_of_3": ["Fentanyl", "Meth", "Lidocaine"],
        "trailing_of_3": ["Fentanyl", "Meth", "Lidocaine"],
        "last_of_2": ["Fentanyl", "Xylazine"],
    }
    credit = tr._credit_from_order(order).set_index(["CaseNumber", "substance"])["credit"]

    assert credit[("solo", "Fentanyl")] == tr.SOLO_CREDIT
    assert credit[("lead_of_3", "Fentanyl")] == tr.NAMED_CREDIT
    assert credit[("trailing_of_3", "Lidocaine")] == tr.TRAILING_CREDIT
    # A 2-substance line can't distinguish "trailing adulterant" from "just
    # not first" -- last-of-2 keeps full credit, not the trailing discount.
    assert credit[("last_of_2", "Xylazine")] == tr.NAMED_CREDIT


def test_credit_from_order_requires_three_to_discount() -> None:
    """Only a 3+-substance line's last slot is treated as the adulterant
    pattern -- a 2-line's last slot is `NAMED_CREDIT`, not `TRAILING_CREDIT`,
    to avoid the short-line bias `polysubstance.profile_table` documents."""
    order = {"c1": ["A", "B"], "c2": ["A", "B", "C"]}
    credit = tr._credit_from_order(order).set_index(["CaseNumber", "substance"])["credit"]
    assert credit[("c1", "B")] != tr.TRAILING_CREDIT
    assert credit[("c2", "C")] == tr.TRAILING_CREDIT


# --------------------------------------------------------------------------
# "Ever independent" -- the xylazine fix
# --------------------------------------------------------------------------

def test_ever_independent_is_solo_or_lead_at_least_once() -> None:
    order = {
        "a": ["Xylazine"],                              # solo
        "b": ["Fentanyl", "Cocaine", "Xylazine"],        # xylazine trailing
        "c": ["Fentanyl", "Meth", "Lidocaine"],          # lidocaine trailing
        "d": ["Fentanyl", "Meth", "Lidocaine"],          # lidocaine trailing
    }
    ind = tr._ever_independent(order)
    assert "Xylazine" in ind      # solo once is enough
    assert "Lidocaine" not in ind  # never solo or lead, anywhere
    assert "Fentanyl" in ind       # lead in b, c, d


def test_ever_independent_reproduces_the_real_lidocaine_xylazine_split() -> None:
    """Mirrors the actual case counts pulled from this repository's data
    (docs/GPS_V2_DESIGN.md): lidocaine never solo/lead across 42 mentions,
    xylazine solo exactly once across 9. Synthetic version of that shape."""
    order = {}
    for i in range(29):
        order[f"lido_trail_{i}"] = ["Fentanyl", "Meth", "Lidocaine"]
    for i in range(13):
        order[f"lido_mid_{i}"] = ["Fentanyl", "Lidocaine", "Meth"]
    for i in range(3):
        order[f"xyl_trail_{i}"] = ["Fentanyl", "Heroin", "Xylazine"]
    for i in range(5):
        order[f"xyl_mid_{i}"] = ["Fentanyl", "Xylazine"]
    order["xyl_solo"] = ["Xylazine"]

    ind = tr._ever_independent(order)
    assert "Lidocaine" not in ind
    assert "Xylazine" in ind


def test_credit_from_order_exempts_independent_substances_from_trailing() -> None:
    order = {"a": ["Fentanyl", "Heroin", "Xylazine"]}  # xylazine trailing
    plain = tr._credit_from_order(order)
    assert plain.set_index("substance").loc["Xylazine", "credit"] == tr.TRAILING_CREDIT

    exempted = tr._credit_from_order(order, independent=frozenset({"Xylazine"}))
    assert exempted.set_index("substance").loc["Xylazine", "credit"] == tr.NAMED_CREDIT


# --------------------------------------------------------------------------
# Dispersion (the bug the user caught)
# --------------------------------------------------------------------------

def test_dispersion_is_one_for_unweighted_credit() -> None:
    """Every mention worth 1 -- the unweighted case -- must be dispersion-
    neutral: phi == 1 exactly, matching `load_quarterly(weighted=False)`."""
    c = pd.Series([1] * 50)
    assert tr._credit_dispersion(c) == pytest.approx(1.0)


def test_dispersion_is_one_for_any_binary_credit_mix() -> None:
    """The property the {0,1} fallback in the design note rests on: for any
    0/1-valued credit, c^2 == c pointwise, so phi == E[c^2]/E[c] == 1
    identically -- not approximately -- however the keep/discard split
    varies. Checked at several mixes, including ones far from 50/50."""
    for p_keep in (0.05, 0.3, 0.5, 0.7, 0.95):
        rng = np.random.default_rng(0)
        c = pd.Series(rng.choice([0, 1], size=20000, p=[1 - p_keep, p_keep]))
        assert tr._credit_dispersion(c) == pytest.approx(1.0, abs=1e-9)


def test_dispersion_exceeds_one_for_the_three_tier_credit() -> None:
    """The {0, 1, 2} scheme actually shipped is *not* dispersion-neutral --
    this is the bug: SOLO_CREDIT == 2 makes c^2 > c whenever c == 2, so phi
    must exceed 1 for any mix that uses the solo tier at all."""
    c = pd.Series([0, 1, 1, 1, 2] * 1000)
    phi = tr._credit_dispersion(c)
    assert phi > 1.0
    # E[c^2]/E[c] by hand for this exact mix: values {0,1,1,1,2}, mean
    # (0+1+1+1+2)/5 = 1.0, mean-of-squares (0+1+1+1+4)/5 = 1.4 -> phi 1.4.
    assert phi == pytest.approx(1.4)


def test_dispersion_matches_hand_computation() -> None:
    c = pd.Series([0, 2, 2, 1])
    # mean = 1.25, mean(c^2) = (0+4+4+1)/4 = 2.25, phi = 2.25/1.25 = 1.8
    assert tr._credit_dispersion(c) == pytest.approx(1.8)


# --------------------------------------------------------------------------
# rank_substances: dual gate (component 2) and the phi correction
# --------------------------------------------------------------------------

def _synthetic_panel():
    """Four quarters (2 baseline, 2 recent), a shrinking denominator (400 ->
    200 OD deaths total), and three kinds of substance:

    - 12 background substances whose own count halves in step with the
      denominator, so their *share* stays flat -- enough non-zero-baseline
      cells for `_fit_gps_prior` to actually fit rather than fall back.
    - "Cocaine": own count exactly flat (40 -> 40) while the total shrinks,
      so its *share* rises credibly (Gate A) with no rise in its own count
      (Gate B) -- the pattern `regime()` currently requires a reader to spot
      by cross-referencing `count_ratio` by hand.
    - "TrueRiser": both count and share rise sharply -- should clear both
      gates.
    """
    quarters = list(pd.date_range("2024-01-01", periods=4, freq="QS"))
    base_q, recent_q = quarters[:2], quarters[2:]
    denom = pd.Series({q: 200 for q in base_q} | {q: 100 for q in recent_q})

    rows = []
    for i in range(12):
        sub = f"Background{i}"
        for q in base_q:
            rows.append((sub, q, 20))
        for q in recent_q:
            rows.append((sub, q, 10))
    for q in base_q:
        rows.append(("Cocaine", q, 20))
    for q in recent_q:
        rows.append(("Cocaine", q, 20))
    for q in base_q:
        rows.append(("TrueRiser", q, 5))
    for q in recent_q:
        rows.append(("TrueRiser", q, 30))

    counts = pd.DataFrame(rows, columns=["substance", "quarter", "n"])
    return counts, denom


def test_gate_a_alone_flags_a_growing_share_of_a_shrinking_total() -> None:
    """Sanity check on the fixture: Cocaine's share-based EB05 alone reads as
    a credible rise, reproducing the real cocaine finding in RESEARCH_LOG.md
    (EB05 1.06 while count_ratio was 0.90) at a scale where it's unambiguous."""
    counts, denom = _synthetic_panel()
    res, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                recent_quarters=2, baseline_quarters=2)
    assert res.loc["Cocaine", "eb05"] > tr.ALARM_THRESHOLD


def test_dual_gate_rejects_share_only_rise() -> None:
    """Cocaine's own count never rose -- Gate B must not clear, and
    `credible_rise` must be False even though Gate A alone fires."""
    counts, denom = _synthetic_panel()
    res, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                recent_quarters=2, baseline_quarters=2)
    assert res.loc["Cocaine", "eb05_own"] < tr.ALARM_THRESHOLD
    assert not bool(res.loc["Cocaine", "credible_rise"])


def test_dual_gate_accepts_a_genuine_rise_in_both() -> None:
    counts, denom = _synthetic_panel()
    res, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                recent_quarters=2, baseline_quarters=2)
    assert res.loc["TrueRiser", "eb05"] > tr.ALARM_THRESHOLD
    assert res.loc["TrueRiser", "eb05_own"] > tr.ALARM_THRESHOLD
    assert bool(res.loc["TrueRiser", "credible_rise"])


def test_phi_defaults_to_one_when_unset() -> None:
    """A `counts` frame built by hand (no `.attrs`), as every caller before
    this change would pass, must score identically to phi == 1 -- the
    unweighted path must be untouched by this change."""
    counts, denom = _synthetic_panel()
    assert "phi" not in counts.attrs
    res, meta = tr.rank_substances(counts, denom, lag_quarters=0,
                                   recent_quarters=2, baseline_quarters=2)
    assert meta["phi"] == pytest.approx(1.0)


def test_phi_correction_is_more_conservative() -> None:
    """A larger phi (more overdispersed credit) must pull EB05 down for a
    substance with a strong raw signal -- the effective sample size shrinks,
    which is the entire point of the correction."""
    counts, denom = _synthetic_panel()
    res_1, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                  recent_quarters=2, baseline_quarters=2)

    counts_phi = counts.copy()
    counts_phi.attrs["phi"] = 2.0
    res_2, meta_2 = tr.rank_substances(counts_phi, denom, lag_quarters=0,
                                       recent_quarters=2, baseline_quarters=2)

    assert meta_2["phi"] == pytest.approx(2.0)
    assert res_2.loc["TrueRiser", "eb05"] < res_1.loc["TrueRiser", "eb05"]


# --------------------------------------------------------------------------
# Conditional-Poisson diagnostic
# --------------------------------------------------------------------------
#
# The trap this diagnostic exists to avoid: a substance whose rate genuinely
# moved fails a *flat* Poisson test for the best possible reason. If `d_trend`
# did not separate that case from real clustering, the check would condemn
# exactly the substances the ranking is built to find --
# `test_a_pure_trend_is_not_clustered` is the one that pins that down.

def _panel_from_series(series: dict[str, list[int]],
                       denom: list[int]) -> tuple[pd.DataFrame, pd.Series]:
    """Build (counts, denom) in `load_quarterly` shape from literal series."""
    qs = pd.period_range("2020Q1", periods=len(denom), freq="Q").to_timestamp()
    rows = [{"substance": s, "quarter": q, "n": n}
            for s, vals in series.items() for q, n in zip(qs, vals)]
    return pd.DataFrame(rows), pd.Series(denom, index=qs, dtype=float)


def test_poisson_data_has_dispersion_near_one() -> None:
    """Counts generated as an honest Poisson around a constant share must land
    near D/df == 1 on both statistics, and be labelled `ok`."""
    rng = np.random.default_rng(0)
    denom = [1000] * 24
    obs = rng.poisson(60, size=24).tolist()
    counts, den = _panel_from_series({"Clean": obs}, denom)
    d = tr.poisson_dispersion(counts, den)
    assert d.loc["Clean", "d_flat"] == pytest.approx(1.0, abs=0.45)
    assert d.loc["Clean", "d_trend"] == pytest.approx(1.0, abs=0.45)
    assert d.loc["Clean", "verdict"] == "ok"


def test_a_pure_trend_is_not_clustered() -> None:
    """The property the whole diagnostic rests on. A substance with a real
    log-linear trend and otherwise perfect Poisson noise must show a large
    `d_flat` (its rate moved) but a `d_trend` near 1 (given the trend, it is
    well behaved) -- and must NOT be labelled `clustered`. Condemning this
    case would condemn every substance the ranking is looking for."""
    rng = np.random.default_rng(1)
    denom = [1000] * 24
    mu = 20 * np.exp(np.linspace(0, 1.6, 24))     # ~5x rise across the window
    obs = rng.poisson(mu).tolist()
    counts, den = _panel_from_series({"Riser": obs}, denom)
    d = tr.poisson_dispersion(counts, den)

    assert d.loc["Riser", "d_flat"] > 5.0          # flat test is badly rejected
    assert d.loc["Riser", "d_trend"] == pytest.approx(1.0, abs=0.5)
    assert d.loc["Riser", "absorbed"] > 4.0        # nearly all of it was trend
    assert d.loc["Riser", "verdict"] != "clustered"


def test_clustered_arrivals_are_caught() -> None:
    """Deaths arriving in batches -- the contaminated-supply mechanism -- must
    be flagged, and must stay flagged after detrending, because no trend
    explains them."""
    rng = np.random.default_rng(2)
    denom = [1000] * 24
    obs = (3 * rng.poisson(8, size=24)).tolist()   # same mean, 3x the variance
    counts, den = _panel_from_series({"Batchy": obs}, denom)
    d = tr.poisson_dispersion(counts, den)

    assert d.loc["Batchy", "d_trend"] > 2.0
    assert d.loc["Batchy", "p_trend"] < 0.05
    assert d.loc["Batchy", "verdict"] == "clustered"
    # The reported rescaling is sqrt of the dispersion, per quasi-likelihood.
    assert d.loc["Batchy", "interval_scale"] == pytest.approx(
        np.sqrt(d.loc["Batchy", "d_trend"]))


def test_underdispersion_reads_as_smooth_not_clustered() -> None:
    """Smoother-than-Poisson is a violation too, but the benign one: it makes
    EB05 conservative rather than overconfident, so it must get its own label
    and never be reported as `clustered`."""
    denom = [1000] * 24
    counts, den = _panel_from_series({"Flat": [50] * 24}, denom)   # zero variance
    d = tr.poisson_dispersion(counts, den)
    assert d.loc["Flat", "d_trend"] < tr.DISPERSION_SMOOTH
    assert d.loc["Flat", "verdict"] == "smooth"


def test_constant_share_against_a_falling_denominator_is_clean() -> None:
    """The offset has to be the denominator, not the quarter index. A
    substance holding a perfectly constant *share* while total OD deaths fall
    22% (LA's actual pattern) has not moved at all, and must not register as
    trending or overdispersed."""
    denom = [int(1000 * 0.98 ** i) for i in range(24)]
    obs = [int(round(0.05 * n)) for n in denom]     # exactly 5% every quarter
    counts, den = _panel_from_series({"Steady": obs}, denom)
    d = tr.poisson_dispersion(counts, den)
    assert d.loc["Steady", "d_flat"] < 0.2
    assert d.loc["Steady", "d_trend"] < 0.2
    assert d.loc["Steady", "verdict"] == "smooth"


def test_sparse_substances_are_excluded_rather_than_guessed() -> None:
    """Below `min_total` the chi-square approximation does not hold, so the
    substance must be dropped from the table entirely rather than given a
    number the sample cannot support."""
    denom = [1000] * 24
    counts, den = _panel_from_series(
        {"Rare": [1] * 24, "Common": [40] * 24}, denom)
    d = tr.poisson_dispersion(counts, den, min_total=tr.MIN_TOTAL_FOR_DISPERSION)
    assert "Rare" not in d.index
    assert "Common" in d.index


def test_trend_fit_recovers_a_known_slope() -> None:
    """`_poisson_trend_mu` is a hand-rolled Nelder-Mead stand-in for a Poisson
    GLM (statsmodels is not a declared dependency). Check it against data with
    a known noiseless log-linear mean."""
    denom = np.full(20, 1000.0)
    mu_true = 30 * np.exp(np.linspace(0, 1.0, 20))
    fitted = tr._poisson_trend_mu(mu_true, denom)
    assert np.allclose(fitted, mu_true, rtol=1e-3)


# --------------------------------------------------------------------------
# Spatial concentration (component 3)
# --------------------------------------------------------------------------
#
# The regression that matters most here is not a wrong number, it is a moved
# one: `spatial_mode="none"` is the default every existing caller and every
# committed result was produced under, so it has to be *bit-identical* to the
# code before this parameter existed.
# `test_spatial_none_is_identical_to_no_spatial_argument` is the guard.

def test_spatial_none_is_identical_to_no_spatial_argument() -> None:
    counts, denom = _synthetic_panel()
    z = pd.Series({"TrueRiser": 4.0, "Cocaine": -2.0})
    plain, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                  recent_quarters=2, baseline_quarters=2)
    with_z, meta = tr.rank_substances(counts, denom, lag_quarters=0,
                                      recent_quarters=2, baseline_quarters=2,
                                      z_spatial=z, spatial_mode="none")
    assert meta["spatial_mode"] == "none"
    for col in ("eb05", "eb05_own", "ebgm", "expected"):
        assert np.allclose(plain[col], with_z.loc[plain.index, col])
    # An ignored covariate must not leak into the reported columns either.
    assert (with_z["z_spatial"] == 0).all()


def test_spatial_discount_raises_eb05_only_for_concentrated_substances() -> None:
    """The discount shrinks `expected`, so the same observed count becomes
    more surprising -- but only where there is spatial evidence. A dispersed
    or unscored substance must be left exactly where it was."""
    counts, denom = _synthetic_panel()
    z = pd.Series({"TrueRiser": 4.0})       # everything else -> 0
    plain, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                  recent_quarters=2, baseline_quarters=2)
    disc, meta = tr.rank_substances(counts, denom, lag_quarters=0,
                                    recent_quarters=2, baseline_quarters=2,
                                    z_spatial=z, spatial_mode="discount")

    assert meta["spatial_mode"] == "discount"
    assert disc.loc["TrueRiser", "expected"] < plain.loc["TrueRiser", "expected"]
    assert disc.loc["TrueRiser", "eb05"] > plain.loc["TrueRiser", "eb05"]
    assert disc.loc["Cocaine", "expected"] == pytest.approx(
        plain.loc["Cocaine", "expected"])


def test_spatial_discount_is_bounded_by_its_weight() -> None:
    """`conc` is in [0, 1], so E can never fall below (1 - w) * E however
    concentrated a substance is -- geography alone cannot manufacture an
    alarm out of a flat count."""
    counts, denom = _synthetic_panel()
    z = pd.Series({s: 50.0 for s in counts["substance"].unique()})
    disc, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                 recent_quarters=2, baseline_quarters=2,
                                 z_spatial=z, spatial_mode="discount",
                                 spatial_weight=0.25)
    plain, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                  recent_quarters=2, baseline_quarters=2)
    ratio = disc["expected"] / plain["expected"]
    assert np.allclose(ratio, 0.75)


def test_covariate_prior_recovers_the_plain_fit_when_z_is_constant() -> None:
    """A covariate with no variation carries no information: the fit must
    drive gamma to the point where it changes nothing, and the posteriors must
    match the unadjusted ones. This is the property that makes the estimated
    covariate safer than the chosen discount -- it can decline to act."""
    counts, denom = _synthetic_panel()
    z = pd.Series({s: 2.0 for s in counts["substance"].unique()})
    plain, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                  recent_quarters=2, baseline_quarters=2)
    prior, meta = tr.rank_substances(counts, denom, lag_quarters=0,
                                     recent_quarters=2, baseline_quarters=2,
                                     z_spatial=z, spatial_mode="prior")
    assert meta["prior"]["z_sd"] == pytest.approx(1.0)   # zero-variance guard
    assert np.allclose(plain["eb05"], prior.loc[plain.index, "eb05"], atol=1e-6)


def test_fit_gps_prior_reports_zero_gamma_without_a_covariate() -> None:
    """`_gps_scores` keys off a truthy gamma to decide whether to build a
    per-substance beta at all, so the no-covariate fit must report 0.0 rather
    than omitting the key."""
    rng = np.random.default_rng(0)
    e = rng.uniform(1, 40, size=120)
    n = rng.poisson(e)
    pr = tr._fit_gps_prior(n.astype(float), e)
    assert pr["gamma"] == 0.0
    assert np.allclose(tr._gps_scores(n.astype(float), e, pr)["eb05"],
                       tr._gps_scores(n.astype(float), e, pr, z=None)["eb05"])


def test_unknown_spatial_mode_is_rejected() -> None:
    counts, denom = _synthetic_panel()
    with pytest.raises(ValueError, match="spatial_mode"):
        tr.rank_substances(counts, denom, lag_quarters=0, recent_quarters=2,
                           baseline_quarters=2, spatial_mode="multiply")


# --------------------------------------------------------------------------
# _apply_treescan_veto -- the default detector's TreeScan cross-check
# --------------------------------------------------------------------------

def _credible_res(rows: dict[str, bool]) -> pd.DataFrame:
    return pd.DataFrame({"credible_rise": rows}).astype({"credible_rise": bool})


def test_treescan_veto_suppresses_only_below_the_floor() -> None:
    res = _credible_res({"low_ri": True, "high_ri": True, "not_rising": True})
    ri = pd.Series({"low_ri": 2.0, "high_ri": 50.0, "not_rising": 2.0})
    out = tr._apply_treescan_veto(res, ri, threshold=10.0)
    assert not bool(out.loc["low_ri", "credible_rise"])
    assert bool(out.loc["high_ri", "credible_rise"])
    assert bool(out.loc["low_ri", "treescan_veto"])
    assert not bool(out.loc["high_ri", "treescan_veto"])
    # not_rising never had credible_rise in the first place -- the veto column
    # must not claim it vetoed something that was never an alarm.
    assert not bool(out.loc["not_rising", "credible_rise"])


def test_treescan_veto_never_promotes_a_non_alarm() -> None:
    """A substance TreeScan loves but EB05 doesn't (`credible_rise` already
    False) must stay False -- this is a veto, not a second way to alarm."""
    res = _credible_res({"eb05_quiet": False})
    ri = pd.Series({"eb05_quiet": 1000.0})
    out = tr._apply_treescan_veto(res, ri, threshold=10.0)
    assert not bool(out.loc["eb05_quiet", "credible_rise"])
    assert not bool(out.loc["eb05_quiet", "treescan_veto"])


def test_treescan_veto_floors_a_substance_absent_from_the_ri_series() -> None:
    """Absent from TreeScan's leaf sweep means no excess at all (RI's own
    floor is 1) -- not a missing observation to be left un-vetoed."""
    res = _credible_res({"never_scanned": True})
    ri = pd.Series({}, dtype=float)
    out = tr._apply_treescan_veto(res, ri, threshold=10.0)
    assert out.loc["never_scanned", "treescan_ri"] == 1.0
    assert not bool(out.loc["never_scanned", "credible_rise"])


# --------------------------------------------------------------------------
# _role_dispersion -- the substance-level role discount
# --------------------------------------------------------------------------

def _lidocaine_shaped_order(n_trailing: int, n_pair: int) -> dict[str, list[str]]:
    """`n_trailing` cases where the substance trails a 3+-line (credit 0),
    `n_pair` cases where it's the second of a plain 2-item line with the same
    partner (credit 1, same shape as Lidocaine's real `[Fentanyl, Lidocaine]`
    pattern) -- enough of each to clear `MIN_CASES_FOR_ROLE_DISCOUNT`."""
    order = {}
    for i in range(n_trailing):
        order[f"trail_{i}"] = ["Fentanyl", "Methamphetamine", "X"]
    for i in range(n_pair):
        order[f"pair_{i}"] = ["Fentanyl", "X"]
    return order


def test_role_dispersion_is_one_for_a_substance_that_is_often_solo() -> None:
    order = {f"solo_{i}": ["X"] for i in range(tr.MIN_CASES_FOR_ROLE_DISCOUNT)}
    phi_s = tr._role_dispersion(order, independent=frozenset())
    assert phi_s["X"] == pytest.approx(1.0)


def test_role_dispersion_is_large_for_an_always_trailing_substance() -> None:
    order = _lidocaine_shaped_order(n_trailing=30, n_pair=0)
    phi_s = tr._role_dispersion(order, independent=frozenset())
    # mean credit -> 0, floored at _ROLE_MEAN_CREDIT_FLOOR -> phi_s at its cap
    assert phi_s["X"] == pytest.approx(tr.NAMED_CREDIT / tr._ROLE_MEAN_CREDIT_FLOOR)


def test_role_dispersion_below_min_cases_is_omitted() -> None:
    order = _lidocaine_shaped_order(n_trailing=tr.MIN_CASES_FOR_ROLE_DISCOUNT - 1,
                                    n_pair=0)
    phi_s = tr._role_dispersion(order, independent=frozenset())
    assert "X" not in phi_s


def test_role_dispersion_matches_the_lidocaine_shape() -> None:
    """22 trailing, 12 short-pair -- this repository's own real Lidocaine mix
    (trends.py's `_role_dispersion` docstring). mean credit = 12/34; phi_s =
    NAMED_CREDIT / mean, hand-computed (the reference point is 1, not
    SOLO_CREDIT's 2 -- see the docstring for why the first cut had this
    wrong and what it broke)."""
    order = _lidocaine_shaped_order(n_trailing=22, n_pair=12)
    phi_s = tr._role_dispersion(order, independent=frozenset())
    mean_credit = 12 / 34
    assert phi_s["X"] == pytest.approx(tr.NAMED_CREDIT / mean_credit)


def test_role_dispersion_never_amplifies_a_well_behaved_substance() -> None:
    """The floor that matters most in practice: a substance exempted by
    `_ever_independent` (so no TRAILING_CREDIT at all) can still have mean
    credit above NAMED_CREDIT purely from occasional solo mentions -- e.g.
    Fentanyl's real mean credit is 1.223. Without the max(1, ...) floor this
    would give phi_s < 1, which *amplifies* evidence instead of discounting
    it, backwards for what this function is for."""
    order = {}
    for i in range(9):
        order[f"named_{i}"] = ["A", "X"]        # credit 1
    for i in range(3):
        order[f"solo_{i}"] = ["X"]               # credit 2
    phi_s = tr._role_dispersion(order, independent=frozenset())
    mean_credit = (9 * 1 + 3 * 2) / 12
    assert mean_credit > tr.NAMED_CREDIT        # confirms the trap is real
    assert phi_s["X"] == pytest.approx(1.0)      # floored, not amplified


def test_role_dispersion_matches_this_repositorys_real_bromazolam_regression() -> None:
    """Bromazolam's real mean credit is 1.000 (33 mentions, all NAMED_CREDIT,
    `_ever_independent`-exempted, never solo). The SOLO_CREDIT-referenced
    first cut gave this phi_s = 2.0 and broke `backtest` across every known
    emergence; pinned down here so it can't regress back to that formula."""
    order = {f"m_{i}": ["Fentanyl", "Bromazolam"] for i in range(33)}
    independent = frozenset({"Bromazolam"})
    phi_s = tr._role_dispersion(order, independent)
    assert phi_s["Bromazolam"] == pytest.approx(1.0)


def test_role_dispersion_is_as_of_safe_like_ever_independent() -> None:
    """Restricting `order` to fewer cases (the as-of discipline every caller
    must apply -- see the docstring) changes phi_s, the same way it changes
    `_ever_independent`; this just pins down that the function actually reads
    from its `order` argument and doesn't cache anything across calls."""
    full = _lidocaine_shaped_order(n_trailing=22, n_pair=12)
    early = dict(list(full.items())[:15])  # an earlier as-of slice
    phi_full = tr._role_dispersion(full, independent=frozenset())
    phi_early = tr._role_dispersion(early, independent=frozenset())
    assert phi_full["X"] != phi_early["X"]


def test_role_discount_leaves_the_raw_ratio_unchanged_but_shrinks_eb05() -> None:
    """The property the whole design rests on: phi_s divides n and E by the
    same amount, so lambda = n/E is exactly invariant, but EB05 (the shrunk
    posterior) is not, because the population prior is fit once and does not
    rescale with one substance's own phi_s."""
    counts, denom = _synthetic_panel()
    counts_disc = counts.copy()
    counts_disc.attrs["phi"] = 1.0
    counts_disc.attrs["phi_s"] = {"TrueRiser": 5.0}

    plain, _ = tr.rank_substances(counts, denom, lag_quarters=0,
                                  recent_quarters=2, baseline_quarters=2)
    disc, _ = tr.rank_substances(counts_disc, denom, lag_quarters=0,
                                 recent_quarters=2, baseline_quarters=2)

    # rank_substances doesn't expose n/expected post-hoc-divided by phi_s
    # directly, but n_recent/expected in the output are the *raw* (pre-phi_s)
    # values either way -- phi_s only touches what goes into the fit, not
    # what's reported -- so the ratio computed from the table is identical
    # by construction. The real claim is on eb05 itself.
    ratio_plain = plain.loc["TrueRiser", "n_recent"] / plain.loc["TrueRiser", "expected"]
    ratio_disc = disc.loc["TrueRiser", "n_recent"] / disc.loc["TrueRiser", "expected"]
    assert ratio_plain == pytest.approx(ratio_disc)
    assert disc.loc["TrueRiser", "eb05"] < plain.loc["TrueRiser", "eb05"]
    # Every other substance's phi_s is still 1 (unnamed in the dict), but its
    # *posterior* still moves a little: alpha/beta are fit jointly across all
    # substances, so shrinking TrueRiser's own (n, E) shifts the population
    # fit slightly, which ripples to everyone's eb05 -- expected, not a bug,
    # and worth pinning down so it isn't mistaken for one later. The shift is
    # small relative to the substance actually being discounted.
    others = plain.index.difference(["TrueRiser"])
    other_shift = (plain.loc[others, "eb05"] - disc.loc[others, "eb05"]).abs().max()
    own_shift = plain.loc["TrueRiser", "eb05"] - disc.loc["TrueRiser", "eb05"]
    assert 0 < other_shift < own_shift
