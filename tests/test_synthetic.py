"""
Tests for the synthetic benchmark generator (`emerging/validation/synthetic.py`,
`docs/SYNTHETIC_BENCHMARK_DESIGN.md`).

Two kinds of failure matter here, and this file is organized around them:

- **A miscalibrated sampler silently changes what a scenario means.** The
  clustered-noise sampler's first cut sampled at 1.42x its intended rate
  because a spawned parent's own children were extra expectation the naive
  split never accounted for -- caught by exactly the moment check below, not
  by inspection. Any noise/shape function that claims a target rate needs a
  test asserting it actually hits that rate, not just that it runs.
- **A generator that shares a detector's own assumptions doesn't test
  anything** (the module docstring's circularity risk) -- so the shape
  functions are tested for behaving like something *other* than a smooth
  log-linear trend where that's the point (step jumps discretely, batch
  concentrates excess in one or two quarters, not spread as a rate).
"""

from __future__ import annotations

import numpy as np
import pytest

from emerging.ingest.extract import MENTIONS_PATH
from emerging.validation import benchmark as bm
from emerging.validation.ground_truth import load_ground_truth
from emerging.validation.synthetic import (GALLERY_PROFILES, Profile,
                                           _best_lag_correlation, _rate_curve,
                                           _sample_counts, generate_panel,
                                           nearest_real_matches)

T = 20


def test_step_shape_jumps_discretely_at_onset() -> None:
    p = Profile("x", baseline_rate=1.0, rises=True, sustained=True,
               shape="step", magnitude=5.0, onset=10, duration=1)
    r = _rate_curve(p, T)
    assert r[9] == pytest.approx(1.0)
    assert r[10] == pytest.approx(5.0)
    assert r[-1] == pytest.approx(5.0)


def test_sigmoid_shape_crosses_the_midpoint_at_its_own_center() -> None:
    p = Profile("x", baseline_rate=1.0, rises=True, sustained=True,
               shape="sigmoid", magnitude=5.0, onset=10, duration=4)
    r = _rate_curve(p, T)
    assert r[9] == pytest.approx(1.0, abs=1e-6)   # flat before onset
    assert r[12] == pytest.approx(3.0, abs=0.01)  # onset + duration/2
    assert r[-1] == pytest.approx(5.0, abs=0.01)  # saturates near the peak


def test_batch_shape_concentrates_excess_rather_than_spreading_a_rate() -> None:
    """The whole point of this shape (module docstring, design note sec 2):
    unlike every other shape, the excess is a fixed lump in 1-2 quarters, not
    a smoothly elevated rate for the rest of the window."""
    p = Profile("x", baseline_rate=1.0, rises=True, sustained=True,
               shape="batch", magnitude=5.0, onset=10, duration=4)
    r = _rate_curve(p, T)
    elevated = np.flatnonzero(r > 1.0 + 1e-9)
    assert set(elevated) == {10, 11}
    assert r[-1] == pytest.approx(1.0)  # back to baseline well after the burst
    expected_excess = 1.0 * (5.0 - 1.0) * 4  # baseline * (magnitude-1) * duration
    assert r.sum() - T * 1.0 == pytest.approx(expected_excess)


def test_batch_shape_n_bursts_repeats_independent_disconnected_spikes() -> None:
    """`n_bursts > 1` is what actually distinguishes "batch" from a
    reverted blip (`Profile.real_positive`'s docstring) -- several separate
    elevated windows, none adjacent, each carrying the same excess a single
    burst would."""
    p = Profile("x", baseline_rate=1.0, rises=True, sustained=True,
               shape="batch", magnitude=5.0, onset=4, duration=4, n_bursts=3)
    r = _rate_curve(p, 40)
    elevated = np.flatnonzero(r > 1.0 + 1e-9)
    # Three separate 2-quarter clusters, not one contiguous elevated block.
    gaps_between_clusters = np.diff(elevated)
    assert (gaps_between_clusters > 1).sum() == 2  # two breaks -> three clusters
    per_burst_excess = 1.0 * (5.0 - 1.0) * 4
    assert r.sum() - 40 * 1.0 == pytest.approx(3 * per_burst_excess)


def test_perturb_at_shocks_only_its_own_quarter_in_the_opposite_direction() -> None:
    """A brief, opposite-direction shock riding on top of an otherwise
    smooth curve -- must not change the curve anywhere else, and must not
    change whether the profile reads as ground-truth positive."""
    base = Profile("base", baseline_rate=1.0, rises=True, sustained=True,
                   shape="log_linear", magnitude=9.0, onset=10, duration=6)
    shocked = Profile("shocked", baseline_rate=1.0, rises=True, sustained=True,
                      shape="log_linear", magnitude=9.0, onset=10, duration=6,
                      perturb_at=25, perturb_factor=0.1)
    r_base = _rate_curve(base, 40)
    r_shocked = _rate_curve(shocked, 40)
    np.testing.assert_allclose(np.delete(r_base, [25, 26]),
                               np.delete(r_shocked, [25, 26]))
    assert r_shocked[25] == pytest.approx(r_base[25] * 0.1)
    assert base.real_positive and shocked.real_positive  # the dip doesn't flip the label


def test_perturb_every_repeats_the_shock_periodically() -> None:
    """`perturb_every` riding on top of a real, sustained rise -- the
    established-substance-with-recurring-batches case -- must still read as
    ground-truth positive (the shocks are noise, not the label) and must
    hit multiple, evenly-spaced quarters, not just the first one."""
    p = Profile("x", baseline_rate=15.0, rises=True, sustained=True,
               shape="log_linear", magnitude=3.0, onset=10, duration=10,
               perturb_at=6, perturb_every=7, perturb_factor=2.2)
    r = _rate_curve(p, 40)
    unperturbed = _rate_curve(Profile("y", baseline_rate=15.0, rises=True,
                                      sustained=True, shape="log_linear",
                                      magnitude=3.0, onset=10, duration=10), 40)
    shocked_quarters = np.flatnonzero(r > unperturbed + 1e-9)
    assert len(shocked_quarters) >= 4  # several shocks, not just one
    gaps = np.diff(shocked_quarters)
    assert (gaps > 1).sum() >= 2  # separated clusters, not one contiguous block
    assert p.real_positive is True


def test_declining_magnitude_reads_as_a_negative_not_an_emergence() -> None:
    """`magnitude <= 1.0` is a decline, not a rise -- `Profile.real_positive`
    must exclude it even with `rises=True, sustained=True`, the same carve-out
    shape="batch" already gets."""
    p = Profile("fading", baseline_rate=15.0, rises=True, sustained=True,
               shape="log_linear", magnitude=0.2, onset=6, duration=8)
    assert p.real_positive is False
    assert p.emergent is None
    r = _rate_curve(p, 40)
    assert r[-1] < r[0]  # genuinely declined, not flat


def test_blip_reverts_to_baseline_after_duration() -> None:
    """The Morphine/Mirtazapine negative shape: a real, sustained-looking
    rise that reverts -- must not be confused with a genuine sustained one."""
    p = Profile("x", baseline_rate=1.0, rises=True, sustained=False,
               shape="step", magnitude=5.0, onset=10, duration=3)
    r = _rate_curve(p, T)
    assert r[10] == pytest.approx(5.0)
    assert r[12] == pytest.approx(5.0)   # still within [onset, onset+duration)
    assert r[13] == pytest.approx(1.0)   # reverted
    assert r[-1] == pytest.approx(1.0)


def test_flat_profile_never_rises_regardless_of_shape_params() -> None:
    p = Profile("x", baseline_rate=4.0, rises=False, shape="sigmoid",
               magnitude=10.0, onset=5)
    r = _rate_curve(p, T)
    assert np.allclose(r, 4.0)


@pytest.mark.parametrize("noise,dispersion,min_var_ratio", [
    ("poisson", 1.0, None),
    ("negbin", 3.0, 1.5),
    ("clustered", 1.0, 1.05),
])
def test_noise_sampler_hits_its_target_mean(noise, dispersion, min_var_ratio) -> None:
    """Every sampler here claims `E[count_t] == rate_t` regardless of noise
    model -- the property the clustered sampler's first cut violated by 42%.
    Checked at a large N so the assertion is about calibration, not luck."""
    rng = np.random.default_rng(0)
    rate = np.full(5000, 10.0)
    draws = _sample_counts(rate, noise, dispersion, rng)
    assert draws.mean() == pytest.approx(10.0, rel=0.05)
    if min_var_ratio is not None:
        assert draws.var() / draws.mean() > min_var_ratio


def test_clustered_noise_is_actually_autocorrelated() -> None:
    """Distinguishes `"clustered"` from `"negbin"`: both overdisperse, but
    only clustering should show positive lag-1 autocorrelation -- otherwise
    it is just negative-binomial noise with extra steps."""
    rng = np.random.default_rng(1)
    rate = np.full(3000, 10.0)
    draws = _sample_counts(rate, "clustered", 1.0, rng)
    ac = np.corrcoef(draws[:-1], draws[1:])[0, 1]
    assert ac > 0.1


def test_generate_panel_ground_truth_matches_profile_labels() -> None:
    """The three-way taxonomy from a single generation call, and the
    `emergent` column real ground truth cannot supply at all."""
    profiles = [
        Profile("flat", baseline_rate=5.0, rises=False),
        Profile("blip", baseline_rate=3.0, rises=True, sustained=False,
                shape="step", magnitude=5, onset=10, duration=2),
        Profile("rare_riser", baseline_rate=1.0, rises=True, sustained=True,
                shape="step", magnitude=8, onset=10, duration=1),
        Profile("established_riser", baseline_rate=20.0, rises=True,
                sustained=True, shape="log_linear", magnitude=2, onset=10,
                duration=4),
    ]
    counts, denom, gt = generate_panel(profiles, T=30, seed=0)
    gt = gt.set_index("substance")

    assert bool(gt.loc["flat", "no_emergence"])
    assert bool(gt.loc["blip", "no_emergence"])       # reverted -- a negative
    assert not bool(gt.loc["rare_riser", "no_emergence"])
    assert bool(gt.loc["rare_riser", "emergent"])
    assert not bool(gt.loc["established_riser", "no_emergence"])
    assert not bool(gt.loc["established_riser", "emergent"])


def test_batch_shape_reads_as_a_negative_even_with_sustained_true() -> None:
    """`shape="batch"` can never actually stay elevated -- `_rate_curve`
    concentrates its excess into 1-2 quarters regardless of `sustained`
    (`test_batch_shape_concentrates_excess_rather_than_spreading_a_rate`).
    Before `Profile.real_positive`'s carve-out, a batch profile with
    `sustained=True` and one with `sustained=False` produced bit-identical
    rate curves and sampled data but opposite ground truth -- caught by the
    synthetic-data profile gallery, where the two were visually
    indistinguishable despite the label disagreeing. `sustained` must not
    be able to flip a batch profile's ground truth at all."""
    always_batch = dict(baseline_rate=3.0, rises=True, shape="batch",
                        magnitude=5, onset=10, duration=2)
    sustained = Profile("batch_sustained", sustained=True, **always_batch)
    reverts = Profile("batch_reverts", sustained=False, **always_batch)
    assert sustained.real_positive is False
    assert sustained.emergent is None

    counts, denom, gt = generate_panel([sustained, reverts], T=30, seed=0)
    gt = gt.set_index("substance")
    assert bool(gt.loc["batch_sustained", "no_emergence"])
    assert bool(gt.loc["batch_reverts", "no_emergence"])
    assert (gt.loc["batch_sustained", ["interval_start", "interval_end"]]
            == gt.loc["batch_reverts", ["interval_start", "interval_end"]]).all()


def test_generate_panel_ground_truth_round_trips_through_load_ground_truth() -> None:
    """The CSV schema (`interval_start`/`interval_end` as `"YYYYQ#"`,
    `ongoing`, `no_emergence`) must be exactly what
    `ground_truth.load_ground_truth` already parses, unmodified -- that is
    the whole point of the injection design (sec 6): reuse the harness,
    don't build a parallel one."""
    profiles = [Profile("riser", baseline_rate=1.0, rises=True, sustained=True,
                        shape="step", magnitude=8, onset=10, duration=1)]
    counts, denom, gt = generate_panel(profiles, T=30, seed=0)
    path = "/tmp/test_synthetic_gt.csv"
    gt.to_csv(path, index=False)
    loaded = load_ground_truth(path, denom.index.max())
    assert loaded.loc[0, "substance"] == "riser"
    assert loaded.loc[0, "interval_start"] < loaded.loc[0, "interval_end"] \
        or loaded.loc[0, "interval_start"] == loaded.loc[0, "interval_end"]


def test_generate_panel_output_scores_through_the_real_model_table() -> None:
    """End-to-end: a synthetic panel must run through `sweep_eb05` and
    `score_model` exactly as real data does -- no separate scoring path."""
    profiles = [
        Profile("flat", baseline_rate=5.0, rises=False),
        Profile("riser", baseline_rate=1.0, rises=True, sustained=True,
                shape="step", magnitude=8, onset=15, duration=1),
    ]
    counts, denom, gt = generate_panel(profiles, T=30, seed=0)
    path = "/tmp/test_synthetic_gt2.csv"
    gt.to_csv(path, index=False)
    loaded = load_ground_truth(path, denom.index.max())

    from emerging.analysis.trends import (BASELINE_QUARTERS, RECENT_QUARTERS,
                                          sweep_eb05)
    sweep = sweep_eb05(counts, denom, RECENT_QUARTERS, BASELINE_QUARTERS)
    row, per_sub = bm.score_model(sweep, bm.MODELS_BY_ID["eb05"], loaded, 1.5)
    assert row["n_labelled"] == 2
    assert set(per_sub["substance"]) == {"flat", "riser"}


def test_build_sweeps_counts_denom_injection_matches_real_load_path() -> None:
    """The `build_sweeps` injection point (sec 6): a phase-1-safe model list
    must run against `counts_denom` without touching `mentions`/`cutoff` at
    all -- checked by passing `None` for both and confirming no crash."""
    profiles = [Profile("riser", baseline_rate=1.0, rises=True, sustained=True,
                        shape="step", magnitude=8, onset=15, duration=1)]
    counts, denom, gt = generate_panel(profiles, T=30, seed=0)
    safe = [bm.MODELS_BY_ID[i] for i in ["eb05", "eb05-dual", "ears", "nb-trend"]]
    sweeps, denom_out = bm.build_sweeps(
        safe, mentions=None, cutoff=None, recent_quarters=4,
        baseline_quarters=8, quiet=True, counts_denom=(counts, denom))
    assert len(sweeps) == 3  # eb05/eb05-dual share one key; ears, nb_trend separate
    assert denom_out.equals(denom)


def test_nb_trend_confirmed_restores_calibration_on_flat_negatives() -> None:
    """The regression test for `docs/findings/synthetic.md` Finding 1 and its
    fix: on a population with no real trend at all, `nb-trend`'s raw
    "ever alarms across its own as-of sweep" rate is wildly above nominal
    (~33% on the full 200-substance run; this smaller N uses a looser bound
    to stay fast without becoming flaky), and `nb-trend-confirmed` must pull
    it back down sharply. This is the property the fix exists to guarantee --
    if a future change to `_confirmed_flags` or `NB_TREND_CONFIRM_QUARTERS`
    silently regressed it, this is what would catch that."""
    rng = np.random.default_rng(11)
    profiles = [Profile(f"flat_{i}", baseline_rate=float(rng.uniform(1, 30)),
                        rises=False) for i in range(80)]
    counts, denom, gt = generate_panel(profiles, T=45, denom_mode="flat",
                                       denom_n0=400, seed=5)

    raw = bm.MODELS_BY_ID["nb-trend"]
    confirmed = bm.MODELS_BY_ID["nb-trend-confirmed"]
    sweeps, _ = bm.build_sweeps(
        [raw, confirmed], mentions=None, cutoff=None, recent_quarters=4,
        baseline_quarters=8, quiet=True, counts_denom=(counts, denom))

    s_raw = bm.model_score(sweeps[raw.sweep_key], raw)
    s_conf = bm.model_score(sweeps[confirmed.sweep_key], confirmed)
    raw_fp_rate = (sweeps[raw.sweep_key].loc[s_raw > raw.fixed_threshold, "substance"]
                  .nunique() / len(profiles))
    conf_fp_rate = (sweeps[confirmed.sweep_key]
                    .loc[s_conf > confirmed.fixed_threshold, "substance"]
                    .nunique() / len(profiles))

    assert raw_fp_rate > 0.15   # the uncorrected defect, reproduced
    assert conf_fp_rate < 0.10  # the fix pulls it back toward nominal (~2.5%)
    assert conf_fp_rate < raw_fp_rate


@pytest.mark.parametrize("model_id", ["treescan", "eb05-weighted", "eb05-spatial"])
def test_build_sweeps_rejects_real_data_models_with_synthetic_counts(model_id) -> None:
    """TreeScan/`--weighted`/spatial all still read real data through
    `mentions` regardless of `counts_denom` -- must raise, not silently mix a
    synthetic panel with the real county's substances."""
    profiles = [Profile("riser", baseline_rate=1.0, rises=True, sustained=True,
                        shape="step", magnitude=8, onset=15, duration=1)]
    counts, denom, gt = generate_panel(profiles, T=30, seed=0)
    with pytest.raises(ValueError, match="cannot be combined"):
        bm.build_sweeps([bm.MODELS_BY_ID[model_id]], mentions=None, cutoff=None,
                        recent_quarters=4, baseline_quarters=8, quiet=True,
                        counts_denom=(counts, denom))


def test_best_lag_correlation_finds_a_shifted_identical_shape() -> None:
    """Two series that are the same shape but offset in time should match
    at max correlation (~1.0) at the lag that undoes the offset -- the whole
    point of searching lags instead of comparing quarter-for-quarter."""
    rng = np.random.default_rng(0)
    shape = rng.uniform(1, 20, size=30)
    a = np.concatenate([np.zeros(5), shape])       # shape starts at index 5
    b = np.concatenate([shape, np.zeros(5)])       # same shape starts at index 0
    result = _best_lag_correlation(a, b, min_overlap=20)
    assert result is not None
    r, lag, n, mean_b = result
    assert r == pytest.approx(1.0, abs=1e-6)


def test_best_lag_correlation_returns_none_below_min_overlap() -> None:
    """Two short, non-overlapping-enough series must not force a spurious
    match just because *some* lag technically has nonzero overlap."""
    a = np.array([1.0, 5.0, 2.0, 8.0])
    b = np.array([3.0, 1.0, 9.0, 2.0])
    assert _best_lag_correlation(a, b, min_overlap=10) is None


def test_best_lag_correlation_ignores_constant_segments() -> None:
    """A constant window has undefined correlation with anything -- must be
    skipped, not silently returned as if it were a real (anti)correlation."""
    a = np.full(30, 5.0)
    b = np.concatenate([np.zeros(10), np.array([1.0, 9.0, 2.0, 7.0, 3.0,
                                                8.0, 1.0, 6.0, 2.0, 9.0,
                                                4.0, 8.0, 3.0, 7.0, 1.0,
                                                6.0, 2.0, 5.0, 9.0, 4.0])])
    assert _best_lag_correlation(a, b, min_overlap=15) is None


@pytest.mark.skipif(not MENTIONS_PATH.exists(),
                    reason="needs the extracted mentions table, outside "
                           "version control")
def test_nearest_real_matches_returns_one_row_per_profile() -> None:
    """Integration check against the real extract -- every gallery profile
    gets a row, and any match found carries a correlation and a real
    aligned window, not just a substance name with no evidence behind it."""
    counts, _, _ = generate_panel(GALLERY_PROFILES, T=40, seed=0)
    matches = nearest_real_matches(counts, GALLERY_PROFILES)
    assert len(matches) == len(GALLERY_PROFILES)
    assert set(matches["profile"]) == {p.name for p in GALLERY_PROFILES}
    found = matches.dropna(subset=["nearest_real_substance"])
    assert len(found) > 0
    assert (found["correlation"] > 0).all()
    assert (found["window_quarters"] >= 16).all()
