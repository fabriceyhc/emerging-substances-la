"""
Tests for `emerging/validation/ground_truth.py`: the interval-overlap metric
and the reference-CSV loader.

The metric itself has to be right independent of any real data -- these are
all synthetic interval sets with a hand-computed answer.
"""

from __future__ import annotations

import pandas as pd
import pytest

from emerging.validation import ground_truth as gt


def _q(s: str) -> pd.Timestamp:
    return pd.Period(s, freq="Q").to_timestamp()


def test_interval_quarters_is_inclusive_of_both_ends() -> None:
    assert gt._interval_quarters(_q("2024Q1"), _q("2024Q1")) == {2024 * 4 + 0}
    assert gt._interval_quarters(_q("2024Q1"), _q("2024Q3")) == {
        2024 * 4 + 0, 2024 * 4 + 1, 2024 * 4 + 2}


def test_overlap_is_one_for_identical_single_intervals() -> None:
    iv = [(_q("2023Q1"), _q("2023Q4"))]
    m = gt.interval_overlap(iv, iv)
    assert m["iou"] == pytest.approx(1.0)
    assert m["recall"] == pytest.approx(1.0)
    assert m["precision"] == pytest.approx(1.0)


def test_overlap_is_zero_for_disjoint_intervals() -> None:
    a = [(_q("2020Q1"), _q("2020Q4"))]
    b = [(_q("2023Q1"), _q("2023Q4"))]
    m = gt.interval_overlap(a, b)
    assert m["iou"] == pytest.approx(0.0)
    assert m["recall"] == pytest.approx(0.0)
    assert m["precision"] == pytest.approx(0.0)


def test_overlap_partial_matches_hand_computation() -> None:
    # GT: 2023Q1-Q4 (4 quarters). Detected: 2023Q3-2024Q2 (4 quarters).
    # Overlap: 2023Q3-Q4 (2 quarters). Union: 2023Q1-2024Q2 (6 quarters).
    gt_iv = [(_q("2023Q1"), _q("2023Q4"))]
    det_iv = [(_q("2023Q3"), _q("2024Q2"))]
    m = gt.interval_overlap(gt_iv, det_iv)
    assert m["overlap_quarters"] == 2
    assert m["gt_quarters"] == 4
    assert m["det_quarters"] == 4
    assert m["iou"] == pytest.approx(2 / 6)
    assert m["recall"] == pytest.approx(2 / 4)      # half the true window caught
    assert m["precision"] == pytest.approx(2 / 4)    # half the alarm time was real


def test_overlap_unions_multiple_intervals_per_side() -> None:
    """A substance with two separate ground-truth episodes and one detected
    episode spanning both should score on the union, not just one pairing."""
    gt_iv = [(_q("2020Q1"), _q("2020Q2")), (_q("2022Q1"), _q("2022Q2"))]
    det_iv = [(_q("2020Q1"), _q("2022Q2"))]  # one long span covering both
    m = gt.interval_overlap(gt_iv, det_iv)
    assert m["gt_quarters"] == 4    # 2 + 2, not 4+4 double counted
    assert m["overlap_quarters"] == 4
    assert m["recall"] == pytest.approx(1.0)


def test_overlap_handles_empty_detected_side() -> None:
    """A substance the detector never flagged at all -- recall 0, precision
    undefined (nothing to be precise about), not a crash."""
    gt_iv = [(_q("2023Q1"), _q("2023Q4"))]
    m = gt.interval_overlap(gt_iv, [])
    assert m["recall"] == pytest.approx(0.0)
    assert m["det_quarters"] == 0
    import math
    assert math.isnan(m["precision"])


def test_load_ground_truth_resolves_ongoing_to_last_q(tmp_path) -> None:
    p = tmp_path / "gt.csv"
    p.write_text(
        "substance,interval_start,interval_end,ongoing,provenance,confidence,"
        "national_context,notes\n"
        "Foo,2023Q1,,True,manual,high,,test row\n"
        "Bar,2022Q1,2022Q4,False,manual,high,,test row\n"
    )
    last_q = _q("2025Q4")
    df = gt.load_ground_truth(p, last_q=last_q)
    foo = df.set_index("substance").loc["Foo"]
    bar = df.set_index("substance").loc["Bar"]
    assert foo["interval_end"] == last_q
    assert bar["interval_end"] == _q("2022Q4")


def test_load_ground_truth_raises_without_last_q_when_ongoing(tmp_path) -> None:
    p = tmp_path / "gt.csv"
    p.write_text(
        "substance,interval_start,interval_end,ongoing,provenance,confidence,"
        "national_context,notes\n"
        "Foo,2023Q1,,True,manual,high,,test row\n"
    )
    with pytest.raises(ValueError, match="ongoing"):
        gt.load_ground_truth(p, last_q=None)


# --------------------------------------------------------------------------
# NFLIS onset -- the sparse-row reindexing bug, and the mitragynine trap
# --------------------------------------------------------------------------

def _halfyear_index(start_year: int, end_year: int) -> list[pd.Timestamp]:
    out = []
    for y in range(start_year, end_year + 1):
        out += [pd.Timestamp(f"{y}-01-01"), pd.Timestamp(f"{y}-07-01")]
    return out


def test_nflis_onset_finds_a_clean_sustained_rise() -> None:
    idx = _halfyear_index(2020, 2023)
    vals = [0, 0, 0, 1, 11, 55, 105, 136]  # onset at index 4 -> 2022-01-01
    s = pd.Series(vals, index=idx)
    assert gt.nflis_onset(s) == pd.Timestamp("2022-01-01")


def test_nflis_onset_requires_two_consecutive_low_periods_first() -> None:
    idx = _halfyear_index(2020, 2022)
    # index 1 (11) clears the floor but has only one preceding low period
    # (there's no index -1) so it can never qualify as an onset regardless;
    # index 2 (8) is itself below the floor; the real onset -- floor-clearing
    # with *both* of the two periods right before it below the floor -- is
    # index 4.
    vals = [9, 11, 8, 8, 12, 30]
    s = pd.Series(vals, index=idx)
    assert gt.nflis_onset(s) == idx[4]


def test_nflis_onset_reindexing_matters_the_way_the_docstring_says() -> None:
    """The bug this function's docstring warns about, reproduced directly:
    a sparse series (the raw DQS export's own shape -- zero-count periods
    are omitted as rows, not written as 0) makes two periods that are years
    apart look adjacent unless reindexed onto a full calendar first."""
    sparse = pd.Series(
        [3, 40],
        index=[pd.Timestamp("2019-01-01"), pd.Timestamp("2020-07-01")],
    )
    # Unreindexed: only two data points, i<2 never holds, so nothing is ever
    # checked -- this is the "silently wrong" failure mode, not a crash.
    assert gt.nflis_onset(sparse) is None

    full = sparse.reindex(_halfyear_index(2019, 2021), fill_value=0)
    # Reindexed: 2019-07-01 and 2020-01-01 are genuine zeros in between, so
    # the >=2-consecutive-low-periods check can actually see them.
    assert gt.nflis_onset(full) == pd.Timestamp("2020-07-01")


def test_nflis_onset_does_not_pretend_to_reject_an_unsustained_spike() -> None:
    """Documented limitation, not a bug: a single spike that collapses right
    back down (mitragynine's real 2015H1 pattern: 33 then 2) still reads as
    an 'onset' by this rule. The function's own docstring says this must be
    read against the series, not trusted blind -- this test pins down that
    the false positive is real, so nobody 'fixes' it by accident later
    without updating the docstring's claim too."""
    idx = _halfyear_index(2013, 2016)
    vals = [1, 1, 1, 1, 33, 2, 1, 4]
    s = pd.Series(vals, index=idx)
    onset = gt.nflis_onset(s)
    assert onset == pd.Timestamp("2015-01-01")   # the spike, not a real trend
    assert s.loc[onset + pd.DateOffset(months=6)] < gt.NFLIS_ONSET_FLOOR  # collapses right after


# --------------------------------------------------------------------------
# no_emergence -- Morphine, the first genuine negative case
# --------------------------------------------------------------------------

def test_overlap_of_empty_ground_truth_against_a_real_detection_is_pure_false_alarm() -> None:
    """The exact shape of the Morphine case: a substance with zero
    ground-truth-positive time that the sweep nonetheless alarmed on. This
    must read as recall undefined (nothing to have covered) and precision
    0.0 (every alarmed quarter was, by definition, a false one) -- the
    reverse of test_overlap_handles_empty_detected_side."""
    det_iv = [(_q("2025Q2"), _q("2025Q2"))]
    m = gt.interval_overlap([], det_iv)
    import math
    assert math.isnan(m["recall"])
    assert m["precision"] == pytest.approx(0.0)
    assert m["gt_quarters"] == 0
    assert m["det_quarters"] == 1


# --------------------------------------------------------------------------
# Wilson CI and the substance-level false-positive rate
# --------------------------------------------------------------------------

def test_wilson_ci_no_data_is_maximum_uncertainty() -> None:
    assert gt.wilson_ci(0, 0) == (0.0, 1.0)


def test_wilson_ci_contains_the_point_estimate() -> None:
    for k, n in [(0, 5), (5, 5), (3, 5), (14, 17), (1, 1)]:
        lo, hi = gt.wilson_ci(k, n)
        p = k / n
        assert lo <= p <= hi


def test_wilson_ci_narrows_as_n_grows_at_the_same_rate() -> None:
    lo5, hi5 = gt.wilson_ci(3, 5)
    lo50, hi50 = gt.wilson_ci(30, 50)
    assert (hi50 - lo50) < (hi5 - lo5)


def test_wilson_ci_never_exceeds_the_unit_interval() -> None:
    lo, hi = gt.wilson_ci(0, 3)
    assert lo == pytest.approx(0.0)
    lo, hi = gt.wilson_ci(3, 3)
    assert hi <= 1.0


def test_wilson_ci_matches_a_hand_computation() -> None:
    """k=14, n=17 -- this repository's own plain-EB05 substance-level
    false-positive rate at the time this was written. Hand-computed via the
    textbook Wilson formula, independent of the implementation."""
    lo, hi = gt.wilson_ci(14, 17)
    assert lo == pytest.approx(0.5897, abs=1e-3)
    assert hi == pytest.approx(0.9380, abs=1e-3)


def test_precision_rate_counts_a_partial_hit_as_a_success() -> None:
    """Benzylfentanyl-shaped case: precision 0.33 (some overlap, not total)
    still counts as a substance-level success -- this metric answers 'did
    the alarm ever touch something real', not 'how much of it was real'."""
    # sub A: partial hit (success). sub B: never alarmed (excluded from n).
    # sub C: pure false alarm (failure).
    per_sub = pd.DataFrame({
        "det_quarters": [3, 0, 4],
        "overlap_quarters": [1, 0, 0],
    })
    r = gt.precision_rate(per_sub)
    assert r["precision_rate_n"] == 2       # only substances with det_quarters > 0
    assert r["precision_rate_k"] == 1       # the partial hit counts as a success
    assert r["precision_rate"] == pytest.approx(0.5)


def test_load_ground_truth_carries_no_emergence_column_through(tmp_path) -> None:
    p = tmp_path / "gt.csv"
    p.write_text(
        "substance,interval_start,interval_end,ongoing,provenance,confidence,"
        "national_context,notes,no_emergence\n"
        "RealOne,2023Q1,2023Q4,False,manual,high,,test row,False\n"
        "FalseAlarm,2024Q1,2024Q1,False,manual,high,,test row,True\n"
    )
    df = gt.load_ground_truth(p, last_q=_q("2025Q4"))
    by_sub = df.set_index("substance")
    assert bool(by_sub.loc["FalseAlarm", "no_emergence"]) is True
    assert bool(by_sub.loc["RealOne", "no_emergence"]) is False
