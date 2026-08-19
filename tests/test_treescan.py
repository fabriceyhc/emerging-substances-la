"""
Tests for the tree-temporal scan.

The scan produces p-values, so the property that matters is *calibration*: on
data generated from the null, it must reject at roughly the nominal rate. A
scan statistic that is subtly miscalibrated still returns plausible-looking
tables, and every downstream recurrence interval is then wrong by an unknown
factor. `test_null_calibration` is the test this file exists for; the rest
guard the specific ways the implementation could break quietly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from emerging.analysis import treescan
from emerging.core import tree


def _toy(n_leaves: int = 12, n_q: int = 12, seed: int = 0):
    """A small overlapping node collection over synthetic leaves."""
    rng = np.random.default_rng(seed)
    p = np.full(n_q, 1.0 / n_q)
    # Three overlapping branches plus the leaves and a root, so the
    # correlation structure the maximum has to absorb is actually present.
    membership = np.zeros((n_leaves + 4, n_leaves))
    membership[0, :] = 1.0                       # root
    membership[1, : n_leaves // 2] = 1.0         # branch A
    membership[2, n_leaves // 2:] = 1.0          # branch B
    membership[3, 2:6] = 1.0                     # overlapping family
    membership[4:, :] = np.eye(n_leaves)
    n_leaf = rng.integers(5, 60, size=n_leaves)
    return p, membership, n_leaf


def _draw(rng, n_leaf, p, membership):
    leaf_counts = rng.multinomial(n_leaf, p).astype(float)
    return membership @ leaf_counts, leaf_counts


def test_llr_matches_hand_computation() -> None:
    # n = 100, p = 0.25 so E = 25; observe 40 in the window.
    c, n, e = np.array([40.0]), np.array([100.0]), np.array([25.0])
    expected = 40 * np.log(40 / 25) + 60 * np.log(60 / 75)
    assert treescan._llr(c, n, e) == pytest.approx(expected)


def test_llr_is_zero_without_excess() -> None:
    """One-sided: a deficit is not a signal, and neither is exactly expected."""
    n, e = np.array([100.0, 100.0]), np.array([25.0, 25.0])
    assert treescan._llr(np.array([10.0, 25.0]), n, e) == pytest.approx([0, 0])


def test_llr_handles_the_whole_node_in_the_window() -> None:
    """c == n makes the second term 0*log(0/x); it must be 0, not NaN."""
    out = treescan._llr(np.array([30.0]), np.array([30.0]), np.array([7.5]))
    assert np.isfinite(out).all()
    assert out[0] == pytest.approx(30 * np.log(30 / 7.5))


def test_window_shares_are_cumulative_from_the_end() -> None:
    p = np.array([0.1, 0.2, 0.3, 0.4])
    got = treescan._window_shares(p, 3)
    assert got == pytest.approx([0.4, 0.7, 0.9])


def test_min_cases_is_applied_to_the_null_as_well() -> None:
    """Filtering only the observed side would bias the null downward.

    With a high `min_cases` almost every simulated node-window is masked, so
    the null maximum collapses toward zero. If the mask were applied to the
    observed data only, the p-values would be tiny for everything — this
    asserts the null is masked too, by checking it responds to the threshold.
    """
    p, membership, n_leaf = _toy()
    rng = np.random.default_rng(1)
    node_counts, leaf_counts = _draw(rng, n_leaf, p, membership)

    lo = treescan.scan_once(node_counts, leaf_counts, p, membership,
                            min_cases=0, n_sim=300, seed=7)[3]
    hi = treescan.scan_once(node_counts, leaf_counts, p, membership,
                            min_cases=10_000, n_sim=300, seed=7)[3]
    assert lo.max() > 0
    assert hi.max() == 0.0, "min_cases did not reach the simulated null"


def test_null_calibration() -> None:
    """Type-I error at the nominal level, on data generated from the null.

    Each replicate draws a fresh dataset in which every leaf's cases are
    distributed over time exactly as the reference says, so no node is truly
    elevated. The fraction of replicates whose best node clears the 10% line
    should be about 10%.
    """
    p, membership, n_leaf = _toy(seed=3)
    rng = np.random.default_rng(11)
    alpha, n_rep, n_sim = 0.10, 240, 199

    rejects = 0
    for r in range(n_rep):
        node_counts, leaf_counts = _draw(rng, n_leaf, p, membership)
        llr, _, _, null_max = treescan.scan_once(
            node_counts, leaf_counts, p, membership,
            min_cases=treescan.MIN_CASES, n_sim=n_sim, seed=1000 + r)
        pmin = treescan._pvalue(np.array([llr.max()]), null_max)[0]
        rejects += pmin <= alpha

    rate = rejects / n_rep
    # Binomial sd at p=0.10, n=240 is ~0.019; +/-3.5 sd is a wide enough band
    # that this does not flake, and narrow enough to catch a real bias.
    assert 0.033 < rate < 0.167, (
        f"rejection rate {rate:.3f} at nominal {alpha}: the scan is "
        f"miscalibrated and every recurrence interval is wrong"
    )


def test_scan_is_deterministic_given_a_seed() -> None:
    p, membership, n_leaf = _toy()
    rng = np.random.default_rng(2)
    node_counts, leaf_counts = _draw(rng, n_leaf, p, membership)
    a = treescan.scan_once(node_counts, leaf_counts, p, membership,
                           n_sim=200, seed=42)[3]
    b = treescan.scan_once(node_counts, leaf_counts, p, membership,
                           n_sim=200, seed=42)[3]
    assert np.array_equal(a, b)


def test_dedupe_keeps_the_more_specific_label() -> None:
    """A category holding one substance is the same test as that substance."""
    nodes = [
        tree.Node("root", "root", frozenset({"A", "B"})),
        tree.Node("Solo cat", "category", frozenset({"A"})),
        tree.Node("A", "substance", frozenset({"A"})),
        tree.Node("B", "substance", frozenset({"B"})),
    ]
    kept = {nd.name for nd in tree._dedupe(nodes)}
    assert "A" in kept and "Solo cat" not in kept


def test_leaf_sweep_keeps_only_substance_nodes_and_renames_to_substance(tmp_path) -> None:
    """`benchmark.py` reads this as an eb05-sweep-shaped table -- a `node`
    column named `substance`, with branch/category/root rows (which are not
    individual-substance detections and have no counterpart in the ground
    truth) dropped rather than silently scored against it."""
    import pandas as pd
    df = pd.DataFrame([
        {"as_of": "2024-01-01", "node": "Fentanyl", "level": "substance",
         "recurrence_interval": 500.0, "p_value": 0.002, "llr": 9.0},
        {"as_of": "2024-01-01", "node": "Fentanyl analogs", "level": "family",
         "recurrence_interval": 200.0, "p_value": 0.005, "llr": 7.0},
        {"as_of": "2024-01-01", "node": "Opioids", "level": "superclass",
         "recurrence_interval": 50.0, "p_value": 0.02, "llr": 4.0},
    ])
    path = tmp_path / "backtest.csv"
    df.to_csv(path, index=False)

    sw = treescan.leaf_sweep(path)
    assert list(sw["substance"]) == ["Fentanyl"]
    assert set(sw.columns) == {"substance", "as_of", "recurrence_interval",
                               "p_value", "llr"}


def test_fentanyl_analogs_excludes_fentanyl_itself() -> None:
    """The regex family would otherwise become a proxy for the parent drug."""
    import pandas as pd
    cats = pd.Series(
        {"Fentanyl": "Fentanyl and Fentanyl-related",
         "Acetyl fentanyl": "Fentanyl and Fentanyl-related",
         "para-Fluorofentanyl": "Fentanyl and Fentanyl-related",
         "Morphine": "Narcotic Analgesics"}
    )
    nodes = {nd.name: nd.leaves for nd in tree.build_nodes(cats)}
    assert "Fentanyl" not in nodes["Fentanyl analogs"]
    assert {"Acetyl fentanyl", "para-Fluorofentanyl"} == nodes["Fentanyl analogs"]


def _branch_fixture(tmp_path, rows):
    """A 12-quarter frame with one family node over three substances.

    Reference is flat, so a 4-quarter window expects exactly a third of each
    substance's total. `A` puts all 12 of its cases in that window (excess
    +8), `B` is flat at 3/quarter (excess 0, present), and `C` appears only
    in the very first quarter (excess -4/3, absent from the window). The
    family's excess is 6.67, so the shares are A 1.20, B 0.00, C -0.20 --
    one substance above any guard, one below it but present, one absent.
    """
    import pandas as pd
    from emerging.core.tree import Node

    qs = pd.date_range("2022-01-01", periods=12, freq="QS")
    counts = []
    for i, q in enumerate(qs):
        counts.append({"substance": "A", "quarter": q, "n": 3 if i >= 8 else 0})
        counts.append({"substance": "B", "quarter": q, "n": 3})
        counts.append({"substance": "C", "quarter": q, "n": 4 if i == 0 else 0})
    counts = pd.DataFrame(counts)
    ref = pd.Series(100.0, index=qs)
    nodes = [
        Node("Fam", "family", frozenset({"A", "B", "C"})),
        Node("Sub", "family", frozenset({"A", "B"})),
        Node("A", "substance", frozenset({"A"})),
        Node("B", "substance", frozenset({"B"})),
        Node("C", "substance", frozenset({"C"})),
    ]
    path = tmp_path / "backtest.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path, (counts, ref, ref, nodes)


def test_branch_sweep_credits_a_family_to_the_members_actually_in_it(tmp_path) -> None:
    """The presence guard is what stops a branch signal from being handed to
    every substance filed under it. `C` is a member of the node that fired but
    has no case inside the window that fired, so it gets nothing -- the
    vectorized form of `backtest`'s "a branch cannot have detected a substance
    that did not exist yet". `B` is present but contributes none of the
    excess, so it survives unguarded and is dropped by the attribution guard.
    """
    path, prep = _branch_fixture(tmp_path, [
        {"as_of": "2024-10-01", "node": "Fam", "level": "family", "n_leaves": 3,
         "window_q": 4, "observed": 24, "expected": 17.3, "n_study": 52,
         "llr": 6.0, "p_value": 0.002, "obs_exp": 1.4,
         "recurrence_interval": 500.0},
    ])
    unguarded = treescan.branch_sweep(path, prepared=prep)
    assert list(unguarded["substance"]) == ["A", "B"]
    assert list(unguarded["source_node"]) == ["Fam", "Fam"]
    assert unguarded.set_index("substance")["excess_share"]["A"] == pytest.approx(1.2)
    assert unguarded.set_index("substance")["excess_share"]["B"] == pytest.approx(0.0)

    guarded = treescan.branch_sweep(path, prepared=prep, min_excess_share=0.15)
    assert list(guarded["substance"]) == ["A"]
    assert list(guarded["recurrence_interval"]) == [500.0]


def test_branch_sweep_share_matches_excess_share(tmp_path) -> None:
    """`branch_sweep` computes every substance's attribution in one array
    pass; `excess_share` computes one at a time from the counts table. They
    are the same quantity and must not drift apart -- the guard's threshold
    is calibrated against the numbers `backtest` prints, which come from
    `excess_share`."""
    path, prep = _branch_fixture(tmp_path, [
        {"as_of": "2024-10-01", "node": "Fam", "level": "family", "n_leaves": 3,
         "window_q": 4, "observed": 24, "expected": 17.3, "n_study": 52,
         "llr": 6.0, "p_value": 0.002, "obs_exp": 1.4,
         "recurrence_interval": 500.0},
    ])
    counts, _, ref, nodes = prep
    got = treescan.branch_sweep(path, prepared=prep).set_index("substance")
    for s in got.index:
        assert got.loc[s, "excess_share"] == pytest.approx(
            treescan.excess_share(counts, ref, {"A", "B", "C"}, s,
                                  pd.Timestamp("2024-10-01"), 4))


def test_branch_sweep_leaf_wins_on_score_and_ties_go_to_the_smaller_node(tmp_path) -> None:
    """With `include_leaves` the whole hierarchy is one detector, so a
    substance takes the best score available to it and a leaf is exempt from
    both guards (it is trivially present in, and wholly responsible for, its
    own excess). Equal scores resolve to the most specific node, so
    `source_node` names the tightest branch that justifies the alarm rather
    than whichever one sorted first."""
    path, prep = _branch_fixture(tmp_path, [
        {"as_of": "2024-10-01", "node": "Fam", "level": "family", "n_leaves": 3,
         "window_q": 4, "observed": 24, "expected": 17.3, "n_study": 52,
         "llr": 6.0, "p_value": 0.002, "obs_exp": 1.4,
         "recurrence_interval": 500.0},
        {"as_of": "2024-10-01", "node": "Sub", "level": "family", "n_leaves": 2,
         "window_q": 4, "observed": 24, "expected": 16.0, "n_study": 48,
         "llr": 6.0, "p_value": 0.002, "obs_exp": 1.5,
         "recurrence_interval": 500.0},
        {"as_of": "2024-10-01", "node": "B", "level": "substance", "n_leaves": 1,
         "window_q": 4, "observed": 12, "expected": 12.0, "n_study": 36,
         "llr": 8.0, "p_value": 0.0005, "obs_exp": 1.0,
         "recurrence_interval": 2000.0},
    ])
    out = treescan.branch_sweep(path, prepared=prep, min_excess_share=0.15,
                                include_leaves=True).set_index("substance")
    # `B` contributes 0% of both branches' excess and is guarded out of them,
    # but its own leaf row stands on its own.
    assert out.loc["B", "recurrence_interval"] == 2000.0
    assert out.loc["B", "source_level"] == "substance"
    # `A` clears the guard on both branches at the same RI; the two-leaf node
    # is the more specific claim.
    assert out.loc["A", "source_node"] == "Sub"
    assert "C" not in out.index
