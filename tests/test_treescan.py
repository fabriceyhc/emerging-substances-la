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
