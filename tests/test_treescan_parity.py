"""
Assert that our scan still agrees with the TreeScan reference binary.

`tests/test_treescan.py` shows the implementation is internally consistent and
correctly calibrated. Neither of those would catch computing the *wrong
statistic* consistently. This file is the external check: same tree, same
counts, same settings, against Martin Kulldorff's TreeScan v2.4.

The binary is not in the repository — it is distributed under its own licence
from treescan.org, and `.gitignore` keeps it out. These tests skip when it is
absent, so they are a no-op on a machine that only has this repository. Run
them wherever the binary has been extracted to `treescan/`:

    pytest tests/test_treescan_parity.py -v

The comparison runs our scan with a *uniform* reference series, which is the
one configuration where the two models coincide — TreeScan's tree-temporal
model spreads each node's cases uniformly over the study period, while ours
takes an arbitrary reference series. See `emerging/treescan_validate.py` for
the full argument and the matched settings.
"""

from __future__ import annotations

import pandas as pd
import pytest

from emerging import treescan_validate as tv
from emerging.extract import MENTIONS_PATH
from emerging.treescan import run_scan

pytestmark = pytest.mark.skipif(
    not tv.TREESCAN_BIN.exists() or not MENTIONS_PATH.exists(),
    reason=f"needs the TreeScan binary at {tv.TREESCAN_BIN} and the extracted "
           f"mentions table; both are outside version control",
)

STUDY_Q, MAX_W = 12, 6


@pytest.fixture(scope="module", params=["strict", "dag"])
def paired(request, tmp_path_factory):
    """Run both implementations once per tree mode."""
    mode = request.param
    work = tmp_path_factory.mktemp(f"treescan_{mode}")
    counts, as_of, quarters, wide, sel = tv._setup(MENTIONS_PATH, mode, STUDY_Q)
    tre, cas, ids = tv.write_inputs(work, sel, wide, quarters)
    back = {v: k for k, v in ids.items()}

    # 0 replications: TreeScan lists every cut it can form, which is what makes
    # the deterministic comparison cover the whole node set rather than the
    # handful that reach significance.
    ref, _, n_nodes = tv.run_reference(work, tre, cas, len(quarters), MAX_W, 0,
                                       mode, 12345678, 1, tag="det")
    ref["node"] = ref["node"].map(back)
    uniform = pd.Series(1.0, index=pd.Index(quarters, name="quarter"))
    ours = run_scan(counts, uniform, sel, as_of, STUDY_Q, MAX_W,
                    min_cases=tv.MIN_CUT_CASES, n_sim=99, seed=12345678)
    merged = ours.merge(ref, on="node", how="inner", suffixes=("_ours", "_ref"))
    return {"mode": mode, "merged": merged, "ref": ref, "ours": ours,
            "n_nodes_ref": n_nodes, "n_nodes_ours": len(sel), "work": work,
            "tre": tre, "cas": cas, "quarters": quarters, "counts": counts,
            "sel": sel, "as_of": as_of}


def test_same_search_space(paired) -> None:
    """Both must scan the same nodes, or the multiplicity correction differs."""
    assert paired["n_nodes_ref"] == paired["n_nodes_ours"]


def test_every_treescan_cut_is_reproduced(paired) -> None:
    missing = set(paired["ref"]["node"]) - set(paired["ours"]["node"])
    assert not missing, f"TreeScan forms cuts we never consider: {sorted(missing)[:5]}"
    assert len(paired["merged"]) >= 40, "comparison covers implausibly few nodes"


def test_log_likelihood_ratios_match(paired) -> None:
    """The deterministic core. Any real difference shows up here."""
    m = paired["merged"]
    diff = (m["llr_ours"] - m["llr_ref"]).abs()
    worst = m.loc[diff.idxmax()]
    assert diff.max() < tv.LLR_TOL, (
        f"[{paired['mode']}] largest LLR gap {diff.max():.2e} at "
        f"{worst['node']}: ours {worst['llr_ours']:.6f} vs TreeScan "
        f"{worst['llr_ref']:.6f}"
    )


def test_best_window_and_counts_match(paired) -> None:
    m = paired["merged"]
    bad = m[m["window_q_ours"] != m["window_q_ref"]]
    assert bad.empty, f"different best window for {list(bad['node'])[:5]}"
    bad = m[m["observed_ours"] != m["observed_ref"]]
    assert bad.empty, f"different window counts for {list(bad['node'])[:5]}"


def test_expected_counts_match(paired) -> None:
    m = paired["merged"]
    diff = (m["expected_ours"] - m["expected_ref"]).abs()
    assert diff.max() < tv.EXPECTED_TOL, (
        f"expected counts differ by {diff.max():.2e}; TreeScan prints 2 dp, "
        f"so anything above {tv.EXPECTED_TOL} is a real disagreement")


@pytest.mark.slow
def test_pvalues_agree_within_monte_carlo_error(paired) -> None:
    """Simulated p-values, which agree only up to simulation noise."""
    p = paired
    ref, _, _ = tv.run_reference(p["work"], p["tre"], p["cas"],
                                 len(p["quarters"]), MAX_W, 9999, p["mode"],
                                 12345678, 1, tag="mc")
    nm = pd.read_csv(p["work"] / "names.csv")
    ref["node"] = ref["node"].map(dict(zip(nm["id"], nm["name"])))
    uniform = pd.Series(1.0, index=pd.Index(p["quarters"], name="quarter"))
    ours = run_scan(p["counts"], uniform, p["sel"], p["as_of"], STUDY_Q, MAX_W,
                    min_cases=tv.MIN_CUT_CASES, n_sim=9999, seed=12345678)
    m = ours.merge(ref, on="node", how="inner", suffixes=("_ours", "_ref"))
    diff = (m["p_value_ours"] - m["p_value_ref"]).abs()
    # Two independent sets of 9,999 replicates; the sd of a p-value difference
    # near 0.5 is about 0.007, so 0.05 is roughly 7 sd and will not flake.
    assert diff.max() < 0.05, (
        f"p-values differ by up to {diff.max():.4f}, which is too much for "
        f"simulation error — the null distributions disagree")
    assert (((m["p_value_ours"] <= 0.01) == (m["p_value_ref"] <= 0.01)).all())
