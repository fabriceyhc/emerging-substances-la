"""
Tree-temporal scan statistic over the substance hierarchy (game plan P1).

Replaces one-substance-at-a-time EB05 testing with a scan over every node of
the hierarchy in `tree.py` and every recent time window, with a Monte Carlo
p-value that is adjusted for all of it at once.

**What it fixes, and what it does not.**

  - *Aggregation power.* The nitazenes are individually unscoreable (1 and 2
    mentions) and previously needed a hand-coded family workaround in
    `trends watch`. The scan tests the branch, so the family is a first-class
    hypothesis instead of a special case.
  - *Multiple testing.* EB05 ranks 205 substances across 45 as-of quarters with
    no adjustment at all, and the 1.5 line is read off a backtest rather than
    derived. Here the null distribution is the *maximum* LLR over every node
    and every window, so a p-value already accounts for the whole search.
  - *Window length.* EB05 fixes the recent window at 4 quarters, and the
    research log records that this matters a great deal — lidocaine scores 1.04
    on 2025Q1–Q4, 2.20 on 2024Q4–2025Q3 and peaks at 4.29 across the sweep.
    Scanning window length 1..6 removes the choice, and the multiplicity
    correction charges us for having made it.
  - *Not fixed: repeated looks.* The recurrence interval here is the
    multiplicity across nodes and windows **within one analysis**. Running the
    scan every quarter is a further multiplicity it does not cover. The game
    plan called MaxSPRT "redundant with TreeScan's recurrence interval"; that
    is only true of the within-analysis part. Read `RI` as "how often a signal
    this strong appears in a single quarter's scan by chance", and divide by
    the number of quarters you have scanned if you want the sequential rate.

**The statistic.** Condition on each substance's total over the study period
and ask whether those cases pile up in a recent window more than the reference
time series says they should. For node G (a set of substances) and window W:

    n = cases of G over the whole study period
    p = reference share of the study period falling in W
    E = n * p                        expected under the null
    c = observed cases of G in W

    LLR = c*log(c/E) + (n-c)*log((n-c)/(n-E))      when c > E, else 0

This is the conditional Bernoulli form used by TreeScan's tree-temporal model.
Under the null the test statistic is max LLR over all (G, W); it is simulated
by redistributing each *substance's* cases across quarters as
Multinomial(n_s, p). Conditioning on the leaf totals automatically conditions
on every internal node's total, and simulating at the leaf preserves the
correlation between a node and its parents — which is exactly what makes the
maximum the right adjustment for overlapping nodes.

**The reference time series** is the choice that decides what "rising" means:

    --reference deaths    p from quarterly OD deaths. Direct analogue of EB05:
                          rising = a growing share of overdose deaths.
    --reference mentions  p from quarterly total substance mentions. The
                          self-controlled form. Immune to panel drift, because
                          a quarter in which the ME reports more substances per
                          death also gets a larger p.

The two disagree exactly when the denominator is moving, which the research log
records as the single most common source of false findings here. Running both
is the point, not a hedge.

**Not the TreeScan binary.** This is a from-scratch implementation of the
method, in Python, so the pipeline stays a single language. Two deliberate
departures: nodes may overlap (`tree.py` explains why), and the study period
rolls with the as-of quarter rather than being fixed, so the comparison against
EB05's rolling 12-quarter frame is like-for-like.

Usage:

    emerging treescan scan
    emerging treescan scan --reference mentions
    emerging treescan backtest
    emerging treescan plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import results_dir
from emerging.core.tree import Node, build
from emerging.config import CUTOFF
# ALARM_THRESHOLD and KNOWN_EMERGENCES are EB05 concepts and the head-to-head
# genuinely depends on them -- this import is the real thing, not plumbing.
from emerging.analysis.trends import (ALARM_THRESHOLD, KNOWN_EMERGENCES,
                             load_quarterly, rank_substances)
from emerging.viz import AQUA, BLUE, INK, INK_2, MUTED, ORANGE, YELLOW

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("treescan")


# Conditioning frame, in quarters. 12 = EB05's 4 recent + 8 baseline, so the
# two methods see the same data and the head-to-head is honest.
STUDY_QUARTERS = 12
# Longest window scanned, in quarters. Must leave a baseline: at 6 of 12 the
# shortest possible comparison period is half the frame.
MAX_WINDOW = 6
# Minimum observed cases in the window for a node to be scannable. Applied
# identically to the observed data and to every simulation — filtering only
# the observed side would bias the null downward and manufacture significance.
MIN_CASES = 3
N_SIM = 9999
# Signal line, as a recurrence interval in analyses: RI 100 means a signal this
# strong turns up in one scan out of 100 by chance (p = 0.01). Quarterly, that
# is roughly once per 25 years of scanning — but see the caveat above about
# repeated looks.
RI_THRESHOLD = 100.0
# Floor below which an LLR is floating-point noise rather than an excess.
# Real signals here are >= 1e-3; see `_llr`.
_LLR_EPS = 1e-9

SCAN_PATH = RESULTS_DIR / "signals.csv"
BACKTEST_PATH = RESULTS_DIR / "backtest.csv"
HEADTOHEAD_PATH = RESULTS_DIR / "vs_eb05.csv"

LEVEL_COLOR = {"root": MUTED, "superclass": INK_2, "category": AQUA,
               "family": ORANGE, "substance": BLUE}


def _llr(c: np.ndarray, n: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Conditional Bernoulli log-likelihood ratio, one-sided (excess only).

    Zero wherever c <= e, so the scan only ever reports elevation. Declines are
    a real question but a different one, and EB95 in `trends rank` already
    answers it.
    """
    c, n, e = np.broadcast_arrays(c, n, e)
    rest_c, rest_e = n - c, n - e
    # e == 0 only when the node has no cases at all, in which case c == 0 too;
    # the `out=ones` default makes log(1) = 0 there rather than leaving the
    # uninitialised memory that a bare `where=` would.
    ratio1 = np.divide(c, e, out=np.ones_like(e), where=e > 0)
    ratio2 = np.divide(rest_c, rest_e, out=np.ones_like(e), where=rest_e > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        term1 = np.where(c > 0, c * np.log(ratio1), 0.0)
        term2 = np.where(rest_c > 0, rest_c * np.log(ratio2), 0.0)
    out = term1 + term2
    # The `out > _LLR_EPS` clause is not cosmetic. When c equals e exactly --
    # a node whose window count lands on its expectation, e.g. metoprolol at
    # 3 observed and 3.0 expected -- `c > e` can still be True on the last bit
    # and the LLR comes out at 3e-16. That never wins a maximum, but it puts a
    # node with literally no excess into the "positive LLR" list, which is how
    # the comparison against the TreeScan binary first surfaced it.
    return np.where((c > e) & (out > _LLR_EPS) & np.isfinite(out), out, 0.0)


def _membership(nodes: list[Node], leaves: list[str]) -> np.ndarray:
    """(n_nodes, n_leaves) 0/1 matrix."""
    idx = {s: i for i, s in enumerate(leaves)}
    m = np.zeros((len(nodes), len(leaves)), dtype=np.float64)
    for r, nd in enumerate(nodes):
        for s in nd.leaves:
            j = idx.get(s)
            if j is not None:
                m[r, j] = 1.0
    return m


def _window_shares(p: np.ndarray, max_window: int) -> np.ndarray:
    """Reference share of the study period in each trailing window, length 1..W."""
    rev = np.cumsum(p[::-1])
    return rev[:max_window]


def _trailing_counts(x: np.ndarray, max_window: int) -> np.ndarray:
    """Trailing-window sums of a (..., n_q) array -> (..., max_window)."""
    rev = np.cumsum(x[..., ::-1], axis=-1)
    return rev[..., :max_window]


def scan_once(
    node_counts: np.ndarray,
    leaf_counts: np.ndarray,
    p: np.ndarray,
    membership: np.ndarray,
    max_window: int = MAX_WINDOW,
    min_cases: int = MIN_CASES,
    n_sim: int = N_SIM,
    seed: int = 20260810,
    batch: int = 500,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run the scan and its Monte Carlo null.

    Returns (llr, obs, expected, null_max) where the first three are
    (n_nodes, max_window) and null_max is (n_sim,).
    """
    n_g = node_counts.sum(axis=1)                       # (n_nodes,)
    pw = _window_shares(p, max_window)                  # (max_window,)
    obs = _trailing_counts(node_counts, max_window)     # (n_nodes, W)
    exp = n_g[:, None] * pw[None, :]
    llr = _llr(obs, n_g[:, None], exp)
    llr = np.where(obs >= min_cases, llr, 0.0)

    rng = np.random.default_rng(seed)
    n_leaf = leaf_counts.sum(axis=1).astype(np.int64)
    null_max = np.empty(n_sim, dtype=np.float64)
    done = 0
    while done < n_sim:
        b = min(batch, n_sim - done)
        # (b, n_leaves, n_q): each substance's cases redistributed over the
        # study period according to the reference time series.
        sim_leaf = rng.multinomial(np.tile(n_leaf, (b, 1)), p)
        sim_node = np.einsum("nl,blq->bnq", membership, sim_leaf.astype(float))
        s_obs = _trailing_counts(sim_node, max_window)
        s_llr = _llr(s_obs, n_g[None, :, None], exp[None, :, :])
        s_llr = np.where(s_obs >= min_cases, s_llr, 0.0)
        null_max[done:done + b] = s_llr.reshape(b, -1).max(axis=1)
        done += b
    return llr, obs, exp, null_max


def _pvalue(stat: np.ndarray, null_max: np.ndarray) -> np.ndarray:
    """Monte Carlo p from the null distribution of the maximum."""
    n_sim = len(null_max)
    ge = (null_max[None, :] >= stat[:, None]).sum(axis=1)
    return (1.0 + ge) / (n_sim + 1.0)


def prepare(
    mentions_path: Path = MENTIONS_PATH,
    cutoff: str | None = CUTOFF,
    reference: str = "deaths",
) -> tuple[pd.DataFrame, pd.Series, pd.Series, list[Node]]:
    """Load counts, the reference series and the node collection."""
    counts, denom = load_quarterly(mentions_path, cutoff=cutoff)
    nodes, _ = build(mentions_path)
    if reference == "deaths":
        ref = denom.astype(float)
    elif reference == "mentions":
        ref = counts.groupby("quarter")["n"].sum().reindex(
            denom.index, fill_value=0).astype(float)
    else:
        raise ValueError("reference must be 'deaths' or 'mentions'")
    return counts, denom, ref, nodes


def run_scan(
    counts: pd.DataFrame,
    ref: pd.Series,
    nodes: list[Node],
    as_of: pd.Timestamp,
    study_quarters: int = STUDY_QUARTERS,
    max_window: int = MAX_WINDOW,
    min_cases: int = MIN_CASES,
    n_sim: int = N_SIM,
    seed: int = 20260810,
) -> pd.DataFrame:
    """Scan as of one quarter. Returns one row per (node, best window)."""
    quarters = [q for q in sorted(ref.index) if q <= as_of][-study_quarters:]
    if len(quarters) < study_quarters:
        raise ValueError(f"need {study_quarters} quarters, have {len(quarters)}")

    p = ref.loc[quarters].to_numpy(dtype=float)
    if p.sum() <= 0:
        raise ValueError("reference series is empty over the study period")
    p = p / p.sum()

    sub = counts[counts["quarter"].isin(quarters)]
    wide = (sub.pivot_table(index="substance", columns="quarter", values="n",
                            aggfunc="sum", fill_value=0)
              .reindex(columns=quarters, fill_value=0))
    leaves = list(wide.index)
    leaf_counts = wide.to_numpy(dtype=float)

    # Only nodes with at least one leaf present in this window are testable.
    live = [nd for nd in nodes if nd.leaves & set(leaves)]
    membership = _membership(live, leaves)
    node_counts = membership @ leaf_counts

    llr, obs, exp, null_max = scan_once(
        node_counts, leaf_counts, p, membership, max_window=max_window,
        min_cases=min_cases, n_sim=n_sim, seed=seed)

    best_w = llr.argmax(axis=1)
    rows = np.arange(len(live))
    best_llr = llr[rows, best_w]
    pvals = _pvalue(best_llr, null_max)

    out = pd.DataFrame({
        "as_of": as_of,
        "node": [nd.name for nd in live],
        "level": [nd.level for nd in live],
        "n_leaves": [nd.n_leaves for nd in live],
        "window_q": best_w + 1,
        "observed": obs[rows, best_w].astype(int),
        "expected": exp[rows, best_w],
        "n_study": node_counts.sum(axis=1).astype(int),
        "llr": best_llr,
        "p_value": pvals,
    })
    out["obs_exp"] = out["observed"] / out["expected"].replace(0, np.nan)
    out["recurrence_interval"] = 1.0 / out["p_value"]
    out.loc[out["llr"] <= 0, ["p_value", "recurrence_interval", "window_q"]] = \
        [1.0, 1.0, np.nan]
    return out.sort_values("llr", ascending=False).reset_index(drop=True)


def leaf_sweep(path: Path = BACKTEST_PATH) -> pd.DataFrame:
    """Substance-leaf rows from a cached `treescan backtest` run, reshaped to
    the (substance, as_of, score) shape `emerging.validation.benchmark` reads
    every other detector's sweep from.

    Re-running the scan here (45 as-of quarters x N_SIM Monte Carlo replicates
    each) is the expensive part of `backtest`, and it already computes every
    node including the substance leaves -- `backtest` just doesn't write out
    the zero-LLR rows. So this reads the cache rather than re-scanning; if the
    mentions data or scan parameters have moved since it was written, re-run
    `emerging treescan backtest` first.
    """
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run "
                                f"`emerging treescan backtest` first")
    df = pd.read_csv(path, parse_dates=["as_of"])
    leaf = df[df["level"] == "substance"].rename(columns={"node": "substance"})
    return leaf[["substance", "as_of", "recurrence_interval", "p_value", "llr"]]


def _frame_excess(
    counts: pd.DataFrame, ref: pd.Series, as_of: pd.Timestamp,
    study_quarters: int = STUDY_QUARTERS, max_window: int = MAX_WINDOW,
) -> tuple[pd.Index, np.ndarray, np.ndarray]:
    """Per-substance observed count and excess, for every trailing window.

    Returns `(substances, observed, excess)` where the two arrays are
    `(n_substances, max_window)` and column `w` is the trailing window of
    length `w + 1` ending at `as_of`. `excess = observed - total * p_w` is the
    same quantity `excess_share` apportions, computed once for every substance
    and window instead of once per (node, substance) pair -- the node-level
    excess is just the sum over the node's leaves, because both the observed
    count and the expectation are sums over leaves.

    Empty index when the study frame is not yet full, matching `run_scan`'s
    refusal to scan a short frame.
    """
    frame = [q for q in sorted(ref.index) if q <= as_of][-study_quarters:]
    if len(frame) < study_quarters:
        return pd.Index([], name="substance"), np.zeros((0, max_window)), \
               np.zeros((0, max_window))
    p = ref.loc[frame].to_numpy(dtype=float)
    p = p / p.sum()

    sub = counts[counts["quarter"].isin(frame)]
    wide = (sub.pivot_table(index="substance", columns="quarter", values="n",
                            aggfunc="sum", fill_value=0)
              .reindex(columns=frame, fill_value=0))
    x = wide.to_numpy(dtype=float)
    obs = _trailing_counts(x, max_window)
    total = x.sum(axis=1)
    pw = _window_shares(p, max_window)
    return wide.index, obs, obs - total[:, None] * pw[None, :]


def branch_sweep(
    path: Path = BACKTEST_PATH,
    mentions_path: Path = MENTIONS_PATH,
    cutoff: str | None = CUTOFF,
    reference: str = "deaths",
    min_excess_share: float | None = None,
    include_leaves: bool = False,
    study_quarters: int = STUDY_QUARTERS,
    max_window: int = MAX_WINDOW,
    prepared: tuple | None = None,
) -> pd.DataFrame:
    """The scan's *internal* nodes, read as a per-substance detector.

    `leaf_sweep` reads only the substance-level nodes, which is the one part
    of the scan that has an EB05 counterpart -- and therefore throws away the
    entire reason a tree scan exists. The nitazenes are 1 and 2 mentions
    apiece and can never clear a leaf test; the family node holding them is a
    single hypothesis with their pooled count. This function propagates every
    branch's signal back down to the substances inside it, so that mechanism
    can be scored on the same (substance, as_of) grid every other detector in
    `emerging.validation.benchmark` is scored on.

    A branch score is credited to a substance only when both hold:

    - **presence** -- the substance has at least one mention inside the
      node's own best window. Without this, `root` firing hands its score to
      all 211 substances, including ones that will not exist for six years.
      This is the vectorized form of `backtest`'s `as_of >= first_seen`
      guard, tightened from "has ever appeared" to "is in the window that
      fired".
    - **attribution**, when `min_excess_share` is set -- the substance
      contributes at least that fraction of the node's excess, exactly
      `excess_share`'s quantity. `Fentanyl and Fentanyl-related` fired in
      2017Q2 on fentanyl's own explosion, in the same quarter carfentanil
      appeared on one case; unguarded, that reads as a detection of
      carfentanil. `backtest` already reports this number and reads 0.5 as
      "attributable" and 0.15 as "partly attributable"; passing `None` here
      keeps the unguarded behaviour so the guard's cost can be measured
      rather than assumed.

    With `include_leaves`, a substance's own leaf row competes with its
    branches and the score is the maximum -- the whole tree read as one
    detector, which is what TreeScan actually is. A leaf is exempt from both
    guards: it is trivially present in and wholly responsible for its own
    excess.

    Ties on the score go to the *smaller* node, so `source_node` names the
    most specific branch that justifies it -- the same rule `backtest` uses
    when it picks which branch to report.

    Reads the same cached `treescan backtest` table `leaf_sweep` does, and
    inherits its assumption that the cache was written with these scan
    parameters and this `--reference`; the attribution weights are recomputed
    here from `mentions_path`, so a stale cache mixes two vintages silently.
    Re-run `emerging treescan backtest` if either has moved.
    """
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run "
                                f"`emerging treescan backtest` first")
    bt = pd.read_csv(path, parse_dates=["as_of"])
    counts, _, ref, nodes = prepared or prepare(mentions_path, cutoff, reference)
    leaves_of = {nd.name: nd.leaves for nd in nodes}
    size_of = {nd.name: nd.n_leaves for nd in nodes}

    rows = []
    for as_of, g in bt.groupby("as_of", sort=True):
        subs, obs, excess = _frame_excess(counts, ref, as_of, study_quarters,
                                          max_window)
        if not len(subs):
            continue
        pos = {s: i for i, s in enumerate(subs)}
        for r in g.itertuples(index=False):
            if r.level == "substance":
                if include_leaves:
                    rows.append((r.node, as_of, r.recurrence_interval,
                                 r.p_value, r.llr, r.node, r.level, 1, 1.0))
                continue
            w = r.window_q
            if pd.isna(w) or r.node not in leaves_of:
                continue
            idx = [pos[s] for s in leaves_of[r.node] if s in pos]
            if not idx:
                continue
            w = int(w) - 1
            e = excess[idx, w]
            node_excess = e.sum()
            # No excess to apportion -- the node is at or below expectation,
            # so no member can be said to have driven it. `_llr` already
            # zeroes those, but a cached row from a different scan setting
            # would otherwise divide by zero here.
            if node_excess <= 0:
                continue
            share = e / node_excess
            present = obs[idx, w] > 0
            keep = present if min_excess_share is None else (
                present & (share >= min_excess_share))
            for j in np.flatnonzero(keep):
                rows.append((subs[idx[j]], as_of, r.recurrence_interval,
                             r.p_value, r.llr, r.node, r.level,
                             size_of.get(r.node, len(idx)), float(share[j])))

    out = pd.DataFrame(rows, columns=[
        "substance", "as_of", "recurrence_interval", "p_value", "llr",
        "source_node", "source_level", "source_n_leaves", "excess_share"])
    if out.empty:
        return out
    out = (out.sort_values(["substance", "as_of", "recurrence_interval",
                            "source_n_leaves", "source_node"],
                           ascending=[True, True, False, True, True],
                           kind="stable")
              .drop_duplicates(subset=["substance", "as_of"], keep="first")
              .reset_index(drop=True))
    return out


def _q(t: pd.Timestamp) -> str:
    return f"{t.year}Q{t.quarter}"


def _qord(s: pd.Series) -> pd.Series:
    """Quarter index, so lags can be differenced. NaT propagates as NA."""
    return (s.dt.year * 4 + s.dt.quarter).astype("Int64")


def excess_share(
    counts: pd.DataFrame, ref: pd.Series, leaves: set[str], substance: str,
    as_of: pd.Timestamp, window_q: int, study_quarters: int = STUDY_QUARTERS,
) -> float:
    """Fraction of a node's excess in its window contributed by one substance.

    A branch signalling while a substance happens to be present is not the same
    as a branch signalling *because of* it. `Fentanyl and Fentanyl-related`
    fired in 2017Q2, the quarter carfentanil first appeared in LA — on 1 case,
    against fentanyl's own explosion in the same node. Scored here that reads
    as 0.00, which is the honest answer.

    Returns NaN when the node has no excess to apportion.
    """
    qs = [q for q in sorted(ref.index) if q <= as_of][-study_quarters:]
    if len(qs) < study_quarters or not window_q or pd.isna(window_q):
        return np.nan
    p = ref.loc[qs].to_numpy(dtype=float)
    pw = p[-int(window_q):].sum() / p.sum()

    sub = counts[counts["quarter"].isin(qs)]
    win = set(qs[-int(window_q):])

    def _oe(members: set[str]) -> tuple[float, float]:
        g = sub[sub["substance"].isin(members)]
        total = float(g["n"].sum())
        obs = float(g.loc[g["quarter"].isin(win), "n"].sum())
        return obs, total * pw

    node_obs, node_exp = _oe(leaves)
    leaf_obs, leaf_exp = _oe({substance})
    denom = node_obs - node_exp
    return (leaf_obs - leaf_exp) / denom if denom > 0 else np.nan


@app.command("scan")
def scan(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    reference: str = typer.Option("deaths", "--reference",
                                  help="deaths | mentions"),
    study_quarters: int = typer.Option(STUDY_QUARTERS, "--study-quarters"),
    max_window: int = typer.Option(MAX_WINDOW, "--max-window"),
    min_cases: int = typer.Option(MIN_CASES, "--min-cases"),
    n_sim: int = typer.Option(N_SIM, "--n-sim"),
    top: int = typer.Option(25, "--top"),
) -> None:
    """Scan as of the most recent settled quarter and write the signal table."""
    counts, denom, ref, nodes = prepare(mentions, cutoff, reference)
    as_of = max(ref.index)
    res = run_scan(counts, ref, nodes, as_of, study_quarters, max_window,
                   min_cases, n_sim)
    res["reference"] = reference

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / SCAN_PATH.name
    if reference != "deaths":
        path = out_dir / f"signals_{reference}.csv"
    res.to_csv(path, index=False)

    quarters = [q for q in sorted(ref.index) if q <= as_of][-study_quarters:]
    typer.echo(
        f"as of {_q(as_of)}; study period {_q(quarters[0])}–{_q(quarters[-1])} "
        f"({study_quarters}q), windows 1–{max_window}q, reference={reference}"
    )
    typer.echo(f"{len(res)} nodes tested x {max_window} windows, "
               f"{n_sim} Monte Carlo replicates; min {min_cases} cases")

    sig = res[res["recurrence_interval"] >= RI_THRESHOLD]
    typer.echo(f"\n{len(sig)} node(s) at recurrence interval >= "
               f"{RI_THRESHOLD:.0f} analyses (p <= {1 / RI_THRESHOLD:.3f})\n")

    disp = res.head(top)[["node", "level", "window_q", "observed", "expected",
                          "obs_exp", "llr", "p_value", "recurrence_interval"]]
    with pd.option_context("display.max_colwidth", 34):
        typer.echo(disp.to_string(index=False,
                                  float_format=lambda v: f"{v:.3f}"))
    typer.echo(f"\nwrote {path}")


@app.command("backtest")
def backtest(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    reference: str = typer.Option("deaths", "--reference"),
    study_quarters: int = typer.Option(STUDY_QUARTERS, "--study-quarters"),
    max_window: int = typer.Option(MAX_WINDOW, "--max-window"),
    n_sim: int = typer.Option(N_SIM, "--n-sim",
                              help="45 scans at this many replicates; the "
                                   "matched-budget cut needs the resolution, "
                                   "because 999 replicates floor RI at 1,000"),
    substances: str = typer.Option(",".join(KNOWN_EMERGENCES), "--substances"),
) -> None:
    """As-of sweep, then head-to-head against EB05 on the known emergences.

    Same ground truth and the same 12-quarter frame as `trends backtest`, so
    the only thing that differs is the statistic. For each known emergence,
    report the first quarter at which each method would have signalled and the
    delay from the substance's first appearance in the data.
    """
    counts, denom, ref, nodes = prepare(mentions, cutoff, reference)
    targets = [s.strip() for s in substances.split(",") if s.strip()]
    all_q = sorted(ref.index)

    # Non-leaf nodes containing each target, smallest first: a branch signal is
    # only interesting if it is the *most specific* branch that fired.
    size = {nd.name: nd.n_leaves for nd in nodes}
    leaf_sets = {nd.name: nd.leaves for nd in nodes}
    containing = {s: sorted((nd.name for nd in nodes
                             if s in nd.leaves and nd.level != "substance"),
                            key=lambda nm: (size[nm], nm))
                  for s in targets}

    frames = []
    for i in range(study_quarters, len(all_q) + 1):
        as_of = all_q[i - 1]
        try:
            res = run_scan(counts, ref, nodes, as_of, study_quarters,
                           max_window, MIN_CASES, n_sim)
        except ValueError:
            continue
        frames.append(res[res["llr"] > 0])
    bt = pd.concat(frames, ignore_index=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    bt.to_csv(out_dir / BACKTEST_PATH.name, index=False)
    swept = sorted(bt["as_of"].unique())

    # EB05 as-of sweep on the same frame and the same quarters. Every
    # substance is scored, not just the targets, because the alarm *budget*
    # below needs the full alarm count, not the targets' share of it.
    eb_rows = []
    for as_of in swept:
        i = all_q.index(as_of) + 1
        sub_denom = denom.loc[all_q[:i]]
        sub_counts = counts[counts["quarter"] <= as_of]
        try:
            r, _ = rank_substances(sub_counts, sub_denom, lag_quarters=0,
                                   recent_quarters=4,
                                   baseline_quarters=study_quarters - 4)
        except ValueError:
            continue
        r = r[r["n_recent"] > 0]
        eb_rows.append(pd.DataFrame({"as_of": as_of, "substance": r.index,
                                     "eb05": r["eb05"].to_numpy(),
                                     "n_recent": r["n_recent"].to_numpy()}))
    eb = pd.concat(eb_rows, ignore_index=True)

    # --- matched alarm budget --------------------------------------------
    # EB05 > 1.5 is unadjusted; RI >= 100 is adjusted for every node and
    # window. Comparing them directly rewards whichever is simply more
    # trigger-happy. So also cut the tree scan at the threshold that raises
    # the *same number of alarms* over the same sweep, which fixes the false-
    # alarm budget and lets detection timing be the only thing that differs.
    n_eb_alarms = int((eb["eb05"] >= ALARM_THRESHOLD).sum())
    ranked = bt.sort_values(["p_value", "llr"], ascending=[True, False],
                            kind="stable").reset_index()
    matched_idx = set(ranked.head(n_eb_alarms)["index"])
    bt["matched_alarm"] = bt.index.isin(matched_idx)
    matched_ri = (float(ranked.iloc[min(n_eb_alarms, len(ranked)) - 1]
                        ["recurrence_interval"]) if len(ranked) else np.nan)

    seen = counts[counts["n"] > 0]
    first_seen = seen.groupby("substance")["quarter"].min()

    def _first(df: pd.DataFrame) -> pd.Timestamp:
        return df["as_of"].min() if not df.empty else pd.NaT

    rows = []
    for s in targets:
        f0 = first_seen.get(s)
        is_leaf = bt["node"] == s
        leaf = bt[is_leaf & (bt["recurrence_interval"] >= RI_THRESHOLD)]
        leaf_m = bt[is_leaf & bt["matched_alarm"]]
        ebs = eb[(eb["substance"] == s) & (eb["eb05"] >= ALARM_THRESHOLD)]
        # A branch cannot have detected a substance that did not exist yet.
        # Without this, "Fentanyl and Fentanyl-related" signalling in 2016Q3
        # scores as an 18-quarter head start on para-fluorofentanyl, which
        # first appears in 2021Q1 — it was fentanyl itself all along.
        branch = bt[bt["node"].isin(containing[s])
                    & (bt["recurrence_interval"] >= RI_THRESHOLD)
                    & (bt["as_of"] >= f0)]
        best_branch, attribution = None, np.nan
        if not branch.empty:
            first = branch[branch["as_of"] == branch["as_of"].min()]
            best_branch = min(first["node"], key=lambda nm: (size[nm], nm))
            row = first[first["node"] == best_branch].iloc[0]
            attribution = excess_share(
                counts, ref, set(leaf_sets[best_branch]), s,
                row["as_of"], row["window_q"], study_quarters)
        rows.append({
            "substance": s,
            "first_seen": f0,
            "eb05_signal": _first(ebs),
            "tree_leaf_signal": _first(leaf),
            "tree_leaf_matched_signal": _first(leaf_m),
            "tree_branch_signal": _first(branch),
            "branch_node": best_branch,
            "branch_excess_share": attribution,
            "peak_llr": bt.loc[is_leaf, "llr"].max() if is_leaf.any() else np.nan,
            "best_ri": (bt.loc[is_leaf, "recurrence_interval"].max()
                        if is_leaf.any() else np.nan),
        })
    h2h = pd.DataFrame(rows)
    for c in ["eb05_signal", "tree_leaf_signal", "tree_leaf_matched_signal",
              "tree_branch_signal"]:
        h2h[c.replace("_signal", "_lag_q")] = _qord(h2h[c]) - _qord(
            h2h["first_seen"])
    h2h.to_csv(out_dir / HEADTOHEAD_PATH.name, index=False)

    typer.echo(f"as-of sweep: {len(swept)} quarters "
               f"({_q(swept[0])}–{_q(swept[-1])}), reference={reference}, "
               f"{n_sim} replicates per quarter\n")
    typer.echo(f"alarm budget over the sweep: EB05 > {ALARM_THRESHOLD} raises "
               f"{n_eb_alarms} substance-quarter alarms "
               f"({eb.loc[eb['eb05'] >= ALARM_THRESHOLD, 'substance'].nunique()}"
               f" distinct substances).")
    typer.echo(f"  tree scan at RI >= {RI_THRESHOLD:.0f}: "
               f"{int((bt['recurrence_interval'] >= RI_THRESHOLD).sum())} "
               f"node-quarter alarms "
               f"({bt.loc[bt['recurrence_interval'] >= RI_THRESHOLD, 'node'].nunique()}"
               f" distinct nodes)")
    typer.echo(f"  tree scan cut to the same {n_eb_alarms} alarms: "
               f"threshold RI {matched_ri:,.0f} "
               f"({bt.loc[bt['matched_alarm'], 'node'].nunique()} "
               f"distinct nodes)\n")

    typer.echo("head-to-head on the known emergences "
               "(quarter of first signal; lag in quarters from first appearance)")
    for _, r in h2h.iterrows():
        fmt = lambda t: "never" if pd.isna(t) else _q(t)  # noqa: E731
        lagf = lambda v: "" if pd.isna(v) else f" +{int(v)}q"  # noqa: E731
        typer.echo(
            f"  {r['substance']:22s} seen {fmt(r['first_seen'])}"
            f" | EB05 {fmt(r['eb05_signal'])}{lagf(r['eb05_lag_q'])}"
            f" | leaf {fmt(r['tree_leaf_signal'])}{lagf(r['tree_leaf_lag_q'])}"
            f" | leaf@budget {fmt(r['tree_leaf_matched_signal'])}"
            f"{lagf(r['tree_leaf_matched_lag_q'])}"
            f" | branch {fmt(r['tree_branch_signal'])}"
            f"{lagf(r['tree_branch_lag_q'])}"
            + (f" via {r['branch_node']}" if r['branch_node'] else "")
        )
        if r["branch_node"]:
            sh = r["branch_excess_share"]
            verdict = ("attributable" if sh >= 0.5 else
                       "partly attributable" if sh >= 0.15 else
                       "NOT attributable — the branch is signalling on "
                       "something else")
            typer.echo(f"  {'':22s}   branch excess from {r['substance']}: "
                       f"{sh:.0%} — {verdict}")
        if pd.isna(r["tree_leaf_signal"]):
            typer.echo(f"  {'':22s}   (leaf peaks at RI {r['best_ri']:,.0f}, "
                       f"below the {RI_THRESHOLD:.0f} line)")

    sig = bt[bt["recurrence_interval"] >= RI_THRESHOLD]
    typer.echo(f"\n{sig['node'].nunique()} distinct node(s) ever signalled "
               f"across the sweep; {len(sig)} node-quarters")
    top = (sig.groupby(["node", "level"])["llr"].max()
              .sort_values(ascending=False).head(15))
    for (node, level), v in top.items():
        n_q = int((sig["node"] == node).sum())
        typer.echo(f"  {node:44s} [{level:10s}] peak LLR {v:6.2f}  "
                   f"{n_q:2d} quarter(s)")
    typer.echo(f"\nwrote {out_dir / BACKTEST_PATH.name}")
    typer.echo(f"wrote {out_dir / HEADTOHEAD_PATH.name}")


@app.command("plot")
def plot(
    results_dir: Path = typer.Option(RESULTS_DIR, "--results-dir"),
    fig_dir: Path = typer.Option(RESULTS_DIR, "--fig-dir"),
    top: int = typer.Option(14, "--top"),
) -> None:
    """Two panels: head-to-head detection timing, and the current signal table."""
    scan_res = pd.read_csv(results_dir / SCAN_PATH.name,
                           parse_dates=["as_of"])
    h2h = pd.read_csv(results_dir / HEADTOHEAD_PATH.name,
                      parse_dates=["first_seen", "eb05_signal",
                                   "tree_leaf_signal",
                                   "tree_leaf_matched_signal",
                                   "tree_branch_signal"])

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(15.5, 7.0), gridspec_kw={"width_ratios": [1.0, 1.12]})

    # -- Panel A: when each method would have signalled -------------------
    h2h = h2h.sort_values("first_seen").reset_index(drop=True)
    y = np.arange(len(h2h))
    series = [("first_seen", "first appearance", MUTED, "o", 46),
              ("eb05_signal", "EB05 > 1.5", INK_2, "s", 44),
              ("tree_leaf_signal", "leaf, RI ≥ 100", BLUE, "D", 42),
              ("tree_leaf_matched_signal", "leaf @ matched budget", YELLOW,
               "v", 52),
              ("tree_branch_signal", "branch", ORANGE, "^", 58)]
    for i in y:
        pts = [h2h.loc[i, c] for c, *_ in series if pd.notna(h2h.loc[i, c])]
        if len(pts) > 1:
            ax1.plot([min(pts), max(pts)], [i, i], color="#e4e2dd", lw=2.4,
                     zorder=1, solid_capstyle="round")
    # Methods that signal in the same quarter land on the same point, so the
    # last one drawn hides the rest — three of five rows have at least one
    # exact tie. Offset each series by a fixed fraction of a row instead.
    offsets = np.linspace(-0.19, 0.19, len(series))
    for (col, label, color, marker, size), dy in zip(series, offsets):
        ok = h2h[col].notna()
        ax1.scatter(h2h.loc[ok, col], y[ok] + dy, s=size, marker=marker,
                    color=color, label=label, zorder=3,
                    edgecolor="white", linewidth=1.1)
    edge = pd.Timestamp("2026-07-01")
    for i in y:
        bits = []
        if pd.isna(h2h.loc[i, "tree_leaf_signal"]):
            bits.append("leaf never signals")
        sh = h2h.loc[i, "branch_excess_share"]
        if pd.notna(sh):
            pct = "≈0%" if abs(sh) < 0.005 else f"{sh:.0%}"
            bits.append(f"branch {pct} attributable")
        if bits:
            ax1.text(edge, i, "  ·  ".join(bits), va="center", ha="left",
                     fontsize=7.6, color=MUTED)

    ax1.set_yticks(y)
    ax1.set_yticklabels(h2h["substance"], fontsize=9.5)
    ax1.set_ylim(-0.7, len(h2h) - 0.3)
    ax1.set_xlim(pd.Timestamp("2015-06-01"), pd.Timestamp("2030-06-01"))
    ax1.set_xticks(pd.to_datetime([f"{yr}-01-01" for yr in range(2016, 2027, 2)]))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax1.invert_yaxis()
    ax1.set_title("Would it have caught the known emergences?",
                  fontsize=12.5, color=INK, loc="left", pad=10)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.42, -0.07), ncol=5,
               frameon=False, fontsize=8.2, columnspacing=1.2,
               handletextpad=0.35)
    ax1.grid(axis="x", color="#eceae5", lw=0.9)
    ax1.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax1.spines[s].set_visible(False)
    ax1.spines["bottom"].set_color(MUTED)
    ax1.tick_params(colors=INK_2, labelsize=9)

    # -- Panel B: current scan ---------------------------------------------
    cur = scan_res[scan_res["llr"] > 0].head(top).iloc[::-1]
    yb = np.arange(len(cur))
    colors = [LEVEL_COLOR.get(lv, BLUE) for lv in cur["level"]]
    ax2.barh(yb, cur["obs_exp"], height=0.62, color=colors, zorder=3)
    ax2.axvline(1.0, color=INK_2, lw=1.2, zorder=4)
    for i, (_, r) in enumerate(cur.iterrows()):
        ri = r["recurrence_interval"]
        star = " *" if ri >= RI_THRESHOLD else ""
        ax2.text(r["obs_exp"] + 0.03, i,
                 f"{int(r['observed'])} vs {r['expected']:.0f} exp · "
                 f"{int(r['window_q'])}q · RI {ri:,.0f}{star}",
                 va="center", ha="left", fontsize=8, color=INK_2)
    ax2.set_yticks(yb)
    ax2.set_yticklabels(
        [f"{n[:38]}  [{lv[:3]}]" for n, lv in zip(cur["node"], cur["level"])],
        fontsize=9)
    ax2.set_xlim(0, max(2.0, float(cur["obs_exp"].max()) * 1.62))
    ax2.set_ylim(-0.7, len(cur) - 0.3)
    ax2.set_xlabel("observed / expected in the scanning window", fontsize=9.5,
                   color=INK_2)
    as_of = pd.Timestamp(scan_res["as_of"].iloc[0])
    ax2.set_title(f"Strongest nodes as of {_q(as_of)}",
                  fontsize=12.5, color=INK, loc="left", pad=10)
    ax2.text(1.0, 1.028, f"*  recurrence interval ≥ {RI_THRESHOLD:.0f}",
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=8.5,
             color=INK_2)
    ax2.grid(axis="x", color="#eceae5", lw=0.9)
    ax2.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)
    ax2.spines["bottom"].set_color(MUTED)
    ax2.tick_params(colors=INK_2, labelsize=9)

    handles = [plt.Line2D([], [], marker="s", ls="", color=LEVEL_COLOR[lv],
                          label=lv, markersize=7)
               for lv in ["substance", "family", "category", "superclass"]]
    ax2.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.10),
               ncol=4, frameon=False, fontsize=8.5, columnspacing=1.4,
               handletextpad=0.35)

    fig.text(0.008, 0.015,
             f"Tree-temporal scan over {len(scan_res)} overlapping nodes x "
             f"{MAX_WINDOW} trailing windows in a {STUDY_QUARTERS}-quarter "
             "frame; p-values come from the Monte Carlo distribution of the "
             "maximum LLR, so they are already adjusted for the whole search. "
             "Left: EB05's 1.5 line is unadjusted, so 'leaf @ budget' re-cuts "
             "the scan to raise the same number of alarms over the sweep. "
             "Recurrence interval is per single analysis and is not corrected "
             "for repeated quarterly looks.",
             fontsize=7.6, color=MUTED, va="bottom", wrap=True)

    fig.tight_layout(rect=(0, 0.115, 1, 1))
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = fig_dir / "scan.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
