"""
Validate `treescan.py` against Martin Kulldorff's TreeScan binary.

`treescan.py` is a from-scratch implementation of the conditional tree-temporal
scan. That is a liability until it is checked against the reference: a scan
statistic that is subtly wrong still returns plausible tables, and every
recurrence interval downstream is then wrong by an unknown factor. The unit
tests in `tests/test_treescan.py` show the implementation is *self*-consistent
and correctly calibrated; only this comparison shows it computes the same
statistic the literature does.

**What is compared, and what cannot be.** TreeScan's tree-temporal model
distributes each node's cases uniformly over the study period (user guide,
"Tree-Temporal Scan Statistic": the analysis is conditioned on the cases in
each node, and "the time of each case ... under the null hypothesis is assumed
to be uniform across the study time period"). Our implementation generalises
that one step, replacing the uniform with an arbitrary reference series — OD
deaths or total mentions per quarter — because a quarter with 30% fewer
overdose deaths is not a quarter with 30% less time.

So the comparison runs our scan with `reference="uniform"`, which reduces it
to exactly TreeScan's model. That validates the tree/DAG construction, the node
aggregation, the window enumeration, the LLR, and the Monte Carlo machinery —
everything except the generalisation, which enters through the single line
`_window_shares` and is covered by `test_window_shares_are_cumulative_from_the_end`.

Matched settings, all of which had to be turned off or overridden because
TreeScan's defaults differ from ours:

    scan-type=1 (tree-temporal)        conditional-type=2 (NODE, not NODEANDTIME)
    probability-model=2 (uniform)      prospective-analysis=y (trailing windows)
    maximum-window-type=1 (fixed)      maximum-window-fixed=MAX_WINDOW
    minimum-window-fixed=1             apply-risk-window-restriction=n
    data-only-on-leaves=y              include-identical-parent-cuts=y

`apply-risk-window-restriction` defaults to y at 20%, which silently drops the
short windows; `include-identical-parent-cuts` defaults to n, which hides
exactly the duplicate-leaf-set nodes our `_dedupe` removes, so leaving it off
would make the two node lists disagree for a reason unrelated to the statistic.

Usage:

    emerging treescan_validate compare
    emerging treescan_validate compare --tree-mode strict
    emerging treescan_validate compare --n-sim 9999
"""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

import pandas as pd
import typer

from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import ROOT, results_dir
from emerging.core.tree import Node, _dedupe
from emerging.analysis.treescan import (MAX_WINDOW, STUDY_QUARTERS, prepare, run_scan)

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("treescan")

TREESCAN_BIN = ROOT / "vendor" / "treescan" / "treescan64"
WORK_DIR = ROOT / "vendor" / "treescan" / "work"
REPORT_PATH = RESULTS_DIR / "vs_reference.csv"

# Tolerances. The LLR is a deterministic function of integer counts and a
# rational expectation, so anything above printf noise is a real disagreement;
# TreeScan writes 6 decimal places.
LLR_TOL = 1e-5
EXPECTED_TOL = 5e-3      # TreeScan prints expected counts to 2 decimals
# TreeScan forms no cut from fewer than 2 cases, in the node or in the window.
# Our `min_cases` is set to the same so both search the same space.
MIN_CUT_CASES = 2


def _leaf_closures(edges: list[tuple[str, str]], leaves: set[str]) -> dict:
    """Descendant-leaf set of every node, from the emitted edge list.

    The point of computing this from the *file* rather than from our objects
    is that it checks the export, which is where a DAG is most likely to be
    mistranslated.
    """
    children: dict[str, list[str]] = {}
    for child, parent in edges:
        if parent:
            children.setdefault(parent, []).append(child)

    memo: dict[str, frozenset[str]] = {}

    def closure(node: str) -> frozenset[str]:
        if node in memo:
            return memo[node]
        memo[node] = frozenset()          # cycle guard
        own = frozenset({node}) if node in leaves else frozenset()
        for c in children.get(node, []):
            own |= closure(c)
        memo[node] = own
        return own

    for n in set(list(children) + [c for c, _ in edges]):
        closure(n)
    return memo


def select_nodes(nodes: list[Node], tree_mode: str, live: set[str]) -> list[Node]:
    """Filter to the requested structure and to leaves present in the window."""
    if tree_mode == "strict":
        nodes = [nd for nd in nodes if nd.level != "family"]
    elif tree_mode != "dag":
        raise ValueError("tree-mode must be 'strict' or 'dag'")
    trimmed = [Node(nd.name, nd.level, nd.leaves & live, nd.parents)
               for nd in nodes]
    return _dedupe([nd for nd in trimmed if nd.leaves])


def build_edges(nodes: list[Node]) -> list[tuple[str, str]]:
    """`(child, parent)` edges reproducing each node's leaf set as its closure.

    Parents are resolved to the nearest *surviving* ancestor, because
    `_dedupe` removes nodes whose leaf set duplicates a more specific one — a
    category holding a single substance disappears, and its leaf then has to
    attach to the superclass instead.
    """
    by_name = {nd.name: nd for nd in nodes}
    by_level: dict[str, list[Node]] = {}
    for nd in nodes:
        by_level.setdefault(nd.level, []).append(nd)

    def nearest(candidates: list[str]) -> str:
        for c in candidates:
            if c in by_name:
                return c
        return "Root"

    edges: list[tuple[str, str]] = [("Root", "")]
    for nd in nodes:
        if nd.level == "root":
            continue
        if nd.level in ("superclass", "family"):
            edges.append((nd.name, "Root"))
        elif nd.level == "category":
            sup = [s.name for s in by_level.get("superclass", [])
                   if nd.leaves <= s.leaves]
            edges.append((nd.name, nearest(sup)))
        elif nd.level == "substance":
            cat = [c.name for c in by_level.get("category", [])
                   if nd.leaves <= c.leaves]
            sup = [s.name for s in by_level.get("superclass", [])
                   if nd.leaves <= s.leaves]
            edges.append((nd.name, nearest(cat + sup)))
            # Every family containing the substance is an additional parent.
            # This is what makes the exported structure a DAG.
            for fam in by_level.get("family", []):
                if nd.leaves <= fam.leaves:
                    edges.append((nd.name, fam.name))
    # `root` in our collection is TreeScan's `Root`; drop the duplicate name.
    return [(c, p) for c, p in edges if c != "root"]


def verify_export(edges, nodes: list[Node], leaves: set[str]) -> list[str]:
    """Every node's closure in the exported file must equal its leaf set."""
    closures = _leaf_closures(edges, leaves)
    problems = []
    for nd in nodes:
        if nd.level == "root":
            continue
        got = closures.get(nd.name, frozenset())
        if got != nd.leaves:
            problems.append(
                f"{nd.name}: exported closure has {len(got)} leaves, "
                f"our node has {len(nd.leaves)} "
                f"(missing {sorted(nd.leaves - got)[:3]}, "
                f"extra {sorted(got - nd.leaves)[:3]})")
    root = closures.get("Root", frozenset())
    if root != leaves:
        problems.append(f"Root closure {len(root)} != {len(leaves)} leaves")
    return problems


def make_ids(nodes: list[Node], leaves: set[str]) -> dict[str, str]:
    """Map every node and leaf name to a comma-free synthetic ID.

    TreeScan's input files are comma-delimited and its reader does not accept
    quoting, so `1,1-Difluoroethane` and `BTMPS (Bis(2,2,6,6-tetramethyl-...))`
    are parsed as extra columns — the first attempt failed with "record 17 has
    node referencing self as parent". Rather than sanitise (which can collide),
    every name gets an opaque ID and the results are mapped back.
    """
    names = sorted({nd.name for nd in nodes} | leaves)
    ids = {name: f"n{i:04d}" for i, name in enumerate(names)}
    ids["Root"] = "Root"
    return ids


def write_inputs(work: Path, nodes: list[Node], wide: pd.DataFrame,
                 quarters: list[pd.Timestamp]) -> tuple[Path, Path, dict]:
    """Write `tree.tre` and `cases.cas`; abort if the closures disagree."""
    work.mkdir(parents=True, exist_ok=True)
    leaves = set(wide.index)
    edges = build_edges(nodes)

    problems = verify_export(edges, nodes, leaves)
    if problems:
        raise typer.BadParameter(
            "the exported tree does not reproduce our node sets — the "
            "comparison would be meaningless:\n  " + "\n  ".join(problems[:8]))

    ids = make_ids(nodes, leaves)
    tre = work / "tree.tre"
    tre.write_text("".join(f"{ids[c]},{ids[p] if p else ''}\n"
                           for c, p in edges))

    rows = []
    for t, q in enumerate(quarters, start=1):
        col = wide[q]
        for leaf, n in col[col > 0].items():
            rows.append(f"{ids[leaf]},{int(n)},{t}\n")
    cas = work / "cases.cas"
    cas.write_text("".join(rows))

    pd.DataFrame({"id": list(ids.values()), "name": list(ids)}).to_csv(
        work / "names.csv", index=False)
    return tre, cas, ids


def run_reference(work: Path, tre: Path, cas: Path, n_q: int, max_window: int,
                  n_sim: int, tree_mode: str, seed: int, processes: int,
                  tag: str = "reference") -> tuple[pd.DataFrame, float, int]:
    """Run the TreeScan binary and return (cuts, wall seconds, nodes in tree)."""
    out = work / f"{tag}.txt"
    cmd = [
        str(TREESCAN_BIN),
        "--tree-filename", str(tre), "--count-filename", str(cas),
        "--scan-type", "1",            # tree-temporal
        "--conditional-type", "2",     # condition on the node only
        "--probability-model", "2",    # uniform over time
        "--scan-rate-type", "0",       # high rate only, as we do
        "--date-precision", "1",
        "--data-time-range", f"[1,{n_q}]",
        "--data-only-on-leaves", "y",
        "--allow-multi-parent-nodes", "y" if tree_mode == "dag" else "n",
        "--allow-multiple-roots", "n",
        "--maximum-window-type", "1",
        "--maximum-window-fixed", str(max_window),
        "--minimum-window-fixed", "1",
        "--apply-risk-window-restriction", "n",
        "--prospective-analysis", "y",
        "--include-identical-parent-cuts", "y",
        # TreeScan will not form a cut from fewer than 2 cases. Discovered by
        # diffing the node lists: every node it declined had either 1 case in
        # total or 1 case in the window. This is a property of the search
        # space, not of reporting, so our `min_cases` is set to match --
        # otherwise our null maximum ranges over nodes theirs never sees.
        "--minimum-node-cases", str(MIN_CUT_CASES),
        "--monte-carlo-replications", str(n_sim),
        "--randomization-seed", str(seed),
        "--parallel-processes", str(processes),
        "--results-csv", "y",
        "--results-filename", str(out),
    ]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"TreeScan failed:\n{proc.stdout[-3000:]}\n"
                           f"{proc.stderr[-2000:]}")

    n_nodes = None
    for line in out.read_text().splitlines():
        if line.startswith("Number of Nodes"):
            n_nodes = int(line.split(":")[1])
            break

    csv = out.with_suffix(".csv")
    if not csv.exists():
        raise RuntimeError(f"TreeScan wrote no CSV at {csv}")
    df = pd.read_csv(csv, dtype={"Cut No.": str})
    # Rows like `2_1` are the contributing children of cut 2, not cuts.
    df = df[~df["Cut No."].str.contains("_", na=False)].copy()
    df = df.rename(columns={
        "Node Identifier": "node", "Node Cases": "n_study",
        "Time Window Start": "win_start", "Time Window End": "win_end",
        "Cases in Window": "observed", "Expected Cases": "expected",
        "Log Likelihood Ratio": "llr", "P-value": "p_value",
        "Recurrence Interval": "recurrence_interval"})
    for c in ["llr", "expected", "p_value"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["window_q"] = df["win_end"] - df["win_start"] + 1
    return df, elapsed, n_nodes


def _setup(mentions: Path, tree_mode: str, study_quarters: int):
    """Shared loading: counts in the study window, live leaves, node set."""
    counts, denom, ref, nodes = prepare(mentions, reference="deaths")
    as_of = max(ref.index)
    quarters = [q for q in sorted(ref.index) if q <= as_of][-study_quarters:]
    sub = counts[counts["quarter"].isin(quarters)]
    wide = (sub.pivot_table(index="substance", columns="quarter", values="n",
                            aggfunc="sum", fill_value=0)
              .reindex(columns=quarters, fill_value=0))
    sel = select_nodes(nodes, tree_mode, set(wide.index))
    return counts, as_of, quarters, wide, sel


@app.command("export")
def export(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    work: Path = typer.Option(WORK_DIR, "--work-dir"),
    tree_mode: str = typer.Option("dag", "--tree-mode"),
    study_quarters: int = typer.Option(STUDY_QUARTERS, "--study-quarters"),
) -> None:
    """Write TreeScan input files and verify they reproduce our node sets."""
    _, as_of, quarters, wide, sel = _setup(mentions, tree_mode, study_quarters)
    tre, cas, _ = write_inputs(work, sel, wide, quarters)
    typer.echo(f"{len(wide)} leaves, {len(sel)} nodes ({tree_mode})")
    typer.echo(f"wrote {tre}\nwrote {cas}")
    typer.echo("every node's closure in the exported tree matches our leaf set")


@app.command("compare")
def compare(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    work: Path = typer.Option(WORK_DIR, "--work-dir"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    tree_mode: str = typer.Option("dag", "--tree-mode",
                                  help="dag (our full node set) | strict "
                                       "(no family layer, a real tree)"),
    study_quarters: int = typer.Option(STUDY_QUARTERS, "--study-quarters"),
    max_window: int = typer.Option(MAX_WINDOW, "--max-window"),
    n_sim: int = typer.Option(9999, "--n-sim"),
    seed: int = typer.Option(12345678, "--seed"),
    processes: int = typer.Option(1, "--processes",
                                  help="TreeScan worker processes; 1 makes "
                                       "the timing comparison single-threaded "
                                       "like ours"),
) -> None:
    """Run both implementations on the same data and report every difference."""
    if not TREESCAN_BIN.exists():
        raise typer.BadParameter(
            f"TreeScan binary not found at {TREESCAN_BIN}. Download it from "
            f"treescan.org (registration required) and extract it there.")

    counts, as_of, quarters, wide, sel = _setup(mentions, tree_mode,
                                                study_quarters)
    typer.echo(f"as of {as_of.year}Q{as_of.quarter}; {study_quarters} quarters, "
               f"{len(wide)} leaves, {len(sel)} nodes ({tree_mode}), "
               f"windows 1–{max_window}q, {n_sim} replicates")

    tre, cas, ids = write_inputs(work, sel, wide, quarters)
    typer.echo(f"exported {tre.name} ({sum(1 for _ in tre.open())} edges) and "
               f"{cas.name}; every node's closure matches our leaf set")

    back = {v: k for k, v in ids.items()}
    # Two runs. With 0 replications TreeScan lists every cut it can form
    # rather than only those that beat a replicate, which widens the exact
    # LLR comparison from a handful of nodes to all of them; the second run
    # is the one that produces p-values.
    det_df, det_secs, n_nodes = run_reference(
        work, tre, cas, len(quarters), max_window, 0, tree_mode, seed,
        processes, tag="reference_det")
    det_df["node"] = det_df["node"].map(back)
    ref_df, ref_secs, _ = run_reference(
        work, tre, cas, len(quarters), max_window, n_sim, tree_mode, seed,
        processes)
    ref_df["node"] = ref_df["node"].map(back)
    if n_nodes is not None and n_nodes != len(sel):
        typer.echo(f"  WARNING: TreeScan built {n_nodes} nodes, we have "
                   f"{len(sel)} — the search spaces differ")
    else:
        typer.echo(f"TreeScan built {n_nodes} nodes; same search space")

    # Ours, with a uniform reference series so the two models coincide.
    uniform = pd.Series(1.0, index=pd.Index(quarters, name="quarter"))
    t0 = time.perf_counter()
    ours = run_scan(counts, uniform, sel, as_of, study_quarters, max_window,
                    min_cases=MIN_CUT_CASES, n_sim=n_sim, seed=seed)
    our_secs = time.perf_counter() - t0

    merged = ours.merge(ref_df, on="node", how="outer",
                        suffixes=("_ours", "_ref"), indicator=True)
    both = merged[merged["_merge"] == "both"].copy()
    positive = both[(both["llr_ours"] > 0) | (both["llr_ref"] > 0)].copy()

    positive["llr_diff"] = (positive["llr_ours"] - positive["llr_ref"]).abs()
    positive["exp_diff"] = (positive["expected_ours"]
                            - positive["expected_ref"]).abs()
    positive["window_same"] = (positive["window_q_ours"]
                               == positive["window_q_ref"])
    positive["p_diff"] = (positive["p_value_ours"]
                          - positive["p_value_ref"]).abs()

    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["node", "level", "n_study_ours", "n_study_ref", "window_q_ours",
            "window_q_ref", "observed_ours", "observed_ref", "expected_ours",
            "expected_ref", "llr_ours", "llr_ref", "llr_diff",
            "p_value_ours", "p_value_ref", "recurrence_interval_ours",
            "recurrence_interval_ref"]
    positive.sort_values("llr_ref", ascending=False)[cols].to_csv(
        out_dir / REPORT_PATH.name, index=False)

    # ---- report ---------------------------------------------------------
    only_ours = merged[(merged["_merge"] == "left_only")
                       & (merged["llr_ours"] > 0)]
    only_ref = merged[merged["_merge"] == "right_only"]
    # TreeScan stops listing cuts once they beat no simulated replicate at all.
    # Ours reports them anyway, so a node only we list is a disagreement only
    # if it carries evidence; at p = 1 there is nothing to disagree about.
    unexplained = only_ours[only_ours["p_value_ours"] < 1.0]

    typer.echo(f"\nnodes with a positive LLR in either: {len(positive)}")
    if len(only_ours):
        typer.echo(f"  listed only by us:         {len(only_ours)}, of which "
                   f"{len(unexplained)} with p < 1 "
                   f"(TreeScan omits cuts that beat no replicate)")
    if len(only_ref):
        typer.echo(f"  listed only by TreeScan:   {len(only_ref)} "
                   f"({', '.join(map(str, only_ref['node'].head(5)))})")

    n_bad_window = int((~positive["window_same"]).sum())
    det = ours.merge(det_df, on="node", how="inner", suffixes=("_ours", "_ref"))
    det["llr_diff"] = (det["llr_ours"] - det["llr_ref"]).abs()
    det["exp_diff"] = (det["expected_ours"] - det["expected_ref"]).abs()
    det_bad_window = int((det["window_q_ours"] != det["window_q_ref"]).sum())
    det_obs = int((det["observed_ours"] == det["observed_ref"]).sum())
    det_only_ref = set(det_df["node"]) - set(ours["node"])

    typer.echo(f"\nDETERMINISTIC PART (must match exactly) — "
               f"{len(det)} nodes, every cut TreeScan can form")
    typer.echo(f"  best window identical:  {len(det) - det_bad_window}/{len(det)}")
    typer.echo(f"  max |LLR difference|:   {det['llr_diff'].max():.3e}"
               f"   (tolerance {LLR_TOL:.0e})")
    typer.echo(f"  max |expected diff|:    {det['exp_diff'].max():.3e}"
               f"   (tolerance {EXPECTED_TOL:.0e}, TreeScan prints 2 dp)")
    typer.echo(f"  cases in window match:  {det_obs}/{len(det)}")
    if det_only_ref:
        typer.echo(f"  cuts TreeScan forms and we do not: {len(det_only_ref)} "
                   f"({', '.join(sorted(det_only_ref)[:4])})")
    obs_same = int((positive["observed_ours"] == positive["observed_ref"]).sum())

    typer.echo("\nMONTE CARLO PART (must agree within simulation error)")
    typer.echo(f"  max |p difference|:     {positive['p_diff'].max():.4f}")
    typer.echo(f"  mean |p difference|:    {positive['p_diff'].mean():.4f}")
    agree = int(((positive["p_value_ours"] <= 0.01)
                 == (positive["p_value_ref"] <= 0.01)).sum())
    typer.echo(f"  same verdict at p=.01:  {agree}/{len(positive)}")

    typer.echo("\nPERFORMANCE")
    typer.echo(f"  TreeScan (C++, {processes} proc), 0 reps:{det_secs:7.2f}s"
               f"   (fixed cost: parse, tree build, one scan)")
    typer.echo(f"  TreeScan (C++, {processes} proc), {n_sim} reps:"
               f"{ref_secs:7.2f}s")
    typer.echo(f"  ours (Python/numpy),      {n_sim} reps:{our_secs:7.2f}s")
    if ref_secs > det_secs and our_secs > 0:
        typer.echo(f"  simulation only:  TreeScan {ref_secs - det_secs:.2f}s "
                   f"vs ours {our_secs:.2f}s "
                   f"({(ref_secs - det_secs) / our_secs:.1f}x)")

    worst = positive.nlargest(6, "llr_diff")[
        ["node", "llr_ours", "llr_ref", "llr_diff", "window_q_ours",
         "window_q_ref", "p_value_ours", "p_value_ref"]]
    typer.echo("\nlargest LLR disagreements")
    typer.echo(worst.to_string(index=False,
                               float_format=lambda v: f"{v:.6f}"))

    ok = (det["llr_diff"].max() < LLR_TOL and det_bad_window == 0
          and det_obs == len(det) and not det_only_ref
          and positive["llr_diff"].max() < LLR_TOL and n_bad_window == 0
          and obs_same == len(positive)
          and len(unexplained) == 0 and len(only_ref) == 0)
    typer.echo("\n" + ("PASS — the two implementations compute the same "
                       "statistic." if ok else
                       "FAIL — see the disagreements above."))
    typer.echo(f"wrote {out_dir / REPORT_PATH.name}")
    raise typer.Exit(0 if ok else 1)


if __name__ == "__main__":
    app()
