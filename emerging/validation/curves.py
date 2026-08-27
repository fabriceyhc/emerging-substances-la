"""
Substance-level ROC and precision-recall curves for every detector in
`benchmark.MODELS`.

**Why this exists.** `benchmark.run` answers "how does this detector do at
one threshold" (twice, at a fixed line and at a matched alarm budget) --
useful for comparing configurations, but it does not show the *shape* of the
trade-off a model offers, which is the thing SAPC/LACME actually needs to see
to pick an operating point: at every possible alarm volume, what recall does
each model buy, and at every recall level, how many of the alarms are real.
This module sweeps every threshold at once instead of reading two fixed
points off it.

**Why substance-level, not substance-quarter-level.** One Bernoulli trial per
substance, the same unit `ground_truth.precision_rate` already uses, for the
same reason: pooling every alarmed *quarter* into one n would inflate the
sample with correlated observations (four quarters of one sustained mistake
are not four independent trials). A substance's score here is the *max*
`model_score` it is ever given -- ROC/PR need one scalar per instance -- but
*where* that max is allowed to come from differs by label, and getting this
wrong is what silently inflated an earlier version of this module's AUCs well
past anything `benchmark.md`'s per-quarter numbers would predict:

- A **positive**'s max is taken only from as-of quarters inside its own
  known interval (`known_emergences.csv`'s `interval_start`/`interval_end`),
  matching `score_model`'s "detected" test. Without this restriction, any
  quarter anywhere in the ~40-quarter sweep counts -- so a substance that
  ever looked unusual even once, years before or after its real emergence,
  reads as a perfect catch. A ranking detector will eventually say that
  about almost any real positive over a long enough history, which is why
  the unrestricted version put recall at 1.0 for `eb05` at its own deployed
  threshold when `benchmark.md`'s `mean_recall` (temporal coverage of the
  true window) is only 0.43 there -- two different questions wearing the
  same name.
- A **negative** (labelled `no_emergence` or unreviewed) has no legitimate
  window to restrict to, so its max is taken across the whole sweep -- any
  alarm on it, whenever it happens, is a genuine false alarm.

A substance a model never scored at all (e.g. `treescan` covers only 95 of
211) is floored the same way `benchmark.ENSEMBLE_FLOOR` already treats an
unscored substance-quarter -- zero evidence, not a missing observation.

**Where the negatives come from.** `known_emergences.csv` hand-labels 22
substances: 15 real local emergences (positive) and 7 documented
`no_emergence` cases. The other ~189 of the panel's 211 substances were never
individually reviewed -- treated here as true negatives too, on the
reasoning (confirmed against this project's own alarm history) that a
substance no detector in this benchmark ever flagged is exactly what an
unreviewed true negative looks like. That is a weaker guarantee than the 7
hand-checked cases, but with 211 substances total there is no realistic path
to hand-labelling the rest, and 196 negatives is enough sample to trace a
real curve rather than the 7-point one the hand-labelled set alone would give.

Usage:
    emerging curves run
    emerging curves run --models eb05,eb05-dual,eb05-v2-role,treescan,nb-trend
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

from emerging.analysis.trends import load_quarterly
from emerging.config import BASELINE_QUARTERS, CUTOFF, RECENT_QUARTERS
from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import KNOWN_EMERGENCES_PATH, results_dir
from emerging.validation.benchmark import (ENSEMBLE_FLOOR, MODELS,
                                           MODELS_BY_ID, Model, build_sweeps,
                                           model_score)
from emerging.validation.ground_truth import (_interval_quarters, _q_ordinal,
                                              load_ground_truth)
from emerging.viz import AQUA, BLUE, INK, INK_2, MUTED, ORANGE, YELLOW

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("curves")

# Six detectors chosen to span the project's real candidates rather than
# dumping all ~40 `benchmark.MODELS` onto one illegible plot: the deployed
# default, its two solo alternative families (docs/findings/benchmark.md),
# and the escalation from no shrinkage that motivates EB05 at all.
DEFAULT_MODELS = ("ratio", "eb05", "eb05-dual", "eb05-v2-role", "treescan",
                  "nb-trend")

# Plotting an ensemble's own score is meaningless as-is: it never sets
# `ensemble` in `ENSEMBLE_FLOOR` because a substance the *combination* never
# scored is already resolved inside `build_ensemble_sweep` (each mode fills
# missing components at their own floor before combining). 0.0 is the right
# "no evidence" value for every combine mode here -- `"threshold"`/`"veto"`
# both read as "no component thinks this is elevated at all" at 0.
_ENSEMBLE_CURVE_FLOOR = 0.0

_PALETTE = [BLUE, ORANGE, AQUA, YELLOW, INK_2, "#8a5fd6"]


def _floor(model: Model) -> float:
    if model.statistic == "ensemble":
        return _ENSEMBLE_CURVE_FLOOR
    return ENSEMBLE_FLOOR[model.statistic]


def substance_universe(mentions: Path, cutoff: str | None) -> list[str]:
    """Every substance the extract ever names, cutoff-filtered -- the full
    panel a curve's negatives are drawn from."""
    counts, _ = load_quarterly(mentions, cutoff=cutoff, quiet=True)
    return sorted(counts["substance"].unique())


def substance_labels(gt: pd.DataFrame, universe: list[str]) -> pd.Series:
    """1 for a labelled real emergence, 0 for everything else in the panel
    (both the 7 hand-labelled `no_emergence` cases and the unreviewed rest --
    see the module docstring for why both count as negatives here)."""
    if "no_emergence" not in gt.columns:
        gt = gt.assign(no_emergence=False)
    positives = set(gt.loc[~gt["no_emergence"], "substance"])
    return pd.Series({s: int(s in positives) for s in universe}, name="label")


def _positive_windows(gt: pd.DataFrame) -> dict[str, set[int]]:
    """Each labelled positive's true interval(s), as a set of quarter
    ordinals -- `no_emergence` cases excluded, since they have no real window
    to restrict to (see `substance_scores`)."""
    no_em = gt["no_emergence"] if "no_emergence" in gt.columns else pd.Series(
        False, index=gt.index)
    pos = gt[~no_em]
    return {
        sub: set().union(*(_interval_quarters(a, b)
                           for a, b in zip(g["interval_start"], g["interval_end"])))
        for sub, g in pos.groupby("substance")
    }


def substance_scores(sweep: pd.DataFrame, model: Model, universe: list[str],
                     gt: pd.DataFrame) -> pd.Series:
    """One score per substance: the max this model ever gave it, floored for
    a substance the sweep never scored at all.

    A labelled positive is scored only inside its own known interval; a
    negative (no such interval) is scored across the whole sweep. See the
    module docstring for why this restriction matters.
    """
    s = model_score(sweep, model)
    sw = pd.DataFrame({"substance": sweep["substance"].to_numpy(),
                       "q_ord": sweep["as_of"].map(_q_ordinal).to_numpy(),
                       "score": s.to_numpy()})

    keep = pd.Series(True, index=sw.index)
    for sub, window in _positive_windows(gt).items():
        sel = sw["substance"] == sub
        keep.loc[sel] = sw.loc[sel, "q_ord"].isin(window)
    sw = sw[keep]

    per_sub = sw.groupby("substance")["score"].max()
    return per_sub.reindex(universe).fillna(_floor(model))


def _rates_at_thresholds(
    scores: pd.Series, labels: pd.Series,
) -> tuple[pd.DataFrame, int, int]:
    """One row per distinct score value in `scores`, descending, with the
    cumulative confusion-matrix counts and rates at that cut -- ties in score
    are collapsed to a single operating point rather than plotted as several,
    which is what would happen sweeping row by row instead of by distinct
    value (common here: `_floor` alone ties every unscored substance).

    Prepends the threshold-at-infinity origin (nothing alarmed: recall 0,
    fpr 0, precision undefined) so every curve starts at the same point.
    """
    df = pd.DataFrame({"score": scores.to_numpy(), "label": labels.to_numpy()})
    n_pos = int(df["label"].sum())
    n_neg = int(len(df) - n_pos)

    grouped = df.groupby("score")["label"].agg(pos="sum", n="count")
    grouped = grouped.sort_index(ascending=False)
    tp = grouped["pos"].cumsum()
    total = grouped["n"].cumsum()
    fp = total - tp

    recall = tp / n_pos if n_pos else pd.Series(np.nan, index=grouped.index)
    fpr = fp / n_neg if n_neg else pd.Series(np.nan, index=grouped.index)
    precision = tp / total

    out = pd.DataFrame({
        "threshold": grouped.index.to_numpy(), "tp": tp.to_numpy(),
        "fp": fp.to_numpy(), "recall": recall.to_numpy(),
        "fpr": fpr.to_numpy(), "precision": precision.to_numpy(),
    })
    origin = pd.DataFrame([{"threshold": np.inf, "tp": 0, "fp": 0,
                            "recall": 0.0, "fpr": 0.0, "precision": np.nan}])
    return pd.concat([origin, out], ignore_index=True), n_pos, n_neg


def _roc_auc(curve: pd.DataFrame) -> float:
    x, y = curve["fpr"].to_numpy(), curve["recall"].to_numpy()
    order = np.argsort(x)
    return float(np.trapezoid(y[order], x[order]))


def _average_precision(curve: pd.DataFrame) -> float:
    """AP = sum over thresholds of (recall_k - recall_{k-1}) * precision_k --
    the step-function estimate (precision held at its actual achieved value
    over each recall interval), not a linear interpolation between precision
    points. PR curves are not concave the way ROC curves are, so linearly
    interpolating between two precision points can claim an achievable
    precision that was never actually observed at any threshold (Davis &
    Goadrich 2006); this avoids that by construction.
    """
    r = curve["recall"].to_numpy()
    p = curve["precision"].to_numpy()
    d_r = np.diff(r, prepend=0.0)
    valid = ~np.isnan(p)
    return float(np.sum(d_r[valid] * p[valid]))


def curves_for_model(
    sweeps: dict[tuple, pd.DataFrame], model: Model, universe: list[str],
    labels: pd.Series, gt: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    scores = substance_scores(sweeps[model.sweep_key], model, universe, gt)
    curve, n_pos, n_neg = _rates_at_thresholds(scores, labels.reindex(universe))
    summary = {
        "model": model.id, "label": model.label, "n_pos": n_pos, "n_neg": n_neg,
        "roc_auc": _roc_auc(curve), "pr_auc": _average_precision(curve),
    }
    curve.insert(0, "model", model.id)
    return curve, summary


def quarter_labels(sweep: pd.DataFrame, gt: pd.DataFrame) -> pd.Series:
    """1 for a (substance, as_of) row inside a labelled positive's known true
    interval, 0 otherwise -- the grain `benchmark.md`'s `mean_recall`/
    `mean_precision` (via `ground_truth.interval_overlap`'s `gt_quarters`/
    `det_quarters`) actually measure, unlike `substance_labels`'s coarser
    per-substance yes/no."""
    windows = _positive_windows(gt)
    q_ord = sweep["as_of"].map(_q_ordinal)
    label = pd.Series(0, index=sweep.index)
    for sub, window in windows.items():
        sel = sweep["substance"] == sub
        label.loc[sel] = q_ord.loc[sel].isin(window).astype(int)
    return label


def quarter_curves_for_model(
    sweeps: dict[tuple, pd.DataFrame], model: Model, gt: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    """The quarter-grain analogue of `curves_for_model`: one instance per
    (substance, as_of) row the model actually scored, rather than one per
    substance -- comparable to `benchmark.run`'s `mean_recall`/
    `mean_precision`, not to `curves_for_model`'s AUCs.

    Two differences from the substance-grain curve worth knowing before
    reading these numbers against each other:

    - **N is not fixed across models.** A sparse sweep (`treescan` covers
      only 95 of 211 substances) contributes fewer rows than a dense one
      (`eb05`'s ~40% substance-quarter panel coverage per `sweep_eb05`) --
      unlike the substance-grain curve, which floors every model onto the
      same fixed 211-substance panel. `n_pos`/`n_neg` in the summary vary by
      model; read them alongside the AUC, not past them.
    - **Quarters within one alarm episode are correlated, not independent.**
      `ground_truth.precision_rate`'s docstring already flags this for a
      Wilson CI specifically (four quarters of one sustained mistake are one
      trial, not four) -- that concern is about the width of an *interval*
      around a proportion, and does not invalidate this curve or its AUC as
      a description of the achievable threshold trade-off, but it means
      these AUCs should not be read with the confidence a genuinely
      independent-trials sample size would imply.
    """
    sweep = sweeps[model.sweep_key]
    scores = model_score(sweep, model)
    labels = quarter_labels(sweep, gt)
    curve, n_pos, n_neg = _rates_at_thresholds(scores, labels)
    summary = {
        "model": model.id, "label": model.label, "n_pos": n_pos, "n_neg": n_neg,
        "roc_auc": _roc_auc(curve), "pr_auc": _average_precision(curve),
    }
    curve.insert(0, "model", model.id)
    return curve, summary


def _plot(curve_tab: pd.DataFrame, summary_tab: pd.DataFrame,
         out_dir: Path, grain: str = "substance") -> None:
    """ROC and PR side by side, one line per model -- the two curves answer
    different questions (`benchmark.py`'s own precision-vs-recall framing
    applies here too) and neither alone is the whole picture, so they are
    drawn together rather than as separate figures a reader has to reconcile
    by hand.
    """
    order = summary_tab["model"].tolist()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    for i, mid in enumerate(order):
        c = curve_tab[curve_tab["model"] == mid]
        s = summary_tab[summary_tab["model"] == mid].iloc[0]
        color = _PALETTE[i % len(_PALETTE)]
        ax1.plot(c["fpr"], c["recall"], color=color, lw=1.6,
                 label=f"{mid} (AUC {s['roc_auc']:.2f})")
        pr = c.dropna(subset=["precision"])
        ax2.plot(pr["recall"], pr["precision"], color=color, lw=1.6,
                 label=f"{mid} (AP {s['pr_auc']:.2f})")

    ax1.plot([0, 1], [0, 1], color=MUTED, lw=1, ls="--")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1.02)
    # Quarter grain's n_pos/n_neg vary per model (uneven sweep coverage --
    # see `quarter_curves_for_model`'s docstring), so the fixed-panel count
    # in the substance-grain label would be misleading here.
    if grain == "substance":
        ax1.set_xlabel(f"false positive rate (of {int(summary_tab['n_neg'].iloc[0])} "
                       "non-emergence substances)")
        ax1.set_ylabel(f"recall (of {int(summary_tab['n_pos'].iloc[0])} known "
                       "emergences)")
    else:
        ax1.set_xlabel("false positive rate (substance-quarters; n varies "
                       "by model, see auc_summary.csv)")
        ax1.set_ylabel("recall (true-emergence substance-quarters)")
    ax1.set_title("ROC", fontsize=10, loc="left", color=INK)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.legend(frameon=False, fontsize=7, loc="lower right")

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1.02)
    ax2.set_xlabel("recall")
    ax2.set_ylabel("precision")
    ax2.set_title("Precision-recall", fontsize=10, loc="left", color=INK)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.legend(frameon=False, fontsize=7, loc="lower left")

    grain_desc = ("all models on the same 211-substance panel" if grain == "substance"
                 else "each model on its own scored (substance, quarter) rows")
    fig.suptitle(f"{grain.capitalize()}-level detection curves -- one point "
                f"per possible alarm threshold, {grain_desc}",
                fontsize=10, color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = out_dir / "curves.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    typer.echo(f"wrote {path}")


@app.command("run")
def run(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    ground_truth: Path = typer.Option(KNOWN_EMERGENCES_PATH, "--ground-truth"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    models: str = typer.Option(",".join(DEFAULT_MODELS), "--models",
                               help="comma-separated model ids"),
    top: int = typer.Option(None, "--top",
                            help="plot only the top N models by pr_auc -- "
                                 "the CSVs always keep every scored model; "
                                 "unset plots all of --models, which is "
                                 "illegible much past ~8 lines"),
    grain: str = typer.Option("substance", "--grain",
                              help="'substance' (default): one Bernoulli "
                                   "trial per substance, coarse recall (ever "
                                   "caught any part of it). 'quarter': one "
                                   "trial per (substance, as_of) row, the "
                                   "grain benchmark.md's mean_recall/"
                                   "mean_precision actually measure."),
) -> None:
    """ROC and PR curves + AUCs for the chosen detectors, at substance or
    quarter grain -- see `--grain`'s help and the module docstring for what
    the two answer differently."""
    if grain not in ("substance", "quarter"):
        raise typer.BadParameter("--grain must be 'substance' or 'quarter'")
    chosen = [MODELS_BY_ID[m.strip()] for m in models.split(",") if m.strip()]
    typer.echo(f"curves ({grain} grain) for {len(chosen)} models "
               f"({recent_quarters}q recent vs {baseline_quarters}q baseline, "
               f"cutoff {cutoff})")
    sweeps, denom = build_sweeps(chosen, mentions, cutoff, recent_quarters,
                                 baseline_quarters)
    gt = load_ground_truth(ground_truth, max(denom.index))

    curves, summaries = [], []
    if grain == "substance":
        universe = substance_universe(mentions, cutoff)
        labels = substance_labels(gt, universe)
        for m in chosen:
            c, s = curves_for_model(sweeps, m, universe, labels, gt)
            curves.append(c)
            summaries.append(s)
    else:
        for m in chosen:
            c, s = quarter_curves_for_model(sweeps, m, gt)
            curves.append(c)
            summaries.append(s)

    curve_tab = pd.concat(curves, ignore_index=True)
    summary_tab = pd.DataFrame(summaries).sort_values("pr_auc", ascending=False)

    out_dir.mkdir(parents=True, exist_ok=True)
    curve_tab.to_csv(out_dir / "curve_points.csv", index=False)
    summary_tab.to_csv(out_dir / "auc_summary.csv", index=False)

    if grain == "substance":
        typer.echo(f"\n{int(labels.sum())} positive substances / "
                   f"{int((labels == 0).sum())} negative, out of {len(universe)} "
                   "in the panel\n")
    else:
        typer.echo(f"\nn_pos/n_neg vary by model coverage -- see "
                   "auc_summary.csv\n")
    typer.echo(summary_tab[["model", "n_pos", "n_neg", "roc_auc", "pr_auc"]]
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    typer.echo(f"\nwrote {out_dir / 'curve_points.csv'} and "
               f"{out_dir / 'auc_summary.csv'}")

    plot_summary = summary_tab.head(top) if top is not None else summary_tab
    if top is not None and top < len(summary_tab):
        typer.echo(f"\nplotting the top {top} of {len(summary_tab)} models "
                   "by pr_auc (full ranking is still in auc_summary.csv)")
    plot_curves = curve_tab[curve_tab["model"].isin(plot_summary["model"])]
    _plot(plot_curves, plot_summary, out_dir, grain=grain)


if __name__ == "__main__":
    app()
