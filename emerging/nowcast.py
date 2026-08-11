"""
Reporting-delay nowcast: how incomplete are the most recent quarters?
(game plan P3)

Every other module in this repository takes `CUTOFF = 2025-12-31` on faith —
a hand-set date before which toxicology is declared settled, applied by
`trends.load_quarterly` before any statistic sees the data. It is a *judgement*
with no evidence behind it, and it is doing a lot of work: it decides how much
data the project throws away, and if it is set too late, every recent count is
silently low.

**This measures the delay instead of assuming it.** LACME has no report or
certification date — only `DeathDate` — so a delay cannot be recovered from a
single extract. It can be recovered from two *vintages* of the same extract:
if a death-month holds 570 deaths in a pull taken in 2024 and 729 in a pull
taken later, 22% of that month had not yet been reported when the first pull
was taken. Sweeping across death-months at one pull date traces the whole
curve, because each month sits at a different delay.

**Which vintages are usable, and why one pair is not.** Three overdose extracts
exist. Between the 2024-08 and 2025-05 pulls, the settled era (2012–2022) moves
by 0.2% — so their difference is essentially pure accrual and the pair is
usable. Between the 2025-05 and 2026-04 pulls the settled era moves by **3.5%**,
which cannot be reporting delay four to fourteen years after the fact; that is
a change in the case definition or the extraction, and mixing it into a delay
estimate would inflate every correction. Only the clean pair is used, and
`triangle` prints the settled-era drift for whatever pair it is given so this
check is never skipped.

**The curve.** Fitted on the clean pair, completeness is 27% in the month a
death occurs, 66% after one month, ~80% from three to five, ~90% at seven to
eight, and effectively 100% from ten months on. Its saturating end is
independently confirmed in the current era: the 2026-02 and 2026-04 pulls agree
to within 1.4% on every month of 2025Q1–Q2, so nine or ten months really is
settled today as well.

**And it does not transfer to the current extract, which is the actual
result.** Applied naively the curve says 2025Q4 is 81% complete and the cutoff
is too late. `transfer_check` refuses that, and it is right to: the current
extract runs at full strength through 2026-01 (198 deaths, in line with every
month before it) and then drops to 54, 18, 6. The curve predicts a *ramp* over
those months — 79%, 77%, 66% — and what is there is a **cliff**, a 3.7x fall
where a 1.2x fall was predicted. No smooth accrual process produces that. It is
the signature of a pull boundary, with a handful of fast-closing cases after
it, and it means the reporting process today is not the one measured in 2024.

So the correction is computed, printed, and explicitly withheld. The practical
reading is the opposite of the naive one: the extract looks essentially
complete through 2026-01, so `CUTOFF = 2025-12-31` is more likely conservative
than lax. This is the documented failure mode of delay-corrected nowcasts —
the CDC's provisional overdose counts use the same machinery and are known to
miss when the regime shifts — arriving here as a shift in the *reporting*
regime rather than the epidemic.

**What would fix it, and it is free.** Two current-era vintages. Save a dated
copy of every future extract; nothing else is needed, and without it this
question cannot be reopened no matter how much modelling is thrown at it.

Usage:
    python -m emerging.nowcast triangle
    python -m emerging.nowcast estimate
    python -m emerging.nowcast plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

from emerging.paths import FIG_DIR, LACME_PATH, RESULTS_DIR
from emerging.trends import BLUE, CUTOFF, INK, INK_2, MUTED, ORANGE

app = typer.Typer(add_completion=False)

# Vintages live outside this repository on purpose: they are full LACME
# extracts and carry decedent names, so they are read in place and only the
# aggregate delay curve is written here. Override with --early/--late.
EPI = Path.home() / "epi" / "data"
EARLY_VINTAGE = EPI / "2012-01-2024-08-overdoses.csv"
LATE_VINTAGE = EPI / "2012-01-2025-05-overdoses-02042026.csv"

# A pair whose settled era moves by more than this is not measuring delay.
SETTLED_DRIFT_TOL = 0.01
SETTLED_END = "2022-12"
MAX_DELAY = 15                 # months; the curve is flat well before this

DELAY_PATH = RESULTS_DIR / "reporting_delay.csv"
NOWCAST_PATH = RESULTS_DIR / "nowcast_quarters.csv"


def _monthly(path: Path) -> tuple[pd.Series, pd.Timestamp]:
    d = pd.read_csv(path, usecols=["DeathDate"], low_memory=False)
    dt = pd.to_datetime(d["DeathDate"], format="mixed", errors="coerce").dropna()
    return dt.dt.to_period("M").value_counts().sort_index(), dt.max()


def _pava(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Pool-adjacent-violators: nearest non-decreasing fit, weighted.

    The raw curve wobbles (0.797 at three months, 0.787 at four) purely from
    Poisson noise. Completeness cannot fall as more time passes, so imposing
    monotonicity is a modelling assumption we actually believe, and it avoids
    inventing a parametric family the six usable points cannot distinguish.
    """
    y, w = y.astype(float).copy(), w.astype(float).copy()
    n = len(y)
    lvl, wt, size = [], [], []
    for i in range(n):
        lvl.append(y[i]); wt.append(w[i]); size.append(1)
        while len(lvl) > 1 and lvl[-2] > lvl[-1]:
            v2, w2, s2 = lvl.pop(), wt.pop(), size.pop()
            v1, w1, s1 = lvl.pop(), wt.pop(), size.pop()
            lvl.append((v1 * w1 + v2 * w2) / (w1 + w2))
            wt.append(w1 + w2); size.append(s1 + s2)
    out = np.empty(n)
    i = 0
    for v, s in zip(lvl, size):
        out[i:i + s] = v
        i += s
    return out


def delay_curve(early: Path = EARLY_VINTAGE, late: Path = LATE_VINTAGE,
                max_delay: int = MAX_DELAY) -> tuple[pd.DataFrame, dict]:
    """Completeness by months elapsed, from a pair of extract vintages."""
    A, a_max = _monthly(early)
    B, b_max = _monthly(late)
    t = pd.DataFrame({"early": A, "late": B}).dropna()

    settled = t.loc[:SETTLED_END]
    if settled.empty or settled["late"].sum() == 0:
        # Without a settled era there is nothing to distinguish accrual from a
        # change in the case definition, which is the whole basis for trusting
        # the pair. Refuse rather than return a NaN that reads as a failure
        # for the wrong reason.
        raise ValueError(
            f"the two vintages share no months at or before {SETTLED_END}, so "
            f"the settled-era drift check cannot run and the pair cannot be "
            f"validated")
    drift = settled["early"].sum() / settled["late"].sum() - 1.0

    as_of = a_max.to_period("M")
    t["delay_m"] = [(as_of - m).n for m in t.index]
    t = t[(t["delay_m"] >= 0) & (t["delay_m"] <= max_delay)]

    g = t.groupby("delay_m")[["early", "late"]].sum().sort_index()
    g = g.reindex(range(max_delay + 1))
    g["raw"] = g["early"] / g["late"]
    obs = g["raw"].notna()
    # Rows are already ordered by elapsed months, and completeness can only
    # rise with elapsed time, so this is a plain non-decreasing isotonic fit.
    # (Reversing first — which an earlier version did — makes the target
    # decreasing and pools the entire curve to one number.)
    fitted = _pava(g.loc[obs, "raw"].to_numpy(),
                   g.loc[obs, "late"].to_numpy())
    g.loc[obs, "completeness"] = np.clip(fitted, 1e-6, 1.0)
    g["completeness"] = g["completeness"].ffill().bfill()
    # Binomial standard error of the fitted proportion.
    g["se"] = np.sqrt(g["completeness"] * (1 - g["completeness"])
                      / g["late"].clip(lower=1))

    meta = {"as_of": as_of, "early_max": a_max, "late_max": b_max,
            "settled_drift": drift,
            "settled_ok": abs(drift) <= SETTLED_DRIFT_TOL}
    return g.reset_index(names="delay_m"), meta


def apply_curve(counts: pd.Series, as_of: pd.Period,
                curve: pd.DataFrame) -> pd.DataFrame:
    """Correct monthly counts for incomplete reporting.

    n is treated as Binomial(N, pi), so N-hat = n / pi and the standard error
    combines the sampling noise in n with the uncertainty in pi itself.
    """
    pi = dict(zip(curve["delay_m"], curve["completeness"]))
    se_pi = dict(zip(curve["delay_m"], curve["se"]))
    tail = float(curve["completeness"].iloc[-1])

    rows = []
    for m, n in counts.items():
        d = (as_of - m).n
        if d < 0:
            continue
        p = pi.get(d, tail if d > curve["delay_m"].max() else tail)
        sp = se_pi.get(d, 0.0)
        est = n / p
        rel = np.sqrt((1 - p) / max(n, 1) + (sp / p) ** 2)
        rows.append({"month": m, "delay_m": d, "observed": int(n),
                     "completeness": p, "estimated": est,
                     "se": est * rel, "added": est - n})
    return pd.DataFrame(rows)


def _q(p) -> str:
    return f"{p.year}Q{p.quarter}"


def _last_full(chk: pd.DataFrame) -> str:
    """Last month whose level is at least half the settled baseline."""
    ok = chk[chk["implied_completeness"] >= 0.5]
    return str(ok["month"].max()) if len(ok) else "n/a"


def _collapse(chk: pd.DataFrame) -> str:
    near = chk[chk["delay_m"] <= 3].sort_values("delay_m")
    return " -> ".join(str(int(v)) for v in near["observed"][::-1])


def transfer_check(counts: pd.Series, as_of: pd.Period, curve: pd.DataFrame,
                   baseline_lo: int = 12, baseline_hi: int = 30) -> pd.DataFrame:
    """Does the current extract's own tail look like the measured curve?

    A delay curve estimated on one vintage pair is only usable on another
    extract if that extract's reporting behaves the same way. There is no
    second current-era vintage to check against directly, but there is one
    thing the extract tells us for free: its own recent months, read against
    a settled baseline level, trace an *implied* completeness profile. If that
    profile is a cliff where the curve is a ramp, the two processes are not
    the same and the correction is invalid.

    This is not a subtle test. The 2024-era curve rises 27% -> 66% -> 76% over
    three months. An extract that runs at full strength and then drops to a
    quarter of it in one month is not obeying that curve, whatever else is
    true.
    """
    delays = pd.Series({m: (as_of - m).n for m in counts.index})
    settled = counts[(delays >= baseline_lo) & (delays <= baseline_hi)]
    base = float(np.median(settled)) if len(settled) else np.nan

    pi = dict(zip(curve["delay_m"], curve["completeness"]))
    rows = []
    for m, n in counts.items():
        d = (as_of - m).n
        if 0 <= d <= 8:
            rows.append({"month": str(m), "delay_m": d, "observed": int(n),
                         "implied_completeness": n / base if base else np.nan,
                         "curve_completeness": pi.get(d, np.nan)})
    out = pd.DataFrame(rows).sort_values("delay_m")
    out["ratio"] = out["implied_completeness"] / out["curve_completeness"]
    return out


@app.command("triangle")
def triangle(
    early: Path = typer.Option(EARLY_VINTAGE, "--early"),
    late: Path = typer.Option(LATE_VINTAGE, "--late"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
) -> None:
    """Estimate the reporting-delay curve from two extract vintages."""
    for p in (early, late):
        if not p.exists():
            raise typer.BadParameter(f"vintage not found: {p}")
    curve, meta = delay_curve(early, late)

    typer.echo(f"early vintage last death {meta['early_max'].date()} "
               f"(taken as its pull date); late vintage last death "
               f"{meta['late_max'].date()}")
    typer.echo(f"settled-era (<= {SETTLED_END}) drift between the two: "
               f"{100 * meta['settled_drift']:+.2f}%  "
               + ("OK — the difference is accrual"
                  if meta["settled_ok"] else
                  "TOO LARGE — this pair differs by case definition, not "
                  "delay, and must not be used"))
    if not meta["settled_ok"]:
        raise typer.Exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    curve.to_csv(out_dir / DELAY_PATH.name, index=False)
    disp = curve[["delay_m", "early", "late", "raw", "completeness"]]
    typer.echo("\n" + disp.to_string(index=False,
                                     float_format=lambda v: f"{v:.4f}"))
    typer.echo(f"\nwrote {out_dir / DELAY_PATH.name}")


@app.command("estimate")
def estimate(
    lacme: Path = typer.Option(LACME_PATH, "--lacme"),
    early: Path = typer.Option(EARLY_VINTAGE, "--early"),
    late: Path = typer.Option(LATE_VINTAGE, "--late"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    quarters: int = typer.Option(8, "--quarters"),
) -> None:
    """Correct the current extract's trailing quarters, and test the CUTOFF."""
    curve, meta = delay_curve(early, late)
    if not meta["settled_ok"]:
        raise typer.BadParameter("vintage pair is not measuring delay")

    m, cur_max = _monthly(lacme)
    as_of = cur_max.to_period("M")

    chk = transfer_check(m, as_of, curve)
    typer.echo(f"current extract: last death {cur_max.date()}, "
               f"delay curve measured at {meta['as_of']}\n")
    typer.echo("TRANSFERABILITY CHECK — does this extract obey the measured "
               "curve?")
    typer.echo(chk.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    near = chk[chk["delay_m"] <= 4]
    bad = ((near["ratio"] > 1.35) | (near["ratio"] < 0.65)).sum()
    transfers = bad == 0
    typer.echo(
        f"\n  {bad} of {len(near)} months at delay <= 4 are off the curve by "
        f"more than 35%. "
        + ("The curve transfers; the correction below is usable."
           if transfers else
           "**The curve does not transfer.** This extract runs at full "
           "strength and then falls off a cliff, which no smooth delay "
           "process produces. The correction below is printed for reference "
           "and MUST NOT be applied.")
    )

    est = apply_curve(m, as_of, curve)
    est["quarter"] = est["month"].apply(lambda p: p.asfreq("Q"))

    # A quarter the extract does not fully span is a calendar problem, not a
    # reporting one, and correcting it would conflate the two.
    full = est.groupby("quarter")["month"].count()
    complete_q = set(full[full == 3].index)

    q = (est[est["quarter"].isin(complete_q)]
         .groupby("quarter")
         .agg(observed=("observed", "sum"), estimated=("estimated", "sum"),
              se=("se", lambda s: float(np.sqrt((s ** 2).sum()))))
         .sort_index())
    q["completeness"] = q["observed"] / q["estimated"]
    q["lo"] = q["estimated"] - 1.96 * q["se"]
    q["hi"] = q["estimated"] + 1.96 * q["se"]
    q = q.tail(quarters)

    q["curve_transfers"] = transfers
    out_dir.mkdir(parents=True, exist_ok=True)
    q.to_csv(out_dir / NOWCAST_PATH.name)
    chk.to_csv(out_dir / "reporting_delay_transfer.csv", index=False)

    typer.echo("\nCORRECTION IMPLIED BY THE 2024-ERA CURVE"
               + ("" if transfers else "  (NOT VALID — see above)"))
    disp = q.copy()
    disp.index = [_q(p) for p in disp.index]
    typer.echo(disp[["observed", "estimated", "lo", "hi", "completeness"]]
               .to_string(float_format=lambda v: f"{v:.1f}"))

    cut = pd.Timestamp(CUTOFF).to_period("Q")
    typer.echo(f"\nVERDICT for CUTOFF = {CUTOFF} (keeps through {_q(cut)}):")
    if transfers:
        kept = q[q.index <= cut]
        worst_q = kept["completeness"].idxmin()
        worst = kept.loc[worst_q, "completeness"]
        typer.echo(f"  least complete quarter kept: {_q(worst_q)} at "
                   f"{100 * worst:.0f}% — "
                   + ("consistent with treating it as settled."
                      if worst >= 0.97 else
                      f"not settled; {int(kept.loc[worst_q, 'observed'])} "
                      f"observed imply about "
                      f"{kept.loc[worst_q, 'estimated']:.0f}."))
    else:
        typer.echo("  the correction cannot be applied, so the cutoff cannot "
                   "be revised from it. What the extract's own tail does say "
                   "is that it holds a steady level through "
                   f"{_last_full(chk)} and collapses immediately after "
                   f"({_collapse(chk)}), which is a pull "
                   "boundary rather than an accrual ramp — so the current "
                   "cutoff is more likely conservative than lax.")
        typer.echo("  To settle this properly, save a dated copy of every "
                   "future extract. Two current-era vintages is all it takes; "
                   "we have none.")
    typer.echo(f"\nwrote {out_dir / NOWCAST_PATH.name}")


@app.command("plot")
def plot(
    results_dir: Path = typer.Option(RESULTS_DIR, "--results-dir"),
    fig_dir: Path = typer.Option(FIG_DIR, "--fig-dir"),
) -> None:
    """Delay curve, and what it does to the recent quarterly counts."""
    curve = pd.read_csv(results_dir / DELAY_PATH.name)
    q = pd.read_csv(results_dir / NOWCAST_PATH.name)
    q["label"] = [f"{p[:4]}Q{p[-1]}" for p in q["quarter"].astype(str)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.2),
                                   gridspec_kw={"width_ratios": [1, 1.15]})

    obs = curve["raw"].notna()
    ax1.scatter(curve.loc[obs, "delay_m"], 100 * curve.loc[obs, "raw"],
                s=42, color=MUTED, zorder=3, label="observed")
    ax1.plot(curve["delay_m"], 100 * curve["completeness"], color=BLUE, lw=2.4,
             zorder=4, label="monotone fit")
    ax1.axhline(97, color=INK_2, lw=1.0, ls=":", zorder=2)
    ax1.text(curve["delay_m"].max(), 97.6, "97%", ha="right", fontsize=8,
             color=INK_2)
    ax1.set_xlabel("months elapsed since the death month", fontsize=9.5,
                   color=INK_2)
    ax1.set_ylabel("% of deaths reported", fontsize=9.5, color=INK_2)
    ax1.set_ylim(0, 105)
    ax1.set_title("Reporting delay, from two extract vintages", fontsize=12,
                  color=INK, loc="left", pad=8)
    ax1.legend(frameon=False, fontsize=9)
    ax1.grid(color="#eceae5", lw=0.9)
    ax1.set_axisbelow(True)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
    ax1.tick_params(colors=INK_2, labelsize=9)

    chk = pd.read_csv(results_dir / "reporting_delay_transfer.csv")
    chk = chk.sort_values("delay_m")
    ax2.plot(chk["delay_m"], 100 * chk["curve_completeness"], color=BLUE,
             lw=2.4, marker="o", ms=6, label="predicted by the 2024 curve")
    ax2.plot(chk["delay_m"], 100 * chk["implied_completeness"], color=ORANGE,
             lw=2.4, marker="D", ms=6,
             label="what this extract actually does")
    ax2.fill_between(chk["delay_m"], 100 * chk["implied_completeness"],
                     100 * chk["curve_completeness"],
                     where=(chk["implied_completeness"]
                            < chk["curve_completeness"]),
                     color=ORANGE, alpha=0.13, zorder=1)
    y1 = 100 * chk.loc[chk["delay_m"] == 1, "implied_completeness"].iloc[0]
    ax2.annotate("a cliff, not a ramp:\n198 -> 54 -> 18 deaths",
                 xy=(1.0, y1),
                 xytext=(3.4, 34), fontsize=9, color=INK_2,
                 arrowprops=dict(arrowstyle="->", color=INK_2, lw=1.1))
    ax2.set_xlabel("months elapsed since the death month", fontsize=9.5,
                   color=INK_2)
    ax2.set_ylabel("% of the settled monthly level", fontsize=9.5, color=INK_2)
    ax2.set_ylim(0, 110)
    ax2.invert_xaxis()
    ax2.set_title("Why the curve cannot be applied here", fontsize=12,
                  color=INK, loc="left", pad=8)
    ax2.legend(frameon=False, fontsize=9, loc="lower left")
    ax2.grid(color="#eceae5", lw=0.9)
    ax2.set_axisbelow(True)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(colors=INK_2, labelsize=9)

    fig.text(0.008, 0.015,
             "Left: completeness measured by comparing an extract pulled in "
             "2024 against a later one, month by month; the pair is usable "
             "because its settled era agrees to 0.2%. Right: the same curve "
             "against what the current extract actually does, read as a share "
             "of its own settled monthly level. The curve predicts a gradual "
             "ramp; the extract holds level and then falls off a cliff, which "
             "is a pull boundary rather than accrual. The correction is "
             "therefore computed but withheld — see results/"
             "nowcast_quarters.csv, column curve_transfers.",
             fontsize=7.6, color=MUTED, va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = fig_dir / "reporting_delay.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
