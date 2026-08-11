"""
Rank and plot substances whose share of LA County overdose deaths is rising.

Detection statistic. The candidates that matter are semi-rare — xylazine has 9
lifetime mentions, nitazenes 3 — and at those counts a raw ratio of recent to
baseline share is worthless: one extra death doubles it. So this uses the
pharmacovigilance answer, a gamma-Poisson empirical-Bayes shrinkage:

    n_s ~ Poisson(lambda_s * E_s),  lambda_s ~ Gamma(alpha, beta)

where `E_s` is the count expected in the recent window if substance `s` held its
baseline share of all overdose deaths. alpha/beta are fitted across substances by
maximum likelihood (the marginal is negative binomial), giving a posterior
Gamma(alpha + n_s, beta + E_s) per substance. Ranking is by **EB05**, its 5th
percentile — the conservative end, so a substance only ranks high when the data
cannot easily be explained by noise. `ebgm`, the posterior *geometric* mean
exp(E[log lambda]), is reported alongside it.

This is a single-component GPS, not full MGPS. DuMouchel's two-component
mixture was implemented first and is not identifiable at this scale — see
`_fit_gps_prior` for what it did instead. It becomes worth revisiting on a
substance x zip x quarter table.

Does it work? `backtest` rolls the as-of date back and re-scores, so the
ranking is falsifiable rather than merely quiet. All five known LA emergences
reach rank 1-4 in the quarter they emerge, xylazine included at n=6.

Right-truncation. Toxicology takes months, so recent quarters are undercounts
and every rising substance looks like it is falling. Two guards, both applied in
`load_quarterly` before any statistic runs: `--cutoff` (default 2025-12-31)
drops quarters whose toxicology is not yet settled, and partial calendar
quarters at the end of the extract are always dropped. `--lag-quarters` is the
older rolling-haircut alternative and defaults to 0, because stacking it on the
cutoff double-corrects. Novel substances confirm *slower* than fentanyl, so even
this is a floor, not a fix — see the README.

Windowing matters more than it looks. The recent/baseline split is a fixed
number of quarters back from the cutoff, so moving the cutoff slides the window.
For a substance whose rise is a short burst rather than a sustained trend, that
can move EB05 across the alarm line — lidocaine does exactly this. Quote the
`backtest` peak alongside the current value when the trajectory is spiky.

Novelty. `date_added` from NFLIS is *catalog* novelty (when labs could report a
substance), left-censored at the 1998-10 inception for 27% of the catalog. Local
epidemiological novelty is `first_seen`, the first LACME quarter — itself
left-censored at 2012, where the extract begins. The scatter plots the two
against each other, because the interesting region is where they disagree: old
chemistry arriving in LA for the first time.

Four figures, written to results/figures/:
    rising_candidates.png  small multiples of the current top-ranked substances
    novelty_vs_growth.png  chemical novelty vs LA arrival
    share_heatmap.png      substance x quarter share of all OD deaths
    known_emergences.png   the backtest set

Usage:
    python -m emerging.trends rank
    python -m emerging.trends rank --recent-quarters 8 --top 30
    python -m emerging.trends backtest
    python -m emerging.trends plot
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
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import gamma as gamma_dist

from emerging.extract import MENTIONS_PATH, load_cohort
from emerging.paths import FIG_DIR, RESULTS_DIR

app = typer.Typer(add_completion=False)


# Validated categorical palette (dataviz skill reference instance; checked with
# scripts/validate_palette.js --mode light: all checks pass, contrast WARN
# relieved by direct labels on every mark).
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#b8b7b2"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f2f6fc", "#2a78d6", "#12365f"])

# Hard completeness cutoff: quarters ending after this date are dropped before
# anything is computed. Toxicology through end-2025 is considered settled;
# 2026 is not. This is a *data judgement*, not a derived quantity — advance it
# by hand when the LACME extract is refreshed and the newer quarters have had
# time to adjudicate.
CUTOFF = "2025-12-31"

# Additional reporting-lag haircut, in quarters, applied on top of the cutoff.
# Defaults to 0 because `CUTOFF` already encodes the settledness judgement —
# these are two ways of expressing the same correction, and stacking them
# discards a year of settled data. Set it only if you drop the cutoff.
LAG_QUARTERS = 0
RECENT_QUARTERS = 4
BASELINE_QUARTERS = 8

# Dated LA-local emergences used as the backtest ground truth. Only
# para-fluorofentanyl is well powered; the rest are shown for honesty about how
# thin the evidence is at LA-county scale.
KNOWN_EMERGENCES = [
    "para-Fluorofentanyl", "Bromazolam", "Mitragynine", "Carfentanil", "Xylazine",
]

# EB05 above which a substance is "in alarm". Read off the backtest: every
# known LA emergence clears it in the quarter it emerges, and it is high enough
# that the bulk of the catalog never does — 24 substances of 219 have ever
# breached it across 45 as-of quarters. It is a calibrated reading line, not a
# significance test; EB05 > 1 already means "credible growth".
ALARM_THRESHOLD = 1.5
ALARM_PATH = RESULTS_DIR / "alarm_history.csv"

# Substances driving the national emerging-drug conversation, checked against
# LA regardless of whether the ranking surfaces them. The ranking can only
# promote what is already in the data; this asks the complementary question,
# "is the thing everyone else is seeing here at all", and a zero is as much a
# result as a hit. Every entry resolves in the lexicon, so a zero means the
# scan looked and found nothing — not that the name was unrecognised.
WATCHLIST: dict[str, str] = {
    "Xylazine": "'tranq' — fentanyl adulterant, widespread in the Northeast",
    "Medetomidine": "vet sedative displacing xylazine in East Coast fentanyl, 2024-",
    "Dexmedetomidine": "medetomidine isomer, same supply signal",
    "BTMPS (Bis(2,2,6,6-tetramethyl-4-piperidinyl) Sebacate)":
        "industrial plastics stabiliser found in fentanyl from mid-2024",
    "Tianeptine": "'gas-station heroin', unscheduled atypical opioid",
    "Bromazolam": "designer benzo, the successor to etizolam/flualprazolam",
    "Carfentanil": "resurgent 2023- after a near-decade absence",
    "Levamisole": "long-standing cocaine adulterant",
    "Desomorphine": "'krokodil' — recurrent media claim, rarely confirmed",
}
# Handled as a family: individually every nitazene is too rare to score, and
# the public-health question is about the class, not the analog.
NITAZENE_RE = r"nitaz"


def load_quarterly(
    mentions_path: Path = MENTIONS_PATH, rollup: bool = True,
    cutoff: str | None = CUTOFF,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (substance x quarter counts, all-OD deaths per quarter).

    Metabolites are rolled up to their parent so NORFENTANYL does not read as
    its own rising substance; class terms are dropped.

    Two completeness filters are applied here, before any statistic sees the
    data. `cutoff` drops quarters ending after a settled-toxicology date (see
    `CUTOFF`). Separately, trailing quarters the extract does not even span are
    always dropped, whatever the cutoff.
    """
    m = pd.read_parquet(mentions_path)
    m = m[~m["is_class_term"].fillna(False)]
    key = "rollup" if rollup else "canonical"
    m[key] = m[key].fillna(m["canonical"])

    # One row per (case, substance) after rollup — a death naming both fentanyl
    # and norfentanyl must count once.
    m = m.drop_duplicates(subset=["CaseNumber", key])
    counts = (m.groupby([key, "quarter"]).size()
                .rename("n").reset_index().rename(columns={key: "substance"}))

    cohort = load_cohort()
    denom = cohort.groupby("quarter").size().rename("N")

    # Drop trailing quarters the extract does not even span. The file ends
    # 2026-04-04, so 2026Q2 is a 4-day stub holding 3 deaths. This is a
    # different defect from reporting lag and has to be removed separately:
    # left in, it silently consumes one of the `lag_quarters` slots, so a
    # 2-quarter lag correction only discards one real quarter of data.
    last_death = cohort["death_date"].max()
    complete = [q for q in denom.index
                if q + pd.offsets.QuarterEnd(0) <= last_death]
    dropped = [q for q in denom.index if q not in set(complete)]
    if dropped:
        typer.echo("dropped partial calendar quarter(s): "
                   + ", ".join(f"{q.year}Q{q.quarter} (n={int(denom[q])})"
                               for q in dropped))
    denom = denom.loc[complete]

    if cutoff is not None:
        cut = pd.Timestamp(cutoff)
        keep = [q for q in denom.index if q + pd.offsets.QuarterEnd(0) <= cut]
        n_cut = len(denom) - len(keep)
        if n_cut:
            typer.echo(f"cutoff {cut.date()}: dropped {n_cut} unsettled "
                       f"quarter(s) after {keep[-1].year}Q{keep[-1].quarter}")
        denom = denom.loc[keep]
        complete = keep

    counts = counts[counts["quarter"].isin(set(complete))]
    return counts, denom


def _nb_logpmf(n: np.ndarray, e: np.ndarray, a: float, b: float) -> np.ndarray:
    """log P(n | E) for n ~ Poisson(lambda*E), lambda ~ Gamma(a, b).

    E is floored at a tiny epsilon so substances absent from the baseline
    (E == 0) score finitely: n == 0 contributes 0, and n > 0 becomes very
    surprising, which is the correct reading of "this appeared from nothing".
    """
    e = np.maximum(e, 1e-9)
    return (gammaln(n + a) - gammaln(a) - gammaln(n + 1)
            + a * np.log(b / (b + e)) + n * np.log(e / (b + e)))


# 0 disables the cap. See `_fit_gps_prior` for when to set one.
MAX_EXPECTED_FOR_FIT = 0.0


def _fit_gps_prior(
    n: np.ndarray, e: np.ndarray, max_expected: float = MAX_EXPECTED_FOR_FIT
) -> dict:
    """ML fit of a Gamma(alpha, beta) prior under a Poisson likelihood.

    The marginal of n given E is negative binomial with size alpha and success
    probability beta/(beta+E), so alpha and beta come straight out of an MLE.

    **Why one component and not DuMouchel's two.** MGPS uses a two-gamma
    mixture, and that was tried here first. It is not identifiable at this
    scale: 203 substances, of which only 111 have a non-zero baseline, cannot
    pin down five parameters. Nelder-Mead converged to a different degenerate
    solution depending on what was in the fit set — a null component of
    Gamma(682, 627) with everything in, and a point mass at 1.08 with the
    largest substances excluded. The single component is stable, has a closed
    form for every posterior quantity, and is honest about what the data
    support. The mixture becomes worth revisiting on a substance x zip x
    quarter table, where the cell count is thousands rather than hundreds.

    Substances with E == 0 are excluded from the fit (they carry no information
    about the prior) but are still scored with it. `max_expected` optionally
    also excludes the largest substances, whose counts dominate the likelihood;
    it is off by default because the standard GPS fit uses all cells and the
    uncapped fit is the more conservative of the two (alpha 7.4 vs 19.9 —
    a more diffuse prior, so less shrinkage).
    """
    keep = e > 0
    if max_expected:
        keep &= e <= max_expected
    nf, ef = n[keep], e[keep]
    if len(nf) < 10:
        return {"alpha": 0.5, "beta": 0.5, "n_fit": int(len(nf)),
                "converged": False}

    def nll(params: np.ndarray) -> float:
        a, b = np.exp(params)
        if min(a, b) < 1e-9:
            return 1e12
        ll = _nb_logpmf(nf, ef, a, b)
        return -float(np.sum(ll)) if np.isfinite(ll).all() else 1e12

    res = minimize(nll, x0=np.array([0.0, 0.0]), method="Nelder-Mead",
                   options={"maxiter": 10000, "maxfev": 10000,
                            "xatol": 1e-8, "fatol": 1e-8})
    a, b = np.exp(res.x)
    return {"alpha": float(a), "beta": float(b), "n_fit": int(len(nf)),
            "converged": bool(res.success), "nll": float(res.fun)}


def _gps_scores(n: np.ndarray, e: np.ndarray, pr: dict) -> dict:
    """EBGM and posterior quantiles from Gamma(alpha + n, beta + E).

    EBGM is the geometric mean exp(E[log lambda | n]), per DuMouchel, rather
    than the arithmetic posterior mean — at these counts the posterior is
    always skewed and the geometric mean is the stable summary.
    """
    from scipy.special import digamma

    post_a = pr["alpha"] + n
    post_b = pr["beta"] + e
    return {
        "ebgm": np.exp(digamma(post_a) - np.log(post_b)),
        "eb05": gamma_dist.ppf(0.05, post_a, scale=1.0 / post_b),
        "eb95": gamma_dist.ppf(0.95, post_a, scale=1.0 / post_b),
    }


def rank_substances(
    counts: pd.DataFrame,
    denom: pd.Series,
    lag_quarters: int = LAG_QUARTERS,
    recent_quarters: int = RECENT_QUARTERS,
    baseline_quarters: int = BASELINE_QUARTERS,
) -> tuple[pd.DataFrame, dict]:
    """Score every substance by EB-shrunken growth in share of OD deaths."""
    quarters = sorted(denom.index)
    if lag_quarters:
        quarters = quarters[:-lag_quarters]
    if len(quarters) < recent_quarters + baseline_quarters:
        raise ValueError(
            f"need {recent_quarters + baseline_quarters} complete quarters, "
            f"have {len(quarters)}"
        )
    recent_q = quarters[-recent_quarters:]
    base_q = quarters[-(recent_quarters + baseline_quarters):-recent_quarters]

    n_recent_all = float(denom.loc[recent_q].sum())
    n_base_all = float(denom.loc[base_q].sum())

    wide = counts.pivot_table(index="substance", columns="quarter",
                              values="n", aggfunc="sum", fill_value=0)
    wide = wide.reindex(columns=sorted(set(wide.columns) | set(quarters)),
                        fill_value=0)

    n_recent = wide[recent_q].sum(axis=1).astype(float)
    n_base = wide[base_q].sum(axis=1).astype(float)
    expected = n_base * (n_recent_all / n_base_all)

    nv, ev = n_recent.to_numpy(), expected.to_numpy()
    prior = _fit_gps_prior(nv, ev)
    scores = _gps_scores(nv, ev, prior)
    ebgm, eb05, eb95 = scores["ebgm"], scores["eb05"], scores["eb95"]

    seen = counts[counts["n"] > 0]
    first_seen = seen.groupby("substance")["quarter"].min()
    last_seen = seen.groupby("substance")["quarter"].max()

    # EB05 measures growth in *share*. When the denominator moves, share and
    # count disagree: LA's total OD deaths fell 22% between the current
    # windows, so a substance with flat counts gains share for free, and
    # cocaine posts credible share growth while its own deaths fall. Carry the
    # count ratio alongside so the two are never confused — `count_ratio < 1`
    # with `eb05 > 1` reads "a growing share of a shrinking total", which is a
    # true statement about supply composition and a false one about deaths.
    per_q = recent_quarters / baseline_quarters
    count_ratio = n_recent / (n_base * per_q).replace(0, np.nan)

    res = pd.DataFrame({
        "n_recent": n_recent.astype(int),
        "n_baseline": n_base.astype(int),
        "expected": expected,
        "share_recent_pct": 100 * n_recent / n_recent_all,
        "share_baseline_pct": 100 * n_base / n_base_all,
        "count_ratio": count_ratio,
        "ebgm": ebgm,
        "eb05": eb05,
        "eb95": eb95,
    })
    res["n_total"] = wide.sum(axis=1).astype(int)
    res["first_seen"] = first_seen
    res["last_seen"] = last_seen
    res = res.sort_values("eb05", ascending=False)

    meta = {
        "prior": prior,
        "recent_window": (recent_q[0], recent_q[-1]),
        "baseline_window": (base_q[0], base_q[-1]),
        "n_recent_all": n_recent_all, "n_base_all": n_base_all,
        "truncated_from": (sorted(denom.index)[-lag_quarters]
                           if lag_quarters else None),
        "quarters": quarters,
    }
    return res, meta


def _attach_novelty(res: pd.DataFrame, mentions_path: Path) -> pd.DataFrame:
    """Join NFLIS catalog novelty onto the ranking."""
    m = pd.read_parquet(mentions_path)
    key = m["rollup"].fillna(m["canonical"])
    nov = (m.assign(_k=key)
             .drop_duplicates(subset="_k")
             .set_index("_k")[["date_added", "date_added_censored", "category"]])
    return res.join(nov, how="left")


def _fmt_prior(pr: dict) -> str:
    return (f"GPS prior: Gamma(alpha={pr['alpha']:.3f}, beta={pr['beta']:.3f}), "
            f"mean {pr['alpha'] / pr['beta']:.3f}, "
            f"sd {np.sqrt(pr['alpha']) / pr['beta']:.3f}, "
            f"fitted on {pr['n_fit']} substances with a non-zero baseline"
            + ("" if pr.get("converged", True) else "  [DID NOT CONVERGE]"))


def _fmt_window(meta: dict) -> str:
    (r0, r1), (b0, b1) = meta["recent_window"], meta["baseline_window"]
    q = lambda t: f"{t.year}Q{t.quarter}"  # noqa: E731
    return (f"recent {q(r0)}–{q(r1)} vs baseline {q(b0)}–{q(b1)}; "
            f"{meta['n_recent_all']:,.0f} vs {meta['n_base_all']:,.0f} OD deaths")


@app.command("rank")
def rank(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    lag_quarters: int = typer.Option(LAG_QUARTERS, "--lag-quarters"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    top: int = typer.Option(25, "--top"),
    min_recent: int = typer.Option(2, "--min-recent",
                                   help="Suppress substances below this recent count"),
) -> None:
    """Rank rising substances by EB05 and write the table."""
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    res, meta = rank_substances(counts, denom, lag_quarters,
                                recent_quarters, baseline_quarters)
    res = _attach_novelty(res, mentions)

    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "rising_substances.csv")

    typer.echo(_fmt_window(meta))
    typer.echo(_fmt_prior(meta["prior"]))
    if meta["truncated_from"] is not None:
        t = meta["truncated_from"]
        typer.echo(f"dropped as incomplete: {t.year}Q{t.quarter} onward "
                   f"({lag_quarters} quarters)")

    show = res[res["n_recent"] >= min_recent].head(top)
    cols = ["n_recent", "n_baseline", "expected", "share_recent_pct",
            "share_baseline_pct", "count_ratio", "ebgm", "eb05", "n_total",
            "first_seen"]
    disp = show[cols].copy()
    disp["first_seen"] = disp["first_seen"].map(
        lambda t: "" if pd.isna(t) else f"{t.year}Q{t.quarter}"
    )
    typer.echo("\n" + disp.to_string(float_format=lambda v: f"{v:.2f}"))
    typer.echo(f"\nwrote {out_dir / 'rising_substances.csv'}")


@app.command("backtest")
def backtest(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    substances: str = typer.Option(",".join(KNOWN_EMERGENCES), "--substances"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
) -> None:
    """Re-run the ranking as of each historical quarter and report where the
    known emergences placed.

    Without this the ranking is unfalsifiable: a detector that never fires and
    a detector that is broken produce the same empty top-10. Rolling the
    as-of date back to the quarter each substance actually emerged answers
    "would this have caught it, and how late" — and, for LA County, mostly
    answers "there was never enough signal to catch."
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    targets = [s.strip() for s in substances.split(",") if s.strip()]
    all_q = sorted(denom.index)
    total_q = recent_quarters + baseline_quarters

    rows = []
    # Step the as-of quarter forward; at each step score only data available
    # then, so nothing after the as-of date can leak into the statistic.
    for i in range(total_q, len(all_q) + 1):
        as_of = all_q[i - 1]
        sub_denom = denom.loc[all_q[:i]]
        sub_counts = counts[counts["quarter"] <= as_of]
        try:
            res, _ = rank_substances(sub_counts, sub_denom, lag_quarters=0,
                                     recent_quarters=recent_quarters,
                                     baseline_quarters=baseline_quarters)
        except ValueError:
            continue
        ranked = res[res["n_recent"] > 0].sort_values("eb05", ascending=False)
        order = {name: r for r, name in enumerate(ranked.index, start=1)}
        for s in targets:
            if s not in res.index:
                continue
            r = res.loc[s]
            rows.append({
                "as_of": as_of, "substance": s, "rank": order.get(s),
                "n_recent": int(r["n_recent"]), "n_baseline": int(r["n_baseline"]),
                "ebgm": r["ebgm"], "eb05": r["eb05"],
            })

    bt = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    bt.to_csv(out_dir / "backtest_known_emergences.csv", index=False)

    typer.echo(f"as-of sweep: {len(all_q) - total_q + 1} quarters, "
               f"{recent_quarters}q recent vs {baseline_quarters}q baseline\n")
    for s in targets:
        g = bt[bt["substance"] == s]
        if g.empty or g["n_recent"].max() == 0:
            typer.echo(f"{s:24s} never reaches a scoreable count")
            continue
        best = g.loc[g["eb05"].idxmax()]
        top10 = g[g["rank"] <= 10]
        first_top10 = top10["as_of"].min() if not top10.empty else None
        typer.echo(
            f"{s:24s} peak EB05 {best['eb05']:.2f} (EBGM {best['ebgm']:.2f}, "
            f"n={int(best['n_recent'])}) at "
            f"{best['as_of'].year}Q{best['as_of'].quarter}, best rank "
            f"{int(g['rank'].min())}"
            + (f", first top-10 at {first_top10.year}Q{first_top10.quarter}"
               if first_top10 is not None else ", never reaches top 10")
        )
    typer.echo(f"\nwrote {out_dir / 'backtest_known_emergences.csv'}")


@app.command("regime")
def regime(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
) -> None:
    """Look for changes in how deaths are recorded, rather than in what people
    are dying of, that would move EB05.

    The statistic is a *share*: numerator the substance, denominator all
    overdose deaths, both counted from ME narratives. Three things can move it
    without any change in the drug supply, and all three are measurable here:

    1. **Denominator drift.** Total OD deaths falling makes every flat
       substance gain share. This is the live one — see `count_ratio` in the
       ranking.
    2. **Panel drift.** Substances named per death rising means the ME is
       finding more of everything, which inflates the numerator of the rare
       substances specifically.
    3. **Simultaneous onsets.** Several substances first breaching in the same
       quarter is the signature of a change in recording, not of several
       independent epidemics starting at once.

    What this cannot do is identify the *policy*. It flags the quarters where
    something structural happened; matching those to naloxone distribution,
    scheduling changes, or a lab protocol revision needs sources outside this
    dataset.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    all_q = sorted(denom.index)
    spd = (counts.groupby("quarter")["n"].sum().reindex(all_q, fill_value=0)
           / denom.reindex(all_q))

    tab = pd.DataFrame({
        "deaths": denom.reindex(all_q).astype(int),
        "yoy_pct": 100 * (denom.reindex(all_q) / denom.reindex(all_q).shift(4) - 1),
        "subs_per_death": spd,
        "n_substances": counts[counts["n"] > 0].groupby("quarter").size()
                        .reindex(all_q, fill_value=0),
    })

    onsets = pd.Series(0, index=all_q, dtype=int)
    if ALARM_PATH.exists():
        h = pd.read_csv(ALARM_PATH, parse_dates=["first_breach"])
        c = h["first_breach"].value_counts()
        onsets = onsets.add(c.reindex(all_q, fill_value=0), fill_value=0)
    tab["breach_onsets"] = onsets.astype(int)

    out_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out_dir / "regime_diagnostics.csv")

    typer.echo("Quarterly recording diagnostics "
               f"({all_q[0].year}Q{all_q[0].quarter}–"
               f"{all_q[-1].year}Q{all_q[-1].quarter})\n")
    d = tab.copy()
    d.index = [f"{q.year}Q{q.quarter}" for q in d.index]
    typer.echo(d.tail(16).to_string(float_format=lambda v: f"{v:.2f}"))

    peak_q = denom.idxmax()
    last4 = float(denom.reindex(all_q[-4:]).sum())
    prev4 = float(denom.reindex(all_q[-8:-4]).sum())
    typer.echo(
        f"\n1. Denominator. Peak {int(denom.max())} deaths in "
        f"{peak_q.year}Q{peak_q.quarter}; last 4 quarters {last4:,.0f} vs the "
        f"4 before {prev4:,.0f} ({100 * (last4 / prev4 - 1):+.0f}%). "
        "Every substance holding a flat count gains "
        f"x{prev4 / last4:.2f} share on this alone.")
    typer.echo(
        f"2. Panel. Substances named per death {spd.iloc[0]:.2f} → "
        f"{spd.iloc[-1]:.2f}; last 12 quarters mean "
        f"{spd.iloc[-12:].mean():.2f}, sd {spd.iloc[-12:].std():.3f} — "
        + ("flat, so recent flags are not a widening panel."
           if spd.iloc[-12:].std() < 0.06 else
           "still moving; treat recent flags with caution."))

    clusters = tab[tab["breach_onsets"] >= 3]
    if len(clusters):
        typer.echo(f"3. Simultaneous onsets (>=3 substances first breaching "
                   f"EB05 {threshold:g} in one quarter):")
        for q, r in clusters.iterrows():
            names = ", ".join(sorted(
                h.loc[h["first_breach"] == q, "substance"]))
            typer.echo(f"     {q.year}Q{q.quarter}  {int(r.breach_onsets)} "
                       f"substances, deaths {int(r.deaths)} "
                       f"({r.yoy_pct:+.0f}% YoY): {names}")
        typer.echo("   Several unrelated substances starting at once is more "
                   "consistent with a change in\n   the denominator or in "
                   "recording than with simultaneous independent epidemics.")
    else:
        typer.echo("3. No quarter has >=3 simultaneous first breaches.")
    typer.echo(f"\nwrote {out_dir / 'regime_diagnostics.csv'}")


@app.command("watch")
def watch(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
) -> None:
    """Check the national emerging-substance watchlist against LA County.

    Answers "is what the rest of the country is seeing showing up here", which
    the ranking cannot: a substance with zero LA mentions has no rank. Reports
    every watchlist entry including the zeros, because absence is the finding
    for most of them.

    **A zero is ambiguous and the output says so.** The lexicon recognises all
    of these names, so a zero means the narratives never mention the
    substance. It does not mean the substance is not in the LA supply — an
    assay the ME does not run produces exactly the same zero. Separating the
    two needs the ME's panel, not this data.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    res, meta = rank_substances(counts, denom, 0, recent_quarters,
                                baseline_quarters)
    m = pd.read_parquet(mentions)
    m = m[~m["is_class_term"].fillna(False)]
    key = m["rollup"].fillna(m["canonical"])

    fam = sorted(set(key[key.str.contains(NITAZENE_RE, case=False, na=False)]))
    targets = list(WATCHLIST.items()) + [
        (f"[family] nitazenes ({len(fam)} analogs)",
         "novel synthetic opioids, individually too rare to score")]

    hist = (pd.read_csv(ALARM_PATH, index_col=0, parse_dates=["first_breach"])
            if ALARM_PATH.exists() else None)

    rows = []
    for name, note in targets:
        members = fam if name.startswith("[family]") else [name]
        sub = counts[counts["substance"].isin(members)]
        total = int(sub["n"].sum())
        r = res.loc[[x for x in members if x in res.index]]
        fired = (hist is not None
                 and any(x in hist.index for x in members))
        rows.append({
            "watch": name, "note": note, "n_total": total,
            "first_seen": sub.loc[sub["n"] > 0, "quarter"].min()
                          if total else pd.NaT,
            "last_seen": sub.loc[sub["n"] > 0, "quarter"].max()
                         if total else pd.NaT,
            "n_recent": int(r["n_recent"].sum()) if len(r) else 0,
            "eb05_max": float(r["eb05"].max()) if len(r) else np.nan,
            "has_fired": fired,
        })

    tab = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out_dir / "watchlist.csv", index=False)

    q = lambda t: "—" if pd.isna(t) else f"{t.year}Q{t.quarter}"  # noqa: E731
    typer.echo(_fmt_window(meta) + "\n")
    for r in rows:
        if r["n_total"] == 0:
            typer.echo(f"  {r['watch'][:52]:54s} ABSENT — 0 mentions, "
                       f"2012-{pd.Timestamp(cutoff).year}")
        else:
            typer.echo(
                f"  {r['watch'][:52]:54s} n={r['n_total']:<5d} "
                f"{q(r['first_seen'])}→{q(r['last_seen'])}  "
                f"recent={r['n_recent']:<4d} EB05={r['eb05_max']:.2f}"
                + ("  [has fired]" if r["has_fired"] else ""))
        typer.echo(f"  {'':54s} {r['note']}")
    typer.echo("\nA zero means the narratives never name it. It cannot be "
               "distinguished from\na substance the ME's panel does not "
               "assay for — that needs the panel, not this data.")
    typer.echo(f"\nwrote {out_dir / 'watchlist.csv'}")


@app.command("alarms")
def alarms(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
) -> None:
    """Sweep every substance across every as-of quarter and record when each
    one crossed the alarm line.

    `backtest` asks "would the detector have caught these five". This asks the
    complementary question — "what has it ever fired on at all" — which is the
    only way to see the false-alarm rate and the clustering of onsets. It
    writes both the full sweep (for plotting) and the per-substance summary.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    sw = sweep_eb05(counts, denom, recent_quarters, baseline_quarters)
    seen = counts[counts["n"] > 0].groupby("substance")["quarter"].min()
    hist = alarm_history(sw, seen, max(denom.index), threshold,
                         first_q=min(denom.index))

    out_dir.mkdir(parents=True, exist_ok=True)
    sw.to_csv(out_dir / "eb05_sweep.csv", index=False)
    hist.to_csv(out_dir / ALARM_PATH.name)

    q = lambda t: "" if pd.isna(t) else f"{t.year}Q{t.quarter}"  # noqa: E731
    all_q = sorted(denom.index)
    typer.echo(f"swept {len(all_q) - recent_quarters - baseline_quarters + 1} "
               f"as-of quarters over {sw['substance'].nunique()} substances; "
               f"{len(hist)} ever exceeded EB05 {threshold}\n")
    disp = hist.copy()
    for c in ("first_seen", "first_breach", "last_breach", "peak_as_of"):
        disp[c] = disp[c].map(q)
    cols = ["first_seen", "first_breach", "lag_quarters", "last_breach",
            "quarters_in_alarm", "n_episodes", "peak_eb05", "n_at_peak",
            "still_open"]
    typer.echo(disp[cols].to_string(float_format=lambda v: f"{v:.2f}"))
    typer.echo(f"\nwrote {out_dir / ALARM_PATH.name} "
               f"and {out_dir / 'eb05_sweep.csv'}")


def sweep_eb05(
    counts: pd.DataFrame, denom: pd.Series,
    recent_quarters: int = RECENT_QUARTERS,
    baseline_quarters: int = BASELINE_QUARTERS,
) -> pd.DataFrame:
    """Re-score *every* substance at every as-of quarter.

    Same as-of discipline as `backtest`, but unrestricted to the known set:
    at each step the prior is refitted and the scores computed from data
    available then, so nothing after the as-of date leaks in. Returns a long
    frame (substance, as_of, eb05, ebgm, n_recent).

    This is the expensive command in the module — it refits the GPS prior once
    per quarter — so it is a separate CLI whose output the plots read from
    disk rather than something `plot` recomputes.
    """
    all_q = sorted(denom.index)
    total_q = recent_quarters + baseline_quarters
    frames = []
    for i in range(total_q, len(all_q) + 1):
        as_of = all_q[i - 1]
        res, _ = rank_substances(
            counts[counts["quarter"] <= as_of], denom.loc[all_q[:i]],
            lag_quarters=0, recent_quarters=recent_quarters,
            baseline_quarters=baseline_quarters)
        r = res.loc[res["n_recent"] > 0, ["n_recent", "ebgm", "eb05"]].copy()
        r["as_of"] = as_of
        frames.append(r.reset_index().rename(columns={"index": "substance"}))
    return pd.concat(frames, ignore_index=True)


def _episodes(quarters: list[pd.Timestamp]) -> list[tuple]:
    """Contiguous runs of alarm quarters -> [(start, end), ...].

    Quarters are contiguous when they are one QuarterBegin apart. A substance
    that breaches, subsides, and breaches again has two episodes, and that is
    reported rather than smoothed into one long span — an intermittent alarm
    and a sustained one mean different things.
    """
    if not quarters:
        return []
    runs, start, prev = [], quarters[0], quarters[0]
    for q in quarters[1:]:
        if q == prev + pd.offsets.QuarterBegin(1, startingMonth=1):
            prev = q
            continue
        runs.append((start, prev))
        start = prev = q
    runs.append((start, prev))
    return runs


def alarm_history(
    sweep: pd.DataFrame, first_seen: pd.Series, last_as_of: pd.Timestamp,
    threshold: float = ALARM_THRESHOLD, first_q: pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Per-substance alarm record: onset, first breach, end, peak, duration.

    The four dates a reader needs to judge a flag, in the order they happen:
    `first_seen` (the first LA death naming it), `first_breach` (when the
    detector would have fired), `last_breach` (when concern ended), and
    `peak_as_of`. `lag_quarters` between the first two is the detection delay
    — the honest measure of how much warning this buys.

    `still_open` marks substances whose last breach is the final as-of
    quarter: their alarm has not ended, it just has not ended *yet*.

    **`lag_quarters` is only meaningful for substances that arrived after the
    extract began.** LACME records start 2012Q1, so anything already
    circulating then reports `first_seen == 2012Q1` — a censoring floor, not an
    arrival. Reading morphine's 53 quarters as a detection delay would be
    nonsense: it was never detected late, it was simply always here and only
    recently grew. Those rows are flagged `first_seen_censored` and their lag
    is NaN, not a number.
    """
    rows = []
    for sub, g in sweep.groupby("substance"):
        g = g.sort_values("as_of")
        hot = g[g["eb05"] > threshold]
        if hot.empty:
            continue
        eps = _episodes(list(hot["as_of"]))
        peak = g.loc[g["eb05"].idxmax()]
        fs = first_seen.get(sub, pd.NaT)
        first_breach = hot["as_of"].min()
        censored = bool(first_q is not None and pd.notna(fs) and fs <= first_q)
        lag = ((first_breach.year - fs.year) * 4 + first_breach.quarter
               - fs.quarter) if pd.notna(fs) else np.nan
        rows.append({
            "substance": sub,
            "first_seen": fs,
            "first_seen_censored": censored,
            "first_breach": first_breach,
            "last_breach": hot["as_of"].max(),
            "lag_quarters": np.nan if censored else lag,
            "quarters_in_alarm": int(len(hot)),
            "n_episodes": len(eps),
            "longest_episode_q": max(
                (e - s).n // 3 + 1 for s, e in
                [(pd.Period(a, "Q"), pd.Period(b, "Q")) for a, b in eps]),
            "peak_eb05": float(peak["eb05"]),
            "peak_ebgm": float(peak["ebgm"]),
            "peak_as_of": peak["as_of"],
            "n_at_peak": int(peak["n_recent"]),
            "still_open": bool(hot["as_of"].max() == last_as_of),
        })
    return (pd.DataFrame(rows).set_index("substance")
            .sort_values("first_breach"))


def _episode_spans(sweep: pd.DataFrame, sub: str,
                   threshold: float) -> list[tuple]:
    g = sweep[sweep["substance"] == sub].sort_values("as_of")
    return _episodes(list(g.loc[g["eb05"] > threshold, "as_of"]))


def _quarterly_series(counts: pd.DataFrame, substance: str,
                      quarters: list) -> pd.Series:
    s = (counts[counts["substance"] == substance]
         .set_index("quarter")["n"].reindex(quarters, fill_value=0))
    return s


def _panel(ax, x, y, title: str, color: str, first_seen=None,
           trunc_from=None, annotate_last: bool = True,
           alarm_spans: list | None = None) -> None:
    """One small-multiple panel.

    Complete quarters are drawn solid; the right-truncated tail is drawn muted
    and dashed rather than dropped, so the reader can see that the apparent
    cliff at the right edge is missing toxicology and not a real collapse.

    `alarm_spans` shades the quarters in which EB05 exceeded the alarm line, so
    the detector's verdict is visible against the raw counts that produced it —
    the point being that the shading starts on the *upslope*, not at the peak.
    """
    x = list(x)
    n_complete = (sum(1 for q in x if q < trunc_from) if trunc_from is not None
                  else len(x))
    xc, yc = x[:n_complete], y.iloc[:n_complete]

    ax.plot(xc, yc, color=color, linewidth=2, solid_capstyle="round", zorder=3)
    ax.fill_between(xc, 0, yc, color=color, alpha=0.12, linewidth=0, zorder=2)
    if n_complete < len(x):
        # Bridge from the last complete point so the dashed tail connects.
        xt, yt = x[n_complete - 1:], y.iloc[n_complete - 1:]
        ax.plot(xt, yt, color=MUTED, linewidth=1.5, linestyle="--", zorder=3)
        ax.axvspan(x[n_complete], x[-1], color=MUTED, alpha=0.22, linewidth=0,
                   zorder=0)
    # Alarm band. Deliberately unlabelled in-panel: at this facet width the
    # label collided with the neighbouring panel's title. The breach quarter
    # goes in the panel title and the band is explained in the caption.
    if alarm_spans:
        for k, (s, e) in enumerate(alarm_spans):
            ax.axvspan(s, e + pd.offsets.QuarterEnd(0), color=YELLOW,
                       alpha=0.22, linewidth=0, zorder=1)
            if k == 0:
                ax.axvline(s, color=YELLOW, linewidth=1.6, zorder=2)
    # Only mark first appearance when it happens inside the plotted window —
    # for substances present since 2012 the line is just noise on the axis.
    if (first_seen is not None and pd.notna(first_seen)
            and len(x) and first_seen > x[0]):
        ax.axvline(first_seen, color=ORANGE, linewidth=1.2, linestyle=":",
                   zorder=1)
    if annotate_last and len(yc) and float(np.max(yc)) > 0:
        ax.annotate(f"{int(yc.iloc[-1])}", (xc[-1], yc.iloc[-1]),
                    textcoords="offset points", xytext=(4, 2),
                    fontsize=8, color=INK_2)
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=4)
    ax.tick_params(labelsize=7.5, colors=INK_2, length=2)
    ax.grid(axis="y", color=MUTED, alpha=0.35, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)


@app.command("plot")
def plot(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(FIG_DIR, "--out-dir"),
    lag_quarters: int = typer.Option(LAG_QUARTERS, "--lag-quarters"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    top: int = typer.Option(12, "--top", help="Facets in the rising-candidate grid"),
    min_recent: int = typer.Option(3, "--min-recent"),
    start: str = typer.Option("2015-01-01", "--start", help="Plot x-axis start"),
    heatmap_start: str = typer.Option(
        "2012-01-01", "--heatmap-start",
        help="Heatmap x-axis start; defaults to the beginning of the extract"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
) -> None:
    """Render the figures."""
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    res, meta = rank_substances(counts, denom, lag_quarters,
                                recent_quarters, baseline_quarters)
    res = _attach_novelty(res, mentions)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_q = sorted(denom.index)
    start_ts = pd.Timestamp(start)
    q_plot = [q for q in all_q if q >= start_ts]
    trunc_from = all_q[-lag_quarters] if lag_quarters else None
    window = _fmt_window(meta)

    # Alarm history, if `alarms` has been run. Read from disk rather than
    # recomputed: the sweep refits the prior 45 times and would dominate the
    # runtime of a command whose job is to draw pictures.
    sweep_path = RESULTS_DIR / "eb05_sweep.csv"
    hist, sweep = None, None
    if ALARM_PATH.exists() and sweep_path.exists():
        hist = pd.read_csv(ALARM_PATH, index_col=0,
                           parse_dates=["first_seen", "first_breach",
                                        "last_breach", "peak_as_of"])
        sweep = pd.read_csv(sweep_path, parse_dates=["as_of"])
    else:
        typer.echo("note: no alarm history on disk — run `alarms` first to get "
                   "breach annotations and the alarm-timeline figure")

    def spans(name: str) -> list:
        if sweep is None:
            return []
        return [(s, e) for s, e in _episode_spans(sweep, name, threshold)
                if e >= q_plot[0]]

    # --- Figure 1: rising candidates, small multiples -------------------
    cand = res[res["n_recent"] >= min_recent].head(top)
    ncol = 4
    nrow = int(np.ceil(len(cand) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.05 * nrow),
                             sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, (name, row) in zip(axes, cand.iterrows()):
        s = _quarterly_series(counts, name, q_plot)
        t = f"{name}\nEB05 {row.eb05:.1f}  n={int(row.n_recent)}"
        if hist is not None and name in hist.index:
            b = hist.loc[name, "first_breach"]
            t += f"  · fired {b.year}Q{b.quarter}"
        _panel(ax, q_plot, s, t, BLUE, row["first_seen"], trunc_from,
               alarm_spans=spans(name))
    for ax in axes[len(cand):]:
        ax.set_visible(False)
    fig.suptitle("Substances rising fastest in LA County overdose deaths\n"
                 f"ranked by EB05 (gamma-Poisson shrunken growth); {window}",
                 fontsize=11, color=INK, x=0.01, ha="left")
    caption = (f"Quarterly deaths naming the substance. Dotted orange = first "
               f"LACME appearance. Yellow band = quarters with EB05 > "
               f"{threshold:g}; its left edge is when the detector fired.")
    if trunc_from is not None:
        caption += (" Dashed grey tail = incomplete reporting, excluded from "
                    "the statistic — the drop there is pending toxicology, not "
                    "a real decline.")
    else:
        caption += (f" Series end at {q_plot[-1].year}Q{q_plot[-1].quarter}; "
                    "later quarters are excluded as unsettled toxicology.")
    fig.text(0.01, 0.005, caption, fontsize=7.5, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(out_dir / "rising_candidates.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Figure 2: the two novelty axes ---------------------------------
    # Chemical novelty (when NFLIS could report it) against local novelty (when
    # LA County first recorded a death). Plotting EBGM here instead was tried
    # and is uninformative: no substance in the current window has a real
    # signal, so the y-axis collapses and every point piles on the 1998-10
    # censoring floor. These two axes have genuine spread and separate the two
    # kinds of "new".
    sc = res[res["date_added"].notna() & res["first_seen"].notna()].copy()
    sc = sc[sc["n_total"] >= 5]
    fig, ax = plt.subplots(figsize=(10.5, 7))
    for censored, color, label in [
        (False, BLUE, "NFLIS date added is a real date"),
        (True, ORANGE, "left-censored — in the catalog at its 1998-10 inception"),
    ]:
        g = sc[sc["date_added_censored"].fillna(False) == censored]
        ax.scatter(g["date_added"], g["first_seen"],
                   s=14 + 9 * np.sqrt(g["n_total"]), c=color, alpha=0.70,
                   edgecolors="white", linewidths=0.8, label=label, zorder=3)

    # Ring the substances the detector has ever fired on. Novelty alone says
    # nothing about whether a substance mattered; this overlays the answer, and
    # shows that the fired-on set is scattered across both novelty axes rather
    # than concentrated among the chemically new.
    ever = set(hist.index) if hist is not None else set()
    fired = sc[sc.index.isin(ever)]
    if len(fired):
        ax.scatter(fired["date_added"], fired["first_seen"],
                   s=90 + 9 * np.sqrt(fired["n_total"]), facecolors="none",
                   edgecolors=YELLOW, linewidths=2.0, zorder=4,
                   label=f"has breached EB05 {threshold:g}")

    lo = min(sc["date_added"].min(), sc["first_seen"].min())
    hi = max(sc["date_added"].max(), sc["first_seen"].max())
    ax.plot([lo, hi], [lo, hi], color=INK_2, linewidth=1, linestyle="--",
            zorder=1)
    ax.annotate("arrived in LA as soon as the chemistry existed",
                (hi, hi), textcoords="offset points", xytext=(-8, 8),
                fontsize=8, color=INK_2, ha="right")

    # Below the diagonal: LA recorded a death before NFLIS catalogued the
    # substance. Shaded and always labelled, because it is the one region a
    # reader is likely to over-interpret — see the region caption.
    ax.fill_between([lo, hi], [lo, lo], [lo, hi], color=MUTED, alpha=0.16,
                    linewidth=0, zorder=0)
    below = sc[sc["first_seen"] < sc["date_added"]]

    # LACME records start in 2012, so every substance already circulating then
    # shares a first_seen of 2012Q1 — the y-axis is left-censored too, and that
    # floor row is an artifact of the extract window, not simultaneous arrival.
    floor = sc["first_seen"].min()
    ax.axhline(floor, color=MUTED, linewidth=1, zorder=1)
    ax.annotate("LACME record start — arrival dates at or below this line are "
                "censored, not simultaneous",
                (hi, floor), textcoords="offset points", xytext=(-8, 6),
                fontsize=8, color=INK_2, ha="right")

    # Label only substances that arrived after the censoring floor. Many share
    # the same x (the 1998-10 column), so labels are de-collided greedily in
    # display space: walk top-down and push each label below the previous one
    # whenever it would overlap, drawing a leader when it moves far.
    # Breached substances are labelled only when their arrival is *observed*.
    # Ten of the 24 sit on the 2012Q1 censoring floor, where the y-coordinate
    # is an artifact of the extract window; labelling them stacks ten names in
    # a column that the figure has already told the reader to disregard.
    labelled = sc[(sc["first_seen"] >= pd.Timestamp("2016-01-01"))
                  | sc.index.isin(below.index)
                  | (sc.index.isin(ever) & (sc["first_seen"] > floor))]

    def _label(name: str) -> str:
        if hist is None or name not in hist.index:
            return name
        b = hist.loc[name, "first_breach"]
        return f"{name}  [fired {b.year}Q{b.quarter}]"

    fig.canvas.draw()  # transforms must be current before we measure
    px_to_pt = 72.0 / fig.dpi
    placed = []
    for name, r in labelled.iterrows():
        # transData works in numbers, not Timestamps, on a date axis.
        xy = (mdates.date2num(r["date_added"]), mdates.date2num(r["first_seen"]))
        _, py = ax.transData.transform(xy)
        placed.append((name, xy, py))
    placed.sort(key=lambda t: -t[2])

    min_gap_px, last_y = 11.0 / px_to_pt, None
    for name, xy, py in placed:
        target_y = py if last_y is None else min(py, last_y - min_gap_px)
        last_y = target_y
        dy_pt = (target_y - py) * px_to_pt
        ax.annotate(
            _label(name), xy, xytext=(8, dy_pt), textcoords="offset points",
            fontsize=8, color=INK, va="center",
            arrowprops=(dict(arrowstyle="-", color=MUTED, linewidth=0.6,
                             shrinkA=0, shrinkB=1) if dy_pt < -5 else None),
        )

    ax.set_xlabel("Date added to the NFLIS catalog  —  chemical novelty",
                  fontsize=9.5, color=INK_2)
    ax.set_ylabel("First LA County overdose death  —  local novelty",
                  fontsize=9.5, color=INK_2)
    # Region caption lives in the figure margin, not inside the axes. Anchored
    # in the shaded region it sat on top of the very points it describes, and
    # every attempt to reposition it collided with a different label — the plot
    # interior has no free space at this label density.
    fig.text(
        0.005, -0.012,
        "Shaded region, below the diagonal: the LA death predates the NFLIS "
        "entry. Not early detection — NFLIS catalogues substances as forensic "
        "labs report them in seized drug\nevidence, and these are prescription "
        "antihistamines that are found at autopsy but rarely seized. "
        f"Yellow rings mark substances that have breached EB05 {threshold:g}; "
        "the quarter is the first breach.",
        fontsize=8, color=INK_2, ha="left", va="top")

    ax.set_title(
        "Two kinds of new: chemical novelty vs arrival in LA County\n"
        "points far above the diagonal are established chemistry reaching LA "
        "late — diffusion, not invention. Marker size = lifetime deaths.",
        fontsize=11, color=INK, loc="left")
    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=7,
                          markerfacecolor=c, markeredgecolor="white", label=lab)
               for c, lab in [(BLUE, "NFLIS date added is a real date"),
                              (ORANGE, "left-censored — in the catalog at its "
                                       "1998-10 inception")]]
    if len(fired):
        handles.append(plt.Line2D(
            [], [], marker="o", linestyle="", markersize=10,
            markerfacecolor="none", markeredgecolor=YELLOW, markeredgewidth=2,
            label=f"has breached EB05 {threshold:g} at some point"))
    ax.legend(handles=handles, frameon=False, fontsize=8.5, loc="lower right")
    ax.grid(color=MUTED, alpha=0.35, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "novelty_vs_growth.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Figure 3: share-of-deaths heatmap ------------------------------
    # Runs the full extract by default rather than the 2015 start the line
    # panels use: a share matrix is legible at 56 columns where a 12-facet line
    # grid is not, and the early years are where the panel-drift caveat below
    # actually shows up.
    q_heat = [q for q in all_q if q >= pd.Timestamp(heatmap_start)]
    topn = res[res["n_recent"] >= min_recent].nlargest(25, "share_recent_pct").index
    wide = (counts[counts["substance"].isin(topn)]
            .pivot_table(index="substance", columns="quarter", values="n",
                         aggfunc="sum", fill_value=0)
            .reindex(index=topn, columns=q_heat, fill_value=0))
    share = 100 * wide.div(denom.reindex(q_heat).to_numpy(), axis=1)

    # Detection-practice drift. Substances named per death climbs from 1.46 in
    # 2012 to ~1.9 and plateaus; before that plateau, a low share can mean the
    # ME was not looking rather than the substance was not there. Mark where
    # the climb ends so the extended left edge is not read as absence. The
    # plateau is the mean over the last 12 quarters and the marker is the first
    # quarter reaching 95% of it.
    spd = (counts.groupby("quarter")["n"].sum()
           .reindex(all_q, fill_value=0) / denom.reindex(all_q))
    plateau = float(spd.iloc[-12:].mean())
    reached = spd.rolling(4, min_periods=4).mean() >= 0.95 * plateau
    drift_end = next((q for q in all_q if bool(reached.get(q, False))), None)
    fig, ax = plt.subplots(figsize=(12, 8))
    # Methamphetamine and fentanyl reach ~70% of deaths; on a linear ramp they
    # saturate the scale and every candidate under ~5% — the entire population
    # of interest — renders as blank white. A power norm keeps the axis in real
    # percentage units while making the low end legible.
    im = ax.imshow(share.to_numpy(), aspect="auto", cmap=SEQ,
                   interpolation="nearest",
                   norm=PowerNorm(gamma=0.4, vmin=0, vmax=float(share.to_numpy().max())))
    labels = [f"{s}  ({int(res.loc[s, 'n_total'])})" for s in share.index]
    ax.set_yticks(range(len(share)), labels, fontsize=8.5, color=INK)
    tick = range(0, len(q_heat), 4)
    ax.set_xticks(list(tick),
                  [f"{q_heat[i].year}Q{q_heat[i].quarter}" for i in tick],
                  fontsize=8, color=INK_2)
    if trunc_from is not None and trunc_from in q_heat:
        ax.axvline(q_heat.index(trunc_from) - 0.5, color=ORANGE, linewidth=1.5)
    drift_note = ""
    if drift_end is not None and drift_end in q_heat:
        xd = q_heat.index(drift_end) - 0.5
        ax.axvline(xd, color=YELLOW, linewidth=1.8, zorder=5)
        drift_note = (
            f"Yellow line ({drift_end.year}Q{drift_end.quarter}): substances "
            f"named per death stops climbing ({spd.iloc[0]:.2f} in "
            f"{all_q[0].year} → {plateau:.2f} at plateau). Left of it the "
            "toxicology panel was still widening, so a faint cell can mean "
            "the ME was not testing rather than the substance was not there.")
    cb = fig.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("% of all overdose deaths that quarter", fontsize=9, color=INK_2)
    cb.ax.tick_params(labelsize=8, colors=INK_2)
    sub = ("orange line = start of incomplete reporting"
           if trunc_from is not None and trunc_from in q_heat
           else f"{q_heat[0].year}Q{q_heat[0].quarter}–"
                f"{q_heat[-1].year}Q{q_heat[-1].quarter}, ending at the last "
                "quarter with settled toxicology")
    ax.set_title("Substance share of LA County overdose deaths, by quarter\n"
                 "top 25 by current share (lifetime deaths in parentheses); "
                 + sub, fontsize=11, color=INK, loc="left")
    if drift_note:
        fig.text(0.005, -0.005, drift_note, fontsize=8, color=INK_2,
                 ha="left", va="top")
    fig.tight_layout()
    fig.savefig(out_dir / "share_heatmap.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Figure 4: known emergences (backtest ground truth) -------------
    present = [s for s in KNOWN_EMERGENCES if s in set(counts["substance"])]
    fig, axes = plt.subplots(1, len(present), figsize=(2.9 * len(present), 2.6),
                             sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, present):
        s = _quarterly_series(counts, name, q_plot)
        fs = res.loc[name, "first_seen"] if name in res.index else None
        sub = f"{name}\nn={int(s.sum())}"
        if hist is not None and name in hist.index:
            h = hist.loc[name]
            sub += f" · fired {h.first_breach.year}Q{h.first_breach.quarter}"
            if pd.notna(h.lag_quarters):
                sub += f" (+{int(h.lag_quarters)}q)"
        _panel(ax, q_plot, s, sub, BLUE, fs, trunc_from,
               alarm_spans=spans(name))
    fig.suptitle("Known emerging substances in LA County — the backtest set\n"
                 "counts are small (xylazine: 9 lifetime), but the EB05 sweep "
                 "still ranks all five top-4 at emergence — see `backtest`",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.text(0.01, 0.005,
             f"Yellow band = quarters with EB05 > {threshold:g}; its left edge "
             "is when the detector would have fired. Dotted orange = first LA "
             "death. (+Nq) = quarters between the two.",
             fontsize=7.5, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.86))
    fig.savefig(out_dir / "known_emergences.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    written = ["rising_candidates.png", "novelty_vs_growth.png",
               "share_heatmap.png", "known_emergences.png"]

    # --- Figure 5: alarm history ----------------------------------------
    # The three dates a reviewer asked for, per substance, on one time axis:
    # first LA death, first EB05 breach, end of the alarm. The gap between the
    # first two is the detection delay and is the point of the figure.
    if hist is not None and len(hist):
        h = hist.sort_values("first_breach", ascending=False)
        fig, ax = plt.subplots(figsize=(13, 0.36 * len(h) + 2.6))
        y = np.arange(len(h))
        t_end = max(all_q) + pd.offsets.QuarterEnd(0)

        x_left = pd.Timestamp("2011-07-01")
        any_open = bool(h["still_open"].any())
        for i, (name, r) in enumerate(h.iterrows()):
            # Presence line: from the first LA death to the end of the window.
            # A censored onset gets a left caret at the axis edge instead of a
            # circle, because the line does not start there — the records do.
            if pd.notna(r["first_seen"]):
                cens = bool(r.get("first_seen_censored", False))
                ax.hlines(i, x_left if cens else r["first_seen"], t_end,
                          color=MUTED, alpha=0.55, linewidth=1.4, zorder=2)
                ax.scatter([x_left if cens else r["first_seen"]], [i],
                           marker="<" if cens else "o", s=30 if cens else 26,
                           c=MUTED if cens else "white", edgecolors=INK_2,
                           linewidths=1.1, zorder=4)
            for s, e in _episode_spans(sweep, name, threshold):
                ax.barh(i, e + pd.offsets.QuarterEnd(0) - s, left=s,
                        height=0.52, color=ORANGE if r["still_open"] else BLUE,
                        edgecolor="white", linewidth=0.8, zorder=3)
            ax.scatter([r["peak_as_of"]], [i], marker="D", s=30, c=YELLOW,
                       edgecolors="white", linewidths=0.9, zorder=5)
            ax.annotate(f"peak {r['peak_eb05']:.1f}  (n={int(r['n_at_peak'])})",
                        (t_end, i), textcoords="offset points", xytext=(8, 0),
                        fontsize=7.5, color=INK_2, va="center",
                        annotation_clip=False)
            # Detection delay. Placed above the bar rather than to its left
            # when the gap is short, or the text lands on the onset marker.
            if pd.notna(r["lag_quarters"]):
                short = (r["first_breach"] - r["first_seen"]).days < 550
                ax.annotate(
                    f"{int(r['lag_quarters'])}q", (r["first_breach"], i),
                    textcoords="offset points",
                    xytext=((2, 9) if short else (-5, 0)),
                    fontsize=7, color=INK_2,
                    va=("bottom" if short else "center"),
                    ha=("left" if short else "right"))

        ax.set_yticks(y, h.index, fontsize=8.5, color=INK)
        ax.set_ylim(-0.8, len(h) - 0.2)
        ax.set_xlim(x_left - pd.Timedelta(days=60),
                    t_end + pd.Timedelta(days=200))
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.grid(axis="x", color=MUTED, alpha=0.30, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8, colors=INK_2, length=2)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(MUTED)

        leg = [
            plt.Line2D([], [], marker="o", linestyle="-", color=MUTED,
                       markerfacecolor="white", markeredgecolor=INK_2,
                       markersize=6,
                       label="first LA death naming it, then present"),
            plt.Line2D([], [], marker="<", linestyle="", color=MUTED,
                       markeredgecolor=INK_2, markersize=7,
                       label="already present when LACME records start (2012)"),
            plt.Rectangle((0, 0), 1, 1, facecolor=BLUE,
                          label=f"in alarm (EB05 > {threshold:g})"),
            plt.Line2D([], [], marker="D", linestyle="", color=YELLOW,
                       markeredgecolor="white", markersize=7,
                       label="peak EB05"),
        ]
        if any_open:
            leg.insert(3, plt.Rectangle((0, 0), 1, 1, facecolor=ORANGE,
                                        label="alarm still open at the cutoff"))
        ax.legend(handles=leg, frameon=False, fontsize=8.5, ncol=len(leg),
                  loc="upper left", bbox_to_anchor=(0, -0.045))
        n_swept = len(all_q) - RECENT_QUARTERS - BASELINE_QUARTERS + 1
        # Denominator is the substances the sweep actually scored (non-zero in
        # some recent window), matching what `alarms` reports — not len(res),
        # which counts rows the sweep never scored.
        n_scored = sweep["substance"].nunique()
        ax.set_title(
            "Every substance the detector has ever fired on, and when\n"
            f"{len(h)} of {n_scored} substances have exceeded EB05 "
            f"{threshold:g} across {n_swept} as-of quarters. The number at each "
            "bar is the detection delay — quarters from the first LA death to "
            "the first breach — and is shown only where arrival is observed, "
            "not censored.",
            fontsize=11, color=INK, loc="left")
        fig.tight_layout()
        fig.savefig(out_dir / "alarm_history.png", dpi=150,
                    bbox_inches="tight", facecolor="white")
        plt.close(fig)
        written.append("alarm_history.png")

    for f in written:
        typer.echo(f"wrote {out_dir / f}")


if __name__ == "__main__":
    app()
