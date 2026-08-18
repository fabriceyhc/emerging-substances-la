"""
Spatial variation in temporal trends: where is mortality rising fastest?
(game plan P7 — the method P2 originally named)

Every other scan in this repository asks where a *level* is unusual. This asks
where a *slope* is unusual, and the difference is the whole point. A
level-based scan on a population denominator would mostly rediscover Skid Row,
South LA and the Antelope Valley, because that is where overdose mortality is
high — driven by who lives and uses where, not by anything changing. SVTT
(Kulldorff 2006) allows every circle its own baseline level and tests only
whether its *trend* differs from the rest of the county, so "this area has
always been bad" cannot produce a signal and "this area is getting worse
faster than everywhere else" can.

**The model.** Cases in zip-year cells, with resident population as the
offset. For a candidate circle Z:

    inside   c_t ~ Poisson( P_t^in  · exp(a_in  + b_in ·t) )
    outside  c_t ~ Poisson( P_t^out · exp(a_out + b_out·t) )

    H0: b_in = b_out  (levels a_in, a_out still free — that is what makes
                       this a test of trend rather than of level)
    H1: b_in > b_out  (one-sided: rising faster inside)

`exp(b) - 1` is the annual relative change, so b = 0.07 is "about 7% a year
faster than the rest of the county".

**The fit.** The MLE of a log-linear Poisson with an offset reduces to a
one-dimensional monotone root find. Differentiating the log-likelihood in `a`
gives exp(a) = C / sum(P_t exp(b t)); substituting into the score for `b`
leaves

    sum_t t P_t exp(b t) / sum_t P_t exp(b t)  =  sum_t t c_t / C

— the population-weighted mean time under trend `b` must equal the
case-weighted mean time. The left side is a mean under an exponential tilt, so
its derivative in `b` is a variance and is strictly positive: the root is
unique and Newton converges from anywhere. That is what makes this affordable,
because it vectorises across all ~27,000 circles at once instead of running a
GLM per circle.

**The null.** Condition on the total number of cases and scatter them over
zip-year cells with probability proportional to `population x exp(b_global t)`,
where `b_global` is the county-wide trend. That is exactly H0 — one common
trend, arbitrary spatial variation in level — so the simulation preserves the
county's overall rise and only destroys spatial variation *in the slope*.

**Read these before the results.**

  - Annual, not quarterly. ACS population is annual, so quarterly figures
    would be interpolation rather than information, and a trend does not need
    the resolution.
  - Population is residential. 90013 has 14,891 residents and a much larger
    daytime population, so its per-resident rate is place-based intensity,
    not individual risk. A *trend* is more robust to this than a level, since
    the bias has to be changing to matter — but Skid Row is exactly where that
    assumption is weakest.
  - ACS 5-year estimates are pre-smoothed and end in 2024; 2025 carries
    forward. A zip whose population genuinely moved in 2025 will show that as
    a spurious rate change.

Usage:
    emerging svtt scan
    emerging svtt scan --substance PCP
    emerging svtt plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

from emerging.ingest.census import load_population
from emerging.ingest.extract import MENTIONS_PATH
from emerging.config import CUTOFF
from emerging.core.spatial import UTM11N, load_points
from emerging.paths import ZIPS_PATH, results_dir
from emerging.viz import INK, INK_2, MUTED, ORANGE

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("svtt")

# LA's overdose rate is a hump, not a trend: 5.3 per 100k in 2013, peaking at
# 29.5 in 2022, down to 21.7 in 2025. A log-linear fit across the whole of that
# is misspecified in the worst way -- it returns +14.5%/yr, which describes
# neither the rise nor the fall, and "rising fastest" then conflates "rose
# earliest during the ramp" with "still rising now". So the default window is
# the **decline phase**, where the county trend is monotone and the question
# is the one worth asking: where is mortality falling slowest, or still rising?
# Pass --first-year 2015 for the whole-epidemic version, and read the annual
# rates the scan prints before believing either.
FIRST_YEAR = 2022
MAX_SPATIAL_FRACTION = 0.25
N_PERM = 999
NEWTON_ITERS = 24
MIN_CASES_IN = 25          # a slope through fewer cases, over 4 years, is noise

SCAN_PATH = RESULTS_DIR / "clusters.csv"


def _fit_beta(C: np.ndarray, TC: np.ndarray, P: np.ndarray,
              t: np.ndarray, iters: int = NEWTON_ITERS) -> np.ndarray:
    """Vectorised MLE of the log-linear trend, one per row of P.

    Solves  weighted_mean_t(P_t e^{bt}) = TC/C  by Newton. The derivative of
    the left side is the variance of t under those weights, hence positive,
    hence the iteration is globally convergent.
    """
    beta = np.zeros(len(C))
    ok = C > 0
    target = np.divide(TC, C, out=np.zeros_like(TC), where=ok)
    for _ in range(iters):
        w = P * np.exp(beta[:, None] * t[None, :])
        A = w.sum(1)
        good = ok & (A > 0)
        m = np.divide((w * t).sum(1), A, out=np.zeros_like(A), where=good)
        v = np.divide((w * t * t).sum(1), A, out=np.zeros_like(A),
                      where=good) - m * m
        step = np.divide(target - m, v, out=np.zeros_like(v),
                         where=good & (v > 1e-12))
        beta = beta + np.clip(step, -2.0, 2.0)
    return np.where(ok, beta, 0.0)


def _profile_ll(C: np.ndarray, TC: np.ndarray, P: np.ndarray, t: np.ndarray,
                beta: np.ndarray) -> np.ndarray:
    """Poisson log-likelihood at the profile MLE, dropping terms free of (a,b).

    At the MLE exp(a) sum(P e^{bt}) = C exactly, so the fitted-count term
    collapses to -C and the whole thing is three cheap pieces.
    """
    A = (P * np.exp(beta[:, None] * t[None, :])).sum(1)
    ok = (C > 0) & (A > 0)
    a = np.log(np.divide(C, A, out=np.ones_like(A), where=ok))
    return np.where(ok, C * a + beta * TC - C, 0.0)


def _fit_common_beta(C1, TC1, P1, C2, TC2, P2, t,
                     iters: int = NEWTON_ITERS) -> np.ndarray:
    """One trend shared by inside and outside, levels still free (H0)."""
    beta = np.zeros(len(C1))
    target = TC1 + TC2
    for _ in range(iters):
        g = np.zeros(len(C1))
        gp = np.zeros(len(C1))
        for C, P in ((C1, P1), (C2, P2)):
            w = P * np.exp(beta[:, None] * t[None, :])
            A = w.sum(1)
            good = A > 0
            m = np.divide((w * t).sum(1), A, out=np.zeros_like(A), where=good)
            v = np.divide((w * t * t).sum(1), A, out=np.zeros_like(A),
                          where=good) - m * m
            g += C * m
            gp += C * v
        step = np.divide(target - g, gp, out=np.zeros_like(gp),
                         where=gp > 1e-12)
        beta = beta + np.clip(step, -2.0, 2.0)
    return beta


def svtt_llr(c_in: np.ndarray, p_in: np.ndarray, c_all: np.ndarray,
             p_all: np.ndarray, t: np.ndarray) -> tuple:
    """LLR for a difference in trend, one row per circle. One-sided."""
    c_out = c_all[None, :] - c_in
    p_out = p_all[None, :] - p_in

    C_in, TC_in = c_in.sum(1), (c_in * t).sum(1)
    C_out, TC_out = c_out.sum(1), (c_out * t).sum(1)

    b_in = _fit_beta(C_in, TC_in, p_in, t)
    b_out = _fit_beta(C_out, TC_out, p_out, t)
    b_0 = _fit_common_beta(C_in, TC_in, p_in, C_out, TC_out, p_out, t)

    ll1 = (_profile_ll(C_in, TC_in, p_in, t, b_in)
           + _profile_ll(C_out, TC_out, p_out, t, b_out))
    ll0 = (_profile_ll(C_in, TC_in, p_in, t, b_0)
           + _profile_ll(C_out, TC_out, p_out, t, b_0))

    llr = np.maximum(ll1 - ll0, 0.0)
    # One-sided and count-floored, applied identically to every replicate.
    keep = (b_in > b_out) & (C_in >= MIN_CASES_IN) & (C_out >= MIN_CASES_IN)
    return np.where(keep, llr, 0.0), b_in, b_out


class Circles:
    """Zip centroids ordered by distance, with a size cap on the radius."""

    def __init__(self, cen: pd.DataFrame, weight: np.ndarray,
                 max_fraction: float = MAX_SPATIAL_FRACTION):
        self.zips = list(cen.index)
        xy = cen.to_numpy()
        d = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=2)
        self.order = np.argsort(d, axis=1, kind="stable")
        self.dist_km = np.take_along_axis(d, self.order, axis=1) / 1000.0
        cum = np.cumsum(weight[self.order], axis=1)
        self.k_len = np.maximum((cum <= max_fraction * weight.sum()).sum(1), 1)
        self.k_max = int(self.k_len.max())
        self.valid = np.arange(self.k_max)[None, :] < self.k_len[:, None]

    def accumulate(self, grid: np.ndarray) -> np.ndarray:
        """(n_zips, T) -> (n_circles, T) for every (centre, radius) circle."""
        g = np.cumsum(grid[self.order[:, :self.k_max]], axis=1)
        return g.reshape(-1, grid.shape[1])

    def labels(self) -> tuple[np.ndarray, np.ndarray]:
        ci = np.repeat(np.arange(len(self.zips)), self.k_max)
        ki = np.tile(np.arange(self.k_max), len(self.zips))
        return ci, ki


def build_grids(substance: str | None, mentions: Path, cutoff: str,
                first_year: int) -> tuple:
    """(cases[zip, year], population[zip, year], zips, years, centroids)."""
    bg, sets = load_points(mentions, since=f"{first_year}-01-01", cutoff=cutoff)
    bg = bg.reset_index(drop=True)
    bg["year"] = pd.to_datetime(
        pd.read_parquet(mentions, columns=["CaseNumber", "death_date"])
        .drop_duplicates("CaseNumber").set_index("CaseNumber")["death_date"]
        .reindex(bg["CaseNumber"]).to_numpy()).year
    bg = bg[bg["zip"].notna() & (bg["zip"] != "<NA>")]
    bg = bg.dropna(subset=["year"])

    if substance:
        if substance not in sets:
            raise typer.BadParameter(f"unknown substance: {substance}")
        bg = bg[bg["CaseNumber"].isin(sets[substance])]

    pop = load_population()
    years = sorted(int(y) for y in bg["year"].unique())
    zips = sorted(set(bg["zip"]) & set(pop["zipcode"]))
    bg = bg[bg["zip"].isin(zips)]

    zpos = {z: i for i, z in enumerate(zips)}
    ypos = {y: i for i, y in enumerate(years)}
    cases = np.zeros((len(zips), len(years)))
    np.add.at(cases, (bg["zip"].map(zpos).to_numpy().astype(int),
                      bg["year"].map(ypos).to_numpy().astype(int)), 1.0)

    wide = pop.pivot_table(index="zipcode", columns="year", values="population",
                           aggfunc="first").reindex(zips)
    P = np.zeros((len(zips), len(years)))
    have = sorted(wide.columns)
    for j, y in enumerate(years):
        # Years past the ACS range carry the last vintage forward.
        src = y if y in wide.columns else max(have)
        P[:, j] = wide[src].to_numpy()
    P = np.nan_to_num(P, nan=0.0)

    cen = bg.groupby("zip")[["x", "y"]].mean().reindex(zips)
    return cases, P, zips, years, cen


@app.command("scan")
def scan(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    substance: str = typer.Option("", "--substance",
                                  help="default: all overdose deaths"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    first_year: int = typer.Option(FIRST_YEAR, "--first-year"),
    n_perm: int = typer.Option(N_PERM, "--n-perm"),
    max_fraction: float = typer.Option(MAX_SPATIAL_FRACTION, "--max-fraction"),
    seed: int = typer.Option(20260811, "--seed"),
    top: int = typer.Option(10, "--top"),
) -> None:
    """Scan for the area whose mortality trend most exceeds the county's."""
    cases, P, zips, years, cen = build_grids(substance or None, mentions,
                                             cutoff, first_year)
    t = np.arange(len(years), dtype=float)
    t = t - t.mean()                       # centred, for conditioning
    circ = Circles(cen, P[:, -1], max_fraction)

    c_all, p_all = cases.sum(0), P.sum(0)
    label = substance or "all overdose deaths"
    typer.echo(f"{label}: {int(cases.sum()):,} deaths, {len(zips)} zips, "
               f"{years[0]}–{years[-1]}, {p_all[-1]:,.0f} residents")
    typer.echo(f"{len(zips) * circ.k_max:,} circles (<= {circ.k_max} zips, "
               f"{max_fraction:.0%} of population), {n_perm} replicates")

    # County-wide trend: both the headline number and the null's H0.
    b_g = _fit_beta(np.array([c_all.sum()]), np.array([(c_all * t).sum()]),
                    p_all[None, :], t)[0]
    typer.echo(f"county-wide trend: {100 * (np.exp(b_g) - 1):+.2f}% per year")
    # Print the series the trend is fitted to. A log-linear model is only
    # honest over a monotone stretch, and LA's epidemic is a hump -- this is
    # the line that stops the whole-epidemic fit being quoted by accident.
    typer.echo("  annual rate per 100k: " + "  ".join(
        f"{y} {1e5 * c / p:.1f}" for y, c, p in zip(years, c_all, p_all)))

    c_cyl = circ.accumulate(cases)
    p_cyl = circ.accumulate(P)
    valid = circ.valid.ravel()
    llr, b_in, b_out = svtt_llr(c_cyl, p_cyl, c_all, p_all, t)
    llr = np.where(valid, llr, 0.0)
    obs = float(llr.max())

    # Null: one common trend, spatial variation in level only.
    q = P * np.exp(b_g * t)[None, :]
    q = (q / q.sum()).ravel()
    total = int(cases.sum())
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for b in range(n_perm):
        sim = rng.multinomial(total, q).reshape(P.shape).astype(float)
        s_llr, _, _ = svtt_llr(circ.accumulate(sim), p_cyl, sim.sum(0), p_all, t)
        null[b] = np.where(valid, s_llr, 0.0).max()
    p_value = (1.0 + float((null >= obs).sum())) / (n_perm + 1.0)

    ci, ki = circ.labels()
    order = np.argsort(-llr)[:top]
    rows = []
    for r in order:
        if llr[r] <= 0:
            continue
        members = circ.order[ci[r], :ki[r] + 1]
        rows.append({
            "substance": label,
            "centre_zip": circ.zips[ci[r]],
            "n_zips": int(ki[r] + 1),
            "radius_km": float(circ.dist_km[ci[r], ki[r]]),
            "cases_inside": int(c_cyl[r].sum()),
            "pop_inside": float(p_cyl[r][-1]),
            "trend_inside_pct": 100 * (np.exp(b_in[r]) - 1),
            "trend_outside_pct": 100 * (np.exp(b_out[r]) - 1),
            "llr": float(llr[r]),
            # Only the maximum has a calibrated p-value; the rest are
            # secondary clusters, reported with the same p for reference.
            "p_value": p_value if r == order[0] else np.nan,
            "zips": "|".join(circ.zips[m] for m in members),
        })
    out = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = ("clusters.csv" if not substance
            else f"clusters_{substance.lower().replace(' ', '_')}.csv")
    out.to_csv(out_dir / name, index=False)

    typer.echo(f"\nmost likely cluster: p = {p_value:.4f} "
               f"(adjusted for all {len(zips) * circ.k_max:,} circles)\n")
    disp = out.drop(columns=["zips", "substance"]).head(top)
    typer.echo(disp.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    typer.echo(f"\nwrote {out_dir / name}")


@app.command("profile")
def profile(
    results_dir: Path = typer.Option(RESULTS_DIR, "--results-dir"),
    name: str = typer.Option("clusters.csv", "--name"),
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    first_year: int = typer.Option(FIRST_YEAR, "--first-year"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
) -> None:
    """Characterise the most likely cluster: what is in it, and is it real?

    A trend cluster is a claim about a place, and the first three questions
    anyone will ask are what substances are involved, whether it is an
    artifact of where people are pronounced dead, and who is dying. This
    answers all three against the rest of the county, so the cluster can be
    reported rather than just located.
    """
    from emerging.ingest.cohort import filter_alcohol_only
    from emerging.paths import LACME_PATH

    res = pd.read_csv(results_dir / name)
    if res.empty:
        raise typer.Exit("no cluster to profile")
    r = res.iloc[0]
    members = set(str(r["zips"]).split("|"))

    bg, sets = load_points(mentions, since=f"{first_year}-01-01")
    meta = filter_alcohol_only(pd.read_csv(LACME_PATH, low_memory=False),
                               verbose=False)
    meta["dt"] = pd.to_datetime(meta["DeathDate"], format="mixed",
                                errors="coerce")
    bg = bg.join(meta[["CaseNumber", "Age", "Gender", "Race", "dt"]]
                 .drop_duplicates("CaseNumber").set_index("CaseNumber"),
                 on="CaseNumber")
    bg["year"] = bg["dt"].dt.year
    inside = bg["zip"].isin(members)
    c, o = bg[inside], bg[~inside]
    typer.echo(f"cluster: {int(r['n_zips'])} zips around {r['centre_zip']} "
               f"({', '.join(sorted(members))})")
    typer.echo(f"{len(c)} deaths inside, {len(o)} outside, "
               f"{first_year}–{bg['year'].max():.0f}\n")

    typer.echo("WHICH ZIPS ACTUALLY MOVED")
    t = c.pivot_table(index="zip", columns="year", values="CaseNumber",
                      aggfunc="count", fill_value=0)
    yrs = sorted(t.columns)
    t["change"] = t[yrs[-2]] + t[yrs[-1]] - t[yrs[0]] - t[yrs[1]]
    typer.echo(t.sort_values("change", ascending=False).to_string())

    typer.echo("\nIS THE RISE A BURST OR A LEVEL SHIFT?")
    byq = c.groupby([c["year"], c["dt"].dt.quarter]).size()
    typer.echo("  " + "  ".join(f"{y}Q{q} {n}" for (y, q), n in byq.items()))

    typer.echo("\nARTIFACT CHECKS")
    typer.echo(f"  at a shared (facility) address: cluster "
               f"{100 * c['at_facility'].mean():.1f}%  vs county "
               f"{100 * o['at_facility'].mean():.1f}%")
    for y in yrs:
        cy = c[c["year"] == y]
        typer.echo(f"    {int(y)}: {len(cy):3d} deaths at "
                   f"{cy[['x', 'y']].drop_duplicates().shape[0]:3d} distinct "
                   f"addresses")

    typer.echo("\nSUBSTANCES (share of deaths naming it, cluster vs county)")
    rows = []
    for sub, cases in sets.items():
        n = int(c["CaseNumber"].isin(cases).sum())
        if n < 8:
            continue
        ci = c["CaseNumber"].isin(cases).mean()
        oi = o["CaseNumber"].isin(cases).mean()
        rows.append({"substance": sub, "n": n, "cluster_pct": 100 * ci,
                     "county_pct": 100 * oi,
                     "ratio": ci / oi if oi else np.nan})
    subs = pd.DataFrame(rows).sort_values("ratio", ascending=False)
    typer.echo(subs.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    spread = subs["ratio"].max() - subs["ratio"].min()
    typer.echo(f"  ratios span {subs['ratio'].min():.2f}–"
               f"{subs['ratio'].max():.2f}"
               + ("  → no substance is driving this; it is more of the same "
                  "supply." if spread < 1.0 else
                  "  → one substance stands out; investigate it."))

    typer.echo("\nWHO")
    for col in ("Gender", "Race"):
        a = c[col].value_counts(normalize=True).mul(100)
        b = o[col].value_counts(normalize=True).mul(100)
        typer.echo(f"  {col}: " + ", ".join(
            f"{k} {a.get(k, 0):.0f}% vs {b.get(k, 0):.0f}%"
            for k in a.head(4).index))
    typer.echo(f"  median age: cluster {c['Age'].median():.0f}, "
               f"county {o['Age'].median():.0f}")

    out_dir.mkdir(parents=True, exist_ok=True)
    subs.to_csv(out_dir / "cluster_substances.csv", index=False)
    typer.echo(f"\nwrote {out_dir / 'svtt_cluster_substances.csv'}")


@app.command("plot")
def plot(
    results_dir: Path = typer.Option(RESULTS_DIR, "--results-dir"),
    fig_dir: Path = typer.Option(RESULTS_DIR, "--fig-dir"),
    name: str = typer.Option("clusters.csv", "--name"),
) -> None:
    """Map the most likely trend cluster and plot its rate against the county."""
    import geopandas as gpd

    res = pd.read_csv(results_dir / name)
    if res.empty:
        typer.echo("no cluster to plot")
        raise typer.Exit()
    r = res.iloc[0]
    members = set(str(r["zips"]).split("|"))

    zips = gpd.read_file(ZIPS_PATH).to_crs(UTM11N)
    zips["zip"] = (pd.to_numeric(zips["ZIPCODE"], errors="coerce")
                     .astype("Int64").astype(str))
    cx, cy = zips.geometry.centroid.x, zips.geometry.centroid.y
    x0, x1 = cx.quantile([0.01, 0.99])
    y0, y1 = cy.quantile([0.01, 0.99])

    cases, P, zl, years, _ = build_grids(None, MENTIONS_PATH, CUTOFF,
                                         FIRST_YEAR)
    inside = np.array([z in members for z in zl])
    rate_in = 1e5 * cases[inside].sum(0) / P[inside].sum(0)
    rate_out = 1e5 * cases[~inside].sum(0) / P[~inside].sum(0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.4, 5.6),
                                   gridspec_kw={"width_ratios": [1, 1.15]})
    zips.plot(ax=ax1, color="#f4f2ee", edgecolor="white", linewidth=0.4)
    zips[zips["zip"].isin(members)].plot(ax=ax1, color=ORANGE,
                                         edgecolor="white", linewidth=0.4)
    ax1.set_xlim(x0, x1)
    ax1.set_ylim(y0, y1)
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_axis_off()
    ax1.set_title(f"{int(r['n_zips'])} zips around {r['centre_zip']}",
                  fontsize=12, color=INK, loc="left", pad=8)

    ax2.plot(years, rate_in, color=ORANGE, lw=2.4, marker="o", ms=5,
             label=f"cluster  ({r['trend_inside_pct']:+.1f}%/yr)")
    ax2.plot(years, rate_out, color=INK_2, lw=2.0, marker="s", ms=4,
             label=f"rest of county  ({r['trend_outside_pct']:+.1f}%/yr)")
    ax2.set_ylabel("overdose deaths per 100,000 residents", fontsize=9.5,
                   color=INK_2)
    ax2.set_xticks(list(years))
    ax2.set_xticklabels([str(y) for y in years])
    ax2.legend(frameon=False, fontsize=9.5)
    ax2.grid(axis="y", color="#eceae5", lw=0.9)
    ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    ax2.tick_params(colors=INK_2, labelsize=9)
    pv = r["p_value"]
    ax2.set_title(f"Trend differs from the rest of the county "
                  f"(p = {pv:.3f})" if pd.notna(pv) else "Trend",
                  fontsize=12, color=INK, loc="left", pad=8)

    fig.text(0.008, 0.015,
             "Spatial variation in temporal trends: each circle is allowed its "
             "own baseline rate and is tested only on its slope, so an area "
             "that has always been high cannot signal. Denominator is ACS "
             "resident population, which for downtown zips understates the "
             "population actually present. p-value is from the maximum over "
             "every circle, under a null of one county-wide trend. Note that "
             "the cluster's rise is concentrated in a single year -- a log-"
             "linear slope is a summary of this shape, not a description of "
             "it.",
             fontsize=7.6, color=MUTED, va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = fig_dir / "trends.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
