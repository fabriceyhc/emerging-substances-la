"""
Space-time scan: *where and when* is a substance over-represented (game plan P2).

`geo.py` answers "is this substance more concentrated than overdose death
itself", county-wide and over one fixed window. It cannot say *where*, and it
cannot say *when*. This scans circles of zipcodes crossed with recent time
windows and returns the single most unlikely cylinder, with a p-value adjusted
for having searched all of them.

**Which SaTScan method this is, and why not the two the game plan named.**
P2 asked for "spatial variation in temporal trends", and mentioned the
space-time permutation model in the same breath. Those are two different
things, and this is a third:

  - *Spatial variation in temporal trends* fits a Poisson trend inside versus
    outside each circle. It needs a population denominator per zip per
    quarter, which is exactly what this project does not have — the reason
    every method here is compositional.
  - *Space-time permutation* (Kulldorff 2005) needs no denominator, which is
    why it was attractive. But it conditions on one case series' own margins,
    so for substance X it asks "do X's deaths cluster in space-time relative
    to X's own spatial and temporal distribution". With 9 xylazine deaths
    there is nothing there, and for a common substance it mostly rediscovers
    where overdose deaths are.
  - **What this uses instead: the Bernoulli space-time scan, case-control.**
    Cases are deaths naming the substance; controls are the *other* overdose
    deaths in the same zip and quarter. The denominator is therefore overdose
    death itself, which is the same convention as EB05, `treescan.py` and
    `geo.py`, and it keeps the whole pipeline answering one question:
    composition of the supply, not incidence.

The null is random labelling — which C of these N deaths name the substance —
identical to `geo.py`'s, and defended there. The p-value comes from the
distribution of the *maximum* LLR over every cylinder, so it is adjusted for
the search in the same way `treescan.py`'s is.

**"Growing" versus "concentrated".** Windows are trailing (they end at the
as-of quarter), so a cluster means "this substance's share of overdose deaths
is elevated in these zips over the last k quarters, relative to its share
across the whole county and the whole study period". That is the operational
form of *growing here*, and it is what distinguishes this from `geo.py`, which
would report the same answer whether the concentration was new or a decade old.

**Read the power warning before reading the results.** The signal-ceiling work
and P1 both landed on the same wall: three of the five known emergences have
fewer than 10 deaths. Splitting those across 280 zips and 12 quarters cannot
create information. This is expected to return negatives for everything except
the few substances with hundreds of deaths, and it is worth running mainly
because that expectation should be tested rather than assumed.

Usage:
    emerging spacetime scan
    emerging spacetime scan --substances PCP,Lidocaine --n-perm 9999
    emerging spacetime plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from scipy.special import xlogy

from emerging.ingest.extract import MENTIONS_PATH
from emerging.config import CUTOFF, SINCE
from emerging.core.spatial import MAX_CASE_FRACTION, UTM11N, load_points
from emerging.paths import ZIPS_PATH, results_dir
from emerging.core.tree import build
from emerging.viz import AQUA, INK, MUTED, ORANGE

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("spacetime")

# Largest circle, as a share of all overdose deaths in the study window. The
# SaTScan convention is 50%; 25% is used here for the reason the TreeScan user
# guide gives for temporal windows — a "cluster" covering most of the map is
# more honestly described as a deficit in the remainder.
MAX_SPATIAL_FRACTION = 0.25
MAX_WINDOW = 6           # trailing quarters, matching treescan.MAX_WINDOW
STUDY_QUARTERS = 12
N_PERM = 999
# Minimum cases for a substance to be scanned at all. Higher than `geo.py`'s 8
# because a cylinder splits those cases across space *and* time.
MIN_CASES = 20

SCAN_PATH = RESULTS_DIR / "clusters.csv"


def _bernoulli_llr(c: np.ndarray, n: np.ndarray, C: float,
                   N: float) -> np.ndarray:
    """Case-control scan LLR for a cylinder holding c cases of n deaths.

    `xlogy` gives 0 log 0 = 0, which is the right value for an empty cell and
    the reason this is not written with a bare `np.log`.
    """
    out_c, out_n = C - c, N - n
    # Empty cylinders and full-county cylinders are masked by `valid3` in the
    # caller, but the division happens first, so silence it here rather than
    # letting a RuntimeWarning fire 160k times per replicate.
    with np.errstate(divide="ignore", invalid="ignore"):
        inside = xlogy(c, c / n) + xlogy(n - c, (n - c) / n)
        outside = xlogy(out_c, out_c / out_n) + xlogy(
            out_n - out_c, (out_n - out_c) / out_n)
        base = xlogy(C, C / N) + xlogy(N - C, (N - C) / N)
        llr = inside + outside - base
    # High-rate only: the cylinder's case share must exceed the county's.
    # Cross-multiplied to avoid dividing by an empty complement.
    elevated = c * out_n > out_c * n
    return np.where(elevated & np.isfinite(llr), llr, 0.0)


def _poisson_llr(c: np.ndarray, e: np.ndarray, C: float) -> np.ndarray:
    """Space-time permutation LLR, conditioned on both of the case margins."""
    with np.errstate(divide="ignore", invalid="ignore"):
        llr = xlogy(c, c / e) + xlogy(C - c, (C - c) / (C - e))
    return np.where((c > e) & np.isfinite(llr), llr, 0.0)


def _trailing(x: np.ndarray, max_window: int) -> np.ndarray:
    """(..., n_q) -> (..., max_window) trailing-window sums.

    Applied *before* the spatial accumulation, not after. The two commute, and
    doing time first keeps the big intermediate array at (centers, k, windows)
    instead of (centers, k, quarters).
    """
    rev = np.cumsum(x[..., ::-1], axis=-1)
    return rev[..., :max_window]


class Geometry:
    """Zip centroids, neighbour ordering, and the fixed denominator cylinders."""

    def __init__(self, bg: pd.DataFrame, quarters: list[pd.Timestamp],
                 max_fraction: float = MAX_SPATIAL_FRACTION,
                 max_window: int = MAX_WINDOW):
        self.quarters = quarters
        self.max_window = max_window

        # Centre each zip on the mean location of its own overdose deaths
        # rather than the polygon centroid: it is the centre of the at-risk
        # set, and it needs no join against the boundary file.
        cen = bg.groupby("zip")[["x", "y"]].mean()
        self.zips = list(cen.index)
        xy = cen.to_numpy()

        d = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=2)
        self.order = np.argsort(d, axis=1, kind="stable")
        self.dist_km = np.take_along_axis(d, self.order, axis=1) / 1000.0

        tot = self._grid(bg, bg.index)
        self.N = float(tot.sum())
        cum_deaths = np.cumsum(tot.sum(axis=1)[self.order], axis=1)
        # k_max[i] = how many zips the circle around i may reach.
        self.k_len = np.maximum((cum_deaths <= max_fraction * self.N).sum(1), 1)
        self.k_max = int(self.k_len.max())
        # Ragged circles in a rectangular array: mask the rows past each
        # centre's own limit so they can never win the maximum.
        self.valid = (np.arange(self.k_max)[None, :] < self.k_len[:, None])

        self.n_cyl = self.accumulate(tot)          # totals: fixed, computed once
        self.valid3 = (self.valid[:, :, None]
                       & (self.n_cyl > 0) & (self.n_cyl < self.N))

    def _grid(self, bg: pd.DataFrame, idx) -> np.ndarray:
        """(n_zips, n_quarters) death counts for a subset of rows."""
        zpos = {z: i for i, z in enumerate(self.zips)}
        qpos = {q: i for i, q in enumerate(self.quarters)}
        g = np.zeros((len(self.zips), len(self.quarters)))
        sub = bg.loc[idx]
        zi = sub["zip"].map(zpos).to_numpy()
        qi = sub["quarter"].map(qpos).to_numpy()
        ok = ~(pd.isna(zi) | pd.isna(qi))
        np.add.at(g, (zi[ok].astype(int), qi[ok].astype(int)), 1.0)
        return g

    def accumulate(self, grid: np.ndarray) -> np.ndarray:
        """(n_zips, n_q) counts -> (n_centres, k_max, n_windows) cylinder sums."""
        win = _trailing(grid, self.max_window)              # (Z, W)
        gathered = win[self.order[:, :self.k_max]]          # (Z, k_max, W)
        return np.cumsum(gathered, axis=1)

    def accumulate_batch(self, grids: np.ndarray) -> np.ndarray:
        """(b, n_zips, n_q) -> (b, n_centres, k_max, n_windows).

        The same operation as `accumulate` for a stack of replicates. Doing it
        as one gather rather than a Python loop over replicates is most of the
        Monte Carlo runtime.
        """
        win = _trailing(grids, self.max_window)              # (b, Z, W)
        gathered = win[:, self.order[:, :self.k_max]]        # (b, Z, k_max, W)
        return np.cumsum(gathered, axis=2)

    def accumulate_spatial(self, vec: np.ndarray) -> np.ndarray:
        """(n_zips,) -> (n_centres, k_max) running totals out from each centre."""
        return np.cumsum(vec[self.order[:, :self.k_max]], axis=1)


def scan_substance(geo: Geometry, case_grid: np.ndarray, C: float,
                   n_perm: int, rng: np.random.Generator,
                   bg: pd.DataFrame, all_idx: np.ndarray,
                   model: str = "bernoulli", batch: int = 40) -> dict:
    """Most likely cylinder for one substance, with a permutation p-value.

    `model="bernoulli"` — cases against the other overdose deaths in the same
    cylinder. Answers *where is this over-represented*. A concentration that
    has been there for a decade wins as easily as one that appeared last year,
    because the reference is the substance's share pooled over the whole study
    period.

    `model="interaction"` — Kulldorff's space-time permutation, conditioned on
    the substance's **own** spatial margin and its **own** temporal margin, so
    expected = (cases in these zips, all quarters) x (cases in this window, all
    zips) / total. Being a busy zip cannot win, and neither can being a busy
    quarter; only the *combination* can. That is "growing here" rather than
    "concentrated here", and it is the distinction P2 exists to draw.
    """
    n_q = len(geo.quarters)
    n_cells = len(geo.zips) * n_q
    cell = bg["_zi"].to_numpy() * n_q + bg["_qi"].to_numpy()
    c_cyl = geo.accumulate(case_grid)

    if model == "bernoulli":
        exp_cyl = geo.n_cyl * C / geo.N
        valid = geo.valid3
        llr = np.where(valid, _bernoulli_llr(c_cyl, geo.n_cyl, C, geo.N), 0.0)
    elif model == "interaction":
        # Both margins are fixed by the null, so the expectation is computed
        # once and reused for every replicate.
        spatial = geo.accumulate_spatial(case_grid.sum(axis=1))     # (Z, k)
        temporal = _trailing(case_grid.sum(axis=0), geo.max_window)  # (W,)
        exp_cyl = spatial[:, :, None] * temporal[None, None, :] / C
        valid = geo.valid[:, :, None] & (exp_cyl > 0) & (exp_cyl < C)
        llr = np.where(valid, _poisson_llr(c_cyl, exp_cyl, C), 0.0)
    else:
        raise ValueError("model must be 'bernoulli' or 'interaction'")

    ci, ki, wi = np.unravel_index(int(llr.argmax()), llr.shape)
    obs_llr = float(llr[ci, ki, wi])

    null = np.empty(n_perm)
    done = 0
    case_cells = None
    if model == "interaction":
        case_zi = np.repeat(np.arange(len(geo.zips)),
                            case_grid.sum(axis=1).astype(int))
        case_qi = np.repeat(np.tile(np.arange(n_q), len(geo.zips)),
                            case_grid.ravel().astype(int))
        case_cells = (case_zi, case_qi)

    while done < n_perm:
        b = min(batch, n_perm - done)
        grids = np.empty((b, len(geo.zips), n_q))
        for j in range(b):
            if model == "bernoulli":
                # Random labelling: which C of the N deaths name the
                # substance. Same null as `geo.py`, defended there.
                pick = rng.choice(all_idx, size=int(C), replace=False)
                grids[j] = np.bincount(cell[pick], minlength=n_cells).reshape(
                    len(geo.zips), n_q)
            else:
                # Space-time permutation: shuffle *when* each case happened
                # while leaving *where* alone, which holds both margins fixed.
                zi, qi = case_cells
                g = np.zeros((len(geo.zips), n_q))
                np.add.at(g, (zi, rng.permutation(qi)), 1.0)
                grids[j] = g
        sim = geo.accumulate_batch(grids)
        if model == "bernoulli":
            s = _bernoulli_llr(sim, geo.n_cyl[None], C, geo.N)
        else:
            s = _poisson_llr(sim, exp_cyl[None], C)
        null[done:done + b] = np.where(valid[None], s, 0.0).reshape(
            b, -1).max(axis=1)
        done += b

    p = (1.0 + float((null >= obs_llr).sum())) / (n_perm + 1.0)
    members = geo.order[ci, :ki + 1]
    n_in = float(geo.n_cyl[ci, ki, wi])
    c_in = float(c_cyl[ci, ki, wi])
    e_in = float(exp_cyl[ci, ki, wi])
    return {
        "model": model,
        "n_cases": int(C),
        "llr": obs_llr,
        "p_value": p,
        "recurrence_interval": 1.0 / p,
        "centre_zip": geo.zips[ci],
        "n_zips": int(ki + 1),
        "radius_km": float(geo.dist_km[ci, ki]),
        "window_q": int(wi + 1),
        "window_from": geo.quarters[-(wi + 1)],
        "observed": int(c_in),
        "expected": e_in,
        "deaths_in_cylinder": int(n_in),
        "rate_ratio": (c_in / e_in) if e_in else np.nan,
        # Every member, not a preview: the plot shades this list, and a
        # truncated one drew a 12-zip blob for a 69-zip cluster.
        "zips": "|".join(geo.zips[m] for m in members),
    }


def _case_sets(sets: pd.Series, nodes: bool) -> dict[str, set]:
    """Substance -> case numbers, optionally with tree branches folded in."""
    out = {s: set(v) for s, v in sets.items()}
    if nodes:
        for nd in build(MENTIONS_PATH)[0]:
            if nd.level in ("substance", "root"):
                continue
            union: set = set()
            for leaf in nd.leaves:
                union |= out.get(leaf, set())
            if union:
                out[nd.name] = union
    return out


@app.command("scan")
def scan(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    since: str = typer.Option(SINCE, "--since"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    substances: str = typer.Option("", "--substances",
                                   help="comma-separated; default is every "
                                        "eligible substance and tree branch"),
    include_nodes: bool = typer.Option(True, "--nodes/--no-nodes",
                                       help="also scan tree branches"),
    min_cases: int = typer.Option(MIN_CASES, "--min-cases"),
    n_perm: int = typer.Option(N_PERM, "--n-perm"),
    max_window: int = typer.Option(MAX_WINDOW, "--max-window"),
    max_fraction: float = typer.Option(MAX_SPATIAL_FRACTION, "--max-fraction"),
    model: str = typer.Option("both", "--model",
                              help="bernoulli (over-represented where) | "
                                   "interaction (growing where) | both"),
    seed: int = typer.Option(20260811, "--seed"),
) -> None:
    """Scan every eligible substance and write the most likely cluster of each."""
    bg, sets = load_points(mentions, since=since, cutoff=cutoff)
    bg = bg.reset_index(drop=True)

    co_q = pd.read_parquet(mentions, columns=["CaseNumber", "quarter"])
    co_q = co_q.drop_duplicates("CaseNumber").set_index("CaseNumber")["quarter"]
    bg["quarter"] = bg["CaseNumber"].map(co_q)
    bg = bg.dropna(subset=["quarter", "zip"])
    bg = bg[bg["zip"] != "<NA>"]
    quarters = sorted(bg["quarter"].unique())[-STUDY_QUARTERS:]
    bg = bg[bg["quarter"].isin(quarters)].reset_index(drop=True)

    geo = Geometry(bg, quarters, max_fraction, max_window)
    bg["_zi"] = bg["zip"].map({z: i for i, z in enumerate(geo.zips)})
    bg["_qi"] = bg["quarter"].map({q: i for i, q in enumerate(quarters)})

    typer.echo(f"{len(bg):,} OD deaths, {len(geo.zips)} zips, "
               f"{len(quarters)} quarters "
               f"({quarters[0].year}Q{quarters[0].quarter}–"
               f"{quarters[-1].year}Q{quarters[-1].quarter})")
    typer.echo(f"circles up to {max_fraction:.0%} of deaths "
               f"(<= {geo.k_max} zips), trailing windows 1–{max_window}q, "
               f"{n_perm} random-labelling replicates")

    pool = _case_sets(sets, include_nodes)
    live = set(bg["CaseNumber"])
    pool = {k: (v & live) for k, v in pool.items()}
    if substances:
        want = {s.strip() for s in substances.split(",") if s.strip()}
        missing = want - set(pool)
        if missing:
            raise typer.BadParameter(f"unknown: {', '.join(sorted(missing))}")
        pool = {k: v for k, v in pool.items() if k in want}

    idx_of = pd.Series(bg.index.to_numpy(), index=bg["CaseNumber"])
    all_idx = bg.index.to_numpy()
    rng = np.random.default_rng(seed)

    eligible, skipped = {}, []
    for name, cases in pool.items():
        n = len(cases)
        frac = n / len(bg)
        if n < min_cases:
            skipped.append((name, n, "too few cases"))
        elif frac > MAX_CASE_FRACTION:
            # Same guard as geo.py: as n/N approaches 1 the null draw becomes
            # the case set itself and the test compares a substance to itself.
            skipped.append((name, n, f"{frac:.0%} of all deaths — out of scope"))
        else:
            eligible[name] = cases
    typer.echo(f"scanning {len(eligible)} of {len(pool)} "
               f"(min {min_cases} cases, <= {MAX_CASE_FRACTION:.0%} of deaths)\n")

    models = ["bernoulli", "interaction"] if model == "both" else [model]
    rows = []
    for name, cases in sorted(eligible.items()):
        rowidx = idx_of.loc[sorted(cases)].to_numpy()
        grid = np.zeros((len(geo.zips), len(quarters)))
        np.add.at(grid, (bg.loc[rowidx, "_zi"].to_numpy().astype(int),
                         bg.loc[rowidx, "_qi"].to_numpy().astype(int)), 1.0)
        for mdl in models:
            res = scan_substance(geo, grid, float(len(cases)), n_perm, rng, bg,
                                 all_idx, model=mdl)
            res["substance"] = name
            rows.append(res)

    # Substance breaks ties on llr, and the sort is stable. Without it two llr
    # values that differ only in the last bit -- which happens across numpy
    # versions -- reorder the table, so a `git diff` of the committed CSV shows
    # several changed rows where nothing changed. Same defect, same fix, as in
    # `polysubstance.profile_table`.
    out = pd.DataFrame(rows).sort_values(
        ["model", "llr", "substance"], ascending=[True, False, True],
        kind="stable")
    # Each p-value is already adjusted for every cylinder searched *within* one
    # substance. It is not adjusted for having scanned 46 substances, and
    # reporting it as though it were is precisely the mistake the tree scan was
    # built to stop. Bonferroni across the scanned set is conservative here,
    # because the set includes nested nodes (Dissociatives contains PCP), but
    # it is the statement that cannot be accused of over-claiming.
    out["p_across_substances"] = (
        out["p_value"] * out.groupby("model")["p_value"].transform("size")
    ).clip(upper=1.0)
    cols = ["substance", "model", "n_cases", "centre_zip", "n_zips",
            "radius_km", "window_q", "window_from", "observed", "expected",
            "rate_ratio", "deaths_in_cylinder", "llr", "p_value",
            "p_across_substances", "recurrence_interval", "zips"]
    out = out[cols]
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / SCAN_PATH.name, index=False)

    label = {"bernoulli": "OVER-REPRESENTED WHERE (vs other OD deaths)",
             "interaction": "GROWING WHERE (space-time interaction, "
                            "conditioned on the substance's own margins)"}
    for mdl in models:
        g = out[out["model"] == mdl]
        sig = g[g["p_value"] <= 0.05]
        surv = g[g["p_across_substances"] <= 0.05]
        typer.echo(f"\n{label[mdl]}")
        typer.echo(f"  {len(sig)} of {len(g)} at p <= 0.05 (adjusted for every "
                   f"cylinder searched, within a substance)")
        typer.echo(f"  {len(surv)} still at 0.05 after correcting for the "
                   f"{len(g)} substances scanned")
        disp = g.head(12).drop(columns=["zips", "window_from", "model",
                                        "recurrence_interval"])
        typer.echo(disp.to_string(index=False,
                                  float_format=lambda v: f"{v:.3f}"))
    if skipped:
        typer.echo(f"\nnot scanned: {len(skipped)}")
        for name, n, why in sorted(skipped, key=lambda t: -t[1])[:6]:
            typer.echo(f"  {name:34s} n={n:5d}  {why}")
    typer.echo(f"\nwrote {out_dir / SCAN_PATH.name}")


@app.command("plot")
def plot(
    results_dir: Path = typer.Option(RESULTS_DIR, "--results-dir"),
    fig_dir: Path = typer.Option(RESULTS_DIR, "--fig-dir"),
    model: str = typer.Option("bernoulli", "--model"),
    top: int = typer.Option(4, "--top"),
) -> None:
    """Map the most likely cylinder of the strongest few substances."""
    import geopandas as gpd

    res = pd.read_csv(results_dir / SCAN_PATH.name)
    res = res[res["model"] == model].nlargest(top, "llr").reset_index(drop=True)
    zips = gpd.read_file(ZIPS_PATH).to_crs(UTM11N)
    zips["zip"] = (pd.to_numeric(zips["ZIPCODE"], errors="coerce")
                     .astype("Int64").astype(str))
    # Catalina and San Clemente sit ~70 km offshore and stretch every panel to
    # the point where the mainland is a smudge. Frame on the mainland mass via
    # centroid quantiles rather than `total_bounds` -- the same fix `geo.py`
    # needed for the identical reason.
    cx, cy = zips.geometry.centroid.x, zips.geometry.centroid.y
    x0, x1 = cx.quantile([0.01, 0.99])
    y0, y1 = cy.quantile([0.01, 0.99])
    padx, pady = 0.06 * (x1 - x0), 0.06 * (y1 - y0)

    ncol = min(len(res), 2)
    nrow = int(np.ceil(len(res) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.6 * ncol, 4.9 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, r) in zip(axes, res.iterrows()):
        members = set(str(r["zips"]).split("|"))
        zips.plot(ax=ax, color="#f4f2ee", edgecolor="white", linewidth=0.4)
        hit = zips[zips["zip"].isin(members)]
        survives = r["p_across_substances"] <= 0.05
        colour = ORANGE if survives else (AQUA if r["p_value"] <= 0.05 else MUTED)
        hit.plot(ax=ax, color=colour, edgecolor="white", linewidth=0.4,
                 alpha=0.85)
        centre = zips[zips["zip"] == str(r["centre_zip"])]
        if not centre.empty:
            centre.plot(ax=ax, facecolor="none", edgecolor=INK, linewidth=0.9)
        pv = ("p = %.3f" % r["p_value"]) if r["p_value"] > 0.001 else "p < 0.001"
        adj = ("survives correction for all 46 substances"
               if survives else "does not survive correction for 46 substances")
        ax.set_title(
            f"{r['substance']}  (n={int(r['n_cases'])})\n"
            f"{int(r['observed'])} vs {r['expected']:.1f} expected · "
            f"{r['rate_ratio']:.2f}x · {int(r['window_q'])}q · {pv}\n"
            f"{adj}",
            fontsize=10.5, color=INK, loc="left", pad=8)
        ax.set_axis_off()
        ax.set_xlim(x0 - padx, x1 + padx)
        ax.set_ylim(y0 - pady, y1 + pady)
        ax.set_aspect("equal", adjustable="box")

    for ax in axes[len(res):]:
        ax.set_axis_off()

    fig.text(0.008, 0.012,
             "Bernoulli space-time scan, case-control: cases are deaths naming "
             "the substance, controls are the other overdose deaths in the "
             "same zip and quarter. Shaded = the most likely cylinder, outlined "
             "= its centre zip. p-values come from the random-labelling "
             "distribution of the maximum LLR, so each is already adjusted for "
             "every circle and window searched within that substance; orange "
             "additionally survives Bonferroni across the 46 substances "
             "scanned. These are static concentrations -- the space-time "
             "interaction scan, which asks where a substance is *growing*, "
             "returns nothing that survives the same correction.",
             fontsize=7.6, color=MUTED, va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = fig_dir / "clusters.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
