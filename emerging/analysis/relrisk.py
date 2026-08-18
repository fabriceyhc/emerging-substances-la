"""
Kelsall-Diggle relative-risk surface, for the substances big enough to carry one
(game plan P4).

`geo.py` reports that PCP is 1.38x more tightly clustered than the overdose
background and anchors on 90013, and `spacetime.py` puts a 2.81x cylinder
around 90002. Both are *summaries*: one number and one circle. Neither is a
map, and neither says how the excess is shaped between those two points. This
estimates the surface itself.

**The estimand.** With cases (deaths naming the substance) and controls (the
other overdose deaths), the log relative risk at location x is

    r(x) = log( f(x) / g(x) )

for case density f and control density g. r(x) = 0 means the substance is
represented there exactly as it is county-wide; r(x) = log 2 means twice as
often. The denominator is again overdose death itself, so this is the same
compositional question as everything else here, just answered continuously.

**Bandwidth is the whole argument, so it is not hidden.** `geo.py` explicitly
rejected kernel methods as *detectors* because "at n=9-39 that choice would
drive the answer", and it is right that the choice is decisive: for PCP the
peak relative risk runs from 13.4x at 1 km to 2.0x at 6 km. Two things are
therefore done rather than one. The default is chosen by **likelihood
cross-validation on the case/control labels** (Kelsall and Diggle 1995, what
`sparr::LSCV.risk` implements) — which targets the *ratio* rather than either
density, and is the only selector that answers the question being asked. And
`surface` prints the whole sensitivity table next to it, so the reader sees
what the choice bought. Silverman's normal-reference rule is reported for
contrast and deliberately not used: it assumes a Gaussian cloud, and LA's
deaths are a multi-centre sprawl across 100 km.

Cases and controls share one bandwidth, also following Kelsall and Diggle: the
ratio of two equally-smoothed densities is far better behaved than the ratio of
two separately-optimised ones, and most edge bias cancels.

**Tolerance contours** come from the same random-labelling null as `geo.py` —
which C of these N deaths name the substance — evaluated pointwise. These are
*pointwise* p-values over a grid of thousands of cells, so they are not
corrected for multiple comparisons and a scattering of isolated significant
cells is expected. Read contiguous regions, not pixels; `spacetime.py` is the
place to go for a single multiplicity-adjusted verdict.

**Implementation note.** `sparr` in R is the reference tool and the game plan
called for it. This is a from-scratch equivalent, kept in Python for the same
reason `treescan.py` was: one language, no data round-trip. Densities are
computed by binning to a grid and convolving with a Gaussian, which is what
`spatstat` does internally, and which incidentally makes the duplicate
coordinates at hospitals and jails harmless to the estimator rather than fatal
to a cross-validated bandwidth.

Usage:
    emerging relrisk surface
    emerging relrisk surface --substance Cocaine --bandwidth-km 3
    emerging relrisk plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from scipy.ndimage import gaussian_filter

from emerging.ingest.extract import MENTIONS_PATH
from emerging.config import CUTOFF, SINCE
from emerging.core.spatial import MAX_CASE_FRACTION, UTM11N, load_points
from emerging.paths import ZIPS_PATH, results_dir
from emerging.viz import INK, INK_2, MUTED

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("relrisk")

CELL_M = 500.0             # grid resolution
N_PERM = 499
SENSITIVITY_KM = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0)
# Wider and finer than the reporting grid, so the CV optimum is not an
# artefact of where the candidates happen to fall.
CV_CANDIDATES_KM = (0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0)
MIN_CASES = 100            # below this a surface is a bandwidth artefact
# Do not estimate a ratio where fewer than this many deaths inform it. The
# Antelope Valley and the mountains hold almost no overdose deaths, and
# without this the log-ratio there is the ratio of two numbers near zero --
# which produced a first map running from -20 to +20 and hid the real range
# of roughly -1 to +1 entirely.
MIN_LOCAL_DEATHS = 5.0
DEFAULT_SUBSTANCE = "PCP"


def normal_reference_km(xy: np.ndarray) -> float:
    """Silverman's pooled normal-reference bandwidth, in km.

    Reported for contrast, not used as the default. It assumes an
    approximately Gaussian point cloud, and LA's overdose deaths are a
    multi-centre sprawl across 100 km; the rule returns ~4.2 km, which smooths
    the Harbor, Skid Row and the Antelope Valley into one blur.
    """
    n = len(xy)
    sd = np.sqrt(0.5 * (xy[:, 0].var(ddof=1) + xy[:, 1].var(ddof=1)))
    return float(sd * n ** (-1.0 / 6.0) / 1000.0)


def lcv_bandwidth_km(xy: np.ndarray, cases: np.ndarray,
                     candidates=CV_CANDIDATES_KM,
                     chunk: int = 400) -> tuple[float, pd.DataFrame]:
    """Kelsall-Diggle likelihood cross-validation for the log-ratio.

    The quantity that matters is the local case proportion
    p(x) = f(x) / (f(x) + g(x)), so the bandwidth is chosen to maximise the
    leave-one-out Bernoulli log-likelihood of the labels

        sum_cases log p_-i(x_i) + sum_controls log(1 - p_-i(x_i))

    rather than to minimise the error of either density on its own. This is
    the selector `sparr::LSCV.risk` implements, and it targets the ratio,
    which is the thing being estimated. Leave-one-out is exact: the point's
    own kernel contribution is subtracted at its own location.

    Evaluated directly at the data points rather than on the grid, because
    LOO at a binned location is not well defined.
    """
    n = len(xy)
    y = cases.astype(float)
    rows = []
    for hk in candidates:
        h = hk * 1000.0
        ll = 0.0
        for s0 in range(0, n, chunk):
            s1 = min(s0 + chunk, n)
            d2 = ((xy[s0:s1, None, :] - xy[None, :, :]) ** 2).sum(-1)
            w = np.exp(-0.5 * d2 / h ** 2)
            # Exact leave-one-out: drop the self term from both sums.
            np.fill_diagonal(w[:, s0:s1], 0.0)
            f = w @ y
            g = w.sum(1) - f
            pi = np.clip(f / (f + g + 1e-300), 1e-12, 1 - 1e-12)
            yi = y[s0:s1]
            ll += float(np.sum(yi * np.log(pi) + (1 - yi) * np.log(1 - pi)))
        rows.append({"bandwidth_km": hk, "cv_loglik": ll})
    tab = pd.DataFrame(rows)
    return float(tab.loc[tab["cv_loglik"].idxmax(), "bandwidth_km"]), tab


class Grid:
    """Raster over the county, with a mask for everything outside it."""

    def __init__(self, cell_m: float = CELL_M, pad_m: float = 3000.0):
        import geopandas as gpd
        from shapely.geometry import box

        zips = gpd.read_file(ZIPS_PATH).to_crs(UTM11N)
        # Catalina and San Clemente sit ~70 km offshore. Union everything
        # first, then keep only the pieces comparable in size to the mainland
        # -- an area test on the *merged* geometry, which a per-zip filter
        # cannot do because individual island zips are not unusually small.
        merged = zips.geometry.union_all()
        pieces = list(getattr(merged, "geoms", [merged]))
        biggest = max(p.area for p in pieces)
        from shapely.ops import unary_union
        self.region = unary_union([p for p in pieces if p.area >= 0.05 * biggest])
        minx, miny, maxx, maxy = self.region.bounds
        self.x = np.arange(minx - pad_m, maxx + pad_m, cell_m)
        self.y = np.arange(miny - pad_m, maxy + pad_m, cell_m)
        self.cell_m = cell_m
        self.shape = (len(self.y), len(self.x))

        gx, gy = np.meshgrid(self.x + cell_m / 2, self.y + cell_m / 2)
        pts = gpd.GeoSeries(gpd.points_from_xy(gx.ravel(), gy.ravel()),
                            crs=UTM11N)
        self.inside = pts.within(self.region).to_numpy().reshape(self.shape)
        self.extent = (self.x[0], self.x[-1], self.y[0], self.y[-1])

    def cell_index(self, xy: np.ndarray) -> np.ndarray:
        ix = np.clip(np.searchsorted(self.x, xy[:, 0]) - 1, 0, len(self.x) - 1)
        iy = np.clip(np.searchsorted(self.y, xy[:, 1]) - 1, 0, len(self.y) - 1)
        return iy * len(self.x) + ix

    def bin(self, idx: np.ndarray) -> np.ndarray:
        return np.bincount(idx, minlength=self.shape[0] * self.shape[1]) \
                 .reshape(self.shape).astype(float)


def log_rr(case_grid: np.ndarray, all_grid: np.ndarray, sigma_cells: float,
           eps: float = 1e-12,
           min_local: float = MIN_LOCAL_DEATHS) -> np.ndarray:
    """Kelsall-Diggle log relative risk, cases against all deaths.

    Both surfaces are smoothed with the *same* kernel, so the normalising
    constants and most of the edge bias divide out; what is left is the log of
    a locally-estimated case proportion against the global one.
    """
    # scipy's gaussian_filter returns the *input* dtype, so an integer count
    # grid comes back truncated to integers -- a smoothed density of 0.3
    # becomes 0 and the whole surface collapses. Coerce rather than trust the
    # caller; this cost an afternoon once.
    f = gaussian_filter(np.asarray(case_grid, dtype=float), sigma_cells,
                        mode="constant")
    g = gaussian_filter(np.asarray(all_grid, dtype=float), sigma_cells,
                        mode="constant")
    p_global = case_grid.sum() / all_grid.sum()
    # Pseudo-count, centred so that it pulls toward "no difference" rather
    # than toward zero risk. Half a death, expressed in the same density units
    # as f and g. Without it a cell holding 40 deaths and no cases reports
    # log(eps/g) -- about -27 with eps = 1e-12 -- and a handful of such cells
    # set the colour scale for the whole map, hiding the real +/-1 range.
    # Algebra check: when f/g equals p_global the added terms cancel and
    # r is exactly 0, whatever the local density.
    a = 0.5 / (2.0 * np.pi * sigma_cells ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.log((f + a) / (g + a / p_global)) - np.log(p_global)
    # Effective number of deaths within one bandwidth of the cell. A Gaussian
    # filter on counts returns deaths per cell, so multiplying by the kernel's
    # area in cells recovers the local sample size the estimate rests on.
    local_n = g * 2.0 * np.pi * sigma_cells ** 2
    return np.where(local_n >= min_local, r, np.nan)


def tolerance(case_idx: np.ndarray, all_idx: np.ndarray, grid: Grid,
              sigma_cells: float, n_perm: int, rng: np.random.Generator,
              observed: np.ndarray) -> np.ndarray:
    """Pointwise upper-tail p-values under random labelling."""
    all_grid = grid.bin(all_idx)
    n_case = len(case_idx)
    ge = np.zeros(grid.shape, dtype=np.int32)
    for _ in range(n_perm):
        pick = rng.choice(len(all_idx), size=n_case, replace=False)
        sim = log_rr(grid.bin(all_idx[pick]), all_grid, sigma_cells)
        ge += np.nan_to_num(sim >= observed, nan=False)
    return (1.0 + ge) / (n_perm + 1.0)


def _load(substance: str, since: str, cutoff: str):
    bg, sets = load_points(MENTIONS_PATH, since=since, cutoff=cutoff)
    if substance not in sets:
        raise typer.BadParameter(f"unknown substance: {substance}")
    cases = bg["CaseNumber"].isin(sets[substance]).to_numpy()
    frac = cases.sum() / len(bg)
    if frac > MAX_CASE_FRACTION:
        raise typer.BadParameter(
            f"{substance} is {frac:.0%} of all deaths — above the "
            f"{MAX_CASE_FRACTION:.0%} guard, so the controls are mostly the "
            f"substance itself and the surface would be flat by construction")
    return bg, cases


@app.command("surface")
def surface(
    substance: str = typer.Option(DEFAULT_SUBSTANCE, "--substance"),
    since: str = typer.Option(SINCE, "--since"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    bandwidth_km: float = typer.Option(0.0, "--bandwidth-km",
                                       help="0 = pooled normal-reference rule"),
    n_perm: int = typer.Option(N_PERM, "--n-perm"),
    cell_m: float = typer.Option(CELL_M, "--cell-m"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    seed: int = typer.Option(20260811, "--seed"),
) -> None:
    """Estimate the log relative-risk surface and its tolerance contours."""
    bg, cases = _load(substance, since, cutoff)
    n_case = int(cases.sum())
    if n_case < MIN_CASES:
        raise typer.BadParameter(
            f"{substance} has {n_case} deaths; below {MIN_CASES} the surface "
            f"is a picture of the bandwidth, not of the data")

    xy = bg[["x", "y"]].to_numpy()
    grid = Grid(cell_m)
    all_idx = grid.cell_index(xy)
    case_idx = all_idx[cases]
    all_grid, case_grid = grid.bin(all_idx), grid.bin(case_idx)

    h_ref = normal_reference_km(xy)
    if bandwidth_km:
        h, cvtab = bandwidth_km, None
    else:
        h, cvtab = lcv_bandwidth_km(xy, cases)
    typer.echo(f"{substance}: {n_case:,} cases against {len(bg) - n_case:,} "
               f"other overdose deaths, {since}..{cutoff}")
    typer.echo(f"grid {grid.shape[1]}x{grid.shape[0]} at {cell_m:.0f} m, "
               f"{int(grid.inside.sum()):,} cells inside the county")
    typer.echo(f"bandwidth {h:.2f} km"
               + ("  (Kelsall-Diggle likelihood CV)" if not bandwidth_km
                  else "  (set by hand)")
               + f"; Silverman's normal-reference rule would say "
                 f"{h_ref:.2f} km")
    if cvtab is not None:
        typer.echo("  CV log-likelihood: " + "  ".join(
            f"{r.bandwidth_km:g}km {r.cv_loglik:,.1f}"
            for r in cvtab.itertuples()))

    # --- sensitivity across bandwidths, before committing to one -----------
    rows = []
    for hk in sorted(set(SENSITIVITY_KM) | {h}):
        r = log_rr(case_grid, all_grid, hk * 1000.0 / cell_m)
        inside = grid.inside & np.isfinite(r)
        rows.append({"bandwidth_km": hk,
                     "max_log_rr": float(np.nanmax(np.where(inside, r, np.nan))),
                     "max_rr": float(np.exp(np.nanmax(np.where(inside, r, np.nan)))),
                     "pct_area_rr_above_2": 100 * float(
                         np.nansum(inside & (r > np.log(2))) / inside.sum()),
                     "chosen": bool(abs(hk - h) < 1e-9)})
    sens = pd.DataFrame(rows)
    typer.echo("\nbandwidth sensitivity (before tolerance testing)")
    typer.echo(sens.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # --- the chosen surface, with tolerance --------------------------------
    sigma = h * 1000.0 / cell_m
    r = log_rr(case_grid, all_grid, sigma)
    rng = np.random.default_rng(seed)
    p = tolerance(case_idx, all_idx, grid, sigma, n_perm, rng, r)

    inside = grid.inside & np.isfinite(r)
    sig = inside & (p <= 0.05)
    typer.echo(f"\n{100 * sig.sum() / inside.sum():.1f}% of the county is "
               f"above the county-wide rate at pointwise p <= 0.05 "
               f"({n_perm} random-labelling replicates)")
    typer.echo(f"peak relative risk {np.exp(np.nanmax(np.where(inside, r, np.nan))):.2f}x")

    # --- per zip, so it joins to everything else ---------------------------
    import geopandas as gpd
    zips = gpd.read_file(ZIPS_PATH).to_crs(UTM11N)
    zips["zip"] = (pd.to_numeric(zips["ZIPCODE"], errors="coerce")
                     .astype("Int64").astype(str))
    gx, gy = np.meshgrid(grid.x + cell_m / 2, grid.y + cell_m / 2)
    cells = gpd.GeoDataFrame(
        {"r": r.ravel(), "p": p.ravel(), "inside": inside.ravel()},
        geometry=gpd.points_from_xy(gx.ravel(), gy.ravel()), crs=UTM11N)
    cells = cells[cells["inside"]]
    j = gpd.sjoin(cells, zips[["zip", "geometry"]], how="inner",
                  predicate="within")
    by_zip = (j.groupby("zip")
                .agg(mean_log_rr=("r", "mean"), max_log_rr=("r", "max"),
                     min_p=("p", "min"), n_cells=("r", "size"))
                .sort_values("mean_log_rr", ascending=False))
    by_zip["mean_rr"] = np.exp(by_zip["mean_log_rr"])

    out_dir.mkdir(parents=True, exist_ok=True)
    tag = substance.lower().replace(" ", "_")
    sens.to_csv(out_dir / f"{tag}_bandwidth.csv", index=False)
    by_zip.to_csv(out_dir / f"{tag}_zip.csv")
    np.savez_compressed(out_dir / f"{tag}_surface.npz",
                        log_rr=r.astype(np.float32), p=p.astype(np.float32),
                        inside=inside, extent=np.array(grid.extent),
                        bandwidth_km=h, substance=substance)

    typer.echo("\nhighest zips by mean relative risk")
    typer.echo(by_zip.head(10)[["mean_rr", "max_log_rr", "min_p", "n_cells"]]
               .to_string(float_format=lambda v: f"{v:.3f}"))
    typer.echo(f"\nwrote {out_dir / f'relrisk_{tag}_zip.csv'} and "
               f"{out_dir / f'relrisk_{tag}_surface.npz'}")


@app.command("plot")
def plot(
    substance: str = typer.Option(DEFAULT_SUBSTANCE, "--substance"),
    results_dir: Path = typer.Option(RESULTS_DIR, "--results-dir"),
    fig_dir: Path = typer.Option(RESULTS_DIR, "--fig-dir"),
) -> None:
    """Map the surface, its tolerance contour, and the bandwidth sensitivity."""
    import geopandas as gpd
    from matplotlib.colors import TwoSlopeNorm

    tag = substance.lower().replace(" ", "_")
    d = np.load(results_dir / f"{tag}_surface.npz", allow_pickle=False)
    r, p, inside = d["log_rr"], d["p"], d["inside"]
    extent = tuple(d["extent"])
    h = float(d["bandwidth_km"])
    sens = pd.read_csv(results_dir / f"{tag}_bandwidth.csv")

    zips = gpd.read_file(ZIPS_PATH).to_crs(UTM11N)
    shown = np.where(inside, r, np.nan)
    lim = float(np.nanpercentile(np.abs(shown), 99.5))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.8, 5.8),
                                   gridspec_kw={"width_ratios": [1.35, 1]})
    zips.boundary.plot(ax=ax1, color="white", linewidth=0.35, zorder=3)
    im = ax1.imshow(shown, origin="lower", extent=extent, cmap="RdBu_r",
                    norm=TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim),
                    zorder=2, interpolation="bilinear")
    ax1.contour(np.where(inside, p, np.nan), levels=[0.05], colors=[INK],
                linewidths=1.3, extent=extent, origin="lower", zorder=4)
    ax1.set_xlim(extent[0], extent[1])
    ax1.set_ylim(extent[2], extent[3])
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_axis_off()
    ax1.set_title(f"{substance}: log relative risk against other overdose "
                  f"deaths\nbandwidth {h:.2f} km; outline = pointwise p ≤ 0.05",
                  fontsize=11.5, color=INK, loc="left", pad=8)
    cb = fig.colorbar(im, ax=ax1, fraction=0.036, pad=0.02)
    cb.set_label("log relative risk", fontsize=9, color=INK_2)
    cb.ax.tick_params(labelsize=8, colors=INK_2)

    ax2.plot(sens["bandwidth_km"], sens["max_rr"], color=INK_2, lw=2.2,
             marker="o", ms=5)
    ch = sens[sens["chosen"]]
    ax2.scatter(ch["bandwidth_km"], ch["max_rr"], s=110, facecolor="none",
                edgecolor=INK, linewidth=1.8, zorder=5, label="chosen")
    ax2.set_xlabel("bandwidth (km)", fontsize=9.5, color=INK_2)
    ax2.set_ylabel("peak relative risk", fontsize=9.5, color=INK_2)
    ax2.set_title("The peak is a function of the bandwidth", fontsize=11.5,
                  color=INK, loc="left", pad=8)
    ax2.legend(frameon=False, fontsize=9)
    ax2.grid(color="#eceae5", lw=0.9)
    ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    ax2.tick_params(colors=INK_2, labelsize=9)

    fig.text(0.008, 0.015,
             "Cases are deaths naming the substance; controls are the other "
             "overdose deaths, so 0 means 'represented here exactly as it is "
             "county-wide'. Both densities use one bandwidth, so normalising "
             "constants and most edge bias divide out. Contours are POINTWISE "
             "p-values over thousands of cells and are not corrected for "
             "multiplicity — read contiguous regions, not pixels. The right "
             "panel is the reason this is a map and not a number.",
             fontsize=7.6, color=MUTED, va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = fig_dir / f"{tag}.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    typer.echo(f"wrote {out}")


if __name__ == "__main__":
    app()
