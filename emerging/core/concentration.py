"""
Spatial concentration of a substance's *recent* deaths, as one standardized
number per (substance, as-of quarter) — the `z_s` that `docs/GPS_V2_DESIGN.md`
#3 needs and did not have.

**Why not reuse `geo.py` or `spacetime.py` directly.** The design note is
explicit that the lifetime `geo.py` localization ratio is the wrong input:
PCP's Skid Row clustering is real and a decade old, and letting stable
historical geography lower the detection bar for an unrelated event is the
failure mode the whole §3 exercise is meant to avoid. `spacetime.py` scans the
right (trailing) windows but costs 999 permutations per substance and refuses
to run below `MIN_CASES = 20`, which excludes every substance this project
actually cares about — xylazine has 9 lifetime mentions. Feeding a 45-quarter
as-of sweep from it is neither affordable nor defined where it matters.

**What this computes instead.** A Cuzick–Edwards-style case-control clustering
count on the deaths in the trailing window, with *exact* null moments in place
of a permutation p-value:

    S_s = #{ordered pairs of s-deaths that are spatial neighbours}

where "neighbour" is the symmetrized k-nearest-neighbour graph over *all*
overdose deaths in the window, k a fixed fraction of the window's deaths
(`NEIGHBOUR_FRACTION`). Under the same random-labelling null every other
spatial analysis here uses — which C of these N deaths name the substance is
arbitrary — S_s has closed-form mean and variance (`_moments`), so

    z_s = (S_s - E[S_s]) / sd[S_s]

costs one sparse mat-vec per substance instead of a Monte Carlo run. That is
what makes a per-quarter sweep over every substance affordable at all.

**Why k-nearest-neighbour and not a fixed radius or the zipcode cells.** Both
alternatives were built and measured first:

- *Zipcode cells* (neighbour == same zip) is the cheapest and matches
  `spacetime.py`'s geography, but the background co-location probability is
  Q = 0.009, so a 6-death substance expects 0.28 co-located pairs and scores
  S = 0 about three quarters of the time. The statistic has essentially no
  resolution in exactly the low-n regime the project lives in.
- *Fixed radius* fixes that but fixes the wrong thing: at any single radius,
  Skid Row deaths have enormous degree and Antelope Valley deaths have none,
  so the statistic mostly measures "is this substance downtown". The null
  absorbs that (it is still a valid test), but the power goes to density, not
  to concentration.

k-NN at a fixed *fraction* of the window adapts to local density and is
scale-free the same way `spacetime.MAX_SPATIAL_FRACTION` is — "the nearest 5%
of the county's overdose deaths" means the same thing in a quarter with 300
deaths and a quarter with 800.

**The low-n floor is not optional.** At C = 3 the statistic can only take the
values 0, 2, 4, 6, so a single coincidental neighbour pair moves z_s by four
standard deviations — on the 2025 window, diphenhydramine (C = 3, S = 6)
scores z = 6.7, which is discreteness, not geography. `MIN_CASES` zeroes the
score below a count where the null distribution has enough support for a
normal standardization to mean anything. See `docs/findings/concentration.md`
for the calibration that sets it.

The consequence is stated plainly rather than buried: **substances below
`MIN_CASES` get no spatial adjustment at all**, and three of the five known
LA emergences are below it in the quarter they emerge. This component cannot
help them. It is built and measured anyway because that was worth testing
rather than assuming — the same reason `spacetime.py` runs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import norm

from emerging.core.spatial import MAX_CASE_FRACTION, UTM11N
from emerging.ingest.extract import MENTIONS_PATH, load_cohort

# Neighbourhood size, as a fraction of the overdose deaths in the trailing
# window. 5% follows `spacetime.MAX_SPATIAL_FRACTION`'s reasoning at a smaller
# scale: a "concentration" spanning a quarter of the county is better described
# as a deficit elsewhere, and the useful end of that range is much tighter. At
# 1% the background co-location probability is 0.013 and low-n substances are
# back to scoring 0; at 10% the neighbourhood is a third of the way across the
# county and stops being a cluster.
NEIGHBOUR_FRACTION = 0.05

# Below this many geocoded deaths in the window, z_s is set to 0 (no spatial
# evidence either way) rather than reported. See the module docstring.
MIN_CASES = 10

# Above this share of the window's deaths, the case-control contrast
# degenerates — the null draw becomes the case set itself. Same constant and
# same reasoning as every other spatial analysis here; see `core/spatial.py`.
# Methamphetamine and fentanyl are withheld by it, as they are there.
MAX_FRACTION = MAX_CASE_FRACTION


def _falling(x: np.ndarray | float, k: int) -> np.ndarray:
    """x(x-1)...(x-k+1), the falling factorial, elementwise."""
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    for j in range(k):
        out = out * (x - j)
    return out


def _moments(a2: float, p2: float, c: np.ndarray, n: float
             ) -> tuple[np.ndarray, np.ndarray]:
    """Exact E[S] and Var[S] under random labelling of C cases among N deaths.

    `S = sum_{i != j} u_i u_j h_ij` for a symmetric 0/1 neighbour kernel `h`
    with zero diagonal and a case indicator `u` drawn as a simple random
    sample of size C without replacement. Expand S^2 and group the terms by
    how many *distinct* population indices they involve; a set of k distinct
    deaths is all in the sample with probability (C)_k / (N)_k, so only three
    graph summaries are needed:

        a2 = sum_{i != j} h_ij              = 2 * (number of edges)
        p2 = sum_v deg_v (deg_v - 1)        = ordered paths of length 2
        a22 = a2^2 - 2*a2 - 4*p2            = the 4-distinct-index remainder

    giving E[S] = a2 (C)_2/(N)_2 and
    E[S^2] = 2 a2 (C)_2/(N)_2 + 4 p2 (C)_3/(N)_3 + a22 (C)_4/(N)_4.

    The `a22` line is where the asymmetric k-NN kernel would have needed a
    separate mutual-neighbour term; symmetrizing the graph first (see
    `window_scores`) is what keeps this to three summaries. Checked against
    2,000-replicate simulation at C in {5, 10, 30, 100, 500} in
    `tests/test_concentration.py` — mean and sd match to within Monte Carlo
    error, which is the real guarantee that this derivation is right.
    """
    a22 = a2 ** 2 - 2 * a2 - 4 * p2
    p_2 = _falling(c, 2) / _falling(n, 2)
    p_3 = _falling(c, 3) / _falling(n, 3)
    p_4 = _falling(c, 4) / _falling(n, 4)
    mean = a2 * p_2
    second = 2 * a2 * p_2 + 4 * p2 * p_3 + a22 * p_4
    return mean, second - mean ** 2


def neighbour_graph(xy: np.ndarray, frac: float = NEIGHBOUR_FRACTION
                    ) -> tuple[sparse.csr_matrix, float, float]:
    """Symmetrized k-NN adjacency over one window's deaths, plus (a2, p2).

    `k` is `ceil(frac * N)`. The graph is symmetrized with OR rather than AND
    so that a death on the edge of a dense cluster still counts as that
    cluster's neighbour — the alternative (mutual k-NN only) prunes exactly
    the boundary cases a growing cluster is made of.

    Ties at identical coordinates (the facility artifact `core/spatial.py`
    documents — hospitals and jails where people are *pronounced*) are left
    in. They inflate degree for everyone at that address, cases and controls
    alike, and the random-labelling null prices that in; what they cannot
    price in is a substance *differentially* concentrated at facilities, which
    is why `geo.py` reports `pct_at_facility` as a caveat and this module
    inherits it rather than pretending it is solved.
    """
    n = len(xy)
    k = max(1, int(np.ceil(frac * n)))
    k = min(k, n - 1)
    _, nn = cKDTree(xy).query(xy, k=k + 1)      # column 0 is the point itself
    rows = np.repeat(np.arange(n), k)
    cols = np.asarray(nn)[:, 1:].ravel()
    h = sparse.coo_matrix((np.ones(len(rows)), (rows, cols)),
                          shape=(n, n)).tocsr()
    h = (h + h.T)
    h.data[:] = 1.0
    h.setdiag(0)
    h.eliminate_zeros()
    deg = np.asarray(h.sum(axis=1)).ravel()
    return h, float(deg.sum()), float((deg * (deg - 1)).sum())


def window_scores(
    bg: pd.DataFrame, cases: pd.DataFrame, frac: float = NEIGHBOUR_FRACTION,
    min_cases: int = MIN_CASES, max_fraction: float = MAX_FRACTION,
) -> pd.DataFrame:
    """z_s for every substance present in one trailing window.

    `bg` is the window's geocoded overdose deaths (`CaseNumber`, `x`, `y`);
    `cases` is the (`CaseNumber`, `substance`) mention rows restricted to the
    same window. Returns one row per substance with `n_cases`, `observed`,
    `expected`, `z`, and `withheld` — the reason a score was zeroed, or "" if
    it was not.
    """
    n = len(bg)
    cols = ["substance", "n_cases", "observed", "expected", "z", "withheld"]
    if n < 3 or cases.empty:
        return pd.DataFrame(columns=cols).set_index("substance")

    h, a2, p2 = neighbour_graph(bg[["x", "y"]].to_numpy(), frac)
    pos = pd.Series(np.arange(n), index=bg["CaseNumber"].to_numpy())

    subs = sorted(cases["substance"].unique())
    sidx = {s: i for i, s in enumerate(subs)}
    u = np.zeros((len(subs), n))
    row = cases["substance"].map(sidx).to_numpy()
    col = cases["CaseNumber"].map(pos).to_numpy()
    ok = ~pd.isna(col)
    u[row[ok].astype(int), col[ok].astype(int)] = 1.0

    c = u.sum(axis=1)
    # diag(U H U^T) without forming the full product.
    s_obs = ((u @ h) * u).sum(axis=1)
    mean, var = _moments(a2, p2, c, float(n))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (s_obs - mean) / np.sqrt(np.maximum(var, 1e-12))
    z = np.where(var > 1e-9, z, 0.0)

    withheld = np.full(len(subs), "", dtype=object)
    too_few = c < min_cases
    too_many = c > max_fraction * n
    withheld[too_few] = "n<min_cases"
    withheld[too_many] = "case_fraction>max"
    z = np.where(too_few | too_many, 0.0, z)

    return pd.DataFrame({
        "substance": subs, "n_cases": c.astype(int), "observed": s_obs,
        "expected": mean, "z": z, "withheld": withheld,
    }).set_index("substance")


def concentration_weight(z: np.ndarray | pd.Series) -> np.ndarray:
    """z_s -> a bounded concentration score in [0, 1].

    `2*Phi(z) - 1`, floored at 0: the one-sided normal confidence that the
    substance is *more* concentrated than overdose death itself. 0 at z = 0,
    0.95 at z = 1.96, and asymptoting to 1 — so the discount it drives (see
    `trends.rank_substances`) is bounded by construction and a spatially
    dispersed substance is never *penalised*, only left alone. Reading it as a
    calibrated probability would be overclaiming at these counts; it is used
    as a monotone squash, which is all the discount needs.
    """
    return np.clip(2 * norm.cdf(np.asarray(z, dtype=float)) - 1, 0.0, 1.0)


def load_spatial_cohort(
    mentions_path: Path = MENTIONS_PATH, quiet: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(geocoded deaths, mention rows) projected to UTM 11N, whole history.

    Deliberately *not* `core.spatial.load_points`, which windows to
    `config.SINCE` (2023) because the case-control analyses only trust the
    ME's coded flags from that era. This needs the full 2012- span, since the
    as-of sweep scores every quarter; coordinates are ~100% complete back to
    2012 (the `SINCE` cut is about flag provenance, not geocoding), so the
    earlier years are usable for a purely geometric statistic that reads no
    coded flags at all.
    """
    import geopandas as gpd

    co = load_cohort(verbose=not quiet)
    co = co[co["in_la_county"].fillna(False).astype(bool)].dropna(
        subset=["lat", "lon"])
    g = gpd.GeoDataFrame(
        co, geometry=gpd.points_from_xy(co["lon"], co["lat"]),
        crs="EPSG:4326").to_crs(UTM11N)
    bg = pd.DataFrame({
        "CaseNumber": co["CaseNumber"].to_numpy(),
        "quarter": co["quarter"].to_numpy(),
        "x": g.geometry.x.to_numpy(), "y": g.geometry.y.to_numpy(),
    })

    m = pd.read_parquet(mentions_path)
    m = m[~m["is_class_term"].fillna(False)].copy()
    m["substance"] = m["rollup"].fillna(m["canonical"])
    m = m[m["CaseNumber"].isin(set(bg["CaseNumber"]))]
    m = m.drop_duplicates(["CaseNumber", "substance"])
    return bg, m[["CaseNumber", "quarter", "substance"]]


def sweep_concentration(
    bg: pd.DataFrame, mentions: pd.DataFrame, quarters: list[pd.Timestamp],
    recent_quarters: int = 4, frac: float = NEIGHBOUR_FRACTION,
    min_cases: int = MIN_CASES,
) -> pd.DataFrame:
    """z_s for every (substance, as-of quarter), over trailing windows.

    The window at each as-of step is the same `recent_quarters` block that
    `trends.rank_substances` scores its numerator over, so the spatial
    evidence and the temporal evidence are about the same deaths. Trailing by
    construction: nothing after `as_of` is in the window, so replaying this
    into an as-of sweep cannot leak the future backward.
    """
    frames = []
    for i in range(recent_quarters, len(quarters) + 1):
        as_of = quarters[i - 1]
        win = set(quarters[i - recent_quarters:i])
        r = window_scores(bg[bg["quarter"].isin(win)],
                          mentions[mentions["quarter"].isin(win)],
                          frac=frac, min_cases=min_cases)
        r = r.reset_index()
        r["as_of"] = as_of
        frames.append(r)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
