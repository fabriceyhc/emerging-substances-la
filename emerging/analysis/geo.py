"""
Test whether a flagged substance is spatially localized within LA County.

This is the place the spatial methods finally earn their keep. SaTScan and
sparr were rejected as *detectors* — at 9 lifetime xylazine deaths there is
nothing to scan for — but once `trends` has named a candidate on temporal
evidence, "and where" is a separate, answerable question, and the answer is
what turns a statistic into an outreach target.

**The design is case-control, not case-only.** Overdose deaths are already
severely clustered: they concentrate in Skid Row, Antelope Valley, and the
Harbor. A substance whose deaths cluster there too has told you nothing. So
the null is not spatial randomness over the county — it is *the same number of
deaths drawn from the other overdose deaths in the same window*. The test asks
whether this substance is more concentrated than overdose death itself, which
is the only version of the question with a public-health answer.

Statistic: mean nearest-neighbour distance among the substance's death
locations, in metres (UTM 11N). Small means clustered. The permutation null
resamples n background deaths 999 times, so the n-dependence of nearest-
neighbour distance cancels — an n=9 substance is compared against n=9 draws.
Raw NND is meaningless across substances for exactly that reason (a random
draw of 9 gives 11.4km, of 367 gives 1.3km); only the ratio is comparable.

**What this is, in the literature.** Drawing a random n-subset of the fixed
point set is precisely a *random labelling* null for a marked point process,
so this is a random-labelling permutation test with a Clark-Evans-style
nearest-neighbour ratio as its summary. Its closest named relative is the
**Cuzick-Edwards k-nearest-neighbour test** (1990), the standard case-control
clustering test, which shares the label-permutation null but counts how many
of each case's k nearest neighbours are cases rather than measuring distance.
The canonical alternatives are **Diggle-Chetwynd** D(s) = K_cases - K_controls
(1991) and **Kelsall-Diggle / sparr** relative-risk surfaces.

Chosen over those because there is no bandwidth or scale radius to pick — at
n=9-39 that choice would drive the answer — it yields one comparable number
per substance across 42 of them, and the p-value is exact rather than
asymptotic. Given up: Cuzick-Edwards has established power properties and this
does not, and D(s) identifies the *scale* of clustering where this collapses
to the nearest-neighbour scale. Promote a substance to sparr/Diggle-Chetwynd
before publishing a spatial claim about it.

Note that including the substance's own deaths in the resampling pool is
required here, not a defect: under random labelling the null is "which n of
these N deaths are cases", so the pool is the union by construction. What it
does mean is that the contrast vanishes as n/N approaches 1 — hence
`MAX_CASE_FRACTION`.

**Three data defects this has to control for, all found in the extract:**

1. *Foreign geocodes.* 46 deaths land at (-22.14, 17.22) — Namibia. The
   `in_la_county` flag catches these; nothing here runs without it.
2. *Facility addresses.* 16% of deaths share a coordinate with at least four
   others: hospitals, jails, treatment facilities. These are nearest-neighbour
   distances of zero and they are an artifact of where people are *pronounced*,
   not where the drug is. The permutation null absorbs the county-wide rate of
   this, but a substance dying disproportionately in hospital would still read
   as clustered, so `pct_at_facility` is reported per substance as the check.
3. *Death vs incident address.* `used_death_address` marks rows where the
   incident location was missing and the death location stood in.

Power is thin and the output says so. At n=10 this test detects only extreme
localization; a non-significant result is not evidence of dispersion.

Usage:
    emerging geo cluster
    emerging geo cluster --since 2022-01-01 --min-n 6
    emerging geo plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from scipy.spatial import cKDTree

from emerging.config import CUTOFF, SINCE
from emerging.ingest.extract import MENTIONS_PATH
from emerging.core.spatial import MAX_CASE_FRACTION, UTM11N, load_points
from emerging.paths import ZIPS_PATH, results_dir
from emerging.viz import INK, INK_2, MUTED, ORANGE

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("geo")


N_PERM = 999
MIN_N = 8


def _mean_nnd(xy: np.ndarray) -> float:
    """Mean distance from each point to its nearest other point, in metres."""
    if len(xy) < 2:
        return np.nan
    d, _ = cKDTree(xy).query(xy, k=2)
    return float(d[:, 1].mean())


def cluster_test(
    bg: pd.DataFrame, cases: set, rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> dict:
    """Permutation case-control test for excess localization."""
    idx = bg.index[bg["CaseNumber"].isin(cases)]
    n = len(idx)
    xy = bg.loc[idx, ["x", "y"]].to_numpy()
    obs = _mean_nnd(xy)

    all_xy = bg[["x", "y"]].to_numpy()
    null = np.empty(n_perm)
    for b in range(n_perm):
        pick = rng.choice(len(all_xy), size=n, replace=False)
        null[b] = _mean_nnd(all_xy[pick])

    # One-sided: clustering means a *smaller* nearest-neighbour distance.
    p = (1 + int((null <= obs).sum())) / (n_perm + 1)
    frac = n / len(bg)
    return {
        "n": n,
        "case_frac": frac,
        "out_of_scope": bool(frac > MAX_CASE_FRACTION),
        "mean_nnd_m": obs,
        "null_nnd_m": float(null.mean()),
        "localization": float(null.mean() / obs) if obs else np.nan,
        "z": float((obs - null.mean()) / null.std(ddof=1)),
        "p_clustered": p,
        "pct_at_facility": 100 * float(bg.loc[idx, "at_facility"].mean()),
        "top_zip": (bg.loc[idx, "zip"].value_counts().index[0]
                    if n else ""),
        "top_zip_pct": (100 * bg.loc[idx, "zip"].value_counts().iloc[0] / n
                        if n else np.nan),
    }


@app.command("cluster")
def cluster(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    since: str = typer.Option(SINCE, "--since"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    min_n: int = typer.Option(MIN_N, "--min-n"),
    n_perm: int = typer.Option(N_PERM, "--n-perm"),
    top: int = typer.Option(30, "--top"),
    seed: int = typer.Option(0, "--seed"),
) -> None:
    """Rank substances by how localized their deaths are, against the
    background geography of overdose death."""
    bg, sets = load_points(mentions, since, cutoff)
    rng = np.random.default_rng(seed)

    rows = []
    for sub, cases in sets.items():
        if len(cases) < min_n:
            continue
        r = cluster_test(bg, cases, rng, n_perm)
        r["substance"] = sub
        rows.append(r)

    res = (pd.DataFrame(rows).set_index("substance")
           .sort_values("p_clustered"))
    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "clustering.csv")

    # County-wide facility rate, as the reference for pct_at_facility.
    fac = 100 * float(bg["at_facility"].mean())
    typer.echo(f"{len(bg):,} OD deaths {since}..{cutoff} in LA County; "
               f"{n_perm} permutations; {fac:.0f}% at a shared facility "
               f"address\n")

    scoped = res[~res["out_of_scope"]]
    cols = ["n", "mean_nnd_m", "null_nnd_m", "localization", "p_clustered",
            "pct_at_facility", "top_zip", "top_zip_pct"]
    typer.echo(scoped.head(top)[cols].to_string(
        float_format=lambda v: f"{v:.3g}"))

    sig = scoped[scoped["p_clustered"] < 0.05]
    typer.echo(f"\n{len(sig)} of {len(scoped)} substances are more localized "
               f"than the overdose background at p<0.05.")
    typer.echo("localization > 1 means tighter than background; p is one-sided "
               "and unadjusted for\nthe number of substances tested — with "
               f"{len(scoped)} tests, expect ~{0.05 * len(scoped):.0f} at "
               "p<0.05 by chance.")

    excluded = res[res["out_of_scope"]]
    if len(excluded):
        typer.echo(
            f"\nWITHHELD — more than {100 * MAX_CASE_FRACTION:.0f}% of the "
            "cohort, so the case-control contrast degenerates:")
        for name, r in excluded.iterrows():
            typer.echo(f"  {name:20s} n={int(r['n']):5d} "
                       f"({100 * r['case_frac']:.0f}% of all OD deaths)")
        typer.echo(
            "  Under random labelling the null draw is mostly these same "
            "deaths, so the test\n  compares the substance against itself and "
            "returns ~1 whatever the geography.\n  Values are in the CSV with "
            "out_of_scope=True; do not quote them. These substances\n  need a "
            "population denominator, not a case-control one.")
    typer.echo(f"\nwrote {out_dir / 'geo_clustering.csv'}")


@app.command("plot")
def plot(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    since: str = typer.Option(SINCE, "--since"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    # Default set pairs the two substances that pass the test (PCP, ketamine)
    # with four flagged candidates that do not, so the figure shows what
    # localized and dispersed actually look like side by side.
    substances: str = typer.Option(
        "PCP,Ketamine,Lidocaine,para-Fluorofentanyl,Carfentanil,Xylazine",
        "--substances"),
    seed: int = typer.Option(0, "--seed"),
) -> None:
    """Map each candidate's deaths over the overdose-death background."""
    import geopandas as gpd

    bg, sets = load_points(mentions, since, cutoff)
    rng = np.random.default_rng(seed)
    names = [s.strip() for s in substances.split(",")
             if s.strip() in sets.index]

    zips = None
    if ZIPS_PATH.exists():
        zips = gpd.read_file(ZIPS_PATH).to_crs(UTM11N)

    ncol = 3
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.9 * nrow))
    axes = np.atleast_1d(axes).ravel()

    # Frame on the populated mainland. Three of 7,581 deaths are on Catalina,
    # ~80km offshore; framing on the raw extent stretches the y-axis to reach
    # them and squeezes the entire county into the top third of every panel.
    # Quantile limits drop those three from view — noted in the caption.
    qx = bg["x"].quantile([0.002, 0.998]).to_numpy()
    qy = bg["y"].quantile([0.002, 0.998]).to_numpy()
    pad = 0.05 * (qx[1] - qx[0])
    xlim = (qx[0] - pad, qx[1] + pad)
    ylim = (qy[0] - pad, qy[1] + pad)
    n_off = int(((bg["x"] < xlim[0]) | (bg["x"] > xlim[1])
                 | (bg["y"] < ylim[0]) | (bg["y"] > ylim[1])).sum())

    for ax, name in zip(axes, names):
        if zips is not None:
            zips.boundary.plot(ax=ax, color=MUTED, linewidth=0.25, alpha=0.55)
        ax.scatter(bg["x"], bg["y"], s=1.6, c=MUTED, alpha=0.30, linewidths=0,
                   zorder=2)
        sel = bg[bg["CaseNumber"].isin(sets[name])]
        # Shrink and soften the marks as n grows, or a few hundred points read
        # as one solid blob and the pattern the test found is invisible.
        size = float(np.clip(1400 / max(len(sel), 1), 9, 46))
        ax.scatter(sel["x"], sel["y"], s=size, c=ORANGE,
                   alpha=0.85 if len(sel) < 60 else 0.6,
                   edgecolors="white", linewidths=0.6 if len(sel) < 200 else 0.3,
                   zorder=4)
        # adjustable="box" shrinks the axes box to honour the aspect ratio.
        # The default, "datalim", instead widens the limits — which would undo
        # the framing above and put the islands back.
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

        r = cluster_test(bg, sets[name], rng, 199)
        if r["out_of_scope"]:
            line2 = (f"{100 * r['case_frac']:.0f}% of all OD deaths — "
                     "case-control contrast degenerates, no score")
            color2 = ORANGE
        else:
            star = "*" if r["p_clustered"] < 0.05 else ""
            line2 = (f"localization {r['localization']:.2f}x  "
                     f"p={r['p_clustered']:.3f}{star}   "
                     f"top zip {r['top_zip']} ({r['top_zip_pct']:.0f}%)")
            color2 = INK
        ax.set_title(f"{name}   n={r['n']}", fontsize=9.5, color=INK,
                     loc="left", pad=14)
        ax.annotate(line2, (0, 1.005), xycoords="axes fraction", fontsize=9,
                    color=color2, ha="left", va="bottom")
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(False)

    for ax in axes[len(names):]:
        ax.set_visible(False)

    fig.suptitle(
        "Where the flagged substances are dying, against the overdose "
        "background\n"
        f"grey = all {len(bg):,} LA County overdose deaths {since[:4]}–"
        f"{cutoff[:4]}; orange = deaths naming the substance. Localization is "
        "the ratio of background to observed mean nearest-neighbour distance.",
        fontsize=11.5, color=INK, x=0.008, ha="left")
    fig.text(0.008, 0.004,
             "p is a one-sided permutation test against n deaths resampled "
             "from the background, so it asks whether the substance is more "
             "concentrated than overdose death already is — not whether it is "
             f"concentrated at all. Underpowered below n~20. Frame is the "
             f"populated county; {n_off} offshore deaths fall outside it.",
             fontsize=8, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.015, 1, 0.94))
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "clusters.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    typer.echo(f"wrote {out_dir / 'geo_clusters.png'}")


if __name__ == "__main__":
    app()
