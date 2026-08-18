"""
Shared spatial primitives: the projection, the geocoded cohort, and the two
guards every spatial analysis needs.

These lived in `geo.py`, which meant `spacetime`, `svtt` and `relrisk` all
imported from the case-control clustering module to get a CRS and a point
loader. The import graph then claimed, falsely, that the space-time scan
depends on the nearest-neighbour permutation test. It does not; they merely
share a way of turning the cohort into projected points.

What stayed behind in `geo.py` is the case-control method itself -- `N_PERM`,
`MIN_N`, `_mean_nnd`, `cluster_test`. Those are meaningful only for that test.
What moved here is what more than one analysis genuinely needs.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer

from emerging.config import CUTOFF, SINCE
from emerging.ingest.extract import MENTIONS_PATH, load_cohort

# LA County in UTM zone 11N. Metres, so nearest-neighbour distances are
# directly interpretable and no great-circle correction is needed at this
# extent.
UTM11N = "EPSG:32611"

# Coordinates shared by at least this many deaths are facilities -- hospitals,
# jails, treatment centres -- not drug geography. They are nearest-neighbour
# distances of zero and an artifact of where people are *pronounced*.
FACILITY_MIN = 5

# Above this share of the cohort, a case-control contrast stops being
# interpretable and the result is withheld rather than printed. Under random
# labelling the null is "which n of these N deaths are cases is arbitrary", so
# as n/N approaches 1 the null draw becomes the case set itself and the
# contrast vanishes. It is not bias — the statistic is still unbiased — it is
# that the question "is this substance more clustered than overdose death?"
# degenerates when the substance *is* most of overdose death. Methamphetamine
# is 61% of the cohort and fentanyl 58%; both need a population denominator,
# not a case-control one.
MAX_CASE_FRACTION = 0.25


def load_points(
    mentions_path: Path = MENTIONS_PATH, since: str = SINCE,
    cutoff: str = CUTOFF,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (projected background OD deaths, case -> substance set)."""
    import geopandas as gpd

    co = load_cohort()
    co = co[(co["death_date"] >= pd.Timestamp(since))
            & (co["death_date"] <= pd.Timestamp(cutoff))]

    n_before = len(co)
    co = co[co["in_la_county"].fillna(False).astype(bool)]
    co = co.dropna(subset=["lat", "lon"])
    if n_before - len(co):
        typer.echo(f"dropped {n_before - len(co)} deaths outside LA County or "
                   "without coordinates")

    g = gpd.GeoDataFrame(
        co, geometry=gpd.points_from_xy(co["lon"], co["lat"]), crs="EPSG:4326"
    ).to_crs(UTM11N)
    bg = pd.DataFrame({
        "CaseNumber": g["CaseNumber"].to_numpy(),
        "x": g.geometry.x.to_numpy(), "y": g.geometry.y.to_numpy(),
        "lon": co["lon"].to_numpy(), "lat": co["lat"].to_numpy(),
        # ZIPCODE arrives as a float, so a plain astype(str) yields "90013.0".
        "zip": pd.to_numeric(co["ZIPCODE"], errors="coerce")
                 .astype("Int64").astype(str).to_numpy(),
    })

    # Facility flag: coordinates shared by many deaths are hospitals, jails and
    # treatment centres, not drug geography.
    dup = bg.groupby(["x", "y"])["CaseNumber"].transform("size")
    bg["at_facility"] = dup >= FACILITY_MIN

    m = pd.read_parquet(mentions_path)
    m = m[~m["is_class_term"].fillna(False)].copy()
    m["s"] = m["rollup"].fillna(m["canonical"])
    m = m[m["CaseNumber"].isin(set(bg["CaseNumber"]))]
    sets = m.drop_duplicates(["CaseNumber", "s"]).groupby("s")["CaseNumber"].apply(set)
    return bg, sets
