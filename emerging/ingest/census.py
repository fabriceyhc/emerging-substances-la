"""
Resident population per zipcode per year, from the ACS 5-year estimates.

Every method in this repository up to now has been **compositional**: it asks
what share of overdose deaths name a substance, never what rate of death a
population experiences. That was a deliberate constraint — the ME extract
carries no denominator — and it shaped `geo.py` and `spacetime.py`, both of
which use other overdose deaths as controls precisely so that no population
count is needed.

The constraint was wrong, or at least stale. ACS table `B01001_001E` (total
population) is public, free, and already being pulled by
`predict_rosla/rosla_nowcast/data_prep/acs_ingest.py` — that script fetches it
and then discards it in `_derive`, keeping only the percentages it is built
from. This module fetches and *keeps* it.

**What a denominator buys, and what it costs.**

  - Buys: rates. "12 deaths per 100,000 residents" is the sentence every
    reviewer and every health officer asks for, and nothing here could say it.
  - Buys: spatial variation in temporal trends (game-plan P2's original
    method), which needs population per zip per period and was ruled out on
    the mistaken belief that we did not have one.
  - Costs: a population-denominated cluster scan mostly rediscovers *where
    overdose mortality is high*, which is driven by who lives and uses where,
    not by supply composition. That is the exact confound the case-control
    design in `geo.py` exists to avoid. Population analyses answer a different
    question; they do not replace the compositional ones.

**Caveats that matter for LA.**

  - ZCTAs are not USPS zipcodes. They are census approximations of them, and
    single-building or PO-box zips have no ZCTA at all. Coverage against the
    zips actually present in the death data is reported by `fetch`.
  - ACS 5-year estimates for year Y pool Y-4..Y, so they are already smoothed
    and lag a fast change. Skid Row is the case to worry about: a 5-year
    average badly describes a population that turns over.
  - **The series ends in 2024** and our study window runs to 2025Q4, so 2025
    is carried forward. ACS 2025 does not publish until late 2026.
  - Residential population is the wrong denominator for a place people travel
    to. 90013 has 14,891 residents and a daytime population several times
    that; a per-resident overdose rate there is an overestimate of individual
    risk and should be read as a place-based intensity instead.

The API key is read from `$CENSUS_API_KEY` or `--key`. It is never written
into this repository; request a free one at
https://api.census.gov/data/key_signup.html.

Usage:
    export CENSUS_API_KEY=...
    emerging census fetch
    emerging census show
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pandas as pd
import typer

from emerging.paths import RAW_DIR, ZIPS_PATH

app = typer.Typer(add_completion=False)

API = "https://api.census.gov/data/{year}/acs/acs5"
POP_VAR = "B01001_001E"          # Sex by Age -> total population
FIRST_YEAR, LAST_YEAR = 2013, 2024
POP_PATH = RAW_DIR / "census" / "acs_zip_population.csv"
# ACS uses large negative sentinels for suppressed cells.
JUMBO = -1e6


def _key(explicit: str | None) -> str:
    key = explicit or os.environ.get("CENSUS_API_KEY")
    if not key:
        # Convenience only: this machine already has one for the modeling repo.
        fallback = Path.home() / "predict_rosla" / "data" / "raw" / "census_api.key"
        if fallback.exists():
            key = fallback.read_text().strip()
    if not key:
        raise typer.BadParameter(
            "no Census API key. Set CENSUS_API_KEY or pass --key; free from "
            "https://api.census.gov/data/key_signup.html")
    return key


def _fetch_year(year: int, key: str) -> pd.DataFrame | None:
    """Total population per ZCTA for one ACS5 vintage."""
    params = {"get": POP_VAR, "for": "zip code tabulation area:*", "key": key}
    # ZCTAs are nested inside states only through 2019; from 2020 the pull is
    # national and has to be filtered afterwards.
    if year <= 2019:
        params["in"] = "state:06"

    for attempt in range(3):
        try:
            import requests
            r = requests.get(API.format(year=year), params=params, timeout=90)
            if r.status_code == 200:
                break
            time.sleep(2 * (attempt + 1))
        except Exception:
            time.sleep(2 * (attempt + 1))
    else:
        typer.echo(f"  {year}: FAILED after 3 attempts", err=True)
        return None

    rows = r.json()
    df = pd.DataFrame(rows[1:], columns=rows[0])
    df = df.rename(columns={"zip code tabulation area": "zipcode"})
    df["population"] = pd.to_numeric(df[POP_VAR], errors="coerce")
    df.loc[df["population"] <= JUMBO, "population"] = pd.NA
    df["year"] = year
    return df[["zipcode", "year", "population"]]


def la_zips() -> set[str]:
    """Zipcodes in the LA County boundary file."""
    import geopandas as gpd
    z = gpd.read_file(ZIPS_PATH)
    return set(pd.to_numeric(z["ZIPCODE"], errors="coerce")
                 .dropna().astype(int).astype(str))


def load_population(path: Path = POP_PATH) -> pd.DataFrame:
    """zipcode x year -> population, as strings/ints ready to join."""
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `emerging census fetch` first")
    d = pd.read_csv(path, dtype={"zipcode": str})
    return d


def population_by_quarter(
    quarters: list[pd.Timestamp], path: Path = POP_PATH
) -> pd.DataFrame:
    """zipcode x quarter population, interpolated within year and held flat
    outside the ACS range.

    Population is interpolated rather than stepped because a step at each
    January would put a discontinuity into every trend estimate at exactly the
    place a temporal scan looks. Years beyond the ACS range are carried
    forward from the last vintage — for the current window that means all of
    2025 uses the 2024 estimate, which is stated in the results rather than
    hidden.
    """
    pop = load_population(path)
    wide = pop.pivot_table(index="zipcode", columns="year",
                           values="population", aggfunc="first")
    # Mid-year reference date for each annual estimate.
    cols = {y: pd.Timestamp(f"{y}-07-01") for y in wide.columns}
    wide = wide.rename(columns=cols)

    out = pd.DataFrame(index=wide.index)
    for q in quarters:
        mid = q + pd.offsets.QuarterBegin(startingMonth=1) - pd.Timedelta(days=1)
        mid = q + (pd.offsets.QuarterEnd(0) - q) / 2
        dates = sorted(wide.columns)
        if mid <= dates[0]:
            out[q] = wide[dates[0]]
        elif mid >= dates[-1]:
            out[q] = wide[dates[-1]]
        else:
            hi = next(d for d in dates if d >= mid)
            lo = max(d for d in dates if d <= mid)
            if hi == lo:
                out[q] = wide[lo]
            else:
                w = (mid - lo) / (hi - lo)
                out[q] = wide[lo] * (1 - w) + wide[hi] * w
    return out


@app.command("fetch")
def fetch(
    out_path: Path = typer.Option(POP_PATH, "--out"),
    first: int = typer.Option(FIRST_YEAR, "--first-year"),
    last: int = typer.Option(LAST_YEAR, "--last-year"),
    key: str = typer.Option("", "--key"),
    la_only: bool = typer.Option(True, "--la-only/--all-ca"),
) -> None:
    """Download total population per zipcode per year and report coverage."""
    api_key = _key(key or None)
    keep = la_zips() if la_only else None

    frames = []
    for year in range(first, last + 1):
        df = _fetch_year(year, api_key)
        if df is None:
            continue
        if keep is not None:
            df = df[df["zipcode"].isin(keep)]
        frames.append(df)
        typer.echo(f"  {year}: {len(df):4d} zips, "
                   f"{df['population'].sum():,.0f} residents")
    if not frames:
        raise typer.Exit(1)

    pop = pd.concat(frames, ignore_index=True).sort_values(["zipcode", "year"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pop.to_csv(out_path, index=False)
    typer.echo(f"\nwrote {out_path}  ({len(pop):,} zip-years)")

    if keep is not None:
        missing = keep - set(pop["zipcode"])
        typer.echo(f"LA boundary file has {len(keep)} zips; ACS returned "
                   f"{pop['zipcode'].nunique()}. {len(missing)} have no ZCTA "
                   f"(PO-box and single-building zips): "
                   f"{', '.join(sorted(missing)[:10])}")


@app.command("show")
def show(
    path: Path = typer.Option(POP_PATH, "--path"),
    mentions_check: bool = typer.Option(True, "--check/--no-check",
                                        help="report coverage against the "
                                             "zips that actually hold deaths"),
) -> None:
    """Summarise the population file and check it covers the death data."""
    pop = load_population(path)
    typer.echo(f"{pop['zipcode'].nunique()} zips x "
               f"{pop['year'].nunique()} years ({pop['year'].min()}–"
               f"{pop['year'].max()}); "
               f"{pop.loc[pop['year'] == pop['year'].max(), 'population'].sum():,.0f}"
               f" residents in the last vintage")

    big = (pop[pop["year"] == pop["year"].max()]
           .nlargest(5, "population")[["zipcode", "population"]])
    typer.echo("\nlargest zips: "
               + ", ".join(f"{r.zipcode} {r.population:,.0f}"
                           for r in big.itertuples()))

    if not mentions_check:
        return
    from emerging.core.spatial import load_points
    from emerging.ingest.extract import MENTIONS_PATH
    bg, _ = load_points(MENTIONS_PATH)
    have = set(pop["zipcode"])
    deaths_by_zip = bg.groupby("zip").size()
    missing = deaths_by_zip[~deaths_by_zip.index.isin(have)]
    typer.echo(f"\ndeath zips without a population: {len(missing)} of "
               f"{len(deaths_by_zip)}, holding {int(missing.sum())} of "
               f"{int(deaths_by_zip.sum())} deaths "
               f"({100 * missing.sum() / deaths_by_zip.sum():.1f}%)")
    if len(missing):
        typer.echo("  worst: " + ", ".join(
            f"{z} ({n})" for z, n in missing.nlargest(8).items()))


if __name__ == "__main__":
    app()
