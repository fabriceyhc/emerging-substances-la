"""
Where the data lives. All paths are relative to the repository root, so every
command is run from there.

    data/raw/lacme/    LACME overdose extract   -- NOT IN GIT, see below
    data/raw/nflis/    DEA NFLIS substance catalog
    data/raw/_geo/     LA County zipcode boundaries
    data/processed/    written by `extract`     -- NOT IN GIT
    results/<analysis>/  tables and figures     -- in git, aggregate only

**The LACME extract is gitignored and must stay that way.** It contains
decedent names. It is copied into `data/raw/lacme/` so the pipeline runs out
of the box on this machine, but `.gitignore` keeps it out of the repository
and nothing here should ever be changed to commit it.

`data/processed/` is also gitignored — the mentions table is one row per
decedent per substance, keyed by case number. It is cheap to regenerate
(`emerging lexicon build` then `emerging extract extract`, about a minute) and does not belong
in version control.

What *is* committed: the NFLIS catalog and the zipcode boundaries, both public,
and everything under `results/`, which is aggregate — counts, rates and test
statistics, no case numbers and no coordinates.
"""

from __future__ import annotations

from pathlib import Path

# Repository root, so paths resolve the same however the module is imported.
ROOT = Path(__file__).resolve().parent.parent

DATA_ROOT = ROOT / "data"
RAW_DIR = DATA_ROOT / "raw"
PROCESSED_DIR = DATA_ROOT / "processed" / "emerging"

LACME_PATH = RAW_DIR / "lacme" / "2012-01-2026-04-overdoses-regex-supplemented.csv"
NFLIS_PATH = RAW_DIR / "nflis" / "NFLIS-Substances_Excel_20260520-222324728.csv"
ZIPS_PATH = RAW_DIR / "_geo" / "zipcodes.geojson"

RESULTS_DIR = ROOT / "results"


def results_dir(analysis: str) -> Path:
    """The output directory for one analysis, created on demand.

    `results/` is the surface this project is reviewed from -- the tables and
    figures are what get checked for reasonableness -- so it is organised one
    directory per module, named after the module. `results/svtt/` pairs with
    `emerging/svtt.py`, `tests/test_svtt.py` and `docs/findings/svtt.md`.

    Figures live beside the tables that produced them rather than in a separate
    `figures/` tree: splitting on medium meant checking a figure against its
    own numbers took two directories. That is why there is no longer a
    `FIG_DIR`.

    The per-substance scans (`relrisk surface --substance`, `svtt scan
    --substance`, `treescan scan --reference`) template their filenames, so the
    old flat directory grew a new file into the same namespace on every run.
    Scoping them under their analysis keeps that growth legible.
    """
    d = RESULTS_DIR / analysis
    d.mkdir(parents=True, exist_ok=True)
    return d
