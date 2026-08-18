"""
Assert that the vendored cohort definition still matches upstream.

`emerging/cohort.py` is a copy of the authoritative definition in the
`predict_rosla` repository. The copy exists so this repository is standalone;
the cost is that upstream can change and nothing here will notice. Counts
would quietly stop being commensurable with `deaths` in
`panel_zip_quarter.parquet` and no error would ever be raised.

That is not hypothetical. The first vendoring transcribed `_CHRONIC_LIVER_RE`
from a terminal `repr()` that had been truncated mid-pattern, silently ending
at `"sequelae of "` and dropping three alternatives. It moved 6 deaths across
the cohort boundary and shifted every downstream denominator. Nothing failed.

These tests skip when `predict_rosla` is not importable, so they are a no-op
on a machine that only has this repository. Run them on a machine that has
both, before publishing any number compared against the modeling work:

    PYTHONPATH=/home/fabrice/predict_rosla pytest tests/ -v
"""

from __future__ import annotations

import inspect

import pytest

from emerging.ingest import cohort

upstream = pytest.importorskip(
    "rosla_nowcast.modeling.zipcode.monthly.signal_check.ccf",
    reason="predict_rosla not on the path; parity cannot be checked here",
)


def _body(fn) -> str:
    """Source with the docstring and all whitespace stripped."""
    parts = inspect.getsource(fn).split('"""')
    return "".join((parts[0] + parts[-1]).split())


@pytest.mark.parametrize("name", ["_ACUTE_INTOX_RE", "_CHRONIC_LIVER_RE"])
def test_regex_patterns_match_upstream(name: str) -> None:
    assert getattr(cohort, name).pattern == getattr(upstream, name).pattern, (
        f"{name} has drifted from predict_rosla. Re-copy it; do not edit by "
        f"hand from a terminal repr(), which is how it broke the first time."
    )


def test_substance_flag_list_matches_upstream() -> None:
    assert cohort.ALCOHOL_FILTER_SUBSTANCES == upstream.ALCOHOL_FILTER_SUBSTANCES


@pytest.mark.parametrize("name", ["filter_alcohol_only", "_norm_zip"])
def test_function_body_matches_upstream(name: str) -> None:
    assert _body(getattr(cohort, name)) == _body(getattr(upstream, name)), (
        f"{name} has drifted from predict_rosla (comparison ignores the "
        f"docstring and whitespace, so this is a real logic difference)."
    )


def test_filter_produces_identical_cohort() -> None:
    """The end-to-end check: same input, same surviving rows.

    Pattern equality is necessary but not sufficient — this is the test that
    would have caught a change in the boolean rule rather than the vocabulary.
    """
    import pandas as pd

    from emerging.paths import LACME_PATH

    if not LACME_PATH.exists():
        pytest.skip(f"LACME extract not present at {LACME_PATH}")

    df = pd.read_csv(LACME_PATH, low_memory=False)
    ours = cohort.filter_alcohol_only(df, verbose=False)
    theirs = upstream.filter_alcohol_only(df, verbose=False)
    assert len(ours) == len(theirs), (
        f"cohort size differs: vendored {len(ours):,} vs upstream "
        f"{len(theirs):,} — every downstream denominator is affected"
    )
    assert set(ours["CaseNumber"]) == set(theirs["CaseNumber"])
