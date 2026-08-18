"""
Overdose-death cohort definition.

**This is a vendored copy, not an original.** The authoritative version lives in
the `predict_rosla` repository at
`rosla_nowcast/modeling/zipcode/monthly/signal_check/ccf.py`, and is itself
kept in sync with the epi project (`epi/scripts/utils.py:filter_alcohol_only`).
It is duplicated here only because this repository is standalone and cannot
import from that tree.

**The duplication is a real liability.** The whole point of the shared
definition is that substance counts here stay commensurable with `deaths` in
`panel_zip_quarter.parquet`; if the upstream rule changes and this copy does
not, every count in this repository silently stops matching the modeling panel
and nothing will fail loudly. Before publishing any number that gets compared
to the modeling work, diff this file against upstream.

Copied 2026-08-10 from `predict_rosla` @ e357baf.

The rule, in words. Among deaths where alcohol is the only substance flagged
(Alcohol=1, every other individual substance flag 0):

  - KEEP if the structured cause fields show acute-intoxication language
    (`intoxic...`, ethanol/alcohol toxicity|poisoning, or `acute` within 5
    characters of alcohol/ethanol/alcoholism), whatever the Mode — this
    catches acute cases the office mislabeled; or
  - KEEP if Mode contains 'ACCIDENT' and there are no chronic-disease markers;
  - REMOVE everything else — chronic alcoholism and liver-disease sequelae,
    and suicide/undetermined/homicide without acute-intoxication language.

Chronic alcoholic liver disease is removed even when coded Accident, because
those are chronic sequelae rather than acute toxicity.
"""

from __future__ import annotations

import re

import pandas as pd
import typer

# The 8 coded flags the ME panel carries. "Alcohol only" means the first is
# set and none of the others are.
ALCOHOL_FILTER_SUBSTANCES = [
    "Methamphetamine", "Heroin", "Cocaine", "Fentanyl", "Alcohol",
    "Prescription.opioids", "Benzodiazepines", "Others",
]

_ACUTE_INTOX_RE = re.compile(
    "intoxic|ethanol toxicit|alcohol toxicit|ethanol poison|alcohol poison|"
    "acute.{0,5}(?:alcohol|ethanol|ethanolism|alcoholism)"
)

_CHRONIC_LIVER_RE = re.compile(
    "cirrhos|laennec|laënnec|fatty liver|steato|hepat|pancreat|"
    "chronic.{0,15}(?:alcohol|ethanol|liver)|alcoholic liver|liver failure|"
    "esophageal var|portal hypertension|ascites|end.stage liver|"
    "sequelae of (?:hepatic|liver)|alcohol use disorder|alcohol abuse"
)
# Verify against upstream rather than trusting the eye:
#   from rosla_nowcast.modeling.zipcode.monthly.signal_check import ccf
#   assert ccf._CHRONIC_LIVER_RE.pattern == _CHRONIC_LIVER_RE.pattern
# The first vendoring of this pattern was transcribed from a truncated
# terminal `repr()` and silently ended at "sequelae of ", dropping the
# hepatic/liver qualifier and the two alcohol-use-disorder alternatives. That
# moved 6 deaths across the cohort boundary and shifted every downstream
# denominator. It failed no test and raised no error — which is the whole
# argument for diffing this file against upstream before publishing.


def _norm_zip(s: pd.Series) -> pd.Series:
    """Normalize to a zero-padded 5-digit zipcode string, NaN for invalid."""
    s = s.astype("string")
    extracted = s.str.extract(r"(\d{5})", expand=False)
    return extracted


def filter_alcohol_only(df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    """Drop chronic alcohol-only deaths; keep acute/accidental cases.

    Cause fields (CauseA..CauseD, CauseOther) are matched case-insensitively.
    Requires `ALCOHOL_FILTER_SUBSTANCES`, a `Mode` column, and at least one
    cause field.
    """
    flags = {c: pd.to_numeric(df[c], errors="coerce").fillna(0)
             for c in ALCOHOL_FILTER_SUBSTANCES}
    non_alc = [c for c in ALCOHOL_FILTER_SUBSTANCES if c != "Alcohol"]
    alcohol_only = (flags["Alcohol"] == 1) & (sum(flags[c] for c in non_alc) == 0)

    mode_norm = df["Mode"].fillna("").astype(str).str.upper()
    is_accidental = mode_norm.str.contains("ACCIDENT", regex=False, na=False)

    cause_candidates = ["CauseA", "CauseB", "CauseC", "CauseD", "CauseOther"]
    col_by_lower = {c.lower(): c for c in df.columns}
    text_cols = [col_by_lower[name.lower()] for name in cause_candidates
                 if name.lower() in col_by_lower]
    combined = df[text_cols[0]].fillna("").astype(str).str.lower()
    for c in text_cols[1:]:
        combined = combined.str.cat(df[c].fillna("").astype(str).str.lower(),
                                    sep=" ")
    is_acute_intox = combined.str.contains(_ACUTE_INTOX_RE, regex=True, na=False)
    is_chronic_liver = combined.str.contains(_CHRONIC_LIVER_RE, regex=True,
                                             na=False)

    # Remove alcohol-only deaths that lack acute-intoxication language and are
    # either non-accidental OR chronic alcoholic liver disease (the latter even
    # when coded Accident — chronic sequelae, not acute toxicity).
    to_remove = alcohol_only & ~is_acute_intox & (~is_accidental | is_chronic_liver)
    if verbose:
        n_chronic_acc = int((alcohol_only & is_accidental & is_chronic_liver
                             & ~is_acute_intox).sum())
        typer.echo(
            f"  Alcohol-only filter: {int(alcohol_only.sum()):,} alcohol-only; "
            f"kept {int((alcohol_only & ~to_remove).sum()):,} acute/accidental; "
            f"removed {int(to_remove.sum()):,} chronic/other "
            f"(incl. {n_chronic_acc} accidental-coded chronic-liver)."
        )
    return df[~to_remove].copy()
