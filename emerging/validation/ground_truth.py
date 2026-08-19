"""
Score the EB05 alarm sweep against a hand-curated reference of known local
emergences, using interval overlap instead of `backtest`'s single-quarter
"first top-10" read.

**Why this exists, and what it adds over `backtest`.** `backtest` answers "did
the ranking notice this substance at all, and how did its rank evolve" — it
has no independent ground truth for *when* an emergence actually started or
ended, so it cannot say whether an alarm episode was well-timed, early, late,
or spuriously extended. This module adds that reference:
`data/raw/ground_truth/known_emergences.csv`, one row per known local
emergence interval (a substance can have more than one, non-contiguous), and
scores the alarm sweep's detected episodes (`trends._episode_spans`) against
it with a temporal-IoU-style metric.

**Labelled substances are not all positives.** Every substance labelled
before this only confirmed a *real* emergence, which meant precision could
only ever be measured on cases already known to be true — never a genuine
test of whether EB05 ever fires on something that isn't real. `no_emergence`
marks a documented negative: a substance that breached the alarm but whose
own raw quarterly count and share trend say it isn't a real local emergence
(Morphine: a substance in clear multi-year structural *decline*, EB05>1.5
traced to a single-quarter blip). `interval_start`/`interval_end` still hold
values for that row for CSV-parsing convenience, but `score` ignores them
and scores against zero ground-truth-positive quarters — so any alarm on it
counts as pure false alarm, the first real precision test in this file.

**Provenance matters, and the CSV records it per row.** Two ways an interval
got there:

- `raw-count-rule`: a fixed, stated, EB05-independent rule applied to the
  substance's raw quarterly mention count — first run of >=2 consecutive
  quarters each with count >=2, following >=2 consecutive quarters each <=1;
  ends at the last quarter of that run before another >=2-consecutive-<=1
  run, or is left open ("ongoing") if the extract ends first. This has to be
  independent of EB05 (or any other fitted statistic) or scoring the
  detector against it would be circular — validating EB05 against labels
  derived from EB05.
- `manual`: the rule above fails on some substances -- carfentanil and
  xylazine are both too sparse and volatile at low n to ever post two
  consecutive qualifying quarters, so their intervals are eyeballed directly
  from the same raw counts and the reasoning is written into `notes`.

**External context lives in a separate file, `data/raw/ground_truth/
health_notices.csv`, and is deliberately not folded into the interval
columns.** It is not used as ground truth for LA's own timing — the point of
local surveillance is to see whether LA lags, leads, or misses the external
picture, so stamping outside dates onto local intervals would beg the
question this project exists to answer. It is there for a reader to compare
by eye against `interval_start`. `known_emergences.csv`'s `national_context`
column carries only a one-line pointer into it now; the full sourcing lives
in `health_notices.csv` so it isn't duplicated (and can't drift) across two
files.

`health_notices.csv` is one row per (substance, jurisdiction) — a substance
can have several, at different government levels (county/state/federal) plus
non-government sources (forensic-toxicology alerts) where relevant. Two
columns carry the honesty this needed after the first pass didn't have it:
`evidence_basis` (did the notice cite actual local detection, or infer it
from a neighboring jurisdiction — these are not the same claim, and
conflating them was a real mistake caught by directly reading LA County's own
xylazine alert, which turned out to cite zero LA-specific data) and
`source_verification` (was this URL actually fetched and read, or is it a
WebSearch tool's synthesized summary of search results — a meaningfully
weaker form of evidence that the first pass didn't distinguish). Coverage is
real and incomplete, and the file says so per row rather than silently
falling back to a weaker source: 3 of the 5 substances have **no** LA County
or California-level notice at all, despite searching. That absence is itself
a finding — the substances *without* a local notice (bromazolam, carfentanil,
para-fluorofentanyl) are just as real in our own ME data as xylazine, which
did get one. A notices-based reference is biased toward whichever substances
get media/political attention, not necessarily whichever are most dangerous
or arrive first — worth knowing before trusting it as a complete record of
anything.

**The metric.** For a substance's ground-truth intervals and detected-alarm
episodes (both lists, since either can have more than one), each side's
quarters are unioned first, then compared as sets:

    iou       = |intersection| / |union|
    recall    = |intersection| / |ground-truth quarters|   -- coverage
    precision = |intersection| / |detected quarters|        -- false-alarm rate

`recall` is the more diagnostic number for a rare-event detector: missing
part of a real emergence is usually worse than an alarm that runs a quarter
long on either side, so a design that trades `precision` for `recall` is
often the right trade here, not a wash.

**A third, stronger external comparator: real NFLIS counts, not notice
dates.** `national_context` (prose from web search) and `health_notices.csv`
(dated alerts, jurisdiction-tiered) both answer "when did officials say
something" — a different question from "when did the substance actually show
up in a forensic count," and the LA County xylazine alert
(`health_notices.csv`) already demonstrated those two things aren't the same
— a notice's date can predate any local evidence for it. `nflis_ca_counts` loads a real
DQS export from DEA's NFLIS-Drug system — California, semi-annual drug
identification counts by substance, the same forensic-lab program this
project's own substance catalog comes from — and `nflis_onset` applies a
rule in the same spirit as `trends.py`'s `raw-count-rule` (a floor, not a
ratio, since half-year state counts run in the tens to hundreds rather than
the single digits our own quarterly LA counts do): first period >= 10,
preceded by two consecutive periods < 10.  This has caught real problems
twice already: the source CSV omits zero-count periods as rows entirely
rather than writing them as 0, so a naive scan of consecutive *rows* mistakes
periods that are actually years apart for adjacent ones (see
`_onset_from_sparse` reindexing to a full half-year calendar before scanning
— para-fluorofentanyl's onset was undetectable without this fix); and the
mechanical rule fires a false positive on mitragynine (a single 2015H1 spike
of 33, collapsing straight back to 2 the next half-year) that has to be read
against the series rather than trusted blind, exactly the same discipline
`docs/GPS_V2_DESIGN.md` had to apply to `_ever_independent`.

Usage:
    emerging ground-truth score
    emerging ground-truth score --threshold 1.0
    emerging ground-truth nflis-lag
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer

from emerging.analysis.trends import (ALARM_THRESHOLD, _episode_spans,
                                      load_quarterly, sweep_eb05)
from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import (KNOWN_EMERGENCES_PATH, NFLIS_CA_COUNTS_PATH,
                            results_dir)

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("ground_truth")
SWEEP_CACHE = results_dir("trends") / "eb05_sweep.csv"


def _parse_quarter(s: str) -> pd.Timestamp:
    return pd.Period(s, freq="Q").to_timestamp()


def load_ground_truth(
    path: Path = KNOWN_EMERGENCES_PATH, last_q: pd.Timestamp | None = None,
) -> pd.DataFrame:
    """Read the reference CSV, resolving `ongoing` rows to `last_q`.

    `last_q` is the last complete quarter the caller is scoring against --
    required whenever any row is `ongoing`, since an open interval has no
    meaning without a "now" to close it at.
    """
    gt = pd.read_csv(path)
    gt["interval_start"] = gt["interval_start"].map(_parse_quarter)

    def _end(row: pd.Series) -> pd.Timestamp:
        if bool(row["ongoing"]) or pd.isna(row["interval_end"]):
            if last_q is None:
                raise ValueError(
                    f"{row['substance']}: an 'ongoing' row needs last_q")
            return last_q
        return _parse_quarter(row["interval_end"])

    gt["interval_end"] = gt.apply(_end, axis=1)
    return gt


def _q_ordinal(t: pd.Timestamp) -> int:
    """Quarters since year 0, so interval arithmetic is plain integer set
    math instead of Timestamp offset juggling."""
    return t.year * 4 + (t.quarter - 1)


def _interval_quarters(start: pd.Timestamp, end: pd.Timestamp) -> set[int]:
    return set(range(_q_ordinal(start), _q_ordinal(end) + 1))


def interval_overlap(
    gt: list[tuple[pd.Timestamp, pd.Timestamp]],
    det: list[tuple[pd.Timestamp, pd.Timestamp]],
) -> dict:
    """Temporal-IoU-style overlap between two sets of (start, end) intervals.

    Each side is unioned across its own intervals before comparison, so a
    substance with two ground-truth episodes and one long detected episode
    (or vice versa) still scores on the actual overlap rather than requiring
    the caller to pre-match individual episodes.
    """
    gt_q: set[int] = set().union(*(_interval_quarters(s, e) for s, e in gt)) if gt else set()
    det_q: set[int] = set().union(*(_interval_quarters(s, e) for s, e in det)) if det else set()
    inter = gt_q & det_q
    union = gt_q | det_q
    return {
        "iou": len(inter) / len(union) if union else float("nan"),
        "recall": len(inter) / len(gt_q) if gt_q else float("nan"),
        "precision": len(inter) / len(det_q) if det_q else float("nan"),
        "gt_quarters": len(gt_q),
        "det_quarters": len(det_q),
        "overlap_quarters": len(inter),
    }


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n, at confidence
    level implied by `z` (1.96 -> ~95%).

    Preferred over the normal approximation (`p +/- z*sqrt(p(1-p)/n)`) at the
    small n this file actually has -- the normal interval can go negative or
    above 1, which a proportion can't, and is a poor approximation exactly
    where it matters most here, n in the single digits. Returns `(0.0, 1.0)`
    for `n == 0` -- no data, maximum uncertainty, not an error.
    """
    if n == 0:
        return 0.0, 1.0
    p = k / n
    denom = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / denom
    margin = (z / denom) * ((p * (1 - p) / n + z ** 2 / (4 * n ** 2)) ** 0.5)
    return max(0.0, centre - margin), min(1.0, centre + margin)


def precision_rate(per_sub: pd.DataFrame) -> dict:
    """Substance-level false-positive rate, with a Wilson CI, computed over
    the labelled substances a model actually alarmed on.

    **Why substance-level, not pooled over alarm-quarters.** Pooling every
    alarmed *quarter* across substances into one big n would inflate the
    apparent sample size with correlated observations -- the four quarters
    of Lidocaine's episode are not four independent trials, they're one
    sustained mistake, and a Wilson interval computed as if they were
    independent would be falsely narrow. One substance, one trial: did any
    part of this model's alarm on it land inside real ground truth
    (`overlap_quarters > 0`) or not (`overlap_quarters == 0`, e.g. Morphine,
    Lidocaine, Levamisole). A partial hit like Benzylfentanyl's 0.33
    precision still counts as a "success" here -- this is answering "did the
    alarm ever touch something real," a different and coarser question than
    `mean_precision`'s "how much of the alarm was real," not a replacement
    for it.

    `per_sub` must have `det_quarters` and `overlap_quarters` columns, e.g.
    `score`'s own `rows` before being written to `overlap.csv`, or
    `benchmark.py`'s `per_sub` from `score_model`.
    """
    alarmed = per_sub[per_sub["det_quarters"] > 0]
    n = len(alarmed)
    k = int((alarmed["overlap_quarters"] > 0).sum())
    lo, hi = wilson_ci(k, n)
    return {"precision_rate": k / n if n else float("nan"),
            "precision_rate_n": n, "precision_rate_k": k,
            "precision_rate_lo": lo, "precision_rate_hi": hi}


# Half-year floor for the NFLIS onset rule. Chosen by inspecting the actual
# series (`nflis_ca_counts.csv` spans single digits to 400+), not tuned to
# produce a particular answer -- see the module docstring.
NFLIS_ONSET_FLOOR = 10

# Half-year period starts a DQS export can cover, used to fill the omitted
# zero-count rows back in before scanning for onset. Widen if a later export
# covers more years.
_NFLIS_FIRST_HALF = pd.Timestamp("2012-07-01")
_NFLIS_LAST_HALF = pd.Timestamp("2025-07-01")


def load_nflis_ca_counts(path: Path = NFLIS_CA_COUNTS_PATH) -> pd.DataFrame:
    """The cleaned NFLIS-Drug California export: one row per (substance,
    half-year), `drug_reports` summed across salt-form variants.

    Built once, by hand, from a DQS query export -- see
    `data/raw/nflis/NFLIS_Drug_DQS_*.csv` for the raw pull and the git history
    of this file for the aggregation script. `BASE_DESCRIPTION` in the raw
    export already merges salt forms (`XYLAZINE` + `XYLAZINE HYDROCHLORIDE` ->
    `Xylazine`); `Fluorofentanyl` was additionally merged into
    `para-Fluorofentanyl` by hand to match this project's own `ISOMER_ROLLUP`
    convention -- NFLIS treats the bare and `para-` forms as different
    catalog entries, but LA's ME office is known to have switched from
    writing one to the other mid-extract (see `docs/findings/trends.md`), so
    treating them as two substances here would undercount whichever name is
    out of fashion in a given period.
    """
    return pd.read_csv(path, parse_dates=["period_start"])


def _full_half_year_index() -> list[pd.Timestamp]:
    out = []
    for y in range(_NFLIS_FIRST_HALF.year, _NFLIS_LAST_HALF.year + 1):
        out += [pd.Timestamp(f"{y}-01-01"), pd.Timestamp(f"{y}-07-01")]
    return [h for h in out if _NFLIS_FIRST_HALF <= h <= _NFLIS_LAST_HALF]


def nflis_onset(series: pd.Series, floor: int = NFLIS_ONSET_FLOOR) -> pd.Timestamp | None:
    """First half-year >= `floor`, preceded by two consecutive half-years
    each < `floor`.

    `series` must already be reindexed onto a *complete* half-year calendar
    (`load_nflis_ca_counts` + a `reindex(_full_half_year_index(),
    fill_value=0)`, done by the caller) -- the raw export omits zero-count
    periods as rows entirely rather than writing them as 0, so scanning an
    un-reindexed series for "two consecutive periods" silently treats periods
    that are actually years apart as adjacent. This function does not
    reindex for you, deliberately: it should fail loudly (wrong answer, easy
    to spot against the printed series) rather than guess at a calendar.

    Does not attempt to reject a single-period spike that immediately
    collapses back below `floor` (mitragynine does exactly this at 2015H1,
    n=33 followed by n=2) -- read the result against the series, the same way
    `trends.py`'s `raw-count-rule` needed manual review on carfentanil and
    xylazine.
    """
    vals = series.tolist()
    for i in range(2, len(vals)):
        if vals[i] >= floor and vals[i - 1] < floor and vals[i - 2] < floor:
            return series.index[i]
    return None


@app.command("nflis-lag")
def nflis_lag(
    nflis_path: Path = typer.Option(NFLIS_CA_COUNTS_PATH, "--nflis-path"),
    ground_truth: Path = typer.Option(KNOWN_EMERGENCES_PATH, "--ground-truth"),
) -> None:
    """Compare each substance's California NFLIS onset to its LA ME-death
    onset -- the strongest external comparator this project has, since it is
    a real forensic count series rather than a notice date or a search
    snippet. Prints the full half-year series so the onset call can be
    checked by eye, not just trusted."""
    nflis = load_nflis_ca_counts(nflis_path)
    gt = pd.read_csv(ground_truth).set_index("substance")["interval_start"]
    halves = _full_half_year_index()

    for sub in sorted(nflis["substance"].unique()):
        s = (nflis[nflis["substance"] == sub]
             .set_index("period_start")["drug_reports"]
             .reindex(halves, fill_value=0))
        onset = nflis_onset(s)
        typer.echo(f"\n=== {sub} ===")
        typer.echo(s.to_string())
        typer.echo(f"CA NFLIS onset (>= {NFLIS_ONSET_FLOOR}, after 2 periods "
                   f"< {NFLIS_ONSET_FLOOR}): "
                   f"{onset.date() if onset is not None else 'none found'}")
        if sub in gt.index:
            typer.echo(f"LA ME-death onset (known_emergences.csv): {gt[sub]}")


def _fmt_intervals(ivs: list[tuple[pd.Timestamp, pd.Timestamp]]) -> str:
    if not ivs:
        return "none"
    return "; ".join(f"{s.year}Q{s.quarter}-{e.year}Q{e.quarter}" for s, e in ivs)


@app.command("score")
def score(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    ground_truth: Path = typer.Option(KNOWN_EMERGENCES_PATH, "--ground-truth"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
    sweep_path: Path = typer.Option(SWEEP_CACHE, "--sweep-path",
                                    help="cached eb05_sweep.csv from `trends "
                                         "alarms`; recomputed if absent"),
) -> None:
    """Score detected EB05 alarm episodes against the reference intervals."""
    counts, denom = load_quarterly(mentions, weighted=False)
    last_q = max(denom.index)

    if sweep_path.exists():
        sweep = pd.read_csv(sweep_path, parse_dates=["as_of"])
    else:
        typer.echo(f"no cached sweep at {sweep_path}; recomputing "
                   "(~1 min) -- run `trends alarms` first to cache it")
        sweep = sweep_eb05(counts, denom)

    gt = load_ground_truth(ground_truth, last_q)
    if "no_emergence" not in gt.columns:
        gt["no_emergence"] = False

    rows = []
    for sub, g in gt.groupby("substance", sort=False):
        # A documented negative case (e.g. Morphine: a structurally declining
        # substance whose one EB05 breach traces to a single-quarter blip) --
        # zero ground-truth-positive time, regardless of whatever placeholder
        # sits in interval_start/interval_end for that row (see the CSV's own
        # notes column for why a placeholder is there at all). Any detected
        # episode is then, correctly, 100% false alarm.
        gt_ivs = [] if bool(g["no_emergence"].iloc[0]) else list(
            zip(g["interval_start"], g["interval_end"]))
        det_ivs = _episode_spans(sweep, sub, threshold)
        m = interval_overlap(gt_ivs, det_ivs)
        rows.append({
            "substance": sub, **m,
            "ground_truth": _fmt_intervals(gt_ivs),
            "detected": _fmt_intervals(det_ivs),
        })

    out = pd.DataFrame(rows).set_index("substance")
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "overlap.csv")

    cols = ["iou", "recall", "precision", "gt_quarters", "det_quarters",
            "overlap_quarters", "ground_truth", "detected"]
    typer.echo(f"EB05 > {threshold:g}, scored against {ground_truth}\n")
    typer.echo(out[cols].to_string(float_format=lambda v: f"{v:.2f}"))
    ok = out.dropna(subset=["iou"])
    typer.echo(f"\nmean IoU {ok['iou'].mean():.3f}  "
               f"mean recall {ok['recall'].mean():.3f}  "
               f"mean precision {ok['precision'].mean():.3f}  "
               f"({len(ok)} substances scored)")

    pr = precision_rate(out.reset_index())
    typer.echo(f"substance-level false-positive rate: {pr['precision_rate_k']}/"
               f"{pr['precision_rate_n']} alarmed substances had a real hit "
               f"({pr['precision_rate']:.3f}, 95% CI "
               f"[{pr['precision_rate_lo']:.3f}, {pr['precision_rate_hi']:.3f}])")
    typer.echo(f"\nwrote {out_dir / 'overlap.csv'}")


if __name__ == "__main__":
    app()
