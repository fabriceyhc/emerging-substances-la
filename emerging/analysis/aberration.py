"""
Temporal aberration detectors: the "model the curve, not the ratio" family.

`trends.py`'s EB05 is squarely a **disproportionality** statistic -- shrink a
ratio of two period counts and threshold the shrunk value. See
`docs/LITERATURE_REVIEW.md` sec 1.5 for the parallel literature that instead
fits the shape of the quarterly series and flags an excessive *derivative*:
CDC's EARS aberration detectors (syndromic surveillance), the EudraVigilance
report-increase algorithm (van Holle & Bauchau, negative-binomial forecast
residual), and Whalen/Hauben/Bate's curve-disturbance detection. None of that
literature proposes replacing disproportionality with a slope test -- every
deployed example runs the two families *alongside* each other -- so these are
built as additional entries in `validation/benchmark.py`'s model table, not a
replacement for `eb05`.

Two candidates, both read off the same `(substance, quarter)` counts and
`quarter -> N` denominator `trends.load_quarterly` already produces, and both
swept over the identical recent/baseline window `sweep_eb05` uses, so a
side-by-side reading isn't also comparing different windows:

- `ears_sweep` -- EARS C2-style: the recent window's share of all OD deaths,
  read as a z-score against the *sample* mean and SD of that same share over
  the baseline window. Deliberately scored on **share**, not raw count, to
  match EB05 Gate A's own definition of "rising" -- the comparison this module
  exists for is "shrinkage vs. a moving-SD z-score", not also "share vs.
  count".
- `nb_trend_sweep` -- fit `log(mu_q) = a + b*t + log(N_q)` (a log-linear
  Poisson trend in *share*, offset by the quarter's own OD-death count) over
  the window, and score by the fitted slope `b`'s Wald z-score, dispersion-
  corrected the same quasi-likelihood way `trends._credit_dispersion` already
  corrects for overdispersion elsewhere in this project. This is the literal
  "model the curve, alarm on a too-steep derivative" detector -- the closest
  deployed precedent is the EudraVigilance algorithm, though that one flags a
  forecast *residual* rather than reading the slope directly.

Both functions return a long frame of `substance`, `as_of`, `n_recent` and the
detector's own score column, exactly the shape `validation/benchmark.py`'s
other sweeps produce, so they slot into `model_score`/`build_sweeps` as two
more statistics rather than a parallel pipeline.

## Porting `GPS_V2_DESIGN.md`'s three extensions

`#1` (position-weighted case definition) and `#2` (dual gate) port over
cleanly and are both implemented here (`weighted`/`role_discount` params;
`nb_trend_own`/`nb_trend_dual` in `validation/benchmark.py`) -- both are
changes to what counts as a death, upstream of whichever statistic is fit on
top, so they apply to a trend slope exactly as they do to a shrinkage ratio.

**`#3` (spatial concentration) does not port the same way, and this is worth
stating rather than skipping silently.** EB05's `"discount"` fusion shrinks
`E` by a constant factor for a spatially concentrated substance
(`E <- E * (1 - w*conc_s)`), which works because EB05 compares two *already-
collapsed* numbers (`n_recent` against one `expected`). `_poisson_trend_slope`
fits `log(mu_t) = a + b*t + log(offset_t)` -- if the same constant discount
were applied to `offset` uniformly across every quarter in the window (the
only way to reuse `concentration_weight`'s single per-(substance, as_of)
score as-is), it would multiply every `offset_t` by the same factor, and a
constant multiplicative rescaling of the offset folds entirely into the
intercept `a`: `log(offset_t * c) = log(offset_t) + log(c)`. **The fitted
slope `b`, and therefore the Wald z this detector actually alarms on, would
be provably unchanged.** Discounting the offset the naive way is not a
weaker version of the spatial term, it is a complete no-op on the statistic
that matters here -- confirmed algebraically, not assumed, before deciding
not to build it. A version that *would* have real effect would need to
discount only the recent-window quarters' offset relative to the baseline
window's (a level shift between the two blocks, closer in spirit to what the
EB05 discount does to the recent/baseline comparison) -- that is a genuinely
different mechanism from a constant rescale and its own design decision, not
a mechanical port of `#3`, and is not built here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import typer

from emerging.analysis.trends import EMERGENT_THRESHOLD, _episodes, load_quarterly
from emerging.config import BASELINE_QUARTERS, CUTOFF, RECENT_QUARTERS
from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import results_dir

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("aberration")

# EARS C1/C2's own literature line: alarm when the current observation is
# more than 3 baseline standard deviations above the moving mean. Kept as the
# fixed reading line for the `ears` model in `validation/benchmark.py`, the
# same role `ALARM_THRESHOLD` (1.5) plays for EB05 -- a different scale
# entirely (a z-score, not a posterior percentile), so it needs its own line,
# not the shared `--threshold`.
EARS_THRESHOLD = 3.0

# Wald z line for the trend-slope detector: the ordinary two-sided 95% line.
# Chosen, not fitted -- like `ALARM_THRESHOLD`, this is a reading convention
# to be checked against the backtest, not a claim that 1.96 is calibrated for
# ~200 simultaneously-scanned substances (no multiple-testing correction is
# applied here, matching EB05's own uncorrected per-substance line). Per-
# quarter, this line is well calibrated (checked directly: `docs/findings/
# synthetic.md` Finding 1, ~1.6% false-alarm rate on 200 genuinely flat
# synthetic substances against a nominal ~2.5%). Reading it across the
# ~34-quarter as-of sweep with no correction is a different story -- see
# `NB_TREND_CONFIRM_QUARTERS` below, which exists because of that finding.
NB_TREND_THRESHOLD = 1.96

# **Repeated-testing correction, found necessary by `docs/findings/
# synthetic.md` Finding 1, not assumed in advance.** `nb_trend_sweep` scores
# every substance at ~34 heavily overlapping as-of quarters (adjacent
# quarters' z-scores run ~0.7 lag-1 autocorrelated -- their windows share
# most of their data), and reading `nb_trend_z > NB_TREND_THRESHOLD` at any
# single quarter, with no correction, inflates the practical "does this
# substance ever alarm across its own sweep" false-positive rate from a
# nominal ~2.5% to ~33% on synthetic data with no real trend at all --
# `eb05`/`ears`, checked against the identical synthetic population, show
# 0%. Requiring `NB_TREND_CONFIRM_QUARTERS` *consecutive* alarming quarters
# before counting an alarm as confirmed brings the same population's rate
# back to ~2.5%, at negligible detection cost on the true positives checked
# there (49/50 vs 50/50 substances ever confirmed) -- 3 was chosen because
# it is the smallest window that actually restores calibration (2 only gets
# to ~15%, still ~6x nominal, because of how correlated adjacent quarters
# are), not because it was assumed. It also happens to match CDC EARS's own
# C3 variant (`docs/LITERATURE_REVIEW.md` sec 1.5), a three-point rule for
# an unrelated reason (a CUSUM-like running sum, not autocorrelation), which
# is a coincidence worth naming rather than a validation of the number --
# the two rules are not the same mechanism.
NB_TREND_CONFIRM_QUARTERS = 3

# Below this many total mentions in the fit window, a slope estimate is
# noise -- `nb_trend_sweep` still computes one (the Wald SE already reflects
# the extra uncertainty), but a substance never seen in the window at all
# gets no row rather than a spuriously confident z from an all-zero fit.
MIN_MENTIONS_FOR_TREND = 1


def _wide(counts: pd.DataFrame, quarters: list[pd.Timestamp]) -> pd.DataFrame:
    """substance x quarter count matrix, zero-filled, columns in `quarters`
    order -- the shared pivot both detectors below window-slice from."""
    wide = counts.pivot_table(index="substance", columns="quarter",
                              values="n", aggfunc="sum", fill_value=0)
    return wide.reindex(columns=quarters, fill_value=0)


def _sliced_counts(
    counts: pd.DataFrame, as_of: pd.Timestamp, weighted: bool,
    role_discount: bool, mentions_path: Path,
) -> pd.DataFrame:
    """The counts table to build this as-of's window from -- a plain slice of
    the already-loaded `counts` when unweighted, or a fresh
    `trends.load_quarterly(weighted=True, ...)` call when not, mirroring
    `trends.sweep_eb05`'s identical branch and for the identical reason:
    `_ever_independent`/`_role_dispersion` are cross-case aggregates, so
    slicing one upfront weighted table by `quarter <= as_of` would leak a
    substance's *later* independence/role evidence back into its earlier
    scores. Recomputing per as-of quarter is the cost of position-weighting
    at all, here as in `sweep_eb05` (~40s per full sweep).
    """
    if not weighted:
        return counts[counts["quarter"] <= as_of]
    sub_counts, _ = load_quarterly(
        mentions_path, cutoff=as_of + pd.offsets.QuarterEnd(0),
        weighted=True, quiet=True, role_discount=role_discount)
    return sub_counts


def ears_sweep(
    counts: pd.DataFrame, denom: pd.Series,
    recent_quarters: int = RECENT_QUARTERS,
    baseline_quarters: int = BASELINE_QUARTERS,
    weighted: bool = False, role_discount: bool = False,
    mentions_path: Path = MENTIONS_PATH,
) -> pd.DataFrame:
    """EARS C2-style z-score, swept over the same as-of grid `sweep_eb05` uses.

    For each as-of quarter: the baseline window's per-quarter *share* series
    (one value per baseline quarter) supplies a sample mean and SD; the recent
    window's aggregate share (summed counts over summed OD deaths, matching
    EB05's own `n_recent / n_recent_all`) is scored against it as a z.

    **The zero-baseline floor.** A substance with a flat, silent baseline has
    `sd_b == 0` by construction -- exactly the "genuinely novel substance"
    case this project cares most about, and the one place a bare z-score
    would divide by zero or blow up to an uninformative `inf`. Floored at the
    Poisson noise a single additional count would already imply at that
    substance's own baseline rate (`sqrt(max(mean_count, 0.25)) / mean_N`),
    so the score reads as "how many baseline-count-scale surprises is the
    recent window worth" rather than claiming more precision than a handful
    of quarterly counts can support -- the same motivation as the epsilon
    floor on `E` in `trends._nb_logpmf`, applied to a variance instead of a
    denominator.

    `weighted`/`role_discount` port `GPS_V2_DESIGN.md` #1 (and its extension)
    onto this detector: the position-weighted case definition changes what
    counts as a "mention" (the input), not the statistic built on top of it,
    so it applies here exactly as it does to `eb05` -- see `_sliced_counts`.
    """
    all_q = sorted(denom.index)
    total_q = recent_quarters + baseline_quarters
    if len(all_q) < total_q:
        raise ValueError(f"need {total_q} complete quarters, have {len(all_q)}")

    frames = []
    for i in range(total_q, len(all_q) + 1):
        as_of = all_q[i - 1]
        window = all_q[i - total_q:i]
        recent_q, base_q = window[-recent_quarters:], window[:-recent_quarters]

        sub_counts = _sliced_counts(counts, as_of, weighted, role_discount,
                                    mentions_path)
        wide = _wide(sub_counts, window)
        denom_w = denom.loc[window]

        share = wide.div(denom_w, axis=1)
        base_share = share[base_q]
        mean_b = base_share.mean(axis=1)
        sd_sample = base_share.std(axis=1, ddof=1).fillna(0.0)

        mean_b_count = wide[base_q].mean(axis=1)
        mean_denom_b = float(denom_w.loc[base_q].mean())
        sd_floor = np.sqrt(np.maximum(mean_b_count, 0.25)) / mean_denom_b
        sd_b = np.maximum(sd_sample, sd_floor)

        n_recent = wide[recent_q].sum(axis=1)
        n_recent_all = float(denom_w.loc[recent_q].sum())
        recent_share = n_recent / n_recent_all

        z = (recent_share - mean_b) / sd_b
        keep = n_recent > 0
        frames.append(pd.DataFrame({
            "substance": wide.index[keep],
            "as_of": as_of,
            "ears_z": z[keep].to_numpy(),
            "n_recent": n_recent[keep].to_numpy(),
        }))
    return pd.concat(frames, ignore_index=True)


def _poisson_trend_slope(
    obs: np.ndarray, offset: np.ndarray,
) -> tuple[float, float, float]:
    """IRLS fit of `log(mu_t) = a + b*t + log(offset_t)`; returns (b, se_b,
    dispersion). t is centred/scaled the same way `trends._poisson_trend_mu`
    centres it, so the two are reading the same kind of trend.

    Standard Poisson-GLM Newton-Raphson: the weighted-least-squares update at
    each step is exact (no line search needed) because the Poisson log-link
    negative log-likelihood is globally convex in `(a, b)`. `se_b` is the
    Wald standard error off the inverse Fisher information
    `(X^T diag(mu) X)^-1`, then quasi-likelihood corrected by the Pearson
    dispersion `sum((obs-mu)^2/mu) / df` -- the same overdispersion correction
    `trends._credit_dispersion`/`_dispersion_index` already use elsewhere in
    this project, applied here to a slope's standard error instead of to
    EB05's effective sample size. Dispersion is floored at 1 (never used to
    *shrink* the SE below the Poisson-implied one): a substance whose counts
    happen to be smoother than Poisson should not get a more confident slope
    than the model already claims.
    """
    n = len(obs)
    t = np.arange(n, dtype=float)
    t = (t - t.mean()) / max(t.std(), 1e-9)
    off = np.maximum(offset, 1e-9)
    X = np.column_stack([np.ones(n), t])

    beta = np.array([np.log(max(obs.sum(), 0.5) / off.sum()), 0.0])
    for _ in range(50):
        eta = X @ beta + np.log(off)
        mu = np.exp(np.minimum(eta, 50.0))
        z = (X @ beta) + (obs - mu) / np.maximum(mu, 1e-9)
        XtW = X.T * mu
        H = XtW @ X + 1e-9 * np.eye(2)
        beta_new = np.linalg.solve(H, XtW @ z)
        if np.max(np.abs(beta_new - beta)) < 1e-10:
            beta = beta_new
            break
        beta = beta_new

    eta = X @ beta + np.log(off)
    mu = np.exp(np.minimum(eta, 50.0))
    H = (X.T * mu) @ X + 1e-9 * np.eye(2)
    cov = np.linalg.inv(H)

    df = max(n - 2, 1)
    dispersion = max(float(np.sum((obs - mu) ** 2 / np.maximum(mu, 1e-9)) / df), 1.0)
    se_b = float(np.sqrt(cov[1, 1] * dispersion))
    return float(beta[1]), se_b, dispersion


def _emergent_flags(
    sweep: pd.DataFrame, score_col: str, threshold: float,
    baseline_quarters: int, emergent_threshold: float = EMERGENT_THRESHOLD,
) -> pd.Series:
    """Per-(substance, as_of) `emergent` flag -- the direct `nb_trend`
    analogue of `trends.alarm_history`'s `emergent` column, and built for the
    identical reason: `nb_trend_z` (and `ears_z`) have no notion of
    "already-established" at all, so a substance whose share has been
    genuinely, significantly rising for a *decade* (Methamphetamine -- see
    `docs/findings/benchmark.md`'s off-target investigation) scores exactly
    like a substance that appeared from nothing last quarter. `emergent`
    answers the separate question `alarm_history` asks for EB05: was this
    alarm's episode previously rare enough that the rise means something new?

    **Anchored at the current episode's onset, not read fresh every
    quarter -- ported from `alarm_history` for the identical reason it exists
    there.** A substance whose rise is real and sustained would otherwise
    "grow into" a higher rolling baseline as its own earlier rise ages into
    the baseline window, and flip from emergent to not partway through the
    very episode the flag exists to describe -- exactly the bug
    `alarm_history`'s docstring names and anchors against. Reading the
    baseline once, at the episode's first alarm quarter, and holding it for
    the episode's whole duration is that same fix, applied per-row here so it
    can gate a live per-quarter statistic (`nb_trend_emergent` in
    `validation/benchmark.py`) rather than only annotate a finished summary
    table the way `alarm_history` does for EB05.

    Requires an `n_baseline` column (`nb_trend_sweep` provides it). Rows
    outside any alarm episode get `False` -- `emergent` is an annotation on
    top of being in alarm, never a standalone signal, matching the EB05
    version exactly.
    """
    out = pd.Series(False, index=sweep.index)
    for _, g in sweep.groupby("substance"):
        g = g.sort_values("as_of")
        hot = g[score_col] > threshold
        if not hot.any():
            continue
        onset_baseline = g.set_index("as_of")["n_baseline"]
        for start, end in _episodes(list(g.loc[hot, "as_of"])):
            avg_baseline_q = float(onset_baseline.loc[start]) / baseline_quarters
            in_episode = (g["as_of"] >= start) & (g["as_of"] <= end) & hot
            out.loc[g.index[in_episode]] = avg_baseline_q < emergent_threshold
    return out


def _confirmed_flags(
    sweep: pd.DataFrame, score_col: str, threshold: float,
    min_consecutive: int = NB_TREND_CONFIRM_QUARTERS,
) -> pd.Series:
    """Per-(substance, as_of) `confirmed` flag: `True` only from the
    `min_consecutive`-th of a run of *consecutive* quarters each individually
    above `threshold` onward -- the repeated-testing correction
    `NB_TREND_CONFIRM_QUARTERS`'s own docstring motivates
    (`docs/findings/synthetic.md` Finding 1).

    **Reindexed to the sweep's full quarter grid before rolling, not rolled
    over adjacent *rows*.** A substance's own sweep is missing a row for any
    as-of quarter where `n_recent == 0` (`nb_trend_sweep`'s own filter) --
    rolling over the filtered rows directly would silently treat two
    alarming quarters either side of a real gap as consecutive. Reindexing
    to every quarter in the sweep and filling the gap with "not hot" first
    is what keeps "consecutive" meaning what it says, the same care
    `trends._episodes` already takes for `emergent`'s own episode logic,
    applied here via a plain rolling window instead of explicit interval
    arithmetic since "consecutive" is the only relationship this needs (no
    onset/duration bookkeeping).
    """
    all_q = sorted(sweep["as_of"].unique())
    out = pd.Series(False, index=sweep.index)
    for sub, g in sweep.groupby("substance"):
        hot = (g.set_index("as_of")[score_col] > threshold).reindex(
            all_q, fill_value=False)
        run = hot.astype(int).rolling(min_consecutive,
                                      min_periods=min_consecutive).sum()
        confirmed_q = set(run.index[run >= min_consecutive])
        out.loc[g.index[g["as_of"].isin(confirmed_q)]] = True
    return out


def nb_trend_sweep(
    counts: pd.DataFrame, denom: pd.Series,
    recent_quarters: int = RECENT_QUARTERS,
    baseline_quarters: int = BASELINE_QUARTERS,
    weighted: bool = False, role_discount: bool = False,
    mentions_path: Path = MENTIONS_PATH,
) -> pd.DataFrame:
    """Log-linear trend slope, Wald z-scored, swept over the same as-of grid
    `sweep_eb05` uses -- "model the curve, alarm on the derivative" read
    literally. See `_poisson_trend_slope` for the fit and
    `docs/LITERATURE_REVIEW.md` sec 1.5 for why a slope rather than a
    forecast residual: the EudraVigilance precedent flags a residual, but the
    slope is the more direct implementation of the question that motivated
    this module, and the two are monotonically related for a two-parameter
    linear trend fit to a fixed window.

    Fits **two** trends per substance per as-of, the same dual-gate move
    `GPS_V2_DESIGN.md` #2 makes for EB05 (`eb05`/`eb05_own`), and for the
    identical reason. `nb_trend_z` (`offset = N_q`, i.e. a trend in *share*)
    inherits EB05 Gate A's own failure mode: a substance whose raw count is
    flat or falling still posts share growth, and therefore a positive slope,
    for free whenever the county-wide total is falling (`docs/RESEARCH_LOG.md`
    / `rank_substances`'s cocaine example). Checked directly on this
    detector's own off-target list before this was added: Methamphetamine
    alarms 26 of the fixed line's 43 off-target alarm-quarters (`docs/
    findings/benchmark.md`), with its raw count in a clean, unbroken decline
    (432 -> 306 over the trailing 8 quarters) over the same window its
    *share* trend reads as significantly rising -- textbook denominator
    drift, not a real signal. `nb_trend_own_z` (`offset = 1`, a trend in the
    raw count itself, mirroring `expected_own`'s "no county-total term at
    all") answers the same question Gate B answers for EB05, and
    `nb_trend_dual` (`min` of the two, in `validation/benchmark.py`) is the
    direct analogue of `eb05_dual`.

    `weighted`/`role_discount` port `GPS_V2_DESIGN.md` #1 (and its extension)
    onto this detector, exactly as they do for `ears_sweep` -- see
    `_sliced_counts`. **No separate correction for the weighted credit
    scheme's overdispersion is applied here**, unlike `rank_substances`'s
    explicit `phi`/`phi_s` division: `_poisson_trend_slope`'s own
    quasi-likelihood dispersion correction is already fit fresh from each
    substance's own window residuals at every as-of step, so it already
    absorbs whatever extra variance the weighted scheme introduces, the same
    role `phi` plays for GPS, just estimated locally instead of pooled
    globally. The trade is power, not correctness: GPS's `phi` is one number
    fit across ~200 substances, while this trend's dispersion is refit from
    as few as `recent_quarters + baseline_quarters` residuals per substance
    per quarter -- noisier, but not double-corrected by construction.

    Also returns `n_baseline` (this window's raw count sum, `phi`-uncorrected
    since it is a diagnostic quantity here, not part of the fit),
    `emergent` (`_emergent_flags`, `nb_trend_z`'s own alarm episodes at
    `NB_TREND_THRESHOLD`) -- the EMERGENT_THRESHOLD analogue this detector
    was missing (see the module docstring and `docs/findings/benchmark.md`)
    -- and `confirmed` (`_confirmed_flags`, `NB_TREND_CONFIRM_QUARTERS`
    consecutive alarming quarters) -- the repeated-testing correction
    `docs/findings/synthetic.md` Finding 1 found necessary.
    """
    all_q = sorted(denom.index)
    total_q = recent_quarters + baseline_quarters
    if len(all_q) < total_q:
        raise ValueError(f"need {total_q} complete quarters, have {len(all_q)}")

    frames = []
    for i in range(total_q, len(all_q) + 1):
        as_of = all_q[i - 1]
        window = all_q[i - total_q:i]
        recent_q, base_q = window[-recent_quarters:], window[:-recent_quarters]

        sub_counts = _sliced_counts(counts, as_of, weighted, role_discount,
                                    mentions_path)
        wide = _wide(sub_counts, window)
        offset = denom.loc[window].to_numpy(dtype=float)
        offset_own = np.ones(len(window))
        n_recent = wide[recent_q].sum(axis=1)
        n_baseline = wide[base_q].sum(axis=1)
        fit_rows = wide.index[wide.sum(axis=1) >= MIN_MENTIONS_FOR_TREND]

        rows = []
        for sub in fit_rows:
            if n_recent.loc[sub] <= 0:
                continue
            obs = wide.loc[sub].to_numpy(dtype=float)
            b, se_b, dispersion = _poisson_trend_slope(obs, offset)
            b_own, se_b_own, _ = _poisson_trend_slope(obs, offset_own)
            rows.append({
                "substance": sub, "as_of": as_of,
                "nb_trend_z": b / se_b if se_b > 0 else 0.0,
                "nb_trend_slope": b, "nb_trend_dispersion": dispersion,
                "nb_trend_own_z": b_own / se_b_own if se_b_own > 0 else 0.0,
                "nb_trend_own_slope": b_own,
                "n_recent": int(n_recent.loc[sub]),
                "n_baseline": int(n_baseline.loc[sub]),
            })
        frames.append(pd.DataFrame(rows))
    out = pd.concat(frames, ignore_index=True)
    out["emergent"] = _emergent_flags(out, "nb_trend_z", NB_TREND_THRESHOLD,
                                      baseline_quarters)
    out["confirmed"] = _confirmed_flags(out, "nb_trend_z", NB_TREND_THRESHOLD)
    return out


@app.command("rank")
def rank(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    top: int = typer.Option(25, "--top"),
    min_recent: int = typer.Option(2, "--min-recent",
                                   help="Suppress substances below this recent count"),
    weighted: bool = typer.Option(
        False, "--weighted/--no-weighted",
        help="Position-weighted case definition (GPS_V2_DESIGN.md #1). Off "
             "by default here, unlike `trends rank` -- the benchmark found a "
             "real but mixed effect for this detector family (best mean IoU "
             "in the whole comparison, but it *adds* off-target substances "
             "rather than removing them; see docs/findings/benchmark.md)"),
    role_discount: bool = typer.Option(False, "--role-discount/--no-role-discount"),
) -> None:
    """Rank substances by the model-the-curve family: nb-trend's Poisson
    trend-slope Wald z and EARS's moving-baseline z, read at the current
    quarter, with the `emergent` annotation alongside them.

    Complementary to `trends rank`'s EB05 ranking, not a replacement for it
    -- see `docs/LITERATURE_REVIEW.md` sec 1.5 and `docs/findings/
    benchmark.md` for the full comparison this command's defaults come from.

    **`emergent` is reported as a column, never used to suppress a row.**
    `_emergent_flags` answers "was this alarm episode previously rare" --
    gating on it throws out real detections (Fentanyl, PCP) along with the
    false one it was built to catch (Methamphetamine), the same trade
    `trends.alarm_history` avoids by keeping its own `emergent` column an
    annotation rather than a filter. A substance alarming with
    `emergent=False` is not a false alarm to dismiss; it reads differently
    from one with `emergent=True` to a reviewer, which is the entire point
    of carrying the column at all.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff, quiet=True)
    nb = nb_trend_sweep(counts, denom, recent_quarters=recent_quarters,
                        baseline_quarters=baseline_quarters, weighted=weighted,
                        role_discount=role_discount, mentions_path=mentions)
    ea = ears_sweep(counts, denom, recent_quarters=recent_quarters,
                    baseline_quarters=baseline_quarters, weighted=weighted,
                    role_discount=role_discount, mentions_path=mentions)

    as_of = nb["as_of"].max()
    cur = nb[nb["as_of"] == as_of].set_index("substance")
    cur_ea = ea[ea["as_of"] == as_of].set_index("substance")["ears_z"]
    res = cur.join(cur_ea, how="outer")

    for col in ["nb_trend_z", "nb_trend_slope", "nb_trend_own_z", "ears_z"]:
        res[col] = res[col].fillna(0.0)
    res["n_recent"] = res["n_recent"].fillna(0).astype(int)
    res["n_baseline"] = res["n_baseline"].fillna(0).astype(int)
    res["emergent"] = res["emergent"].fillna(False)
    res["nb_trend_alarm"] = res["nb_trend_z"] > NB_TREND_THRESHOLD
    res["ears_alarm"] = res["ears_z"] > EARS_THRESHOLD
    res = res.sort_values("nb_trend_z", ascending=False)

    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "aberration_rank.csv")

    typer.echo(f"as-of {as_of.year}Q{as_of.quarter} "
               f"({recent_quarters}q recent vs {baseline_quarters}q baseline)"
               + (" [weighted+role]" if weighted and role_discount else
                  " [weighted]" if weighted else ""))
    typer.echo(f"nb-trend alarm line: z > {NB_TREND_THRESHOLD:g}  "
               f"ears alarm line: z > {EARS_THRESHOLD:g}")

    show = res[res["n_recent"] >= min_recent].head(top)
    cols = ["n_recent", "n_baseline", "nb_trend_z", "nb_trend_slope",
            "nb_trend_own_z", "nb_trend_alarm", "ears_z", "ears_alarm",
            "emergent"]
    typer.echo("\n" + show[cols].to_string(float_format=lambda v: f"{v:.2f}"))

    established = res[res["nb_trend_alarm"] & ~res["emergent"]]
    if len(established):
        typer.echo(f"\n{len(established)} substance(s) alarm on nb-trend but "
                   "read emergent=False -- a real, credible rise in a "
                   "substance that was not previously rare (see "
                   "docs/findings/benchmark.md's Methamphetamine finding for "
                   "why this is reported, not suppressed): "
                   + ", ".join(established.sort_values(
                       "nb_trend_z", ascending=False).head(6).index))
    typer.echo(f"\nwrote {out_dir / 'aberration_rank.csv'}")
