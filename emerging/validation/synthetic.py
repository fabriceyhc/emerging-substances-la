"""
Synthetic benchmark generator for the temporal detectors -- the EB05 family
(`emerging.analysis.trends`) and the model-the-curve family
(`emerging.analysis.aberration`).

See `docs/SYNTHETIC_BENCHMARK_DESIGN.md` for the full design, in particular
the circularity risk this module is built to avoid: a generator whose "in
alarm" cases share a detector's own generative assumptions (Poisson counts
on a log-linear ramp is exactly `nb_trend`'s model) would only prove that
detector can recover its own model, not that it detects real emergence.
`SHAPES`/`NOISE_MODELS` below exist because of that -- `"log_linear"` +
`"poisson"` is the sanity-check scenario, not the only one a fair comparison
can use.

**Phase 1 only** (design note sec 5): every `Profile` here produces
aggregate `(substance, quarter) -> count` and `quarter -> N`, exactly what
`trends.sweep_eb05`/`aberration.ears_sweep`/`aberration.nb_trend_sweep`
already take. No case-level records, so `--weighted`/`--role-discount`
cannot be exercised against synthetic data yet -- that needs phase 2's
case-level generator (design note sec 5), not built here.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

from emerging.analysis.aberration import NB_TREND_THRESHOLD, nb_trend_sweep
from emerging.analysis.trends import load_quarterly
from emerging.config import RECENT_QUARTERS
from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import results_dir
from emerging.viz import AQUA, BLUE, INK, INK_2, MUTED, ORANGE, YELLOW

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("synthetic")

# Reused only as the *default* label boundary when a profile's own
# `emergent` property is read -- never hardcoded into a detector's own
# reading. Design note sec 1: the boundary itself should be swept, not
# fixed at the production value, when testing whether a detector's own
# emergent-flag threshold generalizes rather than just echoing this one.
DEFAULT_EMERGENT_THRESHOLD = 2.0

SHAPES = ("step", "log_linear", "sigmoid", "batch")
NOISE_MODELS = ("poisson", "negbin", "clustered")
DENOM_MODES = ("flat", "declining", "rising")


@dataclass(frozen=True)
class Profile:
    """One synthetic substance's generative recipe (design note sec 1-3).

    `rises`/`sustained`/`baseline_rate` are what fix the ground-truth label
    (see `Profile.emergent` and `generate_panel`'s ground-truth rows) --
    everything else (`shape`, `noise`, `dispersion`) varies which detector
    assumption a given draw exercises, per the module docstring's
    circularity warning.
    """

    name: str
    baseline_rate: float           # avg count/quarter before any rise
    rises: bool = False
    sustained: bool = True         # False = reverts to baseline (a "blip")
    shape: str = "log_linear"      # SHAPES
    magnitude: float = 3.0         # fold-increase at the new plateau
    onset: int = 20                # quarter index the rise starts
    duration: int = 4              # quarters to reach the new level
    noise: str = "poisson"         # NOISE_MODELS
    dispersion: float = 1.0        # NB dispersion phi (Var = phi*mean); noise="negbin" only
    n_bursts: int = 1              # shape="batch" only -- see _rate_curve
    perturb_at: int | None = None  # quarter index of a brief opposite-direction shock
    perturb_factor: float = 1.0    # multiplies rate(t) there; <1 = dip, >1 = spike
    perturb_every: int | None = None  # repeat the shock every N quarters (None = once)

    def __post_init__(self) -> None:
        if self.shape not in SHAPES:
            raise ValueError(f"shape must be one of {SHAPES}, got {self.shape!r}")
        if self.noise not in NOISE_MODELS:
            raise ValueError(f"noise must be one of {NOISE_MODELS}, got {self.noise!r}")

    @property
    def real_positive(self) -> bool:
        """Whether this profile is a genuine ground-truth-positive
        sustained rise: `rises and sustained`, with one carve-out --
        `shape="batch"` can never actually stay elevated. `_rate_curve`'s
        `"batch"` branch concentrates the excess into 1-2 quarters right at
        onset and is flat everywhere else *regardless of `sustained`*
        (`test_batch_shape_concentrates_excess_rather_than_spreading_a_rate`
        pins this down as intentional -- it models a one-off event like
        para-fluorofentanyl's batch arrivals, not a rate change). Before
        this carve-out, `Profile("x", 1.0, rises=True, sustained=True,
        shape="batch", ...)` and the identical profile with `sustained=False`
        produced *bit-identical* rate curves and sampled data, yet opposite
        ground truth (`no_emergence=False` vs `True`, an interval running to
        the series end vs a reverted blip) -- caught by
        `docs/findings/synthetic.md`'s profile-gallery figure, where the two
        were visually indistinguishable despite the opposite label.
        `sustained=True` on a batch profile no longer claims an ongoing
        elevation the data itself does not show. Applies identically
        whether the batch recurs (`n_bursts > 1`) or not -- repeated,
        disconnected events are still not a trend.

        **Second carve-out: `magnitude <= 1.0` is a decline, not a rise, and
        is never a positive.** Every shape formula here is a generic
        multiplicative curve from `baseline_rate` to `baseline_rate *
        magnitude` -- nothing stops `magnitude < 1` from being used to build
        a genuinely declining series (a substance fading out of use), which
        is a real and useful adversarial case (does a detector correctly
        stay quiet on a fading substance, including one with a misleading
        `perturb_at` spike layered on top?) but is definitionally the
        opposite of an emergence. Without this carve-out, `rises=True` on a
        declining profile would have claimed a real, ongoing positive
        exactly the way the unfixed batch carve-out above did.
        """
        return (self.rises and self.sustained and self.shape != "batch"
                and self.magnitude > 1.0)

    @property
    def emergent(self) -> bool | None:
        """`None` when the label doesn't apply (not a real, sustained rise
        -- see `real_positive`); otherwise whether `baseline_rate` reads as
        "was rare" against `DEFAULT_EMERGENT_THRESHOLD` -- the ground-truth
        analogue of `aberration._emergent_flags`'s own reading, computed
        here from the true generative rate rather than an estimated one.
        """
        if not self.real_positive:
            return None
        return self.baseline_rate < DEFAULT_EMERGENT_THRESHOLD


def _rate_curve(profile: Profile, T: int) -> np.ndarray:
    """The true, noise-free `rate(t)` for `t = 0..T-1`, in count/quarter.

    `"batch"` is the one shape that is not a rate change at all -- it adds a
    fixed excess count concentrated in the one or two quarters after onset,
    modeled directly on the real precedent `trends.py` documents for
    para-fluorofentanyl ("arrives in batches rather than at a rate"), and is
    the shape most likely to violate `_poisson_trend_slope`'s smooth-curve
    assumption. Every other shape is a genuine `rate(t)` and gets sampled
    through the same noise model as any other quarter.

    **`n_bursts` is what actually earns "batch" its own class, separate
    from a reverted blip.** A single burst (`n_bursts=1`, the default) is a
    lone spike over an otherwise flat series -- structurally the same shape
    a `sustained=False` blip already produces, just narrower, and not a
    strong enough reason on its own to justify a fourth `SHAPES` entry (see
    `docs/findings/synthetic.md`'s note on this). `n_bursts > 1` repeats the
    burst every `(T - onset) // n_bursts` quarters, independent, same size
    each time -- the actual real-world shape the para-fluorofentanyl
    precedent describes (recurring contaminated batches from a substance
    that never establishes a smooth trend), and the thing a naive
    per-quarter alarm can be fooled by repeatedly while `_confirmed_flags`
    (`NB_TREND_CONFIRM_QUARTERS` *consecutive* quarters) should not be,
    since the bursts are never adjacent by construction.

    **`perturb_at`/`perturb_factor` layer a brief, opposite-direction shock
    on top of whatever the base shape already produced** -- a genuine
    decline with one misleading upward spike partway through, or a genuine
    rise with one misleading downward dip. Applied last, after the
    `sustained` revert, so it modulates the *actual* rate at that quarter
    (elevated if inside a rise, back-to-baseline if reverted) rather than
    the shape's original, unperturbed value. Ground truth (`real_positive`)
    ignores it entirely -- one anomalous quarter should not flip whether an
    otherwise-real trend counts as an emergence, in either direction.

    `perturb_every` repeats the same shock periodically instead of once --
    layered on top of a genuine, sustained rise, this is the same
    "established substance whose deaths arrive in recurring batches rather
    than at a smooth rate" precedent `n_bursts` models for a purely flat
    series (`trends.py`'s para-fluorofentanyl note), just riding on top of a
    real trend instead of a flat one this time.
    """
    t = np.arange(T, dtype=float)
    rate = np.full(T, profile.baseline_rate, dtype=float)
    if not profile.rises:
        return rate

    onset, peak = profile.onset, profile.baseline_rate * profile.magnitude
    if profile.shape == "step":
        rate = np.where(t >= onset, peak, profile.baseline_rate)
    elif profile.shape == "log_linear":
        b = (np.log(max(peak, 1e-9) / max(profile.baseline_rate, 1e-9))
             / max(profile.duration, 1))
        ramp = profile.baseline_rate * np.exp(
            b * np.clip(t - onset, 0, profile.duration))
        rate = np.where(t >= onset, ramp, profile.baseline_rate)
    elif profile.shape == "sigmoid":
        center = onset + profile.duration / 2.0
        width = max(profile.duration / 6.0, 0.5)
        frac = np.where(t >= onset, 1.0 / (1.0 + np.exp(-(t - center) / width)), 0.0)
        rate = profile.baseline_rate + (peak - profile.baseline_rate) * frac
    elif profile.shape == "batch":
        rate = np.full(T, profile.baseline_rate, dtype=float)
        span = max(min(2, T - onset), 1)
        extra_per_burst = profile.baseline_rate * (profile.magnitude - 1.0) * profile.duration
        n_bursts = max(profile.n_bursts, 1)
        gap = max((T - onset) // n_bursts, span) if n_bursts > 1 else 0
        for burst in range(n_bursts):
            start = onset + burst * gap
            for i in range(span):
                if 0 <= start + i < T:
                    rate[start + i] += extra_per_burst / span

    if not profile.sustained:
        revert_at = onset + profile.duration
        rate = np.where(t >= revert_at, profile.baseline_rate, rate)

    if profile.perturb_at is not None:
        period = profile.perturb_every or T  # None -> one shock, never repeats
        for shock_at in range(profile.perturb_at, T, period):
            span = max(min(2, T - shock_at), 1)
            for i in range(span):
                idx = shock_at + i
                if 0 <= idx < T:
                    rate[idx] *= profile.perturb_factor
    return np.maximum(rate, 0.0)


def _sample_counts(
    rate: np.ndarray, noise: str, dispersion: float, rng: np.random.Generator,
) -> np.ndarray:
    """Draw one noisy count series from a true `rate(t)` curve.

    `"negbin"` is parametrized by (mean, dispersion) rather than numpy's own
    (n, p), so `dispersion` reads on the same Var/mean scale
    `trends._credit_dispersion`'s `phi` already uses elsewhere in this
    project -- `dispersion=1` collapses to Poisson exactly.

    `"clustered"` is a light self-exciting process: each quarter's Poisson-
    drawn "parents" independently spawn extra correlated arrivals the
    following quarter, so the series is autocorrelated without being a
    `rate(t)` change -- can be layered onto a flat, non-rising rate, unlike
    the `"batch"` *shape*, which is definitionally a rise.
    """
    if noise == "poisson":
        return rng.poisson(np.maximum(rate, 0.0))
    if noise == "negbin":
        phi = max(dispersion, 1.0 + 1e-6)
        r = np.maximum(rate, 1e-9)
        size = r / (phi - 1.0)
        p = 1.0 / phi
        return rng.negative_binomial(np.maximum(size, 1e-6), p)
    if noise == "clustered":
        # parent_frac + branching*parent_frac + direct_frac must sum to 1 so
        # E[total] == rate; solved for the normalizer rather than hand-tuned,
        # so it stays correct under any (parent_frac, branch) choice. Caught
        # by its own calibration test: the first cut (0.6/0.4 split, no
        # normalizer) sampled at 1.42x the intended rate, since a spawned
        # parent's own children are extra expectation the naive split never
        # accounted for.
        parent_frac, branch = 0.6, 0.7
        norm = parent_frac * (1.0 + branch) + (1.0 - parent_frac)
        r = np.maximum(rate, 0.0) / norm
        parents = rng.poisson(r * parent_frac)
        extra = np.zeros_like(parents)
        for t in range(len(parents) - 1):
            if parents[t] > 0:
                extra[t + 1] += rng.poisson(parents[t] * branch)
        return parents + extra + rng.poisson(r * (1.0 - parent_frac))
    raise ValueError(f"unknown noise {noise!r}")


def _denom_curve(
    T: int, mode: str, n0: float, rng: np.random.Generator,
) -> np.ndarray:
    """The synthetic county-wide total `N_t`, an independent axis from any
    one substance's own rate (design note sec 4). `"declining"`'s -0.5%/
    quarter compounds to ~20% over a 45-quarter window (`(1-0.005)**44 ~=
    0.80`), matching the real ~18-22% decline `regime()` documents --
    calibrated to that number and not to a rounder-looking rate: an earlier
    cut used -3%/quarter, which compounds to a 74% decline over the same
    window and produced an artificially extreme adversarial scenario (every
    flat substance's share exploding, not a realistic stress test of the
    denominator-drift failure mode). Crossing a flat-own-rate profile
    against this mode is what manufactures a Methamphetamine-shaped
    denominator-drift negative on demand -- see `generate_panel`.
    """
    t = np.arange(T, dtype=float)
    if mode == "flat":
        mu = np.full(T, n0)
    elif mode == "declining":
        mu = n0 * np.exp(-0.005 * t)
    elif mode == "rising":
        mu = n0 * np.exp(0.005 * t)  # symmetric with "declining"
    else:
        raise ValueError(f"mode must be one of {DENOM_MODES}, got {mode!r}")
    return np.maximum(rng.poisson(mu).astype(float), 1.0)


def generate_panel(
    profiles: list[Profile], T: int = 45, denom_mode: str = "flat",
    denom_n0: float = 400.0, seed: int = 0, start: str = "2014-01-01",
) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """One Monte Carlo replicate: `(counts, denom, ground_truth)` in exactly
    the shapes `sweep_eb05`/`ears_sweep`/`nb_trend_sweep` and
    `ground_truth.load_ground_truth` already consume, so every existing
    `Model`/`score_model` path in `validation/benchmark.py` runs against
    this unmodified once `build_sweeps` has the injection point design note
    sec 6 describes.

    Every profile gets a ground-truth row: a sustained rise gets a real
    `(interval_start, interval_end)` and `no_emergence=False`; anything else
    (flat, or a reverted blip) gets `no_emergence=True` with a placeholder
    interval at the series midpoint, mirroring the real
    `known_emergences.csv` convention documented in its own negative rows
    (`no_emergence=True` rows carry a placeholder interval "purely so this
    row parses"; `score_model` ignores it). `emergent` is the one column
    real ground truth cannot supply at all -- read here from the true
    generative `baseline_rate`, not an estimate.
    """
    rng = np.random.default_rng(seed)
    quarters = pd.date_range(start, periods=T, freq="QS")
    denom = pd.Series(_denom_curve(T, denom_mode, denom_n0, rng),
                      index=quarters, name="N")
    mid_q = quarters[T // 2]

    rows, gt_rows = [], []
    for p in profiles:
        rate = _rate_curve(p, T)
        counts = _sample_counts(rate, p.noise, p.dispersion, rng)
        for q, n in zip(quarters, counts):
            if n > 0:
                rows.append({"substance": p.name, "quarter": q, "n": int(n)})

        real_positive = p.real_positive
        onset_q = quarters[p.onset] if real_positive else mid_q
        end_q = quarters[-1] if real_positive else mid_q
        gt_rows.append({
            "substance": p.name,
            "interval_start": f"{onset_q.year}Q{onset_q.quarter}",
            "interval_end": f"{end_q.year}Q{end_q.quarter}",
            "ongoing": False,
            "no_emergence": not real_positive,
            "emergent": bool(p.emergent) if real_positive else False,
        })

    counts_df = pd.DataFrame(rows, columns=["substance", "quarter", "n"])
    counts_df.attrs["phi"] = 1.0
    counts_df.attrs["phi_s"] = {}
    gt_df = pd.DataFrame(gt_rows)
    return counts_df, denom, gt_df


# --------------------------------------------------------------------------
# Plots. Purely diagnostic -- these render this module's own generator, never
# real data, so a reader can see what a `Profile` actually produces (the
# design doc's taxonomy made visible) and what the specific synthetic
# populations behind `docs/findings/synthetic.md` Findings 1/2 looked like
# when nb-trend false-alarmed on them.
# --------------------------------------------------------------------------

# One profile per corner of the design space (design note sec 1-3): three
# ground-truth classes (emergent riser, established riser, blip/batch)
# crossed with the four ramp `SHAPES`, plus the three `NOISE_MODELS` --
# including two shown riding a real rise, not only a flat series, since a
# detector's dispersion correction has to hold up in both places. Not
# exhaustive -- chosen to be the smallest set where every axis of the
# taxonomy appears at least once, plus one deliberate near-duplicate pair
# (`single_batch`/`recurring_batch`) kept *because* they are similar: the
# similarity is the point, see the pair's own comment below.
GALLERY_PROFILES = [
    Profile("rare_riser", 1.0, rises=True, sustained=True, shape="log_linear",
            magnitude=8.0, onset=20, duration=4),
    Profile("established_riser", 15.0, rises=True, sustained=True,
            shape="log_linear", magnitude=3.0, onset=20, duration=4),
    Profile("step_riser", 1.5, rises=True, sustained=True, shape="step",
            magnitude=5.0, onset=22),
    Profile("sigmoid_riser", 1.2, rises=True, sustained=True, shape="sigmoid",
            magnitude=7.0, onset=15, duration=8),
    # A single burst (n_bursts=1, "batch"'s default) is just a narrow blip
    # by another name -- not shown separately here for that reason (see the
    # module docstring's note on this). n_bursts>1 is what actually earns
    # "batch" its own class: several independent, disconnected events,
    # visibly different from one blip.
    Profile("recurring_batch", 1.0, rises=True, sustained=True, shape="batch",
            magnitude=6.0, onset=4, duration=2, n_bursts=4),
    Profile("reverting_blip", 3.0, rises=True, sustained=False, shape="step",
            magnitude=5.0, onset=18, duration=3),
    # A real, sustained decline (magnitude<1 -- Profile.real_positive's own
    # carve-out) with one misleading upward spike partway through: does a
    # detector correctly read this as "not an emergence" despite the spike?
    Profile("declining_with_blip", 15.0, rises=True, sustained=True,
            shape="log_linear", magnitude=0.2, onset=6, duration=8,
            perturb_at=26, perturb_factor=5.0),
    # The mirror image: a real, sustained rise with one misleading downward
    # dip partway through -- still a genuine emergence, ground truth ignores
    # the dip (`_rate_curve`'s own docstring on `perturb_at`).
    Profile("rising_with_dip", 1.2, rises=True, sustained=True,
            shape="log_linear", magnitude=9.0, onset=10, duration=6,
            perturb_at=28, perturb_factor=0.15),
    Profile("flat_poisson", 8.0, noise="poisson"),
    # Overdispersion visibly reads on a line plot (occasional much taller
    # spikes than Poisson allows at the same mean) -- unlike "clustered"
    # (autocorrelation), which does not: a clustered series and a plain
    # Poisson one look like the same kind of wiggle to the eye, so it is
    # deliberately not shown here even though the generator supports it and
    # `test_clustered_noise_is_actually_autocorrelated` checks it
    # numerically. A picture earns its place by being visually diagnosable,
    # not just statistically different.
    Profile("flat_negbin", 8.0, noise="negbin", dispersion=4.0),
    Profile("rare_riser_negbin", 1.0, rises=True, sustained=True,
            shape="log_linear", magnitude=8.0, onset=20, duration=4,
            noise="negbin", dispersion=4.0),
    # Closes the loop this whole gallery is built around: an established,
    # genuinely elevated substance (real_positive=True) whose deaths also
    # arrive in periodic batches on top of the real trend -- the actual
    # precedent `trends.py` documents for para-fluorofentanyl, modeled here
    # for the first time as a positive rather than only as a negative
    # (recurring_batch, above, has no underlying trend at all).
    Profile("established_riser_batchy", 15.0, rises=True, sustained=True,
            shape="log_linear", magnitude=3.0, onset=10, duration=10,
            perturb_at=6, perturb_every=7, perturb_factor=2.2),
]


def _caption(*paragraphs: str, width: int = 95) -> str:
    """Hard-wrap each paragraph at `width` chars before handing it to
    `fig.text` -- `fig.text` never wraps on its own, so one long unbroken
    sentence forces `savefig(bbox_inches="tight")` to widen the whole canvas
    to fit it, leaving the actual panel grid stranded in the left portion of
    a much wider image with dead space to its right. `width` should scale
    with the figure's own width (`_caption_width`) -- a fixed constant either
    overflows a narrow figure or under-fills a wide one, wrapping into more
    short lines than the available space needs."""
    return "\n\n".join(textwrap.fill(p, width=width) for p in paragraphs)


def _caption_width(ncol: int, per_col_in: float = 3.3, chars_per_in: float = 13.0) -> int:
    """Chars-per-line that roughly fills a `ncol`-panel-wide figure at the
    caption's own font size, calibrated once against this module's own
    figures rather than assumed -- see `_caption`."""
    return max(int(chars_per_in * per_col_in * ncol), 60)


def _badge_for(p: Profile) -> str:
    """Mirrors `Profile.real_positive`/`generate_panel`'s own ground-truth
    logic exactly, rather than re-deriving it from `rises`/`sustained`
    directly -- the batch/blip mislabeling `real_positive`'s own docstring
    describes was caught by this badge disagreeing with what the line
    visibly does, so this function has to read the same source of truth
    the fix changed, not a copy of the old condition."""
    if p.real_positive:
        return "EMERGENT" if p.emergent else "ESTABLISHED"
    if p.rises and p.magnitude <= 1.0:
        return "DECLINING (not an emergence)"
    if p.rises and p.shape == "batch":
        return (f"BATCH x{p.n_bursts} (recurring, never sustained)"
                if p.n_bursts > 1 else "BATCH x1 (= a narrow blip)")
    if p.rises and not p.sustained:
        return "BLIP (reverts)"
    return "FLAT"


def _badge_color(badge: str) -> tuple[str, str]:
    """Color by badge *kind*, not the exact string -- several badges
    (`BATCH x{n}...`, `STILL ALARMS\\n...`) carry data in their text, so a
    dict keyed on the literal string would silently fall through to the
    default color for any value it hadn't seen before."""
    if badge.startswith("EMERGENT"):
        return ORANGE, "white"
    if badge.startswith("ESTABLISHED"):
        return BLUE, "white"
    if badge.startswith(("BLIP", "BATCH", "SUPPRESSED", "DECLINING")):
        return AQUA, "white"
    if badge.startswith(("NAIVE ALARM", "STILL ALARMS")):
        return YELLOW, INK
    return MUTED, INK


def _panel(
    ax, quarters: list[pd.Timestamp], counts: np.ndarray, rate: np.ndarray,
    title: str, badge: str | None = None, hot_spans: list[tuple] | None = None,
    confirmed_spans: list[tuple] | None = None,
) -> None:
    """One small-multiple: sampled counts (solid) against the true noise-free
    `rate(t)` that generated them (dashed) -- the one comparison a real
    substance's history can never offer, since its true rate is exactly what
    every detector in this project is trying to estimate.

    `hot_spans` shades the recent-quarters window(s) where the trend test
    naively exceeded its alarm line -- a light band, since one look is a
    weak claim. `confirmed_spans` is the subset that also cleared three
    consecutive alarming quarters, drawn as a darker overlay on top of the
    same band, so a reader can see at a glance which alarms were fleeting
    and which held up.
    """
    ax.plot(quarters, counts, color=BLUE, linewidth=1.5, marker="o",
            markersize=2.2, zorder=3)
    ax.plot(quarters, rate, color=INK_2, linewidth=1.1, linestyle="--",
            alpha=0.85, zorder=2)
    for s, e in (hot_spans or []):
        ax.axvspan(s, e, color=YELLOW, alpha=0.28, linewidth=0, zorder=1)
    for s, e in (confirmed_spans or []):
        ax.axvspan(s, e, color=ORANGE, alpha=0.35, linewidth=0, zorder=1)
    ax.set_title(title, fontsize=8.3, color=INK, loc="left")
    if badge:
        color, text_color = _badge_color(badge)
        ax.text(0.98, 0.94, badge, transform=ax.transAxes, fontsize=6.6,
                color=text_color, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.25", facecolor=color,
                          edgecolor="none"))
    ax.tick_params(labelsize=6.3, colors=INK_2)
    ax.set_xlim(quarters[0], quarters[-1])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _plot_gallery(out_dir: Path, seed: int) -> None:
    T = 40
    counts_df, denom, _ = generate_panel(GALLERY_PROFILES, T=T, seed=seed)
    quarters = sorted(denom.index)
    sweep = nb_trend_sweep(counts_df, denom)
    ncol = 4
    nrow = int(np.ceil(len(GALLERY_PROFILES) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.15 * nrow),
                             sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, p in zip(axes, GALLERY_PROFILES):
        rate = _rate_curve(p, T)
        s = (counts_df[counts_df["substance"] == p.name]
             .set_index("quarter")["n"].reindex(quarters, fill_value=0))
        title = f"{p.name}  ({p.shape}, {p.noise})"
        g = sweep[sweep["substance"] == p.name]
        hot_q = g.loc[g["nb_trend_z"] > NB_TREND_THRESHOLD, "as_of"].tolist()
        confirmed_q = g.loc[g["confirmed"], "as_of"].tolist()
        _panel(ax, quarters, s.to_numpy(), rate, title, badge=_badge_for(p),
              hot_spans=_hot_spans(hot_q, quarters),
              confirmed_spans=_hot_spans(confirmed_q, quarters))
    for ax in axes[len(GALLERY_PROFILES):]:
        ax.set_visible(False)
    fig.suptitle("What a synthetic Profile looks like: the design space",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.text(0.01, 0.005, _caption(
             "Solid = sampled quarterly count; dashed = the true "
             "noise-free rate that generated it -- never available for a "
             "real substance, only for a synthetic one. Badge = the "
             "ground-truth label assigned from the generative recipe, "
             "not read off the noisy line. Light yellow = a quarter the "
             "trend test actually alarmed on; darker orange = it also "
             "held for three quarters in a row, the bar an actually-"
             "confirmed detection needs to clear -- these are what the "
             "detector would have flagged, not the ground truth.",
             "declining_with_blip and rising_with_dip each carry one "
             "quarter that goes the opposite way from the rest of the "
             "line -- a real decline with a misleading upward spike, and "
             "a real rise with a misleading downward dip. The label "
             "follows the overall trend, not the one odd quarter: the "
             "decline stays negative despite its spike, the rise stays "
             "positive despite its dip. established_riser_batchy is a "
             "real, ongoing rise whose deaths also arrive in periodic "
             "clumps rather than smoothly -- still labeled a positive "
             "throughout, the clumps are noise riding on a real signal, "
             "not the signal itself.",
             width=_caption_width(ncol)),
             fontsize=7.3, color=INK_2, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.13, 1, 0.93))
    fig.savefig(out_dir / "profile_gallery.png", dpi=150,
               bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _hot_spans(
    hot_as_of: list[pd.Timestamp], quarters: list[pd.Timestamp],
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Turn a list of alarming as-of quarters into (window_start, as_of)
    spans for shading -- each alarm is a property of the whole
    `RECENT_QUARTERS`-wide window ending at that as-of, not a single point."""
    pos = {q: i for i, q in enumerate(quarters)}
    spans = []
    for q in hot_as_of:
        i = pos[q]
        spans.append((quarters[max(0, i - RECENT_QUARTERS + 1)], q))
    return spans


def _plot_findings_close_up(out_dir: Path, seed: int) -> None:
    """The specific synthetic populations behind `docs/findings/
    synthetic.md` Findings 1 and 2, with the naive `nb_trend_z` alarm and its
    fix drawn directly on the raw counts that produced them -- not a
    re-derivation of the findings' headline numbers (those are already
    tabulated there), just making a handful of individual cases visible.
    """
    rng = np.random.default_rng(seed)
    T = 45

    # Row 1 -- Finding 1: repeated-testing false alarms on genuinely flat
    # series. A larger population than the panel needs, so there are enough
    # "hot but never confirmed" cases to choose the clearest examples from.
    flat = [Profile(f"flat_{i}", baseline_rate=float(rng.uniform(2, 20)),
                    rises=False, noise="poisson")
            for i in range(150)]
    counts1, denom1, _ = generate_panel(flat, T=T, seed=seed)
    quarters1 = sorted(denom1.index)
    sweep1 = nb_trend_sweep(counts1, denom1)
    row1 = []
    for sub, g in sweep1.groupby("substance"):
        hot = g.loc[g["nb_trend_z"] > NB_TREND_THRESHOLD, "as_of"].tolist()
        if hot and not g["confirmed"].any():
            row1.append(sub)
    row1 = row1[:4]

    # Row 2 -- Finding 2: denominator-drift false alarms. Same flat-profile
    # recipe, declining county denominator instead of flat.
    drift = [Profile(f"drift_{i}", baseline_rate=float(rng.uniform(2, 20)),
                     rises=False, noise="poisson")
             for i in range(80)]
    counts2, denom2, _ = generate_panel(drift, T=T, denom_mode="declining",
                                        seed=seed + 1)
    quarters2 = sorted(denom2.index)
    sweep2 = nb_trend_sweep(counts2, denom2)
    suppressed, still_alarms = [], []
    for sub, g in sweep2.groupby("substance"):
        hot = g.loc[g["nb_trend_z"] > NB_TREND_THRESHOLD]
        if not len(hot):
            continue
        (still_alarms if hot["emergent"].any() else suppressed).append(sub)
    # At least one "still alarms" case if this draw produced any -- the rarer,
    # more informative outcome (Finding 2's own text: most flat baseline
    # rates here exceed EMERGENT_THRESHOLD, so suppression is the common
    # case; showing only that would hide the caveat the finding itself makes.
    row2 = [(s, True) for s in still_alarms[:1]]
    row2 += [(s, False) for s in suppressed[:4 - len(row2)]]

    ncol = max(len(row1), len(row2), 1)
    fig, axes = plt.subplots(2, ncol, figsize=(3.3 * ncol, 2.3 * 2 + 1.1),
                             sharex=True)
    axes = np.atleast_2d(axes)

    for j in range(ncol):
        ax = axes[0, j]
        if j < len(row1):
            sub = row1[j]
            g = sweep1[sweep1["substance"] == sub]
            hot_q = g.loc[g["nb_trend_z"] > NB_TREND_THRESHOLD, "as_of"].tolist()
            s = (counts1[counts1["substance"] == sub]
                 .set_index("quarter")["n"].reindex(quarters1, fill_value=0))
            rate = np.full(T, [p for p in flat if p.name == sub][0].baseline_rate)
            _panel(ax, quarters1, s.to_numpy(), rate,
                  f"{sub}  (truly flat)",
                  badge="NAIVE ALARM,\nNEVER CONFIRMED",
                  hot_spans=_hot_spans(hot_q, quarters1))
        else:
            ax.set_visible(False)

    for j in range(ncol):
        ax = axes[1, j]
        if j < len(row2):
            sub, was_emergent = row2[j]
            g = sweep2[sweep2["substance"] == sub]
            hot_q = g.loc[g["nb_trend_z"] > NB_TREND_THRESHOLD, "as_of"].tolist()
            s = (counts2[counts2["substance"] == sub]
                 .set_index("quarter")["n"].reindex(quarters2, fill_value=0))
            rate = np.full(T, [p for p in drift if p.name == sub][0].baseline_rate)
            badge = ("STILL ALARMS\n(rare enough to pass)" if was_emergent
                     else "SUPPRESSED\n(too common to flag)")
            _panel(ax, quarters2, s.to_numpy(), rate,
                  f"{sub}  (flat count, falling denom)",
                  badge=badge, hot_spans=_hot_spans(hot_q, quarters2))
        else:
            ax.set_visible(False)

    axes[0, 0].set_ylabel("Finding 1: repeated testing", fontsize=8.5,
                          color=INK, labelpad=10)
    axes[1, 0].set_ylabel("Finding 2: denominator drift", fontsize=8.5,
                          color=INK, labelpad=10)
    fig.suptitle("Why nb-trend false-alarmed on data with no real signal: "
                 "two independent causes, made visible",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.text(0.01, 0.005, _caption(
             "Top row: every count series here is genuinely flat, with no "
             "trend built into it at all. The yellow band marks a quarter "
             "where the trend test alarmed anyway, purely by chance -- a "
             "single look in an unbroken multi-year sweep, not a real "
             "signal. None of these ever alarmed for three quarters in a "
             "row, which is the bar an actually-confirmed detection needs "
             "to clear.",
             "Bottom row: also flat counts, but the countywide "
             "denominator is declining, so each substance's share of all "
             "deaths drifts upward even though its own raw count does "
             "not -- the trend test alarms on that drift, mistaking it "
             "for a real signal. Most of these get suppressed by a "
             "separate rule that only lets already-common substances "
             "alarm, but that rule was not built to catch denominator "
             "drift specifically, and one substance here alarms anyway.",
             width=_caption_width(ncol)),
             fontsize=7.1, color=INK_2, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.17, 1, 0.94))
    fig.savefig(out_dir / "finding_close_up.png", dpi=150,
               bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------
# Nearest real substance. A separate diagnostic from the plots above: not
# "does a detector fire correctly on this shape" but "does any real LA
# County substance actually look like this shape at all", the honesty check
# a hand-picked set of onset/magnitude/duration numbers otherwise has no way
# to answer.
# --------------------------------------------------------------------------

MIN_REAL_TOTAL = 15     # total mentions across the whole extract
MIN_REAL_QUARTERS = 12  # distinct quarters with at least one mention
MIN_MATCH_OVERLAP = 16  # shortest overlap a lag is allowed to score on
THIN_MATCH_MEAN = 3.0   # below this mean count *in the matched window*, flag the match as thin


def _best_lag_correlation(
    a: np.ndarray, b: np.ndarray, min_overlap: int = MIN_MATCH_OVERLAP,
) -> tuple[float, int, int, float] | None:
    """Highest Pearson correlation between `a` and `b` over every time shift
    that leaves at least `min_overlap` quarters of overlap -- shape
    similarity, deliberately not scored on absolute level or on the two
    series happening to share a calendar. A synthetic profile's onset index
    and a real substance's actual history have no reason to line up in
    time, so the comparison that matters is "does the *shape* match at some
    alignment", not "do the raw quarters match". Pearson r already
    normalizes for the two series' very different count scales (a real
    substance can run in the hundreds/quarter, a synthetic one in the
    single digits) -- no separate rescaling needed.

    Returns `(r, lag, overlap_n, mean_b_in_window)` for the winning lag, or
    `None` if no lag has enough overlap or either side is constant there.
    `mean_b_in_window` (the real substance's own mean count over just the
    matched window, not its whole history) exists so a caller can flag a
    match won on a short, thin, mostly-zero window -- Pearson r on a few
    single-digit counts can read deceptively high without meaning the
    shapes genuinely resemble each other the way a high r on a well-
    populated window does.
    """
    la, lb = len(a), len(b)
    best = None
    for lag in range(-(lb - min_overlap), la - min_overlap + 1):
        if lag >= 0:
            n = min(la - lag, lb)
            x, y = a[lag:lag + n], b[:n]
        else:
            n = min(la, lb + lag)
            x, y = a[:n], b[-lag:-lag + n]
        if n < min_overlap or x.std() == 0 or y.std() == 0:
            continue
        r = float(np.corrcoef(x, y)[0, 1])
        if best is None or r > best[0]:
            best = (r, lag, n, float(y.mean()))
    return best


def nearest_real_matches(
    synthetic_counts: pd.DataFrame, profiles: list[Profile],
    mentions_path: Path = MENTIONS_PATH,
) -> pd.DataFrame:
    """For each `profiles` entry, the real LA County substance whose actual
    quarterly mention history correlates best (`_best_lag_correlation`) with
    that profile's own sampled series -- an honesty check on the generator's
    design, not a detector benchmark: if nothing real ever looks like a
    given synthetic shape, that shape may be testing a scenario this
    project has no evidence actually occurs here.

    Real substances are restricted to `MIN_REAL_TOTAL` total mentions and
    `MIN_REAL_QUARTERS` distinct quarters seen, so a substance with one or
    two lifetime mentions (most of the 211-substance catalog) cannot win a
    match purely by being mostly zero, which correlates with nothing
    meaningfully but can score deceptively high noise-to-noise.
    """
    real_counts, _ = load_quarterly(mentions_path, quiet=True)
    real_quarters = sorted(real_counts["quarter"].unique())
    wide = (real_counts.pivot_table(index="substance", columns="quarter",
                                    values="n", aggfunc="sum", fill_value=0)
            .reindex(columns=real_quarters, fill_value=0))
    keep = ((wide.sum(axis=1) >= MIN_REAL_TOTAL)
            & ((wide > 0).sum(axis=1) >= MIN_REAL_QUARTERS))
    wide = wide.loc[keep]

    rows = []
    for p in profiles:
        s = (synthetic_counts[synthetic_counts["substance"] == p.name]
             .set_index("quarter")["n"])
        s = s.reindex(sorted(synthetic_counts["quarter"].unique()),
                      fill_value=0).to_numpy(dtype=float)
        scores = {sub: _best_lag_correlation(s, row.to_numpy(dtype=float))
                  for sub, row in wide.iterrows()}
        scores = {k: v for k, v in scores.items() if v is not None}
        if not scores:
            rows.append({"profile": p.name, "nearest_real_substance": None,
                         "correlation": None, "window_start": None,
                         "window_end": None, "window_quarters": None,
                         "thin_match": None})
            continue
        best_sub = max(scores, key=lambda k: scores[k][0])
        r, lag, n, mean_in_window = scores[best_sub]
        start_idx = 0 if lag >= 0 else -lag
        window = real_quarters[start_idx:start_idx + n]
        rows.append({
            "profile": p.name,
            "nearest_real_substance": best_sub,
            "correlation": round(r, 2),
            "window_start": f"{window[0].year}Q{window[0].quarter}",
            "window_end": f"{window[-1].year}Q{window[-1].quarter}",
            "window_quarters": n,
            "thin_match": mean_in_window < THIN_MATCH_MEAN,
        })
    return pd.DataFrame(rows)


@app.command("match")
def match(
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    seed: int = typer.Option(0, "--seed"),
) -> None:
    """For each `GALLERY_PROFILES` entry, print and save the real LA County
    substance whose actual history looks most like it (`nearest_real_matches`)
    -- do any of these synthetic shapes correspond to something that actually
    happened here, or are some of them purely hypothetical stress tests?
    """
    counts_df, _, _ = generate_panel(GALLERY_PROFILES, T=40, seed=seed)
    matches = nearest_real_matches(counts_df, GALLERY_PROFILES)
    out_dir.mkdir(parents=True, exist_ok=True)
    matches.to_csv(out_dir / "nearest_real_match.csv", index=False)
    typer.echo(matches.to_string(index=False))
    if matches["thin_match"].fillna(False).any():
        thin = matches.loc[matches["thin_match"].fillna(False), "profile"].tolist()
        typer.echo(f"\nthin match (matched window's real mean count < "
                   f"{THIN_MATCH_MEAN:g}/quarter -- treat the correlation "
                   f"as noise-level, not a real resemblance): {', '.join(thin)}")
    typer.echo(f"\nwrote {out_dir / 'nearest_real_match.csv'}")


@app.command("plot")
def plot(
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    seed: int = typer.Option(0, "--seed"),
) -> None:
    """Render this module's own generator, not real data: `profile_gallery.png`
    (what each ground-truth category and generative knob actually produces as
    a line -- `SYNTHETIC_BENCHMARK_DESIGN.md`'s taxonomy made visible) and
    `finding_close_up.png` (the specific kind of synthetic substance
    `docs/findings/synthetic.md` Findings 1 and 2 used to catch nb-trend's
    repeated-testing and denominator-drift false alarms, with the naive
    alarm and its fix drawn on the raw counts that produced them).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_gallery(out_dir, seed)
    _plot_findings_close_up(out_dir, seed)
    typer.echo(f"wrote {out_dir / 'profile_gallery.png'}")
    typer.echo(f"wrote {out_dir / 'finding_close_up.png'}")
