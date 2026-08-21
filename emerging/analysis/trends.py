"""
Rank and plot substances whose share of LA County overdose deaths is rising.

Detection statistic. The candidates that matter are semi-rare — xylazine has 9
lifetime mentions, nitazenes 3 — and at those counts a raw ratio of recent to
baseline share is worthless: one extra death doubles it. So this uses the
pharmacovigilance answer, a gamma-Poisson empirical-Bayes shrinkage:

    n_s ~ Poisson(lambda_s * E_s),  lambda_s ~ Gamma(alpha, beta)

where `E_s` is the count expected in the recent window if substance `s` held its
baseline share of all overdose deaths. alpha/beta are fitted across substances by
maximum likelihood (the marginal is negative binomial), giving a posterior
Gamma(alpha + n_s, beta + E_s) per substance. Ranking is by **EB05**, its 5th
percentile — the conservative end, so a substance only ranks high when the data
cannot easily be explained by noise. `ebgm`, the posterior *geometric* mean
exp(E[log lambda]), is reported alongside it.

This is a single-component GPS, not full MGPS. DuMouchel's two-component
mixture was implemented first and is not identifiable at this scale — see
`_fit_gps_prior` for what it did instead. It becomes worth revisiting on a
substance x zip x quarter table.

Does it work? `backtest` rolls the as-of date back and re-scores, so the
ranking is falsifiable rather than merely quiet. All five known LA emergences
reach rank 1-4 in the quarter they emerge, xylazine included at n=6.

Right-truncation. Toxicology takes months, so recent quarters are undercounts
and every rising substance looks like it is falling. Two guards, both applied in
`load_quarterly` before any statistic runs: `--cutoff` (default 2025-12-31)
drops quarters whose toxicology is not yet settled, and partial calendar
quarters at the end of the extract are always dropped. `--lag-quarters` is the
older rolling-haircut alternative and defaults to 0, because stacking it on the
cutoff double-corrects. Novel substances confirm *slower* than fentanyl, so even
this is a floor, not a fix — see the README.

Windowing matters more than it looks. The recent/baseline split is a fixed
number of quarters back from the cutoff, so moving the cutoff slides the window.
For a substance whose rise is a short burst rather than a sustained trend, that
can move EB05 across the alarm line — lidocaine does exactly this. Quote the
`backtest` peak alongside the current value when the trajectory is spiky.

Novelty. `date_added` from NFLIS is *catalog* novelty (when labs could report a
substance), left-censored at the 1998-10 inception for 27% of the catalog. Local
epidemiological novelty is `first_seen`, the first LACME quarter — itself
left-censored at 2012, where the extract begins. The scatter plots the two
against each other, because the interesting region is where they disagree: old
chemistry arriving in LA for the first time.

Four figures, written to results/figures/:
    rising_candidates.png  small multiples of the current top-ranked substances
    novelty_vs_growth.png  chemical novelty vs LA arrival
    share_heatmap.png      substance x quarter share of all OD deaths
    known_emergences.png   the backtest set

Usage:
    emerging trends rank
    emerging trends rank --recent-quarters 8 --top 30
    emerging trends dispersion
    emerging trends backtest
    emerging trends plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from matplotlib.colors import PowerNorm
from matplotlib.transforms import Bbox
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import chi2, gamma as gamma_dist

from emerging.config import (BASELINE_QUARTERS, CUTOFF, LAG_QUARTERS,
                             RECENT_QUARTERS)
from emerging.core.causeline import causea_order, rollup_map
from emerging.core.concentration import (concentration_weight,
                                         load_spatial_cohort,
                                         sweep_concentration)
from emerging.ingest.extract import MENTIONS_PATH, load_cohort
from emerging.paths import results_dir
from emerging.viz import AQUA, BLUE, INK, INK_2, MUTED, ORANGE, SEQ, YELLOW

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("trends")

# Dated LA-local emergences used as the backtest ground truth. Only
# para-fluorofentanyl is well powered; the rest are shown for honesty about how
# thin the evidence is at LA-county scale.
KNOWN_EMERGENCES = [
    "para-Fluorofentanyl", "Bromazolam", "Mitragynine", "Carfentanil", "Xylazine",
]

# EB05 above which a substance is "in alarm". Read off the backtest: every
# known LA emergence clears it in the quarter it emerges, and it is high enough
# that the bulk of the catalog never does — 24 substances of 219 have ever
# breached it across 45 as-of quarters. It is a calibrated reading line, not a
# significance test; EB05 > 1 already means "credible growth".
ALARM_THRESHOLD = 1.5
ALARM_PATH = RESULTS_DIR / "alarm_history.csv"

# Average pre-rise quarterly death count below which an alarm is also
# annotated `emergent`, not merely `credible_rise`. "Rising" (EB05) and
# "was rare" are separate axes -- fentanyl and PCP clear the alarm line while
# already killing 3.4-3.9/quarter before their rise even starts, which is a
# supply shift in an established killer, not an emergence. Read off the same
# backtest as `ALARM_THRESHOLD`: at first breach the four known LA
# emergences' *baseline* average (`n_baseline` at that quarter, divided by
# `baseline_quarters`) is 0-1.125/quarter, with a clean gap to PCP (3.375)
# and fentanyl (3.875) -- the two alarms in `alarm_history.csv` that are
# obviously not emergences. A chosen reading line, like `ALARM_THRESHOLD`
# itself, not a fitted one.
EMERGENT_THRESHOLD = 2.0

# TreeScan veto floor, as a recurrence interval (`emerging.analysis.treescan`).
# `credible_rise` is suppressed only when TreeScan's own, independent
# substance-leaf scan sees essentially no excess at all this quarter -- far
# below TreeScan's own RI >= 100 alarm line, which would instead re-impose
# TreeScan's detection threshold on top of EB05's and lose real detections
# (see `docs/findings/benchmark.md`'s "mean"/"max" ensembles, both of which
# do exactly that and lose to EB05 alone). Validated on the 22-substance
# backtest as `ensemble-veto-v2role-10`: with `--weighted --role-discount
# --spatial discount` as the case definition, this is a same-or-better
# result than the case definition alone on every column at the matched
# budget, and the best negative-avoidance (6/7) of any configuration tried.
TREESCAN_VETO_THRESHOLD = 10.0

# Substances driving the national emerging-drug conversation, checked against
# LA regardless of whether the ranking surfaces them. The ranking can only
# promote what is already in the data; this asks the complementary question,
# "is the thing everyone else is seeing here at all", and a zero is as much a
# result as a hit. Every entry resolves in the lexicon, so a zero means the
# scan looked and found nothing — not that the name was unrecognised.
WATCHLIST: dict[str, str] = {
    "Xylazine": "'tranq' — fentanyl adulterant, widespread in the Northeast",
    "Medetomidine": "vet sedative displacing xylazine in East Coast fentanyl, 2024-",
    "Dexmedetomidine": "medetomidine isomer, same supply signal",
    "BTMPS (Bis(2,2,6,6-tetramethyl-4-piperidinyl) Sebacate)":
        "industrial plastics stabiliser found in fentanyl from mid-2024",
    "Tianeptine": "'gas-station heroin', unscheduled atypical opioid",
    "Bromazolam": "designer benzo, the successor to etizolam/flualprazolam",
    "Carfentanil": "resurgent 2023- after a near-decade absence",
    "Levamisole": "long-standing cocaine adulterant",
    "Desomorphine": "'krokodil' — recurrent media claim, rarely confirmed",
}
# Handled as a family: individually every nitazene is too rare to score, and
# the public-health question is about the class, not the analog.
NITAZENE_RE = r"nitaz"

# Position-weighted case definition (`docs/GPS_V2_DESIGN.md` #1).
#
# A substance is discounted only when it is named *last* on a line of three or
# more -- the adulterant pattern documented for lidocaine in
# `docs/findings/polysubstance.md` (last on 82% of its lines, first on none).
# A substance that is a necessary contributor but rarely or never lead -- a
# benzodiazepine potentiating opioid toxicity, say -- is not on that pattern
# and keeps full credit; only "habitually relegated to the trailing slot" is
# treated as evidence of a non-causal role. A 2-substance line cannot tell
# "trailing adulterant" from "just not first" apart (the same short-line
# caution `polysubstance.profile_table` documents), so last-of-2 also keeps
# full credit.
#
# **Integer-valued does not mean Poisson-safe.** The first cut of this
# reasoned that keeping credit integer-valued was enough to keep the weighted
# sum a legitimate count for `_nb_logpmf`'s Poisson assumption. It is not: for
# a per-mention credit c with any spread across {0, 1, 2}, c^2 >= c pointwise
# whenever c is 0 or >=1, so E[c^2] >= E[c] and the weighted sum is
# structurally *overdispersed* relative to a Poisson with the same mean --
# integer or not. Fed straight into `_fit_gps_prior`, that overdispersion
# reads as signal: on this repository's data it inflated the false-EB05>1
# rate from 5.1/203 substances (raw counts) to 19.0/203. `SOLO_CREDIT == 2` is
# the offender; a pure {0, 1} scheme has c^2 == c identically, so phi == 1 no
# matter how the keep/discard split varies case to case, which is what makes
# it dispersion-neutral *by construction* rather than approximately so. That
# is kept as {0, 1, 2} anyway -- a solo mention is arguably the strongest
# available evidence of causal weight, not merely of independence, which is
# what `polysubstance.alone_pct` already measures -- and the dispersion is
# corrected explicitly instead: see `_credit_dispersion` and its use in
# `rank_substances`.
# Spatial concentration (`docs/GPS_V2_DESIGN.md` #3), `spatial_mode="discount"`.
# The maximum share of a substance's expected count the concentration term may
# discount. 0.25 is a *chosen* number, and the design note is explicit that
# this being chosen rather than estimated is the fallback's weakness relative
# to the covariate-adjusted prior. It is set where it is because
# E -> 0.75*E raises EB05 by roughly what a 33% rise in the observed count
# would -- large enough to change a borderline ranking, small enough that
# geography alone can never manufacture an alarm out of a flat count. A
# substance at conc == 1 (z_s well past 2) and n_recent == n_baseline still
# scores EB05 below the 1.5 line on its own.
SPATIAL_WEIGHT = 0.25

# Where the per-(substance, as-of) concentration sweep is cached. It costs
# about 5s to recompute, so this is a convenience rather than a dependency --
# `_spatial_sweep` regenerates it silently when absent, unlike
# `ALARM_PATH`, which `plot` refuses to recompute.
SPATIAL_PATH = RESULTS_DIR / "spatial_concentration.csv"

SOLO_CREDIT = 2       # unambiguous: nothing else named
NAMED_CREDIT = 1      # named, and not (last on a line of 3+)
TRAILING_CREDIT = 0   # last on a line of 3+ -- the adulterant pattern
# Mentions the cause-line scan cannot place (it has no fuzzy fallback, see
# `causea_order`) get neutral credit: no ordering evidence is not evidence of
# a trailing role.
UNMATCHED_CREDIT = 1


def _ever_independent(order: dict[str, list[str]]) -> frozenset[str]:
    """Substances that have, at least once, been named alone or first on the
    cause line -- direct evidence the ME considered them capable of
    independently dominating a death, even if they usually ride along with
    something else.

    Fixes the construct-validity problem the design note's validation caught:
    xylazine is a fentanyl adulterant by pattern (habitually trailing,
    same as lidocaine) but is not pharmacologically inert the way lidocaine
    is, and the surveillance question for it is "is this entering the
    supply," not "did the ME consider it the proximate cause here." On this
    repository's data the split is clean: lidocaine and levamisole (known
    inert cutting agents) are named alone or first *zero* times across their
    entire histories (42 and 12 mentions); xylazine clears it once, in the
    same quarter its emergence is later detected.

    Solo and lead collapse to the same check: `seq[0]` is the substance
    itself when `len(seq) == 1` (solo) and the lead substance otherwise, so
    "named alone or first" is exactly "appears at index 0 of some case's
    ordered sequence" -- no separate solo-mention branch needed.
    """
    return frozenset(seq[0] for seq in order.values() if seq)


def _credit_from_order(
    order: dict[str, list[str]], independent: frozenset[str] = frozenset(),
) -> pd.DataFrame:
    """Per (CaseNumber, substance) -> integer credit from cause-line position.

    `independent` is the output of `_ever_independent`: a substance in that
    set never receives `TRAILING_CREDIT`, even in a case where it is named
    last on a 3+-substance line, because other evidence already establishes
    it can be the whole story.

    Pure function of an already-computed case -> ordered-substances map (and
    an already-computed independence set), so it is testable without the
    lexicon or cohort I/O `causea_order` needs.
    """
    rows = []
    for case, seq in order.items():
        L = len(seq)
        for idx, sub in enumerate(seq):
            if L == 1:
                c = SOLO_CREDIT
            elif idx == L - 1 and L >= 3:
                c = NAMED_CREDIT if sub in independent else TRAILING_CREDIT
            else:
                c = NAMED_CREDIT
            rows.append((case, sub, c))
    return pd.DataFrame(rows, columns=["CaseNumber", "substance", "credit"])


def _position_credit(cohort: pd.DataFrame, roll: dict[str, str]) -> pd.DataFrame:
    """Cause-line order, converted to credit. See `_credit_from_order`.

    **As-of caution.** Per-case credit only ever reads that case's own
    cause-line text, so slicing the resulting counts by quarter afterward
    (as `backtest`/`sweep_eb05` do to the *unweighted* counts today) cannot
    leak future information. `_ever_independent` breaks that property: it is
    a cross-case aggregate, so a substance's *later* independence evidence
    would retroactively un-discount its *earlier* trailing mentions if this
    were computed once and reused for every as-of quarter. Callers that
    replay history (`backtest`, `sweep_eb05`, once either supports
    `--weighted`) must call this -- via `load_quarterly` -- with `cohort`
    already restricted to `quarter <= as_of`, once per as-of step, not once
    up front. `load_quarterly` does this by construction: it is always given
    the cohort already filtered to its own `cutoff`.
    """
    order = causea_order(cohort, roll)
    return _credit_from_order(order, _ever_independent(order))


def _credit_dispersion(credit: pd.Series) -> float:
    """phi = E[c^2] / E[c], the dispersion a per-mention credit distribution
    introduces relative to a plain Poisson count (phi == 1 for raw mentions,
    where every credit is exactly 1).

    A compound-Poisson sum W = sum(c_i) over a Poisson-count number of
    mentions has Var(W) = mu * E[c^2] against Mean(W) = mu * E[c], so phi is
    exactly the variance-to-mean ratio the weighted sum carries beyond what a
    same-mean Poisson would -- dividing both the count and the expectation
    that go into `_fit_gps_prior` by phi (see `rank_substances`) rescales the
    effective sample size down to what the credit scheme's own variance
    supports, which is the standard quasi-likelihood correction for
    overdispersed count data. `gammaln` is defined on the reals, so feeding
    the rescaled, non-integer n through `_nb_logpmf` is mechanically fine;
    what changes is that its output is a quasi-likelihood rather than a
    normalized pmf, which is the accepted trade for overdispersion this way.
    """
    c = credit.to_numpy(dtype=float)
    m1 = c.mean()
    return float((c ** 2).mean() / m1) if m1 > 0 else 1.0


# Minimum credited mentions before a substance gets its own role discount.
# `_role_dispersion` is a per-substance aggregate over exactly the deaths
# that name it, the same statistical-power problem `polysubstance.MIN_CASES`
# already guards against -- below this, one or two differently-timed deaths
# would swing phi_s wildly, and the substance keeps phi_s = 1 (no extra
# discount) rather than a number this repository's rarest candidates
# (xylazine, n=9) cannot support.
MIN_CASES_FOR_ROLE_DISCOUNT = 10

# Floor on a substance's own mean credit before computing phi_s, so a
# substance whose *every* credited mention happens to be zero (fully
# trailing, no exceptions at all) gets a large but finite discount rather
# than a division by zero.
_ROLE_MEAN_CREDIT_FLOOR = 0.05


def _role_dispersion(
    order: dict[str, list[str]], independent: frozenset[str],
) -> dict[str, float]:
    """Per-substance phi_s = max(1, NAMED_CREDIT / mean(credit)) -- an
    *additional* effective-sample-size discount, on top of the per-death
    position credit and the global phi correction, for a substance whose own
    aggregate role still reads as adulterant-patterned even where per-death
    position cannot prove it.

    **Why per-death credit alone leaves a gap.** Only a mention *trailing on
    a line of 3+* is discounted (`_credit_from_order`); a 2-substance line
    keeps full credit even when it's "just co-occurring with the same one
    thing every time," because a 2-line can't distinguish that from "really
    was a contributor." Checked directly against real data: 12 of
    Lidocaine's 34 recent mentions are literally `[Fentanyl, Lidocaine]` --
    full per-death credit every time, by design, and enough raw signal on
    their own to keep tripping the alarm after the per-death discount zeroes
    the other 22. `_credit_from_order` cannot fix this without breaking the
    short-line carve-out generally; this function targets it a different
    way, by reading the substance's own *aggregate* mean credit.

    **The reference point is NAMED_CREDIT (1), not SOLO_CREDIT (2) -- this
    was wrong in the first cut and caught by running `backtest`, not by
    inspection.** An `_ever_independent`-exempted substance never receives
    `TRAILING_CREDIT`, so its credit values are all 1 or 2 -- but *solo* is
    genuinely rare for almost any real drug, dangerous or not: checked
    directly on this repository's data, Bromazolam sits at mean credit
    1.000, Carfentanil 1.167, Mitragynine 1.013, PCP 1.001, Fentanyl 1.223 --
    all comfortably at or above 1, none anywhere near 2. Using `SOLO_CREDIT`
    as the reference punished every one of them for not being frequently
    solo, a bar almost nothing clears, and gutted `backtest`'s known
    emergences across the board (Bromazolam's peak EB05 dropped from 2.69 to
    1.00, Xylazine's dual gate stopped clearing entirely) -- a real
    regression the first version of this function shipped with, caught only
    by re-running `backtest`, the same discipline that caught the xylazine
    problem in `_ever_independent` itself. `NAMED_CREDIT` -- "was this
    mention ever relegated to the trailing slot" -- is the right bar: a
    substance sitting at or above it (never trailing) gets `phi_s = 1`,
    exactly the substances above; Lidocaine (mean 0.310) and Levamisole
    (mean 0.667) sit clearly below it and still get a real discount (phi_s
    3.2 and 1.5 respectively). The `max(1, ...)` floor matters on its own:
    without it, any substance with mean credit above 1 (Fentanyl, at 1.223)
    would get `phi_s < 1`, which *amplifies* its evidence rather than
    discounting it -- backwards for what this function is for.

    **A constant, substance-level factor doesn't distort the ranking the way
    it might look like it should.** `phi_s` divides both `n` and `E`
    (`rank_substances`) by the same amount for a substance, so it leaves the
    raw ratio `n/E` exactly unchanged -- checked numerically on Lidocaine:
    the ratio is 1.971 whether n and E are full-scale or both divided down.
    What it changes is the *posterior*: `_fit_gps_prior`'s `alpha`/`beta`
    are fit once, globally, across every substance, so shrinking one
    substance's own `n`/`E` pulls its own posterior harder toward that fixed
    population mean -- the identical mechanism `_credit_dispersion`'s global
    `phi` already uses, just substance-specific rather than pooled. This is
    a deliberate design choice, not a coincidence: a substance-specific
    factor is only informative if it is allowed to change the *how sure are
    we* answer, since it structurally cannot change the *what's the rate*
    answer, and trying to make it change the rate directly (e.g. dividing
    only `n`) would give a substance's full lifetime reputation a permanent,
    asymmetric veto over ever detecting a real change in its own behaviour.

    **As-of caution — identical to `_ever_independent`, for the identical
    reason.** This is a cross-case aggregate over each substance's own
    mentions, so it must be computed from `order` already restricted to
    `quarter <= as_of`, recomputed fresh at every as-of step, never cached
    across a `backtest` sweep and reused for earlier quarters.
    """
    credit = _credit_from_order(order, independent)
    out: dict[str, float] = {}
    for sub, g in credit.groupby("substance"):
        if len(g) < MIN_CASES_FOR_ROLE_DISCOUNT:
            continue
        m = max(float(g["credit"].mean()), _ROLE_MEAN_CREDIT_FLOOR)
        out[sub] = max(1.0, NAMED_CREDIT / m)
    return out


def load_quarterly(
    mentions_path: Path = MENTIONS_PATH, rollup: bool = True,
    cutoff: str | None = CUTOFF, weighted: bool = False, quiet: bool = False,
    role_discount: bool = False,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (substance x quarter counts, all-OD deaths per quarter).

    Metabolites are rolled up to their parent so NORFENTANYL does not read as
    its own rising substance; class terms are dropped.

    Two completeness filters are applied here, before any statistic sees the
    data. `cutoff` drops quarters ending after a settled-toxicology date (see
    `CUTOFF`). Separately, trailing quarters the extract does not even span are
    always dropped, whatever the cutoff.

    `weighted` replaces the raw mention count (every death worth 1) with the
    position-weighted case definition from `docs/GPS_V2_DESIGN.md` #1: a death
    is discounted to 0 only when the substance is named last on a cause line of
    three or more *and* the substance has never been named alone or first
    anywhere in the cohort up to `cutoff` (`_ever_independent`) -- the
    adulterant pattern, unless other evidence already rules that out. It is a
    case-definition change, not a different statistic -- everything
    downstream of `load_quarterly` is unchanged.

    `role_discount` (only meaningful when `weighted=True`) additionally
    applies `_role_dispersion`'s per-substance phi_s -- see that function's
    docstring for the Lidocaine case it targets (a substance whose mentions
    are overwhelmingly short co-occurrence lines, which per-death credit
    alone cannot discount) and why a constant per-substance factor changes
    the posterior rather than the ranked ratio. Off by default: it changes
    `--weighted`'s existing, already-validated behaviour, and should be
    opted into per `backtest` run, not silently folded in.

    `quiet` suppresses the diagnostic `typer.echo` lines below (and
    `load_cohort`'s alcohol-filter line). Off by default; a caller that
    re-invokes this once per as-of quarter (`backtest --weighted`, see
    `_position_credit`) sets it, or the same three lines print on every
    iteration and bury the actual output.
    """
    m = pd.read_parquet(mentions_path)
    m = m[~m["is_class_term"].fillna(False)]
    key = "rollup" if rollup else "canonical"
    m[key] = m[key].fillna(m["canonical"])

    # One row per (case, substance) after rollup — a death naming both fentanyl
    # and norfentanyl must count once.
    m = m.drop_duplicates(subset=["CaseNumber", key])

    cohort = load_cohort(verbose=not quiet)
    phi_s: dict[str, float] = {}

    if weighted:
        roll = rollup_map(m) if rollup else {}
        # `_ever_independent` is a cross-case aggregate (see `_position_credit`
        # for why that matters) -- bound it to `cutoff`, not the full extract,
        # so a caller replaying history one `cutoff` at a time (`backtest`,
        # once it recomputes rather than slices -- see there) never lets
        # evidence from after the point being scored change an earlier score.
        credit_cohort = cohort
        if cutoff is not None:
            credit_cohort = cohort[cohort["quarter"] <= pd.Timestamp(cutoff)]
        order = causea_order(credit_cohort, roll)
        independent = _ever_independent(order)
        credit = _credit_from_order(order, independent).rename(
            columns={"substance": key})
        m = m.merge(credit, how="left", on=["CaseNumber", key])
        m["credit"] = m["credit"].fillna(UNMATCHED_CREDIT).astype(int)
        phi = _credit_dispersion(m["credit"])
        if role_discount:
            phi_s = _role_dispersion(order, independent)
    else:
        m["credit"] = 1
        phi = 1.0  # raw counts are exactly Poisson; no correction to make

    counts = (m.groupby([key, "quarter"])["credit"].sum()
                .rename("n").reset_index().rename(columns={key: "substance"}))

    denom = cohort.groupby("quarter").size().rename("N")

    # Drop trailing quarters the extract does not even span. The file ends
    # 2026-04-04, so 2026Q2 is a 4-day stub holding 3 deaths. This is a
    # different defect from reporting lag and has to be removed separately:
    # left in, it silently consumes one of the `lag_quarters` slots, so a
    # 2-quarter lag correction only discards one real quarter of data.
    last_death = cohort["death_date"].max()
    complete = [q for q in denom.index
                if q + pd.offsets.QuarterEnd(0) <= last_death]
    dropped = [q for q in denom.index if q not in set(complete)]
    if dropped and not quiet:
        typer.echo("dropped partial calendar quarter(s): "
                   + ", ".join(f"{q.year}Q{q.quarter} (n={int(denom[q])})"
                               for q in dropped))
    denom = denom.loc[complete]

    if cutoff is not None:
        cut = pd.Timestamp(cutoff)
        keep = [q for q in denom.index if q + pd.offsets.QuarterEnd(0) <= cut]
        n_cut = len(denom) - len(keep)
        if n_cut and not quiet:
            typer.echo(f"cutoff {cut.date()}: dropped {n_cut} unsettled "
                       f"quarter(s) after {keep[-1].year}Q{keep[-1].quarter}")
        denom = denom.loc[keep]
        complete = keep

    counts = counts[counts["quarter"].isin(set(complete))]
    # Dispersion correction for `rank_substances` -- see `_credit_dispersion`
    # and, for `phi_s`, `_role_dispersion`. Set after the filtering above
    # rather than before: `.attrs` propagation through pandas indexing is not
    # guaranteed, so this is the last write.
    counts.attrs["phi"] = phi
    counts.attrs["phi_s"] = phi_s
    return counts, denom


def _spatial_sweep(
    quarters: list[pd.Timestamp], mentions_path: Path = MENTIONS_PATH,
    recent_quarters: int = RECENT_QUARTERS, cache: Path | None = SPATIAL_PATH,
    refresh: bool = False,
) -> pd.DataFrame:
    """(substance, as_of) -> z_spatial, for the quarters this ranking uses.

    Computed over the *same* trailing `recent_quarters` window the ranking
    scores its numerator over, so the two halves of the statistic are about
    the same deaths. `quarters` is passed in from the caller's already
    cutoff-filtered `denom` index rather than read from the cohort, so the
    concentration sweep inherits the completeness filters exactly instead of
    reimplementing them -- otherwise the 2026Q2 stub `load_quarterly` drops
    would quietly reappear here as a window of its own.

    Cached to CSV because it is the same numbers for every model config in
    `validation/benchmark.py`, which would otherwise recompute it six times.
    The cache carries the `recent_quarters` it was built with and is discarded
    on a mismatch: the as-of *dates* are identical whatever the window length,
    so a set-membership check alone would silently hand 4-quarter windows to a
    caller that asked for 6.
    """
    if cache is not None and cache.exists() and not refresh:
        sw = pd.read_csv(cache, parse_dates=["as_of"])
        fresh = (set(sw["as_of"]) >= set(quarters[recent_quarters - 1:])
                 and "recent_quarters" in sw
                 and set(sw["recent_quarters"]) == {recent_quarters})
        if fresh:
            return sw
    bg, mn = load_spatial_cohort(mentions_path, quiet=True)
    keep = set(quarters)
    sw = sweep_concentration(bg[bg["quarter"].isin(keep)],
                             mn[mn["quarter"].isin(keep)],
                             sorted(keep), recent_quarters=recent_quarters)
    sw["recent_quarters"] = recent_quarters
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        sw.to_csv(cache, index=False)
    return sw


def _z_at(sweep: pd.DataFrame, as_of: pd.Timestamp) -> pd.Series:
    """The `z_spatial` column of one as-of quarter, indexed by substance."""
    g = sweep[sweep["as_of"] == as_of]
    return g.set_index("substance")["z"] if not g.empty else pd.Series(dtype=float)


def _nb_logpmf(n: np.ndarray, e: np.ndarray, a: float, b: float) -> np.ndarray:
    """log P(n | E) for n ~ Poisson(lambda*E), lambda ~ Gamma(a, b).

    E is floored at a tiny epsilon so substances absent from the baseline
    (E == 0) score finitely: n == 0 contributes 0, and n > 0 becomes very
    surprising, which is the correct reading of "this appeared from nothing".
    """
    e = np.maximum(e, 1e-9)
    return (gammaln(n + a) - gammaln(a) - gammaln(n + 1)
            + a * np.log(b / (b + e)) + n * np.log(e / (b + e)))


# 0 disables the cap. See `_fit_gps_prior` for when to set one.
MAX_EXPECTED_FOR_FIT = 0.0

# Bound on the spatial covariate's coefficient in the covariate-adjusted prior
# (`_fit_gps_prior(z=...)`). z is standardized, so gamma = 1 already means a
# one-sd rise in spatial concentration divides the prior rate parameter by e --
# a factor-2.7 shift in the shrinkage target. Anything past 2 is the optimizer
# escaping into a corner of a 3-parameter fit that ~110 informative cells
# cannot pin down, not a finding.
MAX_GAMMA = 2.0


def _fit_gps_prior(
    n: np.ndarray, e: np.ndarray, max_expected: float = MAX_EXPECTED_FOR_FIT,
    z: np.ndarray | None = None,
) -> dict:
    """ML fit of a Gamma(alpha, beta) prior under a Poisson likelihood.

    The marginal of n given E is negative binomial with size alpha and success
    probability beta/(beta+E), so alpha and beta come straight out of an MLE.

    **Why one component and not DuMouchel's two.** MGPS uses a two-gamma
    mixture, and that was tried here first. It is not identifiable at this
    scale: 203 substances, of which only 111 have a non-zero baseline, cannot
    pin down five parameters. Nelder-Mead converged to a different degenerate
    solution depending on what was in the fit set — a null component of
    Gamma(682, 627) with everything in, and a point mass at 1.08 with the
    largest substances excluded. The single component is stable, has a closed
    form for every posterior quantity, and is honest about what the data
    support. The mixture becomes worth revisiting on a substance x zip x
    quarter table, where the cell count is thousands rather than hundreds.

    Substances with E == 0 are excluded from the fit (they carry no information
    about the prior) but are still scored with it. `max_expected` optionally
    also excludes the largest substances, whose counts dominate the likelihood;
    it is off by default because the standard GPS fit uses all cells and the
    uncapped fit is the more conservative of the two (alpha 7.4 vs 19.9 —
    a more diffuse prior, so less shrinkage).

    **`z` is the covariate-adjusted prior of `docs/GPS_V2_DESIGN.md` #3**, the
    design note's *primary* direction for folding spatial concentration into
    the statistic rather than multiplying it in afterwards. Given a
    per-substance spatial-concentration score, the rate parameter becomes

        beta_s = beta_0 * exp(-gamma * z_s)

    with `gamma` estimated by the same MLE that already produces alpha and
    beta_0, so the shrinkage *target* moves with concentration: a lower
    beta_s is a prior expecting a higher lambda, so a spatially concentrated
    substance needs less raw temporal evidence to clear the same EB05 bar.

    The property that makes this preferable to the `E`-discount fallback is
    that **gamma is estimated, not chosen**: if spatial concentration carries
    no information about which substances are genuinely rising, the
    likelihood drives gamma to 0 and the fit reduces exactly to the
    unadjusted one. The discount's `w` cannot do that — it applies whatever
    weight it was handed. `z` is standardized inside the fit set (mean and sd
    are returned in the prior dict) so gamma is on a per-standard-deviation
    scale and the optimizer sees a well-conditioned problem; substances
    excluded from the fit are still scored, using those same two constants.
    gamma is bounded to `MAX_GAMMA` — an unbounded exponential in a
    3-parameter Nelder-Mead fit over ~110 informative cells will happily run
    off to a degenerate corner, the same identifiability wall the two-component
    mixture hit above.
    """
    keep = e > 0
    if max_expected:
        keep &= e <= max_expected
    nf, ef = n[keep], e[keep]
    zf = None if z is None else np.asarray(z, dtype=float)[keep]
    if len(nf) < 10:
        return {"alpha": 0.5, "beta": 0.5, "gamma": 0.0, "z_mean": 0.0,
                "z_sd": 1.0, "n_fit": int(len(nf)), "converged": False}

    z_mean, z_sd = 0.0, 1.0
    if zf is not None:
        z_mean = float(zf.mean())
        sd = float(zf.std())
        z_sd = sd if sd > 1e-9 else 1.0
        zs = (zf - z_mean) / z_sd

    def nll(params: np.ndarray) -> float:
        a, b = np.exp(params[:2])
        if min(a, b) < 1e-9:
            return 1e12
        if zf is None:
            bs = b
        else:
            g = params[2]
            if abs(g) > MAX_GAMMA:
                return 1e12
            bs = b * np.exp(-g * zs)
        ll = _nb_logpmf(nf, ef, a, bs)
        return -float(np.sum(ll)) if np.isfinite(ll).all() else 1e12

    x0 = np.array([0.0, 0.0]) if zf is None else np.array([0.0, 0.0, 0.0])
    res = minimize(nll, x0=x0, method="Nelder-Mead",
                   options={"maxiter": 10000, "maxfev": 10000,
                            "xatol": 1e-8, "fatol": 1e-8})
    a, b = np.exp(res.x[:2])
    return {"alpha": float(a), "beta": float(b),
            "gamma": float(res.x[2]) if zf is not None else 0.0,
            "z_mean": z_mean, "z_sd": z_sd, "n_fit": int(len(nf)),
            "converged": bool(res.success), "nll": float(res.fun)}


def _gps_scores(n: np.ndarray, e: np.ndarray, pr: dict,
                z: np.ndarray | None = None) -> dict:
    """EBGM and posterior quantiles from Gamma(alpha + n, beta_s + E).

    EBGM is the geometric mean exp(E[log lambda | n]), per DuMouchel, rather
    than the arithmetic posterior mean — at these counts the posterior is
    always skewed and the geometric mean is the stable summary.

    `z` reconstructs the per-substance `beta_s = beta_0 * exp(-gamma * z_s)`
    of a covariate-adjusted fit, using the standardization constants the fit
    returned — so substances excluded from the fit (E == 0) are scored on the
    same scale as those inside it. With no covariate, or a fit that drove
    gamma to 0, this collapses to the scalar `beta` and the posteriors are
    identical to the unadjusted ones.
    """
    from scipy.special import digamma

    beta = pr["beta"]
    if z is not None and pr.get("gamma"):
        zs = (np.asarray(z, dtype=float) - pr.get("z_mean", 0.0)) / pr.get(
            "z_sd", 1.0)
        beta = pr["beta"] * np.exp(-pr["gamma"] * zs)

    post_a = pr["alpha"] + n
    post_b = beta + e
    return {
        "ebgm": np.exp(digamma(post_a) - np.log(post_b)),
        "eb05": gamma_dist.ppf(0.05, post_a, scale=1.0 / post_b),
        "eb95": gamma_dist.ppf(0.95, post_a, scale=1.0 / post_b),
    }


def rank_substances(
    counts: pd.DataFrame,
    denom: pd.Series,
    lag_quarters: int = LAG_QUARTERS,
    recent_quarters: int = RECENT_QUARTERS,
    baseline_quarters: int = BASELINE_QUARTERS,
    z_spatial: pd.Series | None = None,
    spatial_mode: str = "none",
    spatial_weight: float = SPATIAL_WEIGHT,
) -> tuple[pd.DataFrame, dict]:
    """Score every substance by EB-shrunken growth in share of OD deaths.

    `z_spatial` is the per-substance spatial-concentration score for *this*
    as-of quarter's recent window (`core.concentration.window_scores`), and
    `spatial_mode` selects which of `docs/GPS_V2_DESIGN.md` #3's two fusions
    it drives. Both leave `lambda = n/E`'s meaning intact; neither is the
    rejected `EBGM * spatial_RR` multiplication, which would silently discard
    the conservative-bound property that is the reason to use EB05 at all.

    - `"none"` (default) — ignore it. Identical output to before this
      parameter existed.
    - `"discount"` — the design note's simpler fallback,
      `E_s <- E_s * (1 - w * conc_s)`, with `conc_s` the bounded [0, 1] squash
      of `z_s` (`concentration_weight`). Applied to Gate A and Gate B alike,
      and applied *before* the prior is fitted rather than after, so the prior
      describes the ensemble under the adjusted expectations instead of being
      fitted to one set of numbers and used on another.
    - `"prior"` — the design note's primary direction, the covariate-adjusted
      prior `beta_s = beta_0 * exp(-gamma * z_s)` with gamma estimated by the
      same MLE. See `_fit_gps_prior`.

    Substances with no spatial score — absent from `z_spatial`, or withheld by
    `concentration.MIN_CASES` / `MAX_FRACTION` — get `z_s = 0`, which is a
    no-op under both fusions. That is the intended reading: no spatial
    evidence either way, so no adjustment either way.
    """
    # Read before `pivot_table`, whose attrs propagation is not guaranteed.
    # 1.0 for unweighted counts (`load_quarterly(weighted=False)` sets it
    # explicitly; a `counts` frame built by hand, as in tests, has none, and
    # "no correction" is the right default for that case too).
    phi = float(counts.attrs.get("phi", 1.0))
    phi_s = dict(counts.attrs.get("phi_s", {}))
    quarters = sorted(denom.index)
    if lag_quarters:
        quarters = quarters[:-lag_quarters]
    if len(quarters) < recent_quarters + baseline_quarters:
        raise ValueError(
            f"need {recent_quarters + baseline_quarters} complete quarters, "
            f"have {len(quarters)}"
        )
    recent_q = quarters[-recent_quarters:]
    base_q = quarters[-(recent_quarters + baseline_quarters):-recent_quarters]

    n_recent_all = float(denom.loc[recent_q].sum())
    n_base_all = float(denom.loc[base_q].sum())

    wide = counts.pivot_table(index="substance", columns="quarter",
                              values="n", aggfunc="sum", fill_value=0)
    wide = wide.reindex(columns=sorted(set(wide.columns) | set(quarters)),
                        fill_value=0)

    n_recent = wide[recent_q].sum(axis=1).astype(float)
    n_base = wide[base_q].sum(axis=1).astype(float)
    per_q = recent_quarters / baseline_quarters

    # Gate A -- share of all OD deaths. Rises when a substance takes a bigger
    # slice of the pie, which a shrinking pie can produce on its own.
    expected = n_base * (n_recent_all / n_base_all)

    # Spatial concentration (`docs/GPS_V2_DESIGN.md` #3). Aligned to the
    # ranking's own index and zero-filled: a substance the concentration sweep
    # never scored is "no spatial evidence", which both fusions treat as no
    # adjustment.
    if spatial_mode not in ("none", "discount", "prior"):
        raise ValueError("spatial_mode must be 'none', 'discount' or 'prior'")
    if z_spatial is None or spatial_mode == "none":
        zs = pd.Series(0.0, index=n_recent.index)
    else:
        zs = z_spatial.reindex(n_recent.index).fillna(0.0).astype(float)
    conc = pd.Series(concentration_weight(zs.to_numpy()), index=zs.index)

    if spatial_mode == "discount":
        # A concentrated substance is expected to have produced *fewer* deaths
        # than its county-wide baseline share implies, so the same observed
        # count is more surprising. Bounded by construction: `conc` is in
        # [0, 1] and `spatial_weight` caps the discount, so E can never be
        # driven to zero and no substance is ever *penalised* for being
        # diffuse -- conc == 0 leaves E untouched.
        expected = expected * (1.0 - spatial_weight * conc)

    nv, ev = n_recent.to_numpy(), expected.to_numpy()
    zv = zs.to_numpy() if spatial_mode == "prior" else None

    # phi > 1 (weighted counts, see `_credit_dispersion`) rescales both n and
    # E down together before fitting or scoring. lambda = n/E, so the point
    # estimate is untouched; what shrinks is the effective sample size behind
    # it, which is the correct response to a count that carries less
    # information per unit than a same-sized raw count would.
    #
    # phi_s (`_role_dispersion`, `docs/GPS_V2_DESIGN.md` #1 role discount)
    # compounds with phi the same way, substance by substance: a substance
    # not in the dict (below `MIN_CASES_FOR_ROLE_DISCOUNT`, or `role_discount`
    # not requested at all) gets phi_s = 1, no additional discount.
    phi_s_arr = np.array([phi_s.get(sub, 1.0) for sub in n_recent.index])
    nv_c, ev_c = nv / phi / phi_s_arr, ev / phi / phi_s_arr
    prior = _fit_gps_prior(nv_c, ev_c, z=zv)
    scores = _gps_scores(nv_c, ev_c, prior, z=zv)
    ebgm, eb05, eb95 = scores["ebgm"], scores["eb05"], scores["eb95"]

    # Gate B -- self-referential, no county-total term at all (`docs/
    # GPS_V2_DESIGN.md` #2). Expected is just this substance's own baseline
    # rate continued flat into the recent window, so a shrinking denominator
    # elsewhere in the county cannot move it. This is `count_ratio`'s
    # denominator, run back through the same GPS shrinkage as Gate A instead
    # of reported as a raw ratio.
    expected_own = n_base * per_q
    if spatial_mode == "discount":
        expected_own = expected_own * (1.0 - spatial_weight * conc)
    ev_own = expected_own.to_numpy()
    ev_own_c = ev_own / phi / phi_s_arr
    prior_own = _fit_gps_prior(nv_c, ev_own_c, z=zv)
    scores_own = _gps_scores(nv_c, ev_own_c, prior_own, z=zv)
    ebgm_own, eb05_own, eb95_own = (scores_own["ebgm"], scores_own["eb05"],
                                    scores_own["eb95"])

    seen = counts[counts["n"] > 0]
    first_seen = seen.groupby("substance")["quarter"].min()
    last_seen = seen.groupby("substance")["quarter"].max()

    # EB05 measures growth in *share*. When the denominator moves, share and
    # count disagree: LA's total OD deaths fell 22% between the current
    # windows, so a substance with flat counts gains share for free, and
    # cocaine posts credible share growth while its own deaths fall. Carry the
    # count ratio alongside so the two are never confused — `count_ratio < 1`
    # with `eb05 > 1` reads "a growing share of a shrinking total", which is a
    # true statement about supply composition and a false one about deaths.
    count_ratio = n_recent / (n_base * per_q).replace(0, np.nan)

    res = pd.DataFrame({
        "n_recent": n_recent.astype(int),
        "n_baseline": n_base.astype(int),
        "expected": expected,
        "share_recent_pct": 100 * n_recent / n_recent_all,
        "share_baseline_pct": 100 * n_base / n_base_all,
        "count_ratio": count_ratio,
        "ebgm": ebgm,
        "eb05": eb05,
        "eb95": eb95,
        "expected_own": expected_own,
        "ebgm_own": ebgm_own,
        "eb05_own": eb05_own,
        "eb95_own": eb95_own,
        "z_spatial": zs,
        "concentration": conc,
    })
    # "Credible rise" needs both gates: share up (Gate A) *and* the
    # substance's own count above its own trend (Gate B). Cocaine-shaped
    # cases -- share credibly up, own count flat or down -- clear A alone and
    # land at `credible_rise == False`, which is the label
    # `regime()` currently requires a reader to work out by hand.
    res["credible_rise"] = (res["eb05"] > ALARM_THRESHOLD) & (
        res["eb05_own"] > ALARM_THRESHOLD)
    res["n_total"] = wide.sum(axis=1).astype(int)
    res["first_seen"] = first_seen
    res["last_seen"] = last_seen
    res = res.sort_values("eb05", ascending=False)

    meta = {
        "prior": prior,
        "prior_own": prior_own,
        "phi": phi,
        "phi_s": phi_s,
        "spatial_mode": spatial_mode,
        "spatial_weight": spatial_weight if spatial_mode == "discount" else None,
        "n_spatial_scored": int((zs != 0).sum()),
        "recent_window": (recent_q[0], recent_q[-1]),
        "baseline_window": (base_q[0], base_q[-1]),
        "n_recent_all": n_recent_all, "n_base_all": n_base_all,
        "truncated_from": (sorted(denom.index)[-lag_quarters]
                           if lag_quarters else None),
        "quarters": quarters,
    }
    return res, meta


# --------------------------------------------------------------------------
# Conditional-Poisson diagnostic
# --------------------------------------------------------------------------
#
# GPS assumes n | lambda ~ Poisson(lambda * E) *within* a substance. The
# Gamma prior models spread *between* substances, so the observed table being
# overdispersed is expected and is not evidence against anything -- it is the
# signal `_fit_gps_prior` exists to measure. What breaks the closed form is
# overdispersion in the conditional: deaths that do not arrive independently.
#
# The cross-section cannot see that (one (n, E) pair per substance carries no
# information about within-substance arrival), but the quarterly series can,
# with one caveat that dominates everything: **a substance whose rate genuinely
# moved also fails a flat Poisson test**, because lambda really did change.
# That is the signal, not a defect. So the statistic that matters is dispersion
# around a *fitted trend*, not around a flat line. On the real data the gap
# between the two is the whole diagnostic -- fentanyl falls 7.68 -> 2.19 once
# its decline is removed and PCP falls 2.05 -> 0.86, while para-fluorofentanyl
# holds at 10.27 -> 10.55 because its counts arrive in batches rather than at a
# rate. See `docs/GAMMA_POISSON_PRIMER.md` step 3.
#
# This reports; it does not correct. Applying a per-substance phi to EB05 is a
# defensible next step (`_credit_dispersion` already does the pooled version
# for credit weighting) but it changes the ranking, so it wants its own
# backtest rather than being switched on inside a diagnostic.

# Below this many mentions across the window, chi-square is too far from its
# asymptotic distribution to read -- most cells would have expected counts
# under 5.
MIN_TOTAL_FOR_DISPERSION = 40

# D/df below this reads as materially smoother than Poisson. Coarse on
# purpose, in the spirit of `polysubstance._verdict`: separate the obvious
# cases, leave the rest unlabelled rather than pretending to resolve them.
DISPERSION_SMOOTH = 0.6


def _poisson_trend_mu(obs: np.ndarray, offset: np.ndarray) -> np.ndarray:
    """Fitted means from ML on log(mu_q) = a + b*t + log(offset_q).

    The offset is the quarter's all-OD death count, so `a + b*t` is a trend in
    the substance's *share*. Fitting the raw count instead would read the
    county's own 22% decline as a substance-level trend and absorb it.

    Nelder-Mead on two parameters, matching `_fit_gps_prior` rather than
    pulling in a GLM dependency; the Poisson log-likelihood is convex here, so
    the method only has to be reliable, not clever. `t` is centred and scaled
    because an uncentred quarter index makes the intercept and slope strongly
    correlated and slows the simplex badly.
    """
    t = np.arange(len(obs), dtype=float)
    t = (t - t.mean()) / max(t.std(), 1e-9)

    def nll(p: np.ndarray) -> float:
        lin = p[0] + p[1] * t
        if not np.isfinite(lin).all() or lin.max() > 50:
            return 1e12
        mu = offset * np.exp(lin)
        if (mu <= 0).any():
            return 1e12
        return float(np.sum(mu - obs * np.log(mu)))

    a0 = np.log(max(obs.sum(), 0.5) / max(offset.sum(), 1e-9))
    res = minimize(nll, np.array([a0, 0.0]), method="Nelder-Mead",
                   options={"maxiter": 5000, "maxfev": 5000,
                            "xatol": 1e-10, "fatol": 1e-10})
    return offset * np.exp(res.x[0] + res.x[1] * t)


def _dispersion_index(obs: np.ndarray, exp: np.ndarray,
                      n_params: int) -> tuple[float, float]:
    """Pearson dispersion D/df, and the chi-square tail probability of D.

    D/df == 1 is exactly Poisson; the p-value is one-sided against the
    overdispersed alternative, which is the direction that costs false alarms.
    """
    df = len(obs) - n_params
    stat = float(((obs - exp) ** 2 / np.maximum(exp, 1e-9)).sum())
    return stat / df, float(chi2.sf(stat, df))


def _dispersion_verdict(d_trend: float, p_trend: float) -> str:
    """Coarse label for whether a substance's own counts obey Poisson."""
    if p_trend < 0.05:
        return "clustered"
    if d_trend < DISPERSION_SMOOTH:
        return "smooth"
    return "ok"


def poisson_dispersion(
    counts: pd.DataFrame, denom: pd.Series,
    quarters: list[pd.Timestamp] | None = None,
    min_total: int = MIN_TOTAL_FOR_DISPERSION,
) -> pd.DataFrame:
    """Per-substance check that quarterly counts vary the way Poisson demands.

    Returns one row per substance with enough mentions to test, carrying:

    `d_flat`   dispersion around a constant share -- inflated by any real
               trend, so high values here mean little on their own.
    `d_trend`  dispersion around a fitted log-linear share trend. This is the
               one to read: it asks whether the counts scatter more than
               Poisson allows *after* granting the substance whatever rise or
               fall it actually had.
    `absorbed` d_flat - d_trend, i.e. how much of the raw excess was trend.
    `verdict`  `clustered` (Poisson rejected -- EB05's interval is too narrow
               for this substance), `smooth` (underdispersed, so EB05 is
               conservative), or `ok`.

    A `clustered` verdict does not invalidate the substance's point estimate;
    lambda = n/E is unaffected. What it invalidates is the *confidence*: the
    interval should be roughly sqrt(d_trend) times wider.
    """
    quarters = sorted(denom.index) if quarters is None else list(quarters)
    wide = counts.pivot_table(index="substance", columns="quarter", values="n",
                              aggfunc="sum", fill_value=0)
    wide = wide.reindex(columns=quarters, fill_value=0)
    offset = denom.loc[quarters].to_numpy(dtype=float)

    rows = []
    for substance, series in wide.iterrows():
        obs = series.to_numpy(dtype=float)
        if obs.sum() < min_total:
            continue
        flat = obs.sum() * offset / offset.sum()
        d_flat, _ = _dispersion_index(obs, flat, 1)
        d_trend, p_trend = _dispersion_index(
            obs, _poisson_trend_mu(obs, offset), 2)
        rows.append({
            "substance": substance,
            "n_window": int(obs.sum()),
            "d_flat": d_flat,
            "d_trend": d_trend,
            "absorbed": d_flat - d_trend,
            "p_trend": p_trend,
            "interval_scale": float(np.sqrt(max(d_trend, 0.0))),
            "verdict": _dispersion_verdict(d_trend, p_trend),
        })

    if not rows:
        return pd.DataFrame(
            columns=["n_window", "d_flat", "d_trend", "absorbed", "p_trend",
                     "interval_scale", "verdict"],
            index=pd.Index([], name="substance"))
    return (pd.DataFrame(rows).set_index("substance")
            .sort_values("d_trend", ascending=False))


def _attach_dispersion(res: pd.DataFrame, counts: pd.DataFrame,
                       denom: pd.Series, meta: dict) -> pd.DataFrame:
    """Join the conditional-Poisson verdict onto a ranking.

    Uses exactly the quarters the ranking scored (baseline start through
    recent end), so the diagnostic describes the same data the EB05 came from.
    Substances below `MIN_TOTAL_FOR_DISPERSION` get NaN/`too_few` rather than a
    number the sample cannot support.
    """
    lo = meta["baseline_window"][0]
    hi = meta["recent_window"][1]
    window = [q for q in sorted(denom.index) if lo <= q <= hi]
    d = poisson_dispersion(counts, denom, window)
    out = res.join(d[["d_trend", "interval_scale", "verdict"]], how="left")
    return out.rename(columns={"verdict": "poisson"}).fillna({"poisson": "too_few"})


def _attach_novelty(res: pd.DataFrame, mentions_path: Path) -> pd.DataFrame:
    """Join NFLIS catalog novelty onto the ranking."""
    m = pd.read_parquet(mentions_path)
    key = m["rollup"].fillna(m["canonical"])
    nov = (m.assign(_k=key)
             .drop_duplicates(subset="_k")
             .set_index("_k")[["date_added", "date_added_censored", "category"]])
    return res.join(nov, how="left")


def _fmt_prior(pr: dict) -> str:
    return (f"GPS prior: Gamma(alpha={pr['alpha']:.3f}, beta={pr['beta']:.3f}), "
            f"mean {pr['alpha'] / pr['beta']:.3f}, "
            f"sd {np.sqrt(pr['alpha']) / pr['beta']:.3f}, "
            f"fitted on {pr['n_fit']} substances with a non-zero baseline"
            + ("" if pr.get("converged", True) else "  [DID NOT CONVERGE]"))


def _fmt_window(meta: dict) -> str:
    (r0, r1), (b0, b1) = meta["recent_window"], meta["baseline_window"]
    q = lambda t: f"{t.year}Q{t.quarter}"  # noqa: E731
    return (f"recent {q(r0)}–{q(r1)} vs baseline {q(b0)}–{q(b1)}; "
            f"{meta['n_recent_all']:,.0f} vs {meta['n_base_all']:,.0f} OD deaths")


def _treescan_ri_live(
    mentions_path: Path, as_of: pd.Timestamp, cutoff: str | None = CUTOFF,
    reference: str = "deaths",
) -> pd.Series:
    """This quarter's TreeScan substance-leaf recurrence interval, live.

    A local import: `treescan.py` imports from this module, so importing it
    back at module load time would be circular. ~5s for one quarter's
    9999-replicate Monte Carlo null -- fine for an interactive `rank`, and
    orders of magnitude cheaper than the 45-quarter sweep `backtest` uses the
    cache for instead (see `_treescan_ri_history`).
    """
    from emerging.analysis import treescan
    counts, denom, ref, nodes = treescan.prepare(mentions_path, cutoff, reference)
    res = treescan.run_scan(counts, ref, nodes, as_of)
    leaf = res[res["level"] == "substance"]
    return leaf.set_index("node")["recurrence_interval"]


def _treescan_ri_history() -> pd.DataFrame:
    """Cached historical TreeScan substance-leaf RI, one row per (substance,
    as_of) -- `treescan backtest`'s output, reused rather than re-scanning 45
    quarters here. Refresh with `emerging treescan backtest` if the mentions
    data has moved since it was last written."""
    from emerging.analysis import treescan
    return treescan.leaf_sweep()


def _apply_treescan_veto(
    res: pd.DataFrame, ri: pd.Series, threshold: float = TREESCAN_VETO_THRESHOLD,
) -> pd.DataFrame:
    """Suppress `credible_rise` where TreeScan sees no corroborating rise.

    `ri`, indexed by substance, is one quarter's TreeScan leaf recurrence
    interval; a substance absent from it gets TreeScan's own floor (RI=1, no
    excess at all under its case definition -- see `treescan._llr`). Adds
    `treescan_ri` and `treescan_veto` columns and never raises `eb05` or
    `credible_rise` -- see `TREESCAN_VETO_THRESHOLD` for why this beats the
    symmetric ensemble combinations tried in `docs/findings/benchmark.md`.
    """
    out = res.copy()
    r = ri.reindex(out.index, fill_value=1.0)
    out["treescan_ri"] = r
    out["treescan_veto"] = out["credible_rise"] & (r < threshold)
    out["credible_rise"] = out["credible_rise"] & ~out["treescan_veto"]
    return out


@app.command("rank")
def rank(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    lag_quarters: int = typer.Option(LAG_QUARTERS, "--lag-quarters"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    top: int = typer.Option(25, "--top"),
    min_recent: int = typer.Option(2, "--min-recent",
                                   help="Suppress substances below this recent count"),
    weighted: bool = typer.Option(True, "--weighted/--no-weighted",
                                  help="Position-weighted case definition "
                                       "(GPS_V2_DESIGN.md #1) instead of raw "
                                       "mention counts"),
    role_discount: bool = typer.Option(True, "--role-discount/--no-role-discount",
                                       help="With --weighted, add each "
                                            "substance's own phi_s "
                                            "(_role_dispersion) -- targets "
                                            "the Lidocaine residual a 2-line "
                                            "co-occurrence leaves undiscounted"),
    spatial: str = typer.Option("discount", "--spatial",
                                help="Spatial-concentration fusion "
                                     "(GPS_V2_DESIGN.md #3): none, discount, "
                                     "or prior"),
    treescan_veto: bool = typer.Option(
        True, "--treescan-veto/--no-treescan-veto",
        help="Suppress credible_rise where a live TreeScan scan of this "
             "quarter sees no corroborating rise (docs/findings/"
             "benchmark.md's ensemble-veto-v2role-10)"),
    veto_threshold: float = typer.Option(
        TREESCAN_VETO_THRESHOLD, "--veto-threshold",
        help="TreeScan recurrence-interval floor below which a rise is vetoed"),
) -> None:
    """Rank rising substances by EB05 and write the table.

    The defaults are the project's current recommendation
    (`docs/findings/benchmark.md`'s `ensemble-veto-v2role-10`): position-
    weighted + role-discounted case definition, dual-gated, spatial E-
    discount, TreeScan veto. Every piece can be switched off individually to
    reproduce an earlier variant from that comparison, including plain EB05
    (`--no-weighted --no-role-discount --spatial none --no-treescan-veto`).
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff, weighted=weighted,
                                   role_discount=role_discount)
    quarters = sorted(denom.index)
    z_now = None
    if spatial != "none":
        z_now = _z_at(_spatial_sweep(quarters, mentions,
                                     recent_quarters=recent_quarters),
                      quarters[-1 - lag_quarters] if lag_quarters
                      else quarters[-1])
    res, meta = rank_substances(counts, denom, lag_quarters,
                                recent_quarters, baseline_quarters,
                                z_spatial=z_now, spatial_mode=spatial)
    res = _attach_novelty(res, mentions)
    res = _attach_dispersion(res, counts, denom, meta)
    if treescan_veto:
        as_of = quarters[-1 - lag_quarters] if lag_quarters else quarters[-1]
        ri = _treescan_ri_live(mentions, as_of, cutoff)
        res = _apply_treescan_veto(res, ri, veto_threshold)
    if spatial == "none":
        # A column of zeros reads as "no substance is concentrated", which is
        # false -- it means the term was switched off. Drop it rather than
        # publish an inactive statistic that looks like a measured one.
        res = res.drop(columns=["z_spatial", "concentration"])

    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "rising_substances.csv")

    typer.echo(_fmt_window(meta)
               + (" [weighted+role]" if weighted and role_discount else
                  " [weighted]" if weighted else "")
               + (f" [spatial: {spatial}]" if spatial != "none" else "")
               + (" [treescan-veto]" if treescan_veto else ""))
    typer.echo(_fmt_prior(meta["prior"]))
    if weighted:
        typer.echo(f"dispersion correction: phi={meta['phi']:.3f} "
                   "(1.0 = no overdispersion from the credit weighting)")
        if role_discount:
            ps = meta["phi_s"]
            typer.echo(f"role discount: {len(ps)} substances scored "
                       f"(rest below MIN_CASES_FOR_ROLE_DISCOUNT), "
                       + (f"max phi_s={max(ps.values()):.2f} "
                          f"({max(ps, key=ps.get)})" if ps else "none scored"))
    if spatial != "none":
        typer.echo(f"spatial: {meta['n_spatial_scored']} substances scored "
                   f"(rest below concentration.MIN_CASES or withheld)"
                   + (f", gamma={meta['prior']['gamma']:+.3f}"
                      if spatial == "prior" else
                      f", w={meta['spatial_weight']:g}"))
    if meta["truncated_from"] is not None:
        t = meta["truncated_from"]
        typer.echo(f"dropped as incomplete: {t.year}Q{t.quarter} onward "
                   f"({lag_quarters} quarters)")
    if treescan_veto:
        vetoed = res[res["treescan_veto"]]
        typer.echo(f"treescan veto: RI < {veto_threshold:g} suppresses "
                   f"credible_rise; {len(vetoed)} substance(s) vetoed this "
                   "quarter"
                   + (": " + ", ".join(
                       f"{n} (RI={r.treescan_ri:.1f})"
                       for n, r in vetoed.iterrows()) if len(vetoed) else ""))

    show = res[res["n_recent"] >= min_recent].head(top)
    cols = ["n_recent", "n_baseline", "expected", "share_recent_pct",
            "share_baseline_pct", "count_ratio", "ebgm", "eb05", "eb05_own",
            "credible_rise", "poisson", "n_total", "first_seen"]
    if treescan_veto:
        cols = cols[:cols.index("credible_rise") + 1] + ["treescan_ri"] + \
               cols[cols.index("credible_rise") + 1:]
    disp = show[cols].copy()
    disp["first_seen"] = disp["first_seen"].map(
        lambda t: "" if pd.isna(t) else f"{t.year}Q{t.quarter}"
    )
    typer.echo("\n" + disp.to_string(float_format=lambda v: f"{v:.2f}"))

    # `poisson` qualifies the interval, not the estimate: `clustered` means
    # this substance's deaths do not arrive independently, so its EB05 is
    # narrower than the data supports. See `poisson_dispersion`.
    bad = res[res["poisson"] == "clustered"]
    if len(bad):
        typer.echo(f"\nconditional Poisson rejected for {len(bad)} substance(s) "
                   "-- EB05 is overconfident there: "
                   + ", ".join(f"{n} (x{r.interval_scale:.1f})"
                               for n, r in bad.head(6).iterrows()))
    typer.echo(f"\nwrote {out_dir / 'rising_substances.csv'}")


@app.command("dispersion")
def dispersion(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    min_total: int = typer.Option(MIN_TOTAL_FOR_DISPERSION, "--min-total"),
) -> None:
    """Check the assumption EB05 rests on: do each substance's deaths arrive
    independently, the way Poisson requires?

    Read `d_trend`, not `d_flat`. A substance that genuinely rose or fell fails
    a flat Poisson test for the best possible reason -- its rate moved -- so
    the honest question is whether it scatters more than Poisson allows *after*
    granting it that trend. `absorbed` is how much of the raw excess the trend
    explained; a large `absorbed` with a small `d_trend` is a clean substance
    that happens to be moving, which is exactly what the ranking looks for.

    `clustered` does not mean the EB05 is wrong, it means it is overconfident:
    the point estimate lambda = n/E is untouched, but the interval should be
    about `interval_scale` times wider. The mechanism is deaths sharing a
    source -- one contaminated batch killing several people -- which the Gamma
    prior cannot absorb, because that layer models differences *between*
    substances, not clustering inside one.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    quarters = sorted(denom.index)[-(recent_quarters + baseline_quarters):]
    d = poisson_dispersion(counts, denom, quarters, min_total=min_total)

    out_dir.mkdir(parents=True, exist_ok=True)
    d.to_csv(out_dir / "poisson_dispersion.csv")

    lo, hi = quarters[0], quarters[-1]
    typer.echo(f"window {lo.year}Q{lo.quarter}-{hi.year}Q{hi.quarter} "
               f"({len(quarters)} quarters); {len(d)} substances with "
               f">= {min_total} mentions")
    typer.echo("D/df == 1 is exactly Poisson; d_trend is the one to read.\n")
    typer.echo(d.to_string(float_format=lambda v: f"{v:.2f}"))

    n_bad = int((d["verdict"] == "clustered").sum())
    typer.echo(f"\n{n_bad} of {len(d)} reject the conditional Poisson "
               f"(p < 0.05): their EB05 intervals are too narrow.")
    typer.echo(f"wrote {out_dir / 'poisson_dispersion.csv'}")


@app.command("backtest")
def backtest(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    substances: str = typer.Option(",".join(KNOWN_EMERGENCES), "--substances"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
    weighted: bool = typer.Option(True, "--weighted/--no-weighted",
                                  help="Position-weighted case definition "
                                       "(GPS_V2_DESIGN.md #1) instead of raw "
                                       "mention counts"),
    role_discount: bool = typer.Option(True, "--role-discount/--no-role-discount",
                                       help="With --weighted, add each "
                                            "substance's own phi_s "
                                            "(_role_dispersion)"),
    spatial: str = typer.Option("discount", "--spatial",
                                help="Spatial-concentration fusion "
                                     "(GPS_V2_DESIGN.md #3): none, discount, "
                                     "or prior"),
    treescan_veto: bool = typer.Option(
        True, "--treescan-veto/--no-treescan-veto",
        help="Suppress credible_rise where TreeScan's cached historical "
             "scan (`emerging treescan backtest`) sees no corroborating "
             "rise that quarter"),
    veto_threshold: float = typer.Option(
        TREESCAN_VETO_THRESHOLD, "--veto-threshold",
        help="TreeScan recurrence-interval floor below which a rise is vetoed"),
) -> None:
    """Re-run the ranking as of each historical quarter and report where the
    known emergences placed.

    Without this the ranking is unfalsifiable: a detector that never fires and
    a detector that is broken produce the same empty top-10. Rolling the
    as-of date back to the quarter each substance actually emerged answers
    "would this have caught it, and how late" — and, for LA County, mostly
    answers "there was never enough signal to catch."

    Also the validation gate for `--weighted` and the dual-gate `eb05_own` /
    `credible_rise` columns (`docs/GPS_V2_DESIGN.md`): a change that moves a
    known emergence to a *later* detection quarter than plain, unweighted
    single-gate EB05 is a regression in this table, not an improvement. The
    defaults are the project's current recommendation (`docs/findings/
    benchmark.md`'s `ensemble-veto-v2role-10`); see `rank`'s docstring.

    **`--weighted` recomputes at every as-of step instead of slicing.**
    `_ever_independent` is a cross-case aggregate (see `_position_credit`),
    so reusing counts built from the full `cutoff`-bounded history at an
    earlier as-of quarter would let a substance's *later* independence
    evidence retroactively un-discount its *earlier* trailing mentions --
    exactly the leak the as-of discipline here exists to prevent. The
    unweighted path has no such hazard (each death's credit only ever reads
    its own cause-line text) and keeps the cheap slice-after-the-fact it
    always used.

    **`--treescan-veto` reads the cached sweep, not a live scan per as-of
    quarter.** A live scan is ~4.5s/quarter (`_treescan_ri_live`); over a
    45-quarter sweep that's minutes, for a check this module's own docstring
    already asks to be re-run often. `emerging treescan backtest` writes the
    cache this reads; refresh it there if the mentions data has moved.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff, weighted=weighted)
    targets = [s.strip() for s in substances.split(",") if s.strip()]
    all_q = sorted(denom.index)
    total_q = recent_quarters + baseline_quarters
    z_sweep = (_spatial_sweep(all_q, mentions, recent_quarters=recent_quarters)
               if spatial != "none" else None)
    ri_hist = (_treescan_ri_history().set_index(["substance", "as_of"])
              ["recurrence_interval"] if treescan_veto else None)

    rows = []
    # Step the as-of quarter forward; at each step score only data available
    # then, so nothing after the as-of date can leak into the statistic.
    for i in range(total_q, len(all_q) + 1):
        as_of = all_q[i - 1]
        sub_denom = denom.loc[all_q[:i]]
        if weighted:
            sub_counts, _ = load_quarterly(
                mentions, cutoff=as_of + pd.offsets.QuarterEnd(0),
                weighted=True, quiet=True, role_discount=role_discount)
        else:
            sub_counts = counts[counts["quarter"] <= as_of]
        try:
            res, _ = rank_substances(
                sub_counts, sub_denom, lag_quarters=0,
                recent_quarters=recent_quarters,
                baseline_quarters=baseline_quarters,
                z_spatial=_z_at(z_sweep, as_of) if z_sweep is not None else None,
                spatial_mode=spatial)
        except ValueError:
            continue
        ranked = res[res["n_recent"] > 0].sort_values("eb05", ascending=False)
        order = {name: r for r, name in enumerate(ranked.index, start=1)}
        for s in targets:
            if s not in res.index:
                continue
            r = res.loc[s]
            credible_rise = bool(r["credible_rise"])
            ri = np.nan
            if treescan_veto:
                ri = float(ri_hist.get((s, as_of), 1.0))
                credible_rise = credible_rise and not (ri < veto_threshold)
            rows.append({
                "as_of": as_of, "substance": s, "rank": order.get(s),
                "n_recent": int(r["n_recent"]), "n_baseline": int(r["n_baseline"]),
                "ebgm": r["ebgm"], "eb05": r["eb05"],
                "eb05_own": r["eb05_own"], "treescan_ri": ri,
                "credible_rise": credible_rise,
            })

    bt = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    bt.to_csv(out_dir / "backtest_known_emergences.csv", index=False)

    typer.echo(f"as-of sweep: {len(all_q) - total_q + 1} quarters, "
               f"{recent_quarters}q recent vs {baseline_quarters}q baseline"
               + (" [weighted+role]" if weighted and role_discount else
                  " [weighted]" if weighted else "")
               + (f" [spatial: {spatial}]" if spatial != "none" else "")
               + (" [treescan-veto]" if treescan_veto else "") + "\n")
    for s in targets:
        g = bt[bt["substance"] == s]
        if g.empty or g["n_recent"].max() == 0:
            typer.echo(f"{s:24s} never reaches a scoreable count")
            continue
        best = g.loc[g["eb05"].idxmax()]
        top10 = g[g["rank"] <= 10]
        first_top10 = top10["as_of"].min() if not top10.empty else None
        cr = g[g["credible_rise"]]
        first_cr = cr["as_of"].min() if not cr.empty else None
        typer.echo(
            f"{s:24s} peak EB05 {best['eb05']:.2f} (EBGM {best['ebgm']:.2f}, "
            f"n={int(best['n_recent'])}) at "
            f"{best['as_of'].year}Q{best['as_of'].quarter}, best rank "
            f"{int(g['rank'].min())}"
            + (f", first top-10 at {first_top10.year}Q{first_top10.quarter}"
               if first_top10 is not None else ", never reaches top 10")
            + (f", dual gate first clears at "
               f"{first_cr.year}Q{first_cr.quarter}"
               if first_cr is not None else ", dual gate never clears")
        )
    typer.echo(f"\nwrote {out_dir / 'backtest_known_emergences.csv'}")


@app.command("regime")
def regime(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
) -> None:
    """Look for changes in how deaths are recorded, rather than in what people
    are dying of, that would move EB05.

    The statistic is a *share*: numerator the substance, denominator all
    overdose deaths, both counted from ME narratives. Three things can move it
    without any change in the drug supply, and all three are measurable here:

    1. **Denominator drift.** Total OD deaths falling makes every flat
       substance gain share. This is the live one — see `count_ratio` in the
       ranking.
    2. **Panel drift.** Substances named per death rising means the ME is
       finding more of everything, which inflates the numerator of the rare
       substances specifically.
    3. **Simultaneous onsets.** Several substances first breaching in the same
       quarter is the signature of a change in recording, not of several
       independent epidemics starting at once.

    What this cannot do is identify the *policy*. It flags the quarters where
    something structural happened; matching those to naloxone distribution,
    scheduling changes, or a lab protocol revision needs sources outside this
    dataset.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    all_q = sorted(denom.index)
    spd = (counts.groupby("quarter")["n"].sum().reindex(all_q, fill_value=0)
           / denom.reindex(all_q))

    tab = pd.DataFrame({
        "deaths": denom.reindex(all_q).astype(int),
        "yoy_pct": 100 * (denom.reindex(all_q) / denom.reindex(all_q).shift(4) - 1),
        "subs_per_death": spd,
        "n_substances": counts[counts["n"] > 0].groupby("quarter").size()
                        .reindex(all_q, fill_value=0),
    })

    onsets = pd.Series(0, index=all_q, dtype=int)
    if ALARM_PATH.exists():
        h = pd.read_csv(ALARM_PATH, parse_dates=["first_breach"])
        c = h["first_breach"].value_counts()
        onsets = onsets.add(c.reindex(all_q, fill_value=0), fill_value=0)
    tab["breach_onsets"] = onsets.astype(int)

    out_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out_dir / "regime_diagnostics.csv")

    typer.echo("Quarterly recording diagnostics "
               f"({all_q[0].year}Q{all_q[0].quarter}–"
               f"{all_q[-1].year}Q{all_q[-1].quarter})\n")
    d = tab.copy()
    d.index = [f"{q.year}Q{q.quarter}" for q in d.index]
    typer.echo(d.tail(16).to_string(float_format=lambda v: f"{v:.2f}"))

    peak_q = denom.idxmax()
    last4 = float(denom.reindex(all_q[-4:]).sum())
    prev4 = float(denom.reindex(all_q[-8:-4]).sum())
    typer.echo(
        f"\n1. Denominator. Peak {int(denom.max())} deaths in "
        f"{peak_q.year}Q{peak_q.quarter}; last 4 quarters {last4:,.0f} vs the "
        f"4 before {prev4:,.0f} ({100 * (last4 / prev4 - 1):+.0f}%). "
        "Every substance holding a flat count gains "
        f"x{prev4 / last4:.2f} share on this alone.")
    typer.echo(
        f"2. Panel. Substances named per death {spd.iloc[0]:.2f} → "
        f"{spd.iloc[-1]:.2f}; last 12 quarters mean "
        f"{spd.iloc[-12:].mean():.2f}, sd {spd.iloc[-12:].std():.3f} — "
        + ("flat, so recent flags are not a widening panel."
           if spd.iloc[-12:].std() < 0.06 else
           "still moving; treat recent flags with caution."))

    clusters = tab[tab["breach_onsets"] >= 3]
    if len(clusters):
        typer.echo(f"3. Simultaneous onsets (>=3 substances first breaching "
                   f"EB05 {threshold:g} in one quarter):")
        for q, r in clusters.iterrows():
            names = ", ".join(sorted(
                h.loc[h["first_breach"] == q, "substance"]))
            typer.echo(f"     {q.year}Q{q.quarter}  {int(r.breach_onsets)} "
                       f"substances, deaths {int(r.deaths)} "
                       f"({r.yoy_pct:+.0f}% YoY): {names}")
        typer.echo("   Several unrelated substances starting at once is more "
                   "consistent with a change in\n   the denominator or in "
                   "recording than with simultaneous independent epidemics.")
    else:
        typer.echo("3. No quarter has >=3 simultaneous first breaches.")
    typer.echo(f"\nwrote {out_dir / 'regime_diagnostics.csv'}")


@app.command("watch")
def watch(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
) -> None:
    """Check the national emerging-substance watchlist against LA County.

    Answers "is what the rest of the country is seeing showing up here", which
    the ranking cannot: a substance with zero LA mentions has no rank. Reports
    every watchlist entry including the zeros, because absence is the finding
    for most of them.

    **A zero is ambiguous and the output says so.** The lexicon recognises all
    of these names, so a zero means the narratives never mention the
    substance. It does not mean the substance is not in the LA supply — an
    assay the ME does not run produces exactly the same zero. Separating the
    two needs the ME's panel, not this data.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff)
    res, meta = rank_substances(counts, denom, 0, recent_quarters,
                                baseline_quarters)
    m = pd.read_parquet(mentions)
    m = m[~m["is_class_term"].fillna(False)]
    key = m["rollup"].fillna(m["canonical"])

    fam = sorted(set(key[key.str.contains(NITAZENE_RE, case=False, na=False)]))
    targets = list(WATCHLIST.items()) + [
        (f"[family] nitazenes ({len(fam)} analogs)",
         "novel synthetic opioids, individually too rare to score")]

    hist = (pd.read_csv(ALARM_PATH, index_col=0, parse_dates=["first_breach"])
            if ALARM_PATH.exists() else None)

    rows = []
    for name, note in targets:
        members = fam if name.startswith("[family]") else [name]
        sub = counts[counts["substance"].isin(members)]
        total = int(sub["n"].sum())
        r = res.loc[[x for x in members if x in res.index]]
        fired = (hist is not None
                 and any(x in hist.index for x in members))
        rows.append({
            "watch": name, "note": note, "n_total": total,
            "first_seen": sub.loc[sub["n"] > 0, "quarter"].min()
                          if total else pd.NaT,
            "last_seen": sub.loc[sub["n"] > 0, "quarter"].max()
                         if total else pd.NaT,
            "n_recent": int(r["n_recent"].sum()) if len(r) else 0,
            "eb05_max": float(r["eb05"].max()) if len(r) else np.nan,
            "has_fired": fired,
        })

    tab = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out_dir / "watchlist.csv", index=False)

    q = lambda t: "—" if pd.isna(t) else f"{t.year}Q{t.quarter}"  # noqa: E731
    typer.echo(_fmt_window(meta) + "\n")
    for r in rows:
        if r["n_total"] == 0:
            typer.echo(f"  {r['watch'][:52]:54s} ABSENT — 0 mentions, "
                       f"2012-{pd.Timestamp(cutoff).year}")
        else:
            typer.echo(
                f"  {r['watch'][:52]:54s} n={r['n_total']:<5d} "
                f"{q(r['first_seen'])}→{q(r['last_seen'])}  "
                f"recent={r['n_recent']:<4d} EB05={r['eb05_max']:.2f}"
                + ("  [has fired]" if r["has_fired"] else ""))
        typer.echo(f"  {'':54s} {r['note']}")
    typer.echo("\nA zero means the narratives never name it. It cannot be "
               "distinguished from\na substance the ME's panel does not "
               "assay for — that needs the panel, not this data.")
    typer.echo(f"\nwrote {out_dir / 'watchlist.csv'}")


@app.command("alarms")
def alarms(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
    weighted: bool = typer.Option(True, "--weighted/--no-weighted",
                                  help="Position-weighted case definition "
                                       "(GPS_V2_DESIGN.md #1) instead of raw "
                                       "mention counts"),
    role_discount: bool = typer.Option(True, "--role-discount/--no-role-discount",
                                       help="With --weighted, add each "
                                            "substance's own phi_s "
                                            "(_role_dispersion)"),
    spatial: str = typer.Option("discount", "--spatial",
                                help="Spatial-concentration fusion "
                                     "(GPS_V2_DESIGN.md #3): none, discount, "
                                     "or prior"),
    treescan_veto: bool = typer.Option(
        True, "--treescan-veto/--no-treescan-veto",
        help="Suppress credible_rise where TreeScan's cached historical "
             "scan (`emerging treescan backtest`) sees no corroborating "
             "rise that quarter"),
    veto_threshold: float = typer.Option(
        TREESCAN_VETO_THRESHOLD, "--veto-threshold",
        help="TreeScan recurrence-interval floor below which a rise is vetoed"),
) -> None:
    """Sweep every substance across every as-of quarter and record when each
    one crossed the alarm line.

    `backtest` asks "would the detector have caught these five". This asks the
    complementary question — "what has it ever fired on at all" — which is the
    only way to see the false-alarm rate and the clustering of onsets. It
    writes both the full sweep (for plotting) and the per-substance summary.

    The defaults match `rank`/`backtest`'s current recommendation, and "in
    alarm" here means `credible_rise` (the dual gate, TreeScan-vetoed) rather
    than the looser single-gate `eb05 > threshold` -- the same episodes the
    deployed detector would actually have raised. `--threshold` still governs
    the dual gate's own two `eb05 > threshold` tests (see `rank_substances`),
    just no longer as a *separate*, single-gate alarm definition.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff, weighted=weighted,
                                   role_discount=role_discount)
    all_q = sorted(denom.index)
    z_sweep = (_spatial_sweep(all_q, mentions, recent_quarters=recent_quarters)
               if spatial != "none" else None)
    sw = sweep_eb05(counts, denom, recent_quarters, baseline_quarters,
                    z_sweep=z_sweep, spatial_mode=spatial,
                    weighted=weighted, role_discount=role_discount,
                    mentions_path=mentions)
    if treescan_veto:
        ri_hist = (_treescan_ri_history()
                  .set_index(["substance", "as_of"])["recurrence_interval"])
        sw = sw.set_index(["substance", "as_of"])
        sw["treescan_ri"] = ri_hist.reindex(sw.index, fill_value=1.0)
        sw["credible_rise"] = sw["credible_rise"] & (sw["treescan_ri"] >= veto_threshold)
        sw = sw.reset_index()

    seen = counts[counts["n"] > 0].groupby("substance")["quarter"].min()
    hist = alarm_history(sw, seen, max(denom.index), threshold,
                         first_q=min(denom.index), gate_col="credible_rise",
                         baseline_quarters=baseline_quarters)

    out_dir.mkdir(parents=True, exist_ok=True)
    sw.to_csv(out_dir / "eb05_sweep.csv", index=False)
    hist.to_csv(out_dir / ALARM_PATH.name)

    q = lambda t: "" if pd.isna(t) else f"{t.year}Q{t.quarter}"  # noqa: E731
    typer.echo(f"swept {len(all_q) - recent_quarters - baseline_quarters + 1} "
               f"as-of quarters over {sw['substance'].nunique()} substances; "
               f"{len(hist)} ever fired EB05+"
               + (" [weighted+role]" if weighted and role_discount else
                  " [weighted]" if weighted else "")
               + (f" [spatial: {spatial}]" if spatial != "none" else "")
               + (" [treescan-veto]" if treescan_veto else "") + "\n")
    disp = hist.copy()
    for c in ("first_seen", "first_breach", "last_breach", "peak_as_of",
              "latest_onset"):
        disp[c] = disp[c].map(q)
    cols = ["first_seen", "first_breach", "lag_quarters", "last_breach",
            "quarters_in_alarm", "n_episodes", "peak_eb05", "n_at_peak",
            "still_open", "latest_onset", "avg_baseline_deaths_q", "emergent"]
    typer.echo(disp[cols].to_string(float_format=lambda v: f"{v:.2f}"))
    typer.echo(f"\nwrote {out_dir / ALARM_PATH.name} "
               f"and {out_dir / 'eb05_sweep.csv'}")


SWEEP_COLUMNS = ["n_recent", "n_baseline", "expected", "ebgm", "eb05",
                 "eb05_own", "credible_rise", "z_spatial"]


def sweep_eb05(
    counts: pd.DataFrame, denom: pd.Series,
    recent_quarters: int = RECENT_QUARTERS,
    baseline_quarters: int = BASELINE_QUARTERS,
    z_sweep: pd.DataFrame | None = None,
    spatial_mode: str = "none",
    spatial_weight: float = SPATIAL_WEIGHT,
    weighted: bool = False,
    role_discount: bool = False,
    mentions_path: Path = MENTIONS_PATH,
) -> pd.DataFrame:
    """Re-score *every* substance at every as-of quarter.

    Same as-of discipline as `backtest`, but unrestricted to the known set:
    at each step the prior is refitted and the scores computed from data
    available then, so nothing after the as-of date leaks in. Returns a long
    frame of `SWEEP_COLUMNS` plus `substance` and `as_of` — enough for every
    model variant `validation/benchmark.py` scores to be read off one sweep
    without re-deriving anything: `expected` carries the unshrunken n/E, the
    two `eb05` columns carry the single and dual gates.

    `spatial_mode` / `z_sweep` apply `docs/GPS_V2_DESIGN.md` #3, and
    `weighted` (with `role_discount`) applies #1 -- with the same
    recompute-don't-slice discipline `backtest` documents, because both
    `_ever_independent` and `_role_dispersion` are cross-case aggregates and
    slicing either would leak a substance's later evidence backward into its
    earlier scores.

    This is the expensive command in the module — it refits the GPS prior once
    per quarter — so it is a separate CLI whose output the plots read from
    disk rather than something `plot` recomputes.
    """
    all_q = sorted(denom.index)
    total_q = recent_quarters + baseline_quarters
    frames = []
    for i in range(total_q, len(all_q) + 1):
        as_of = all_q[i - 1]
        if weighted:
            sub_counts, _ = load_quarterly(
                mentions_path, cutoff=as_of + pd.offsets.QuarterEnd(0),
                weighted=True, quiet=True, role_discount=role_discount)
        else:
            sub_counts = counts[counts["quarter"] <= as_of]
        res, _ = rank_substances(
            sub_counts, denom.loc[all_q[:i]],
            lag_quarters=0, recent_quarters=recent_quarters,
            baseline_quarters=baseline_quarters,
            z_spatial=_z_at(z_sweep, as_of) if z_sweep is not None else None,
            spatial_mode=spatial_mode, spatial_weight=spatial_weight)
        r = res.loc[res["n_recent"] > 0, SWEEP_COLUMNS].copy()
        r["as_of"] = as_of
        frames.append(r.reset_index().rename(columns={"index": "substance"}))
    return pd.concat(frames, ignore_index=True)


def _episodes(quarters: list[pd.Timestamp]) -> list[tuple]:
    """Contiguous runs of alarm quarters -> [(start, end), ...].

    Quarters are contiguous when they are one QuarterBegin apart. A substance
    that breaches, subsides, and breaches again has two episodes, and that is
    reported rather than smoothed into one long span — an intermittent alarm
    and a sustained one mean different things.
    """
    if not quarters:
        return []
    runs, start, prev = [], quarters[0], quarters[0]
    for q in quarters[1:]:
        if q == prev + pd.offsets.QuarterBegin(1, startingMonth=1):
            prev = q
            continue
        runs.append((start, prev))
        start = prev = q
    runs.append((start, prev))
    return runs


def alarm_history(
    sweep: pd.DataFrame, first_seen: pd.Series, last_as_of: pd.Timestamp,
    threshold: float = ALARM_THRESHOLD, first_q: pd.Timestamp | None = None,
    gate_col: str | None = None, baseline_quarters: int = BASELINE_QUARTERS,
    emergent_threshold: float = EMERGENT_THRESHOLD,
) -> pd.DataFrame:
    """Per-substance alarm record: onset, first breach, end, peak, duration.

    The four dates a reader needs to judge a flag, in the order they happen:
    `first_seen` (the first LA death naming it), `first_breach` (when the
    detector would have fired), `last_breach` (when concern ended), and
    `peak_as_of`. `lag_quarters` between the first two is the detection delay
    — the honest measure of how much warning this buys.

    `still_open` marks substances whose last breach is the final as-of
    quarter: their alarm has not ended, it just has not ended *yet*.

    **`lag_quarters` is only meaningful for substances that arrived after the
    extract began.** LACME records start 2012Q1, so anything already
    circulating then reports `first_seen == 2012Q1` — a censoring floor, not an
    arrival. Reading morphine's 53 quarters as a detection delay would be
    nonsense: it was never detected late, it was simply always here and only
    recently grew. Those rows are flagged `first_seen_censored` and their lag
    is NaN, not a number.

    `gate_col`, when given, names a boolean column in `sweep` (e.g.
    `credible_rise`, already dual-gated and optionally TreeScan-vetoed by the
    caller) that defines "in alarm" instead of the legacy single-gate
    `eb05 > threshold`. `peak_eb05`/`peak_as_of` still read off `eb05` either
    way -- the gate changes which quarters count, not what "peak" means.

    **`emergent` is a second, independent axis from "in alarm".** `eb05`
    answers "is this substance's share credibly rising"; `emergent` answers
    the different question "was it previously rare enough that a rise means
    something new", read off `avg_baseline_deaths_q` -- `n_baseline` at
    `latest_onset` (the `baseline_quarters` immediately before the *most
    recent* alarm episode started) divided by `baseline_quarters`.

    Anchored on `latest_onset` rather than `first_breach` on purpose. A
    substance can alarm, subside, and alarm again years later (`n_episodes`
    already tracks this) -- fentanyl's 2016-2021 rise and a hypothetical
    fresh 2025 one would be two episodes of the same row. Reading the
    baseline off the *first* episode would judge a 2025 rise by how rare
    fentanyl was in 2014, which is not the question "is this rise, now,
    something new" is asking. `latest_onset` is the start of `_episodes`'
    last run, so a substance with one episode reads the same baseline either
    way and only a repeat-offender's classification can change between
    episodes.

    Also anchored on a *baseline* window rather than the recent or blended
    count: gating on a window that includes the rise itself would let a
    substance's own growth push it past the count threshold and silently
    un-label it as emergent right when the rise is most newsworthy. Every
    substance can be both rising and established (fentanyl, PCP -- already
    killing 3+/quarter before their most recent rise) or rising and emergent
    (xylazine, para-fluorofentanyl -- at or near zero beforehand); `emergent`
    is an annotation on top of `credible_rise`/`eb05 > threshold`, not a
    separate detector, and rows that never breach never get one.
    """
    rows = []
    for sub, g in sweep.groupby("substance"):
        g = g.sort_values("as_of")
        hot = g[g[gate_col]] if gate_col else g[g["eb05"] > threshold]
        if hot.empty:
            continue
        eps = _episodes(list(hot["as_of"]))
        peak = g.loc[g["eb05"].idxmax()]
        fs = first_seen.get(sub, pd.NaT)
        first_breach = hot["as_of"].min()
        latest_onset = eps[-1][0]
        censored = bool(first_q is not None and pd.notna(fs) and fs <= first_q)
        lag = ((first_breach.year - fs.year) * 4 + first_breach.quarter
               - fs.quarter) if pd.notna(fs) else np.nan
        avg_baseline_q = (float(g.loc[g["as_of"] == latest_onset, "n_baseline"]
                                .iloc[0]) / baseline_quarters)
        rows.append({
            "substance": sub,
            "first_seen": fs,
            "first_seen_censored": censored,
            "first_breach": first_breach,
            "last_breach": hot["as_of"].max(),
            "lag_quarters": np.nan if censored else lag,
            "quarters_in_alarm": int(len(hot)),
            "n_episodes": len(eps),
            "longest_episode_q": max(
                (e - s).n // 3 + 1 for s, e in
                [(pd.Period(a, "Q"), pd.Period(b, "Q")) for a, b in eps]),
            "peak_eb05": float(peak["eb05"]),
            "peak_ebgm": float(peak["ebgm"]),
            "peak_as_of": peak["as_of"],
            "n_at_peak": int(peak["n_recent"]),
            "still_open": bool(hot["as_of"].max() == last_as_of),
            "latest_onset": latest_onset,
            "avg_baseline_deaths_q": avg_baseline_q,
            "emergent": avg_baseline_q < emergent_threshold,
        })
    return (pd.DataFrame(rows).set_index("substance")
            .sort_values("first_breach"))


def _episode_spans(sweep: pd.DataFrame, sub: str,
                   threshold: float, gate_col: str | None = None) -> list[tuple]:
    g = sweep[sweep["substance"] == sub].sort_values("as_of")
    mask = g[gate_col] if gate_col else (g["eb05"] > threshold)
    return _episodes(list(g.loc[mask, "as_of"]))


def _quarterly_series(counts: pd.DataFrame, substance: str,
                      quarters: list) -> pd.Series:
    s = (counts[counts["substance"] == substance]
         .set_index("quarter")["n"].reindex(quarters, fill_value=0))
    return s


def _bold_title(name: str) -> str:
    """Bold a substance name for a matplotlib panel title via inline
    mathtext (`ax.set_title` renders a `$...$` span as math and leaves the
    rest of the string as plain text, so this can sit mid-title). A bare
    space or hyphen inside `\\mathbf` gets mathtext's math spacing -- a
    dropped gap, or a wide padded minus sign -- rather than reading as
    ordinary text, so both are escaped."""
    escaped = name.replace(" ", r"\ ").replace("-", r"\text{-}")
    return r"$\mathbf{" + escaped + "}$"


def _panel(ax, x, y, title: str, color: str, first_seen=None,
           trunc_from=None, annotate_last: bool = True,
           alarm_spans: list | None = None, badge: str | None = None) -> None:
    """One small-multiple panel.

    Complete quarters are drawn solid; the right-truncated tail is drawn muted
    and dashed rather than dropped, so the reader can see that the apparent
    cliff at the right edge is missing toxicology and not a real collapse.

    `alarm_spans` shades the quarters in which EB05 exceeded the alarm line, so
    the detector's verdict is visible against the raw counts that produced it —
    the point being that the shading starts on the *upslope*, not at the peak.

    `badge`, when given, draws a small filled rounded-rect label in the top
    right corner -- a CSS-tag look, `emergent` read at a glance rather than
    parsed out of the title text.
    """
    x = list(x)
    n_complete = (sum(1 for q in x if q < trunc_from) if trunc_from is not None
                  else len(x))
    xc, yc = x[:n_complete], y.iloc[:n_complete]

    ax.plot(xc, yc, color=color, linewidth=2, solid_capstyle="round", zorder=3)
    ax.fill_between(xc, 0, yc, color=color, alpha=0.12, linewidth=0, zorder=2)
    if n_complete < len(x):
        # Bridge from the last complete point so the dashed tail connects.
        xt, yt = x[n_complete - 1:], y.iloc[n_complete - 1:]
        ax.plot(xt, yt, color=MUTED, linewidth=1.5, linestyle="--", zorder=3)
        ax.axvspan(x[n_complete], x[-1], color=MUTED, alpha=0.22, linewidth=0,
                   zorder=0)
    # Alarm band. Deliberately unlabelled in-panel: at this facet width the
    # label collided with the neighbouring panel's title. The breach quarter
    # goes in the panel title and the band is explained in the caption.
    if alarm_spans:
        for k, (s, e) in enumerate(alarm_spans):
            ax.axvspan(s, e + pd.offsets.QuarterEnd(0), color=YELLOW,
                       alpha=0.22, linewidth=0, zorder=1)
            if k == 0:
                ax.axvline(s, color=YELLOW, linewidth=1.6, zorder=2)
    # Only mark first appearance when it happens inside the plotted window —
    # for substances present since 2012 the line is just noise on the axis.
    if (first_seen is not None and pd.notna(first_seen)
            and len(x) and first_seen > x[0]):
        ax.axvline(first_seen, color=ORANGE, linewidth=1.2, linestyle=":",
                   zorder=1)
    if annotate_last and len(yc) and float(np.max(yc)) > 0:
        ax.annotate(f"{int(yc.iloc[-1])}", (xc[-1], yc.iloc[-1]),
                    textcoords="offset points", xytext=(4, 2),
                    fontsize=8, color=INK_2)
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=4)
    if badge:
        ax.text(1.0, 1.30, badge, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=6.5, color="white", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.32,rounding_size=0.5",
                          facecolor=ORANGE, edgecolor="none"), zorder=5)
    ax.tick_params(labelsize=7.5, colors=INK_2, length=2)
    ax.grid(axis="y", color=MUTED, alpha=0.35, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)


@app.command("plot")
def plot(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    lag_quarters: int = typer.Option(LAG_QUARTERS, "--lag-quarters"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff",
                               help="Drop quarters ending after this date"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    top: int = typer.Option(12, "--top", help="Facets in the rising-candidate grid"),
    min_recent: int = typer.Option(3, "--min-recent"),
    start: str = typer.Option("2015-01-01", "--start", help="Plot x-axis start"),
    heatmap_start: str = typer.Option(
        "2012-01-01", "--heatmap-start",
        help="Heatmap x-axis start; defaults to the beginning of the extract"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold"),
    weighted: bool = typer.Option(True, "--weighted/--no-weighted",
                                  help="Position-weighted case definition "
                                       "(GPS_V2_DESIGN.md #1) instead of raw "
                                       "mention counts"),
    role_discount: bool = typer.Option(True, "--role-discount/--no-role-discount",
                                       help="With --weighted, add each "
                                            "substance's own phi_s "
                                            "(_role_dispersion)"),
    spatial: str = typer.Option("discount", "--spatial",
                                help="Spatial-concentration fusion "
                                     "(GPS_V2_DESIGN.md #3): none, discount, "
                                     "or prior"),
    treescan_veto: bool = typer.Option(
        True, "--treescan-veto/--no-treescan-veto",
        help="Suppress credible_rise (today's ranking) where a live "
             "TreeScan scan sees no corroborating rise"),
    veto_threshold: float = typer.Option(
        TREESCAN_VETO_THRESHOLD, "--veto-threshold",
        help="TreeScan recurrence-interval floor below which a rise is vetoed"),
) -> None:
    """Render the figures.

    `res` (today's ranking, used for labels in figures 1/2/4) uses the same
    defaults as `rank` -- current recommendation, live TreeScan veto. The
    alarm bands and figure 5's episodes are read from `eb05_sweep.csv` /
    `alarm_history.csv` on disk (`credible_rise`, already dual-gated and
    TreeScan-vetoed by whatever flags `alarms` was last run with) rather than
    recomputed here -- see `alarms`.
    """
    counts, denom = load_quarterly(mentions, cutoff=cutoff, weighted=weighted,
                                   role_discount=role_discount)
    quarters = sorted(denom.index)
    z_now = None
    if spatial != "none":
        z_now = _z_at(_spatial_sweep(quarters, mentions,
                                     recent_quarters=recent_quarters),
                      quarters[-1 - lag_quarters] if lag_quarters
                      else quarters[-1])
    res, meta = rank_substances(counts, denom, lag_quarters,
                                recent_quarters, baseline_quarters,
                                z_spatial=z_now, spatial_mode=spatial)
    if treescan_veto:
        as_of = quarters[-1 - lag_quarters] if lag_quarters else quarters[-1]
        ri = _treescan_ri_live(mentions, as_of, cutoff)
        res = _apply_treescan_veto(res, ri, veto_threshold)
    res = _attach_novelty(res, mentions)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_q = sorted(denom.index)
    start_ts = pd.Timestamp(start)
    q_plot = [q for q in all_q if q >= start_ts]
    trunc_from = all_q[-lag_quarters] if lag_quarters else None
    window = _fmt_window(meta)

    # Alarm history, if `alarms` has been run. Read from disk rather than
    # recomputed: the sweep refits the prior 45 times and would dominate the
    # runtime of a command whose job is to draw pictures.
    sweep_path = RESULTS_DIR / "eb05_sweep.csv"
    hist, sweep = None, None
    if ALARM_PATH.exists() and sweep_path.exists():
        hist = pd.read_csv(ALARM_PATH, index_col=0,
                           parse_dates=["first_seen", "first_breach",
                                        "last_breach", "peak_as_of"])
        sweep = pd.read_csv(sweep_path, parse_dates=["as_of"])
    else:
        typer.echo("note: no alarm history on disk — run `alarms` first to get "
                   "breach annotations and the alarm-timeline figure")

    def spans(name: str) -> list:
        if sweep is None:
            return []
        return [(s, e) for s, e in _episode_spans(
                    sweep, name, threshold, gate_col="credible_rise")
                if e >= q_plot[0]]

    # --- Figure 1: rising candidates, small multiples -------------------
    cand = res[res["n_recent"] >= min_recent].head(top)
    ncol = 4
    nrow = int(np.ceil(len(cand) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.05 * nrow),
                             sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, (name, row) in zip(axes, cand.iterrows()):
        s = _quarterly_series(counts, name, q_plot)
        is_emergent = bool(hist is not None and name in hist.index
                           and hist.loc[name, "emergent"])
        t = f"{_bold_title(name)}  (n={int(row.n_recent)})"
        t += f"\nEB05+ {row.eb05:.1f}"
        if hist is not None and name in hist.index:
            b = hist.loc[name, "first_breach"]
            t += f"  · fired {b.year}Q{b.quarter}"
        _panel(ax, q_plot, s, t, BLUE, row["first_seen"], trunc_from,
               alarm_spans=spans(name),
               badge="EMERGENT" if is_emergent else None)
    for ax in axes[len(cand):]:
        ax.set_visible(False)
    fig.suptitle("Substances rising fastest in LA County overdose deaths\n"
                 f"ranked by EB05+ (gamma-Poisson shrunken growth); {window}",
                 fontsize=11, color=INK, x=0.01, ha="left")
    caption = (f"Quarterly deaths naming the substance. Orange EMERGENT tag "
               f"= averaged under {EMERGENT_THRESHOLD:g}/quarter before its "
               f"most recent rise; untagged panels are a credible rise in an "
               f"already-established substance. Dotted orange = first LACME "
               f"appearance.\nYellow band = quarters EB05+ actually fired; "
               f"its left edge is when it first fired.")
    if trunc_from is not None:
        caption += (" Dashed grey tail = incomplete reporting, excluded from "
                    "the statistic — the drop there is pending toxicology, not "
                    "a real decline.")
    else:
        caption += (f" Series end at {q_plot[-1].year}Q{q_plot[-1].quarter}; "
                    "later quarters are excluded as unsettled toxicology.")
    fig.text(0.01, 0.005, caption, fontsize=7.5, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.045, 1, 0.93))
    fig.savefig(out_dir / "rising_candidates.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Figure 2: the two novelty axes ---------------------------------
    # Chemical novelty (when NFLIS could report it) against local novelty (when
    # LA County first recorded a death). Plotting EBGM here instead was tried
    # and is uninformative: no substance in the current window has a real
    # signal, so the y-axis collapses and every point piles on the 1998-10
    # censoring floor. These two axes have genuine spread and separate the two
    # kinds of "new".
    sc = res[res["date_added"].notna() & res["first_seen"].notna()].copy()
    sc = sc[sc["n_total"] >= 5]
    fig, ax = plt.subplots(figsize=(10.5, 7))
    for censored, color, label in [
        (False, BLUE, "NFLIS date added is a real date"),
        (True, AQUA, "left-censored — in the catalog at its 1998-10 inception"),
    ]:
        g = sc[sc["date_added_censored"].fillna(False) == censored]
        ax.scatter(g["first_seen"], g["date_added"],
                   s=14 + 9 * np.sqrt(g["n_total"]), c=color, alpha=0.70,
                   edgecolors="white", linewidths=0.8, label=label, zorder=3)

    # Ring the substances the detector has ever fired on. Novelty alone says
    # nothing about whether a substance mattered; this overlays the answer, and
    # shows that the fired-on set is scattered across both novelty axes rather
    # than concentrated among the chemically new.
    ever = set(hist.index) if hist is not None else set()
    fired = sc[sc.index.isin(ever)]
    if len(fired):
        ax.scatter(fired["first_seen"], fired["date_added"],
                   s=90 + 9 * np.sqrt(fired["n_total"]), facecolors="none",
                   edgecolors=ORANGE, linewidths=2.0, alpha=0.55, zorder=4,
                   label="has fired EB05+")

    # Zoom past the dead 1998-2010 stretch (no arrivals are recorded there
    # anyway) so the crowded post-2012 cluster gets more horizontal room.
    # Set before any pixel-space measurement below, since transData needs to
    # reflect the final view for the label de-collision to be accurate.
    ax.set_xlim(left=mdates.date2num(pd.Timestamp("2010-01-01")))

    lo = min(sc["date_added"].min(), sc["first_seen"].min())
    hi = max(sc["date_added"].max(), sc["first_seen"].max())
    ax.plot([lo, hi], [lo, hi], color=INK_2, linewidth=1, linestyle="--",
            zorder=1)
    diag_ann = ax.annotate("arrived in LA as soon as the chemistry existed",
                           (hi, hi), textcoords="offset points",
                           xytext=(-8, 8), fontsize=8, color=INK_2, ha="right")

    # Above the diagonal: LA recorded a death before NFLIS catalogued the
    # substance. Shaded and always labelled, because it is the one region a
    # reader is likely to over-interpret — see the region caption.
    ax.fill_between([lo, hi], [lo, hi], [hi, hi], color=MUTED, alpha=0.16,
                    linewidth=0, zorder=0)
    below = sc[sc["first_seen"] < sc["date_added"]]

    # LACME records start in 2012, so every substance already circulating then
    # shares a first_seen of 2012Q1 — the x-axis is left-censored too, and that
    # floor column is an artifact of the extract window, not simultaneous
    # arrival.
    floor = sc["first_seen"].min()
    ax.axvline(floor, color=MUTED, linewidth=1, zorder=1)
    # The 2010 x-axis start leaves too little room left of the floor line for
    # a caption, so it grows rightward instead, into the shaded above-diagonal
    # band -- the one part of the plot with almost no circles in it.
    floor_ann = ax.annotate(
        "LACME record start\narrival dates left of this line\n"
        "are censored, not simultaneous",
        xy=(mdates.date2num(floor), mdates.date2num(pd.Timestamp("2018-01-01"))),
        textcoords="offset points", xytext=(8, 0),
        fontsize=8, color=INK_2, ha="left", va="center")

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=7,
                          markerfacecolor=c, markeredgecolor="white", label=lab)
               for c, lab in [(BLUE, "First LA Death"),
                              (AQUA, "NFLIS Add Date")]]
    if len(fired):
        handles.append(plt.Line2D(
            [], [], marker="o", linestyle="", markersize=10,
            markerfacecolor="none", markeredgecolor=ORANGE, markeredgewidth=2,
            alpha=0.55, label="EB05++TS Alarmed"))
    # (2024, 2012) in data space, converted to axes fraction rather than
    # anchored via transData directly -- a data-space bbox_to_anchor fights
    # tight_layout's own bbox computation and blows up the saved figure.
    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()
    lx = mdates.date2num(pd.Timestamp("2024-01-01"))
    ly = mdates.date2num(pd.Timestamp("2012-01-01"))
    legend = ax.legend(handles=handles, frameon=False, fontsize=8.5,
                       loc="center right",
                       bbox_to_anchor=((lx - xlo) / (xhi - xlo),
                                       (ly - ylo) / (yhi - ylo)),
                       bbox_transform=ax.transAxes, markerfirst=False)

    # Label only substances that arrived after the censoring floor, plus the
    # two dominant killers (for scale) even though neither is novel by either
    # axis. Breached substances are labelled only when their arrival is
    # *observed*. Ten of the labelled sit on the 2012Q1 censoring floor, where
    # the y-coordinate is an artifact of the extract window.
    labelled = sc[(sc["first_seen"] >= pd.Timestamp("2016-01-01"))
                  | sc.index.isin(below.index)
                  | (sc.index.isin(ever) & (sc["first_seen"] > floor))
                  | sc.index.isin(["Fentanyl", "Methamphetamine"])]

    # Many labelled substances share an x (the floor column) or a y (the
    # 1998-10 chemical-catalog floor), so de-collision has to work in both
    # directions at once rather than just walking down a single column: try
    # each label at its natural position, then alternating steps above and
    # below, and keep the first spot that neither overlaps an already-placed
    # label nor falls outside the axes -- so a crowded row pushes labels up
    # into open space instead of stacking off the bottom of the plot.
    fig.canvas.draw()  # transforms must be current before we measure
    renderer = fig.canvas.get_renderer()
    ax_bbox = ax.get_window_extent(renderer)
    px_per_pt = fig.dpi / 72.0

    pts = []
    for name, r in labelled.iterrows():
        # transData works in numbers, not Timestamps, on a date axis.
        xy = (mdates.date2num(r["first_seen"]), mdates.date2num(r["date_added"]))
        px, py = ax.transData.transform(xy)
        # Clear the marker's own edge, not its center — a fixed offset is
        # fine for the small dots but gets swallowed by the large
        # lifetime-death bubbles (and their orange rings, drawn even bigger)
        # for a name like Carfentanil or Fentanyl.
        s = 14 + 9 * np.sqrt(r["n_total"])
        if name in ever:
            s = max(s, 90 + 9 * np.sqrt(r["n_total"]))
        clear_pt = np.sqrt(s / np.pi) + 4.0
        probe = ax.annotate(name, xy, xytext=(clear_pt, 0),
                            textcoords="offset points", fontsize=8,
                            color=INK, va="center")
        fig.canvas.draw()
        bb = probe.get_window_extent(renderer)
        probe.remove()
        pts.append((name, xy, px, py, clear_pt, bb.width + 4, bb.height + 2))

    # The aqua (left-censored) substances all share the identical 1998-10
    # floor date_added, so their markers form one dense row directly above
    # a wide-open band. Stack those labels straight up instead of running
    # them through the diagonal search below: constant x, varying only the
    # height needed to clear whichever neighbor is already there.
    censored = set(sc.index[sc["date_added_censored"].fillna(False)])
    up_pts = sorted((p for p in pts if p[0] in censored), key=lambda t: t[2])
    pts = [p for p in pts if p[0] not in censored]
    # Top of the plot first, then left to right for anything sharing a row.
    pts.sort(key=lambda t: (-t[3], t[2]))

    # Labelled markers are fixed obstacles too, not just other labels --
    # otherwise a label happily lands on top of a *different* labelled point
    # a short distance away (Furanyl fentanyl's label sitting on U-47700's
    # dot, e.g.) instead of using open space a few points further out. Limited
    # to labelled markers, not all ~200 plotted points: avoiding every small
    # unlabelled dot in the dense floor cluster forces detours through empty
    # space that reads worse than the odd close pass would have.
    marker_boxes = {}
    for name2, r2 in labelled.iterrows():
        xy2 = (mdates.date2num(r2["first_seen"]), mdates.date2num(r2["date_added"]))
        px2, py2 = ax.transData.transform(xy2)
        s2 = 14 + 9 * np.sqrt(r2["n_total"])
        if name2 in ever:
            s2 = max(s2, 90 + 9 * np.sqrt(r2["n_total"]))
        rad2 = np.sqrt(s2 / np.pi) * px_per_pt
        marker_boxes[name2] = Bbox.from_bounds(px2 - rad2, py2 - rad2,
                                               2 * rad2, 2 * rad2)

    def box_for(px, py, clear_pt, w, h, dy, side):
        cy = py + dy * px_per_pt
        clear_px = clear_pt * px_per_pt
        x0 = px + clear_px if side == "right" else px - clear_px - w
        return Bbox.from_bounds(x0, cy - h / 2, w, h)

    def line_clear(px, py, box, side, obstacles, own_box, start_clear_px):
        # The connecting line from marker to text must not cut through a
        # label placed earlier -- two text boxes can be individually clear
        # of each other while the line reaching a *third*, taller label still
        # pierces straight through one of them (a short label like "MDA"
        # sitting right in the x-span of "Ephedrine"'s text, e.g.).
        anchor_x = box.x0 if side == "right" else box.x1
        anchor_y = box.y0 + box.height / 2
        length = np.hypot(anchor_x - px, anchor_y - py)
        if length <= start_clear_px:
            return True
        # Start past this marker's own clearance radius, not at its center --
        # markers can sit concentrically (Fentanyl and Methamphetamine share
        # one point), and every line leaving a shared center necessarily
        # starts inside a sibling's own circle. That's fine; it's only a real
        # block once the line is out in open space and still hits something.
        t0 = start_clear_px / length
        for t in np.linspace(t0, 1.0, 6):
            pt = (px + t * (anchor_x - px), py + t * (anchor_y - py))
            for b in obstacles:
                if b is own_box:
                    continue
                if b.x0 <= pt[0] <= b.x1 and b.y0 <= pt[1] <= b.y1:
                    return False
        return True

    step_pt = 9.0
    placed_boxes = [diag_ann.get_window_extent(renderer),
                    floor_ann.get_window_extent(renderer),
                    legend.get_window_extent(renderer)] + list(marker_boxes.values())

    # Straight up, 3/4" (54pt) above the marker, centered on it. Each label
    # only climbs past that base height if a narrower gap would collide with
    # a neighbor already placed -- most stay at the same height, the few in
    # a tight patch of the row step up just enough to clear it.
    base_up_pt, step_up_pt = 54.0, 12.0
    for name, xy, px, py, clear_pt, w, h in up_pts:
        offset_pt = base_up_pt
        while True:
            offset_px = offset_pt * px_per_pt
            box = Bbox.from_bounds(px - w / 2, py + offset_px, w, h)
            if box.y1 <= ax_bbox.y1 and not any(box.overlaps(b) for b in placed_boxes):
                break
            offset_pt += step_up_pt
        placed_boxes.append(box)
        ax.annotate(
            name, xy, xytext=(0, offset_pt), textcoords="offset points",
            fontsize=8, color=INK, ha="center", va="bottom", zorder=6,
            bbox=dict(boxstyle="square,pad=0.1", facecolor="white",
                     edgecolor="none", alpha=0.7),
            arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=0.6,
                            shrinkA=0, shrinkB=1),
        )

    for name, xy, px, py, clear_pt, w, h in pts:
        own_box = marker_boxes.get(name)
        best = None  # (side, dy, box) -- smallest |dy| wins, right breaks ties
        fallback = None
        for side in ("right", "left"):
            for k in range(16):
                for dy in ([0.0] if k == 0 else [step_pt * k, -step_pt * k]):
                    box = box_for(px, py, clear_pt, w, h, dy, side)
                    if box.y0 < ax_bbox.y0 or box.y1 > ax_bbox.y1:
                        continue
                    if fallback is None:
                        fallback = (side, dy, box)
                    if any(box.overlaps(b) for b in placed_boxes):
                        continue
                    if not line_clear(px, py, box, side, placed_boxes, own_box,
                                      clear_pt * px_per_pt):
                        continue
                    if best is None or abs(dy) < abs(best[1]):
                        best = (side, dy, box)
                    break
                if best is not None and abs(best[1]) == (0.0 if k == 0 else step_pt * k):
                    break
            if best is not None and best[1] == 0.0:
                break  # can't beat a clean, undisplaced placement
        side, chosen_dy, chosen_box = best if best is not None else fallback
        placed_boxes.append(chosen_box)
        dx0 = clear_pt if side == "right" else -clear_pt
        ax.annotate(
            name, xy, xytext=(dx0, chosen_dy), textcoords="offset points",
            fontsize=8, color=INK, va="center",
            ha=("left" if side == "right" else "right"), zorder=6,
            bbox=dict(boxstyle="square,pad=0.1", facecolor="white",
                     edgecolor="none", alpha=0.7),
            arrowprops=(dict(arrowstyle="-", color=MUTED, linewidth=0.6,
                             shrinkA=0, shrinkB=1) if abs(chosen_dy) > 4 else None),
        )

    ax.set_xlabel("First LA County overdose death  —  local novelty",
                  fontsize=9.5, color=INK_2)
    ax.set_ylabel("Date added to the NFLIS catalog  —  chemical novelty",
                  fontsize=9.5, color=INK_2)
    # Region caption lives in the figure margin, not inside the axes. Anchored
    # in the shaded region it sat on top of the very points it describes, and
    # every attempt to reposition it collided with a different label — the plot
    # interior has no free space at this label density.
    fig.text(
        0.005, -0.012,
        "Shaded region, above the diagonal: the LA death predates the NFLIS "
        "entry. Not early detection — NFLIS catalogues substances as forensic "
        "labs report them in seized drug\nevidence, and these are prescription "
        "antihistamines that are found at autopsy but rarely seized. "
        "Orange rings mark substances that have gone EB05++TS Alarmed; "
        "the quarter is the first time.",
        fontsize=8, color=INK_2, ha="left", va="top")

    ax.set_title(
        "Two kinds of new: arrival in LA County vs chemical novelty\n"
        "points far below the diagonal are established chemistry reaching LA "
        "late — diffusion, not invention. Marker size = lifetime deaths.",
        fontsize=11, color=INK, loc="left")
    ax.grid(color=MUTED, alpha=0.35, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "novelty_vs_growth.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Figure 3: share-of-deaths heatmap ------------------------------
    # Runs the full extract by default rather than the 2015 start the line
    # panels use: a share matrix is legible at 56 columns where a 12-facet line
    # grid is not, and the early years are where the panel-drift caveat below
    # actually shows up.
    q_heat = [q for q in all_q if q >= pd.Timestamp(heatmap_start)]
    topn = res[res["n_recent"] >= min_recent].nlargest(25, "share_recent_pct").index
    wide = (counts[counts["substance"].isin(topn)]
            .pivot_table(index="substance", columns="quarter", values="n",
                         aggfunc="sum", fill_value=0)
            .reindex(index=topn, columns=q_heat, fill_value=0))
    share = 100 * wide.div(denom.reindex(q_heat).to_numpy(), axis=1)

    # Detection-practice drift. Substances named per death climbs from 1.46 in
    # 2012 to ~1.9 and plateaus; before that plateau, a low share can mean the
    # ME was not looking rather than the substance was not there. Mark where
    # the climb ends so the extended left edge is not read as absence. The
    # plateau is the mean over the last 12 quarters and the marker is the first
    # quarter reaching 95% of it.
    spd = (counts.groupby("quarter")["n"].sum()
           .reindex(all_q, fill_value=0) / denom.reindex(all_q))
    plateau = float(spd.iloc[-12:].mean())
    reached = spd.rolling(4, min_periods=4).mean() >= 0.95 * plateau
    drift_end = next((q for q in all_q if bool(reached.get(q, False))), None)
    fig, ax = plt.subplots(figsize=(12, 8))
    # Methamphetamine and fentanyl reach ~70% of deaths; on a linear ramp they
    # saturate the scale and every candidate under ~5% — the entire population
    # of interest — renders as blank white. A power norm keeps the axis in real
    # percentage units while making the low end legible.
    im = ax.imshow(share.to_numpy(), aspect="auto", cmap=SEQ,
                   interpolation="nearest",
                   norm=PowerNorm(gamma=0.4, vmin=0, vmax=float(share.to_numpy().max())))
    labels = [f"{s}  ({int(res.loc[s, 'n_total'])})" for s in share.index]
    ax.set_yticks(range(len(share)), labels, fontsize=8.5, color=INK)
    tick = range(0, len(q_heat), 4)
    ax.set_xticks(list(tick),
                  [f"{q_heat[i].year}Q{q_heat[i].quarter}" for i in tick],
                  fontsize=8, color=INK_2)
    if trunc_from is not None and trunc_from in q_heat:
        ax.axvline(q_heat.index(trunc_from) - 0.5, color=ORANGE, linewidth=1.5)
    drift_note = ""
    if drift_end is not None and drift_end in q_heat:
        xd = q_heat.index(drift_end) - 0.5
        ax.axvline(xd, color=YELLOW, linewidth=1.8, zorder=5)
        drift_note = (
            f"Yellow line ({drift_end.year}Q{drift_end.quarter}): substances "
            f"named per death stops climbing ({spd.iloc[0]:.2f} in "
            f"{all_q[0].year} → {plateau:.2f} at plateau). Left of it the "
            "toxicology panel was still widening, so a faint cell can mean "
            "the ME was not testing rather than the substance was not there.")
    cb = fig.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("% of all overdose deaths that quarter", fontsize=9, color=INK_2)
    cb.ax.tick_params(labelsize=8, colors=INK_2)
    sub = ("orange line = start of incomplete reporting"
           if trunc_from is not None and trunc_from in q_heat
           else f"{q_heat[0].year}Q{q_heat[0].quarter}–"
                f"{q_heat[-1].year}Q{q_heat[-1].quarter}, ending at the last "
                "quarter with settled toxicology")
    ax.set_title("Substance share of LA County overdose deaths, by quarter\n"
                 "top 25 by current share (lifetime deaths in parentheses); "
                 + sub, fontsize=11, color=INK, loc="left")
    if drift_note:
        fig.text(0.005, -0.005, drift_note, fontsize=8, color=INK_2,
                 ha="left", va="top")
    fig.tight_layout()
    fig.savefig(out_dir / "share_heatmap.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    # --- Figure 4: known emergences (backtest ground truth) -------------
    present = [s for s in KNOWN_EMERGENCES if s in set(counts["substance"])]
    fig, axes = plt.subplots(1, len(present), figsize=(2.9 * len(present), 2.6),
                             sharex=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, present):
        s = _quarterly_series(counts, name, q_plot)
        fs = res.loc[name, "first_seen"] if name in res.index else None
        sub = f"{_bold_title(name)}  n={int(s.sum())}"
        if hist is not None and name in hist.index:
            h = hist.loc[name]
            sub += f"\nfired {h.first_breach.year}Q{h.first_breach.quarter}"
            if pd.notna(h.lag_quarters):
                sub += f" (+{int(h.lag_quarters)}q)"
            sub += f", peak {h.peak_eb05:.1f}"
        elif sweep is not None:
            g = sweep[sweep["substance"] == name]
            if len(g):
                sub += f"\nnever fired, peak {g['eb05'].max():.1f}"
        _panel(ax, q_plot, s, sub, BLUE, fs, trunc_from,
               alarm_spans=spans(name))
    fig.suptitle("Known emerging substances in LA County\n"
                 "quarterly deaths naming each substance, with EB05+'s alarm history",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.text(0.01, 0.005,
             "Yellow band = quarters EB05+ actually fired; its left edge is "
             "when it first fired. Dotted orange = first LA death. (+Nq) = "
             "quarters between the two.",
             fontsize=7.5, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.86))
    fig.savefig(out_dir / "known_emergences.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)

    written = ["rising_candidates.png", "novelty_vs_growth.png",
               "share_heatmap.png", "known_emergences.png"]

    # --- Figure 5: alarm history ----------------------------------------
    # The three dates a reviewer asked for, per substance, on one time axis:
    # first LA death, first EB05 breach, end of the alarm. The gap between the
    # first two is the detection delay and is the point of the figure.
    if hist is not None and len(hist):
        h = hist.sort_values("first_breach", ascending=False)
        fig, ax = plt.subplots(figsize=(13, 0.36 * len(h) + 2.6))
        y = np.arange(len(h))
        t_end = max(all_q) + pd.offsets.QuarterEnd(0)

        x_left = pd.Timestamp("2011-07-01")
        any_open = bool(h["still_open"].any())
        for i, (name, r) in enumerate(h.iterrows()):
            # Presence line: from the first LA death to the end of the window.
            # A censored onset gets a left caret at the axis edge instead of a
            # circle, because the line does not start there — the records do.
            if pd.notna(r["first_seen"]):
                cens = bool(r.get("first_seen_censored", False))
                ax.hlines(i, x_left if cens else r["first_seen"], t_end,
                          color=MUTED, alpha=0.55, linewidth=1.4, zorder=2)
                ax.scatter([x_left if cens else r["first_seen"]], [i],
                           marker="<" if cens else "o", s=30 if cens else 26,
                           c=MUTED if cens else "white", edgecolors=INK_2,
                           linewidths=1.1, zorder=4)
            for s, e in _episode_spans(sweep, name, threshold,
                                       gate_col="credible_rise"):
                ax.barh(i, e + pd.offsets.QuarterEnd(0) - s, left=s,
                        height=0.52, color=ORANGE if r["still_open"] else BLUE,
                        edgecolor="white", linewidth=0.8, zorder=3)
            ax.scatter([r["peak_as_of"]], [i], marker="D", s=30, c=YELLOW,
                       edgecolors="white", linewidths=0.9, zorder=5)
            ax.annotate(f"peak {r['peak_eb05']:.1f}  (n={int(r['n_at_peak'])})",
                        (t_end, i), textcoords="offset points", xytext=(8, 0),
                        fontsize=7.5, color=INK_2, va="center",
                        annotation_clip=False)
            # Detection delay. Placed above the bar rather than to its left
            # when the gap is short, or the text lands on the onset marker.
            if pd.notna(r["lag_quarters"]):
                short = (r["first_breach"] - r["first_seen"]).days < 550
                ax.annotate(
                    f"{int(r['lag_quarters'])}q", (r["first_breach"], i),
                    textcoords="offset points",
                    xytext=((2, 9) if short else (-5, 0)),
                    fontsize=7, color=INK_2,
                    va=("bottom" if short else "center"),
                    ha=("left" if short else "right"))

        ax.set_yticks(y, h.index, fontsize=8.5, color=INK)
        ax.set_ylim(-0.8, len(h) - 0.2)
        ax.set_xlim(x_left - pd.Timedelta(days=60),
                    t_end + pd.Timedelta(days=200))
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.grid(axis="x", color=MUTED, alpha=0.30, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8, colors=INK_2, length=2)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(MUTED)

        leg = [
            plt.Line2D([], [], marker="o", linestyle="-", color=MUTED,
                       markerfacecolor="white", markeredgecolor=INK_2,
                       markersize=6,
                       label="first LA death naming it, then present"),
            plt.Line2D([], [], marker="<", linestyle="", color=MUTED,
                       markeredgecolor=INK_2, markersize=7,
                       label="already present when LACME records start (2012)"),
            plt.Rectangle((0, 0), 1, 1, facecolor=BLUE,
                          label="in alarm (EB05+)"),
            plt.Line2D([], [], marker="D", linestyle="", color=YELLOW,
                       markeredgecolor="white", markersize=7,
                       label="peak EB05"),
        ]
        if any_open:
            leg.insert(3, plt.Rectangle((0, 0), 1, 1, facecolor=ORANGE,
                                        label="alarm still open at the cutoff"))
        ax.legend(handles=leg, frameon=False, fontsize=8.5, ncol=len(leg),
                  loc="upper left", bbox_to_anchor=(0, -0.045))
        n_swept = len(all_q) - RECENT_QUARTERS - BASELINE_QUARTERS + 1
        # Denominator is the substances the sweep actually scored (non-zero in
        # some recent window), matching what `alarms` reports — not len(res),
        # which counts rows the sweep never scored.
        n_scored = sweep["substance"].nunique()
        ax.set_title(
            "Every substance the detector has ever fired on, and when\n"
            f"{len(h)} of {n_scored} substances have fired EB05+ "
            f"across {n_swept} as-of quarters. The number at each "
            "bar is the detection delay — quarters from the first LA death to "
            "the first breach — and is shown only where arrival is observed, "
            "not censored.",
            fontsize=11, color=INK, loc="left")
        fig.tight_layout()
        fig.savefig(out_dir / "alarm_history.png", dpi=150,
                    bbox_inches="tight", facecolor="white")
        plt.close(fig)
        written.append("alarm_history.png")

    for f in written:
        typer.echo(f"wrote {out_dir / f}")


if __name__ == "__main__":
    app()
