"""
One table, one row per detector, scored against the same backtest: how
reliably does each variant of the ranking statistic identify the known LA
emergences, and what does it cost to run it that hot.

**Why this exists.** `trends backtest` answers "did *this* configuration catch
the five known emergences, and when", and `ground-truth score` answers "did
its alarm episodes line up with the reference intervals". Neither compares
configurations. `docs/GPS_V2_DESIGN.md` accumulated three proposed extensions
to the base GPS/EB05 statistic, each validated only against plain EB05 one at
a time, and each with its own separate verdict written up in its own section.
That is not enough to say which detector to actually run. This scores all of
them side by side on identical data, identical as-of discipline, and two
identical reading rules.

**The two reading rules, and why one alone would be misleading.**

1. *Fixed line* — every model read at `ALARM_THRESHOLD` (1.5), the line the
   backtest calibrated for plain EB05. This is what a reader would actually
   do, but it flatters any model that inflates its statistic: the spatial
   discount raises `E`, which raises EB05 for everything scored, so it buys
   recall by firing more often rather than by discriminating better.
2. *Matched alarm budget* — each model's threshold is re-cut so it raises the
   **same total number of substance-quarter alarms** as plain EB05 at 1.5.
   Now a model can only improve recall by spending its fixed alarm budget on
   the right substance-quarters. This is the comparison that means something,
   and the precedent is `RESEARCH_LOG.md`'s tree-scan head-to-head, which
   found the scan's apparent earlier detection disappeared once its budget was
   matched to EB05's.

A model that improves on the fixed line and not on the matched budget has not
improved the detector; it has lowered the threshold.

**What the precision column is and is not.** `data/raw/ground_truth/
known_emergences.csv` labels five substances. It does not assert that every
*other* alarm is false — lidocaine and PCP both alarm for documented reasons
and neither is in the file. So `precision` here is computed the way
`ground_truth.score` computes it, within the labelled substances only, and the
county-wide cost of running hot is reported separately as `off_target_*`:
alarms on unlabelled substances, which are a *load* on a reviewer, not a
proven error rate. Reading `off_target_quarters` as a false-positive count
would be overclaiming.

Usage:
    emerging benchmark run
    emerging benchmark run --threshold 1.0
    emerging benchmark run --models eb05,eb05-v2
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import typer

from emerging.analysis.aberration import EARS_THRESHOLD, NB_TREND_THRESHOLD
from emerging.analysis.aberration import ears_sweep as aberration_ears_sweep
from emerging.analysis.aberration import \
    nb_trend_sweep as aberration_nb_trend_sweep
from emerging.analysis.trends import (ALARM_THRESHOLD, SPATIAL_WEIGHT,
                                      _episodes, _spatial_sweep,
                                      load_quarterly, sweep_eb05)
from emerging.analysis.treescan import RI_THRESHOLD as TREESCAN_RI_THRESHOLD
from emerging.analysis.treescan import branch_sweep as treescan_branch_sweep
from emerging.analysis.treescan import leaf_sweep as treescan_leaf_sweep
from emerging.analysis.treescan import prepare as treescan_prepare
from emerging.config import BASELINE_QUARTERS, CUTOFF, RECENT_QUARTERS
from emerging.ingest.extract import MENTIONS_PATH
from emerging.paths import KNOWN_EMERGENCES_PATH, results_dir
from emerging.validation.ground_truth import (interval_overlap,
                                              load_ground_truth, precision_rate)
from emerging.viz import BLUE, INK, MUTED, ORANGE

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("benchmark")

# The three readings of one tree-temporal scan: all read off
# `recurrence_interval` and all cut at TreeScan's own RI line, but each is a
# different table and so a different sweep key. Slot 0 of that key is the
# statistic name, which is what keeps it from colliding with the EB05
# family's `(weighted, role_discount, spatial)`.
TREESCAN_STATISTICS = frozenset(
    {"treescan", "treescan_branch", "treescan_tree"})
# The "model the curve, not the ratio" family (docs/LITERATURE_REVIEW.md sec
# 1.5): a moving-baseline z-score (EARS C2) and a log-linear trend slope's
# Wald z (the direct disproportionality-free analogue of "large derivative").
# Both take the same `(weighted, role_discount)` axis GPS_V2_DESIGN.md #1
# gives the EB05 family (`Model.weighted`/`.role_discount`, reused rather than
# duplicated) -- `aberration.ears_sweep`/`nb_trend_sweep` accept the identical
# flags and recompute `load_quarterly` per as-of the same way `sweep_eb05`
# does. `nb_trend`/`nb_trend_own`/`nb_trend_dual`/`nb_trend_emergent`/
# `nb_trend_confirmed` are five readings of one sweep (`aberration.
# nb_trend_sweep` fits both trends and the `emergent`/`confirmed` flags
# together, same reason `eb05`/`eb05_dual` share `sweep_eb05`), so they
# collapse to the same key -- see `Model.sweep_key`.
NB_TREND_SWEEP_STATISTICS = frozenset(
    {"nb_trend", "nb_trend_own", "nb_trend_dual", "nb_trend_emergent",
     "nb_trend_confirmed", "nb_trend_confirmed_emergent"})
# The attribution line the branch models are guarded at, and the band
# `treescan.backtest` already calls "partly attributable" -- below it, a
# branch is signalling on something other than the substance being credited.
MIN_EXCESS_SHARE = 0.15


@dataclass(frozen=True)
class Model:
    """One detector: a case definition, a gate, and a spatial fusion.

    `statistic` picks what the alarm is read off:

    - `"ratio"` — the unshrunken `n_recent / expected`, the thing GPS exists
      to replace. Included as the floor, not as a candidate: at n = 6 one
      extra death doubles it, which is the whole argument for shrinkage and is
      worth showing rather than asserting.
    - `"eb05"` — the posterior 5th percentile, single gate (Gate A, share of
      all OD deaths).
    - `"eb05_dual"` — `min(eb05, eb05_own)`, which is exactly the
      `credible_rise` conjunction of `docs/GPS_V2_DESIGN.md` #2 written as one
      monotone number. Writing it as a min rather than as two comparisons is
      what lets the dual gate be threshold-swept for the matched-budget
      reading at all; at the fixed line the two forms are identical by
      construction.
    - `"eb05_own"` — Gate B alone (the self-referential, own-count-only
      posterior 5th percentile), with no share gate at all. Not a fallback or
      an approximation of `eb05_dual` -- it answers a different question ("is
      this substance's own count elevated against its own history") with no
      reference to the county total, so it never catches the denominator-
      drift failure Gate A exists for (see `docs/findings/benchmark.md`'s
      cocaine example). Included to check whether Gate A is pulling its
      weight at all once `--weighted`/`--role-discount` is already in the
      case definition, not because it's expected to win.
    - `"treescan"` — the substance-leaf node's recurrence interval from
      `emerging.analysis.treescan`, P1's tree-temporal scan run solo (no
      EB05 underneath it at all). Its native scale (RI, roughly 1..10000) has
      nothing in common with EB05's (a posterior percentile, roughly 0..3), so
      the shared `--threshold` at the fixed reading would be meaningless for
      it -- see `fixed_threshold`.
    - `"treescan_branch"` — the same scan read at its *internal* nodes: a
      substance scores whatever the best branch containing it scored. This is
      the mechanism a tree scan exists for and the one the leaf reading
      cannot exercise at all -- the nitazenes are 1 and 2 mentions apiece and
      can never clear a leaf test, but the family node holding them is a
      single hypothesis with their pooled count. `min_excess_share` decides
      whether the branch's score is credited to every member present or only
      to the members actually driving it.
    - `"treescan_tree"` — `max(leaf, branch)`, the whole hierarchy read as
      one detector, which is what TreeScan is. Included because neither half
      alone is the method: the leaf reading is the method with its
      aggregation removed, and the branch reading is the method with its
      most-specific hypotheses removed.
    - `"ears"` — EARS C2-style moving-baseline z-score
      (`emerging.analysis.aberration.ears_sweep`), on the same share-of-OD-
      deaths metric as EB05 Gate A. Its native scale is a z-score, not a
      posterior percentile, so it needs `fixed_threshold` the same way
      `"treescan"` does. `weighted`/`role_discount` apply
      `GPS_V2_DESIGN.md` #1 to the case definition underneath it, exactly as
      they do for `eb05` -- see `aberration.ears_sweep`.
    - `"nb_trend"` — log-linear Poisson trend slope, Wald z-scored
      (`emerging.analysis.aberration.nb_trend_sweep`) -- "model the curve,
      alarm on the derivative" read literally, rather than a shrunk ratio of
      two periods. Also its own scale, also needs `fixed_threshold`.
      `weighted`/`role_discount` apply the same way as for `"ears"`.
    - `"nb_trend_own"` — the same slope fit with no county-total offset, so
      it reads a trend in the substance's *raw count* rather than its share.
      The `nb_trend` analogue of Gate B (`eb05_own`).
    - `"nb_trend_dual"` — `min(nb_trend, nb_trend_own)`, the `nb_trend`
      analogue of `eb05_dual` -- see `nb_trend_sweep`'s docstring for the
      Methamphetamine/Cocaine denominator-drift case this exists to catch.
    - `"nb_trend_emergent"` — `nb_trend_z` gated by `aberration._emergent_flags`
      (`sweep["emergent"]`): floored to `-inf` wherever the substance's
      current alarm episode was not previously rare, exactly the way
      `credible_rise` is a conjunction rather than a blend. The `nb_trend`
      analogue of `trends.alarm_history`'s `emergent` annotation, but wired
      as a hard gate here rather than a label, since `EMERGENT_THRESHOLD` is
      never used to *suppress* an EB05 alarm and this benchmark wants to ask
      whether it should for `nb_trend` specifically (Methamphetamine).
    - `"nb_trend_confirmed"` — `nb_trend_z` gated by
      `aberration._confirmed_flags` (`sweep["confirmed"]`): floored to
      `-inf` unless the current quarter is the `NB_TREND_CONFIRM_QUARTERS`-th
      or later of a run of *consecutive* alarming quarters. Not an EB05
      analogue -- built from `docs/findings/synthetic.md` Finding 1, which
      found `nb_trend`'s single-quarter read has no correction for the
      ~34-quarter as-of sweep it's read across, inflating its real
      false-alarm rate from a nominal ~2.5% to ~33% on synthetic data with
      no real trend at all.
    - `"nb_trend_confirmed_emergent"` — `nb_trend_z` gated by *both*
      `confirmed` and `emergent` at once (a conjunction, not a blend, the
      same `min`-as-AND move `eb05_dual` makes). Checked against real data
      before this was added: `nb-trend-confirmed` and `nb-trend-emergent`
      are non-redundant -- `confirmed` cuts the genuinely noise-driven
      off-target cases (a single quarter clearing the bar by chance) and
      barely touches Methamphetamine (a real, sustained trend that easily
      clears 3 consecutive quarters); `emergent` fully zeroes Methamphetamine
      and does nothing for one-off noise. This asks whether the two combine
      to close both failure modes, or whether stacking a second precision
      mechanism spends more recall than it's worth (`docs/findings/
      benchmark.md`'s recall/precision frontier discussion applies here too).
    - `"ensemble"` — `components` (other model ids) combined per `combine`;
      see `build_ensemble_sweep` for the four modes (`"mean"`, `"max"`,
      `"threshold"`, `"veto"`) and what each one actually does to the
      backtest. `"mean"`/`"max"`/`"threshold"` rank or normalize every
      component onto a common scale before combining, because the candidates
      otherwise share no common unit (EB05's posterior percentile,
      TreeScan's Monte Carlo recurrence interval, the raw ratio's unbounded
      n/E); `"veto"` does not need to, since it never combines two scores
      into one number -- it reads `components[0]`'s own raw score and only
      asks yes/no whether `components[1]` clears `veto_threshold`.

    `fixed_threshold` lets a model be read at its own natively-calibrated line
    for the *fixed* reading instead of the `--threshold` shared by the EB05
    family; the matched-budget reading always re-cuts every model's own score
    column regardless, so it needs no override. A `"mean"`/`"max"` ensemble's
    score is a rank in [0, 1] and has no natural non-budget line either, so it
    sets one too (0.5 = "above the middle of its components' own rankings, on
    average"). `"threshold"`'s score is each component's own ratio to *its*
    calibrated line, so 1.0 is the natural line. `"veto"` stays in its
    proposer's native units, so it needs no override at all -- `None` falls
    back to the shared `--threshold`, exactly like the proposer read alone.
    """

    id: str
    label: str
    statistic: str          # "ratio" | "eb05" | "eb05_dual" | "treescan*" | "ensemble"
    weighted: bool = False  # GPS_V2_DESIGN.md #1
    role_discount: bool = False  # #1's substance-level phi_s, only if weighted
    spatial: str = "none"   # GPS_V2_DESIGN.md #3: none | discount | prior
    fixed_threshold: float | None = None  # override --threshold for this model
    # Only for the treescan branch/tree statistics: the fraction of a node's
    # excess a substance must contribute before that node's score is credited
    # to it. `None` credits every member present in the node's window, which
    # is what lets the guard's cost be measured rather than assumed.
    min_excess_share: float | None = None
    components: tuple[str, ...] = ()  # model ids, only if statistic == "ensemble"
    combine: str = "mean"  # "mean" | "max" | "threshold" | "veto", ensemble only
    veto_threshold: float | None = None  # combine == "veto" only; see below

    @property
    def sweep_key(self) -> tuple:
        """Models differing only in how they *read* one sweep share it."""
        if self.statistic in TREESCAN_STATISTICS:
            return (self.statistic, self.min_excess_share)
        if self.statistic in NB_TREND_SWEEP_STATISTICS:
            return ("nb_trend", self.weighted, self.role_discount)
        if self.statistic == "ears":
            return ("ears", self.weighted, self.role_discount)
        if self.statistic == "ensemble":
            return ("ensemble", self.combine, self.veto_threshold) + self.components
        return (self.weighted, self.role_discount, self.spatial)


MODELS: list[Model] = [
    Model("ratio", "raw n/E, no shrinkage", "ratio"),
    Model("eb05", "classic EB05 (single gate, raw counts)", "eb05"),
    Model("eb05-weighted", "EB05, position-weighted case definition (#1)",
          "eb05", weighted=True),
    Model("eb05-dual", "EB05 dual-gated, share AND own-count (#2)",
          "eb05_dual"),
    Model("eb05-spatial", "EB05, spatial-concentration E-discount (#3)",
          "eb05", spatial="discount"),
    Model("eb05-spatial-prior", "EB05, spatial covariate-adjusted prior (#3)",
          "eb05", spatial="prior"),
    Model("eb05-v2", "GPS v2: weighted + dual-gated + spatial (#1+#2+#3)",
          "eb05_dual", weighted=True, spatial="discount"),
    Model("eb05-weighted-role",
          "EB05, position-weighted + substance role discount (#1, extended)",
          "eb05", weighted=True, role_discount=True),
    Model("eb05-v2-role",
          "GPS v2 + role discount: weighted+role, dual-gated, spatial "
          "(#1 extended+#2+#3) -- candidate to replace eb05-v2",
          "eb05_dual", weighted=True, role_discount=True, spatial="discount"),
    Model("eb05-dual-role",
          "weighted+role, dual-gated, no spatial (#1 extended+#2) -- "
          "eb05-v2-role without #3",
          "eb05_dual", weighted=True, role_discount=True),
    Model("eb05-own-role",
          "weighted+role, Gate B (self-referential) alone -- no share gate, "
          "no dual AND, no spatial",
          "eb05_own", weighted=True, role_discount=True),
    Model("treescan",
          "TreeScan, substance-leaf node, deaths reference (P1, solo)",
          "treescan", fixed_threshold=TREESCAN_RI_THRESHOLD),
    Model("treescan-branch",
          "TreeScan, branch nodes only, credited to every member present -- "
          "the unguarded form, kept to price the attribution guard",
          "treescan_branch", fixed_threshold=TREESCAN_RI_THRESHOLD),
    Model("treescan-branch-attr",
          "TreeScan, branch nodes only, credited only where the substance "
          f"drives >= {MIN_EXCESS_SHARE:.0%} of the branch's excess",
          "treescan_branch", fixed_threshold=TREESCAN_RI_THRESHOLD,
          min_excess_share=MIN_EXCESS_SHARE),
    Model("treescan-tree",
          "TreeScan as designed: max(leaf, attributed branch) over the whole "
          "hierarchy -- the aggregation the leaf-only row cannot show",
          "treescan_tree", fixed_threshold=TREESCAN_RI_THRESHOLD,
          min_excess_share=MIN_EXCESS_SHARE),
    Model("ears", "EARS C2-style moving-baseline z-score on share of OD "
          "deaths -- CDC syndromic-surveillance aberration detection, no "
          "EB05/GPS underneath it at all (docs/LITERATURE_REVIEW.md sec 1.5)",
          "ears", fixed_threshold=EARS_THRESHOLD),
    Model("nb-trend", "log-linear Poisson trend, Wald z of the fitted slope "
          "over the same window EB05 sees -- 'model the curve, alarm on the "
          "derivative' read literally (docs/LITERATURE_REVIEW.md sec 1.5)",
          "nb_trend", fixed_threshold=NB_TREND_THRESHOLD),
    Model("nb-trend-own", "nb-trend's Gate B analogue: trend slope on the "
          "substance's raw count, no county-total offset at all",
          "nb_trend_own", fixed_threshold=NB_TREND_THRESHOLD),
    Model("nb-trend-dual", "nb-trend's eb05-dual analogue: min(nb-trend, "
          "nb-trend-own) -- built to catch the same share-vs-count "
          "denominator drift Gate B catches for EB05 (Methamphetamine, "
          "Cocaine on nb-trend's own off-target list, see "
          "docs/findings/benchmark.md)",
          "nb_trend_dual", fixed_threshold=NB_TREND_THRESHOLD),
    Model("nb-trend-emergent", "nb-trend gated by the EMERGENT_THRESHOLD "
          "analogue (aberration._emergent_flags): alarms only where the "
          "substance's own alarm episode reads as previously rare, "
          "mirroring trends.alarm_history's `emergent` column but wired as "
          "a hard gate -- built to test whether it, not the dual gate, is "
          "what actually resolves the Methamphetamine finding",
          "nb_trend_emergent", fixed_threshold=NB_TREND_THRESHOLD),
    Model("nb-trend-confirmed", "nb-trend gated by requiring "
          "NB_TREND_CONFIRM_QUARTERS consecutive alarming quarters "
          "(aberration._confirmed_flags) -- the repeated-testing "
          "correction docs/findings/synthetic.md Finding 1 found necessary, "
          "not an EMERGENT_THRESHOLD analogue",
          "nb_trend_confirmed", fixed_threshold=NB_TREND_THRESHOLD),
    Model("nb-trend-confirmed-emergent", "nb-trend gated by both confirmed "
          "AND emergent at once -- does stacking the repeated-testing fix "
          "on top of the established-vs-emerging fix close both failure "
          "modes, given they're independently non-redundant on real data?",
          "nb_trend_confirmed_emergent", fixed_threshold=NB_TREND_THRESHOLD),
    Model("ears-weighted-role", "ears, position-weighted + role-discounted "
          "case definition (GPS_V2_DESIGN.md #1, extended) -- does the same "
          "case-definition fix that helps eb05-weighted-role transfer to a "
          "detector with no shrinkage in it at all?",
          "ears", weighted=True, role_discount=True,
          fixed_threshold=EARS_THRESHOLD),
    Model("nb-trend-weighted-role", "nb-trend, position-weighted + "
          "role-discounted case definition (#1, extended)",
          "nb_trend", weighted=True, role_discount=True,
          fixed_threshold=NB_TREND_THRESHOLD),
    Model("nb-trend-dual-weighted-role", "nb-trend-dual, position-weighted "
          "+ role-discounted case definition (#1, extended) -- the "
          "position-weighted analogue of eb05-v2-role's own #1-extended "
          "layer, minus the spatial term (see aberration.py's docstring for "
          "why #3 doesn't port the same way)",
          "nb_trend_dual", weighted=True, role_discount=True,
          fixed_threshold=NB_TREND_THRESHOLD),
    Model("ensemble-mean", "ensemble, soft AND: mean rank of "
          "eb05-dual + treescan -- kept as a negative comparison, see "
          "build_ensemble_sweep's docstring",
          "ensemble", fixed_threshold=0.5, combine="mean",
          components=("eb05-dual", "treescan")),
    Model("ensemble-max", "ensemble, soft OR by rank: max rank of "
          "eb05-dual + treescan -- also loses to eb05-dual alone, see "
          "build_ensemble_sweep's docstring",
          "ensemble", fixed_threshold=0.5, combine="max",
          components=("eb05-dual", "treescan")),
    Model("ensemble-threshold", "ensemble, soft OR by native threshold: "
          "eb05-dual OR treescan, each at its own calibrated line",
          "ensemble", fixed_threshold=1.0, combine="threshold",
          components=("eb05-dual", "treescan")),
    Model("ensemble-threshold-ratio", "ensemble, soft OR by native "
          "threshold: eb05-dual + treescan + raw ratio, each at its own line",
          "ensemble", fixed_threshold=1.0, combine="threshold",
          components=("eb05-dual", "treescan", "ratio")),
    Model("ensemble-veto-2", "ensemble: eb05-dual proposes, treescan vetoes "
          "below RI 2 -- matches or slightly beats eb05-dual alone at both "
          "readings, at zero-to-one fewer off-target substances",
          "ensemble", combine="veto", veto_threshold=2.0,
          components=("eb05-dual", "treescan")),
    Model("ensemble-veto-10", "ensemble: eb05-dual proposes, treescan vetoes "
          "below RI 10 -- best IoU and fully clean off-target of any model "
          "in the table, at one fewer detected emergence",
          "ensemble", combine="veto", veto_threshold=10.0,
          components=("eb05-dual", "treescan")),
    Model("ensemble-veto-v2role-2", "ensemble: eb05-v2-role (the best solo "
          "GPS detector) proposes, treescan vetoes below RI 2",
          "ensemble", combine="veto", veto_threshold=2.0,
          components=("eb05-v2-role", "treescan")),
    Model("ensemble-veto-v2role-10", "ensemble: eb05-v2-role proposes, "
          "treescan vetoes below RI 10",
          "ensemble", combine="veto", veto_threshold=10.0,
          components=("eb05-v2-role", "treescan")),
    # The veto's failure mode is the vetoer's sparsity: `treescan` scores only
    # 95 of 211 substances anywhere in the sweep, and a substance it never
    # scored is floored to RI 1 and vetoed unconditionally. Reading the same
    # scan at its branches roughly doubles that coverage without lowering its
    # own line, so these ask whether aggregation is what the veto was missing.
    Model("ensemble-veto-tree-2", "ensemble: eb05-dual proposes, the full "
          "tree reading (treescan-tree) vetoes below RI 2",
          "ensemble", combine="veto", veto_threshold=2.0,
          components=("eb05-dual", "treescan-tree")),
    Model("ensemble-veto-tree-10", "ensemble: eb05-dual proposes, "
          "treescan-tree vetoes below RI 10",
          "ensemble", combine="veto", veto_threshold=10.0,
          components=("eb05-dual", "treescan-tree")),
    Model("ensemble-veto-v2role-tree-10", "ensemble: eb05-v2-role proposes, "
          "treescan-tree vetoes below RI 10 -- the deployed default's veto "
          "with the leaf-only vetoer swapped for the whole hierarchy",
          "ensemble", combine="veto", veto_threshold=10.0,
          components=("eb05-v2-role", "treescan-tree")),
    # -- nb-trend/ears ensembles, following the benchmark.md verdict that
    # `docs/LITERATURE_REVIEW.md` sec 1.5's model-the-curve family is a
    # genuinely independent axis worth the same veto/threshold-union
    # machinery already validated on treescan, not a fourth family competing
    # for the same slot. `nb-trend`'s own alarm line is a Wald z (1.96), so
    # its veto floors are chosen on that scale rather than reused from
    # treescan's RI floors -- 0.0 ("the slope is at least pointing the right
    # way") is the loose analogue of RI 2, 1.0 ("within reach of its own
    # significance line") the tighter analogue of RI 10.
    Model("ensemble-threshold-nbtrend", "ensemble, soft OR by native "
          "threshold: eb05-dual OR nb-trend, each at its own calibrated "
          "line -- the union design that worked for treescan, tried on the "
          "other independent axis",
          "ensemble", fixed_threshold=1.0, combine="threshold",
          components=("eb05-dual", "nb-trend")),
    Model("ensemble-veto-eb05dual-nbtrend-0", "ensemble: eb05-dual proposes, "
          "nb-trend vetoes below z 0 -- does an independently-significant "
          "trend corroborate eb05-dual's alarms the way treescan's veto "
          "does?",
          "ensemble", combine="veto", veto_threshold=0.0,
          components=("eb05-dual", "nb-trend")),
    Model("ensemble-veto-eb05dual-nbtrend-1", "ensemble: eb05-dual proposes, "
          "nb-trend vetoes below z 1 -- tighter floor than the 0-line above",
          "ensemble", combine="veto", veto_threshold=1.0,
          components=("eb05-dual", "nb-trend")),
    Model("ensemble-veto-nbtrend-ears-0", "ensemble: nb-trend proposes, ears "
          "vetoes below z 0 -- does the cleanest off-target detector in the "
          "table clean up nb-trend's own off-target load (Methamphetamine "
          "et al., see docs/findings/benchmark.md)?",
          "ensemble", combine="veto", veto_threshold=0.0,
          fixed_threshold=NB_TREND_THRESHOLD,
          components=("nb-trend", "ears")),
    Model("ensemble-veto-nbtrend-ears-1", "ensemble: nb-trend proposes, ears "
          "vetoes below z 1 -- tighter floor than the 0-line above",
          "ensemble", combine="veto", veto_threshold=1.0,
          fixed_threshold=NB_TREND_THRESHOLD,
          components=("nb-trend", "ears")),
    Model("ensemble-veto-v2role-nbtrend-1", "ensemble: eb05-v2-role "
          "proposes, nb-trend vetoes below z 1 -- directly comparable to "
          "ensemble-veto-v2role-10 (same proposer, treescan swapped for "
          "nb-trend as the vetoer)",
          "ensemble", combine="veto", veto_threshold=1.0,
          components=("eb05-v2-role", "nb-trend")),
]
MODELS_BY_ID = {m.id: m for m in MODELS}

# Each component's "no evidence" floor, used only to fill a (substance, as_of)
# an ensemble needs but the component's own sweep never scored -- every sweep
# here is sparse (eb05's own sweep only has ~40% substance-quarter coverage;
# see `sweep_eb05`), and an absent row means the component saw nothing worth
# reporting, not a missing observation. Ranking with it excluded instead of
# floored would let a component's sparsest quarters inflate the few
# substances it does score.
ENSEMBLE_FLOOR = {"ratio": 0.0, "eb05": 0.0, "eb05_dual": 0.0, "eb05_own": 0.0,
                  "treescan": 1.0, "treescan_branch": 1.0, "treescan_tree": 1.0,
                  "ears": 0.0, "nb_trend": 0.0, "nb_trend_own": 0.0,
                  "nb_trend_dual": 0.0, "nb_trend_emergent": 0.0,
                  "nb_trend_confirmed": 0.0,
                  "nb_trend_confirmed_emergent": 0.0}


def model_score(sweep: pd.DataFrame, model: Model) -> pd.Series:
    """The one number each model's alarm is read off, per (substance, as_of).

    `ratio` divides by an `expected` that is 0 for a substance with no
    baseline at all. That is a real and common case — it is what a genuinely
    novel substance looks like — and it is precisely where the unshrunken
    ratio has nothing to say: 3/0 and 30/0 are both infinity. Mapped to `inf`
    rather than dropped, so the model is scored on the behaviour it actually
    has (every from-nothing appearance ties at the top of its ranking)
    instead of being quietly excused from its worst regime.
    """
    if model.statistic == "ratio":
        with np.errstate(divide="ignore", invalid="ignore"):
            return sweep["n_recent"] / sweep["expected"].replace(0, np.nan)
    if model.statistic == "eb05":
        return sweep["eb05"]
    if model.statistic == "eb05_dual":
        return sweep[["eb05", "eb05_own"]].min(axis=1)
    if model.statistic == "eb05_own":
        return sweep["eb05_own"]
    if model.statistic in TREESCAN_STATISTICS:
        return sweep["recurrence_interval"]
    if model.statistic == "ears":
        return sweep["ears_z"]
    if model.statistic == "nb_trend":
        return sweep["nb_trend_z"]
    if model.statistic == "nb_trend_own":
        return sweep["nb_trend_own_z"]
    if model.statistic == "nb_trend_dual":
        return sweep[["nb_trend_z", "nb_trend_own_z"]].min(axis=1)
    if model.statistic == "nb_trend_emergent":
        return sweep["nb_trend_z"].where(sweep["emergent"], -np.inf)
    if model.statistic == "nb_trend_confirmed":
        return sweep["nb_trend_z"].where(sweep["confirmed"], -np.inf)
    if model.statistic == "nb_trend_confirmed_emergent":
        return sweep["nb_trend_z"].where(
            sweep["confirmed"] & sweep["emergent"], -np.inf)
    if model.statistic == "ensemble":
        return sweep["ensemble_score"]
    raise ValueError(f"unknown statistic {model.statistic!r}")


def build_ensemble_sweep(component_ids: tuple[str, ...],
                         sweeps: dict[tuple, pd.DataFrame],
                         combine: str = "mean",
                         veto_threshold: float | None = None) -> pd.DataFrame:
    """Combine several sweeps into one percentile-rank ensemble sweep.

    Each component's raw score is incomparable across components (EB05's
    posterior percentile, TreeScan's Monte Carlo recurrence interval, the raw
    ratio's unbounded n/E) but comparable *within* a component across
    substances in the same quarter, so the combination ranks each component's
    score among the substances it scored that quarter (0..1) before combining.
    A (substance, as_of) a component never scored is filled at that
    component's `ENSEMBLE_FLOOR` before ranking -- zero evidence, not a
    missing observation -- so sparse components (every sweep here covers well
    under the full substance panel each quarter; see `sweep_eb05`) can't have
    their few scored rows inflated by a smaller ranking denominator.

    `combine`:

    - `"mean"` — a soft AND, by percentile rank. A substance only scores high
      if every component agrees it's unusual *for that component*. Checked
      directly against the backtest: this loses to `eb05-dual` alone at a
      matched budget (9/22 vs 15/22 detected combining with `treescan`, whose
      own solo coverage is only 8/22) — averaging in a component that is
      silent on 14 of 22 known emergences drags every one of those down by
      its `ENSEMBLE_FLOOR`, even when the other component is confident. Kept
      as the negative comparison, not recommended.
    - `"max"` — a soft OR, by percentile rank. A substance scores high if
      *any* component ranks it highly *among what that component scored this
      quarter*. Checked against the backtest: still loses to `eb05-dual`
      alone at a matched budget (11/22 vs 15/22) — a within-quarter
      percentile is relative to how many substances a sparse component
      bothered to score at all, so `treescan`'s rank-1-of-3 on a quiet
      quarter counts the same as `eb05-dual`'s rank-1-of-80 on a busy one,
      which is not a fair trade against the component that is actually
      driving detection.
    - `"threshold"` — a soft OR, by each component's own already-calibrated
      line (`model.fixed_threshold`, or `ALARM_THRESHOLD` for a component
      with none): score = max over components of `raw_score / that line`, so
      1.0 is exactly where the component would fire on its own. At the fixed
      reading this is a literal alarm union, and it is the one combination
      that actually beats a solo model on the backtest: `eb05-dual` alone
      gets 0.36 mean recall / 1 off-target substance at the fixed line;
      `eb05-dual OR treescan` gets 0.42 mean recall on 2 off-target
      substances — because a union can only ever *add* alarms, a sparse
      component can never dilute what a confident one already caught, unlike
      `"mean"` and `"max"`, which both re-derive every component's opinion
      from scratch each quarter and so are exposed to how much that
      component happened to score that quarter.
    - `"veto"` — not a combination at all, an asymmetric gate: `components`
      must be exactly `(proposer, vetoer)`. The score is the proposer's own
      raw value, unchanged, wherever the vetoer's raw value is at or above
      `veto_threshold`; everywhere else it is floored to the vetoer's
      `ENSEMBLE_FLOOR` (never an alarm, however high the proposer scored).
      Unlike every mode above, the proposer's own scale is never touched, so
      it stays interpretable at the proposer's own fixed line -- this is the
      one design where "fixed reading" is not an artifact. Whether it beats
      the proposer running alone depends entirely on `veto_threshold`: set to
      the vetoer's own alarm line, this is `"threshold"`'s AND rather than
      its OR, and inherits `"mean"`'s failure mode (a sparse vetoer kills
      almost everything); set low enough, it approaches a no-op. See
      `docs/findings/benchmark.md` for where the real threshold lands.

    One further asymmetry, inherited rather than introduced here: `ratio`
    scores a substance with no baseline at all (`expected == 0`) as NaN, not
    a rank-topping `inf` (see `model_score`'s docstring vs. its actual
    behaviour) -- `rank`, `mean` and `max` all skip NaN by default, so
    `ratio` simply abstains on that substance-quarter rather than casting the
    loudest possible vote either way.
    """
    parts: dict[str, pd.Series] = {}
    for cid in component_ids:
        m = MODELS_BY_ID[cid]
        sw = sweeps[m.sweep_key]
        s = model_score(sw, m)
        parts[cid] = pd.Series(
            s.to_numpy(), index=pd.MultiIndex.from_arrays(
                [sw["substance"].to_numpy(), sw["as_of"].to_numpy()],
                names=["substance", "as_of"]))

    idx = parts[component_ids[0]].index
    for cid in component_ids[1:]:
        idx = idx.union(parts[cid].index)

    if combine == "threshold":
        normalized = {}
        for cid, s in parts.items():
            m = MODELS_BY_ID[cid]
            line = m.fixed_threshold if m.fixed_threshold is not None else ALARM_THRESHOLD
            filled = s.reindex(idx, fill_value=ENSEMBLE_FLOOR[m.statistic])
            normalized[cid] = filled / line
        agg = pd.DataFrame(normalized).max(axis=1)
        out = agg.rename("ensemble_score").reset_index()
        return out

    if combine == "veto":
        if len(component_ids) != 2:
            raise ValueError("combine='veto' takes exactly (proposer, vetoer)")
        if veto_threshold is None:
            raise ValueError("combine='veto' requires veto_threshold")
        proposer_id, vetoer_id = component_ids
        proposer_m, vetoer_m = MODELS_BY_ID[proposer_id], MODELS_BY_ID[vetoer_id]
        proposer_floor = ENSEMBLE_FLOOR[proposer_m.statistic]
        proposer = parts[proposer_id].reindex(idx, fill_value=proposer_floor)
        vetoer = parts[vetoer_id].reindex(
            idx, fill_value=ENSEMBLE_FLOOR[vetoer_m.statistic])
        vetoed = vetoer < veto_threshold
        agg = proposer.where(~vetoed, proposer_floor)
        out = agg.rename("ensemble_score").reset_index()
        return out

    ranks = {}
    for cid, s in parts.items():
        filled = s.reindex(idx, fill_value=ENSEMBLE_FLOOR[MODELS_BY_ID[cid].statistic])
        ranks[cid] = filled.groupby(level="as_of").rank(pct=True)

    combined = pd.DataFrame(ranks)
    if combine == "mean":
        agg = combined.mean(axis=1)
    elif combine == "max":
        agg = combined.max(axis=1)
    else:
        raise ValueError(f"unknown combine {combine!r}")
    out = agg.rename("ensemble_score").reset_index()
    return out


def _spans(alarms: pd.DataFrame, sub: str) -> list[tuple]:
    g = alarms[alarms["substance"] == sub].sort_values("as_of")
    return _episodes(list(g["as_of"]))


def _q_gap(a: pd.Timestamp, b: pd.Timestamp) -> int:
    """Quarters from a to b, signed."""
    return (b.year - a.year) * 4 + (b.quarter - a.quarter)


def score_model(
    sweep: pd.DataFrame, model: Model, gt: pd.DataFrame, threshold: float,
) -> tuple[dict, pd.DataFrame]:
    """Every metric for one model at one threshold.

    Detection lag is measured from each reference interval's `interval_start`
    to the first alarm quarter *inside* that interval — not to the first alarm
    anywhere, which would let a spurious alarm years earlier be scored as a
    prescient one. A model that never alarms inside the interval contributes
    no lag and drops out of `mean_lag_q`, so `detected` must be read alongside
    it: a model that detects one substance perfectly and misses four has an
    excellent mean lag.
    """
    s = model_score(sweep, model)
    alarms = sweep.loc[s > threshold, ["substance", "as_of"]]
    labelled = set(gt["substance"])
    if "no_emergence" not in gt.columns:
        gt = gt.assign(no_emergence=False)

    rows, lags = [], []
    for sub, g in gt.groupby("substance", sort=False):
        # A documented negative (e.g. Morphine -- see ground_truth.py's
        # module docstring): interval_start/interval_end hold a placeholder
        # for CSV-parsing convenience only, and must not be scored as a real
        # interval, or a genuine false alarm reads as a hit.
        gt_ivs = [] if bool(g["no_emergence"].iloc[0]) else list(
            zip(g["interval_start"], g["interval_end"]))
        det_ivs = _spans(alarms, sub)
        m = interval_overlap(gt_ivs, det_ivs)
        hits = alarms.loc[alarms["substance"] == sub, "as_of"]
        inside = [q for q in hits
                  if any(a <= q <= b for a, b in gt_ivs)]
        lag = _q_gap(min(a for a, _ in gt_ivs), min(inside)) if inside else np.nan
        lags.append(lag)
        rows.append({"substance": sub, **m, "lag_q": lag})

    per_sub = pd.DataFrame(rows)
    off = alarms[~alarms["substance"].isin(labelled)]
    pr = precision_rate(per_sub)
    return {
        "model": model.id,
        "label": model.label,
        "threshold": threshold,
        "detected": int(per_sub["lag_q"].notna().sum()),
        "n_labelled": len(per_sub),
        "mean_iou": float(per_sub["iou"].mean()),
        "mean_recall": float(per_sub["recall"].mean()),
        "mean_precision": float(per_sub["precision"].mean()),
        "precision_rate": pr["precision_rate"],
        "precision_rate_n": pr["precision_rate_n"],
        "precision_rate_ci": f"[{pr['precision_rate_lo']:.2f}, "
                             f"{pr['precision_rate_hi']:.2f}]",
        "mean_lag_q": float(np.nanmean(lags)) if per_sub["lag_q"].notna().any()
                      else np.nan,
        "alarm_quarters": int(len(alarms)),
        "off_target_quarters": int(len(off)),
        "off_target_subs": int(off["substance"].nunique()),
    }, per_sub


def matched_threshold(sweep: pd.DataFrame, model: Model, budget: int) -> float:
    """The threshold at which this model raises exactly `budget` alarms.

    Returns the largest score that still fits inside the budget, nudged down
    so the strict `>` comparison in `score_model` admits it. Ties at the
    budget boundary (common for `ratio`, whose from-nothing cases all sit at
    `inf`) are resolved in the model's favour — it gets *at least* the budget,
    never less — so no model is scored on a smaller alarm allowance than the
    baseline it is being compared against.
    """
    s = model_score(sweep, model).dropna()
    s = s[np.isfinite(s)]
    if len(s) <= budget:
        return float(s.min()) - 1e-9 if len(s) else 0.0
    cut = float(np.sort(s.to_numpy())[::-1][budget - 1])
    return np.nextafter(cut, -np.inf)


def build_sweeps(
    models: list[Model], mentions: Path, cutoff: str,
    recent_quarters: int, baseline_quarters: int, quiet: bool = False,
    counts_denom: tuple[pd.DataFrame, pd.Series] | None = None,
) -> tuple[dict[tuple, pd.DataFrame], pd.Series]:
    """One as-of sweep per distinct (case definition, role discount, spatial
    fusion), plus one per distinct set of ensemble components, and the
    quarterly denominator they were built from.

    The eight EB05-family models collapse to five sweeps -- `ratio`, `eb05`
    and `eb05-dual` are three readings of the same numbers -- which matters
    because the position-weighted sweeps cannot slice a precomputed table
    (`_ever_independent` and, when `role_discount` is on, `_role_dispersion`
    are both cross-case aggregates; see `trends._position_credit`) and so
    re-run `load_quarterly` once per as-of quarter, at ~40s each. An
    ensemble's components are built even if the caller did not select them
    directly -- `--models ensemble-eb-tree` alone still needs `eb05-dual` and
    `treescan`'s sweeps to combine.

    The three TreeScan readings are all reshapings of one cached
    `treescan backtest` table and cost no Monte Carlo work here; the two
    branch readings additionally need the mentions table and the tree to
    recompute each node's excess attribution, which is loaded at most once
    however many of them were selected.

    `counts_denom`, when given, is used in place of
    `load_quarterly(mentions, cutoff=cutoff)` -- the injection point
    `docs/SYNTHETIC_BENCHMARK_DESIGN.md` sec 6 designs, so
    `emerging.validation.synthetic.generate_panel`'s output scores through
    every existing `Model`/`score_model` path without a parallel pipeline.
    Restricted to phase-1-compatible models: TreeScan (needs the real
    substance tree) and any `weighted`/`role_discount`/spatial sweep (both
    still read real case-level or geocoded records straight from `mentions`
    regardless of `counts_denom`) would silently mix a synthetic panel with
    real county data, so selecting one alongside `counts_denom` raises
    instead of quietly doing that.
    """
    all_ids = {m.id for m in models}
    for m in models:
        all_ids.update(m.components)
    all_models = [MODELS_BY_ID[i] for i in all_ids]

    if counts_denom is not None:
        needs_real = sorted(m.id for m in all_models
                            if m.statistic in TREESCAN_STATISTICS
                            or m.weighted or m.spatial != "none")
        if needs_real:
            raise ValueError(
                "counts_denom (synthetic data) cannot be combined with "
                f"models that still read real data through `mentions`: "
                f"{needs_real} -- TreeScan, --weighted/--role-discount and "
                "spatial fusion all need real case-level or geocoded "
                "records a phase-1 synthetic panel does not have "
                "(docs/SYNTHETIC_BENCHMARK_DESIGN.md sec 5, phase 2)")
        counts, denom = counts_denom
    else:
        counts, denom = load_quarterly(mentions, cutoff=cutoff, quiet=quiet)
    quarters = sorted(denom.index)
    need_spatial = {m.spatial for m in all_models} - {"none"}
    z_sweep = (_spatial_sweep(quarters, mentions,
                              recent_quarters=recent_quarters)
               if need_spatial else None)

    out: dict[tuple, pd.DataFrame] = {}
    ensemble_keys = []
    prepared: tuple | None = None  # treescan's counts/tree, built at most once
    for key in dict.fromkeys(m.sweep_key for m in all_models):
        if key[0] == "ensemble":
            ensemble_keys.append(key)  # needs other sweeps built first
            continue
        if key[0] in TREESCAN_STATISTICS:
            stat, min_share = key
            if not quiet:
                detail = {"treescan": "substance leaves only",
                          "treescan_branch": "branch nodes only",
                          "treescan_tree": "leaves + branches"}[stat]
                guard = ("every member present" if min_share is None
                         else f"excess share >= {min_share:.0%}")
                typer.echo(f"  sweep: {stat} ({detail}, {guard}; cached "
                           "backtest, `emerging treescan backtest` to refresh)")
            if stat == "treescan":
                out[key] = treescan_leaf_sweep()
            else:
                # `prepare` re-reads the mentions table and rebuilds the tree
                # to recompute the attribution weights, so it is done once and
                # shared across however many branch readings were selected.
                if prepared is None:
                    prepared = treescan_prepare(mentions, cutoff)
                out[key] = treescan_branch_sweep(
                    mentions_path=mentions, cutoff=cutoff,
                    min_excess_share=min_share,
                    include_leaves=(stat == "treescan_tree"),
                    prepared=prepared)
            continue
        if key[0] in ("ears", "nb_trend"):
            stat, weighted, role_discount = key
            if not quiet:
                detail = ("moving-baseline z-score" if stat == "ears" else
                          "trend slope + own-count trend + emergent flag")
                typer.echo(f"  sweep: {stat} (model-the-curve family, {detail}, "
                           f"{recent_quarters}q recent / {baseline_quarters}q "
                           f"baseline, same window as eb05, weighted="
                           f"{weighted} role_discount={role_discount})")
            fn = aberration_ears_sweep if stat == "ears" else aberration_nb_trend_sweep
            out[key] = fn(counts, denom, recent_quarters=recent_quarters,
                          baseline_quarters=baseline_quarters,
                          weighted=weighted, role_discount=role_discount,
                          mentions_path=mentions)
            continue
        weighted, role_discount, spatial = key
        if not quiet:
            typer.echo(f"  sweep: weighted={weighted} role_discount="
                       f"{role_discount} spatial={spatial}"
                       + (" (recomputes per quarter, ~40s)" if weighted else ""))
        out[key] = sweep_eb05(
            counts, denom, recent_quarters, baseline_quarters,
            z_sweep=z_sweep if spatial != "none" else None,
            spatial_mode=spatial, weighted=weighted,
            role_discount=role_discount, mentions_path=mentions)

    for key in ensemble_keys:
        combine, veto_threshold, components = key[1], key[2], key[3:]
        if not quiet:
            detail = (f"{components[0]} proposes, {components[1]} vetoes "
                      f"below {veto_threshold}" if combine == "veto"
                      else f"({combine}) of {', '.join(components)}")
            typer.echo(f"  sweep: ensemble {detail}")
        out[key] = build_ensemble_sweep(components, out, combine=combine,
                                        veto_threshold=veto_threshold)
    return out, denom


@app.command("run")
def run(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    ground_truth: Path = typer.Option(KNOWN_EMERGENCES_PATH, "--ground-truth"),
    threshold: float = typer.Option(ALARM_THRESHOLD, "--threshold",
                                    help="the fixed reading line"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    recent_quarters: int = typer.Option(RECENT_QUARTERS, "--recent-quarters"),
    baseline_quarters: int = typer.Option(BASELINE_QUARTERS, "--baseline-quarters"),
    models: str = typer.Option("", "--models",
                               help="comma-separated model ids; default all"),
    baseline: str = typer.Option("eb05", "--baseline",
                                 help="model whose alarm budget the matched "
                                      "reading is cut to"),
) -> None:
    """Score every detector against the backtest, twice: at a fixed reading
    line and at a matched alarm budget."""
    chosen = ([MODELS_BY_ID[m.strip()] for m in models.split(",") if m.strip()]
              if models else MODELS)
    if baseline not in {m.id for m in chosen}:
        raise typer.BadParameter(f"--baseline {baseline} is not among the "
                                 "selected models")

    typer.echo(f"benchmarking {len(chosen)} models "
               f"({recent_quarters}q recent vs {baseline_quarters}q baseline, "
               f"cutoff {cutoff})")
    sweeps, denom = build_sweeps(chosen, mentions, cutoff, recent_quarters,
                                 baseline_quarters)
    gt = load_ground_truth(ground_truth, max(denom.index))

    fixed, per_sub = [], []
    for m in chosen:
        t = threshold if m.fixed_threshold is None else m.fixed_threshold
        row, sub = score_model(sweeps[m.sweep_key], m, gt, t)
        row["reading"] = "fixed"
        fixed.append(row)
        sub.insert(0, "model", m.id)
        sub.insert(1, "reading", "fixed")
        per_sub.append(sub)

    budget = next(r["alarm_quarters"] for r in fixed if r["model"] == baseline)
    matched = []
    for m in chosen:
        t = matched_threshold(sweeps[m.sweep_key], m, budget)
        row, sub = score_model(sweeps[m.sweep_key], m, gt, t)
        row["reading"] = "matched"
        matched.append(row)
        sub.insert(0, "model", m.id)
        sub.insert(1, "reading", "matched")
        per_sub.append(sub)

    tab = pd.DataFrame(fixed + matched)
    subs = pd.concat(per_sub, ignore_index=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out_dir / "models.csv", index=False)
    subs.to_csv(out_dir / "per_substance.csv", index=False)

    cols = ["model", "threshold", "detected", "mean_recall", "mean_iou",
            "mean_precision", "precision_rate", "precision_rate_ci",
            "mean_lag_q", "alarm_quarters",
            "off_target_quarters", "off_target_subs"]
    fmt = lambda v: f"{v:.2f}"  # noqa: E731
    typer.echo(f"\n=== fixed reading line: score > {threshold:g} "
              "(per-model override in the threshold column, e.g. treescan "
              f"at its own RI >= {TREESCAN_RI_THRESHOLD:.0f}) ===")
    typer.echo(tab[tab["reading"] == "fixed"][cols].to_string(
        index=False, float_format=fmt))
    typer.echo(f"\n=== matched alarm budget: {budget} substance-quarters "
               f"(= {baseline} at {threshold:g}) ===")
    typer.echo(tab[tab["reading"] == "matched"][cols].to_string(
        index=False, float_format=fmt))
    typer.echo(f"\ndetected = of {len(gt['substance'].unique())} labelled "
               "emergences, how many get an alarm inside their reference "
               "interval; off_target = alarms on unlabelled substances, a "
               "reviewer load, not a proven error rate")
    typer.echo(f"\nwrote {out_dir / 'models.csv'} and "
               f"{out_dir / 'per_substance.csv'}")


@app.command("plot")
def plot(
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    table: Path = typer.Option(None, "--table",
                               help="models.csv from `benchmark run`"),
) -> None:
    """Recall at the fixed line against recall at a matched alarm budget.

    The two bars per model are the entire argument of this module, so they are
    drawn side by side rather than in separate panels: the gap between them is
    the share of a model's apparent advantage that is threshold inflation
    rather than discrimination. A model whose orange bar is much shorter than
    its blue one did not detect better, it fired more.
    """
    table = table or (out_dir / "models.csv")
    if not table.exists():
        raise typer.BadParameter(f"{table} not found -- run `benchmark run`")
    tab = pd.read_csv(table)

    order = [m.id for m in MODELS if m.id in set(tab["model"])]
    fixed = tab[tab["reading"] == "fixed"].set_index("model").reindex(order)
    matched = tab[tab["reading"] == "matched"].set_index("model").reindex(order)

    y = np.arange(len(order))
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(11, 0.62 * len(order) + 2.2), sharey=True,
        gridspec_kw={"width_ratios": [2.1, 1]})

    # y is inverted below, so the *smaller* offset draws on top: the fixed
    # reading sits above the matched one for every model, matching the order
    # the two are described in.
    ax.barh(y - 0.19, fixed["mean_recall"], height=0.36, color=BLUE,
            label=f"fixed line (> {fixed['threshold'].iloc[0]:g})")
    ax.barh(y + 0.19, matched["mean_recall"], height=0.36, color=ORANGE,
            label=f"matched budget ({int(matched['alarm_quarters'].iloc[0])} "
                  "alarm-quarters)")
    for i, (f, m) in enumerate(zip(fixed["mean_recall"], matched["mean_recall"])):
        ax.text(f + 0.008, i - 0.19, f"{f:.2f}", va="center", fontsize=8,
                color=INK)
        ax.text(m + 0.008, i + 0.19, f"{m:.2f}", va="center", fontsize=8,
                color=INK)
    ax.set_yticks(y, order, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, max(fixed["mean_recall"].max(),
                       matched["mean_recall"].max()) * 1.14)
    ax.set_xlabel("mean recall over the labelled emergences")
    ax.set_title("How much of each detector's coverage survives a matched "
                 "alarm budget", fontsize=10, loc="left", color=INK)
    # Above the axes, not inside: at seven models the lower-right corner is
    # occupied by the last model's own bars.
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower left",
              bbox_to_anchor=(0, 1.02))
    ax.spines[["top", "right"]].set_visible(False)

    # The cost side, on a log axis: `ratio` runs an order of magnitude hotter
    # than everything else at the fixed line, so a linear scale would collapse
    # the six models that matter into one indistinguishable stub.
    ax2.barh(y, fixed["alarm_quarters"], height=0.55, color=MUTED)
    for i, v in enumerate(fixed["alarm_quarters"]):
        ax2.text(v * 1.08, i, f"{int(v)}", va="center", fontsize=8, color=INK)
    ax2.set_xscale("log")
    ax2.set_xlim(right=fixed["alarm_quarters"].max() * 2.2)
    ax2.set_xticks([50, 100, 200, 500, 1000])
    ax2.get_xaxis().set_major_formatter(ScalarFormatter())
    ax2.set_xlabel("alarm-quarters raised at the fixed line")
    ax2.set_title("what it costs to run it (log scale)", fontsize=10,
                  loc="left", color=INK)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    path = out_dir / "benchmark.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    typer.echo(f"wrote {path}")


if __name__ == "__main__":
    app()
