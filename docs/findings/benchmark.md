# Every detector, scored side by side

Findings for `emerging/validation/benchmark.py`. Results in
[`results/benchmark/`](../../results/benchmark/). The variants being compared
are specified in [`GPS_V2_DESIGN.md`](../GPS_V2_DESIGN.md); the reference
intervals and the overlap metric come from
[ground_truth.md](ground_truth.md), including its metric definitions. Back to
the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

**`ground_truth.md`'s own results tables show only plain EB05** — their
"Detected" column is that one model's alarm episode
(`trends.sweep_eb05`/`_episode_spans`, single-gate, unweighted, no spatial
term). This is the doc that answers "how does *every* variant do," side by
side, on the same reference.

## Why the fixed reading line is not enough on its own

`GPS_V2_DESIGN.md` accumulated three extensions to plain EB05, each validated
against plain EB05 one at a time. That does not say which detector to run. This
scores all of them on identical data with identical as-of discipline — and
**twice**, because reading them all at 1.5 turns out to be misleading.

The reason is mechanical: several of the variants change the *scale* of the
statistic, not only its ordering. The spatial discount shrinks `E`, so EB05
rises for every substance with a spatial score, and the model buys recall by
firing more often rather than by discriminating better. The matched-budget
reading re-cuts each model's threshold so it raises the **same 85
substance-quarter alarms** plain EB05 raises at 1.5. Now a model can only
improve by spending a fixed alarm budget on better substance-quarters.

The precedent is this repository's own tree-scan head-to-head
(`RESEARCH_LOG.md`, 2026-08-10), where an apparent detection advantage
vanished once the alarm budget was matched. The same thing happens here.

**Rebuilt against the expanded 22-substance reference** (was 5, then 15,
then 18, now 22 — see `ground_truth.md`'s selection-bias finding and its
false-positive-rate section). This matters here specifically: the
5-substance version had no documented false alarm anywhere in the file, so
no model's precision could ever be told apart from another's. Seven
`no_emergence` cases now exist (Morphine, Lidocaine, Levamisole, Alprazolam,
Mirtazapine, Hydromorphone, Citalopram), which is what makes `mean_precision`
and `precision_rate` below actual rates rather than one anecdote.

## Three columns `ground_truth.md` doesn't have

`mean_iou`/`mean_recall`/`mean_precision` are defined in
[ground_truth.md's Definitions](ground_truth.md#definitions) — those are
per-substance metrics, averaged over the 15 labelled ones. The three columns
below are different in kind: they're computed over **every substance in the
~205-substance catalog**, not just the labelled 15, because a model's cost
isn't only about how well it does on the substances we've checked — it's
also about how much unverified output it hands a reviewer.

| Column | Computed as | Reads as |
|---|---|---|
| `alarm_quarters` | count of `(substance, quarter)` pairs where the model's score exceeds the threshold, across the full catalog | total alarm volume — how much output this model produces |
| `off_target_quarters` | of those, how many are on a substance **not in the 22 labelled ones** | alarm-quarters with no verdict either way |
| `off_target_subs` | the count of *distinct* unlabelled substances that got at least one alarm | how many different substances a reviewer would have to go check by hand |

**These are a workload measure, not a false-positive count.** An off-target
alarm might be a real emergence that simply isn't labelled yet — Amphetamine
and 1,1-Difluoroethane remain investigated and unlabelled on exactly those
grounds (`ground_truth.md`). It's not, however, a blanket excuse: six
off-target substances from earlier passes (Lidocaine, Levamisole, Alprazolam,
Mirtazapine, Hydromorphone, Citalopram) turned out to be legitimate
`no_emergence` cases once checked directly, and are labelled negatives now,
not off-target anymore — see below. One (Codeine) turned out to be a real
positive on closer inspection. Reading `off_target_subs` as "this many wrong
alarms" overstates what the number means; reading it as "this many alarms
nobody has checked yet" is accurate.

## Results

**`precision_rate` and its 95% Wilson CI are defined in
[ground_truth.md's Definitions](ground_truth.md#definitions)** — one
Bernoulli trial per alarmed substance (did the alarm ever touch real ground
truth), not pooled over alarm-quarters, which would treat one substance's
correlated quarters as independent trials and understate the uncertainty.
It is a different, coarser question than `mean_precision` (which measures
how *much* of the alarm was real, not just whether any of it was) — reported
alongside it, not as a replacement.

**Fixed reading line, score > 1.5:**

| model | detected | mean recall | mean precision | precision rate (95% CI) | mean IoU | lag (q) | alarm-q | off-target subs |
|---|---|---|---|---|---|---|---|---|
| ratio (no shrinkage) | 15/22 | 0.64 | 0.46 | 0.71 [0.50, 0.86] | 0.30 | 1.20 | **872** | 93 |
| eb05 | 15/22 | 0.43 | 0.72 | 0.75 [0.53, 0.89] | 0.31 | 1.67 | 85 | 4 |
| eb05-weighted (#1) | **12/22** | 0.30 | **0.80** | **0.80 [0.55, 0.93]** | 0.25 | 1.67 | 62 | 0 |
| eb05-dual (#2) | 14/22 | 0.36 | 0.74 | 0.74 [0.51, 0.88] | 0.28 | 1.79 | 75 | 1 |
| eb05-spatial (#3, discount) | 15/22 | **0.48** | 0.65 | 0.68 [0.47, 0.84] | **0.31** | 1.60 | 106 | 7 |
| eb05-spatial-prior (#3, covariate) | 15/22 | 0.43 | 0.75 | 0.79 [0.57, 0.91] | **0.33** | 1.67 | 84 | 3 |
| eb05-v2 (#1+#2+#3) | **11/22** | 0.34 | 0.79 | 0.79 [0.52, 0.92] | 0.28 | **1.27** | 71 | 0 |
| eb05-v2-role (#1 ext.+#2+#3, candidate) | **11/22** | 0.34 | **0.85** | **0.85 [0.58, 0.96]** | 0.30 | 1.36 | 69 | 0 |
| eb05-weighted-role (#1, extended) | **12/22** | 0.28 | **0.86** | **0.86 [0.60, 0.96]** | 0.24 | 1.67 | 59 | 0 |
| treescan (P1, solo, leaves only) | 8/22 | 0.32 | 0.73 | 0.73 [0.43, 0.90] | 0.27 | 1.12 | 79 | 1 |
| treescan-branch (branches, unguarded) | 12/22 | 0.39 | 0.40 | 0.63 [0.41, 0.81] | 0.23 | 2.92 | **389** | 56 |
| treescan-branch-attr (branches, share ≥ 15%) | 10/22 | 0.34 | 0.70 | 0.77 [0.50, 0.92] | 0.27 | 2.10 | 105 | 7 |
| **treescan-tree (leaves + attributed branches)** | 11/22 | **0.44** | 0.70 | 0.73 [0.48, 0.89] | **0.34** | 1.64 | 121 | 7 |
| ears (EARS C2, moving-baseline z-score) | 11/22 | 0.23 | 0.79 | 0.79 [0.52, 0.92] | 0.19 | 2.09 | 43 | **0** |
| nb-trend (Poisson trend slope, Wald z) | 14/22 | **0.49** | 0.73 | 0.78 [0.55, 0.91] | **0.37** | 1.79 | 161 | 11 |
| nb-trend-own (own-count trend, no denominator) | 14/22 | 0.53 | 0.64 | 0.70 [0.48, 0.85] | 0.36 | 1.43 | 299 | 25 |
| nb-trend-dual (min of the two above) | 14/22 | 0.45 | 0.74 | 0.78 [0.55, 0.91] | 0.35 | 1.86 | 146 | 9 |
| nb-trend-emergent (gated by the EMERGENT_THRESHOLD analogue) | 10/22 | 0.31 | 0.67 | 0.71 [0.45, 0.88] | 0.25 | 1.50 | 69 | **7** |
| nb-trend-confirmed (>=3 consecutive alarming quarters, synthetic-benchmark fix) | 12/22 | 0.26 | 0.85 | 0.86 [0.60, 0.96] | 0.24 | 3.75 | 87 | **4** |
| nb-trend-confirmed-emergent (both gates at once) | 8/22 | 0.15 | 0.80 | 0.80 [0.49, 0.94] | 0.14 | 3.38 | 28 | **2** |
| ears-weighted-role (#1 extended, ported) | 11/22 | 0.22 | 0.79 | 0.79 [0.52, 0.92] | 0.18 | 2.09 | 40 | **0** |
| nb-trend-weighted-role (#1 extended, ported) | 14/22 | 0.51 | 0.73 | 0.78 [0.55, 0.91] | **0.39** | 1.71 | 157 | 16 |
| nb-trend-dual-weighted-role (#1 extended, ported) | 14/22 | 0.49 | 0.75 | 0.78 [0.55, 0.91] | 0.38 | 1.71 | 137 | 13 |
| ensemble-threshold-nbtrend (eb05-dual OR nb-trend) | 14/22 | **0.52** | 0.63 | 0.67 [0.45, 0.83] | 0.36 | 1.50 | 170 | 12 |
| ensemble-veto-eb05dual-nbtrend-0 (eb05-dual, nb-trend vetoes < z 0) | 14/22 | 0.36 | 0.74 | 0.74 [0.51, 0.88] | 0.28 | 1.79 | 75 | 1 |
| ensemble-veto-eb05dual-nbtrend-1 (eb05-dual, nb-trend vetoes < z 1) | 14/22 | 0.36 | 0.78 | 0.78 [0.55, 0.91] | 0.28 | 1.79 | 73 | **0** |
| ensemble-veto-nbtrend-ears-0 (nb-trend, ears vetoes < z 0) | 14/22 | 0.49 | 0.73 | 0.78 [0.55, 0.91] | 0.37 | 1.79 | 161 | 11 |
| ensemble-veto-nbtrend-ears-1 (nb-trend, ears vetoes < z 1) | 14/22 | 0.45 | 0.73 | 0.78 [0.55, 0.91] | 0.34 | 1.79 | 143 | 8 |
| ensemble-veto-v2role-nbtrend-1 (eb05-v2-role, nb-trend vetoes < z 1) | 11/22 | 0.34 | 0.85 | 0.85 [0.58, 0.96] | 0.30 | 1.36 | 69 | **0** |
| ensemble-mean (eb05-dual + treescan, rank) | 15/22 | 0.94 | 0.43 | 0.68 [0.47, 0.84] | 0.40 | 0.13 | 1585 | 94 |
| ensemble-max (eb05-dual + treescan, rank) | 15/22 | 0.98 | 0.41 | 0.68 [0.47, 0.84] | 0.40 | 0.00 | 1811 | 127 |
| ensemble-threshold (eb05-dual OR treescan) | 14/22 | 0.42 | 0.70 | 0.70 [0.48, 0.85] | 0.31 | 1.71 | 96 | 2 |
| ensemble-threshold-ratio (+ ratio) | 15/22 | 0.71 | 0.44 | 0.68 [0.47, 0.84] | 0.34 | 0.73 | 894 | 94 |
| ensemble-veto-2 (eb05-dual, treescan vetoes < RI 2) | 14/22 | 0.36 | **0.78** | **0.78 [0.55, 0.91]** | 0.28 | 1.79 | 73 | **0** |
| ensemble-veto-10 (eb05-dual, treescan vetoes < RI 10) | 12/22 | 0.33 | 0.80 | 0.80 [0.55, 0.93] | 0.27 | 1.75 | 68 | **0** |
| ensemble-veto-tree-2 (eb05-dual, treescan-tree vetoes < RI 2) | 14/22 | 0.36 | **0.78** | **0.78 [0.55, 0.91]** | 0.28 | 1.79 | 73 | **0** |
| ensemble-veto-tree-10 (eb05-dual, treescan-tree vetoes < RI 10) | 12/22 | 0.33 | 0.80 | 0.80 [0.55, 0.93] | 0.28 | 1.67 | 69 | **0** |

**The rank-based ensembles' fixed-line numbers are not meaningful, and the
table above is what demonstrates it, not just a caveat about it** — mean/max
of a percentile rank clusters near 0.5 (mean of `k` uniforms) or well above
it (max of `k` uniforms) *by construction*, regardless of any real signal, so
"score > 1.5" reads as "score > 0.5" for these — see the ensembles section
below for what these rows actually show.

**Matched alarm budget, 85 substance-quarters:**

| model | threshold | detected | mean recall | mean precision | precision rate (95% CI) | mean IoU | lag (q) | off-target subs |
|---|---|---|---|---|---|---|---|---|
| ratio | 4.78 | 10/22 | 0.26 | 0.69 | 0.77 [0.50, 0.92] | 0.19 | 2.20 | 28 |
| eb05 | 1.50 | 15/22 | 0.43 | 0.72 | 0.75 [0.53, 0.89] | 0.31 | 1.67 | 4 |
| eb05-weighted | 1.28 | 14/22 | 0.39 | 0.78 | 0.78 [0.55, 0.91] | 0.31 | 1.71 | 4 |
| eb05-dual | 1.43 | **15/22** | **0.43** | 0.73 | 0.75 [0.53, 0.89] | 0.32 | **1.60** | 3 |
| eb05-spatial | 1.70 | 14/22 | 0.39 | 0.74 | 0.74 [0.51, 0.88] | 0.31 | 1.71 | 4 |
| eb05-spatial-prior | 1.49 | 15/22 | 0.43 | 0.75 | 0.79 [0.57, 0.91] | 0.33 | 1.67 | 4 |
| eb05-v2 | 1.39 | 12/22 | 0.37 | 0.74 | 0.75 [0.51, 0.90] | 0.29 | 1.42 | 4 |
| eb05-v2-role (candidate) | 1.39 | 12/22 | 0.37 | 0.79 | 0.80 [0.55, 0.93] | 0.31 | 1.42 | 4 |
| eb05-weighted-role | 1.24 | 13/22 | 0.37 | 0.76 | 0.76 [0.53, 0.90] | 0.30 | 1.69 | 6 |
| treescan | RI 60.2 | 8/22 | 0.36 | 0.67 | 0.67 [0.39, 0.86] | 0.28 | 1.00 | 1 |
| treescan-branch | RI 10,000 (uncuttable) | 12/22 | 0.31 | 0.50 | 0.63 [0.41, 0.81] | 0.19 | 3.17 | 30 |
| treescan-branch-attr | RI 500 | 9/22 | 0.30 | 0.67 | 0.75 [0.47, 0.91] | 0.25 | 2.11 | 3 |
| **treescan-tree** | RI 3,333 | 10/22 | 0.36 | 0.77 | 0.77 [0.50, 0.92] | 0.30 | 1.60 | 2 |
| ears | 2.08 | 14/22 | 0.39 | 0.73 | 0.78 [0.55, 0.91] | 0.30 | 1.43 | 5 |
| nb-trend | 2.61 | 11/22 | 0.25 | 0.91 | **0.92 [0.65, 0.99]** | 0.23 | 2.09 | 4 |
| nb-trend-own | 4.32 | 4/22 | 0.09 | 0.66 | 0.67 [0.30, 0.90] | 0.08 | 1.75 | 4 |
| nb-trend-dual | 2.48 | 11/22 | 0.25 | 0.91 | **0.92 [0.65, 0.99]** | 0.23 | 2.00 | 4 |
| nb-trend-emergent | 1.96 | 10/22 | 0.31 | 0.67 | 0.71 [0.45, 0.88] | 0.25 | 1.50 | 7 |
| nb-trend-confirmed | 2.01 | 12/22 | 0.26 | 0.85 | **0.92 [0.67, 0.99]** | 0.24 | 3.75 | 4 |
| nb-trend-confirmed-emergent | 1.96 | 8/22 | 0.15 | 0.80 | 0.80 [0.49, 0.94] | 0.14 | 3.38 | **2** |
| ears-weighted-role | 2.17 | 14/22 | 0.38 | 0.73 | 0.78 [0.55, 0.91] | 0.28 | 1.50 | 7 |
| nb-trend-weighted-role | 2.54 | 12/22 | 0.29 | 0.80 | 0.80 [0.55, 0.93] | 0.24 | 2.17 | 3 |
| nb-trend-dual-weighted-role | 2.37 | **13/22** | 0.30 | 0.81 | 0.81 [0.57, 0.93] | 0.25 | 2.08 | 3 |
| ensemble-threshold-nbtrend | 1.44 | 11/22 | 0.31 | 0.91 | **0.92 [0.65, 0.99]** | 0.29 | 1.73 | **2** |
| ensemble-veto-eb05dual-nbtrend-0 | 1.43 | **15/22** | 0.43 | 0.72 | 0.75 [0.53, 0.89] | 0.32 | 1.60 | 3 |
| ensemble-veto-eb05dual-nbtrend-1 | 1.37 | **15/22** | **0.44** | 0.72 | 0.75 [0.53, 0.89] | 0.32 | 1.67 | 2 |
| ensemble-veto-nbtrend-ears-0 | 2.61 | 11/22 | 0.25 | 0.91 | **0.92 [0.65, 0.99]** | 0.23 | 2.09 | 4 |
| ensemble-veto-nbtrend-ears-1 | 2.60 | 11/22 | 0.25 | 0.91 | **0.92 [0.65, 0.99]** | 0.23 | 2.09 | 4 |
| ensemble-veto-v2role-nbtrend-1 | 1.38 | 12/22 | 0.37 | 0.79 | 0.80 [0.55, 0.93] | 0.31 | 1.42 | 3 |
| ensemble-mean | rank 0.98 | 9/22 | 0.30 | 0.67 | 0.69 [0.42, 0.87] | 0.24 | 1.33 | 2 |
| ensemble-max | rank 0.99 | 11/22 | 0.33 | 0.69 | 0.73 [0.48, 0.89] | 0.26 | 1.73 | 3 |
| ensemble-threshold | 1.25× | 12/22 | 0.38 | 0.80 | 0.80 [0.55, 0.93] | 0.32 | 1.58 | 1 |
| ensemble-threshold-ratio | 4.80× | 11/22 | 0.37 | 0.74 | 0.79 [0.52, 0.92] | 0.31 | 1.45 | 7 |
| ensemble-veto-2 | 1.37 | **15/22** | **0.44** | 0.72 | 0.75 [0.53, 0.89] | 0.32 | 1.60 | 2 |
| ensemble-veto-10 | 1.12 | 14/22 | 0.42 | 0.77 | 0.78 [0.55, 0.91] | **0.33** | 1.79 | **0** |
| ensemble-veto-tree-2 | 1.37 | **15/22** | **0.44** | 0.72 | 0.75 [0.53, 0.89] | 0.32 | 1.60 | 2 |
| ensemble-veto-tree-10 | 1.23 | 14/22 | 0.41 | 0.73 | 0.74 [0.51, 0.88] | 0.30 | 1.71 | 1 |

**TreeScan run solo (`model_score`'s `"treescan"` statistic — the substance-leaf
node's recurrence interval, no EB05 underneath it at all) is the cleanest
detector in the table by off-target load — 1 unlabelled substance ever
alarmed, at either reading, against 3-7 for every EB05 variant — but the
lowest `detected` (8/22) and, distinctively, the *lowest lag* on what it does
catch (1.00-1.12q, best in the table). Read together with `off_target`: this
is a genuinely different failure shape from the EB05 family, not a strictly
worse version of it — fast and quiet when it fires, but conservative about
firing. Missed entirely: Mitragynine, Xylazine, Benzylfentanyl, Clonazepam,
Ketamine, Ephedrine, Codeine, the same low-n substances the EB05 family
already struggles with, plus a few more — consistent with the leaf-only test
losing TreeScan's actual selling point (branch aggregation power for exactly
these thin substances; see `treescan.py`'s module docstring and the
nitazenes example) by construction, since this row only reads the
substance-level node. **That is no longer a caveat but a measurement**: the
three `treescan-branch*`/`treescan-tree` rows read the same scan at its
internal nodes and recover Xylazine, Clonazepam and Codeine from exactly
this missed list — see [Reading the scan at its
branches](#reading-the-scan-at-its-branches--the-aggregation-this-table-used-to-omit)
for what that costs. **It does not avoid Lidocaine** — 5 alarmed quarters at
the fixed line, more than any EB05 variant including plain `eb05` (4q) —
because nothing in the scan's case definition discounts a mention's position
on the cause line; `--role-discount` (#1, extended) is solving a problem
TreeScan's own multiplicity correction does not touch at all.

## Is `eb05-v2` still the best-of-everything combination?

No — `--role-discount` (#1, extended) postdates `eb05-v2`'s definition, and
`eb05-v2-role` (same #1+#2+#3 recipe, with `--role-discount` folded into #1)
**Pareto-dominates it: equal or better on every column, at both readings.**
Fixed line: same 11/22 detected, same 0.34 mean recall, but precision rises
0.79 → 0.85 and mean IoU 0.28 → 0.30. Matched budget: same 12/22 detected,
same 0.37 mean recall, precision 0.75 → 0.80 and mean IoU 0.29 → 0.31. The
only column that moves the wrong way is mean lag (1.27 → 1.36 at the fixed
line), and that comparison is over only the 11 substances each one catches —
too thin a base to read as a real cost. `eb05-v2` is kept in the table as the
`GPS_V2_DESIGN.md` #1+#2+#3 checkpoint it was defined as; `eb05-v2-role` is
the current recommendation whenever that combination is wanted, and the
checkpoint's own documentation should be read as superseded on this point.

## Is `eb05-v2-role`'s spatial term (#3) pulling its weight, and is the dual gate itself an ensemble?

Two follow-on questions, both checked directly rather than argued from
intuition. **`eb05-v2-role` is not a post-hoc ensemble the way the TreeScan
combinations above are.** `--weighted`/`--role-discount` are applied once, to
the substance's own case counts, in `load_quarterly` -- before either GPS
gate is fit. Both `eb05` (Gate A, share) and `eb05_own` (Gate B,
self-referential) are then computed from that *same* already-reweighted
`counts` table, and only combined at the very end via `min()`. There is one
case definition, run through GPS twice with two different `expected`
formulas, not two independently-scored detectors glued together after the
fact.

**Is #3 (spatial) worth keeping once #1-extended and #2 are both in, or is it
safe to drop as the more speculative mechanism?** `eb05-dual-role` (`#1
extended + #2`, no spatial) against `eb05-v2-role` (same, `+ #3`):

| reading | model | detected | mean recall | mean IoU | precision rate | off-target subs |
|---|---|---|---|---|---|---|
| fixed | eb05-dual-role (no #3) | 11/22 | 0.26 | 0.23 | 0.85 [0.58, 0.96] | 0 |
| fixed | eb05-v2-role (+ #3) | 11/22 | **0.34** | **0.30** | 0.85 [0.58, 0.96] | 0 |
| matched | eb05-dual-role (no #3) | **13/22** | 0.37 | 0.29 | 0.76 [0.53, 0.90] | 5 |
| matched | eb05-v2-role (+ #3) | 12/22 | 0.37 | **0.31** | **0.80 [0.55, 0.93]** | **4** |

Not a clean Pareto win either way. Dropping spatial costs real fixed-line
recall and IoU (0.34→0.26, 0.30→0.23 — the biggest single-column gap in this
section) and matched-budget precision/off-target, but *gains* one extra
detected substance at the matched budget (13 vs 12). Both `precision_rate`
CIs overlap the same way every pair in this file does at n=22, so "spatial
adds real recall" is a defensible point-estimate reading, not a proven one —
but the numbers don't support the intuition that #3 is inert or safely
droppable; if anything it's still the axis with the largest single measured
effect in this comparison. Keep it unless the matched-budget detection count
specifically is what a deployment cares about.

**Is Gate A (the share comparison) still pulling weight once role/position
weighting is already in the case definition — could the detector just be the
self-referential gate alone, `eb05_own`, weighted and role-discounted?**
Tested directly as `eb05-own-role`: **no, clearly not.** Dropping Gate A
does not simplify the detector for free — it makes it *worse* on every axis
that matters: precision rate falls to 0.73 (fixed) from 0.85, and off-target
load explodes to 10 substances / 105 alarm-quarters (fixed line), the worst
in this entire section by a wide margin, worse even than plain `eb05-dual`'s
1 off-target substance. Consistent with why Gate A exists at all
(`docs/findings/benchmark.md`'s cocaine example, and the RI writeup
above): it is not redundant with a role-discounted own-count comparison — a
substance can have a genuinely elevated *absolute* count relative to its own
history without that being surveillance-relevant, and Gate A's share
comparison is an independent check the self-referential gate has no way to
provide. The dual gate's AND is doing real, non-redundant work; "just the
self-tracked amount" is a strictly worse detector, not a leaner version of
the same one.

## Ensembling `eb05-dual` and TreeScan

Three ways of combining `eb05-dual` (the recall/IoU leader in the original
six) and `treescan` (the cleanest off-target load, and a genuinely
independent method — no EB05 underneath it) into one detector, built in
`build_ensemble_sweep` — full mechanism and the failure diagnosis for each in
its docstring, summarised here against the actual backtest numbers.

**`"mean"` and `"max"` (percentile rank, soft AND / soft OR) both lose to
`eb05-dual` running alone, and the fixed-line numbers above (15/22 detected,
recall 0.94–0.98) are an artifact, not a result** — ranking each component's
score among only the substances *that component scored that quarter* means a
sparse quarter for `treescan` (which might score 3 substances) and a busy one
for `eb05-dual` (which might score 80) count on the same 0–1 scale, so the
mean/max of a rank clusters near or above 0.5 regardless of whether anything
real happened. At the matched 85-alarm budget — the only reading that means
anything here — `"mean"` gets 9/22 detected and `"max"` gets 11/22, both
below `eb05-dual` alone's 15/22. `"mean"`'s specific failure: a substance
`treescan` is silent on gets floored and dragged down even when `eb05-dual`
is confident, since 14 of the 22 known emergences are outside `treescan`'s
own solo coverage. `"max"` fixes that dilution but not the denominator
problem — a substance that happens to be `treescan`'s top-ranked-of-3 that
quarter competes evenly against `eb05-dual`'s top-ranked-of-80.

**`"threshold"` — each component read at its own already-calibrated line
(`eb05-dual` > 1.5, `treescan` RI > 100), combined as a literal OR — is the
one combination that adds real value, and it is a precision trade, not a
free win.** At the fixed line it nearly matches plain `eb05` (14/22 detected
vs 15/22, 0.42 vs 0.43 mean recall) while cutting off-target substances in
half (2 vs 4) — a union can only ever *add* alarms, so `treescan`'s
contribution never dilutes what `eb05-dual` already caught, unlike the two
rank-based combinations above. At the matched budget it does **not** beat
`eb05-dual` alone on detection (12/22 vs 15/22, 0.38 vs 0.43 mean recall) —
some of the fixed 85-alarm budget is necessarily spent on `treescan`'s
picks instead of `eb05-dual`'s own best ones — but it posts the **best
precision rate (0.80) and cleanest off-target load (1 substance) of every
model in the entire table**, EB05 family included, tied only by
`eb05-v2-role`. Read together with `eb05-weighted-role`'s trade-off above:
this project now has two independent ways to buy precision at a real,
measured recall cost, from two structurally different mechanisms (a
per-substance case-definition discount, vs. an independent second detection
method), not one technique restated twice.

**Adding `ratio` as a third component (`ensemble-threshold-ratio`) makes
every reading worse**, confirming the "floor, not a candidate" framing
`Model`'s own docstring already gives it: fixed-line off-target jumps to 94
substances (`ratio`'s own indiscriminate 872-alarm-quarter pathology floods
back in through the union), and at the matched budget off-target rises from
1 to 7 substances for no recall gain (0.37 vs `"threshold"`'s 0.38). Not
included going forward.

### A fourth design: EB05 proposes, TreeScan vetoes

Every design above *combines* two scores into one, which is exactly what
costs `"mean"`/`"max"` their recall (a sparse component's opinion counts as
much as a dense one's) and costs `"threshold"` some of `eb05-dual`'s own best
picks (the union's alarm budget has to be shared). `combine="veto"` in
`build_ensemble_sweep` does not combine at all: it keeps `eb05-dual`'s own
raw score exactly as-is, and only asks a yes/no question of `treescan` —
does this substance clear a *veto floor*, well below `treescan`'s own RI 100
alarm line — before letting the alarm stand. A substance `eb05-dual` likes
that `treescan` sees literally no corroborating rise for at all gets
suppressed; everything else passes through untouched.

**This is the one ensemble design that matches or beats `eb05-dual` running
alone, at both readings, not just one.** At `veto_threshold = 2` (RI ≥ 2 —
twice as likely as chance, a very low bar, far below `treescan`'s own 100):

| reading | model | detected | mean recall | mean IoU | precision rate | off-target subs |
|---|---|---|---|---|---|---|
| fixed | eb05-dual alone | 14/22 | 0.36 | 0.28 | 0.74 [0.51, 0.88] | 1 |
| fixed | ensemble-veto-2 | 14/22 | 0.36 | 0.28 | **0.78 [0.55, 0.91]** | **0** |
| matched | eb05-dual alone | 15/22 | 0.43 | 0.32 | 0.75 [0.53, 0.89] | 3 |
| matched | ensemble-veto-2 | **15/22** | **0.44** | 0.32 | 0.75 [0.53, 0.89] | **2** |

Detected, recall and IoU are tied or fractionally better at both readings —
within-noise, not a real recall gain — while off-target drops to fully clean
at the fixed line and by one substance at matched budget. Because the veto
only ever *removes* alarms `treescan` flatly disagrees with, it can't dilute
`eb05-dual`'s own best picks the way `"threshold"`'s union has to when the
alarm budget is shared; and because the floor is a fixed absolute value
rather than a rank, it doesn't inherit `"mean"`'s problem either. Checked
directly: 2 of `eb05-dual`'s 75 fixed-line alarms get vetoed, both on
substances outside the 22-substance reference — real off-target alarms
`treescan` never independently saw anything happening on. At the higher
`veto_threshold = 10`, the trade becomes explicit rather than free: 12/22 →
14/22 detected (down from `eb05-dual` alone) but **off-target hits zero at
both readings** and matched-budget mean IoU (0.33) is the best of any model
in the entire table, EB05 family included. `ensemble-veto-2` is the
recommendation if the goal is "clean up `eb05-dual` for nothing"; either
veto variant avoids 3-4 of the 7 labelled negatives (Morphine, Alprazolam,
Hydromorphone, and Mirtazapine at `veto_threshold = 10`) without touching
Lidocaine at all — `treescan` independently agrees Lidocaine is elevated (its
own solo RI clears even the RI ≥ 10 bar easily), so the veto correctly leaves
it alone; that specific case stays `--role-discount`'s to solve, not this
one's, and the two remain complementary rather than redundant fixes.

**Stacking the veto on top of `eb05-v2-role` (the best solo detector, not
plain `eb05-dual`) instead — `ensemble-veto-v2role-*` — answers "is
`eb05-v2-role` alone still the ceiling, or does the independent cross-check
add anything on top of it too."** It does, and at `veto_threshold = 10` the
gain is close to unambiguous at the matched budget: **13/22 detected (up
from 12), mean recall 0.39 (up from 0.37), mean IoU 0.32 (up from 0.31),
precision rate 0.81 (up from 0.80), off-target 1 substance (down from 4)** —
every column ties or improves against `eb05-v2-role` running alone. At the
fixed line it also posts **6 of 7 negatives avoided, the best count anywhere
in this investigation** (only Citalopram still alarms) — Lidocaine included,
because `eb05-v2-role` already suppresses it upstream via `--role-discount`,
so there's nothing left for the veto to have to catch there. This is *not*,
however, evidence that layering the veto is a strictly-better move in
general: compared against `ensemble-veto-10` (the plain-`eb05-dual` version),
it trades away real detections for the extra precision (13/22 vs 14/22
detected, 0.39 vs 0.42 recall) in exchange for one more avoided negative and
a slightly better precision rate. Compounding every precision-favouring
mechanism at once (role discount + spatial + dual + veto) keeps paying for
itself here, but each layer is still buying precision with detections,
same as every other result in this file — there is no configuration in this
table that is simply "best," only points further along the same recall/
precision frontier every other result on this page sits on.

**`eb05-weighted-role` has the highest precision in the entire table at the
fixed line (0.86), and the highest `detected` among the position-weighted
family (tied with plain `eb05-weighted` at 12/22) — the first #1 variant
that improves precision without also being the worst detector.** It's an
extension of `--weighted`, not a new axis: `_role_dispersion` adds a
per-substance discount, `phi_s`, computed from the substance's own aggregate
mean credit (not just each individual death's own position), specifically to
catch the pattern `--weighted` alone cannot — see `GPS_V2_DESIGN.md` #1 for
the full mechanism and the real bug (a wrong reference point that broke
Bromazolam and every other known emergence) caught by running `backtest`
before this shipped. Checked directly in `per_substance.csv`: Lidocaine goes
from 1 alarmed quarter under `eb05-weighted` to **0** under
`eb05-weighted-role` — full suppression, not a reduction — while every other
negative is untouched (they're different failure shapes this mechanism
wasn't built for: Alprazolam is flat rather than adulterant-patterned,
Citalopram is just noisy).

**The matched-budget row shows the real cost.** At a fixed threshold,
`eb05-weighted-role` is a clean improvement on `eb05-weighted`. Forced to
match an 85-alarm budget, its threshold has to drop to 1.24 — the lowest in
the whole table — because a role-discounted substance needs more raw
evidence to clear any given bar, so hitting the same alarm count means
lowering the bar for everything else too. `detected` and mean recall both
end up slightly below plain `eb05-weighted` there. Read together: this
extension makes the *fixed-threshold* deployment strictly better on
Lidocaine specifically, without pretending the trade-off disappears once the
alarm budget is held fixed instead.

**Every one of these confidence intervals overlaps every other one.** That
is itself the headline result of adding the Wilson column: at n=22 (17-20
alarmed substances per model), this benchmark cannot statistically
distinguish any of the seven models' false-positive rates from each other.
The point estimates still have a consistent ordering worth reading (below),
but none of the gaps in this table would survive being asked "is this
difference real."

## Reading the scan at its branches — the aggregation this table used to omit

Every TreeScan row above this section is the **substance-leaf** node's
recurrence interval, which is the one node level that has an EB05
counterpart and therefore the only one the 22-substance reference can score
directly. It is also the one level at which a tree scan has no advantage over
testing substances one at a time — the whole argument of
`emerging/analysis/treescan.py` is that a *family* of thin substances is a
single hypothesis with their pooled count. Earlier versions of this file
flagged that as an untested caveat. `treescan.branch_sweep` now closes it by
propagating each internal node's score back down to the substances inside it,
onto the same `(substance, as_of)` grid every other model is scored on.

**Two guards decide whether a branch's score may be credited to a member**,
and the three new rows exist to price them:

- **Presence** — the substance must have at least one case inside the window
  that fired. This is the vectorized form of `treescan.backtest`'s
  `as_of >= first_seen` rule, and without it `root` alone hands its score to
  all 211 substances.
- **Attribution** — the substance must supply at least 15% of the node's
  excess (`excess_share`, the same quantity `backtest` prints, and the same
  15% it already calls "partly attributable"). `treescan-branch` runs with
  presence only; `treescan-branch-attr` and `treescan-tree` add attribution.

### What aggregation catches that the leaf reading cannot

Three labelled emergences the leaf-only row misses outright, each via a real
pharmacological family and each with the attribution to back it:

| substance | leaf | branch that catches it | excess share | quarters |
|---|---|---|---|---|
| **Xylazine** (n=9 lifetime) | never | `Cutting agents and adulterants` (5 leaves) | 21%, 19% | 2024Q4–2025Q1 |
| **Clonazepam** | never | `Benzodiazepines` (19), `Prescription benzodiazepines` (11) | 16–22% | 2019Q4–2021Q1, 2025Q1–Q2 |
| **Codeine** | never | `Prescription opioids` (13) | 18% | 2025Q1 |

Xylazine is the case this whole mechanism exists for, and it is the same
substance the README names as beyond reach of the *spatial* term too: at 9
lifetime mentions it scores identically zero in all 45 quarters of
`core/concentration.py` and never clears a leaf test in any quarter of the
scan. Pooled into a five-member adulterant family it clears RI 10,000, with
21% of that family's excess attributable to xylazine specifically — the
"partly attributable" band, not a clean catch, and worth quoting with that
qualifier attached.

The sharper demonstration is the **2020 designer-benzodiazepine wave**, where
aggregation does not just add a detection but moves one earlier. Etizolam and
Flualprazolam both detect at lag **0q** through `Designer benzodiazepines`
(8 leaves) against 1q at the leaf, at 59–64% and 27–46% attribution
respectively — neither substance needs the other's count to clear on its own,
but the family node sees them a quarter before either leaf does, and the
attribution says both are genuinely driving it. Bromazolam's leaf goes quiet
after 2024Q4 and the same `Benzodiazepines` node carries it through 2025Q2.

### What it costs, named

**Alprazolam — a labelled negative — alarms 9 times, and the attribution
guard cannot stop it.** It rides `Benzodiazepines` and then `Depressants` at
34–46% excess share, comfortably above any threshold that keeps Clonazepam
(16–22%) or Xylazine (19–21%). This is the honest limit of the guard:
`excess_share` measures *contribution*, not causation, and a large stable
member of a genuinely rising family will always look attributable. The same
mechanism that produces the Clonazepam catch produces the Alprazolam false
alarm; there is no threshold that separates them, which is why
`treescan-tree` avoids only 3 of the 7 negatives against the leaf reading's
4.

**Unguarded branch propagation is not a detector at all.** `treescan-branch`
raises 389 alarm-quarters across 56 off-target substances at the fixed line —
mean precision 0.40, the worst of anything in this file except the raw ratio
— and it **cannot be read at a matched budget at all**: 234 of its rows tie
at RI 10,000, the Monte Carlo ceiling at 9,999 replicates, because one node
firing at p = 1e-4 stamps that same score onto every member present. The
matched-budget row reports 234 alarms against an 85 budget for that reason,
and its numbers should be read as "uncuttable", not as a threshold choice.
The guard removes 284 of those 389 alarms and 49 of those 56 off-target
substances while costing 2 of 12 detections. That trade is the entire case
for computing attribution rather than trusting the tree.

### The whole tree, read as one detector

`treescan-tree` is `max(leaf, attributed branch)` — TreeScan as it is
actually specified, rather than either half of it. Against the leaf-only row:

| | leaf only | whole tree | |
|---|---|---|---|
| detected (fixed / matched) | 8 / 8 | **11 / 10** | aggregation reaches three more |
| mean recall (fixed / matched) | 0.32 / 0.36 | **0.44** / 0.36 | **no matched-budget gain** |
| mean IoU (fixed / matched) | 0.27 / 0.28 | **0.34** / **0.30** | |
| precision rate (matched) | 0.67 | **0.77** | |
| off-target subs (fixed / matched) | 1 / 1 | 7 / 2 | |
| mean lag (fixed / matched) | **1.12** / **1.00** | 1.64 / 1.60 | branches arrive later than leaves |

**The result that matters is the third row.** Mean recall at a matched
budget is identical (0.36 both ways) while `detected` rises 8 → 10. At a
fixed alarm budget, aggregation does not deepen coverage of the emergences
the scan already found — it *spreads the same budget across more of them*,
and the matched-budget threshold has to climb from RI 60 to RI 3,333 to pay
for it. Bromazolam (1q → 2q) and Flualprazolam (0q → 1q) each give back a
quarter of lead time at the matched budget to fund Xylazine and Clonazepam.
Whether that is an improvement depends entirely on whether breadth or depth
is what the deployment needs, and this benchmark does not decide that.

What did *not* change: **the nitazenes still never signal.** `Novel synthetic
opioids` peaks at RI 1.16 across the entire 45-quarter sweep, and the
`Nitazenes` node never fires either. That family is the case
`core/tree.py`'s docstring gives for building a family layer at all, and
aggregation does not rescue it — three mentions pooled is still three
mentions. Branch aggregation helps families whose members already have real
counts (benzodiazepines, adulterants, prescription opioids); it does not
manufacture power where there is none.

### Does a less sparse vetoer improve the veto?

The veto's known failure mode is the vetoer's sparsity: `treescan` scores
only 95 of 211 substances anywhere in the sweep, and anything it never scored
is floored to RI 1 and vetoed unconditionally. Reading the same scan at its
branches nearly doubles that coverage (162 substances) without lowering its
own line, so the obvious question is whether that is what the veto was
missing. **It is not.**

| | vetoer = `treescan` | vetoer = `treescan-tree` |
|---|---|---|
| `eb05-dual`, veto < RI 2 (matched) | 15/22, recall 0.44, IoU 0.32 | 15/22, recall 0.44, IoU 0.32 — *identical* |
| `eb05-dual`, veto < RI 10 (matched) | 14/22, recall 0.42, **IoU 0.33**, **0 off-target** | 14/22, recall 0.41, IoU 0.30, 1 off-target |
| `eb05-v2-role`, veto < RI 10 (fixed) | **6 of 7 negatives avoided** | 5 of 7 — Alprazolam gets through |

At RI 2 the two vetoers are indistinguishable, because a threshold that low
is nearly a no-op either way. At RI 10 the tree vetoer is consistently
*slightly worse*: a vetoer that scores more substances vetoes fewer alarms,
and the alarms it now lets through are the ones the veto existed to kill —
including Alprazolam on the deployed default. The extra coverage the branch
reading buys is coverage of substances riding large families, which is
exactly the population a precision-oriented veto wants to keep blocking.
`ensemble-veto-v2role-10` keeps the leaf-only vetoer, and that choice is now
checked rather than inherited.

## A different family entirely: EARS and NB-trend, from `docs/LITERATURE_REVIEW.md` sec 1.5

Every model above this section is a **disproportionality** statistic — a
ratio of two period counts, shrunk or not. `docs/LITERATURE_REVIEW.md` sec
1.5 surveys a parallel pharmacovigilance/syndromic-surveillance literature
that instead fits the *shape* of the quarterly curve and flags an excessive
derivative — CDC's EARS aberration detectors, the EudraVigilance
report-increase algorithm, MaxSPRT read as a slope test. Two candidates cheap
enough to prototype without new infrastructure, both implemented in
`emerging/analysis/aberration.py` and swept over the identical
recent/baseline window every EB05 variant uses, so this is a fair swap of
statistic, not also a change of window:

- **`ears`** — EARS C2's own mechanism: the recent window's share of all OD
  deaths, read as a z-score against the sample mean and SD of that share over
  the baseline window. A classical, pre-empirical-Bayes technique.
- **`nb-trend`** — a log-linear Poisson trend fit to the window, scored by
  the fitted slope's Wald z (quasi-likelihood dispersion corrected). The
  literal "model the curve, alarm on the derivative" detector.
- **`nb-trend-own`/`nb-trend-dual`** — `nb-trend`'s own Gate-B analogue (the
  same slope fit with no county-total offset, so it reads a trend in the
  *raw count*) and the `min()` of the two, mirroring `eb05`/`eb05_own`/
  `eb05_dual` exactly. Added after checking `nb-trend` alone against its own
  off-target list — see below.

### `ears` is the cleanest model in the entire table by off-target load, at a real recall cost

**Zero off-target substances at the fixed line — the only model in this
entire file to post that at any reading, matched budget included** (5
substances there, still among the lowest). It avoids every one of the 7
labelled negatives outright, ties `nb-trend`/`nb-trend-dual` for the
highest fixed-line precision rate (0.79) of any single (non-ensemble) EB05-
family or TreeScan variant, and it does this with no gamma-Poisson shrinkage
at all — a 1954-vintage moving-average z-score run on the same window EB05
sees. That is itself informative: shrinkage is not what is buying EB05's own
precision, since `ears` gets similar or better precision from window
discipline and a classical z-test alone.

The cost is real and shows up exactly where the fixed/matched gap is
designed to expose it: fixed-line `detected` (11/22) and mean recall (0.23)
are the lowest of any non-TreeScan-leaf model in the table. `sd_b`'s sample
estimate over only `baseline_quarters` (8) points is noisy for a substance
whose baseline is thin and lumpy — an inherent cost of not pooling
information across substances the way the Gamma prior does, which is
precisely the problem GPS shrinkage was built to fix. At the matched budget
`ears` closes most of the gap (14/22, mean recall 0.39, within 0.04 of plain
`eb05`) while still running the second-cleanest off-target load (5
substances) after `nb-trend`/`nb-trend-dual` (4) — a genuinely different
point on the frontier from anything above it: buy `eb05`-comparable recall
at a matched budget, for a lower off-target load, using a method with no
Bayesian machinery in it at all.

### `nb-trend` posts the best mean recall and mean IoU of any single detector in the whole table

**At the fixed line, `nb-trend`'s mean recall (0.49) and mean IoU (0.37) are
each the highest of any non-`ratio` model here** — ahead of every EB05
variant, every TreeScan reading, and every ensemble except the two
rank-based ones already flagged as threshold-scale artifacts. Checked
substance by substance (`per_substance.csv`): `nb-trend` matches or beats
plain `eb05`'s recall on every one of the 15 detected positives, most
sharply on the ones EB05's shrinkage discounts hardest for being thin —
Xylazine (0.43 vs `eb05`'s 0.14), PCP (0.48 vs 0.17), Ketamine (0.50 vs
0.14), Clonazepam (0.38 vs 0.23), Codeine (0.44 vs 0.22) — and on Fentanyl
specifically, recall 0.97 / IoU 0.90 against `eb05`'s 0.69 / 0.69, essentially
the whole reference interval. The mechanism is not shrinkage doing something
cleverer; it is what a Wald test naturally does that a fixed EB05>1.5 line
does not — a substance with enough accumulated evidence keeps testing
significant for as long as its trend holds, rather than needing a fresh
period-over-period ratio to clear a fixed bar every quarter, which is a
better match to `recall`'s own construct (overlap with a *presence*
interval, not just the initial rise — see `ground_truth.md`'s construct-
mismatch note, also flagged in this file's own caveats section below).

**That same mechanism is what costs it on off-target load — 11 substances,
45 alarm-quarters at the fixed line, the second-worst in the table after
`treescan-branch`'s unguarded form.** Checked directly: Methamphetamine alone
accounts for 26 of those 45 alarm-quarters. This is *not* the denominator-
drift artifact `eb05_own`/Gate B was built to catch — Methamphetamine's own
raw count really was climbing for most of the sweep's history (176 in
2014Q4 to 1,802 by 2023Q1, a genuine, well-documented decade-long rise), and
`nb-trend-own` (the no-denominator, own-count-only reading) confirms it,
scoring *more* confidently than the share-based reading at every one of
those quarters (z up to 8.2, against `nb-trend`'s own z up to 4.7 there).
**This is the same category of finding `trends.py`'s `EMERGENT_THRESHOLD`
exists to name for EB05 (Fentanyl, PCP: real, credible growth in an already-
large, already-known killer, not a novel emergence)** — `nb-trend` has no
equivalent concept, and a slope test has no natural way to distinguish "this
just started happening" from "this has been climbing for a decade and is
already the county's largest cause of overdose death," so it flags both
identically. The 22-substance reference is a list of *emergences*
specifically, so Methamphetamine (and, on 4 alarm-quarters, Amphetamine —
already separately flagged as investigated-but-unlabelled in
`ground_truth.md`) reads as reviewer load here, but would not obviously read
as a modeling error to a deployment asking the plainer question "what is
still rising."

**Read against *live* deployment rather than the full historical sweep, most
of that load evaporates.** `nb-trend`'s score for Methamphetamine at every
2024–2025 as-of quarter is comfortably below the 1.96 line (0.55–1.34,
against a 2014–2023 peak of 4.7), and `nb-trend-own` has been strongly
*negative* since 2025Q1 (-3.5 to -7.4) as the county's meth deaths have
themselves been declining (1,586 → 1,345 recent-window count, 2024Q3 →
2025Q4). The 45 off-target alarm-quarters are a backtest-wide total across
45 replayed as-of dates, most of them from the real, already-over meth
epidemic's actual rise; the number a live monthly run would actually surface
today is close to zero. This is a caveat about how to read
`off_target_quarters` for a *slope* detector specifically — it is even less
of a live false-positive count than the same column already is for the
disproportionality family (see "What this benchmark cannot tell you"),
because a sustained real trend keeps testing significant for its whole
duration, not just its onset.

### The dual gate helps, but not for the reason it helps EB05

`nb-trend-dual` (`min(nb-trend, nb-trend-own)`) was built expecting to
reproduce Gate A/Gate B's win against denominator drift. Checked directly:
**it does not remove Methamphetamine from the off-target list at all** —
still 26 alarm-quarters, because both the share trend and the own-count
trend really were significantly positive through most of the meth
epidemic's actual rise, so the `min()` has nothing to suppress there. What it
*does* remove is thinner, noisier corroboration: total off-target
alarm-quarters fall 45 → 42 and off-target substances 11 → 9 (Oxycodone and
one of Hydroxyzine's two quarters drop out), mean recall gives back a little
(0.49 → 0.45) and mean IoU barely moves (0.37 → 0.35) — a small, real
precision-for-recall trade, not the clean win the EB05 dual gate posted.
**`nb-trend-own` alone is a poor detector** — matched-budget `detected` falls
to 4/22 (the worst of anything in the whole table) because a raw-count trend
with no share normalization is dominated by exactly the large, already-
established substances (Methamphetamine, Fentanyl) the dual gate exists to
downweight, and has almost nothing to say about a substance too thin to
move its own count trend on its own. It is included only as the component
that makes the dual gate legible, the same role `eb05-own-role` plays for
the EB05 family.

### Ensembling `nb-trend`/`ears` with the EB05/TreeScan family

Six designs, all following the same veto/threshold machinery validated on
`eb05-dual` + TreeScan above, none of them nesting one ensemble inside
another (`build_sweeps` does not resolve nested ensemble dependencies, so
every component here is a base detector). **One caught a real bug before it
caught anything about the data**: `combine="veto"` reads the *proposer's*
own score unchanged (`build_ensemble_sweep`'s docstring), so a veto
ensemble's fixed reading line must equal the proposer's own native line —
`ensemble-veto-nbtrend-ears-*` first shipped without that override and
silently read `nb-trend` (native line 1.96, a Wald z) at EB05's 1.5 instead,
which is a looser bar on that scale and inflated the fixed-line numbers
(recall 0.69, the highest in the whole table, off-target 40 substances, the
second-worst). Caught by the same discipline this file already applies
throughout — checking a surprising number rather than reporting it — fixed,
and now covered by a regression test
(`test_veto_ensemble_fixed_threshold_matches_its_proposers_own_line`,
`tests/test_benchmark.py`). The numbers below are post-fix.

**`ensemble-threshold-nbtrend` (`eb05-dual` OR `nb-trend`, each at its own
line) is the one design that improves on `nb-trend` alone, not just matches
it.** At the matched budget it Pareto-dominates plain `nb-trend`: identical
`detected` (11/22) and precision rate (0.92), but higher mean recall (0.31
vs 0.25), higher mean IoU (0.29 vs 0.23), and half the off-target substances
(2 vs 4). The mechanism is the same one that made `eb05-dual OR treescan`
work: a union can only add alarms, never dilute, so spending part of the
shared 85-alarm budget on `eb05-dual`'s own surgical picks instead of
`nb-trend`'s noisier ones buys precision for free. At the fixed line it
posts the best mean recall of any non-veto-corrupted single or ensemble
model in the table (0.52) at a real off-target cost (12 substances) — a
predictable result given `nb-trend` alone already dominates `eb05-dual` on
recall for every one of the 14 positives both detect (see above), so the
union's fixed-line behaviour is close to `nb-trend` alone with a little
extra noise, not two complementary coverages meeting.

**`nb-trend` and TreeScan are interchangeable as `eb05-dual`'s vetoer, which
is itself a finding.** `ensemble-veto-eb05dual-nbtrend-1` (floor z > 1) posts
matched-budget numbers (15/22, recall 0.44, IoU 0.32, 2 off-target
substances) that are statistically indistinguishable from `ensemble-veto-2`'s
(treescan vetoer: 15/22, 0.44, 0.32, 2 off-target substances) — two
structurally unrelated independent-corroboration signals (a multiplicity-
corrected tree scan vs. a Poisson trend slope's Wald test) agreeing on which
of `eb05-dual`'s own alarms deserve to survive is convergent evidence that
the veto mechanism itself is sound, not an artifact of TreeScan
specifically. The looser floor (`-0`, z > 0) is a near no-op in both
readings — identical to `eb05-dual` alone at the fixed line, meaning nothing
`eb05-dual` alarms on ever has a *negative* `nb-trend` slope, which is itself
informative about how rarely the two disagree in direction even when they
disagree on significance.

**`ears` does not clean up `nb-trend`'s Methamphetamine problem, at either
floor, and checking why is the important result here.** Floor z > 0 is a
total no-op (identical numbers to `nb-trend` alone at both readings — every
one of `nb-trend`'s alarms already has `ears_z >= 0`). Floor z > 1 removes
some real noise (off-target substances 11 → 8, Oxycodone/Chlordiazepoxide/
Quetiapine drop out) at a real recall/IoU cost (0.49 → 0.45, 0.37 → 0.34),
but **Methamphetamine keeps 22 of its 26 off-target alarm-quarters even at
the tighter floor** — checked directly (`ears_z >= 1` alongside
`nb_trend_z > 1.96`). This is not a weak-floor problem the way the `z > 0`
no-op is; `ears` and `nb-trend` are both share-based statistics, so a
substance whose *share* was genuinely, sustainedly elevated for a decade
(which Methamphetamine's was — see above) clears both bars by construction.
**No combination of these detectors resolves this, because it is not a
detection error to resolve.** Methamphetamine's decade-long rise is real by
every measure in this file; what none of `ears`/`nb-trend`/their ensembles
have is `trends.EMERGENT_THRESHOLD`'s distinction between "genuinely new"
and "already the county's largest cause of overdose death, still climbing."
Fixing this would mean building an EMERGENT_THRESHOLD equivalent for the
curve family (e.g. gating on the fitted intercept, or on `n_baseline`, the
same input EB05's own version already uses) rather than layering another
vetoer — a different, targeted piece of work, not an ensembling question.

**`nb-trend` as a second vetoer stacked on `eb05-v2-role` does not beat
TreeScan there.** `ensemble-veto-v2role-nbtrend-1` posts matched numbers
(12/22, recall 0.37, IoU 0.31, 1 fewer off-target substance than
`eb05-v2-role` alone) that are a small, real improvement on the unvetoed
baseline, but clearly behind `ensemble-veto-v2role-10` (treescan vetoer:
13/22, recall 0.39, IoU 0.32, only 1 off-target substance). `eb05-v2-role`'s
alarms are already precision-filtered by `--role-discount` and the spatial
term, and TreeScan's multiplicity correction is evidently better matched to
what's left to filter there than a trend slope is.

**Verdict.** Consistent with `docs/LITERATURE_REVIEW.md` sec 1.5's own
reading of the literature — nobody proposes replacing disproportionality
with a curve/derivative test — none of these six ensembles is proposed as a
new default. Two are worth keeping in the model table going forward:
`ensemble-threshold-nbtrend` as a genuine, if modest, Pareto improvement on
`nb-trend` alone at a matched budget, and `ensemble-veto-eb05dual-nbtrend-1`
as confirmation that `eb05-dual`'s veto mechanism generalizes to a second,
independent corroborating signal — worth having in reserve if TreeScan's own
scan ever needs to be swapped out or is unavailable, since the two currently
perform equivalently as `eb05-dual`'s vetoer. The Methamphetamine finding is
the one result here that changes the project's model of the problem rather
than just its scoreboard: it is evidence that the curve family's actual gap
against EB05 is not off-target noise fixable by ensembling, but a missing
concept (established-vs-emerging) that both a future curve-family threshold
and, in the meantime, this repository's own EMERGENT_THRESHOLD framing on
the EB05 side already exist to name.

## Closing the curve family's gap: an EMERGENT_THRESHOLD analogue, and porting #1/#3

The Methamphetamine finding above identified a real, named gap: `nb-trend`
(and `ears`) have no notion of "already-established," so a decade-long real
rise in the county's largest cause of overdose death scores identically to a
substance appearing from nothing. `trends.EMERGENT_THRESHOLD` names the
identical distinction for EB05 — but only ever as an *annotation*
(`alarm_history`'s `emergent` column), never a filter that suppresses an
alarm. Building the `nb-trend` analogue as an actual gate, rather than just
another annotation, tests directly whether that restraint was arbitrary or
load-bearing.

### `nb_trend_emergent` (`emerging.analysis.aberration._emergent_flags`)

Ported faithfully, including the one subtlety that matters:
**anchored at the current alarm episode's onset, not read fresh every
quarter** — the identical fix `alarm_history` needed, and for the identical
reason. A substance whose own rise is real and sustained would otherwise
grow its rolling baseline as the earlier part of its own rise ages into the
baseline window, and flip from emergent to not partway through the very
episode the flag exists to describe. Reading the baseline once, at first
alarm, and holding it for the episode's duration avoids that.

**It fully closes the Methamphetamine case, on the metric that matters
most.** Gated at `EMERGENT_THRESHOLD = 2.0` (reused directly from
`trends.py`, since `n_baseline` is the identical quantity computed the
identical way), Methamphetamine and Cocaine both drop to **zero**
alarm-quarters — checked directly, not inferred from the aggregate off-target
count. Off-target substances at the fixed line fall 11 → 7 (Amphetamine,
1,1-Difluoroethane, Hydroxyzine, Chlordiazepoxide, Gabapentin, Quetiapine,
Temazepam remain — all thin, none the sustained-epidemic shape Methamphetamine
was).

**The cost is real, and it is the same cost EB05 was avoiding by keeping
`EMERGENT_THRESHOLD` an annotation, not a gate.** Checked directly
(`per_substance.csv`): four labelled positives drop out entirely —
**Fentanyl, PCP, Acetyl fentanyl, Codeine.** Fentanyl and PCP are the exact
two substances `trends.py`'s own `EMERGENT_THRESHOLD` calibration note names
as "real, credible growth in an already-large, already-known killer" — this
project's reference set still counts their documented rises as labelled
emergences, so a hard gate throws those detections away along with
Methamphetamine's false one. `detected` falls 14/22 → 10/22 and mean recall
0.49 → 0.31 at the fixed line. **This is the finding that explains
`alarm_history`'s design choice rather than just replicating it**: annotating
"was this previously rare" is safe and adds real interpretive value
(Methamphetamine reads differently from Xylazine to a reviewer even while
both stay flagged), but gating on it is a real, asymmetric trade that throws
out true positives this project's own reference set still wants to see.
**Recommendation: expose `emergent` as the annotation `nb_trend_sweep`
already returns, read alongside `nb_trend_z` the way `alarm_history` reads it
alongside `eb05` — not as `nb_trend_emergent`'s hard-gated statistic**, which
is kept in the model table as the check that proved this rather than as a
candidate for adoption. **Built**: `emerging aberration rank` — the model-
the-curve family's own current-quarter report, `nb_trend`/`ears` alongside
`emergent`, mirroring `trends rank` without duplicating its EB05 machinery.
`--weighted`/`--role-discount` are off by default here (unlike `trends
rank`), since the effect is real but mixed for this family rather than a
clean win -- see the next section.

### Porting `GPS_V2_DESIGN.md` #1 (position-weighted case definition) and #3 (spatial)

**#1 ports mechanically and cleanly** — position-weighting is a change to
what counts as a death, upstream of whichever statistic runs on top of the
resulting counts, so `ears_sweep`/`nb_trend_sweep` now take the identical
`weighted`/`role_discount` flags `sweep_eb05` does and recompute
`trends.load_quarterly` fresh per as-of quarter, the identical
cross-case-aggregate discipline `sweep_eb05` already needs for the same
reason (`_ever_independent`/`_role_dispersion` must never leak later
evidence into an earlier as-of's score). Results are a real but *mixed* win,
not the clean precision gain #1 gave EB05:

- **`ears-weighted-role`** barely moves — recall/IoU both drift down by
  ~0.01-0.05 at both readings, off-target stays clean at the fixed line (0
  substances, same as unweighted) but gets *worse* at the matched budget (7
  vs 5). `ears` was already the cleanest detector in the table with nothing
  for position-weighting to fix.
- **`nb-trend-weighted-role`** posts the **best mean IoU of any model in the
  entire table (0.39)**, edging out plain `nb-trend`'s own 0.37 — a real
  gain. But fixed-line off-target substances rise 11 → **16**, not fall, the
  opposite direction from `eb05-weighted-role`'s signature win. Checked
  directly: none of plain `nb-trend`'s 11 off-target substances are removed
  by weighting (Methamphetamine and Cocaine both remain); five *new* ones
  appear (Ethanol, Hydrocodone, MDMA, Olanzapine, Venlafaxine). This makes
  sense once named rather than assumed: `--role-discount` was built and
  validated against EB05's specific failure mode (a trailing adulterant
  inflating a period-over-period *ratio*, e.g. Lidocaine), a different
  mechanism from what puts Methamphetamine on `nb-trend`'s list
  (denominator drift / already-established, nothing to do with cause-line
  position). Reweighting a count series changes its quarter-to-quarter
  *shape*, which a ratio-based statistic barely notices but a slope fit can
  react to directly — for a handful of substances that reshaping apparently
  creates a new spurious trend rather than removing an old one.
- **`nb-trend-dual-weighted-role`** is the best-behaved of the three:
  matched-budget `detected` rises to **13/22**, the best any nb-trend dual
  variant has posted (up from 11/22 unweighted), with mean IoU up
  (0.23 → 0.25) and off-target substances essentially unchanged (3 vs 4).
  Combining the two axes that each individually help (`#1` for IoU, `#2` for
  precision) does compound here, just less cleanly than it does for EB05.

**#3 was checked and found not to port at all, algebraically, before
anything was built.** EB05's spatial discount shrinks the single collapsed
`expected` a ratio is read against; `nb_trend`'s statistic is a fitted slope,
and discounting the trend's offset by the same constant factor in every
quarter of the window — the only way to reuse `concentration_weight`'s
existing per-(substance, as_of) score as-is — folds entirely into the
fitted intercept and leaves the slope, and therefore the Wald z this
detector alarms on, *provably unchanged* (`log(offset*c) = log(offset) +
log(c)`, and only the intercept absorbs an additive constant). See
`emerging/analysis/aberration.py`'s module docstring for the full argument.
A version that discounted only the recent-window quarters differently from
the baseline window's would have real effect, but that is a new mechanism to
design, not a port of #3 as it exists — not built here.

## A second, independent gap: repeated-testing inflation, found by synthetic data

The 22-substance benchmark could not have found this — it needs a
population with a *known* null (no real trend at all) at a scale no real
substance count provides. `docs/findings/synthetic.md`: `nb-trend`'s
per-quarter Wald test is calibrated in isolation, but reading it across the
~34-quarter as-of sweep with no repeated-testing correction inflates the
real false-alarm rate from a nominal ~2.5% to ~33% on synthetic data with no
trend at all — `eb05`/`ears`, checked against the identical population, show
0%. The fix, `nb-trend-confirmed` (requiring
`aberration.NB_TREND_CONFIRM_QUARTERS = 3` consecutive alarming quarters,
chosen empirically off the synthetic calibration check, not assumed), is in
the table above at both readings.

**Checked against real data, this is a second, additive cause of `nb-trend`'s
off-target load, not a restatement of the Methamphetamine finding.**
`nb-trend-confirmed` cuts the genuinely noise-driven off-target cases hard —
1,1-Difluoroethane, Amphetamine and Cocaine each fall from 3-4 alarm-quarters
to 1 — while Methamphetamine keeps 19 of its 26 alarm-quarters, because that
alarm was never noise to begin with. `nb-trend-emergent` does the reverse:
fully zeroes Methamphetamine, does nothing for a single noisy quarter
clearing the bar by chance. Two gates, each closing the failure mode it was
built for, each leaving the other's failure mode almost untouched — the
clearest evidence yet that `nb-trend`'s off-target load is not one problem,
and that treating it as one would have left half of it unaddressed either
way. The cost is real and large: mean lag roughly doubles (1.79 → 3.75
quarters at the fixed line) — a confirmation rule cannot fire before its
confirmation window elapses, by construction.

**Combined (`nb-trend-confirmed-emergent`), the two gates get about as
clean as `nb-trend` can get — 2 off-target substances, neither a documented
negative, and Methamphetamine and Cocaine both fully gone for the first
time — but the recall cost does not stay additive.** `detected` falls to
8/22, below *either* gate alone (10/22, 12/22): the confirmation gate's lag
cost lands hardest on exactly the thin, short-lived signals
(Acetyl fentanyl, Codeine, Ephedrine, Mitragynine) this benchmark keeps
finding are the most fragile category for every precision mechanism tried.
Full numbers and the substance-level detail in `docs/findings/synthetic.md`.
Read as the precision extreme of the frontier — a candidate secondary
filter alongside a higher-recall primary, the same proposer/vetoer shape
already validated for `eb05-dual` + TreeScan — not a new default.

## What the table says

**1. Shrinkage is still unambiguously worth its complexity.** The raw `n/E`
ratio has the best recall at the fixed line (0.64) and craters once its alarm
budget is matched (0.26, mean precision drops too) — at the fixed line it is
raising **872** alarm-quarters, ten times any shrinkage model, spread across
93 substances. That is not a detector, it is a list. Unchanged from every
earlier version of this table.

**2. Seven negatives, checked per model, per substance
(`per_substance.csv`) — no variant avoids more than 4 of 7, and each has a
different profile of which it gets right:**

| | eb05 | eb05-weighted | eb05-dual | eb05-spatial | eb05-spatial-prior | eb05-v2 | eb05-weighted-role | treescan |
|---|---|---|---|---|---|---|---|---|
| Morphine (blip) | alarms (1q) | **avoids** | **avoids** | alarms (1q) | alarms (1q) | **avoids** | **avoids** | alarms (1q) |
| Lidocaine (adulterant) | alarms (4q) | alarms, less (1q) | alarms (4q) | alarms (4q) | alarms (4q) | alarms, less (1q) | **avoids** | alarms, more (5q) |
| Levamisole (adulterant) | alarms (1q) | **avoids** | alarms (1q) | alarms (1q) | **avoids** | **avoids** | **avoids** | alarms (1q) |
| Alprazolam (flat) | **avoids** | alarms (1q) | **avoids** | alarms (2q) | **avoids** | alarms (2q) | alarms (1q) | **avoids** |
| Mirtazapine (blip) | alarms (1q) | **avoids** | alarms (1q) | alarms (1q) | alarms (1q) | **avoids** | **avoids** | **avoids** |
| Hydromorphone (decline) | **avoids** | **avoids** | **avoids** | alarms (1q) | **avoids** | **avoids** | **avoids** | **avoids** |
| Citalopram (noise) | alarms (2q) | alarms (1q) | alarms (1q) | alarms (2q) | alarms (2q) | alarms (1q) | alarms (1q) | **avoids** |
| **avoided / 7** | **2** | **4** | **3** | **0** | **3** | **4** | **5** | **4** |

**`eb05-weighted-role` is the only variant that fully avoids Lidocaine**, and
posts the best avoided-count in the table (5/7) doing it — the position-
weighted case definition, `GPS_V2_DESIGN.md` #1, was built specifically
*because of* Lidocaine, and `--role-discount` is what finally closes that
case rather than just reducing it (4 quarters → 1 → 0). `eb05-weighted` and
`eb05-v2` tie for second (4/7: Hydromorphone, Levamisole, Mirtazapine,
Morphine, still alarming on Lidocaine at reduced volume). `eb05-spatial` is
worst: it alarms on all seven negatives, avoiding none — direct,
substance-level confirmation of the scale-inflation finding in point 5.
`treescan` ties `eb05-weighted`/`eb05-v2` at 4/7, but with a different
profile: it's the *only* model in the table to avoid Citalopram — every
EB05 variant alarms on it, including `eb05-weighted-role` — yet it's also
the only model that alarms *harder* on Lidocaine than plain `eb05` (5q vs
4q). Its window-scan case definition has no equivalent of the position
discount, so a steady low-share adulterant is exactly the shape it's built
to catch; Citalopram, by contrast, is described elsewhere as noise
(scattered, no sustained window), which is exactly what the scan's windowed
excess test is built to *not* fire on. A genuinely different result from the
EB05-family pattern in the rest of this table, not just a weaker one.

**3. The precision gain from `--weighted` is not free, and the mean-precision
column alone would overstate it.** `eb05-weighted` and `eb05-v2` post the
lowest `detected` count in the table (12/22 and 11/22) — the substances they
miss entirely (`det_quarters = 0` where `gt_quarters > 0`) are Xylazine,
Benzylfentanyl, and Ephedrine, the three lowest-n positives in the reference
set. These are not correctly avoided false alarms; they are real emergences
the model never flagged. Both variants are **more conservative across the
board** — that's what suppresses 4 of 7 negatives and dents Lidocaine, and
it's also what costs the three thin real detections. Not a clean win, a real
trade-off, and now with a Wilson CI that says the precision gain itself
isn't even statistically distinguishable from plain EB05's at this n.

**4. `eb05-dual` is still the one variant that improves without an obvious
new cost.** It avoids 3 of 7 negatives (Morphine, Levamisole, Hydromorphone)
while still detecting 14/22 at the fixed line and 15/22 at matched budget —
only Benzylfentanyl is missed entirely among positives, and even that is a
near-miss (it does alarm, just outside the reference window: precision 0.0
there, not a total miss). At matched budget it ties for the best mean recall
(0.43) *and* posts the best mean IoU (0.32) of every model in the table. It
does nothing for Lidocaine or Citalopram — a real limitation, not a rounding
error.

**5. The spatial discount's fixed-line recall advantage is still, at least
partly, threshold inflation** — its mean recall falls at matched budget (0.49
→ 0.41), the same direction as in the 5-substance version of this finding.
Worth re-checking against `concentration.md` now that PCP and several other
substances are labelled, since the earlier version of this claim was checked
against a much smaller substance set.

**6. eb05-spatial-prior continues to decline to act** — within 0.01 of plain
EB05 on both readings, consistent with the covariate γ self-calibrating to
near-zero when spatial concentration carries no information. See
`concentration.md` for the fitted γ.

## Checking whether off-target substances are legitimate, not just unlabelled

`off_target_*` is a workload measure, deliberately not an error count (see
the definitions above) — but that's a reason to check some of them, not a
reason not to. This is the investigation that produced six of the seven
`no_emergence` cases and one new positive (Codeine); it's kept here for the
substances that *didn't* turn into labelled rows, which is as important a
result as the ones that did.

Plain `eb05`'s original off-target substances were the pool this whole
second-round investigation started from (`ground_truth.md`): four —
Lidocaine, Levamisole, Mirtazapine, Citalopram — turned into labelled
negatives, one — Codeine — turned into a labelled positive on a closer,
quarter-level look, and two — Amphetamine, 1,1-Difluoroethane — are still
unresolved either way.

**`eb05-spatial` separately surfaced 5 substances no other model alarmed
on** — Acetaminophen, Alprazolam, Hydromorphone, Hydroxyzine, Lamotrigine —
checked directly against raw quarterly counts and share. Three (Acetaminophen,
Hydroxyzine, Lamotrigine) are thin (n=22–52 lifetime), scattered, no
sustained pattern — not threats. The other two earned labelled rows for two
different reasons: Alprazolam is substantial (n=922, 3–7% share every
quarter since 2020) but flat, not rising — `eb05-spatial` alarmed on real,
non-noise data that just isn't an emergence. Hydromorphone is a real opioid
in genuine decline (share ~0.6–1.6% in 2018 down to ~0.1–0.4% since) —
Morphine's mechanism, on a substance only `eb05-spatial` was sensitive
enough to alarm on at all.

**All of this points the same direction as finding 5 above.** `eb05-spatial`
alarms on every one of the 7 negatives (0 avoided, the worst of any model —
see point 2's table) while also being the model most likely to surface a
*new*, previously-uninvestigated false positive (Hydromorphone and
Alprazolam both came from its off-target list, not any other model's).
Independent, now substance-level evidence for the same conclusion:
`eb05-spatial`'s recall advantage is bought with alarm volume, not
discrimination.

## What this benchmark cannot tell you

- **22 labelled substances is still a small denominator for anything but the
  point estimates** — the Wilson CI column exists precisely so this doesn't
  have to be estimated by eye. One quarter of overlap on a thin substance
  (Ephedrine, Benzylfentanyl) moves mean recall noticeably more than the same
  quarter would on Fentanyl. Treat point-estimate gaps under ~0.05 as
  unresolved, and remember every CI in the precision-rate table overlaps
  every other one. **This is the specific limit
  [`SYNTHETIC_BENCHMARK_DESIGN.md`](../SYNTHETIC_BENCHMARK_DESIGN.md) exists
  to lift** — known ground truth and Monte Carlo replication in place of one
  real draw per substance, so a gap like this one can be resolved rather than
  only flagged.
- **`off_target_*` is a reviewer load, not an error rate.** The reference file
  labels 22 substances now (was 5, then 15, then 18, then +4 more) — it does
  not assert every remaining unlabelled alarm is false. Amphetamine and
  1,1-Difluoroethane, still unlabelled in `ground_truth.md`, alarm for
  reasons that simply haven't been resolved either way, not confirmed errors.
- **Recall is measured against *presence* intervals**, the construct mismatch
  `ground_truth.md` documents in detail — reference intervals track whether a
  substance is still around, EB05 tracks whether it is still *rising*. This
  caps recall for every model alike; it is a fair comparison between models
  and a pessimistic absolute number for all of them.
- **`mean_lag_q` is averaged over different substances per model.** A model
  that detects fewer emergences but detects them fast posts a better mean lag
  than one that catches more slowly — `eb05-v2`'s 1.27q lag looks best in the
  table and is partly an artifact of it detecting the easiest 11.
- **Branch aggregation is now measured, but it is measured against a
  substance-level reference, which is not the unit it operates on.** The
  `treescan-branch*`/`treescan-tree` rows credit a branch's score to
  individual members so the 22-substance ground truth can score it at all.
  That is the only way to put aggregation in this table, and it imports the
  question `excess_share` cannot settle: a branch signal is evidence about
  the *family*, and apportioning it to members is an attribution model, not
  a measurement. Alprazolam is the standing proof — 34–46% of a genuinely
  rising benzodiazepine node's excess, and a labelled negative. A reference
  that labelled *families* rather than substances would score this mechanism
  on its own terms; this one cannot, so read the branch rows as a lower
  bound on what aggregation is worth and an upper bound on how cleanly it can
  be attributed.
- **The branch rows inherit the cached scan's `--reference deaths` and its
  9,999-replicate ceiling.** `branch_sweep` reshapes `results/treescan/
  backtest.csv` and recomputes attribution from the mentions table; nothing
  checks that the two were built from the same vintage. The ceiling is the
  more visible constraint: at 9,999 replicates the smallest p is 1e-4, so
  every strongly-signalling node-quarter ties at RI 10,000, which is what
  makes `treescan-branch`'s matched budget uncuttable below 234 alarms.

## Verdict

**There is no single "best" configuration — there is a recall/precision
frontier, and every result in this file is a point on it, not a winner.**
For a reader who wants one name: if recall/coverage matters most,
`eb05-dual` or `ensemble-veto-2` (plain `eb05-dual`, treescan veto at RI 2)
— highest detected counts in the whole table, close to free relative to
`eb05-dual` alone. If precision and reviewer load matter most,
`ensemble-veto-v2role-10` (`eb05-v2-role`, treescan veto at RI 10) — 6 of 7
negatives avoided and the best off-target load of anything tested, at a real
cost in detected count relative to the recall end of the frontier. `eb05-v2-
role` alone sits in between, and stacking the veto on top of it is a small
but genuine improvement *on eb05-v2-role specifically* — not evidence that
either end of the frontier is more correct than the other. Which point to
deploy at is a product decision (how much reviewer time exists to check
alarms) more than a modeling one, and every comparison in this file carries
the same caveat: at n=22 labelled substances, none of these differences have
been shown to be statistically real, only consistently ordered.

**`ensemble-veto-v2role-10` is now `emerging trends rank`/`backtest`'s
default** (`emerging/analysis/trends.py`'s `_apply_treescan_veto`, threaded
through `--weighted --role-discount --spatial discount --treescan-veto`, all
four on by default) -- picked as the precision end of the frontier over the
recall end because a live surveillance deployment pays the reviewer-load cost
of every alarm on every run, while a missed detection is one bad quarter
among the ones this file's own `mean_lag_q` column already shows most
variants catch within 1-2 quarters of each other. Runtime was checked before
committing to this as a live default, not assumed: the TreeScan veto costs
~4.5s per `rank` call (2.3s -> 6.8s total; the weighted/role/spatial pieces
together add only ~0.9s) -- real but tolerable for an interactive command run
occasionally, not in a hot loop. `backtest`'s historical veto reuses
`treescan.leaf_sweep()`'s cache rather than live-scanning all 45 quarters
(which would cost minutes, not seconds), so the full backtest is unaffected
runtime-wise (~38s, matching the pre-veto weighted sweep cost). Every
component remains independently switchable back off via its own flag, down
to plain single-gate EB05. **The vetoer stays leaf-only on purpose**:
swapping it for the whole-hierarchy reading (`ensemble-veto-v2role-tree-10`)
lets Alprazolam back through and drops the default from 6 of 7 negatives
avoided to 5, for no gain anywhere else — see the branch-aggregation section
above.

**Branch aggregation is a candidate for the *recall* end of the frontier, not
the precision end.** `treescan-tree` reaches three emergences no leaf-level
detector in this file reaches (Xylazine, Clonazepam, Codeine) and detects the
2020 designer-benzodiazepine wave a quarter earlier than any of them, but at
a matched budget it buys that breadth by spending lead time on substances it
already had, and it carries a negative (Alprazolam) that no attribution
threshold can separate from its true catches. It is not proposed as a default
here. Its clearest standalone use is the one `treescan.py`'s own `backtest`
already serves — asking which *family* is moving, where the answer is a
family and does not have to be apportioned to a member at all.

**The point estimates still favour `eb05-dual`; the confidence intervals say
this benchmark can't prove it yet.** `eb05-dual` ties for the best mean
recall and posts the best mean IoU at matched budget, avoids 3 of 7
negatives at no clear detection cost, and remains the variant with the most
coherent case for adoption. But its `precision_rate` 95% CI ([0.51, 0.88] at
the fixed line) overlaps plain EB05's ([0.53, 0.89]) almost entirely — at
n=22, none of the seven models' false-positive rates are statistically
distinguishable from each other yet. Read the point estimates as the
current best guess, not as a settled ranking; more labelled substances
(or more quarters of real data accumulating) narrows this, nothing else
does. `eb05-weighted`/`eb05-v2` avoid 4 of 7 negatives at the fixed line, but
that's bought with three real misses on the thinnest labelled substances, a
genuine trade-off that should stay a judgment call rather than a default.
**Lidocaine, the finding flagged here as the next concrete target, is now
closed at the fixed threshold**: `eb05-weighted-role` (`--role-discount`,
`GPS_V2_DESIGN.md` #1's extension) cuts it from 1 alarmed quarter to 0,
avoids 5 of 7 negatives — the best count of any variant — and posts the
highest precision in the whole table (0.86) without giving up detection
count relative to plain `eb05-weighted`. It is not a free win: at matched
budget its threshold has to drop to 1.24, the lowest in the table, and its
matched-budget recall and precision both land a little below plain
`eb05-weighted`'s. Whether that trade is worth it depends on whether the
deployment reads at a fixed line or a fixed alarm budget — worth a real
choice, not a default either way. `eb05-spatial` avoids zero of seven
negatives, the only model to do so — revisit it against `concentration.md`
now that the labelled set includes substances (PCP, several others) it
didn't before.

**`treescan` run solo is not a contender for best single detector here — 8/22
detected is the lowest in the table — but it is not measuring the same thing
as any EB05 variant, and the project's own game plan already points at
ensembling rather than picking one.** It is the cleanest model in the table
on off-target load (1 unlabelled substance ever, vs 3-7 elsewhere) and posts
the best lag on what it does catch, but it does not solve Lidocaine at all
(worse than plain `eb05`) — a concrete illustration that `--role-discount`
and TreeScan's multiplicity/window correction are answers to different
problems, not competing solutions to the same one. And this table only
exercises TreeScan's leaf-level reading; its actual selling point (branch
aggregation for substances too thin to score alone) is untested here by
construction — see the caveat above. Worth carrying forward into an ensemble
design as a detector with a distinct, complementary failure profile, not as
a ranked-and-rejected alternative to the EB05 family.

**Ensembling `eb05-dual` and TreeScan has a real answer now, and it depends
entirely on *how* the two are combined.** Three symmetric designs —
`"mean"`, `"max"` (both percentile rank) and `"threshold"` (each component's
own line, unioned) — are all a trade, not a free improvement, and the first
two actively lose to `eb05-dual` running alone (their fixed-line numbers in
the results table are a threshold-scale artifact, flagged there for exactly
that reason). The asymmetric fourth design — `eb05-dual` proposes, `treescan`
can only *veto* below a low floor — is different in kind: at
`veto_threshold = 2`, `ensemble-veto-2` **matches or fractionally beats
`eb05-dual` alone on detected, recall and IoU at both readings, while cutting
off-target substances to zero at the fixed line.** The mechanism is why: a
veto can only remove alarms, never dilute an alarm budget or drag a score
down by a sparse component's absence, which is exactly what cost the other
three designs their recall. At `veto_threshold = 10` the same design becomes
an explicit, real trade (12/22 vs 14-15/22 detected) in exchange for zero
off-target at both readings and the best matched-budget mean IoU (0.33) of
every model in the table. Net: this project now has both a validated way to
trade recall for precision that doesn't touch the case definition at all
(the veto, independent of `--role-discount`), *and*, distinct from every
other finding here, a specific combination (`veto_threshold = 2`) that
appears to be close to a free win rather than a trade — worth treating as
provisional at n=22 (see the confidence-interval caveat throughout this
file) rather than banked, but the strongest single result in this
investigation.

**`eb05-v2` is superseded by `eb05-v2-role` as the current best-of-everything
recommendation** — see the dedicated section above; Pareto-dominant at both
readings, kept as a new model id rather than a silent redefinition so
`GPS_V2_DESIGN.md`'s own `eb05-v2` checkpoint stays reproducible as written.

**The model-the-curve family (`docs/LITERATURE_REVIEW.md` sec 1.5) adds two
new points to the frontier, not a replacement for anything on it.** `ears`
is the cleanest single detector in the entire table on off-target load (0
substances at the fixed line, at a real recall cost); `nb-trend` posts the
best mean recall (0.49) and mean IoU (0.37) of any single detector here,
built from nothing more than a Poisson trend slope and a Wald test — no
gamma-Poisson shrinkage, no TreeScan multiplicity correction. Both come from
a statistically unrelated literature to everything else in this file, which
is exactly what makes them worth carrying into the same veto/ensemble
machinery validated on `eb05-dual` + TreeScan above rather than scoring them
as a fourth family competing for the same slot. **That ensembling is now
done, not just proposed**: `ensemble-threshold-nbtrend` (`eb05-dual` OR
`nb-trend`) Pareto-dominates plain `nb-trend` at the matched budget, and
`nb-trend` swapped in for TreeScan as `eb05-dual`'s vetoer
(`ensemble-veto-eb05dual-nbtrend-1`) lands statistically indistinguishable
from the TreeScan-vetoed version — two unrelated independent-corroboration
signals agreeing is itself evidence the veto mechanism is sound, not a
TreeScan-specific trick. **Methamphetamine's off-target alarm survived every
ensemble tried, but not the fix that actually targets it**: none of the
ensembles above are what closes it, because it was never an ensembling
problem — `nb_trend_emergent` (the `EMERGENT_THRESHOLD` analogue,
"Closing the curve family's gap" below) reduces Methamphetamine and Cocaine
to zero alarm-quarters directly, at the cost of also gating out Fentanyl,
PCP, Acetyl fentanyl and Codeine — the same trade `trends.alarm_history`
avoids by keeping `emergent` an annotation rather than a filter, now
demonstrated rather than assumed. `#1` (position-weighting) also ports, with
mixed results (`nb-trend-weighted-role` posts the best mean IoU in the whole
table); `#3` (spatial) does not port at all, and the reason is algebraic, not
empirical — see that section and `emerging/analysis/aberration.py`'s
docstring. See the dedicated sections above for the full results, including
the `fixed_threshold` bug the veto ensembles caught and the regression test
that now guards it.
