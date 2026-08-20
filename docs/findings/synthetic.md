# Synthetic benchmark: first findings

Design in [`SYNTHETIC_BENCHMARK_DESIGN.md`](../SYNTHETIC_BENCHMARK_DESIGN.md);
generator in `emerging/validation/synthetic.py`; tests in
`tests/test_synthetic.py`. This is the phase-1 (aggregate-count) generator's
first real output, not a finished benchmark — two scenarios, run once each,
reported here because the first one surfaced a real, previously-invisible
defect, not because the numbers below are meant to be read as final.

**What this is for, stated precisely, because it is easy to overclaim.**
This is a diagnostic tool for checking whether a detector's statistical
*properties* (calibration, power under a known-true rate) hold up, using
ground truth no real substance can supply. It is not a parameter-fitting
exercise — a threshold retuned to look optimal against this synthetic
generator's own noise assumptions would not be trusted to transfer to real
LA County data without being re-checked against the real 22-substance
benchmark, the same discipline every other change in this project already
follows. The one class of finding synthetic data *can* motivate a direct fix
for is a **generic statistical defect** — one whose correctness follows from
probability theory, not from matching this generator's specific assumptions.
The finding below is exactly that kind.

## Finding 1: `nb-trend`'s per-quarter test is calibrated; the swept sequence of quarters is not

**Setup.** 200 synthetic substances, `rises=False` (truly flat rate, no
generative trend at all), Poisson noise, flat county denominator, 45
quarters — the simplest possible true-negative population, with no
adversarial element at all.

| | per-(substance, quarter) exceedance rate | "ever alarms across its own 34-quarter sweep" rate |
|---|---|---|
| `eb05` (> 1.5) | 0.0% | 0.0% |
| `ears` (> 3.0) | 0.0% | 0.0% |
| `nb-trend` (z > 1.96) | 1.6% | **33%** |

`nb-trend`'s **per-quarter** Wald test is fine — 1.6% is close to (if
anything, slightly under) the nominal ~2.5% one-sided rate a z > 1.96 line
implies, confirming `_poisson_trend_slope`'s standard error is not
miscalibrated at the single-test level. What is not fine: `nb_trend_sweep`
scores every substance at ~34 overlapping as-of quarters, with **no
correction for the repeated look at all**, and the practical, "would this
substance ever have tripped the alarm" rate a live monthly/quarterly
deployment actually cares about is 33% — roughly 20x the nominal rate. A
naive independent-looks calculation (`1 - 0.975^34 ~= 58%`) overstates it
somewhat, since consecutive quarters' overlapping windows are correlated,
not independent tests, but the direction and scale of the inflation is real
and large.

**`eb05` and `ears` show zero inflation on the identical population.** This
is not because they correct for repeated testing explicitly — neither does,
and the real 22-substance benchmark's own literature review already flags
this as an open EB05 weakness (`LITERATURE_REVIEW.md`: "205 substances x 45
quarters, entirely unadjusted"). It is because gamma-Poisson shrinkage pools
information across the whole substance panel each quarter, which turns out
to be dramatically more conservative against pure per-substance noise than
an unpooled per-substance frequentist test, at least at this synthetic
population's scale. That asymmetry is itself new information — the real
benchmark's off-target counts could never cleanly separate "detector is
noisy" from "detector correctly saw something real about this specific
substance," because every real substance's history is entangled with
whatever else was actually happening to it.

**This directly extends the Methamphetamine finding from
`docs/findings/benchmark.md`, rather than replacing it.** That investigation
found `nb-trend`'s off-target alarms were mostly a real, sustained historical
share rise with no `EMERGENT_THRESHOLD` equivalent to name it as
"established, not emerging." This finding says a second, independent,
purely statistical mechanism is *also* inflating the off-target rate on top
of that: some fraction of `nb-trend`'s real-data alarms are not attributable
to any real signal at all, just the cost of testing the same substance at
34 correlated but non-identical windows with a threshold calibrated for one
look.

## Finding 2: denominator drift alone (no repeated-testing artifact removed) adds a second, compounding false-alarm source

**Setup.** 60 flat-rate substances (`rises=False`), Poisson noise, county
denominator declining ~0.5%/quarter (`(1-0.005)^44 ~= 0.80`, calibrated to
match the real ~18-22% decline `regime()` documents over a comparable span
— see the note in `synthetic._denom_curve`'s docstring about an earlier,
uncalibrated -3%/quarter draft that produced a 74% decline and an
unrealistically extreme scenario, caught before being reported).

| model | denominator-drift false positives |
|---|---|
| `eb05` | 0/60 |
| `eb05-dual` | 0/60 |
| `ears` | 0/60 |
| `nb-trend` | 33/60 |
| `nb-trend-own` | 27/60 |
| `nb-trend-dual` | 26/60 |
| `nb-trend-emergent` | 0/60 |

**A genuinely calibrated ~18-22% county-wide decline alone is not enough to
fool `eb05` even without Gate B** — 0/60, matching `eb05-dual`'s own 0/60
exactly. That is a real, useful negative result: it says the real
Methamphetamine case was not explainable by denominator drift of this
magnitude alone, consistent with the earlier finding that Methamphetamine's
own raw count really was elevated for most of its alarmed history, not only
its share.

**`nb-trend-own`'s 27/60 (which reads *no* denominator at all — `offset=1`)
confirms this is mostly Finding 1's repeated-testing artifact riding along,
not a second denominator-specific failure.** A trend fit with no county-total
term cannot be fooled by the county total; its false-alarm rate here (45%)
is in the same range as Finding 1's flat-Poisson "ever" rate (33%, on a
different random draw), which is the expected signature of the same
underlying cause. **The dual gate (`nb-trend-dual`, 26/60) barely improves on
`nb-trend-own` alone (27/60)** — a real, useful negative result of its own:
`min(nb_trend, nb_trend_own)` only cancels two error sources that are
largely uncorrelated, and here both components are independently inflated by
the *same* repeated-testing mechanism, riding the same count series, so
their false alarms are correlated rather than independent and the `min()`
does not buy as much as it did against the real, cleanly-separable
Methamphetamine case.

**`nb-trend-emergent`'s 0/60 here should not be over-read as a validation of
the emergent gate's precision specifically** — this run's flat profiles had
`baseline_rate` drawn from `Uniform(2, 25)`, mostly above
`EMERGENT_THRESHOLD` (2.0), so most of them would be gated out by
`baseline_rate` alone regardless of whether the gate is doing anything
subtle. A fair test of the gate's precision needs baseline rates
deliberately straddling the threshold — not yet run.

## The fix: `nb-trend-confirmed`, built and validated

`aberration._confirmed_flags`, wired in as a fourth `nb_trend_sweep` column
(`confirmed`) and a new `validation/benchmark.py` statistic
(`nb_trend_confirmed`): floor `nb_trend_z` to `-inf` unless the current
quarter is the `NB_TREND_CONFIRM_QUARTERS`-th (3) or later of a run of
*consecutive* alarming quarters. 3 was chosen empirically, not assumed: on
the same 200-substance flat-Poisson population as Finding 1, requiring 2
consecutive quarters only pulls the false-alarm rate down to ~15% (adjacent
quarters run ~0.7 lag-1 autocorrelated, so 2 in a row is a much weaker bar
than independence would suggest); 3 pulls it to ~2.5%, right at nominal.
Recall cost on the true-positive population checked alongside it: 49/50 vs
50/50 substances "ever confirmed" — negligible. `_confirmed_flags` reindexes
to the sweep's full quarter grid before rolling, specifically so a
substance's own missing rows (`n_recent == 0` quarters) can't be silently
treated as adjacent to a real gap — the same care `_episodes` already takes
for `emergent`'s own logic, applied here to a plain consecutive-run check
instead of interval arithmetic. Regression-tested directly:
`tests/test_synthetic.py::test_nb_trend_confirmed_restores_calibration_on_flat_negatives`
reproduces Finding 1's inflation and asserts the fix pulls it back down, so
a future change that silently breaks this is caught by the test suite, not
by rerunning this investigation by hand.

**Checked against the real 22-substance benchmark, not only synthetic
data**, and the result confirms the two fixes are genuinely complementary,
not two versions of the same thing:

| model | fixed-line off-target subs | fixed-line detected | mean recall | mean lag (q) |
|---|---|---|---|---|
| `nb-trend` | 11 | 14/22 | 0.49 | 1.79 |
| `nb-trend-emergent` | 7 | 10/22 | 0.31 | 1.50 |
| `nb-trend-confirmed` | **4** | 12/22 | 0.26 | **3.75** |

Checked directly which off-target substances survive each gate:
`nb-trend-confirmed` cuts the genuinely noise-driven cases hard —
1,1-Difluoroethane, Amphetamine and Cocaine each drop from 3-4 alarm-quarters
to 1 — but **Methamphetamine keeps 19 of its 26 alarm-quarters**, because
that alarm was never noise; it is a real, sustained multi-year trend that
trivially clears 3 consecutive quarters. `nb-trend-emergent` is the
opposite: it fully zeroes Methamphetamine (an established-vs-emerging
problem) but does nothing about a single noisy quarter clearing the bar by
chance (a repeated-testing problem). **Each gate closes the failure mode it
was built for and leaves the other one's failure mode almost untouched** —
concrete, real-data confirmation that `nb-trend`'s off-target load has (at
least) two additive, independent causes, not one problem with two
descriptions.

**The lag cost is the real, honest price, and it is large.** Mean lag
roughly doubles (1.79 -> 3.75 quarters) at the fixed line — a confirmation
rule cannot fire before its confirmation window has elapsed, by
construction. Whether that trade is worth it is the same kind of
deployment-specific judgment call `docs/findings/benchmark.md`'s own
recall/precision frontier discussion already makes for every other
precision-buying mechanism in this project; not resolved here.

## Combining both gates: `nb-trend-confirmed-emergent`

Built and checked against the real backtest (`nb_trend_confirmed_emergent`
in `validation/benchmark.py`, a conjunction — `confirmed & emergent`, the
same `min`-as-AND move `eb05_dual` makes for its own two gates):

| model | fixed-line off-target subs | detected | mean recall | mean lag (q) |
|---|---|---|---|---|
| `nb-trend` | 11 | 14/22 | 0.49 | 1.79 |
| `nb-trend-emergent` | 7 | 10/22 | 0.31 | 1.50 |
| `nb-trend-confirmed` | 4 | 12/22 | 0.26 | 3.75 |
| `nb-trend-confirmed-emergent` | **2** | **8/22** | **0.15** | 3.38 |

**The off-target load is now essentially as clean as it gets** — 2
alarm-quarters total, on two substances (1,1-Difluoroethane, Amphetamine)
neither confirmed nor ruled out as real by `ground_truth.md`, not a single
one of the 7 documented negatives. Checked directly: Methamphetamine and
Cocaine are both **fully** gone, not just reduced — the first nb-trend
reading in this whole investigation to clear both of them at once.

**The recall cost does not stay additive; it compounds.** `detected` falls
to 8/22, below either individual gate (10/22, 12/22) — the lowest of any
nb-trend variant tried. Checked directly which true positives are lost on
top of `nb-trend-emergent`'s own known cost (Fentanyl, PCP): **Acetyl
fentanyl, Codeine, Ephedrine and Mitragynine also drop out.** These are not
established-vs-emerging losses (they are all thin, genuinely emergent
substances) — they are the `confirmed` gate's lag cost landing on exactly
the low-n, short-lived signals this benchmark has repeatedly found to be the
most fragile category for every precision-buying mechanism tried so far
(the same substances `--weighted`'s own cost fell on in `docs/findings/
benchmark.md`). A substance whose real elevation only lasts 1-2 quarters
before the reference interval's own construct (`ground_truth.md`) considers
it "over" may never accumulate 3 consecutive alarming quarters at all.

**Verdict: this is the precision extreme of the frontier, not a new
default.** Worth having for a specific use case — a very-high-confidence
secondary filter layered alongside a higher-recall primary signal (the same
proposer/vetoer structure already validated for `eb05-dual` + TreeScan in
`docs/findings/benchmark.md`), not as `nb-trend`'s own replacement. Neither
gate should be silently assumed to compose for free; this is the concrete
demonstration that they don't.

## A generator bug the visualization itself caught: batch profiles could claim an ongoing rise they never showed

`emerging synthetic plot` (`profile_gallery.png`, `finding_close_up.png`) —
small multiples of individual synthetic substances against the true `rate(t)`
that generated them, and against the naive/fixed `nb_trend_z` alarm history
for the specific substances behind Findings 1/2 above. Built to see what a
`Profile` actually looks like, not just its ground-truth number.

The gallery immediately surfaced something Finding 1/2's aggregate tables
could not: `batch_riser` (`shape="batch", sustained=True`) and
`reverting_blip` (`shape="step", sustained=False`) were visually
indistinguishable lines — both a short spike reverting to baseline — yet
`generate_panel` labelled one a real, ongoing positive (interval running to
the series end) and the other a negative. Root cause: `_rate_curve`'s
`"batch"` branch concentrates its excess into 1-2 quarters and is flat
everywhere else *regardless of `sustained`* (a deliberate, tested property —
it models a one-off event like para-fluorofentanyl's batch arrivals, not a
rate change) — so `sustained` was a complete no-op for that shape, and two
profiles differing only in that flag produced bit-identical data with
opposite labels. Fixed with `Profile.real_positive` (`rises and sustained and
shape != "batch"`), used consistently by `generate_panel`'s ground truth and
the gallery's own badge logic; regression-tested
(`tests/test_synthetic.py::test_batch_shape_reads_as_a_negative_even_with_sustained_true`).
Did not affect any number reported above — Findings 1/2 only used flat
(`rises=False`) profiles — but would have silently mislabeled any future use
of `shape="batch", sustained=True` as a true positive. Worth naming as a
concrete instance of this project's own stated purpose for synthetic data:
a generic defect, caught by looking at the data rather than trusting the
recipe, not a parameter tuned to a specific benchmark reading.

**Second look, on review: is "batch" worth its own shape at all, if it now
always reads as a negative like a blip?** Fair challenge — a single batch
(`n_bursts=1`, the module's original and only mode) *is* just a narrower
blip, and the gallery said so honestly rather than asserting a difference
that wasn't there. What actually earns "batch" a class of its own is
recurrence: `n_bursts` (new `Profile` field) repeats the burst every
`(T - onset) // n_bursts` quarters, independent occurrences, none adjacent
— the literal reading of the para-fluorofentanyl precedent ("arrives in
batches", plural) rather than a single one-off event. `single_batch` and
`recurring_batch` are now shown side by side in the gallery on purpose: the
first still looks like a blip (that's the honest baseline), the second
visibly doesn't — several disconnected spikes, still ground-truth-negative
throughout (`Profile.real_positive`'s carve-out applies regardless of
`n_bursts` — repeated, disconnected events are still not a trend), and a
direct real-world stress test of `_confirmed_flags`: intermittent bursts
should never accumulate 3 *consecutive* alarming quarters by construction,
which is exactly the property `NB_TREND_CONFIRM_QUARTERS` is supposed to
buy. Not yet run as a dedicated finding (would need a population of
recurring-batch substances swept the way Findings 1/2's flat populations
were) — noted here as the natural next check, not run speculatively.

## A third finding from the gallery: a blip lasting exactly as long as the confirmation window can confirm anyway

Adding the actual alarm history on top of the profile-gallery figure (light
band = a quarter the trend test alarmed on; darker band = it also held for
three consecutive quarters, the confirmation rule's own bar) surfaced
something neither Finding 1 nor Finding 2 predicted: `reverting_blip` — a
temporary elevation, 3 quarters long, that reverts to baseline right after —
gets `confirmed=True` for two as-of quarters (checked directly:
`nb_trend_z` = 3.24 and 2.63, both above the 1.96 line, at the two as-of
quarters whose 4-quarter recent window happens to catch the blip's full
3-quarter plateau).

This is not the repeated-testing artifact Finding 1 describes — that
mechanism produces one-off, single-quarter false alarms on data with no
elevation at all. This is different: a real, temporary elevation whose
*duration* happens to be long enough to satisfy a rule built to distinguish
"real trend" from "one noisy quarter." The confirmation window can only ever
tell "sustained across several quarters" from "one quarter," not "sustained
indefinitely" from "sustained for about as long as the window itself" —
those two are genuinely indistinguishable from inside a fixed-width window,
by construction, not a bug in how the window is checked. `recurring_batch`
right next to it in the same figure, whose individual bursts last only 1-2
quarters, correctly never alarms at all (checked directly: zero quarters
above the alarm line across its whole sweep) — the confirmation window does
distinguish a short burst from a real trend, just not a medium one from a
real trend. Not fixed here — flagged as a real, specific limit of a fixed
consecutive-quarters rule, found by looking at where a picture disagreed
with what the label said, the same way the earlier findings on this page
were found.

## Does anything real actually look like these shapes? Checking, not assuming

`emerging synthetic match` correlates each gallery profile's sampled series
against every real LA County substance with at least 15 lifetime mentions
across at least 12 distinct quarters (60 of the 211 in the catalog clear
that bar), searching every possible time alignment and keeping the best
Pearson correlation — an honesty check on the twelve hand-picked recipes
above, not a detector benchmark: if nothing real ever resembles a given
synthetic shape, that shape is a hypothetical stress test, not a reflection
of something LA County's data has actually shown.

| profile | nearest real substance | correlation | matched window | thin? |
|---|---|---|---|---|
| rare_riser | Fentanyl | 0.88 | 2014Q1–2023Q4 (40q) | no |
| established_riser | PCP | 0.90 | 2014Q3–2024Q2 (40q) | no |
| step_riser | Ethanol | 0.80 | 2015Q1–2024Q4 (40q) | no |
| sigmoid_riser | Methamphetamine | 0.86 | 2015Q3–2025Q2 (40q) | no |
| rising_with_dip | Fentanyl | 0.82 | 2017Q2–2025Q4 (35q) | no |
| declining_with_blip | Fentanyl | 0.86 | 2022Q1–2025Q4 (16q) | no |
| rare_riser_negbin | Acetyl fentanyl | 0.81 | 2015Q4–2025Q3 (40q) | no |
| recurring_batch | Clozapine | 0.85 | 2012Q1–2015Q4 (16q) | **yes** |
| reverting_blip | Lidocaine | 0.77 | 2020Q1–2025Q4 (24q) | **yes** |
| flat_poisson | Zolpidem | 0.75 | 2012Q1–2015Q4 (16q) | **yes** |
| flat_negbin | Mirtazapine | 0.78 | 2021Q3–2025Q4 (18q) | **yes** |
| established_riser_batchy | Lidocaine | 0.84 | 2019Q4–2025Q4 (25q) | **yes** |

"Thin" (flagged automatically, `THIN_MATCH_MEAN`) means the matched window's
real mean count is under 3/quarter — Pearson r on a few single-digit counts
can read deceptively high without the shapes genuinely resembling each
other the way a high r on a well-populated window does. Read plainly: **7 of
12 profiles have a credible, well-populated real analogue; 5 do not.** The 5
are exactly the negative/noise-flavored classes (batch, blip, flat, and the
batchy positive) — LA County's real high-volume substances don't replicate
those micro-patterns convincingly at any real scale, which is itself useful
information: those shapes are legitimate adversarial stress tests precisely
*because* nothing this project has actually observed looks like them yet,
not because the recipes are unrealistic.

**Fentanyl matches three different profiles at three different windows of
its own history** (the full-series rare-to-huge rise, most of the series
with a mid-life dip, and just its last 16 quarters reading as a decline) —
not a methodology quirk, a real substance's own trajectory passing through
multiple archetypes over its lifetime. **The last of those three is worth a
specific caveat rather than a confident readout**: `declining_with_blip`'s
best match is Fentanyl's own 2022Q1–2025Q4 window, whose raw counts do fall
substantially over that span (peaking near 520/quarter in 2022Q3, down to
roughly 225–270/quarter by late 2025) — but this project's own reporting-
delay finding documents that the most recent several quarters of *any*
substance are under-reported (toxicology confirmation lags), and that
correction is deliberately withheld from the default view. Part of this
apparent decline is very likely real and part of it is very likely still-
pending toxicology; this analysis cannot tell those apart, and does not
try to. Worth checking against the reporting-delay correction directly
before treating Fentanyl's recent trajectory as a settled finding — flagged
here, not resolved.

## Recall, for context (not the headline result here)

Same synthetic population as Finding 1's true positives (30 "rare, sustained
riser" profiles, 20 "already-large, sustained riser" profiles, mixed ramp
shapes, Poisson noise, flat denominator): `nb-trend`'s mean recall (0.37-0.46
by category) is consistently higher than `eb05`'s (0.18-0.28) here too —
replicating the real 22-substance benchmark's own finding
(`docs/findings/benchmark.md`: best mean recall/IoU of any single detector)
under controlled conditions, evidence that result is a real structural
property of the slope-test approach and not an artifact of which 22 real
substances happen to exist in LA County's data.
