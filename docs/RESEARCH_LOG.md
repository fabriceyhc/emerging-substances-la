# Research log

Reverse-chronological. Each entry records what changed, what it showed, and
what it cost. Findings that later turned out to be artifacts are kept, marked,
not deleted — the artifact is usually the more useful record.

---

## 2026-08-19 — `emergent` is a second axis on `alarm_history`, not a second detector

A question about the project's framing surfaced a gap: `alarm_history.csv`
tracks every substance that ever clears `credible_rise`, and the write-up
calls a subset of those "emergent" by hand, but nothing in the table said
which. Fentanyl and PCP breach the alarm line the same as xylazine or
para-fluorofentanyl, and treating all four as one category erases the
distinction that actually matters for a reader — "growing share of deaths"
(what EB05 measures) and "was this rare before it grew" are independent
questions, and a substance already killing hundreds a quarter can have a
statistically real, alarm-worthy rise without being a new threat.

**Added `avg_baseline_deaths_q`, `latest_onset` and `emergent` to
`alarm_history()` (`emerging/analysis/trends.py`).** `n_baseline` from the
sweep (raw count over the `baseline_quarters` immediately *before* the rise)
divided by `baseline_quarters` gives a pre-rise average death rate;
`emergent = avg_baseline_deaths_q < EMERGENT_THRESHOLD` (2.0/quarter, a
chosen reading line calibrated the same way `ALARM_THRESHOLD` was — off the
backtest, not fitted). On the current sweep this cleanly separates the two
groups: Fentanyl (3.88/quarter) and PCP (3.38–10.38, below) are the only
breaching substances that fail it; every other alarm — Citalopram, Etizolam,
Clonazepam, Flualprazolam, para-Fluorofentanyl, Acetyl fentanyl, Mitragynine,
Ketamine, Bromazolam, Carfentanil — sits at 0-1.5/quarter and passes.

**Anchored on the baseline window specifically, not recent or blended,
because the obvious alternative self-cancels.** Gating on a window that
includes the rise itself would let a substance's own growth push its average
past the threshold and un-label it as emergent right when the rise is most
newsworthy — the more alarming the trajectory, the less likely the flag
would survive it.

**And anchored on the *latest* episode's onset, not the first.** A substance
can alarm, subside for years, and alarm again — `n_episodes`/`_episodes`
already track this, and the first version of this change read the baseline
off `first_breach` regardless. That is wrong for a repeat offender: fentanyl
first breached in 2016Q3 off a 2014-2016 baseline of 3.88/quarter, but if it
alarmed again today the relevant question is how rare it was *before this
rise*, not before the one nine years ago. Fixed by reading `n_baseline` at
`latest_onset = eps[-1][0]` (the start of `_episodes`' last run) instead.
Single-episode substances are unaffected — `latest_onset == first_breach`
whenever `n_episodes == 1`, which is every current row but PCP. PCP's second
episode (2021Q4) reads baseline 10.38/quarter, off a window that already
includes its own first episode — a stronger "not emergent" read than the
first-breach version's 3.38, for the same underlying reason `_role_dispersion`
and the position-weighted case definition keep re-learning: a substance's own
recent history is never a neutral yardstick for it once it has already been
elevated.

**This is deliberately not a new detection line.** `emergent` never fires on
its own — only annotates rows `credible_rise`/`eb05 > threshold` already
selected — for the same reason `alarms`' own docstring keeps it a distinct
diagnostic from `rank`/`backtest`'s dual gate: conflating "is this real" with
"is this new" would change what the table measures, not just add a column to
it. A hard count cutoff will still flicker at its own boundary the way raw
ratios did before EB05 shrinkage was added for the growth question — a
substance sitting at 1.9 vs 2.1 baseline deaths/quarter isn't epidemiologically
different — so `avg_baseline_deaths_q` is kept alongside the boolean rather
than only the collapsed flag, and `EMERGENT_THRESHOLD` is documented as a
convention, not a significance test, exactly like `ALARM_THRESHOLD`.

**Surfaced on `rising_candidates.png` too, not just the CSV.** `_panel`
(`trends.py`) takes an optional `badge` and draws it as a small filled
rounded-rect in the panel's top-right corner — an orange `EMERGENT` tag,
CSS-chip style, rather than folding it into the title text where it would
have to compete with `_bold_title`'s existing mathtext escaping. Panel
ordering stays plain EB05+ (unchanged) rather than grouping emergent
candidates first — sorting by tag was tried and reverted, since the ranking
axis and the annotation axis are supposed to stay visibly separate, the same
point the CSV-side design makes. The caption line was long enough to need a
manual break before "Yellow band =" — `fig.text` does not reflow on its own.

---

## 2026-08-19 — Reading TreeScan at its branches: the aggregation the benchmark had never exercised

`benchmark.py`'s `treescan` row was always the substance-**leaf** node's
recurrence interval — the one node level a tree scan has no advantage at, and
the only one the substance-level ground truth can score directly. The entry
below (2026-08-18) flagged that explicitly: *"This comparison only exercises
TreeScan's substance-leaf reading; its actual case (branch-level aggregation
for substances too thin to score alone) is untested by construction."* This
closes that.

**New: `treescan.branch_sweep`**, which propagates each internal node's score
back down to its member substances onto the same `(substance, as_of)` grid
every other detector is scored on, under two guards — *presence* (the
substance has a case inside the window that fired; the vectorized form of
`backtest`'s `as_of >= first_seen`, without which `root` alone hands its score
to all 211 substances) and *attribution* (`excess_share` ≥ 15%, the band
`backtest` already calls "partly attributable"). Attribution is recomputed in
one array pass per as-of quarter rather than per (node, substance) pair, and
`tests/test_treescan.py` pins it against `excess_share`'s own value so the
two cannot drift. Four new benchmark models: `treescan-branch` (unguarded),
`treescan-branch-attr`, `treescan-tree` (`max(leaf, attributed branch)` —
TreeScan as actually specified), and three ensembles swapping the tree
reading in as the vetoer. No new Monte Carlo work: all of it reshapes the
cached `treescan backtest` table.

**Aggregation reaches three labelled emergences no leaf-level detector in
the file reaches** — Xylazine via `Cutting agents and adulterants` (21% of
that family's excess, RI 10,000; the same substance the README names as
beyond reach of the spatial term too, at 9 lifetime mentions), Clonazepam via
`Benzodiazepines` (16–22%), Codeine via `Prescription opioids` (18%). The
cleaner demonstration is the **2020 designer-benzodiazepine wave**: Etizolam
and Flualprazolam both detect at **lag 0q** through `Designer
benzodiazepines` against 1q at the leaf, at 59–64% and 27–46% attribution —
the family node sees them a quarter before either leaf does and the
attribution says both are genuinely driving it.

**And the recall gain does not survive a matched budget.** `treescan-tree`
vs. leaf-only: `detected` 8 → 11 fixed, 8 → 10 matched, mean IoU 0.27 → 0.34
/ 0.28 → 0.30, matched precision rate 0.67 → 0.77 — but **mean recall at a
matched budget is identical, 0.36 both ways**, with the threshold climbing
from RI 60 to RI 3,333 to pay for the extra substances. Aggregation spreads a
fixed alarm budget across *more* emergences rather than covering the ones it
already had more completely; Bromazolam (1q → 2q) and Flualprazolam (0q → 1q)
each give back a quarter of lead time to fund Xylazine and Clonazepam. Same
lesson as this log's original tree-scan head-to-head, arrived at from the
other direction: the fixed line flatters, the budget tells.

**The attribution guard is load-bearing, and it still cannot save
Alprazolam.** Unguarded branch propagation raises 389 alarm-quarters over 56
off-target substances (mean precision 0.40) and **cannot be read at a matched
budget at all** — 234 of its rows tie at RI 10,000, the Monte Carlo ceiling
at 9,999 replicates, because one node firing at p = 1e-4 stamps that score
onto every member present. The guard removes 284 of those alarms and 49 of
those substances for 2 of 12 detections. What it cannot remove: Alprazolam, a
labelled negative, riding `Benzodiazepines`/`Depressants` at 34–46% excess
share — *above* Clonazepam's 16–22% and Xylazine's 19–21%. `excess_share`
measures contribution, not causation, and a large stable member of a rising
family will always look attributable. There is no threshold separating the
true catch from the false one; `treescan-tree` avoids 3 of 7 negatives
against the leaf reading's 4.

**Negative result: a less sparse vetoer does not improve the veto.** The
hypothesis was that `ensemble-veto-*`'s limitation is `treescan`'s sparsity
(95 of 211 substances ever scored; anything unscored is floored to RI 1 and
vetoed unconditionally), and that the branch reading's 162 substances would
fix it. At RI 2 the two vetoers are byte-identical (the threshold is nearly a
no-op either way); at RI 10 the tree vetoer is consistently *slightly worse*,
and on the deployed default it drops 6-of-7 negatives avoided to 5 by letting
Alprazolam through. A vetoer that scores more substances vetoes fewer alarms,
and the alarms it lets through are substances riding large families —
precisely the population a precision-oriented veto wants to keep blocking.
`ensemble-veto-v2role-10` keeps the leaf-only vetoer; that is now checked
rather than inherited.

**Unchanged: the nitazenes still never signal.** `Novel synthetic opioids`
peaks at RI 1.16 across all 45 quarters and the `Nitazenes` node never fires.
That family is `core/tree.py`'s stated reason for building a family layer at
all, and aggregation does not rescue it — three mentions pooled is still
three mentions. Branch aggregation helps families whose members already carry
real counts; it does not manufacture power where there is none. The 2026-08-10
entry's "**Aggregation power — failed**" verdict stands as written for that
case, and is now bounded rather than general: it failed for the family it was
built for, and works for benzodiazepines, adulterants and prescription
opioids.

**What this still cannot settle.** Aggregation operates on families, and the
22-substance reference labels substances, so scoring it at all requires
apportioning a family's evidence to members — an attribution model, not a
measurement. Alprazolam is the standing proof. Read the branch rows as a
lower bound on what aggregation is worth and an upper bound on how cleanly it
can be attributed. Full tables, the per-branch attribution and the vetoer
comparison: `docs/findings/benchmark.md`.

---

## 2026-08-18 — Wired `ensemble-veto-v2role-10` in as `rank`/`backtest`'s default

Moved the strongest result from the comparison framework into the live
pipeline: `emerging trends rank` and `backtest` now default to `--weighted
--role-discount --spatial discount --treescan-veto` (all four on), i.e.
`eb05-v2-role` with TreeScan able to veto (never add) an alarm it sees no
corroborating rise for at all.

**New: `_apply_treescan_veto` (`emerging/analysis/trends.py`)**, plus two RI
sources -- `_treescan_ri_live` (one quarter, ~4.5s, a local import to avoid
the circular dependency with `treescan.py`, which already imports from
`trends.py`) for `rank`, and `_treescan_ri_history` (the cached
`treescan.leaf_sweep()`, effectively free) for `backtest`'s 45-quarter sweep.
Deliberately not derived from `benchmark.py`'s `build_ensemble_sweep` --
production code depending on the comparison/validation module would be a
layering inversion, so the veto logic is small enough to duplicate cleanly
rather than share.

**Runtime was checked, not assumed, before committing to this as a live
default** (a concern raised directly): `rank` goes from 2.3s to 6.8s, almost
entirely the TreeScan veto's Monte Carlo scan (the weighted/role/spatial
pieces add ~0.9s combined). Real but tolerable for an interactive command.
The historical path in `backtest` avoids the expensive route entirely --
using the cache instead of a live scan per as-of quarter keeps the full
45-quarter run at ~38s, matching the pre-veto weighted-sweep cost, instead of
the minutes a live per-quarter scan would cost.

**Deliberately left `alarms`/`alarm_history` on plain single-gate `eb05`,
not touched.** Its own docstring frames it as a distinct diagnostic --
raw-`eb05`'s breach history, for false-alarm-rate and clustering analysis --
not "the alarm system"; rewiring its comparison basis to the dual gate plus
veto would change what it measures, not just which config it uses, so that
was scoped out rather than done partially. Every piece of the new default
remains independently switchable: `--no-weighted --no-role-discount
--spatial none --no-treescan-veto` reproduces plain EB05 exactly.

---

## 2026-08-18 — Is `eb05-v2-role` "the best configuration"? No -- there's a frontier, not a winner

Direct follow-on: does stacking the TreeScan veto (the one ensemble design
that helped, above) on top of `eb05-v2-role` — the best *solo* GPS
detector — improve on `eb05-v2-role` itself, or was the earlier veto result
specific to plain `eb05-dual`?

**Genuine, close-to-unambiguous improvement at `veto_threshold = 10`.**
`ensemble-veto-v2role-10` at the matched budget: 13/22 detected (up from
`eb05-v2-role` alone's 12), mean recall 0.39 (up from 0.37), mean IoU 0.32
(up from 0.31), precision rate 0.81 (up from 0.80), off-target down to 1
substance (from 4) — every column ties or improves. At the fixed line it
avoids **6 of 7 labelled negatives, the best count anywhere in this whole
investigation** (only Citalopram still alarms), Lidocaine included this time
— unlike the plain-`eb05-dual` veto, which never touches Lidocaine because
treescan agrees it's really elevated, `eb05-v2-role` already suppresses
Lidocaine upstream via `--role-discount`, so the veto has nothing left to
catch there.

**Not proof that stacking is generally correct, though — it's still on the
same frontier, just further along it.** Against `ensemble-veto-10` (built on
plain `eb05-dual` instead), the `eb05-v2-role`-based version trades away real
detections for the extra precision (13/22 vs 14/22, 0.39 vs 0.42 recall).
Every precision-favouring mechanism compounded so far (role discount +
spatial + dual + veto) has kept paying for itself in this backtest, but each
layer still buys precision by spending detections, exactly like every other
result in this file. Landed on stating this explicitly in the Verdict rather
than naming a single winner: `docs/findings/benchmark.md` now names two
points on the frontier (`eb05-dual`/`ensemble-veto-2` for recall,
`ensemble-veto-v2role-10` for precision) instead of one champion.

---

## 2026-08-18 — Checked two intuitions about the dual gate: drop spatial, or drop Gate A

Two follow-on questions about `eb05-v2-role`, both checked against the
backtest rather than argued from priors.

**Is `eb05-v2-role` actually an ensemble, or a genuinely integrated
detector?** The latter — `--weighted`/`--role-discount` reweight the
substance's case counts once, in `load_quarterly`, before either GPS gate is
fit; `eb05` and `eb05_own` are then both computed from that same reweighted
`counts` table and combined at the end via `min()`. Unlike the TreeScan
combinations, there's one case definition run through GPS twice, not two
independently-scored detectors glued together post hoc.

**Is #3 (spatial) safe to drop as the more speculative mechanism, once #1
extended + #2 are already in?** Checked (`eb05-dual-role`, no spatial, vs
`eb05-v2-role`): not a clean win either way. Dropping spatial costs the
largest single-column gap in the whole benchmark (fixed-line recall
0.34→0.26, IoU 0.30→0.23) but gains one extra detected substance at the
matched budget (13/22 vs 12/22). The numbers don't support the intuition
that spatial is inert — it's still the largest measured effect in this
comparison — though at n=22 the precision-rate CIs overlap the same way
every pair in this file does, so it isn't a proven effect either.

**Is Gate A (share) still doing anything once role/position weighting is
already in the case definition, or would the self-referential gate alone
(`eb05_own`, weighted + role-discounted) be a leaner equivalent?** Tested as
`eb05-own-role`: clearly worse, not leaner — precision rate falls to 0.73
from 0.85 and off-target load explodes to 10 substances (worst in the
section), because a substance can have a genuinely elevated absolute count
against its own history without that being surveillance-relevant, which is
exactly the check Gate A's share comparison provides and a self-referential
gate structurally cannot. The dual gate's AND is non-redundant work, not
belt-and-suspenders. Full tables: `docs/findings/benchmark.md`.

---

## 2026-08-18 — A fourth ensemble design: EB05 proposes, TreeScan vetoes

Follow-on to the entry below, prompted by a specific proposal: instead of
combining `eb05-dual` and `treescan` into one score (`"mean"`, `"max"`,
`"threshold"`, all below), let `eb05-dual` propose its own alarms unchanged
and have `treescan` veto only the ones it sees no corroborating rise for at
all, below a low floor well under `treescan`'s own RI 100 alarm line.
Implemented as `combine="veto"` in `build_ensemble_sweep` — the score stays
in the proposer's own native units except where the vetoer's raw score falls
below `veto_threshold`, where it is floored, so it's the one combination
where the *fixed* reading is not a scale artifact (see the "mean"/"max"
caveat below).

**This is the first ensemble design that matches or beats `eb05-dual` alone
at both readings, rather than trading against it.** At `veto_threshold = 2`
(RI ≥ 2 — a very low bar): detected, mean recall and mean IoU are tied or
fractionally better than `eb05-dual` solo at both the fixed line and the
matched budget, while off-target substances drop from 1 → 0 (fixed) and
3 → 2 (matched). The mechanism explains why the earlier three designs
couldn't do this: a veto only ever *removes* an alarm the vetoer flatly
disagrees with, so it can't dilute the proposer's own best picks the way a
shared alarm budget (`"threshold"`) or a rank average (`"mean"`, silently
punished for the vetoer's sparse coverage) both do. At `veto_threshold = 10`
the trade becomes explicit rather than free — 12/22 vs 14-15/22 detected —
in exchange for zero off-target substances at both readings and the best
matched-budget mean IoU (0.33) in the entire table. Neither veto variant
touches Lidocaine (`treescan`'s own solo signal there clears even RI ≥ 10
easily, so the veto correctly does not fire) — that stays `--role-discount`'s
problem to solve, confirming the two mechanisms are complementary rather
than redundant. Full numbers and the worked example:
`docs/findings/benchmark.md`.

---

## 2026-08-18 — `eb05-v2-role`, and three ways to ensemble EB05 with TreeScan

Follows directly from the TreeScan integration below: once TreeScan was in
the same comparison as the EB05 family, the natural next questions were
whether `eb05-v2` should absorb `--role-discount`, and whether combining
`eb05-dual` with `treescan` beats either running alone.

**`eb05-v2-role` (`weighted + role_discount + dual-gated + spatial`) Pareto-
dominates `eb05-v2`** at both readings — equal-or-better on every column
(precision 0.79 → 0.85 fixed, 0.75 → 0.80 matched; mean IoU up at both) with
no real cost (the one column that moves the wrong way, mean lag, is averaged
over only 11 substances and too thin to read as real). `eb05-v2` predates
`--role-discount` and was never rebuilt with it. Added as a new model id
rather than redefining `eb05-v2` in place, so `GPS_V2_DESIGN.md`'s own
checkpoint stays reproducible as documented; `eb05-v2-role` is now the
recommendation whenever that combination is wanted.

**Three ensemble designs for `eb05-dual + treescan`, in `build_ensemble_sweep`:**
percentile rank averaged (`"mean"`), percentile rank maxed (`"max"`), and each
component's own already-calibrated threshold maxed (`"threshold"` — a literal
alarm union). The first two were the obvious approach and both **lose to
`eb05-dual` running alone** at a matched alarm budget (9/22 and 11/22 detected
vs. 15/22) — diagnosed, not just observed: ranking a component's score among
only the substances *that component scored that quarter* means TreeScan's
"top-of-3-scored" on a quiet quarter and EB05-dual's "top-of-80-scored" on a
busy one land on the same 0–1 scale, so a sparse component's opinion counts
for as much as a dense one's regardless of how much evidence backed it up.
Their fixed-line numbers (15/22 detected, 0.94–0.98 recall) are worthless —
mean/max of `k` percentile ranks clusters at or above 0.5 by construction,
independent of any real signal, which is documented prominently in
`docs/findings/benchmark.md` next to the table that would otherwise look like
a win. `"threshold"` (each component read at its own line, unioned) is the
one design that adds real value: best precision rate (0.80) and cleanest
off-target load (1 substance) of every model in the entire comparison, EB05
family included — but it is a precision trade, not a free win, since it still
doesn't beat `eb05-dual` alone on matched-budget recall (0.38 vs 0.43).
Adding `ratio` as a third component made every reading worse (off-target
1 → 7 substances at matched budget, no recall gain), which is the "floor, not
a candidate" framing `Model`'s docstring already gives it, now checked rather
than assumed. Full tables: `docs/findings/benchmark.md`.

---

## 2026-08-18 — Ground truth to 22 substances, a role-discount extension to #1, TreeScan added to the benchmark

Three follow-ons to the head-to-head below, each triggered by a question about
the benchmark itself rather than a new detector.

**Ground truth grew from 5 to 22 labelled substances**
(`data/raw/ground_truth/known_emergences.csv`), sourced from
`alarm_history.csv`'s unlabelled breaching substances — every one that is ever
a *primary cause of death*, deliberately excluding inert adulterants like
lidocaine that never are, to keep the label LA-death-count-based rather than
drug-supply-presence-based (presence necessarily precedes death by an unknown
lag; scoring against presence would validate a different, laggily-related
construct). 15 positive, 7 `no_emergence` negatives — the first version of
this file had none, so no model's precision could ever be told apart from
another's. Added Wilson score CIs (`ground_truth.wilson_ci`,
`precision_rate`) over normal-approximation, since n is small per substance;
substance-level (not quarter-pooled) to avoid treating one substance's
correlated alarm-quarters as independent Bernoulli trials.

**`--role-discount` extends GPS v2 #1** (`trends._role_dispersion`): a
substance-level `phi_s`, computed from that substance's own aggregate mean
credit across every mention (not just each death's own position), stacked on
top of `--weighted`'s per-death discount. Answers a question `--weighted`
structurally cannot: lidocaine's typical cause line is a plain two-item
`[Fentanyl, Lidocaine]`, which keeps full per-death credit under position-only
weighting no matter how often it happens. First version used
`SOLO_CREDIT / mean(credit)` as the discount and broke `backtest` across the
board (Bromazolam's peak EB05 2.69 → 1.00) — caught by running the actual
backtest, not by inspection. Root cause: `_ever_independent`-exempted
substances still rarely average near `SOLO_CREDIT = 2`, since solo mentions
are rare for almost any real drug (checked directly: Bromazolam mean 1.000,
Fentanyl 1.223, PCP 1.001). Fixed by rebasing on `NAMED_CREDIT = 1` with a
`max(1.0, ...)` floor so no substance gets *amplified*. Re-validated: known
emergences barely move, lidocaine goes from 1 alarmed quarter under
`eb05-weighted` to 0.

**TreeScan added to `emerging/validation/benchmark.py` as a ninth model**,
run solo (`Model(statistic="treescan")`, reading its own leaf node's
recurrence interval, no EB05 underneath it) rather than only compared
head-to-head the way `treescan.py`'s own `backtest` command already does it.
Reuses that command's cached Monte Carlo output (`treescan.leaf_sweep`) rather
than re-running 45 quarters × 9999 replicates inside the benchmark. Needed a
`fixed_threshold` override on `Model`, since TreeScan's native scale
(recurrence interval, ~1-10,000) has nothing in common with EB05's (a
posterior percentile, ~0-3) — the shared `--threshold` would be meaningless
applied to it. Result: lowest `detected` in the table (8/22) but cleanest by
far on off-target load (1 unlabelled substance ever alarmed, vs 3-7 for every
EB05 variant) and best lag on what it does catch — a different failure shape,
not a strictly worse detector, and it does *not* solve lidocaine (5 alarmed
quarters, worse than plain EB05's 4) since nothing in its case definition
discounts cause-line position. This comparison only exercises TreeScan's
substance-leaf reading; its actual case (branch-level aggregation for
substances too thin to score alone) is untested by construction — see
`docs/findings/benchmark.md`'s caveats section. *(Closed 2026-08-19 by
`treescan.branch_sweep` and the `treescan-branch*`/`treescan-tree` models —
see the entry above.)* Full tables and the
per-negative breakdown: `docs/findings/benchmark.md`; `docs/findings/
ground_truth.md` for the 22-substance reference and its definitions.

---

## 2026-08-18 — GPS v2 §3 (spatial concentration), and a head-to-head of every detector

Closes the last open component of `docs/GPS_V2_DESIGN.md` and then, because
three extensions each validated separately against plain EB05 still do not say
which detector to run, scores all seven side by side
(`emerging/validation/benchmark.py`, `docs/findings/benchmark.md`).

**The `z_s` §3 asked for did not exist, and the obvious source could not be
made to work.** The design note is right that `geo.py`'s lifetime localization
ratio is the wrong input — PCP's Skid Row cluster is real and a decade old,
and letting stable historical geography lower the bar for an unrelated event
is the failure this whole component is meant to avoid — and it points at
`spacetime.py`'s trailing cylinders instead. That turns out to be unusable
for the purpose on two counts: 999 permutations per substance cannot be run 45
times over an as-of sweep, and `MIN_CASES = 20` leaves it undefined for four
of the five labelled emergences.

**What shipped: exact null moments instead of a permutation test**
(`emerging/core/concentration.py`). A Cuzick–Edwards-style case-control
clustering count — ordered pairs of a substance's deaths that are spatial
neighbours, over the symmetrized k-NN graph on all overdose deaths in the
trailing window, k = 5% of the window. Under the same random-labelling null
`geo.py` defends, that statistic has closed-form mean and variance, so
`z_s` is one sparse mat-vec per substance and the whole 45-quarter ×
211-substance sweep costs **4.6 seconds**. The derivation is the load-bearing
claim, so it is checked against 4,000 draws from the null it describes at
C ∈ {5, 10, 30, 100} rather than against intuition — mean and sd match
throughout.

*Two designs measured and rejected first, both on real data:* zipcode cells
(background co-location probability 0.009, so a 6-death substance scores S = 0
three quarters of the time — no resolution where it is needed), and a fixed
radius (valid, but its power goes to "is this substance downtown" rather than
to concentration, because degree tracks density). k-NN at a fixed *fraction*
is density-adaptive and scale-free the way `spacetime.MAX_SPATIAL_FRACTION`
is.

*The floor is not decorative.* At C = 3 the statistic can only be 0, 2, 4 or
6, so one coincidental pair is four standard deviations — diphenhydramine
scored z = 6.7 on three deaths before `MIN_CASES = 10` was added. The
simulated tail rate at C = 5 is 10.6% against a nominal 5%, and at C = 10 it
is 6.1%; that table is what sets the floor.

**Positive result: it tracks onset, and the decay is the interesting part.**
para-Fluorofentanyl scores z = 3.18 in 2021Q2, the quarter after its first LA
death, and decays to −1.83 by 2022Q3 as it saturates the county-wide fentanyl
supply. Mitragynine scores 4.63 at its onset, carfentanil 3.47 at its 2024
cluster. A substance appearing in a few places and then diffusing is the shape
a point-source hypothesis predicts, and it is now a number rather than an
impression from a map. PCP (3.84) and cocaine (5.19) corroborate `geo.py` and
`relrisk.py` from a different statistic on a different window.

**Negative result: the E-discount fallback fails, and it fails in the way the
design note warned about for a different mechanism.** Read at the fixed 1.5
line it posts the best recall in the benchmark (0.44 vs plain EB05's 0.38).
Re-cut so it raises the *same 85 substance-quarter alarms*, it drops to 0.40 —
and the substance it appeared to gain, xylazine, has `z_s` identically zero in
all 45 quarters, so the discount cannot have moved xylazine at all. It moved
the threshold, by inflating every scored substance's EB05 at once. §3 rejects
`EBGM × spatial_RR` for silently discarding the conservative-bound property;
the discount arrives at the same place through the back door.

**The covariate-adjusted prior is the one to keep, because it can decline to
act.** `beta_s = beta_0 · exp(−γ · z_s)` with γ estimated by the same MLE:
fitted γ at the current quarter is **+0.061**, and the model lands within 0.01
of plain EB05 on both readings. That is the property `w = 0.25` structurally
cannot have — a chosen weight applies itself whatever the data say. The design
note called this "less principled"; the benchmark shows it is not an
aesthetic objection but the mechanism of the fallback's failure.

**The benchmark, and why it is read twice.** Seven models — raw `n/E`, classic
EB05, and the five variants — scored against the ground-truth intervals at a
fixed line and at a matched alarm budget. The matched reading exists because
this repository already learned the lesson once: the tree-scan head-to-head
(2026-08-10) found an apparent earlier-detection advantage that disappeared
once its budget was matched to EB05's.

  - *Shrinkage is the one thing here that unambiguously earns its complexity.*
    Raw `n/E` looks best at the fixed line (recall 0.54, all five found) and
    is worst by a wide margin at matched budget (0.24, one emergence lost, 85
    alarms spread over 37 substances against EB05's 19). At the fixed line it
    raises **872** alarm-quarters — ten times every other model. That is not a
    detector, it is a list. `trends.py` has asserted this since it was
    written; it is now measured.
  - *None of the three GPS v2 extensions improves on plain EB05.* At matched
    budget the six shrinkage models span 0.35–0.40 mean recall. Five labelled
    substances contribute 57 reference quarters, so that whole spread is about
    four quarters of overlap — inside what a couple of differently-timed
    deaths would move. The honest reading is "unresolved", not "eb05-dual
    wins".
  - *Where the spatial term does act it acts correctly, and the limit is power
    not modelling.* Its matched-budget gains are on para-fluorofentanyl
    (recall 0.40 → 0.45, detection lag 1 quarter → 0) and mitragynine (0.13 →
    0.20) — exactly the two substances that clear `MIN_CASES` through their
    emergence, in exactly the quarters where z is 3.18 and 4.63.
  - *`eb05-v2` trades a substance for a quarter of lead time.* Best detection
    lag in the table (1.50 vs 2.20 quarters) and one of two models that miss
    an emergence outright — xylazine, in both readings. `GPS_V2_DESIGN.md` §1
    predicted this ("xylazine never clears the §2 dual gate under
    `--weighted`"); the benchmark confirms it end to end. Two independent 1.5
    crossings at n = 7 is a materially harder bar than one.

**Decision: plain EB05 stays the default.** `--weighted`, `--spatial` and the
dual gate remain available as what they always were — diagnostics answering
specific questions (adulterant or driver, denominator drift, localized or
diffuse) — not as a better ranking.

*The binding constraint is now the reference set, not the statistic.* Five
labelled substances cannot resolve a 0.05 difference in mean recall, and most
of the table sits inside that. A better detector cannot be *demonstrated* on
this many labels however good it is, which makes growing
`known_emergences.csv` the highest-value next move for this line of work.

*Recorded so it is not re-litigated:* `off_target_quarters` in the benchmark
counts alarms on unlabelled substances. It is a reviewer load, not a
false-positive rate — lidocaine and PCP both alarm for documented reasons and
neither is in the reference file. It is reported because a detector firing on
109 substances costs something real, not because those 109 are wrong.

---

## 2026-08-18 — GPS v2: position weighting, a dual gate, and a ground-truth reference for `backtest`

Prompted by a literature-review discussion of the "keep GPS, triangulate
around it" verdict: three formalizations of the corroborating checks
(`polysubstance.py`'s position analysis, `regime()`'s denominator-drift
diagnostic, spatial concentration) into the ranking statistic itself, rather
than downstream cross-references a reader has to reconcile by hand. Full
design record, including two rejected variants and open questions, in
`docs/GPS_V2_DESIGN.md`. Only the first two of the three are built; the third
(spatial concentration in the GPS prior) is still a proposal.

**1. Position-weighted case definition — failed its own validation gate, was
fixed, and the fix taught a general lesson.** `trends rank --weighted` /
`trends backtest --weighted` replace the raw mention count with an integer
credit from cause-line position (`_credit_from_order`): full credit for a
solo mention, discounted to zero for a substance named last on a 3+-substance
line — the adulterant pattern documented for lidocaine (last on 82% of its
lines, first on none). Two problems surfaced in that order, both caught by
running the thing rather than reasoning about it:

  - *The scheme was overdispersed, not merely integer-valued.* The design note
    first argued integer credit was enough to keep `_nb_logpmf`'s Poisson
    assumption intact. It is not — dispersion, not integer-vs-float, is the
    property that matters, and any credit spread across `{0, 1, 2}` has
    `E[c²] > E[c]`, i.e. more variance than a same-mean Poisson. On this
    repository's data it inflated the false-`EB05>1` rate from 5.1/203
    substances (raw counts) to 19.0/203. Fixed with a quasi-likelihood
    dispersion correction (`_credit_dispersion`, φ = E[c²]/E[c] ≈ 1.42
    measured directly off the credit distribution): divide both `n` and every
    `expected` array by φ before fitting — restores the false-alarm rate to
    4.9/203, back in line with raw counts.
  - *Even dispersion-corrected, it regressed a known emergence.* `backtest
    --weighted` — the validation the design note specified before trusting
    this — showed xylazine's peak EB05 dropping from 1.93 to 1.05 and its best
    rank from 4 to 12, never reaching top-10. Cause: xylazine is a fentanyl
    adulterant by *pattern* (habitually trailing, same surface shape as
    lidocaine) but not by *pharmacology* — it is exactly the kind of emerging
    contaminant this pipeline exists to flag, not an inert cutting agent.
    Cause-line position measures "was this compound clinically dominant in
    this death," a different construct from "should surveillance care about
    it," and the credit scheme conflated them. Fixed with `_ever_independent`:
    a substance that has ever been named alone or first *anywhere in its
    history* is exempted from the trailing discount everywhere, even in cases
    where it is itself trailing. Clean separation on real data — lidocaine and
    levamisole (known inert adulterants) are solo-or-lead **zero** times
    across 42 and 12 lifetime mentions; xylazine clears it once, the quarter
    before its emergence is detected. Xylazine's peak EB05 recovers to 1.46,
    best rank to 5, top-10 restored at 2024Q3 (one quarter behind unweighted).

  `_ever_independent` is a cross-case aggregate, unlike per-death credit
  (which only ever reads that one death's own text and so cannot leak
  future information forward). Computing it once from the full history and
  reusing it for every earlier `backtest` quarter would let a substance's
  *later* proof of independence retroactively un-discount its *earlier*
  trailing mentions. Fixed by bounding `load_quarterly`'s credit computation
  to its own `cutoff` (previously unbounded — latent, harmless for per-death
  credit, would not have been for this) and having `backtest --weighted`
  recompute fresh at every as-of step instead of slicing one precomputed
  table. Cost: ~40s for the full 45-quarter sweep, still cheap.

  **Still open:** xylazine never clears the dual gate below under
  `--weighted` — at its peak (n=7) EB05 is 1.46 and `eb05_own` is 1.32, both
  just under 1.5. Signal ceiling at low n (`RESEARCH_LOG.md`'s prior finding:
  three of five known emergences sit at n < 10, where Poisson noise binds, not
  the estimator), not a defect in the fix. `--weighted` stays opt-in, not a
  default.

**2. Joint dual gate for denominator drift — validated cleanly, no caveats.**
`eb05_own` and `credible_rise` columns on `rank`/`backtest`: a second GPS fit
whose `expected` is the substance's own baseline rate continued flat
(`n_base * recent_quarters/baseline_quarters`, no county-total term), run
alongside the existing share-based fit. "Credible rise" now requires both to
clear the alarm line. Formalizes what `regime()`'s `count_ratio` diagnostic
already showed by hand (cocaine: EB05 1.06 on a shrinking total, count itself
flat-to-down) into a labelled cell in the ranking output. All five known
emergences still clear the joint gate in `backtest`, at a cost of at most one
extra quarter of lag.

**Housekeeping this needed.** `causea_order`/`rollup_map` moved out of
`polysubstance.py` into `emerging/core/causeline.py`, once `trends.py` needed
the same cause-line ordering — same reasoning as the `core/spatial.py`
extraction in the reorg above, applied to a second case.

**3. A ground-truth reference for scoring `backtest`, and the metric it
enables.** `backtest` had no independent record of *when* an emergence
actually started or ended — only whether the ranking noticed it eventually.
Added `data/raw/ground_truth/known_emergences.csv`: one row per known local
emergence interval (a substance may have more than one), each carrying its
provenance. To avoid validating EB05 against labels derived from EB05, the
rule is fixed and reads only *raw* quarterly mention counts: a run of ≥2
consecutive quarters each ≥2, following ≥2 consecutive quarters each ≤1.
That rule fails outright on carfentanil and xylazine — both too sparse and
volatile at low n to ever post two qualifying consecutive quarters — so
those two intervals are hand-annotated, with the rule's failure and the
reasoning both written into the CSV. `national_context` carries cited
external evidence (DEA alerts, CDC MMWR, trade press, pulled this session)
for when each substance became a recognised threat *nationally*; it is
context for a reader, not ground truth for LA's own timing, since the
question this project asks is whether LA's local picture lags or leads the
national one.

`emerging/validation/ground_truth.py` (`interval_overlap`, temporal-IoU-style,
tested independent of any real data) scores the EB05 alarm sweep's detected
episodes against this reference. Result on the current five: **precision
1.00 across the board** — every alarm quarter these five substances ever
raised fell inside a real emergence window, no false alarms against this
reference. Recall is low for three of five (para-Fluorofentanyl 0.40,
mitragynine 0.13, xylazine 0.14), and **this is a construct mismatch, not a
detector failure** — the reference measures *presence* (still showing up in
LA deaths at all) while EB05 measures *growth* (still rising relative to its
own recent baseline), and a substance can stay present for years after its
growth phase ends. Para-Fluorofentanyl is the clean case: raw counts never
return to near-zero through 2025Q4, so the presence-based reference stays
open for 20 quarters, while EB05 correctly stops alarming once growth levels
off around 2023Q1 — which is also exactly what this log already reports
independently as a credible *decline* (EB95 0.37). Full table and discussion
in `docs/findings/ground_truth.md`. Open question, not resolved: whether the
reference schema should separately annotate an "actively rising"
sub-interval so recall doesn't need this caveat attached every time.

---

## 2026-08-18 — Reorganisation, and three reproducibility defects it exposed

Structural only: no method changed, and 24 result tables were regenerated
byte-identical to the committed ones to prove it (including the TreeScan binary
parity table and the 261x253 relrisk surface at max|diff| = 0). Tests held at
54 passed / 1 skipped throughout.

**What moved.** `CUTOFF` was defined three times (`trends`, `geo`,
`polysubstance`) and the palette three times, with five more modules importing
colours *from* `trends.py` — so rendering the relative-risk map pulled in the
empirical-Bayes ranking machinery. Both now live in `config.py` and `viz.py`,
which import nothing from the package. `load_points`, `UTM11N` and
`MAX_CASE_FRACTION` moved out of `geo.py` into `core/spatial.py`, since three
spatial analyses were importing them from the case-control clustering module.
The import graph now has exactly one cross-analysis edge — `treescan -> trends`
for `ALARM_THRESHOLD` and `KNOWN_EMERGENCES`, which is real.

Modules are grouped `ingest/ core/ analysis/ validation/`; `results/` is one
directory per module with figures beside their tables (`results/figures/` and
`FIG_DIR` are gone); findings prose moved from a 1,048-line README into
`docs/findings/<analysis>.md`. `emerging <module> <command>` is now a single
CLI, and `emerging pipeline` encodes the `alarms`-before-`plot` ordering that
was previously only a README paragraph.

**Three defects, all pre-existing, found by diffing regenerated output against
committed output:**

1. *`polysubstance_profile.csv` was not reproducible.* `sort_values` used
   pandas' default unstable quicksort on a tied column, fed by sets iterating
   in string-hash order, so tied rows reshuffled every run — two runs of
   unmodified code disagreed. A `git diff` of the committed CSV could not
   distinguish a real change from noise. Fixed with an explicit name tiebreak
   plus `kind="stable"`; verified across three `PYTHONHASHSEED` values, values
   unchanged.

2. *Same defect in `spacetime.py`.* Fixed the same way. It is why a last-bit
   float difference in one `llr` presented as four changed rows. That last-bit
   difference is itself real but environmental, not ours: running the scan on
   pre-reorganisation code at `b4ed99b` gives output identical to the current
   code on all 92 rows (max|llr diff| = 0), while *both* differ from the
   committed CSV by 4.44e-16 on the same single row — i.e. the committed file
   was written under a different numpy/scipy build. Harmless to conclusions,
   and no longer able to churn the table now that the sort is deterministic.

3. **The committed svtt per-substance p-values do not reproduce, and one of
   them crosses 0.05.** `svtt_clusters_fentanyl.csv` holds p = 0.63 where the
   default invocation gives 0.657; `svtt_clusters_methamphetamine.csv` holds
   **p = 0.02 where the default gives 0.054**. The test statistics are
   bit-identical (8.180210406528204), so it is purely the Monte Carlo draw:
   those files were generated with a seed or `--n-perm` other than the
   defaults, and the invocation was never recorded. Confirmed not to be a
   consequence of this reorganisation by running the scan on pre-reorganisation
   code in a separate worktree at `b4ed99b`, which reproduces 0.657 exactly.
   **Unresolved and deliberately not papered over** — regenerating would turn a
   significant cluster non-significant, which is a research decision, not a
   cleanup. Either re-run at defaults and accept 0.054, or recover the original
   invocation and record it in the file.

The general lesson is the one this log keeps relearning: a results directory
that is read to check the work has to be diff-stable, or the diff stops being
evidence.

---

## Current state (2026-08-10)

**What exists.** A working pipeline from LACME cause-of-death narratives to a
ranked, falsifiable, spatially- and pharmacologically-annotated list of
substances rising in LA County overdose deaths.

| Stage | Module | Output |
|---|---|---|
| Alias map | `lexicon.py` | 4,089 aliases → 3,127 canonicals |
| Extraction | `extract.py` | 219 substances, 99.8% decedent coverage, F1 0.962–0.997 vs the ME's own flags |
| Detection | `trends.py` | GPS/EB05 ranking, as-of backtest, alarm history, watchlist, regime diagnostics |
| Hierarchy | `tree.py` | 211 substances → 246 testable nodes (NFLIS categories + a pharmacological family layer) |
| Scan | `treescan.py` | tree-temporal scan, multiplicity-adjusted; head-to-head vs EB05 — **validated against the TreeScan C++ binary, LLRs identical to 5e-07** |
| Interpretation | `polysubstance.py` | cause-of-death vs passenger verdicts |
| Localization | `geo.py` | case-control random-labelling permutation test |
| Space-time | `spacetime.py` | Bernoulli + space-time-permutation scan over zip circles x recent quarters |
| Denominator | `census.py` | ACS resident population, 301 LA zips x 2013–2024 |
| Trend scan | `svtt.py` | spatial variation in temporal trends — where a *slope* differs from the county |
| Completeness | `nowcast.py` | reporting-delay curve from extract vintages, with a transferability guard |
| Risk surface | `relrisk.py` | Kelsall–Diggle log relative risk with CV bandwidth and tolerance contours |

**Headline results, settled window (recent 2025Q1–Q4 vs baseline 2023Q1–2024Q4;
2,153 vs 5,491 OD deaths):**

- GPS prior Gamma(α=7.267, β=7.547), fitted on 113 substances.
- **Nothing currently clears the EB05 > 1.5 alarm line.** Top: PCP 1.27,
  Alprazolam 1.11, Cocaine 1.06, Lidocaine 1.04.
- 24 of 205 substances have breached 1.5 at some point across 45 as-of quarters.
- Credible *declines*: Fentanyl (EB95 0.88), para-Fluorofentanyl (EB95 0.37).
- **The tree-temporal scan disagrees with the first line above.** After
  adjusting for all 158 nodes × 6 windows, four nodes clear RI ≥ 100:
  Lidocaine (RI 1,111), Dissociatives (455), Depressants and Tranquilizers
  (270), PCP (185). EB05 missed lidocaine only because its burst straddles the
  fixed 4-quarter boundary. Both hold under either reference series.

**The four findings that changed what we can claim:**

1. **Lidocaine is an adulterant, not a killer.** All 17 recent deaths name
   fentanyl and none name it alone. On the ME's combined cause line — which is
   ordered by clinical significance, not alphabetically — lidocaine is written
   **last in 14 of 17 deaths (82%) and first in none**. Fentanyl is the
   mirror image: first on 82% of multi-substance lines, last on 9%. The correct
   headline is "the fentanyl supply is being cut differently".
2. **Share growth ≠ death growth.** Total OD deaths fell 18% between windows,
   so a flat-count substance gains ×1.22 share for free. Cocaine posts EB05
   1.06 with `count_ratio` 0.90 — *a growing share of a shrinking total*.
   "Cocaine is rising" is the wrong sentence.
3. **Seven substances first breached in 2024Q3 simultaneously**, in a quarter
   when deaths fell 19% YoY, with panel density flat (sd 0.049). Seven
   independent epidemics is not the parsimonious reading; the denominator is.
4. **Medetomidine is absent from LA.** Zero mentions, as are dexmedetomidine,
   BTMPS and desomorphine — all four resolve in the lexicon, so the scan looked.
   LA has not had the wave that hit Philadelphia and Chicago. Caveat that
   cannot be resolved from this data: an assay the ME does not run produces an
   identical zero.

**Spatial findings.** PCP is 1.38× tighter than the overdose background
(p=0.005, n=367, anchored on 90013, Skid Row); the flagged *novel* substances
are not localized — they track the county-wide fentanyl supply. The space-time
scan sharpens this: **five substances are credibly concentrated somewhere
(cocaine 2.07× in South LA, PCP 2.81× around Watts, dissociatives, depressants
and tranquilizers, and cutting agents 3.42×), and none is credibly growing
anywhere** after correcting for the 46 scanned. The geography is real and it is
static.

**The one place that is changing** turns up only once a population denominator
is used and the test is on the *slope* rather than the level: against a county
falling at 9.3%/yr, the **Pomona Valley (91765, 91766, 91768, 91748, 91789) is
rising at 25.2%/yr**, p = 0.016 adjusted for all 23,406 circles. It started at
a third of the county rate and has converged upward, so it is a catch-up rather
than a new epicentre — and invisible to any level-based scan. Read the caveats
in the P7 entry before quoting it: the rise is concentrated in 2024, and the
cluster *boundary* is not stable to dropping 2025 even though the region is.

### How EB05 and the tree scan divide the work

Read this before quoting either one. They are **complementary, but not in the
direction the game plan assumed** — P1 was written expecting the scan to be a
better detector. It is not. It is a better adjudicator.

| | EB05 | Tree scan |
|---|---|---|
| Detects earlier (matched 85-alarm budget) | — | never |
| Known emergences found | 5/5 | 3/5 |
| Effect size with a credible interval | yes | no |
| Detects credible *declines* | yes | no |
| Multiplicity-adjusted | no | yes |
| Free of a window choice | no | yes |
| Can test a family | no | yes |
| Cost of a 45-quarter sweep | ~1 s | ~2.5 min |

**What only EB05 can do.** It *estimates*; the scan only *tests*. EBGM with
EB05/EB95 is a shrunken effect size with an interval, so 205 substances can be
ranked; p-values cannot rank, because half of them tie at 1.0. EB95 < 1 is the
only way we detect declines (fentanyl 0.88, para-fluorofentanyl 0.37) — the
scan is high-rate only. And it stays sensitive at the counts that matter:
xylazine reaching rank 4 on n=6 is the headline EB05 result, where the scan
peaks at RI 13.

**What only the scan can do.** Adjust for having looked at 205 substances × 45
quarters. Search window length, which the log already shows is load-bearing
(lidocaine reads 1.04 / 2.20 / 4.29 depending on where the 4-quarter boundary
falls). Test a branch. And decline the shrinking-denominator artifact
unprompted — cocaine is non-significant under both reference series without
being told about `count_ratio`.

**Why neither wins, and what that means.** This lands where the signal-ceiling
work already did: three of the five known emergences sit at n < 10, where
Poisson noise and not the estimator is binding. Family-wise error control costs
exactly the power needed for a substance with nine deaths. The scan does not
lower that ceiling — it makes it visible. Its real contribution is honesty
about how much of EB05's sensitivity is bought with an uncontrolled
false-positive rate.

**Operating procedure: EB05 screens and ranks, the scan adjudicates.** Report
the pair, and read the three cells:

| | RI ≥ 100 | RI < 100 |
|---|---|---|
| **EB05 > 1.5** | a real signal; defensible to a reviewer | a candidate to watch, not to publish — this is where EB05's unadjusted threshold was quietly doing the work (xylazine, mitragynine) |
| **EB05 < 1.5** | EB05's fixed window missed it (lidocaine, PCP) | nothing |

The two off-diagonal cells are the new information, and the bottom-left one
matters most right now: the current headline is "nothing clears the line",
and lidocaine and PCP contradict it.

*Correction to an earlier framing in this log:* the scan raising 165
node-quarter alarms against EB05's 85 was reported as "twice the noise for no
detection gain". True at RI ≥ 100, but that is a threshold choice, not a
property of the method — at a matched budget it is 11 nodes. The finding that
survives is the first row of the table above: never earlier, and 3 of 5.

---

## Game plan

Priorities from `LITERATURE_REVIEW.md`. Ordered by expected value, not effort.

### P1 — TreeScan, prospective tree-temporal `done 2026-08-10`

Built as `tree.py` + `treescan.py`; see the entry below for what it showed. Of
the three things it was supposed to fix, **one delivered, one delivered
partially, one failed**:

- **Multiple testing — delivered.** p-values now come from the Monte Carlo
  distribution of the maximum LLR over 158 nodes × 6 windows, calibrated to
  within Monte Carlo error. This is the whole return on the work.
- **Repeated looks — partial.** The recurrence interval covers multiplicity
  *within* an analysis, not across quarterly looks. The claim below that
  MaxSPRT is redundant with it was too strong.
- **Aggregation power — failed.** The nitazene branch never signals (peak LLR
  3.38, RI 1). Aggregation helps for individually *marginal* members; the
  nitazenes are individually *absent*.

It **runs alongside EB05 rather than replacing it** (step 4 of the original
plan): at a matched alarm budget it never detects earlier and misses two of the
five known emergences. EB05 screens and ranks, the scan adjudicates — see "How
EB05 and the tree scan divide the work" above for the operating procedure and
the three-cell rule for reading them together. The risk noted here — "standalone
C++ binary, breaks the pure-Python pipeline" — was avoided by implementing the
method directly; it is ~400 lines and runs the full 45-quarter sweep in 2.5
minutes, and has since been validated against that binary.

### P2 — space-time scan `done 2026-08-11`

Built as `spacetime.py`; see the entry below. **Answered, and the answer is
no.** Five substances are strongly and credibly *concentrated* somewhere;
**none is credibly growing anywhere**. The distinction only exists because the
first model could not draw it — see the log entry for why a second one was
needed.

The space-time permutation applied to a substance's own case series is the
*growth* test and was built, but as the secondary model, because on its own it
mostly rediscovers where overdose deaths are. The primary is a Bernoulli
case-control space-time scan, which keeps the denominator as overdose death
itself and so stays commensurable with EB05, `treescan.py` and `geo.py`.

**Correction, 2026-08-11.** This entry originally also ruled out spatial
variation in temporal trends on the grounds that it "needs a population
denominator per zip per quarter, which is the thing this project has never
had". **That was wrong.** ACS `B01001_001E` is public and free, and
`predict_rosla/rosla_nowcast/data_prep/acs_ingest.py` was already fetching it
and discarding it in `_derive`. `census.py` now keeps it: 301 LA zips x
2013–2024, and **every one of the 280 zips holding deaths has a population**.
SVTT is therefore back on the table — see P7.

### P3 — Delay-distribution nowcast `done 2026-08-11`

Built as `nowcast.py`. **The delay curve is measurable and it does not
transfer**, so the cutoff stays where it is — see the entry below. The warning
in the original note turned out to be the finding: the regime that shifted was
the *reporting* process, not the epidemic.

### P4 — relative-risk surface for PCP `done 2026-08-11`

Built as `relrisk.py` — in Python rather than the planned R/`sparr` detour, for
the same reason `treescan.py` was. Turns "1.38x tighter" into a map, and the
map disagrees with the summary in an informative way. See the entry below.

### P5 — READUS-PV reporting conventions in the manuscript `not started`

### P7 — spatial variation in temporal trends, on a real denominator `done 2026-08-11`

Built as `census.py` + `svtt.py`; see the entry below. **The first
population-denominated finding the project has produced:** against a county
falling at −9.3%/yr, the Pomona Valley is rising at +25.2%/yr, p = 0.016
adjusted for all 23,406 circles.

Original notes, kept because the caveats still bind:

Unblocked 2026-08-11 by `census.py`. Every method in this repository so far is
**compositional** — what share of overdose deaths names a substance — because
the ME extract carries no denominator. With population per zip per year, two
things become possible that were not:

1. **Rates.** "N deaths per 100,000 residents" is the sentence every reviewer
   and every health officer asks for, and nothing here can currently say it.
2. **SVTT** (Kulldorff 2006): fit a log-linear Poisson trend inside and
   outside each circle and scan for where the *trend* differs. Unlike a
   level-based scan it is not confounded by where overdose mortality is
   simply high, which is the objection that made `geo.py` case-control in the
   first place.

The single highest-value run is **all overdose deaths, not per substance**:
"where in LA County is overdose mortality rising fastest per capita" is
directly actionable and nothing in the pipeline answers it. Per-substance SVTT
is a follow-on for the handful with enough counts.

Cost estimate: the MLE for a log-linear Poisson trend reduces to a
one-dimensional monotone root find, so it vectorises across all ~27k circles
at once; ~20 Newton iterations per replicate. Feasible, but it is real work
and needs the same calibration discipline as P1 and P2.

Caveats to carry into the results, all recorded in `census.py`: ACS 5-year
estimates are already smoothed and lag fast change (Skid Row is the case to
worry about); the series ends in 2024 so 2025 is carried forward; ZCTAs spill
across the county line, which is why the total reads 10.1M against the
county's ~9.7M; and **residential population is the wrong denominator for a
place people travel to** — 90013 has 14,891 residents and a daytime
population several times that, so a per-resident rate there is an
overestimate of individual risk and should be read as place-based intensity.

### P6 — release the scan as a Python package `deferred 2026-08-11`

**Come back to this.** There appears to be no Python implementation of the
tree-based scan statistic: nothing on PyPI under `treescan`, `pytreescan`,
`scanstatistics` or the other obvious names, and R's `TreeMineR` covers only
the tree-*only* variants (Kulldorff 2003/2013), not the tree-temporal one. The
official tool is a registration-walled C++ binary. `tree.py` + `treescan.py`
are ~700 lines with no dependency beyond numpy/pandas, and now match the
reference binary exactly (see the 2026-08-11 entry), so the raw material is
there.

Deferred rather than dropped, because two things need doing first and neither
is a side quest off the LA analysis:

1. **Scope honesty.** We implement one corner of TreeScan's surface —
   conditional-on-node, prospective, high-rate. TreeScan also does
   unconditional / total-cases / node-and-time conditioning, Bernoulli and
   Poisson models, retrospective windows, low-rate scans and sequential
   analysis. A release either says so plainly or grows to cover more.
2. **The two things that make ours worth having are extensions, not a
   reimplementation** — the arbitrary reference series and the overlapping-node
   DAG. Neither is in the literature. They would have to be presented as such,
   and the reference-series generalisation validated against a simulated ground
   truth, because it is precisely the part the binary cannot check.

Licensing is not an obstacle: the method is published science and our code was
written before the binary was ever downloaded. Do not vendor any part of their
distribution.

### Open questions needing outside input

- **Which policies?** `regime` identifies *when* something structural happened
  (2024Q3, and the pre-2020 panel widening). Matching those to naloxone
  saturation, scheduling changes, or an ME protocol revision needs sources
  outside this dataset.
- **Does the ME assay for medetomidine and BTMPS?** Determines whether four
  zeros mean absence or blindness. One question to the office resolves it.
- **Is the PCP/Skid Row concentration actionable**, or already known to the
  outreach teams working that area?

### Explicitly not doing

- **MaxSPRT** — was "redundant with TreeScan's recurrence interval"; that is
  **only true within a single analysis**, and building P1 made the distinction
  concrete. The RI does not adjust for looking every quarter. Still not worth
  doing — the sequential inflation is bounded and we report RI, not a decision
  — but the reason has changed and the old one was wrong.
- **LGCP / inlabru** — prior-dominated at n=9–39; revisit only if a single
  substance becomes a funded investigation.
- **Two-component MGPS** — not identifiable at this scale; documented in
  `_fit_gps_prior`. Revisit on a substance × zip × quarter table.
- **Spatial methods as *detectors*** — they stay downstream of the temporal
  detector. At n=9 there is nothing to scan.

---

## 2026-08-11 — P4, Kelsall-Diggle relative-risk surface

`relrisk.py`. Cases are deaths naming the substance, controls the other
overdose deaths, so r(x) = log(f/g) is 0 where the substance is represented
exactly as it is county-wide. Same compositional question as everything else,
answered continuously instead of as one number (`geo.py`) or one circle
(`spacetime.py`). Written in Python rather than driving R's `sparr`.

**Bandwidth was the objection and it is now measured, not assumed.** `geo.py`
rejected kernel methods partly because "at n=9-39 that choice would drive the
answer", and it does: PCP's peak relative risk runs 3.7x at 1 km to 2.0x at
6 km. The default is therefore chosen by **likelihood cross-validation on the
case/control labels** (Kelsall-Diggle 1995, what `sparr::LSCV.risk` does),
which targets the *ratio* rather than either density. It picks **4.0 km**, with
a flat interior optimum — 3.5, 4 and 5 km are within 2 log-likelihood units of
each other, which is itself the finding: PCP's excess is a broad regional
gradient, not a hotspot, and the data will not support fine structure.
Silverman's normal-reference rule says 4.19 km but is not used; it assumes a
Gaussian cloud and LA is a 100 km multi-centre sprawl.

**The map disagrees with `geo.py`'s headline zip, and both are right.** `geo.py`
reports 90013 (Skid Row) as PCP's top zip — that is the *modal* zip by raw
count (23 deaths). By *share of local overdose deaths* the peak is elsewhere:
90044 South LA at 15.2%, 90744 Wilmington 14.3%, 90806 Long Beach 13.8%,
against 4.84% county-wide. The surface accordingly puts the significant excess
in a contiguous **South LA to Harbor corridor** (90731, 90744, 90732, 90802,
90813, 90044) at roughly 2.1x, with the valleys and the Antelope Valley
depleted. Skid Row is inside an elevated downtown area but is not the peak —
it has the most PCP deaths because it has the most deaths.

12.6% of the county is above the county-wide rate at pointwise p <= 0.05.
These are **pointwise** contours over 41,000 cells with no multiplicity
correction; read contiguous regions, and go to `spacetime.py` for a single
adjusted verdict.

*Three defects found while building, all of which produced a plausible-looking
wrong map:*
- `eps = 1e-12` floored the log-ratio at about -27 wherever a cell had zero
  cases, so a handful of empty-desert cells set the colour scale from -20 to
  +20 and hid the real +/-1 range. Replaced with a half-death pseudo-count
  centred so that f/g = p_global gives exactly 0 at any local density. This
  also regularised the upper tail: bandwidth sensitivity fell from 13.4x-2.0x
  to 3.7x-2.0x, so a chunk of the "bandwidth drives everything" problem was
  really this bug.
- No low-density mask, so the ratio was being estimated where almost no deaths
  informed it. Now requires >= 5 deaths within one bandwidth.
- `scipy.ndimage.gaussian_filter` returns the *input* dtype, so an integer
  count grid comes back truncated — a smoothed density of 0.3 becomes 0. Caught
  by the calibration test firing at 84% instead of 5%. `log_rr` now coerces.

---

## 2026-08-11 — P3, reporting-delay nowcast

`nowcast.py`. `CUTOFF = 2025-12-31` has always been a judgement with nothing
behind it, and it decides how much data the project throws away. This measures
the thing it was guessing at.

**LACME has no report date**, only `DeathDate`, so delay cannot be recovered
from one extract. It can from two *vintages*: sweeping death-months at a single
pull date traces the whole curve, because each month sits at a different delay.
Three overdose extracts exist. The 2024-08 / 2025-05 pair is usable — its
settled era (2012–2022) agrees to **0.2%**, so the difference is accrual. The
2025-05 / 2026-04 pair is **not**: its settled era moves **3.5%**, which cannot
be reporting delay four to fourteen years later. That is a case-definition or
extraction change, and folding it into a delay estimate would inflate every
correction. `triangle` prints this drift for whatever pair it is handed.

**The curve** (clean pair, isotonic): 27% complete in the month of death, 66%
at one month, ~80% at three to five, ~90% at seven to eight, ~100% from ten.
Its saturating end is independently confirmed in the *current* era — the
2026-02 and 2026-04 pulls agree to within 1.4% on every month of 2025Q1–Q2.

**The result I nearly shipped, and why it was wrong.** Applied to the current
extract the curve says 2025Q4 is 81% complete (560 observed → ~690) and 2025Q3
is 93%, which would mean `CUTOFF` is too late, the recent window is undercounted
by ~7%, and the monotone 2025 decline (681, 628, 602, 560) is really a decline
and rebound (681, 629, 650, 690). That is a headline-changing claim and it does
not survive its own check.

**`transfer_check` refuses it.** The current extract runs at full strength
through 2026-01 — 198 deaths, in line with every month before it — and then
drops to **54, 18, 6**. The curve predicts a ramp over those months (79%, 77%,
66%); what is there is a **cliff**, a 3.7x fall where 1.2x was predicted. No
smooth accrual process produces that. It is a pull boundary with a few
fast-closing cases trailing it, and it means today's reporting is not the
process measured in 2024. The correction is computed, printed, and withheld.

The practical reading inverts: the extract looks essentially complete through
2026-01, so **`CUTOFF = 2025-12-31` is more likely conservative than lax** and
stays put. Nothing downstream changes.

*What would settle it, and it is free:* two current-era vintages. Save a dated
copy of every future extract. Without that this question cannot be reopened,
however much modelling is thrown at it.

*Defect found by the tests:* an empty settled window made the drift check
return NaN, which happened to reject — for the wrong reason, and silently.
It now raises.

*Also fixed:* the isotonic fit was applied to a reversed array, pooling the
entire curve to a single constant 0.8636. It looked like a plausible flat
completeness and would have produced a uniform 16% inflation of every quarter.

---

## 2026-08-11 — P7, spatial variation in temporal trends

`svtt.py`, on the denominator `census.py` unlocked. Every other scan here asks
where a *level* is unusual; this asks where a *slope* is. Each circle keeps its
own baseline rate and is tested only on its trend, so "this area has always
been bad" cannot signal — which is the objection that made `geo.py`
case-control in the first place, now handled inside the model rather than by
avoiding the denominator.

The MLE of a log-linear Poisson with a population offset reduces to a
one-dimensional monotone root find: the population-weighted mean time under
trend *b* must equal the case-weighted mean time. Its derivative is a variance,
so the root is unique and Newton converges from anywhere, which is what makes
23,406 circles affordable without a GLM per circle.

**Caught on the first run: the window was misspecified.** With the whole
2015–2025 series the county trend came out at **+14.5%/yr**. LA's overdose rate
is a hump, not a trend — 5.3 per 100k in 2013, peaking at 29.5 in 2022, down to
21.7 in 2025 — so a log-linear fit across all of it describes neither phase,
and "rising fastest" would have conflated *rose earliest during the ramp* with
*still rising now*. The default is now the decline phase (2022–2025, county
−9.3%/yr, monotone), and `scan` prints the annual rate series it fitted so this
cannot recur silently.

**The finding.** A cluster of 5 zips — 91765 Diamond Bar, 91766/91768 Pomona,
91748 Rowland Heights, 91789 Walnut — rising **+25.2%/yr while the county
falls 9.8%/yr**. p = 0.016, adjusted for every circle searched. Verified by
hand, independent of the scan code:

| Year | Cluster /100k | Rest of county |
|---|---|---|
| 2022 | 11.5 | 30.0 |
| 2023 | 12.5 | 28.8 |
| 2024 | 23.0 | 26.3 |
| 2025 | 20.4 | 21.7 |

**It is a catch-up, not a new epicentre.** The cluster started at a third of
the county rate and converged upward; it is still slightly *below* county
average in level. A level-based scan — including a population-denominated one
— would miss it completely, which is exactly the case SVTT exists for and the
reason P7 was worth doing rather than just computing rates.

*Three things to say alongside it, not after being asked:*
- **The rise is a single year.** 2023→2024 is the whole step; 2025 falls back.
  A log-linear slope summarises that shape, it does not describe it. Quote the
  table, not just the +25%.
- **The boundary is not stable.** Dropping 2025 moves the top cluster from
  this tight 5-zip core to a broad 58-zip eastern-county region (+7.9% vs
  −7.9%, LLR 8.50 against 10.15). Every Pomona zip is inside both, so the
  *region* is robust and the *extent* is not. Report "eastern San Gabriel /
  Pomona Valley", not five zipcodes.
- **156 deaths.** Enough for a slope, not enough to subdivide.

*Calibration.* `tests/test_svtt.py` checks the trend MLE against the score
equation it claims to solve and against a planted slope, that the LLR is
exactly zero when two circles share a slope but differ 10× in level, that it is
one-sided, and that a planted cluster is recovered. The calibration test runs
the null with a deliberate **5× spread in baseline level across zips** — a
level-based scan would light up on that, and SVTT must not. 8 tests pass.

*Follow-up, same day — the cluster profiled and the substance identified.*
`svtt profile` characterises the most likely cluster against the rest of the
county, because "a place is rising" is not reportable until you know what is in
it and whether it is an artifact.

- **It is Pomona, not five zips.** 91768 (+19) and 91766 (+16) supply 35 of the
  cluster's 44-death increase; 91789, 91765 and 91748 contribute 9 between them
  and are along for the ride via the circle geometry. Report "Pomona", with the
  neighbours as context.
- **A level shift, not a burst.** Quarterly: 2022–23 runs 2–9 deaths, 2024–25
  runs 6–17, with the step at 2024Q1 and no reversion. One incident or one
  batch of late reporting would not look like this.
- **Not an artifact of where people are pronounced.** 3.8% of cluster deaths
  are at shared facility addresses against 18.6% county-wide, and the deaths
  are at 26/27/44/47 distinct addresses across the four years. These are homes.
- **No substance is driving it.** Cluster-vs-county shares: heroin 1.46x, PCP
  1.44x, methamphetamine 1.16x, fentanyl 0.89x, cocaine 0.74x. The spread is
  narrow; this is more of the same supply, not a new drug.
- **Who:** 55% LATINE against 38% county-wide, median age 48 against 44, sex
  identical. Pomona is a heavily Latino city so this partly reflects who lives
  there, but it has obvious implications for who any outreach has to reach.

*Per-substance SVTT, and it separates the two epidemics cleanly:*

| | county trend | most likely cluster | p |
|---|---|---|---|
| Fentanyl | **−16.25%/yr** | 91765, 7 zips, +3.8%/yr | **0.63** |
| Methamphetamine | **−7.56%/yr** | 91741, 42 zips (Glendora–Azusa–Claremont–Pomona), **+5.98%/yr** | **0.02** |

**Fentanyl is falling uniformly everywhere** — no significant spatial variation
in its trend anywhere in the county. **Methamphetamine is not**: it falls half
as fast countywide and is *rising* across a broad eastern San Gabriel Valley
region that contains the all-cause cluster. So the eastern-county rise is a
methamphetamine phenomenon, and the county-wide decline is fentanyl-led.

That distinction changes what a response would be. A naloxone-and-fentanyl
-test-strip package is aimed at the part of the epidemic that is already
receding fastest; the part that is growing, in the place it is growing, is
stimulant-involved and has no equivalent pharmacological answer.

---

## 2026-08-11 — P2, space-time scan

`spacetime.py`. Circles of zipcodes crossed with trailing quarters; 7,572 OD
deaths, 280 zips, 12 quarters (2023Q1–2025Q4); circles capped at 25% of deaths
(≤ 97 zips), windows 1–6q, 999 replicates. 46 of 157 substances and tree
branches were eligible (≥ 20 cases, ≤ 25% of all deaths).

**Two models, because the obvious one cannot answer the question P2 asked.**
The Bernoulli case-control scan — cases are deaths naming the substance,
controls are the other OD deaths in the same zip and quarter — measures the
substance's share against its share pooled over the *whole* study period. A
concentration that has been there for a decade therefore wins as easily as one
that appeared last year. It flagged cocaine at 2.07× over a 5-quarter window in
South LA, which reads as an emerging hotspot and is not one. So the second
model is Kulldorff's space-time permutation, conditioned on the substance's own
spatial margin *and* its own temporal margin: a busy zip cannot win, a busy
quarter cannot win, only the interaction can. That is "growing here" as opposed
to "concentrated here".

**The result is a clean split.**

| | concentrated (Bernoulli) | growing (interaction) |
|---|---|---|
| p ≤ 0.05 | 21 of 46 | 1 of 46 |
| survives Bonferroni across the 46 | **5** | **0** |

The five: Cocaine 2.07× (137 vs 66.1, 12 zips around 90011, 5q), PCP 2.81×
(68 vs 24.2, 37 zips around 90002, 4q), Dissociatives 2.25×, Depressants and
Tranquilizers 2.12×, and Cutting agents and adulterants 3.42× (22 vs 6.4, 42
zips around 90043, 6q). **The spatial structure of the LA supply is real and
strong, but static — these substances have a geography, not a moving front.**
PCP corroborates `geo.py` on an independent statistic, and the adulterant
family is the spatial echo of the lidocaine finding.

*Had only the Bernoulli model been built, this entry would have claimed the
opposite.* That is the transferable lesson: with trailing windows, "recent
elevation" and "long-standing concentration" are the same picture unless you
condition them apart.

**Multiplicity, again.** Each p-value is adjusted for every cylinder searched
*within* a substance, and not at all for having scanned 46 of them. Reported
per-substance, the interaction model has a hit — diphenhydramine, p = 0.003,
4 deaths against 0.44 expected across the affluent Westside (90024, 90049,
90210, 90272 …) over 3 quarters. Corrected across the 46, it is 0.138. It is
not a facility artifact (11% of its deaths at shared addresses, against 15%
county-wide), and it is not a finding either. `p_across_substances` is now a
column so nobody has to remember this.

*Two defects found and fixed while building:*
- The observed grid was built with `np.add.at` and the simulated grids with
  `np.bincount`. They agreed, but only by luck of construction — if they had
  not, the p-values would have been wrong with nothing failing. Now asserted.
- The reported member-zip list was truncated to 12, so the map drew a 12-zip
  blob for a 69-zip cluster. The statistics were right; the picture was not.

*Also recorded:* methamphetamine (61% of deaths), fentanyl (58%) and their
parent nodes are withheld by the same `MAX_CASE_FRACTION` guard as in
`geo.py` — 111 of 157 nodes were skipped, most for having too few cases.

*Calibration.* Synthetic nulls for both models are in
`tests/test_spacetime.py`, along with a planted-cluster recovery test. On the
real background, rejection rates are 0.075 / 0.113 at nominal 0.05 / 0.10 for
C=40 and 0.037 / 0.087 for C=400 — within Monte Carlo error at the thresholds
that matter (80 replicate datasets, so sd ≈ 0.024 and 0.034). The body of the
distribution runs mildly conservative, which errs toward under-claiming.

---

## 2026-08-11 — Validated the scan against the TreeScan binary

Downloaded TreeScan v2.4 (treescan.org, registration required) and ran it
against ours on the same tree and the same counts. **They compute the same
statistic.**

| | strict tree | DAG (production) |
|---|---|---|
| Nodes in the search space | 141 both | 155 both |
| Cuts compared exactly | 57 | 65 |
| Max abs LLR difference | 5.0e-07 | 5.0e-07 |
| Best window / window counts identical | 57/57 | 65/65 |
| Max abs p difference (9,999 reps each) | 0.0041 | 0.0039 |

5e-07 is TreeScan's printf precision — the LLRs agree to the last digit either
program prints. This closes the main open liability from the P1 entry below:
the unit tests showed we were self-consistent and calibrated, neither of which
would have caught computing the wrong statistic consistently.

**The comparison is only possible in one configuration, and that is worth
recording.** TreeScan's tree-temporal model spreads each node's cases
*uniformly* over the study period. Ours generalises that to an arbitrary
reference series, because a quarter with 18% fewer overdose deaths is not a
quarter with 18% less time — which is the entire point of `--reference deaths`.
So the validation runs ours with a uniform reference, where the two coincide.
Everything except the generalisation is checked; the generalisation itself is
one line (`_window_shares`) with its own unit test.

**One real bug found.** Metoprolol: 3 observed against 3.000 expected, LLR
3.3e-16, because `c > e` came out true on the last bit of a float. It could
never win a maximum, so no published number moves, but it put a node with no
excess into the positive-LLR list. `_llr` now floors at 1e-9.

**Four settings had to be matched**, each of which would have looked like a
disagreement: `conditional-type=2` (node, not node-and-time);
`apply-risk-window-restriction=n` (defaults *on* at 20%, silently dropping
short windows); `include-identical-parent-cuts=y`; and `minimum-node-cases=2`,
because TreeScan forms no cut from fewer than 2 cases. That last one is a
property of the *search space*, not of reporting — left unmatched, our null
maximum would range over nodes TreeScan never sees.

*Practical note for anyone re-running this:* TreeScan's input files are
comma-delimited with no quoting, so `1,1-Difluoroethane` parses as an extra
column and the run dies with "record 17 has node referencing self as parent".
The exporter writes opaque IDs and maps names back. It also recomputes every
node's leaf closure *from the emitted file* and refuses to run if it disagrees
with ours, so a mistranslated DAG fails loudly rather than quietly comparing
two different trees.

**Performance** (155 nodes, 12 quarters, single-threaded): 999 reps — TreeScan
2.13s, ours 0.43s; 9,999 — 2.25s vs 3.5s; 99,999 — 12.4s vs 34s. Different
shapes: TreeScan carries a ~2s constant then scales well and uses all cores
(99,999 in 2.5s multi-threaded); ours is near-linear at ~0.35 ms/replication
and single-threaded. Neither is a constraint at the counts we need.

**Decision: keep ours, drop the binary from the workflow.** It stays in
`treescan/` (gitignored — it is Kulldorff's, under its own licence) purely so
`tests/test_treescan_parity.py` can re-run the check. Those tests skip when it
is absent, so the repository still stands alone.

*Open, not yet acted on:* no Python implementation of the tree-based scan
statistic appears to exist — nothing on PyPI under the obvious names, and R's
`TreeMineR` covers the tree-only variants (Kulldorff 2003/2013), not the
tree-temporal one. Extracting `tree.py` + `treescan.py` as a standalone package
looks worthwhile, but the reference-series generalisation is ours and not in
the literature, so it would need to be presented as an extension rather than a
reimplementation, and validated on more than one dataset first.

---

## 2026-08-10 — P1, tree-temporal scan

Implemented the method rather than driving the TreeScan binary, so the pipeline
stays one language. ~400 lines; the 45-quarter sweep runs in 2.5 minutes.

**Calibration first, because everything else depends on it.** On null data the
rejection rates are 0.050 / 0.096 / 0.183 / 0.496 against nominal 0.05 / 0.10 /
0.20 / 0.50. Asserted in `tests/test_treescan.py`, together with the specific
way it would break silently: applying the `min_cases` floor to the observed
data but not to the simulations, which would bias the null downward and
manufacture significance everywhere.

**Positive result: it changes the current answer.** EB05 says nothing clears
1.5. After adjusting for the entire search, four nodes clear RI ≥ 100 —
Lidocaine 1,111, Dissociatives 455, Depressants and Tranquilizers 270, PCP 185.
Lidocaine was missed only because its burst straddles the fixed 4-quarter
boundary; scanning window length 1–6 finds it at 6q, 31 observed vs 17.4
expected. PCP corroborates the `geo.py` Skid Row finding on an independent
axis. Both survive under either reference series, and **cocaine — the EB05
near-miss that `count_ratio` exposed as a growing share of a shrinking total —
is nowhere near significant under either.** The scan declines the artifact
without being told about it.

**Negative result: it does not detect earlier.** EB05 > 1.5 is unadjusted and
RI ≥ 100 is adjusted for the whole search, so comparing them directly rewards
whichever is more trigger-happy. The sweep therefore also re-cuts the scan at
the threshold raising the *same* 85 alarms over the sweep. At that matched
budget the scan never detects earlier (para-fluorofentanyl +1q, bromazolam
+12q vs 11, carfentanil +29q) and **misses mitragynine and xylazine entirely**
(leaf peaks RI 48 and 13).

**Negative result: the aggregation argument, which is what motivated P1, did
not survive.** The nitazene branch never signals — peak LLR 3.38, RI 1.
Aggregation pools *marginal* members; the nitazenes are *absent*, 1 and 2
mentions in fourteen years, and no pooling makes signal from three cases.
Branch-only detection — a branch signalling in a quarter when no member leaf
did — occurs in **3 node-quarters out of 165** across the whole sweep.

**Defect found and fixed: branch signals are not attributable by default.**
The first head-to-head reported carfentanil detected at +0q via `Fentanyl and
Fentanyl-related` in 2017Q2, and para-fluorofentanyl at −18q — a branch
"detecting" a substance six years before it existed. Two corrections: branch
signals now require `as_of >= first_seen`, and `excess_share()` reports what
fraction of the branch's excess the substance actually contributes. That turns
carfentanil's apparent 0-quarter detection into **1% attributable** (it was
fentanyl's own explosion) and mitragynine's into ≈0%. Only para-fluorofentanyl
earns a genuine branch detection, 74% attributable via `Fentanyl analogs`.
Without this the tree scan manufactures detections for any rare substance
filed under a large, busy branch.

**The tree turned out to be the hard part, not the statistic.** The plan called
for `root → NFLIS category → substance` because the field is already populated.
The NFLIS categories are a forensic-chemistry filing system: `Depressants and
Tranquilizers` puts PCP (724) beside xylazine (9), `Narcotic Analgesics` puts
morphine (1,124) beside the nitazenes (3), and `Other substances` is a
71-member junk drawer holding lidocaine and levamisole next to atorvastatin.
Aggregating the nitazenes to their NFLIS parent *loses* the signal. So
`tree.py` keeps the categories and adds a pharmacological family layer beside
them. Families cross category boundaries, which makes it a **DAG, not a tree** —
fine, because the LLR is per node and the Monte Carlo maximum absorbs the
overlap; the strict-tree requirement belongs to TreeScan's implementation, not
to the method.

*Recorded so it is not re-litigated:* `Veterinary sedatives` collapses to
`{Xylazine}` and is dropped as a duplicate of the leaf, because medetomidine
and dexmedetomidine are absent from LA. Absent members stay in the constant on
purpose — the branch then already exists in the quarter one is first coded.

*Also worth knowing:* `Fentanyl and Fentanyl-related` is in signal in 29 of 45
quarters and `Fentanyl` in 28. A node that is almost always in alarm carries no
information; the fentanyl branches should be read as background, not events.

---

## 2026-08-10 — Reviewer feedback round

Six items from a presentation. All six built.

**Polysubstance / principal cause** (`polysubstance.py`, new). Two independent
discriminators. Co-occurrence (`alone_pct`, `with_fentanyl_pct`) and — the
better one — **position on the ME's combined cause line**, which is ordered by
clinical significance, not alphabetically. Reported three ways: `lead_pct`
(named first), `last_pct` (named last), and `mean_rank_pct`, the average
position on a 0 = always first to 1 = always last scale. Fentanyl 0.13,
lidocaine 0.93.

*Why the position is normalized rather than a mean ordinal (1…N).* Line length
varies with substance type, and a raw mean position conflates "named late" with
"appears on a long line" — a substance always named last scores 2 on a
two-substance line and 5 on a five-substance line. The bias runs the wrong way
for us: adulterants co-occur with more substances, sit on longer lines, and so
would be inflated by a raw ordinal precisely where we most need the measure to
be trustworthy. Ranking all 26 substances both ways gives **Spearman 0.79**,
with large disagreements — methamphetamine is 20th of 26 raw but 7th
normalized (short lines, so position 2 of 2 is *last*), and
para-fluorofentanyl is 6th raw but 15th normalized (long lines, so position
2.83 is the middle).

*Cost of normalizing, recorded so nobody re-derives it:* on a two-substance
line the score can only be 0 or 1, so short lines push toward the extremes;
single-substance lines are undefined (0/0) and excluded from both `lead_pct`
and `mean_rank_pct`; and it is an unweighted mean of ratios, so a 2-substance
line counts as much as a 5-substance one. `last_pct` was added alongside
because the normalized mean is the right thing to *rank* on but hard to read —
"0% first, 82% last" states lidocaine's case more plainly than 0.93 does.

*Tried and rejected:* which **field** the substance appeared in. Lidocaine is
100% `CauseA`, identical to fentanyl, because the ME puts every substance on
one line. Position within the line is signal; presence in it is not.

*Defect found and fixed:* `lead_pct` initially counted solo mentions as "named
first" while `mean_rank_pct` excluded them, so methamphetamine read as both 49%
first and near-last. Both statistics now use the same multi-substance
denominator.

**Alarm history** (`trends alarms`, new). Full sweep of every substance across
all 45 as-of quarters; 24 have ever breached EB05 1.5. Peaks: pFF 7.00,
Carfentanil 5.02, Flualprazolam 4.84, Bromazolam 4.78, Etizolam 4.49,
Lidocaine 4.29.

*Correction built in:* `first_seen` is left-censored at 2012Q1, so morphine's
apparent "53-quarter detection delay" is the censoring floor, not a late catch.
Censored rows now carry a caret marker and a NaN lag rather than a fake number.

**Regime diagnostics** (`trends regime`, new). Tested and **rejected** the
obvious hypothesis: panel expansion is *not* driving recent flags —
substances-per-death has been flat at ~1.89 since 2020 (sd 0.049). The real
artifact is the denominator (finding 2 above). Panel drift is historical:
1.45 (2012) → 1.89, plateauing 2020Q1, which is why the heatmap marks that line.

**Watchlist** (`trends watch`, new). Generalized a one-off medetomidine query
into a standing check of nine nationally-tracked substances plus the nitazene
family. Reviewer's "metatodanine" was medetomidine.

**Heatmap extended to 2012.** Reveals the heroin→fentanyl transition and the
methadone decline, neither visible from 2015.

**Geo clustering** (`geo.py`, new). Case-control random-labelling permutation
test. Three data defects found while building it, all now handled: **46 deaths
geocoded to Namibia** (caught by `in_la_county`); 16% of deaths at addresses
shared by ≥5 others — hospitals and jails, where people are *pronounced*
(reported as `pct_at_facility`; amphetamine is 44% and should be discounted);
3 Catalina deaths stretching the map frame.

*Guard added after review:* methamphetamine (61% of the cohort) and fentanyl
(58%) are **withheld, not scored**. Under random labelling the null draw
becomes the case set itself as n/N → 1, so the test compares a substance
against itself and returns ≈1 regardless of geography. `MAX_CASE_FRACTION =
0.25`.

*Reasoning corrected mid-review:* including a substance's own deaths in the
resampling pool is **required** by the random-labelling null, not
contamination as first described. The practical conclusion was unchanged; the
justification was wrong.

---

## Earlier — build-out

**Detection statistic.** Two-component MGPS implemented first and abandoned:
not identifiable at 205 substances / 111 with non-zero baselines. Nelder-Mead
converged to different degenerate solutions depending on the fit set —
Gamma(682, 627) with everything in, a point mass at 1.08 with large substances
excluded. Single-component GPS is stable and has closed-form posteriors.

**Backtest established the alarm line.** All five known LA emergences reach
rank 1–4 in the quarter they emerge, **including xylazine at n=6**. Real
emergences score EB05 1.93–7.00, which is what makes 1.5 defensible rather than
arbitrary.

*Correction:* an earlier claim that only para-fluorofentanyl was well enough
powered to support a detection curve was **wrong** — the backtest disproved it.

**Extraction defects found and fixed** (99.4% → 99.8% coverage):
- The hyphen rule was all-or-nothing, discarding known substances whenever a
  sibling token was a modifier (`METHAMPHETAMINE-INDUCED`).
- `ALCOHOLIC INTOXICATION` — the adjective form was the single largest gap, 67
  deaths.
- Brand names and non-catalog poisons (Xanax, arsenic, oleandrin): NFLIS
  catalogues trafficked *drugs*, so poisons are simply absent from it.

**Do not fix the remainder by loosening the fuzzy cutoff.** At 0.85 it admits
six false positives, worst being bare `ALCOHOLIC` → Ethanol at 104 mentions,
which would tag alcoholic cirrhosis as acute alcohol involvement. Those tokens
are in `FUZZY_BLOCK` so the mistake cannot be reintroduced by changing one
number.

**Two artifacts that would have manufactured findings:**

1. **Fluorofentanyl isomer split.** `PARA-FLUOROFENTANYL` and bare
   `FLUOROFENTANYL` ran as separate series; no case names both; the ME switched
   convention in 2023Q2. Unrolled, that reads as one substance collapsing and
   another emerging. Fixed by `ISOMER_ROLLUP`. **Assume more of these exist.**
2. **A 4-day 2026Q2 stub** holding 3 deaths was silently consuming one of two
   `lag_quarters` slots, so a 2-quarter correction discarded only one real
   quarter. Fixed by dropping calendar-incomplete quarters separately from the
   settledness cutoff.

**One extraction fix changed a result:** `CARENTANYL EFFECTS` (2017-06) is a
real carfentanil death, moving carfentanil's first LA appearance from 2024Q1 to
2017Q2 — an isolated case during the national 2016–17 wave, then a seven-year
gap before the 2024 cluster.

**Windowing is not innocuous.** Lidocaine scores EB05 1.04 on 2025Q1–Q4, 2.20
on 2024Q4–2025Q3, and peaks at 4.29 in the sweep. For spiky trajectories,
quote the backtest peak alongside the current value.

**Methods rejected at design time**, with reasons that still hold:
SaTScan/sparr as *detectors* (need n ≥ 50–100); substance co-occurrence graph
as a detector (two-hub star — meth 54%, fentanyl 54%, everything else < 6%;
retained as an interpretive layer in `polysubstance.py`); zero-inflated Poisson
and hurdle models (counts run 0–5); occupancy-detection models to separate
"absent" from "not assayed" (needs detection probability to vary across units,
and there is one ME office running one evolving panel).
