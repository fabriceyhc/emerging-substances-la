# GPS v2: a design note

Status: **all three built, plus one extension to §1, plus a TreeScan
cross-check -- and this combination is now `emerging trends rank`'s
default.** §1 behind `--weighted` (and `--role-discount` on top of it), §2 as
the `eb05_own`/`credible_rise` columns, §3 behind `--spatial discount|prior`
with the score itself in `emerging/core/concentration.py`. A dozen-plus
variants are scored head to head, against a reference that has grown from 5
to 22 labelled substances (including 7 documented negatives, with a Wilson
confidence interval on precision) in `emerging/validation/benchmark.py`
([findings/benchmark.md](findings/benchmark.md)) — read that doc for the
current verdict, since it has been revised more than once as the reference
grew and should not be summarized here where it will go stale again. **The
production default** (`rank`/`backtest`'s own defaults as of this revision):
`--weighted --role-discount --spatial discount --treescan-veto`, i.e.
`eb05-v2-role` (§1 extended + §2 + §3) with `emerging.analysis.treescan`'s
substance-leaf scan able to veto (not add) an alarm it sees no corroborating
rise for at all (`docs/findings/benchmark.md`'s `ensemble-veto-v2role-10`) --
every EB05-family alarm the veto can suppress lives in
`emerging.analysis.trends._apply_treescan_veto`, not in the comparison
framework, so it's live in `rank`'s output, not only in a backtest table.
This is a point on a recall/precision frontier, not a proven optimum (see
that doc's Verdict) -- every piece is independently switchable back off with
its own flag, including to plain single-gate EB05
(`--no-weighted --no-role-discount --spatial none --no-treescan-veto`).
Section by section: §2 validates cleanly against `backtest`. §1 failed its
first backtest gate (the xylazine finding), was fixed with
`_ever_independent`, and now passes at the single-gate level — see its
Validation result section for what's still marginal at very low n. The
failure and the fix are both kept in the record rather than cleaned up after
the fact, per this project's convention of keeping findings that turn out to
be artifacts. Prior context in [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md)
and [RESEARCH_LOG.md](RESEARCH_LOG.md); this note predates the repo
reorganisation into `emerging/analysis/`, so earlier path references below
may say `trends.py` where the file is now `emerging/analysis/trends.py`.

## Motivation

`trends.py`'s GPS/EB05 ranking has three known gaps, each currently patched by
a *separate* downstream module rather than fixed in the statistic itself:

1. It counts every mention equally, so an adulterant that always rides along
   with the real driver (lidocaine on fentanyl) ranks as a credible riser.
   Caught after the fact by `polysubstance.py`'s cause-line position analysis.
2. It measures *share* of all OD deaths, so a shrinking denominator inflates
   every flat-count substance's score (cocaine: EB05 1.06, `count_ratio`
   0.90). Caught after the fact by `regime()`'s diagnostics.
3. It has no spatial term at all, so a temporally marginal but tightly
   localized signal — plausibly a single contaminated batch or dealer — gets
   no extra weight over a diffuse, county-wide one of the same size. Caught,
   if at all, by manually cross-referencing `geo.py`/`spacetime.py`.

This note designs three extensions that build the fixes into the ranking
statistic itself, in the same closed-form GPS family, rather than adding a
fourth downstream cross-reference. The precedent for combining detectors
formally rather than reading them separately by eye is the existing
EB05×TreeScan 2×2 table in `RESEARCH_LOG.md` ("How EB05 and the tree scan
divide the work") — this generalizes that pattern.

---

## 1. Position-weighted case definition

### Problem with the naive fix

A hard restriction to "count only when the substance is lead on the cause
line" throws away a real category: a substance that is a necessary
contributor to a lethal combination (e.g. a benzodiazepine potentiating
opioid toxicity) but is conventionally listed second, never first. That
substance should still register a spike in weighted count if its rate rises —
restricting to lead-only would zero it out entirely.

### Design: graded integer credit, not a binary filter

For death *i* naming substance *s* on a cause line of length *L_i*:

| Position | Credit `c_i` |
|---|---|
| solo mention (`L_i == 1`) | 2 |
| named, not last (`L_i > 1`) | 1 |
| named last, `L_i >= 3` | 0 |
| named last, `L_i == 2` | 1 *(see below)* |

The `L_i == 2` carve-out matters for the same reason `polysubstance.py`
already flags: on a two-substance line, "last" and "not first" are the same
fact, and a substance's line-length distribution is itself confounded with
its role (adulterants co-occur with more substances, so sit on longer lines).
Zeroing "last on a 2-line" would reproduce that bias. Only "last on a line of
3+" is treated as the adulterant signature, matching the empirical pattern
that grounded this design: lidocaine is last in 82% of (mostly 3+) lines and
first in 0% — see [`docs/findings/polysubstance.md`](findings/polysubstance.md).

`n_recent` / `n_baseline` / `expected` in `rank_substances` (`trends.py:265`)
are recomputed as sums of `c_i` instead of raw mention counts. Nothing else
in `_fit_gps_prior` or `_gps_scores` changes.

### Why integers, not a continuous weight — and why that turned out to be the
### wrong axis to worry about

GPS's closed form depends on `n ~ Poisson(λE)` being exactly true — that is
what makes the negative-binomial marginal, and therefore the closed-form
`Gamma(α+n, β+E)` posterior, exact. A continuous-weighted sum is not
Poisson-distributed in general (its variance does not equal its mean once
weights vary), so feeding floats into `_nb_logpmf` would fit a
negative-binomial likelihood to data that violates the assumption the
derivation depends on. Options considered, most to least expensive:

- **Tweedie / compound Poisson-Gamma.** The statistically "correct" model if
  credit is treated as a continuous mark on a Poisson arrival process (the
  standard actuarial model for Poisson-count-of-claims × Gamma-ish claim
  size). No closed-form empirical-Bayes shrinkage as clean as GPS's — loses
  the property that makes the 45-quarter backtest sweep cost ~1s.
- **Quasi-Poisson / weighted NB GLM.** Valid point estimates off exact
  Poisson, but the closed-form posterior percentile that EB05 *is* does not
  survive; recovering an EB-flavored interval needs a numeric hierarchical
  GLMM.
- **Small-integer credit (originally chosen, and wrong on its own reasoning).**
  The draft argued the weighted sum staying a literal nonnegative integer was
  enough to keep it "approximately Poisson to the same degree the existing
  raw-count model already is." **That claim doesn't survive contact with the
  actual scheme.** Integer-valued is not the property that matters — dispersion
  is. For any per-mention credit `c` with `SOLO_CREDIT = 2` in play, `c² ≥ c`
  pointwise, so `E[c²] ≥ E[c]` and the weighted sum is *structurally
  overdispersed* relative to a Poisson with the same mean, integer-valued or
  not. Caught on this repository's real data: feeding the raw {0,1,2}-weighted
  counts straight into `_fit_gps_prior` inflated the false-`EB05>1` rate from
  5.1/203 substances (plain counts) to 19.0/203 — the scheme was finding
  spurious signal by construction, not by design.

  The actual dispersion-safe property is `φ = E[c²]/E[c] = 1`, which holds
  **identically**, not approximately, for any `{0,1}`-valued credit — `c² ≡ c`
  pointwise for a binary variable, so this holds regardless of how the
  keep/discard split varies case to case. That is a materially stronger
  argument than "approximately Poisson," and it is the reason a `{0,1}` scheme
  (drop `SOLO_CREDIT`, keep only the trailing-adulterant discount) would need
  *zero* changes to `_fit_gps_prior` / `_gps_scores`. It was kept as `{0,1,2}`
  anyway — a solo mention is arguably the strongest available evidence of
  causal weight, a different claim than the independence question
  `polysubstance.alone_pct` already answers — and the dispersion is corrected
  explicitly instead, below.

### The dispersion correction actually shipped

`φ = E[c²]/E[c]`, computed once over the full weighted-mention population in
`load_quarterly` (`_credit_dispersion`, one line of pandas) and carried on
`counts.attrs["phi"]` (1.0 for the unweighted path, so nothing about existing
callers changes). `rank_substances` divides both `n` and every `expected`
array by `φ` before fitting or scoring, for both Gate A and Gate B (§2):
`λ = n/E` is unchanged, since both are scaled by the same factor, but the
posterior's effective sample size shrinks — precisely the right response,
since a compound-Poisson sum `W = Σc_i` has `Var(W) = μ·E[c²]` against
`Mean(W) = μ·E[c]`, so `φ` is exactly the excess variance-to-mean ratio the
weighted sum carries beyond what a same-mean Poisson would. This is the
standard quasi-likelihood move for overdispersed count data; `gammaln` is
defined on the reals, so the now-fractional `n/φ` runs through `_nb_logpmf`
mechanically fine, at the cost of its output being a quasi-likelihood rather
than a normalized pmf. On this repository's data: `φ ≈ 1.42`, and the
corrected false-`EB05>1` rate returns to 4.9/203 — back in line with plain
counts' 5.1/203.

### Validation result: §1 fails its own gate

The original validation plan was to check `backtest` on weighted counts and
confirm all five known emergences still reach rank 1-4. It was checked, and it
fails for a real reason, not a calibration nit:

**Xylazine — one of the five known emergences — regresses badly under
`--weighted`.** Peak EB05 drops from 1.93 to 1.05, best rank from 4 to 12, and
it never reaches top-10 (unweighted: top-10 at 2024Q2). Cause, confirmed
directly against the credit table: of xylazine's 9 lifetime mentions, 3 are
`TRAILING_CREDIT`-discounted to zero — a third of an already n=9 substance's
total evidence, discarded.

The reason is not a calibration bug, it is a **construct-validity problem**.
Xylazine is a fentanyl adulterant — the `WATCHLIST` entry in `trends.py`
describes it exactly that way — so it is *frequently coded last on the cause
line alongside fentanyl*, the same surface pattern lidocaine has. But the two
substances are not the same kind of "adulterant" for this project's purpose:
lidocaine is pharmacologically inert and the right headline is "not a killer";
xylazine is exactly the kind of emerging contaminant this whole pipeline
exists to flag, regardless of whether the ME's narrative treats it as the
dominant cause of a given death. Cause-line position measures "was this
compound clinically dominant in this specific death," which is a different
construct from "is this compound entering the supply and worth tracking" —
and §1's design conflates them. Bromazolam (peak EB05 4.78 → 1.75) and
Carfentanil (5.02 → 2.24) show the same direction of effect, more mildly.

**§1 is not ready to be a default.** It stays available behind `--weighted`
for further work, but should not be treated as a validated improvement over
plain EB05 until this is resolved — candidates worth discussing: exempting
known-adulterant-pattern substances (`WATCHLIST`) from the trailing discount,
using a milder discount than `TRAILING_CREDIT = 0`, or accepting that position
on the cause line is the wrong signal for this half of the problem entirely
and this design should be shelved in favour of something that reads
adulterant-vs-driver status from `polysubstance.py`'s output as a *label*
rather than folding it into the case definition.

### The fix: `_ever_independent` — read the substance's full history, not
### just each death's own line

Chosen over the alternatives above. For substance *s*, check whether it has
ever been named *alone* or *first* anywhere in the cohort up to the point
being scored. If so, it never receives `TRAILING_CREDIT`, even in a case
where it is named last on a 3+-line — other evidence already establishes it
can be the whole story. Solo and lead collapse to one check: the first
element of a case's ordered sequence is the lead substance if the line has
more than one entry, and the substance itself if it doesn't — so
"named alone or first" is exactly "appears at index 0 of some case's ordered
sequence," `_ever_independent` in `trends.py`.

**Clean separation on real data.** Lidocaine and levamisole (the two known
inert cutting agents) are named alone or first **zero** times across their
entire histories (42 and 12 mentions respectively) — not a close call.
Xylazine clears it once, in 2024Q2, the quarter immediately before its
backtest emergence is detected. Carfentanil and PCP clear it dozens of times;
they were never really at risk from this problem.

**This is a cross-case aggregate, and that changes the as-of story.**
Per-death credit only ever reads that one death's own cause-line text, so
slicing the resulting counts by quarter after the fact — what `backtest` does
to the unweighted path today — cannot leak the future backward. `_ever_
independent` is different: it looks across *every* case up to some boundary,
so computing it once from the full `cutoff`-bounded cohort and reusing it for
every earlier as-of quarter would let a substance's *later* proof of
independence retroactively un-discount its *earlier* trailing mentions —
using information from after the point being scored. Fixed by making
`load_quarterly`'s weighted path bound `_ever_independent`'s input to its own
`cutoff` (previously it saw the whole unfiltered extract regardless of
`cutoff`, a latent gap that didn't matter for per-death credit but would have
for this), and by having `backtest --weighted` call `load_quarterly` fresh at
every as-of step instead of slicing one precomputed table. Cost: `causea_order`
re-runs 45 times instead of once — a full `backtest --weighted` run takes
~40s, still cheap enough to be a normal part of iterating on this.

**Result.** Xylazine: peak EB05 1.05 → 1.46, best rank 12 → 5, and it reaches
top-10 again (2024Q3, one quarter later than plain EB05's 2024Q2, versus never
under the unfixed `--weighted`). Bromazolam: 1.75 → 2.69. Carfentanil:
2.24 → 2.83. All still short of the unweighted numbers, which is expected —
`--weighted` is a stricter case definition by design — but no longer a
regression against the *unweighted* backtest's rank/top-10 outcome.

**What's still unresolved: xylazine never clears the §2 dual gate under
`--weighted`.** At its peak quarter (2024Q3, n_recent=7, n_baseline=1), EB05
is 1.46 and `eb05_own` is 1.32 — both just under the 1.5 line. This is the
signal-ceiling problem `RESEARCH_LOG.md` already documents (three of the five
known emergences sit at n < 10, where Poisson noise is what's binding, not the
estimator) rather than a defect in `_ever_independent` itself: requiring *two*
independent 1.5-crossings at n=7 is a materially harder bar than one, and
xylazine was already a marginal case for §2 alone even on unweighted counts.
Worth knowing about before reading `credible_rise` as the final word on a
substance this rare.

### `--role-discount`: the Lidocaine residual, and a real bug caught before
### it shipped

`_ever_independent` fixed xylazine but left a scoped gap: Lidocaine still
alarmed once under `--weighted` (down from 4 quarters to 1, not 0). Checked
directly — of 34 recent Lidocaine mentions, 12 are literally
`[Fentanyl, Lidocaine]`, a plain 2-item line, which keeps full per-death
credit by design (the 2-line carve-out that stops a genuinely short-lined
substance from being punished for its own line length). That's real,
undiscountable signal at the per-death level.

`--role-discount` adds `_role_dispersion`: a per-substance, as-of-safe
`phi_s` computed from the substance's own *aggregate* mean credit, divided
into `n` and `E` the same way the global `phi` correction already is — a
constant per-substance factor that leaves the raw ratio `n/E` unchanged
(checked numerically) but pulls a substance's *posterior* harder toward the
population prior, since `alpha`/`beta` are fit once, globally, and don't
rescale with one substance's own discount.

**The first version of this used `SOLO_CREDIT` (2) as the reference point
and broke every known emergence in `backtest`** — Bromazolam's peak EB05
dropped from 2.69 to 1.00, Xylazine's dual gate stopped clearing at all.
Cause: `_ever_independent`-exempted substances never receive
`TRAILING_CREDIT`, but *solo* is rare for nearly any real drug — checked
directly, Bromazolam's own mean credit is 1.000, Carfentanil's 1.167,
PCP's 1.001, Fentanyl's 1.223, none within reach of `SOLO_CREDIT`. Using it
as the reference punished every one of them for not being frequently solo,
a bar almost nothing clears. Fixed by using `NAMED_CREDIT` (1) as the
reference instead — "was this mention ever relegated to trailing" — with a
`max(1, ...)` floor so a substance sitting above 1 (Fentanyl) never gets
*amplified*. Re-validated: known emergences are within a few hundredths of
their `--weighted`-only EB05 (Bromazolam 2.69→2.55, Carfentanil 2.83→2.59),
and **Lidocaine's alarm is fully suppressed — 0 quarters, down from 1 —
across its entire real alarm history (2024Q4–2025Q3), not just the current
window.** The bug and the fix are both kept in the record, the same
convention as the xylazine finding above: caught by running `backtest`,
not by inspection.

Scored properly as `eb05-weighted-role` in `emerging/validation/benchmark.py`
— see `docs/findings/benchmark.md` for the full result, including the
trade-off it makes at a matched alarm budget (best precision in the whole
table at a fixed threshold; a lower re-cut threshold at matched budget lets
in some weaker alarms, since a role-discounted substance needs more raw
evidence to clear any given bar).

---

## 2. Joint dual-gate for denominator drift

**Implemented and validated.** All five known emergences still clear the
joint gate in `backtest`, at a cost of at most one extra quarter of lag versus
plain top-10 (para-Fluorofentanyl, Bromazolam, Carfentanil: +1 quarter;
Mitragynine: same quarter; Xylazine: +1 quarter, still clears — contrast with
§1, where it never clears once `--weighted` is added). Unlike §1, this holds
up on its own.

### Design

Two GPS fits, differing only in `expected`, both reusing `rank_substances`
unchanged:

- **Gate A — share (current).**
  `E_s = n_base * (n_recent_all / n_base_all)` — is this substance's share of
  all OD deaths rising.
- **Gate B — self-referential.**
  `E_s = n_base * (recent_quarters / baseline_quarters)` — no county-total
  term. Is this substance's own count above where its own baseline rate would
  put it. (This is exactly `count_ratio`'s denominator, `trends.py:313`, run
  back through GPS/shrinkage instead of reported as a raw ratio.)

### Decision rule

|  | EB05^own > 1.5 | EB05^own ≤ 1.5 |
|---|---|---|
| **EB05^share > 1.5** | credible rise | share-only — denominator-driven (cocaine's pattern) |
| **EB05^share ≤ 1.5** | count-only — rare; total deaths rising faster elsewhere | no signal |

"Rising" requires both cells to clear the line. This formalizes what
`regime()` currently does as a manual post-hoc diagnostic into a labeled cell
in the ranking output itself.

### Implementation cost

Minimal — one more call to the existing `rank_substances` with a different
`expected` formula, plus a join on the two result frames. No new fitting
machinery.

---

## 3. Spatial concentration in the core formula

**Built, and it settles its own open question in favour of the covariate
prior — for a reason that is not "it scored better".** Both fusions below were
implemented; the score `z_s` they consume is
`emerging/core/concentration.py`, written up in
[findings/concentration.md](findings/concentration.md), and the head-to-head
is [findings/benchmark.md](findings/benchmark.md). Result summary before the
design text, which is left as drafted:

- **The `z_s` this section asked for did not exist in any usable form.** The
  section correctly rules out `geo.py`'s lifetime ratio and points at
  `spacetime.py`'s trailing windows — but `spacetime.py` costs 999
  permutations per substance and refuses to run below `MIN_CASES = 20`, so it
  cannot feed a 45-quarter sweep and is undefined for four of the five
  labelled emergences. What shipped instead is a Cuzick–Edwards-style
  case-control clustering count with **closed-form null moments** in place of
  a permutation p-value: a 45-quarter × 211-substance sweep in 4.6s. The open
  question below ("raw LLR, its RI-derived p-value, or the cylinder relative
  risk") is answered by not using the cylinder scan at all.
- **The fallback (`E`-discount) fails the matched-budget test and the
  covariate prior passes it.** At the fixed 1.5 line the discount posts the
  best recall in the whole benchmark (0.44 vs 0.38). At a matched alarm budget
  it drops to 0.40, and the substance it appeared to gain — xylazine — has
  `z_s` identically zero in all 45 quarters, so the discount could not have
  moved it directly. It moved the threshold. That is precisely the "silently
  discards the conservative-bound property" hazard this section warns about
  under multiplicative fusion, arriving through a different door.
- **The covariate prior's virtue is that it can decline.** Fitted γ at the
  current quarter is **+0.061**, effectively zero, and the model is within
  0.01 of plain EB05 on both readings. Because γ is *estimated*, a covariate
  carrying no information reduces the fit exactly to the unadjusted one. `w`
  cannot do that. The draft below calls the discount "less principled because
  `w` is chosen, not estimated"; the benchmark shows that is not an aesthetic
  objection, it is the mechanism of the failure.
- **Where it does act, it acts correctly, and the limit is power not design.**
  para-Fluorofentanyl scores z = 3.18 in 2021Q2 — the quarter after its first
  LA death — decaying to −1.83 by 2022Q3 as it saturates the county-wide
  supply. Mitragynine scores 4.63 at its onset, carfentanil 3.47 at its 2024
  cluster. Those are the two-and-a-bit substances with ≥10 deaths in a
  trailing window; xylazine, at 9 lifetime mentions, is structurally invisible
  to any spatial statistic and `MIN_CASES` says so rather than returning a
  number one coincidental neighbour pair would have produced.

### Rejected: multiplicative fusion

Ranking by `EBGM * spatial_RR` (or similar) is tempting but statistically
incoherent: multiplying a shrunk lower-percentile estimate by an
independently-computed relative risk does not correspond to a percentile of
anything, so it silently discards the "conservative bound" property that is
the reason to use EB05 over a raw ratio in the first place.

### Primary direction: covariate-adjusted prior

`_fit_gps_prior` (`trends.py:198`) currently fits one global `Gamma(α, β)`
shared across all substances. Extend the rate parameter to depend on a
spatial-concentration covariate:

    beta_s = beta_0 * exp(-gamma * z_s)

where `z_s` is a standardized concentration score and `gamma` is one
additional parameter fit by the same MLE that already estimates `alpha` and
`beta_0`. A higher `z_s` lowers the effective shrinkage target for that
substance — a spatially concentrated substance needs less raw temporal
evidence to clear the same EB05 bar. This changes the shrinkage target
itself, rather than adjusting the output score after the fact.

Closest available precedent: DuMouchel & Pregibon's covariate-adjusted
extensions to GPS for FDA AERS screening. **Flagged as needing citation
verification** — the exact mechanism has not been checked against the
original source and should not be cited in a manuscript until it is.

### Simpler fallback: discount the expected count

    E_s_spatial = E_s * (1 - w * concentration_s)

for a tunable `w`. No new fitting step — one line in `rank_substances`. Less
principled than the covariate-prior approach (`w` is chosen, not estimated
from the data), but useful as a first prototype to sanity-check whether the
covariate-prior approach's added complexity is worth it.

### Which spatial signal to use — this matters

`geo.py`'s localization ratio is whole-history and, per `RESEARCH_LOG.md`,
"the geography is real and it is static" (PCP's clustering reflects a stable
local market, not a recent event). The point-source-contamination hypothesis
this design is meant to catch — one dealer spiking a batch — is about a *new*
cluster in a recent window, which is what `spacetime.py`'s trailing cylinders
test, not what `geo.py`'s fixed-window test does. `z_s` / `concentration_s`
should be built from the **trailing-window space-time signal**
(`spacetime.py`), not the lifetime `geo.py` ratio — using the lifetime metric
would let a substance's stable historical geography (e.g. PCP in Skid Row)
lower the bar for an unrelated, unrolling event.

### Keep the decision-layer column too

The covariate-prior fusion and a plain third decision-layer column (space-time
cylinder significant / not, alongside the Gate A / Gate B table from §2) are
not redundant. Composition and geography are largely orthogonal in this data
already — flagged novel substances mostly "track the county-wide fentanyl
supply" rather than localize, per `RESEARCH_LOG.md`. The prior fusion answers
"how suspicious should a given temporal signal be"; the column answers "is
there also an independently significant space-time cluster right now."
Reporting both lets each fire independently, which is itself informative
rather than double-counting.

---

## Prototyping order

1. ~~Integer-credit position weighting (§1)~~ — built, dispersion-corrected,
   **failed its first backtest gate** on xylazine, then fixed with
   `_ever_independent` (read the substance's full history, not just each
   death's own line) and **re-passes at the single-gate level**. `eb05_own`
   dual-gate still doesn't clear for xylazine specifically — a low-n signal
   ceiling, not a defect in the fix. Behind `--weighted`; still not promoted
   to a default, but no longer blocked on a known regression.
2. ~~Dual gate (§2)~~ — built and **validated**: all five known emergences
   still clear it, at a cost of at most one extra quarter of lag. Exposed as
   `eb05_own` / `credible_rise` in `rank` and `backtest`.
3. ~~Spatial fallback (§3, `E_s` discount) first, covariate-prior second~~ —
   both built, on a `z_s` (`core/concentration.py`) that had to be designed
   from scratch because neither `geo.py` nor `spacetime.py` could supply one.
   The ordering advice turned out to be backwards: the fallback is the one
   that *fails*, by inflating the whole EB05 distribution rather than
   reordering it, and the added fitting complexity of the covariate prior is
   what buys the ability to return "no effect" — which is what it returns
   (γ = +0.061). Neither improves on plain EB05 at a matched alarm budget.

Each is validated against the existing `backtest` set before being read as an
improvement — a change that moves a known emergence to a *later* detection
quarter than plain EB05, or stops it from clearing at all, is a regression,
not a refinement. That is exactly what happened to §1 on the first cut, and
exactly what stopped happening once `_ever_independent` was added.

## Open questions

- **Is a single solo-or-lead occurrence enough evidence, or too fragile?**
  `_ever_independent` flips on exactly one case for xylazine. The separation
  from lidocaine/levamisole (0 occurrences across 42 and 12 mentions) is clean
  enough on today's data that this wasn't tightened, but a substance whose
  *only* solo/lead case is a data-entry idiosyncrasy would be permanently
  exempted from the discount on that single point. Worth a minimum-occurrence
  or rate-based version if a future case looks fragile the way xylazine's
  clean split doesn't.
- Does the position-credit scheme need per-substance-class tuning (opioids vs.
  adulterants vs. stimulants), or is `_ever_independent` enough separation on
  its own now?
- ~~What is a defensible source for `z_s` beyond the space-time cylinder's
  LLR — raw LLR, its RI-derived p-value, or the cylinder relative risk
  itself?~~ **Answered: none of the three.** The cylinder scan cannot be run
  45 times per substance and is undefined below n = 20. `z_s` is a
  Cuzick–Edwards case-control clustering count with exact null moments
  (`core/concentration.py`), which needs no permutation at all. Two
  alternatives measured and rejected on the way — zipcode cells (no resolution
  below n ≈ 50) and a fixed radius (power goes to "is this downtown" rather
  than to concentration) — are recorded in
  [findings/concentration.md](findings/concentration.md).
- **New, from the benchmark:** the five-substance reference set cannot resolve
  differences below ~0.05 mean recall, which is most of the spread between the
  six shrinkage models. Growing the reference set is now the binding
  constraint on comparing detectors at all — a better statistic cannot be
  demonstrated on this many labels.
- Should Gate B's baseline be a fixed two-window ratio (as drafted) or a
  smoothed secular trend (borrowing `nowcast.py`'s regression machinery),
  the way `svtt.py` already fits trends rather than window ratios?
- Verify the DuMouchel & Pregibon citation before it appears anywhere public.
