# Ground-truth interval overlap

Findings for `emerging/validation/ground_truth.py`. Reference data in
[`data/raw/ground_truth/known_emergences.csv`](../../data/raw/ground_truth/known_emergences.csv),
results in [`results/ground_truth/`](../../results/ground_truth/). Back to the
[README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## Why `backtest` alone wasn't enough

`backtest` answers "did the ranking notice this substance, and how did its
rank evolve" — it has no independent record of *when* an emergence actually
started or ended, so it can't say whether an alarm episode was well-timed,
early, late, or ran on too long. This adds that reference: one row per known
local emergence interval, scored against the EB05 alarm sweep's detected
episodes with a temporal-IoU-style metric (`interval_overlap`).

**Avoiding circularity was the design constraint.** A reference derived from
EB05 itself would validate the detector against its own output. Two
provenances instead: `raw-count-rule` — a fixed, stated rule applied to the
substance's *raw* quarterly mention count, independent of any fitted
statistic — and `manual`, for the two substances (carfentanil, xylazine) too
sparse and volatile at low n to ever satisfy that rule, annotated by eye with
the reasoning written into the CSV's `notes` column.

## Definitions

Everything downstream reduces two sets of quarters to six numbers.
**Ground-truth quarters** are the reference interval(s) for a substance from
`known_emergences.csv` (unioned if it has more than one, non-contiguous
episode). **Detected quarters** are the EB05 alarm sweep's episode(s) for
that substance (`trends._episode_spans`, also unioned if more than one).
Everything else is built from comparing those two sets.

| Metric | Formula | Reads as |
|---|---|---|
| `gt_quarters` | `\|ground-truth quarters\|` | how long the real emergence lasted |
| `det_quarters` | `\|detected quarters\|` | how long EB05 stayed in alarm |
| `overlap_quarters` | `\|ground-truth ∩ detected\|` | quarters both sides agree on — the number every other metric is built from |
| **precision** | `overlap / det_quarters` | of the quarters EB05 alarmed, what fraction were a real emergence — the **false-alarm** read: low precision means EB05 was in alarm when it shouldn't have been |
| **recall** | `overlap / gt_quarters` | of the quarters that were a real emergence, what fraction EB05 caught — the **coverage** read: low recall means EB05 missed part of a real emergence |
| **IoU** (Intersection over Union) | `overlap / \|ground-truth ∪ detected\|` | one number capturing both directions of mismatch at once — an alarm that runs short *and* one that runs long both pull IoU down, which is why it's reported alongside precision/recall rather than instead of them |

**Worked example, Bromazolam:** ground truth 2023Q3–2025Q1 (7 quarters),
detected 2024Q1–2025Q1 (5 quarters), and the detected span sits entirely
inside the ground-truth one, so overlap = 5, union = 7.
`precision = 5/5 = 1.00` (every alarmed quarter was real), `recall = 5/7 ≈
0.71` (the alarm covered 5 of the 7 real quarters — it started late),
`IoU = 5/7 ≈ 0.71`. Reading precision and recall together tells you *which
way* the mismatch runs; IoU alone would only tell you that one exists.

**Two ways a metric comes out undefined (`NaN`), not zero, and why that
distinction matters.** If a substance never alarmed at all, `det_quarters =
0` and precision is `0/0` — undefined, not 0, because there's nothing to be
imprecise *about*. If a substance is a documented `no_emergence` case
(Morphine), `gt_quarters = 0` and recall is `0/0` — undefined, because there
was no real emergence to have covered. Both are reported as `NaN` and
excluded from the `mean_*` summary figures (`out.dropna(subset=["iou"])` in
`score`), rather than silently counted as 0 and dragging the average down
for a reason that has nothing to do with detector quality.

**`mean_precision` is an average of ratios, one per substance — it is not,
by itself, a rate with a meaningful confidence interval.** `precision_rate`
(`emerging.validation.ground_truth.precision_rate`) is a different
quantity, built to have one: one Bernoulli trial per *alarmed* substance
(`det_quarters > 0`), "success" if any part of the alarm overlapped real
ground truth (`overlap_quarters > 0`, so Benzylfentanyl's partial 0.33
precision still counts as a success — this is answering "did the alarm
ever touch something real," a coarser question than `mean_precision`'s "how
much of the alarm was real"). Each substance contributes exactly one trial,
regardless of how many quarters it was alarmed for, specifically so that
Lidocaine's 4 correlated quarters of one sustained mistake don't get
counted as 4 independent data points — that would inflate the apparent
sample size and produce a falsely narrow confidence interval. The interval
itself is a **Wilson score interval** (`wilson_ci`), preferred over the
normal approximation at the single-digit-to-low-double-digit n this file
actually has, where the normal interval can extend below 0 or above 1.

`emerging/validation/benchmark.py` adds one more metric on top of these,
specific to comparing several detector variants: **`mean_lag_q`**, the
number of quarters from a ground-truth interval's start to the first alarm
*inside* that interval (not the first alarm anywhere — an earlier spurious
alarm elsewhere shouldn't score as prescient). A substance the model never
alarms on inside its interval contributes no lag value and drops out of the
mean, so `mean_lag_q` has to be read alongside `detected`, the count of
substances that got any qualifying alarm at all — see
`docs/findings/benchmark.md`.

## Result, original five known emergences

| Substance | IoU | Recall | Precision | Ground truth | Detected |
|---|---|---|---|---|---|
| para-Fluorofentanyl | 0.40 | 0.40 | 1.00 | 2021Q1–2025Q4 | 2021Q2–2023Q1 |
| Bromazolam | 0.71 | 0.71 | 1.00 | 2023Q3–2025Q1 | 2024Q1–2025Q1 |
| Mitragynine | 0.13 | 0.13 | 1.00 | 2022Q2–2025Q4 | 2022Q4–2023Q1 |
| Carfentanil | 0.50 | 0.50 | 1.00 | 2024Q1–2025Q4 | 2024Q3–2025Q2 |
| Xylazine | 0.14 | 0.14 | 1.00 | 2023Q3–2025Q1 | 2024Q3–2024Q3 |

**Precision is 1.00 across all five.** That looked reassuring, and it was
the wrong number to trust — see below.

## The 1.00 precision figure was measured on a selection-biased sample

24 substances have ever breached EB05 > 1.5 (`alarm_history.csv`). Only 5
were labelled. **Precision of 1.00 measured on 5 substances we already knew
were real emergencies says nothing about whether EB05 ever cries wolf** —
there was no case in the file where it *could* have been wrong. Expanded
the labelled set to cover 19 of the 24 (excluding lidocaine and other
known-inert adulterants — `polysubstance.py` already established those
aren't real emergencies, and re-litigating that isn't what this file is
for), using the same non-circular methodology throughout: raw quarterly
LACME mention counts, never EB05 or any fitted statistic. Not presence in
the drug supply either — that's a different, laggily-related construct
(see the NFLIS section above); the ground truth for a *death*-based
detector has to be built from *death* data.

**Two shapes of "real emergence," not one rule.** The original raw-count-rule
(arrival from near-absence — `≥2 consecutive quarters ≥2` after `≥2
consecutive quarters ≤1`) assumes the substance is new. Four of the
candidates aren't: Fentanyl, PCP, Codeine and Morphine are censored at
`first_seen = 2012Q1` — present since before the extract even starts, so
there's no arrival to find. For those, onset/end were read directly off the
share trajectory instead (same principle, different shape: where does the
*rate* credibly step up or down).

**Seven substances are documented negatives, not positives.** Each is a
different shape of "EB05 alarmed, but this isn't a real local emergence,"
which matters more than the count: a false positive can happen for
genuinely different reasons, and conflating them into one number would hide
which mechanism is responsible.

- **Morphine, Mirtazapine** — *spurious blips*. Morphine: share fell from
  18–27% (2014–2015) to 1–2% (2022–2025), a clean multi-year structural
  decline; the one breach (peak 1.53, 2025Q2) traces to a single-quarter
  raw-count blip (2024Q3 = 20 against a surrounding range of 2–12).
  Mirtazapine: same mechanism on much sparser data (n=40 lifetime) — one
  standout quarter (2023Q4, 0.75% share) against a background of near-zero.
- **Hydromorphone** — a *real decline*, Morphine-shape but without the blip:
  share falls cleanly from ~0.6–1.6% (2018) to ~0.1–0.4% (2019–2023), no
  reversal worth calling real.
- **Lidocaine, Levamisole** — *real signal, wrong interpretation*.
  Lidocaine's rise is real, not noise (1 → 10 → 11 mentions,
  2024Q3–2025Q1) — `polysubstance.py` already established why it still
  isn't a threat: named last on 82% of its cause lines, first on 0%, 100%
  co-occurring with fentanyl. Levamisole is the same adulterant pattern,
  thinner data (n=12).
- **Citalopram** — *no clean story, lowest confidence of the seven*. Share
  oscillates irregularly (0.0–1.2%) with no sustained run anywhere and a
  late partial resurgence that muddies a simple decline narrative. Included
  because there's no real emergence to find either, not because the
  negative case is as clean as the others — revisit if more data changes
  this read.
- **Alprazolam** — *no real change at all*. Substantial and continuously
  present (n=922, 3–7% share every quarter since 2020) but flat, not rising
  — noise around a stable baseline. Doesn't breach plain EB05; only
  `eb05-spatial` alarms on it (see `benchmark.md`).

`no_emergence=True` scores each against zero ground-truth-positive quarters,
so any alarm is scored as pure false alarm rather than silently excluded.

**Codeine moved the other direction — from "considered and excluded" to a
labelled positive — on a second, closer look.** The first pass used 3-year
block averages and read it as flat noise (peak share ~1.3%). A quarter-level
look shows something different: low and noisy through 2019–2023 (mostly
0.1–0.6%), then a real step-up sustained from 2023Q4 (0.90%) through 2025Q4
(consistently 0.6–1.3%) — a modest but genuine elevation relative to its own
prior baseline, same "read the rate off share" method as Fentanyl/PCP/
Ketamine. **The lesson generalizes: block-level summaries can hide a real
signal a quarter-level read catches, in both directions** — that's exactly
what almost mislabelled Codeine as noise and what caught it before the
mistake shipped.

**Two candidates remain deliberately unlabelled** after the same scrutiny:
Amphetamine (later-quarter readings hint at a possible real uptick, same way
Codeine's did, but the pattern is noisier and less resolved — a genuine "not
enough evidence yet" rather than a confident call either way) and
1,1-Difluoroethane (scattered singles throughout, no sustained pattern — the
same failure mode `nflis_onset` documents for mitragynine's spike).
Buprenorphine, Methocarbamol, Acetaminophen, Hydroxyzine, and Lamotrigine
were also checked and are all too thin (lifetime n < 30, no sustained
pattern in either direction) to responsibly call anything.

## Result, all 22 labelled substances

| Substance | IoU | Recall | Precision | Ground truth | Detected |
|---|---|---|---|---|---|
| Acetyl fentanyl | 0.80 | 0.80 | 1.00 | 2022Q1–2024Q2 | 2022Q2–2024Q1 |
| Flualprazolam | 0.80 | 0.80 | 1.00 | 2020Q1–2021Q1 | 2020Q2–2021Q1 |
| Fentanyl | 0.69 | 0.69 | 1.00 | 2016Q2–2023Q2 | 2016Q3–2021Q2 |
| Bromazolam | 0.71 | 0.71 | 1.00 | 2023Q3–2025Q1 | 2024Q1–2025Q1 |
| Ephedrine | 0.50 | 0.50 | 1.00 | 2024Q1–2024Q4 | 2024Q3–2024Q4 |
| Carfentanil | 0.50 | 0.50 | 1.00 | 2024Q1–2025Q4 | 2024Q3–2025Q2 |
| Etizolam | 0.44 | 0.44 | 1.00 | 2019Q4–2021Q4 | 2020Q1–2020Q4 |
| para-Fluorofentanyl | 0.40 | 0.40 | 1.00 | 2021Q1–2025Q4 | 2021Q2–2023Q1 |
| Benzylfentanyl | 0.25 | 0.50 | **0.33** | 2024Q2–2024Q3 | 2024Q3–2025Q1 |
| Clonazepam | 0.23 | 0.23 | 1.00 | 2019Q3–2022Q3 | 2020Q1–2020Q3 |
| Mitragynine | 0.13 | 0.13 | 1.00 | 2022Q2–2025Q4 | 2022Q4–2023Q1 |
| Xylazine | 0.14 | 0.14 | 1.00 | 2023Q3–2025Q1 | 2024Q3–2024Q3 |
| PCP | 0.17 | 0.17 | 1.00 | 2020Q2–2025Q4 | 2020Q3–2021Q2 |
| Ketamine | 0.14 | 0.14 | 1.00 | 2022Q3–2025Q4 | 2022Q4–2023Q1 |
| Codeine | 0.22 | 0.22 | 1.00 | 2023Q4–2025Q4 | 2024Q3–2024Q4 |
| **Morphine** | 0.00 | n/a | **0.00** | *(none — negative case)* | 2025Q2–2025Q2 |
| **Lidocaine** | 0.00 | n/a | **0.00** | *(none — negative case)* | 2024Q4–2025Q3 |
| **Levamisole** | 0.00 | n/a | **0.00** | *(none — negative case)* | 2022Q1–2022Q1 |
| **Mirtazapine** | 0.00 | n/a | **0.00** | *(none — negative case)* | 2024Q3–2024Q3 |
| **Citalopram** | 0.00 | n/a | **0.00** | *(none — negative case)* | 2017Q4; 2024Q4 |
| **Alprazolam** | n/a | n/a | n/a | *(none — negative case)* | *(never alarms)* |
| **Hydromorphone** | n/a | n/a | n/a | *(none — negative case)* | *(never alarms)* |

**mean IoU 0.307, mean recall 0.426, mean precision 0.717** across the 20
substances plain EB05 actually alarms on (Alprazolam and Hydromorphone
contribute nothing to plain EB05's numbers — both are `eb05-spatial`-only
false positives, never triggering plain EB05 at all; see `benchmark.md`).

**Substance-level false-positive rate, with a Wilson 95% confidence
interval** (`emerging.validation.ground_truth.precision_rate` — one
independent trial per alarmed substance: did the alarm ever touch real
ground truth, yes/no; see the function's docstring for why this is
computed per-substance rather than pooled over alarm-quarters, which would
treat correlated quarters within one episode as independent and understate
the true uncertainty): **15 of 20 alarmed substances had a real hit — 0.75,
95% CI [0.53, 0.89].** That interval is wide — with n=20 it has to be — but
it's now a real interval instead of a bare point estimate implying more
precision than the data supports. Citalopram alarms twice, in two separate
episodes (2017Q4 and 2024Q4) that don't overlap either — both count as
false alarms.

Lidocaine alarms for four straight quarters (2024Q4–2025Q3) — the single
worst false-positive episode in the file, and the one this project's
`--weighted` case definition was specifically built to fix (see below).
Benzylfentanyl's 0.33 is a different kind of imperfection — a real,
correctly-detected positive whose interval is just imprecise, not a false
alarm; don't conflate the two when reading this table.

**Does `--weighted` actually fix Lidocaine? Partially.** Checked directly
against `results/benchmark/per_substance.csv`: plain `eb05` alarms on
Lidocaine for 4 quarters; `eb05-weighted` and `eb05-v2` cut that to 1 —
a real, large reduction, and a direct validation of the position-weighted
case definition's motivating case (`GPS_V2_DESIGN.md` #1 was built
*because of* Lidocaine) — but not a complete fix. The same two models also
correctly avoid Morphine and Levamisole entirely (0 alarmed quarters, not
just fewer). `eb05-dual` does not help with Lidocaine at all (still 4
quarters) but does avoid Morphine. No single variant is clean across all
seven negatives — see `benchmark.md` for the full breakdown, now including
Wilson CIs per model.

`emerging/validation/benchmark.py` (built separately, scoring several
detector variants against this same reference) needed the identical
`no_emergence` fix — it read `known_emergences.csv` directly and would have
silently scored Morphine's placeholder interval as a real one otherwise.

## The recall numbers need a caveat before they're read as detector quality

Recall is low for three of five, and the reason is not that EB05 missed an
emergence — it's a **construct mismatch between what the reference measures
and what EB05 measures.** The reference interval is *presence*: is the
substance still showing up in LA deaths at all, above a near-absence floor.
EB05 measures *growth*: is its share of overdose deaths still rising relative
to its own recent baseline. A substance can stay present for years after its
growth phase ends — para-Fluorofentanyl is the clean case: raw counts never
return to near-zero through 2025Q4 (reference stays "ongoing" for 20
quarters), but EB05 correctly stops alarming once the rate of increase
levels off around 2023Q1, which is also exactly what `RESEARCH_LOG.md`
already reports independently as a credible *decline* (EB95 0.37). The
detector is doing its job; the reference interval is answering a related but
different question for the back half of that window. Mitragynine and
xylazine show the same pattern for the same reason — both have long
presence-based reference windows built on `ongoing: True` while their actual
EB05-alarm episodes are much shorter.

**What this means for reading the metric.** Precision here is the trustworthy
number: it says the detector didn't cry wolf, and it isn't vulnerable to the
presence-vs-growth conflation since a detected quarter falling inside a
presence window is unambiguous regardless of which construct produced the
window. Recall against a presence-based reference should be read as "coverage
of the substance's active local history," not "coverage of its emergence" —
conflating the two would make EB05 look worse than it is on exactly the
substances (para-Fluorofentanyl, mitragynine) where it is doing what it's
designed to do: stop alarming once growth stops, not once the substance
disappears.

**Open question, not resolved here:** should the reference schema
distinguish an "actively rising" sub-interval from the outer "present"
interval, so recall can be computed against the former without needing this
caveat every time? That would need its own annotation pass and its own stated
rule — deferred rather than guessed at.

## Does LA lag the nation, and is "the nation" even the right comparator?

`national_context` in the CSV was added as context, explicitly not as ground
truth for LA's own timing. Checked directly, substance by substance:

| Substance | National/regional marker | Apparent LA lag | Driven by |
|---|---|---|---|
| para-Fluorofentanyl | First national deaths 2020Q3; peak 2021Q2 (CDC MMWR 71(39), SUDORS) | **~0** — LA's onset (2021Q1) precedes the national peak | Rides LA's *already-established* fentanyl network |
| Bromazolam | National designer-benzo supply share 4%→73%, 2021→mid-2023; surge "late 2022" | **~2–3 quarters** | Genuinely new to LA's market |
| Carfentanil | Near-zero nationally through 2022; 3400% detection rise 2023; deaths triple 2024 | **~1 year** behind the initial 2023 signal, but synced with the 2024 escalation | Responded to scale, not the earliest signal |
| Xylazine | DEA national alert + CA DPH alert, both ~2023Q1 | **~2 quarters**, and *real* — LA fentanyl-supply positivity was 0% the same quarter the alerts fired (independent drug-checking data, not this project's mortality data) | Different regional market — xylazine is Northeast-concentrated; the "national" signal is substantially a Philadelphia/Boston/NYC signal |
| Mitragynine | No comparable marker exists | n/a | Not an outbreak — legal supplement, gradual use growth |

**A single "LA lags the nation by N quarters" figure would be actively
misleading — the spread across these five (0 to ~4 quarters, plus one
substance the framework doesn't even apply to) is the finding, not noise
around a stable number.** Two different mechanisms produce it, pointing
opposite directions: a new fentanyl analog riding an already-established
local trafficking network shows no lag at all, while a substance requiring a
new adulteration practice or supply-chain change to take hold shows a real,
multi-quarter one.

**The xylazine row is the clearest case that "national" is often the wrong
comparator entirely, not just an imprecise one.** Xylazine's documented
concentration in specific Northeast cities means its national signal is not
a spatially uniform epidemic curve LA is a lagging point on — it's
substantially a different region's local epidemic, reported at national
scope. An independent LA-specific drug-checking study (DART-MS on fentanyl
samples from a community program, medRxiv 10.1101/2025.05.13.25327478 —
unrelated to LACME mortality data) puts LA fentanyl-supply xylazine
positivity at 0% in 2023Q1, the same quarter the DEA and CA DPH alerts
fired, rising to 29.5% by 2025Q1. That confirms two things at once: the
2-quarter gap between the alerts and LA's first xylazine deaths is a real
local-market lag, not a reporting-delay artifact of the ME data — and the
national alert genuinely predated any LA supply presence, so it was never
going to be a good predictor of *when*, only a warning that it eventually
would. Xylazine's row confidence was raised from `low` to `medium` on the
strength of this independent corroboration.

## `national_context` wasn't good enough, and the fix is `health_notices.csv`

The table above was built entirely from WebSearch tool summaries — never a
directly-read source. Asked to do better, the fix wasn't a bigger dataset; it
was **checking the actually-comparable jurisdiction first.** LA County's own
health department, not "the nation," is the peer to LA County ME data, so
`data/raw/ground_truth/health_notices.csv` now records one row per
(substance, jurisdiction), county through federal, each with an
`evidence_basis` and a `source_verification` field, populated from what was
actually fetched and read rather than summarized.

**Reading LA County DPH's actual xylazine alert (2023-03-08) directly changed
the finding.** It cites **zero LA County detection data.** It warns xylazine
is "likely present" in LA's supply *because San Francisco and San Diego had
found it*, not because LA had. The county's own alert was anticipatory, not
evidence — a distinction the previous, search-summary-only version of this
table could not have caught, because it never read past the headline claim
"CA DPH issued an alert." This is now the strongest single citation in the
file, specifically because it was fetched directly instead of summarized.

**No LA County or California-level notice exists for bromazolam, carfentanil,
or para-fluorofentanyl**, despite real searching — and that absence is a
finding, not a gap to eventually fill. All three are genuine local
emergences in this project's own ME data; none generated a local government
alert the way xylazine's "zombie drug" coverage did. **A notices-based
reference is structurally biased toward whichever substances become
media-salient, not toward whichever are actually most dangerous or arrive
first** — a real limitation of the whole "compose it from DOH notices"
approach, stated in the file rather than papered over by quietly falling
back to a federal citation that reads the same as a verified local one.

A genuine count dataset (not annual-report PDFs, which were set aside as too
much extraction work for what they'd buy) turned out to be available a
different way: the user has direct access to DEA's **NFLIS-Drug Data Query
System**, and pulled a real export — California, semi-annual drug
identification counts, 2012–2025, for the designer-benzodiazepine and
fentanyl-analog categories. This is the strongest external comparator in the
file by a wide margin: a real forensic-lab count series, state-level (a
better-matched jurisdiction than "national" to begin with), not a notice
date and not a search snippet.

## The NFLIS cross-check revises two of the earlier lag estimates, substantially

`data/raw/ground_truth/nflis_ca_counts.csv` (cleaned from the raw DQS export
in `data/raw/nflis/`; `BASE_DESCRIPTION` already merges salt-form variants,
and `Fluorofentanyl`/`para-Fluorofentanyl` were merged by hand to match this
project's own `ISOMER_ROLLUP`). `emerging.validation.ground_truth.nflis_onset`
applies a rule in the same spirit as `trends.py`'s `raw-count-rule` — first
half-year >= 10, after two consecutive half-years < 10 — run via
`emerging ground-truth nflis-lag`.

| Substance | CA NFLIS onset | LA ME-death onset | Lag |
|---|---|---|---|
| Xylazine | 2019H1 | 2023Q3 | **~4.5 years** |
| Bromazolam | 2021H2 | 2023Q3 | **~2 years** |
| Carfentanil | 2024H1 | 2024Q1 | ~0 |
| para-Fluorofentanyl | 2020H2 | 2021Q1 | ~0.5 years |
| Mitragynine | none (2015H1 reads as one, but collapses right back to n=2 the next half-year — noise, not a trend) | 2022Q2 | not comparable |

**Xylazine's lag is roughly triple the earlier estimate, and it's not the
same claim as before.** The earlier "~2 quarters" figure came from comparing
LA's ME-death onset to the county alert (2023-03-08) and the LA-specific
drug-checking study's 0%-in-2023Q1 reading — both of which measure *LA's own
supply*. CA NFLIS shows xylazine substantial statewide (n=24–61 per
half-year) as early as **2019–2020**, years before it was detectable in LA's
own fentanyl supply at all. Read together rather than as a contradiction,
this is a two-stage story: xylazine took roughly four years to migrate from
"present somewhere in California" to "present in LA's own drug supply," and
then only a couple of quarters more to show up in LA's fatal-overdose data
once it was locally present. The earlier "~2 quarters" number was correct
for the question it was answering; it just wasn't the question that matters
for "does LA lag California."

**Bromazolam's lag roughly doubled too** — 2021H2 statewide vs. the
secondhand "late 2022" national-surge citation the earlier version relied
on. The real CA data shows bromazolam taking off a full year earlier than
the best available secondhand source suggested.

**Carfentanil and para-fluorofentanyl are confirmed, not revised** — both
were already read as near-zero lag from weaker sourcing, and the real CA
series lines up almost exactly with LA's own onset for both. Mitragynine's
"no coherent trend" read is also confirmed, now against real forensic data
rather than just the absence of a comparable marker.

**A methodological note worth keeping, because it nearly produced a wrong
answer:** the raw DQS export omits zero-count half-years as rows entirely
rather than writing them as 0. A first pass at `nflis_onset` that scanned
consecutive *rows* of the sparse series missed para-fluorofentanyl's onset
completely, because its first nonzero row had no rows before it to check —
two periods that were actually years apart elsewhere in the data would have
silently read as adjacent. Fixed by reindexing onto a complete half-year
calendar before scanning; both the bug and the fix are pinned down in
`tests/test_ground_truth.py`.
