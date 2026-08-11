# Research log

Reverse-chronological. Each entry records what changed, what it showed, and
what it cost. Findings that later turned out to be artifacts are kept, marked,
not deleted — the artifact is usually the more useful record.

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

**One substantive spatial finding:** PCP, 1.38× tighter than the overdose
background, p=0.005, n=367, anchored on 90013 (Skid Row). The flagged *novel*
substances are not localized — they track the county-wide fentanyl supply.

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
five known emergences. The risk noted here — "standalone C++ binary, breaks the
pure-Python pipeline" — was avoided by implementing the method directly; it is
~400 lines and runs the full 45-quarter sweep in 2.5 minutes.

### P2 — SaTScan spatial variation in temporal trends `not started`

Answers *where is it growing*, which is the question the project actually asks;
`geo.py` only answers *where is it concentrated*. Fuses the temporal and
spatial tests into one statistic. The space-time permutation variant needs no
population denominator, which matches our constraint exactly.

### P3 — Delay-distribution nowcast `not started`

Replace the hard `CUTOFF` (which deletes 2026 outright) with a fitted
reporting-delay distribution, so recent quarters are usable with honest
uncertainty. Note the CDC precedent *and* its documented failure when the
epidemic changes regime — which ours is doing.

### P4 — `sparr` relative-risk surface for PCP only `not started`

n=367 clears the n ≥ 50–100 threshold that ruled kernel methods out. Turns
"1.38× tighter" into an actual map. R detour; PCP is the only substance that
qualifies.

### P5 — READUS-PV reporting conventions in the manuscript `not started`

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
