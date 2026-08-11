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
| Interpretation | `polysubstance.py` | cause-of-death vs passenger verdicts |
| Localization | `geo.py` | case-control random-labelling permutation test |

**Headline results, settled window (recent 2025Q1–Q4 vs baseline 2023Q1–2024Q4;
2,153 vs 5,491 OD deaths):**

- GPS prior Gamma(α=7.267, β=7.547), fitted on 113 substances.
- **Nothing currently clears the EB05 > 1.5 alarm line.** Top: PCP 1.27,
  Alprazolam 1.11, Cocaine 1.06, Lidocaine 1.04.
- 24 of 205 substances have breached 1.5 at some point across 45 as-of quarters.
- Credible *declines*: Fentanyl (EB95 0.88), para-Fluorofentanyl (EB95 0.37).

**The four findings that changed what we can claim:**

1. **Lidocaine is an adulterant, not a killer.** 17 of 17 recent deaths name
   fentanyl, 0 name it alone, and it sits at 0.93 on the ME's cause line
   (fentanyl: 0.13). The correct headline is "the fentanyl supply is being cut
   differently".
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

### P1 — TreeScan, prospective tree-temporal `not started`

The single highest-value change. Replaces one-substance-at-a-time testing with
a scan over the NFLIS category hierarchy, which is already in the data (23
categories, fully populated).

Fixes three things at once:
- **Aggregation power** — nitazenes are currently "individually too rare to
  score" and needed a hand-coded family workaround. TreeScan tests the branch.
- **Multiple testing** — 205 substances × 45 quarters is presently unadjusted.
- **Repeated looks** — replaces the hand-picked EB05 > 1.5 line with a
  recurrence interval.

Steps:
1. Export a `root → NFLIS category → substance` tree plus quarterly counts in
   TreeScan's input format.
2. Run the tree-temporal model with the same as-of discipline as
   `trends backtest`.
3. **Head-to-head against EB05 on the five known emergences** — same
   ground truth, so the comparison is honest. Report detection delay for each.
4. Decide whether TreeScan replaces EB05 or runs alongside it.

Risk: standalone C++ binary, not a Python library — breaks the pure-Python
pipeline. Nomenclature lag means a genuinely novel analog may sit uncategorised
at the root.

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

- **MaxSPRT** — redundant with TreeScan's recurrence interval for our use.
- **LGCP / inlabru** — prior-dominated at n=9–39; revisit only if a single
  substance becomes a funded investigation.
- **Two-component MGPS** — not identifiable at this scale; documented in
  `_fit_gps_prior`. Revisit on a substance × zip × quarter table.
- **Spatial methods as *detectors*** — they stay downstream of the temporal
  detector. At n=9 there is nothing to scan.

---

## 2026-08-10 — Reviewer feedback round

Six items from a presentation. All six built.

**Polysubstance / principal cause** (`polysubstance.py`, new). Two independent
discriminators. Co-occurrence (`alone_pct`, `with_fentanyl_pct`) and — the
better one — **position on the ME's combined cause line**, which is ordered by
clinical significance, not alphabetically. Fentanyl sits at mean relative
position 0.13, lidocaine at 0.93.

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
