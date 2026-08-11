# Emerging-substance detection (LA County)

Recovers named substances from LACME cause-of-death narratives and ranks which
ones are rising. Scoped to LA County only.

## Why this exists

The modeling panel carries 8 coded substance flags and collapses everything else
into `Others` (4.8% of post-2022 deaths, with no indication of what the "other"
was). The narrative text is far richer — `EFFECTS OF FENTANYL,
METHAMPHETAMINE, AND LIDOCAINE` — and names xylazine, bromazolam, carfentanil,
mitragynine, para-fluorofentanyl and ~200 more. Extraction turns 8 columns into
**219 substances at 99.8% decedent coverage**, which is the input every method
in this directory depends on.

## Pipeline

```bash
python -m emerging.lexicon      build     # alias -> canonical map
python -m emerging.extract      extract   # decedent x substance table
python -m emerging.extract      validate  # score vs the 8 coded flags
python -m emerging.trends       rank      # EB05 ranking
python -m emerging.trends       backtest  # does the detector fire?
python -m emerging.trends       alarms    # every breach, ever  (run before plot)
python -m emerging.trends       watch     # national watchlist vs LA
python -m emerging.trends       regime    # recording artifacts that move EB05
python -m emerging.trends       plot      # five figures
python -m emerging.polysubstance profile  # cause vs passenger
python -m emerging.polysubstance plot
python -m emerging.geo          cluster   # is it localized?
python -m emerging.geo          plot
```

`alarms` must run before `plot`: it writes `alarm_history.csv` and
`eb05_sweep.csv`, which supply the breach annotations and the alarm-timeline
figure. `plot` warns and skips those rather than recomputing a 45-quarter
sweep.

| Module | Role |
|---|---|
| `lexicon.py` | 4,089 aliases → 3,127 canonicals from the DEA NFLIS catalog + a curated supplement (ethanol, ME misspellings, metabolite/isomer rollups) |
| `extract.py` | greedy longest-n-gram matcher over the narrative fields, with an audited fuzzy fallback for ME typos |
| `trends.py` | gamma-Poisson empirical-Bayes ranking, as-of backtest, alarm history, watchlist, recording diagnostics, figures |
| `polysubstance.py` | is a flagged substance a cause of death or a passenger — co-occurrence plus position on the ME's cause line |
| `geo.py` | case-control permutation test for spatial localization against the overdose-death background |

| `cohort.py` | overdose-death cohort definition — **a vendored copy**, see below |
| `paths.py` | where the data lives |

## Setup

```bash
pip install -e .
```

All commands are run from the repository root. Paths are fixed (`emerging/paths.py`):

| Path | In git? | Notes |
|---|---|---|
| `data/raw/lacme/` | **no** | LACME extract — **contains decedent names** |
| `data/processed/` | **no** | one row per decedent per substance; regenerate in ~1 min |
| `data/raw/nflis/` | yes | DEA NFLIS catalog, public |
| `data/raw/_geo/` | yes | LA County zipcode boundaries, public |
| `results/` | yes | aggregate only — no case numbers, no coordinates |

**Do not add `.gitignore` exceptions for the first two.** Regenerate the
processed tables with `lexicon build` then `extract extract`.

### The cohort definition is duplicated, and that is a liability

`cohort.py` (`filter_alcohol_only`, `_norm_zip`) is a **vendored copy**. The
authoritative version lives in the `predict_rosla` repository at
`rosla_nowcast/modeling/zipcode/monthly/signal_check/ccf.py` and is itself kept
in sync with the epi project. It is duplicated here only because this
repository is standalone.

The point of the shared definition is that counts here stay commensurable with
`deaths` in `panel_zip_quarter.parquet`. If upstream changes and this copy does
not, every count silently stops matching and **nothing fails loudly**. Diff
before publishing any number that gets compared to the modeling work.
Copied 2026-08-10 from `predict_rosla` @ e357baf.

## Documents

- [`RESEARCH_LOG.md`](RESEARCH_LOG.md) — current state, findings, game plan, and
  the artifacts that would have manufactured false results
- [`LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md) — rare-event detection and
  spatial clustering methods; what is better than this and what is not

## Validation

Extraction scored against the ME's own coded flags — the only ground truth
available (`results/flag_validation.csv`):

| Flag | Precision | Recall | F1 |
|---|---|---|---|
| Fentanyl (incl. analogs) | 0.997 | 0.998 | 0.997 |
| Heroin | 0.995 | 0.998 | 0.996 |
| Methamphetamine | 0.993 | 0.989 | 0.991 |
| Cocaine | 0.985 | 0.997 | 0.991 |
| Alcohol | 0.986 | 0.962 | 0.974 |
| Benzodiazepines | 0.990 | 0.936 | 0.962 |

Disagreement is not automatically extractor error — the coded flags are
themselves regex-derived — so both directions are reported.

### Why coverage is 99.8% and not 100%

55 of 22,004 deaths (0.25%) get no named substance. Almost none of it is a
fixable extractor defect:

| Cause | n | Fixable? |
|---|---|---|
| The examiner named a **class, not a drug** — "opioid intoxication", "effects of opioids" | 34 | No. The record does not say which opioid. |
| **Chronic alcohol** phrasing — alcoholic hepatitis / pancreatitis / cardiomyopathy, chronic alcoholism | 14 | No, by design. These are chronic conditions contributing to a death, not acute substance involvement. |
| **Mid-word OCR damage** — `ALCOHOLIC INTOX ICATION`, `100COLCHI CINE` | 4 | Not without intra-word repair, which risks more than it recovers. |
| **Not overdoses at all** — anaphylaxis to ceftriaxone / Zosyn, pulmonary coccidioidomycosis | 3 | No — these are cohort-definition edge cases. |

Getting from 99.4% to 99.8% came from three fixes, all worth knowing about:

1. **The hyphen rule was too strict.** It only split a hyphenated token when
   *every* part resolved, so `METHAMPHETAMINE-INDUCED`, `COCAINE-ASSOCIATED`,
   `ETHANOL-OLEANDRIN` and `ALPRAZOLAM--HELIUM` all discarded a known substance
   because one sibling token was a modifier or an uncatalogued compound. It now
   keeps whichever parts resolve. This is safe because whole aliases are matched
   and consumed first, so `PARA-FLUOROFENTANYL` never reaches the split branch.
2. **`ALCOHOLIC INTOXICATION`** — the adjective form was the single largest gap
   (67 deaths). Only the two-word acute phrasings are mapped; a bare `ALCOHOLIC`
   would also fire on alcoholic cirrhosis and hepatitis.
3. **Brand names and non-catalog poisons** — Xanax, OxyContin, arsenic, zinc,
   colchicine, oleandrin. NFLIS catalogues trafficked *drugs*, so poisons and
   plant toxins that appear in overdose casework are simply absent from it.

**Do not fix the rest by loosening the fuzzy cutoff.** Tested at 0.85, it
recovers a few real typos and admits six false positives, the worst being bare
`ALCOHOLIC` → Ethanol at 104 mentions — which would have tagged alcoholic
cirrhosis as acute alcohol involvement and corrupted the ethanol series. Those
tokens are now in `FUZZY_BLOCK` so the mistake cannot be reintroduced by
changing one number. Verified typos go in `SUPPLEMENT` by name instead.

One of those by-name additions changed a result: `CARENTANYL EFFECTS` (2017-06)
is a real carfentanil death, moving carfentanil's first LA appearance from
2024Q1 to 2017Q2 — an isolated case during the national 2016–17 carfentanil
wave, then a seven-year gap before the 2024 cluster.

Two review surfaces, both worth reading before trusting a result:
`extract unmatched` ranks tokens the lexicon never consumed (where a genuinely
novel substance would first surface), and `fuzzy_resolutions.csv` logs every
typo substitution the matcher guessed at.

## Does the ranking work?

`backtest` re-runs the ranking as of each historical quarter, scoring only data
available then. Without it a detector that never fires and a detector that is
broken look identical.

| Substance | Peak EB05 | Best rank | First top-10 |
|---|---|---|---|
| para-Fluorofentanyl | 7.00 | 1 | 2021Q1 |
| Carfentanil | 5.02 | 1 | 2024Q2 |
| Bromazolam | 4.78 | 1 | 2023Q4 |
| Mitragynine | 2.36 | 2 | 2022Q4 |
| Xylazine | 1.93 | 4 | 2024Q2 |

All five reach the top 4 in the quarter they emerge — **including xylazine at
n=6**, which is a stronger result than the raw counts suggest. The shrinkage is
what makes this possible: it measures a share shift against ~600 deaths per
quarter, not the count in isolation.

This calibrates a threshold: real emergences score EB05 1.93–7.00, so
**EB05 > 1.5** is a defensible alarm line.

### Nothing currently clears the line — but watch lidocaine

On the settled window (recent 2025Q1–2025Q4), the highest EB05 across all 219
substances is **1.27** (PCP). Nothing reaches 1.5.

**Lidocaine is the one to keep an eye on, and it is a lesson in windowing.**
Its counts run 2 → 4 → **10 → 11** → 2 → 2 → 2 across 2024Q2–2025Q4: a sharp
burst straddling the 2024/2025 boundary, then a return to baseline. Where the
4-quarter window falls decides what you see:

| Recent window | Lidocaine n | EB05 | Verdict |
|---|---|---|---|
| 2025Q1–2025Q4 (settled, current default) | 17 | 1.04 | below the line |
| 2024Q4–2025Q3 (window shifted one quarter) | 25 | 2.20 | above the line |
| Peak over the as-of sweep (2025Q1) | 27 | **4.29** | rank 1 |

All three are correct; they answer slightly different questions. The substance
is real — 43 lifetime cases, all but one named in the primary `CauseA` field
rather than stray narrative text, and 38 of 43 co-occurring with fentanyl, the
adulterant-attaching-to-the-hub pattern. What the table shows is that a
**short burst is window-sensitive in a way a sustained trend is not.**

Practical rule: for any spiky trajectory, quote the `backtest` peak alongside
the current value. A single-window EB05 will under- or over-state a burst
depending on where the boundary lands.

## How the expected count is computed

`E_s` = the substance's **share of all overdose deaths in the baseline window**,
applied to the **recent window's death total**:

```
E_s = n_baseline_s x (N_recent / N_baseline)
```

Two fixed adjacent windows, not a fitted trend: recent = the last 4 quarters,
baseline = the 8 immediately before it. The baseline rolls forward with the
data, which is what makes this an *emergence* detector — a substance that has
been high for years does not flag, because its own recent past is the yardstick.

Worked example (lidocaine, current window):

| | |
|---|---|
| baseline 2023Q1–2024Q4 | 22 lidocaine deaths / 5,491 OD deaths = **0.401%** |
| recent 2025Q1–Q4 | 2,153 OD deaths |
| expected | 0.401% x 2,153 = **8.63** |
| observed | **17** → raw ratio 1.97 → EBGM 1.47, EB05 1.04 |

### The denominator is deaths, not people — and that matters

This is a **compositional** measure: share of overdose deaths, not a population
rate. LA's total overdose deaths fell **22%** between the two windows (686/qtr →
538/qtr), so a substance with perfectly flat absolute counts gains ~1.28x share
automatically.

The consequence is that "rising" in this table can mean "a larger slice of a
shrinking pie":

| Substance | deaths/qtr, baseline → recent | absolute change | share ratio | EB05 |
|---|---|---|---|---|
| Lidocaine | 2.8 → 4.3 | **1.55x** | 1.97 | 1.04 |
| PCP | 28.9 → 34.3 | **1.19x** | 1.51 | 1.27 |
| Cocaine | 126.1 → 113.8 | **0.90x** | 1.15 | 1.06 |
| Methamphetamine | 416.0 → 336.3 | **0.81x** | 1.03 | 0.98 |

Cocaine clears the "credible growth" bar on share while its actual death count
*fell*. That is not a bug — share is the right frame for "is the drug supply
changing" — but **"cocaine is rising" is the wrong sentence.** Say "cocaine is a
growing share of a shrinking total," and check the absolute column before
describing any substance as increasing.

## Reading the output — four caveats that change conclusions

**1. Right-truncation is the biggest practical hazard.** Toxicology takes
months, so the last two quarters are undercounts and every rising substance
looks like it is falling. `--lag-quarters` (default 2) drops them before any
statistic is computed, and the figures draw the tail dashed rather than hiding
it. This is a *floor, not a fix*: novel substances need send-out confirmation
and so confirm more slowly than fentanyl, meaning delay correlates with exactly
the thing being detected. Properly handling this needs a reporting triangle —
LACME would have to supply an adjudication or report date alongside death date.
**Worth requesting; it is cheap to ask for and currently unavailable.**

Completeness is handled in `load_quarterly`, before any statistic runs, by two
separate guards:

- **`--cutoff` (default `2025-12-31`)** drops quarters whose toxicology is not
  yet settled. This is a *data judgement, not a derived quantity* — **advance it
  by hand when the LACME extract is refreshed** and the newer quarters have had
  time to adjudicate. Leaving it stale silently shrinks the analysis window.
- **Partial calendar quarters are always dropped.** The extract ends 2026-04-04,
  so 2026Q2 is a 4-day stub holding 3 deaths. Left in, it looks like an ordinary
  quarter to the older `--lag-quarters` haircut and silently absorbs one of its
  slots, so a "2-quarter" correction discards only one real quarter.

`--lag-quarters` is the older rolling alternative to the cutoff and now defaults
to **0**: the two express the same correction, and stacking them throws away a
year of settled data. Use one or the other.

**2. `date_added` is catalog novelty, not street novelty.** It records when
NFLIS could report a substance. 852 of 3,101 entries (27%) carry the catalog's
1998-10 inception stamp, including para-fluorofentanyl and carfentanil — both of
which only reached LA deaths in 2021 and 2024. `date_added_censored` marks these.
Local arrival (`first_seen`) is the more useful axis, and is itself censored at
2012 where the extract starts.

The corollary is the shaded region in `novelty_vs_growth.png`. Three substances
sit below the diagonal, meaning the LA death predates the NFLIS entry:
**chlorcyclizine** (LA 2012Q3, NFLIS 2023-10), **desloratadine** (LA 2014Q2,
NFLIS 2023-10), and **U-47700** (a one-month gap, effectively on the line).
This is not LA detecting something before the DEA. NFLIS catalogues substances
as forensic labs report them **in seized drug evidence**; the first two are
prescription antihistamines that turn up in postmortem blood but are not
trafficked, so they entered the catalog late. The region is a statement about
NFLIS's sampling frame, not about local surveillance being ahead.

**3. Naming variants can manufacture a fake emergence.** LACME writes both
`PARA-FLUOROFENTANYL` and a bare `FLUOROFENTANYL`; the second form takes over in
2023Q2 when the office changed convention. Unrolled, that reads as one substance
collapsing and another appearing. `ISOMER_ROLLUP` collapses them. Assume more
cases like this exist and check `unmatched_tokens.csv` when a new substance
appears to jump.

**4. The prior is single-component, not MGPS.** DuMouchel's two-gamma mixture
was tried first and is not identifiable here — 203 substances, only 111 with a
non-zero baseline, cannot pin down five parameters, and the optimizer converged
to different degenerate solutions depending on the fit set. Details in
`_fit_gps_prior`. The mixture is worth revisiting on a substance × zip × quarter
table.

## What this does not do

- **Detect a substance nobody assayed.** If the ME's panel does not test for it
  and the narrative does not name it, no count model recovers it. That bound is
  chemistry and lab scope, not statistics.
- **Spatial *detection*.** Case-control relative-risk surfaces (`sparr`) and
  Kulldorff scans need roughly n≥50–100 cases; only para-fluorofentanyl and
  mitragynine clear that. Spatial remains a characterization tool for
  substances the temporal detector already named — which is exactly what
  `geo.py` does, and it is deliberately downstream of `trends`, never a
  detector in its own right.
- **Occupancy-detection modeling.** Separating "absent" from "not tested for"
  needs detection probability to vary across units. One ME office running one
  evolving panel gives no such variation — every zip shares the same lab. This
  needs multi-jurisdiction data to be identifiable.
- **Substance co-occurrence *graph*.** Post-2022 prevalence is methamphetamine
  54% / fentanyl 54% / everything else under 6% — a two-hub star graph where
  centrality trajectories carry little discriminating information, so it is
  not a detector. Pairwise co-occurrence is still highly informative as an
  *interpretive* layer, and that is what `polysubstance.py` uses it for.
- **Identify the policy behind a regime change.** `regime` flags the quarters
  where the denominator, the panel, or the onset pattern shifted. Matching
  those to naloxone saturation, a scheduling change, or an ME protocol
  revision needs sources outside this dataset.
- **Distinguish "absent" from "not assayed".** `watch` reports zeros for
  medetomidine, dexmedetomidine, BTMPS and desomorphine. Every one of those
  names resolves in the lexicon, so the scan looked and found nothing in the
  narratives — but an assay the ME does not run produces an identical zero.
  Only the ME's panel can separate them.

## Is it a cause of death, or a passenger?

`polysubstance.py`. EB05 measures growth in share; it cannot say whether
anyone died *of* the substance. Two independent discriminators, over the last
4 settled quarters:

- **Company kept.** `alone_pct` (deaths naming nothing else) and
  `with_fentanyl_pct`.
- **Position on the cause line.** The ME writes one combined string —
  `EFFECTS OF FENTANYL, METHAMPHETAMINE, AND LIDOCAINE` — ordered by clinical
  significance, not alphabetically. `lead_pct` is how often the substance is
  named first, `last_pct` how often it is named last, and `mean_rank_pct` its
  average position from 0 (always first) to 1 (always last). All three are
  computed only on multi-substance lines, so a solo mention cannot count as
  "named first".

| | n (2025) | alone | with fentanyl | 1st | last | mean rank | verdict |
|---|---|---|---|---|---|---|---|
| Fentanyl | 1101 | 15% | — | 82% | 9% | 0.13 | principal driver |
| Cocaine | 455 | 26% | 49% | 44% | 31% | 0.43 | principal driver |
| Methamphetamine | 1345 | 37% | 49% | 16% | 66% | 0.75 | independent |
| Alprazolam | 110 | 1% | 75% | 21% | 46% | 0.64 | adulterant/marker |
| PCP | 137 | 0% | 51% | 2% | 81% | 0.90 | secondary/incidental |
| **Lidocaine** | **17** | **0%** | **100%** | **0%** | **82%** | **0.93** | **adulterant/marker** |

Lidocaine is a cutting agent. Every one of its 17 recent deaths names
fentanyl, and the ME writes it **last on the cause line in 14 of 17 and first
in none** — the mirror image of fentanyl. The right headline is "the fentanyl
supply is being cut differently", not "lidocaine is killing people".

### Why position is normalized rather than a mean ordinal

A raw mean position (1…N) conflates "named late" with "appears on a long
line": a substance always named last scores 2 on a two-substance line and 5 on
a five-substance line. The bias runs the wrong way here — adulterants
co-occur with more substances and so sit on longer lines, inflating a raw
ordinal exactly for the class we are trying to identify.

Ranking all 26 substances both ways gives **Spearman 0.79**. Methamphetamine
is 20th of 26 raw but 7th normalized (short lines, so position 2 of 2 is
*last*); para-fluorofentanyl is 6th raw but 15th normalized (long lines, so
2.83 is the middle).

The cost: on a two-substance line the score can only be 0 or 1, so short lines
push toward the extremes; single-substance lines are undefined and excluded;
and it is an unweighted mean of ratios. `last_pct` is reported alongside
because the normalized mean is the right thing to rank on but harder to read.

*A discriminator that was tried and does not work:* which field the substance
was found in. Lidocaine is 100% `CauseA`, identical to fentanyl — the ME puts
everything on one line. Position within that line is the signal; presence in
it is not.

## Alarm history — when did each flag actually fire?

`trends alarms` sweeps every substance across all 45 as-of quarters and
records the four dates that matter: first LA death, first EB05 > 1.5 breach,
end of the alarm, and the peak. **24 of 205 substances have ever breached.**

| First breach | Substance | Peak EB05 | Quarters in alarm | Detection delay |
|---|---|---|---|---|
| 2020Q1 | Etizolam | 4.49 | 4 | 28q |
| 2020Q2 | Flualprazolam | 4.84 | 4 | **3q** |
| 2021Q2 | para-Fluorofentanyl | **7.00** | 8 | **1q** |
| 2022Q2 | Acetyl fentanyl | 2.82 | 8 | 27q |
| 2024Q1 | Bromazolam | 4.78 | 5 | 11q |
| 2024Q3 | Carfentanil | 5.02 | 4 | 29q |
| 2024Q3 | Xylazine | 1.93 | 1 | 4q |
| 2024Q4 | Lidocaine | 4.29 | 4 | 12q |

Detection delay is reported **only where arrival is observed.** LACME records
start 2012Q1, so anything already circulating then has a censored
`first_seen`; morphine's apparent "53-quarter delay" is the censoring floor,
not a late catch. Those rows carry `first_seen_censored` and a NaN lag.

## Recording changes that move EB05, without any change in the drug supply

`trends regime`. The statistic is a share, so three things move it for reasons
that have nothing to do with drugs — and all three are measurable:

1. **Denominator drift — this one is live.** Total OD deaths fell 18% between
   the current windows (2,619 → 2,153). Every substance holding a flat count
   gains ×1.22 share for free. This is why `rank` now carries `count_ratio`
   alongside `eb05`: `count_ratio < 1` with `eb05 > 1` means *a growing share
   of a shrinking total*, which is true about supply composition and false
   about deaths. Cocaine is the case in point — EB05 1.06, count ratio 0.90.
2. **Panel drift — historical, not current.** Substances named per death climbs
   1.45 (2012) → 1.89 and plateaus around 2020Q1; the last 12 quarters have
   sd 0.049. Recent flags are not a widening panel, but pre-2020 shares are
   depressed by narrower testing, which is why the heatmap marks that line.
3. **Simultaneous onsets.** **7 substances first breached in 2024Q3** —
   benzylfentanyl, carfentanil, codeine, ephedrine, methocarbamol,
   mirtazapine, xylazine — in a quarter when deaths fell 19% year over year.
   Seven unrelated epidemics starting at once is not the parsimonious reading;
   the denominator is.

## Is it localized?

`geo.py`. A **random-labelling permutation test** on a marked point pattern:
holding the 7,581 death locations fixed, the null asks which n of them are
cases. Overdose deaths already concentrate in Skid Row, the Antelope Valley
and the Harbor, so the comparison is against *n deaths drawn from all overdose
deaths in the same window*, not against the population — the question is
whether a substance is more concentrated than overdose death itself. Statistic
is mean nearest-neighbour distance in metres (UTM 11N), 999 permutations,
2023–2025.

`localization = null mean NND / observed mean NND`. Raw NND is not comparable
across substances — it is dominated by n (a random draw of 9 gives 11.4 km, of
367 gives 1.3 km) — and the ratio cancels that because the null uses the same
n. Under a Poisson idealization, mean NND scales as 1/√density, so 1.39× closer
is roughly 1.9× denser.

**3 of 42 substances beat the background at p<0.05** (about 2 expected by
chance across 42 unadjusted tests, so read the top of the list, not the tail):

| Substance | n | localization | p | top zip |
|---|---|---|---|---|
| PCP | 367 | 1.38× | 0.005 | 90013 (Skid Row) |
| Ketamine | 58 | 1.40× | 0.015 | 91754 |
| Cocaine | 1454 | 1.09× | 0.015 | 90013 |
| Lidocaine | 39 | 1.22× | 0.115 | 90057 |
| Carfentanil | 22 | 0.98× | 0.540 | 91601 |
| Xylazine | 9 | 0.99× | 0.560 | 90061 |

PCP is the one substantive spatial finding: well powered and tightly bound to
Skid Row. The flagged novel substances are *not* localized — they track the
fentanyl supply, which is county-wide.

**Methamphetamine (61% of the cohort) and fentanyl (58%) are withheld, not
scored.** As the case fraction approaches 1 the random-labelling null draw
becomes the case set itself, so the test compares the substance against itself
and returns ≈1 regardless of geography. That is a degenerate question, not a
null result. `MAX_CASE_FRACTION = 0.25` gates it; the values are still in the
CSV under `out_of_scope=True` and should not be quoted. Those two need a
population denominator.

**Power, not just significance.** The smallest localization detectable at
p<0.05 is ~1.97× at n=9, 1.52× at n=22, 1.25× at n=39, 1.09× at n=367. So
xylazine's 0.99× is not evidence of dispersion — it is no evidence at all, and
carfentanil is little better. Below n≈20 a non-significant result is
uninformative.

### Provenance of the method

This is a **random-labelling test** in the marked-point-process sense, using a
Clark–Evans-style nearest-neighbour ratio as the summary statistic. The
closest named ancestor is the **Cuzick–Edwards k-nearest-neighbour test**
(1990), which is the standard case-control clustering test and uses the same
label-permutation null; it counts how many of each case's k nearest neighbours
are also cases, where this uses mean distance to the nearest case. The
canonical alternative is **Diggle–Chetwynd's** D(s) = K_cases(s) −
K_controls(s) (1991), and **Kelsall–Diggle / `sparr`** kernel relative-risk
surfaces are the estimation counterpart.

Why this variant rather than one of those off the shelf: no bandwidth or
scale-radius to choose (Diggle–Chetwynd needs a range of s, `sparr` needs a
kernel bandwidth, and at n=9–39 those choices drive the answer), one number
per substance so 42 of them tabulate, and an exact permutation p-value with no
asymptotics. What it gives up: Cuzick–Edwards has known power properties and
this does not, and D(s) reports the *scale* at which clustering occurs where
this collapses everything to the nearest-neighbour scale. If a substance ever
becomes the object of a real spatial investigation rather than a triage
column, use Diggle–Chetwynd or `sparr` on it.

Three data defects this had to control for, all found while building it:
46 deaths geocoded to **Namibia** (handled by `in_la_county`); **16%** of
deaths at coordinates shared by ≥5 others — hospitals and jails, where people
are pronounced rather than where the drug is, reported per substance as
`pct_at_facility` (amphetamine is 44%, and should be discounted accordingly);
and 3 Catalina deaths that stretch the map frame.
