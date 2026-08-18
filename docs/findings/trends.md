# EB05 ranking

Findings for `emerging/trends.py`. Results in [`results/trends/`](../../results/trends/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

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
