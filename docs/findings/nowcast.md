# Reporting delay and completeness

Findings for `emerging/nowcast.py` (tests: `tests/test_nowcast.py`). Results in [`results/nowcast/`](../../results/nowcast/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## How complete are the most recent quarters?

`CUTOFF = 2025-12-31` in `trends.py` is a hand-set date before which toxicology
is declared settled, applied before any statistic sees the data. It was a
judgement with nothing behind it. `nowcast.py` measures the thing it guesses at.

LACME carries no report date, only `DeathDate`, so delay cannot be recovered
from one extract — but it can from two **vintages**. Sweeping death-months at a
single pull date traces the whole curve, because each month sits at a different
delay. Of the three overdose extracts, the 2024-08 / 2025-05 pair is usable
(settled era agrees to 0.2%); the 2025-05 / 2026-04 pair is not (settled era
moves 3.5%, which is a case-definition change, not delay). `triangle` prints
that drift for any pair, so the check is never skipped.

| Months elapsed | 0 | 1 | 2 | 3–5 | 6 | 7–8 | 9 | 10+ |
|---|---|---|---|---|---|---|---|---|
| % reported | 27 | 66 | 76 | ~79 | 85 | ~90 | 99.5 | ~100 |

The saturating end is independently confirmed in the current era: the 2026-02
and 2026-04 pulls agree to within 1.4% on every month of 2025Q1–Q2.

### The correction is computed and then withheld

Applied naively, the curve says 2025Q4 is 81% complete (560 observed → ~690),
2025Q3 is 93%, the cutoff is too late, and the monotone 2025 decline
(681, 628, 602, 560) is really a decline and rebound (681, 629, 650, 690).
That would change the project's headline.

It does not survive its own check. The current extract runs at full strength
through 2026-01 — 198 deaths, in line with every month before it — and then
drops to **54, 18, 6**. The curve predicts a *ramp* across those months
(79%, 77%, 66%); what is there is a **cliff**, a 3.7× fall where 1.2× was
predicted. No smooth accrual process produces that. It is a pull boundary with
a few fast-closing cases trailing it, which means today's reporting is not the
process measured in 2024.

`transfer_check` detects this and `estimate` prints the correction marked
**NOT VALID**, with `curve_transfers=False` in `results/nowcast/quarters.csv`.

**The practical reading inverts.** The extract looks essentially complete
through 2026-01, so `CUTOFF = 2025-12-31` is more likely conservative than lax.
It stays, and nothing downstream changes.

This is the documented failure mode of delay-corrected nowcasting — the CDC's
provisional overdose counts use the same machinery and are known to miss when
the regime shifts. Here the regime that shifted was the reporting process, not
the epidemic.

### What would settle it costs nothing

Two current-era vintages. **Save a dated copy of every future extract.** That
is the entire fix, and without it this question cannot be reopened no matter
how much modelling is applied to it.
