# Tree-temporal scan

Findings for `emerging/treescan.py` (tests: `tests/test_treescan.py`). Results in [`results/treescan/`](../../results/treescan/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## Scanning the hierarchy instead of one substance at a time

`treescan.py` implements a tree-temporal scan statistic over the substance
hierarchy in `tree.py`: every node — leaf, pharmacological family, NFLIS
category, superclass, root — against every trailing window of 1–6 quarters
inside the same 12-quarter frame EB05 uses. The p-value comes from the Monte
Carlo distribution of the **maximum** log-likelihood ratio over the whole
search, so it is already adjusted for all 158 nodes × 6 windows. Reported as a
recurrence interval: RI 100 means a signal this strong turns up in one scan in
100 by chance.

It removes the two choices the section above just showed to be load-bearing —
which window, and which of 205 substances you looked at — and the multiplicity
correction charges for having made them.

**It changes the current answer.** EB05 says nothing clears 1.5. The scan,
after full adjustment, says four nodes do:

| Node | Level | Window | Observed | Expected | RI |
|---|---|---|---|---|---|
| Lidocaine | substance | 6q | 31 | 17.4 | 1,111 |
| Dissociatives | family | 5q | 203 | 158.3 | 455 |
| Depressants and Tranquilizers | category | 5q | 225 | 178.4 | 270 |
| PCP | substance | 5q | 174 | 134.5 | 185 |

Lidocaine's burst is real; EB05 missed it only because the burst straddles the
4-quarter boundary. PCP corroborates the spatial finding from `geo.py` on an
independent axis. Both hold under either reference series (`--reference
deaths`, the EB05 analogue, and `--reference mentions`, which is immune to
panel drift) — and cocaine, the EB05 near-miss that `count_ratio` exposed as a
growing share of a shrinking total, is nowhere near significant under either.

### But it does not detect earlier than EB05, and the aggregation argument failed

`treescan backtest` runs the same as-of sweep as `trends backtest` on the same
ground truth, so the only thing that differs is the statistic. EB05 > 1.5 is
unadjusted while RI ≥ 100 is adjusted for the entire search, so the sweep also
re-cuts the scan at the threshold that raises **the same number of alarms** —
85 over the sweep — and that matched-budget column is the fair comparison.

| Substance | First seen | EB05 | Leaf, RI ≥ 100 | Leaf @ matched budget |
|---|---|---|---|---|
| para-Fluorofentanyl | 2021Q1 | +1q | **+0q** | +1q |
| Bromazolam | 2021Q2 | **+11q** | +11q | +12q |
| Carfentanil | 2017Q2 | +29q | +29q | +29q |
| Mitragynine | 2017Q4 | **+20q** | never (peak RI 48) | never |
| Xylazine | 2023Q3 | **+4q** | never (peak RI 13) | never |

At a matched false-alarm budget the tree scan never detects earlier, and it
misses two of five outright. **EB05's apparent advantage on mitragynine and
xylazine is bought with an unadjusted threshold** — but the scan does not buy
anything back with its adjustment either.

Two further negative results worth recording:

- **The aggregation-power argument, which is what motivated the whole thing,
  did not survive contact with the data.** The nitazene branch never signals:
  peak LLR 3.38, RI 1. Aggregation helps when members are individually
  *marginal*; the nitazenes are individually *absent* (1 and 2 mentions in
  fourteen years), and no amount of pooling creates signal from three cases.
  Across the whole sweep, branch-only detection — a branch signalling in a
  quarter when no member leaf did — happens in **3 node-quarters out of 165**.
- **Branch signals are not attributable by default.** `Fentanyl and
  Fentanyl-related` fires in 2017Q2, the exact quarter carfentanil first
  appears in LA, which scores as a 0-quarter detection until you ask how much
  of the branch's excess carfentanil actually contributes: **1%**. It was
  fentanyl's own explosion. `excess_share()` computes this and the head-to-head
  reports it, because without it a large branch containing a rare substance
  manufactures detections. Only para-fluorofentanyl gets a genuine branch
  detection (74% attributable, via `Fentanyl analogs`).

**Verdict: it runs alongside EB05, it does not replace it.** What it adds is
multiplicity-honest inference and window-free detection, which is exactly what
promoted lidocaine and PCP from "below the line" to "signal". What it does not
add is earlier warning.

### Reading the two together

EB05 *estimates* and the scan *tests*, so they answer different questions and
only EB05 can rank (p-values tie at 1.0 for half the catalog) or detect
declines (the scan is high-rate only). EB05 screens, the scan adjudicates:

| | RI ≥ 100 | RI < 100 |
|---|---|---|
| **EB05 > 1.5** | a real signal, defensible to a reviewer | a candidate to watch, not publish — where EB05's unadjusted threshold was doing the work (xylazine, mitragynine) |
| **EB05 < 1.5** | EB05's fixed window missed it (lidocaine, PCP) | nothing |

Both off-diagonal cells are new information. The bottom-left one is the
current headline: `trends rank` says nothing clears 1.5, and lidocaine and PCP
contradict it.

### The tree is where the findings are made or lost

The game plan called for `root → NFLIS category → substance`, on the grounds
that the category field is already in the data. Building it showed the NFLIS
categories are a forensic-chemistry filing system, not a pharmacological one:

- `Depressants and Tranquilizers` holds PCP (724 mentions) next to xylazine (9)
  and zolpidem — PCP is 60% of the node, so any tranq signal drowns.
- `Narcotic Analgesics` mixes morphine (1,124) with the nitazenes (3 total).
  The family that motivated a tree scan is 0.1% of its own parent, so
  aggregating to the NFLIS parent *loses* the signal rather than pooling it.
- `Other substances` is a 71-member junk drawer holding lidocaine and
  levamisole alongside atorvastatin and vancomycin.

So `tree.py` keeps the NFLIS categories — they are the published taxonomy and a
reader can check them — and adds a family layer beside them. Families cross
category boundaries, which makes the structure a **DAG rather than a tree**;
nothing in the statistic requires otherwise, since the LLR is computed per node
and the Monte Carlo maximum absorbs however much the nodes overlap. Nodes whose
leaf set duplicates another's are dropped, keeping the more specific label.

`emerging tree check` lists every family member absent from LA, so
the constant stays a claim about pharmacology rather than a wish list. Absent
members are kept on purpose: medetomidine is in `Veterinary sedatives` today,
so the branch already exists the quarter it is first coded.

### Validated against the TreeScan binary

`treescan.py` is a from-scratch implementation, which is a liability until it
is checked against the reference. `emerging treescan-validate compare`
runs both on the same tree and the same counts and reports every difference.

| | strict tree | DAG (production) |
|---|---|---|
| Nodes in the search space | 141 both | **155 both** |
| Cuts compared exactly | 57 | **65** |
| Max abs LLR difference | 5.0e-07 | **5.0e-07** |
| Best window identical | 57/57 | **65/65** |
| Cases in window identical | 57/57 | **65/65** |
| Max abs p-value difference (9,999 reps each) | 0.0041 | **0.0039** |
| Same verdict at p = 0.01 | all | all |

5e-07 is TreeScan's printf precision — it writes 6 decimal places, so the LLRs
are identical to the last digit either program prints. The p-values are two
independent Monte Carlo runs and agree to within simulation error.

Getting there required matching four settings whose TreeScan defaults differ
from ours, and each one is a way the comparison could have looked like a
disagreement while being nothing of the kind:

- `conditional-type=2` (condition on the node), not the shipped example's
  `NODEANDTIME`.
- `apply-risk-window-restriction=n`. It defaults to *on* at 20%, silently
  dropping short windows.
- `include-identical-parent-cuts=y`, so TreeScan reports the duplicate
  leaf-set nodes our `_dedupe` removes.
- `minimum-node-cases=2`. TreeScan will not form a cut from fewer than 2
  cases — in the node or in the window. This is a property of the *search
  space*, not of reporting, so leaving ours at 0 would have had our null
  maximum range over nodes TreeScan never sees.

Two things the comparison found:

- **A real bug in ours.** Metoprolol showed 3 observed against 3.000 expected
  and an LLR of 3.3e-16 — `c > e` came out true on the last bit of a float. It
  could never win a maximum, but it put a node with no excess into the
  positive-LLR list. `_llr` now floors at 1e-9.
- Substance names break TreeScan's input format. It is comma-delimited with no
  quoting, so `1,1-Difluoroethane` parses as an extra column; the first run
  died with "record 17 has node referencing self as parent". The exporter now
  writes opaque IDs and maps names back.

Also worth recording: the export is checked before it is used. `verify_export`
recomputes every node's descendant-leaf closure *from the emitted file* and
requires it to equal our leaf set, so a mistranslated DAG fails loudly instead
of quietly comparing two different trees.

**Performance** (155 nodes, 125 leaves, 12 quarters, single-threaded):

| Replications | TreeScan (C++) | ours (Python/numpy) |
|---|---|---|
| 999 | 2.13 s | 0.43 s |
| 9,999 | 2.25 s | 3.5 s |
| 99,999 | 12.4 s | 34 s |

Comparable, with different shapes: TreeScan carries a ~2 s constant and then
scales well, ours is near-linear at ~0.35 ms per replication. TreeScan also
uses all cores (99,999 replications in 2.5 s with `--parallel-processes 0`)
where ours is single-threaded. At the replication counts this project needs,
neither is a constraint.

### Calibration

The scan produces p-values, so the property that matters is that they mean what
they say. On null data — every leaf's cases distributed over time exactly as
the reference series says — the observed rejection rates are 0.050, 0.096,
0.183 and 0.496 against nominal 0.05, 0.10, 0.20 and 0.50. `tests/
test_treescan.py::test_null_calibration` asserts this, along with the specific
way it would break quietly: applying the `min_cases` floor to the observed data
but not to the simulations.

`tests/test_treescan_parity.py` re-runs the binary comparison as a test. It
skips when the binary is absent, so it is a no-op for anyone who only has this
repository, and it is the guard against our implementation drifting away from
the reference on some later edit.

**One thing the recurrence interval does not cover: repeated looks.** It is the
multiplicity across nodes and windows *within one analysis*. Running the scan
every quarter is a further multiplicity, and the game plan's claim that
TreeScan's RI makes MaxSPRT redundant is only true of the within-analysis part.
