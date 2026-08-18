# Space-time scan

Findings for `emerging/spacetime.py` (tests: `tests/test_spacetime.py`). Results in [`results/spacetime/`](../../results/spacetime/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## Where, and is it growing there?

`geo.py` asks whether a substance is more concentrated than overdose death
itself, county-wide, over one fixed window. `spacetime.py` asks *where* and
*when*: circles of zipcodes crossed with trailing quarters, returning the
single most unlikely cylinder with a p-value adjusted for every cylinder
searched. 7,572 OD deaths, 280 zips, 12 quarters (2023Q1–2025Q4); circles
capped at 25% of deaths (≤ 97 zips), windows 1–6 quarters, 999 replicates.
46 of 157 substances and tree branches are eligible (≥ 20 cases, and ≤ 25% of
all deaths, the same `MAX_CASE_FRACTION` guard `geo.py` uses).

### Two models, because one cannot answer the question

The **Bernoulli case-control** scan takes cases as deaths naming the substance
and controls as the other overdose deaths in the same zip and quarter, so the
denominator stays overdose death itself — the same convention as EB05,
`treescan.py` and `geo.py`. It measures a substance's local share against its
share pooled over the whole study period, which means **a concentration that
has been there for a decade wins as easily as one that appeared last year.**

That is not a subtle problem. It flagged cocaine at 2.07× over a 5-quarter
window in South LA, which reads as an emerging hotspot and is not one. So the
second model is Kulldorff's **space-time permutation**, conditioned on the
substance's own spatial margin *and* its own temporal margin: a busy zip cannot
win, a busy quarter cannot win, only the interaction can. That is "growing
here" rather than "concentrated here".

### Concentrated: yes. Growing: no.

| | concentrated (Bernoulli) | growing (interaction) |
|---|---|---|
| p ≤ 0.05 | 21 of 46 | 1 of 46 |
| survives Bonferroni across the 46 scanned | **5** | **0** |

The five that survive everything:

| Substance | n | Cluster | Window | Observed | Expected | Ratio |
|---|---|---|---|---|---|---|
| Cocaine | 1,454 | 12 zips around 90011 (South LA + downtown) | 5q | 137 | 66.1 | 2.07× |
| PCP | 367 | 37 zips around 90002 (Watts) | 4q | 68 | 24.2 | 2.81× |
| Dissociatives | 429 | 69 zips around 90274 (South Bay) | 5q | 90 | 40.0 | 2.25× |
| Depressants and Tranquilizers | 481 | 66 zips around 90277 | 5q | 94 | 44.3 | 2.12× |
| Cutting agents and adulterants | 59 | 42 zips around 90043 | 6q | 22 | 6.4 | 3.42× |

**The spatial structure of the LA supply is real and strong, but static —
these substances have a geography, not a moving front.** PCP corroborates the
`geo.py` finding on an independent statistic, and the adulterant family is the
spatial echo of the lidocaine result.

Had only the Bernoulli model been built, this section would claim the opposite.
With trailing windows, "recent elevation" and "long-standing concentration"
produce the same picture unless you condition them apart.

### Multiplicity, again

Each p-value is adjusted for every cylinder searched **within** a substance,
and not at all for having scanned 46 of them. Read per-substance, the
interaction model does have a hit: diphenhydramine, p = 0.003, 4 deaths against
0.44 expected across the affluent Westside (90024, 90049, 90210, 90272 …) over
3 quarters. Corrected across the 46, that is 0.138. It is not a facility
artifact — 11% of its deaths are at shared addresses against 15% county-wide —
and it is not a finding either. `p_across_substances` is a column in
`results/spacetime/clusters.csv` so this does not have to be remembered.

Bonferroni is conservative here, because the scanned set contains nested nodes
(`Dissociatives` contains `PCP`), but it is the statement that cannot be
accused of over-claiming.

### Calibration and the checks that mattered

Both models have synthetic null calibration tests plus a planted-cluster
recovery test in `tests/test_spacetime.py`. On the real background, rejection
rates are 0.075 / 0.113 against nominal 0.05 / 0.10 at C=40, and 0.037 / 0.087
at C=400 — within Monte Carlo error at the thresholds that drive decisions
(80 replicate datasets, so sd ≈ 0.024 and 0.034). The body of the distribution
runs mildly conservative, which errs toward under-claiming.

Two defects found while building, both of the kind that fail silently:

- The observed grid was built with `np.add.at` and the simulated grids with
  `np.bincount`. They happened to agree, but nothing checked it, and a
  mismatch would have made every p-value wrong with no error raised. Now
  asserted before the results are believed.
- The reported member-zip list was truncated to 12, so the map drew a 12-zip
  blob for a 69-zip cluster. The statistics were right; the picture was not.
