# Spatial concentration as a component of the ranking statistic

Findings for `emerging/core/concentration.py` — the `z_s` that
[`GPS_V2_DESIGN.md` §3](../GPS_V2_DESIGN.md) specified and did not have.
Results in [`results/trends/spatial_concentration.csv`](../../results/trends/spatial_concentration.csv);
the head-to-head that decides whether it earns its place is
[benchmark.md](benchmark.md). Back to the [README](../../README.md);
chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## The statistic, and why it is not a permutation test

Every other spatial analysis here (`geo.py`, `spacetime.py`, `relrisk.py`)
buys its p-value with 999 random-labelling replicates. That is affordable once
per substance. It is not affordable **45 times per substance**, which is what
feeding an as-of sweep requires, and `spacetime.py` additionally refuses to
run below `MIN_CASES = 20` — excluding every substance this project exists to
find.

So the score is a Cuzick–Edwards-style case-control clustering count with
*closed-form* null moments instead:

    S_s = #{ordered pairs of s-deaths that are spatial neighbours}

over the symmetrized k-nearest-neighbour graph on all overdose deaths in the
trailing window, k = 5% of the window. Under random labelling — which C of
these N deaths name the substance, the same null `geo.py` defends — `S_s` has
exact mean and variance, so `z_s = (S_s − E[S_s]) / sd[S_s]` costs one sparse
mat-vec rather than a Monte Carlo run. The full 45-quarter × 211-substance
sweep takes **4.6 seconds**.

**The exactness is the load-bearing claim, so it is tested against the null it
describes**, not against intuition. `tests/test_concentration.py::test_moments_match_simulation`
draws 4,000 replicates at C ∈ {5, 10, 30, 100} and checks mean and sd:

| C | sim mean | exact mean | sim sd | exact sd | P(z > 1.645) |
|---|---|---|---|---|---|
| 5 | 1.29 | 1.25 | 1.53 | 1.54 | 0.106 |
| 10 | 5.56 | 5.63 | 3.21 | 3.28 | 0.061 |
| 20 | 23.73 | 23.78 | 6.78 | 6.80 | 0.054 |
| 50 | 153.33 | 153.29 | 18.05 | 17.78 | 0.059 |
| 150 | 1399.74 | 1398.38 | 57.40 | 58.24 | 0.057 |
| 400 | 9988.79 | 9985.74 | 179.57 | 177.53 | 0.063 |

The moments are right everywhere. **The normal approximation to them is not,
below about C = 10** — at C = 5 the nominal 5% tail is 10.6%, because the
statistic is skewed and lumpy at low counts. That table is what sets
`MIN_CASES = 10`, rather than a round number chosen for looking careful.

## Two designs measured and rejected before this one

Both were built and run on the real data, not reasoned about:

- **Zipcode cells** (neighbour ≡ same zip). Cheapest, and it matches
  `spacetime.py`'s geography. Background co-location probability is
  Q = 0.009, so a 6-death substance expects 0.28 co-located pairs and scores
  S = 0 about three quarters of the time. On the 2025 window it separated
  nothing below C ≈ 50. **No resolution in the regime that matters.**
- **Fixed radius.** Fixes the resolution and breaks something else: at any one
  radius, Skid Row deaths have enormous degree and Antelope Valley deaths have
  almost none, so the statistic's power goes to "is this substance downtown"
  rather than to concentration. Still a *valid* test — the null absorbs the
  density structure — but it spends its power on the wrong contrast.

k-NN at a fixed *fraction* of the window is density-adaptive and scale-free in
the same sense `spacetime.MAX_SPATIAL_FRACTION` is: "the nearest 5% of the
county's overdose deaths" means the same thing in a 300-death quarter and an
800-death one.

## What it finds

It reproduces the geography the pipeline already knows about, on an axis
nothing else in the ranking sees. On the final window, PCP scores z = 3.84 and
cocaine z = 5.19 — PCP's Skid Row concentration is the finding `geo.py` and
`relrisk.py` both report independently, arrived at here from a different
statistic on a different window.

More usefully, it tracks the **onset** of the labelled emergences:

| Substance | as-of | n in window | S | E[S] | z |
|---|---|---|---|---|---|
| para-Fluorofentanyl | 2021Q2 | 24 | 60 | 34.1 | **3.18** |
| para-Fluorofentanyl | 2021Q3 | 41 | 134 | 101.2 | 2.29 |
| para-Fluorofentanyl | 2022Q3 | 109 | 658 | 734.7 | −1.83 |
| Mitragynine | 2023Q1 | 21 | 60 | 26.5 | **4.63** |
| Carfentanil | 2024Q3 | 12 | 22 | 8.2 | **3.47** |

**para-Fluorofentanyl's trajectory is the signature the design was after and
did not know it would get.** It enters *concentrated* (z = 3.18 in 2021Q2, the
quarter after its first LA death), and the concentration decays to
*negative* by 2022Q3 as it saturates the county-wide fentanyl supply. A new
substance appearing in a few places and then diffusing is exactly the shape a
point-source contamination hypothesis predicts, and it is visible here as a
number rather than by eye on a map.

## The limit, stated plainly

**Three of the five labelled emergences never clear `MIN_CASES` at all.**
Xylazine's z is identically 0 in every one of the 45 quarters — it has 9
lifetime mentions and never puts 10 of them in a 4-quarter window. Bromazolam
is scored in 6 quarters, carfentanil in 5; only para-fluorofentanyl (19) and
mitragynine (13) are scored through their emergence.

So this component **cannot help the substances the project most wants help
with**, and no tuning fixes that: it is the same signal ceiling
`RESEARCH_LOG.md` records for the space-time scan and for P1's aggregation
argument. Splitting 9 deaths across 290 zips does not create information. What
the component *can* do is sharpen the two substances that do clear the floor,
and [benchmark.md](benchmark.md) measures whether that is worth anything.

## Caveats

- **Facility clustering is priced in for the background, not for a
  differential.** Deaths at hospitals and jails share coordinates
  (`core/spatial.py`'s `FACILITY_MIN`), which inflates degree for cases and
  controls alike, and the random-labelling null handles that. What it cannot
  handle is a substance *differentially* concentrated at facilities —
  amphetamine is 44% at-facility per `geo.py`. Read a high z for such a
  substance as "concentrated where it is pronounced", not "concentrated where
  it is used".
- **`z_s` is not a p-value.** `concentration_weight` squashes it with a normal
  CDF, and above the `MIN_CASES` floor that is a reasonable monotone
  transform, but reading `2Φ(z) − 1 = 0.95` as 95% confidence would be
  overclaiming at these counts. It is used as a bounded weight, which is all
  the discount needs.
- **It reads geography, not novelty.** A substance with a decade-old stable
  cluster scores high every quarter — PCP does. The trailing window keeps
  `geo.py`'s lifetime ratio from leaking in, but it does not by itself
  distinguish a *new* cluster from an old one. The temporal half of the
  statistic is what is supposed to do that, which is the argument for fusing
  them rather than reading the concentration score alone.
