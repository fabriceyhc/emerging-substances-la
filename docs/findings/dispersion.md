# Does the Poisson assumption hold?

Findings for `poisson_dispersion` in [`emerging/analysis/trends.py`](../../emerging/analysis/trends.py).
Results in [`results/trends/poisson_dispersion.csv`](../../results/trends/poisson_dispersion.csv).
Back to the [README](../../README.md); the statistics are built up from scratch
in [GAMMA_POISSON_PRIMER.md](../GAMMA_POISSON_PRIMER.md), chronology in
[RESEARCH_LOG.md](../RESEARCH_LOG.md).

## The question

EB05 is a posterior percentile from `Gamma(alpha + n, beta + E)`, and that
closed form is exact only because `n | lambda ~ Poisson(lambda * E)` — which
forces variance to equal mean. Every other module here tests the *inputs* to
that statistic. Nothing tested the assumption underneath it.

The distinction that makes the test possible at all:

- **Between substances**, spread is expected and is the thing being fitted.
  The `Gamma` prior *is* that spread, and the observed table is consequently
  overdispersed by construction. Measured against a plain Poisson the 113
  scored substances give `phi = 2.38`; measured against the negative-binomial
  marginal the model actually claims, `phi = 1.06`. **Overdispersion in the
  cross-section is not evidence of anything wrong.**
- **Within a substance**, spread above Poisson is a real violation, and the
  cross-section cannot see it — one `(n, E)` pair per substance carries no
  information about how that substance's deaths arrived.

The quarterly series can see it. That is what this module does.

## The trap, and why `d_flat` is not the answer

A substance whose rate genuinely moved *also* fails a flat Poisson test,
because `lambda` really did change. That is the signal, not a defect. Testing
against a flat line would condemn precisely the substances the ranking exists
to find.

So the statistic is Pearson dispersion around a **fitted log-linear trend in
the substance's share** (offset by the quarter's all-OD deaths, so the county's
own 22% decline is not read as a substance-level trend). `d_flat` is reported
beside it only to show how much the trend absorbed.

## Result, 2023Q1–2025Q4

| substance | n | `d_flat` | `d_trend` | absorbed | verdict |
|---|---|---|---|---|---|
| para-Fluorofentanyl | 181 | 10.27 | **10.55** | −0.28 | clustered |
| Acetyl fentanyl | 108 | 6.26 | **5.33** | 0.93 | clustered |
| Morphine | 89 | 4.29 | **4.26** | 0.03 | clustered |
| Fentanyl | 4461 | 7.68 | 2.19 | **5.50** | clustered |
| Cocaine | 1464 | 2.07 | 2.01 | 0.05 | clustered |
| Alprazolam | 320 | 2.70 | 2.05 | 0.65 | clustered |
| PCP | 368 | 2.05 | **0.86** | **1.19** | ok |
| Ethanol | 1059 | 0.84 | 0.66 | 0.17 | ok |
| Methamphetamine | 4673 | 0.57 | **0.51** | 0.05 | smooth |

**9 of 21** substances reject the conditional Poisson at p < 0.05.

Read the `absorbed` column against `d_trend`:

**Fentanyl (7.68 → 2.19).** Most of its excess was the trend — a monotone
slide from 509 deaths/quarter to 225. Its rate genuinely fell. A flat test
would have called this the worst violation in the table; it is the second
best-explained.

**PCP (2.05 → 0.86).** *All* of it was trend. Around its rising line PCP is
cleaner than Poisson requires. The top-ranked substance in the ranking is also
the best-behaved one, so its `eb05` of 1.27 is earned rather than borrowed.

**para-Fluorofentanyl (10.27 → 10.55).** Detrending changes nothing, and the
series shows why: `36, 4, 15, 18, 22, 20, 42, 12, 5, 3, 1, 3`. That is not a
rate, it is **batches** — deaths sharing a contaminated supply arrive together.
This is the textbook violation, and the Gamma prior cannot absorb it, because
that layer models differences *between* substances rather than clustering
inside one.

**Methamphetamine (0.57 → 0.51).** Underdispersed — smoother than random. A
violation too, but the benign direction: it makes `eb05` conservative rather
than overconfident.

## What it costs, and what it does not

A `clustered` verdict does not touch the point estimate. `lambda = n/E` is
unaffected. What it invalidates is the *confidence*: the interval should be
roughly `sqrt(d_trend)` times wider.

```
para-Fluorofentanyl  n=  12 E= 66.26   as scored: EB05 0.17 EB95 0.37   corrected: 0.31 / 0.99   x3.5
Morphine             n=  29 E= 23.53   as scored: EB05 0.87 EB95 1.50   corrected: 0.65 / 1.59   x1.5
PCP                  n= 137 E= 90.57   as scored: EB05 1.27 EB95 1.68   corrected: 1.29 / 1.67   x0.9
Methamphetamine      n=1345 E=1304.90  as scored: EB05 0.98 EB95 1.08   corrected: 1.00 / 1.06   x0.7
```

**The current ranking survives this, but by luck rather than design.** The
badly-violated substances (para-fluorofentanyl, acetyl fentanyl) are
*collapsing*, so their overconfident intervals sit at the bottom of the table
where nothing depends on them. The substances near the decision boundary — PCP,
methamphetamine, ethanol — are clean or conservative. That would not hold if a
batch-driven adulterant were on the way *up*, which is exactly the scenario
the tool exists to catch.

## Deliberately not applied

This reports; it does not correct. Dividing `n` and `E` by a per-substance
`d_trend` is the natural next step — `_credit_dispersion` already does the
pooled version of exactly this arithmetic for credit weighting — but it
changes the ranking, so it needs its own `backtest` run against the five known
emergences rather than being switched on inside a diagnostic. Two reasons for
caution beyond the usual:

1. `d_trend` is itself estimated from 12 points, so it is noisy; shrinking it
   toward 1 would probably be wiser than using it raw.
2. A substance can be clustered *because* it is emerging — a new adulterant
   enters via a small number of batches — so a correction risks damping the
   signal at exactly the moment it matters. This wants checking against
   `backtest` before it is trusted.

## Implementation note

The trend fit is Nelder-Mead on the two-parameter Poisson log-likelihood
rather than a GLM, because `statsmodels` is not a declared dependency. Checked
against `statsmodels.GLM(family=Poisson)` across all 21 substances: worst
disagreement `1.4e-07`.

## Usage

```
emerging trends dispersion
emerging trends dispersion --min-total 25
emerging trends rank            # carries a `poisson` column from the same test
```
