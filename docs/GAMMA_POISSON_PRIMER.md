# Gamma–Poisson shrinkage, from the ground up

A primer on the statistics behind `trends.py`'s EB05 ranking, assuming no
background beyond "I know what an average is."

Every numbered block below is a **complete, runnable script**. Each seeds its
own random number generator, so pasting any one block on its own reproduces the
output printed under it exactly. Requirements are what the repo already has:
`numpy`, `scipy`.

The destination is one line of code — `trends.py:358`:

```python
post_a = pr["alpha"] + n
post_b = pr["beta"] + e
...
"eb05": gamma_dist.ppf(0.05, post_a, scale=1.0 / post_b),
```

By the end you should be able to say what every symbol in it is, where it came
from, and what it promises. Twelve steps.

**Contents**

- [0. The problem this solves](#0-the-problem-this-solves)
- [1. Mean and variance](#1-mean-and-variance)
- [2. The Poisson distribution](#2-the-poisson-distribution)
- [3. Variance = mean, and why that's the whole game](#3-variance-mean-and-why-thats-the-whole-game)
- [4. Rates, expected counts, and why `n/E` fails](#4-rates-expected-counts-and-why-ne-fails)
- [5. Bayes' rule](#5-bayes-rule)
- [6. The Gamma distribution](#6-the-gamma-distribution)
- [7. Conjugacy: the trick that makes it fast](#7-conjugacy-the-trick-that-makes-it-fast)
- [8. Shrinkage, seen directly](#8-shrinkage-seen-directly)
- [9. Where does the prior come from? The negative binomial](#9-where-does-the-prior-come-from-the-negative-binomial)
- [10. Empirical Bayes: fitting the prior from the data](#10-empirical-bayes-fitting-the-prior-from-the-data)
- [11. EBGM and EB05: what the ranking promises](#11-ebgm-and-eb05-what-the-ranking-promises)
- [12. Why weighted counts break it, and the φ fix](#12-why-weighted-counts-break-it-and-the-φ-fix)
- [Cheat sheet](#cheat-sheet)

---

## 0. The problem this solves

You have 203 substances, and you want to rank them by "which is rising."

For each one you have two numbers. **`n`** is the observed count: how many
overdose deaths named it in the recent window. **`E`** is the *expected* count,
and it is worth being concrete about where it comes from, because everything
downstream is a comparison against it.

`rank_substances` splits the settled quarters into two windows — a **recent**
one (2025Q1–2025Q4, holding 2,153 OD deaths) and the **baseline** before it
(2023Q1–2024Q4, holding 5,491). Then, at `trends.py:414`:

```python
expected = n_base * (n_recent_all / n_base_all)
```

`E` is the substance's *own baseline count*, scaled down by how much smaller the
recent window is. In words: **if this substance had held its baseline share of
all OD deaths, this is how many we'd have seen.** Nothing external is assumed —
`E` is the substance's own past, restated at the recent window's size.

The scale factor here is 2153/5491 = **0.392096**, the same for every substance.
So:

```
substance             n_base   x 0.392096   E     n     n/E
Nitrous oxide              1       0.3921  0.39    4   10.20
Chlordiazepoxide           5       1.9605  1.96    6    3.06
PCP                      231      90.5742 90.57  137    1.51
```

One subtlety to note now and pick up in step 12: that factor is 0.392, not 0.5,
even though the recent window is 4 quarters against a baseline of 8. Those would
match if total OD deaths had held steady; they fell ~22% instead. So this `E`
measures **share** — "rising" means taking a larger slice of a shrinking pie.
(`trends.py:433` computes a second expectation, `expected_own = n_base * per_q`,
which uses the 0.5 and asks about the substance's own count instead. Both are run
through the machinery below; that's the `eb05` / `eb05_own` dual gate.)

The obvious statistic is now the ratio `n/E`. It does not work, and the reason it
does not work is not subtle:

| substance | n | E | n/E | *(EB05)* |
|---|---|---|---|---|
| Nitrous oxide | 4 | 0.39 | **10.20** | *0.80* |
| Chlordiazepoxide | 6 | 1.96 | 3.06 | *0.83* |
| PCP | 137 | 90.57 | 1.51 | ***1.27*** |

(Real rows from `results/trends/rising_substances.csv`. Ignore the last column
for now — it's where we're going.)

Note what `E` is made of in each row. Nitrous oxide's 0.39 is *one death across
two years*, scaled down. PCP's 90.57 is 231 deaths, scaled down. Identical
arithmetic, wildly different amounts of evidence underneath — and the ratio
`n/E` shows no trace of that difference.

Nitrous oxide went up tenfold. It went from 1 death to 4. Four deaths in a county
with 2,153 overdose deaths in the window is a number you would get from noise
alone, routinely. Meanwhile PCP's unremarkable-looking 1.51 rests on 137 deaths
and is the most solid signal in the table.

Ranking by `n/E` puts the noise on top, every time. Everything that follows is
machinery for fixing that in a principled way — one that *automatically* knows
that 4-vs-0.39 is weak evidence and 137-vs-90.6 is strong evidence, without you
hand-tuning a minimum-count threshold.

---

## 1. Mean and variance

Two numbers summarize any batch of random outcomes.

The **mean** is the average — where the outcomes sit.

The **variance** is the average of *squared* distances from the mean — how
spread out they are. Squaring is what makes it work: distances above and below
the mean would otherwise cancel to zero.

The **standard deviation** (`sd`) is the square root of the variance, which puts
it back into the original units. Read it as "the typical distance an outcome
lands from the mean."

```python
# --- 1 ---
import numpy as np
rng = np.random.default_rng(0)

rolls = rng.integers(1, 7, size=100_000)
print(f"mean     {rolls.mean():.3f}")
print(f"variance {rolls.var():.3f}      = mean of (roll - mean)^2")
print(f"std dev  {rolls.std():.3f}      = typical distance from the mean")
print(f"check    {((rolls - rolls.mean())**2).mean():.3f}")
```

```
mean     3.498
variance 2.910      = mean of (roll - mean)^2
std dev  1.706      = typical distance from the mean
check    2.910
```

A die averages 3.5 and typically lands about 1.7 away from that. The last line
is just the definition of variance computed by hand, to show `np.var` isn't
doing anything mysterious.

**Why you care:** "is this substance rising" is the question "is this count
further above its expectation than the spread would ordinarily produce." You
cannot answer it without knowing the spread. The rest of this document is
largely about where the spread comes from.

---

## 2. The Poisson distribution

A **distribution** is a full accounting of how likely each outcome is — not just
the average, the whole shape.

Consider how a substance count is actually generated. There were about 30,000
deaths in the cohort. Each one either named substance X or didn't. Suppose each
death independently has some small probability *p* of naming it.

Counting successes across N independent tries is the **binomial** distribution.
But when N is huge and p is tiny, something convenient happens: the binomial
stops caring about N and p separately, and depends only on their product
`μ = N·p`. That limit is the **Poisson distribution**.

```python
# --- 2 ---
import numpy as np
from scipy.stats import binom, poisson
rng = np.random.default_rng(0)

N, p = 30_000, 0.0005            # 30k deaths, each 0.05% likely to name substance X
b = rng.binomial(N, p, 200_000)  # exact: count them one death at a time
q = rng.poisson(N * p, 200_000)  # approximation: only the product N*p matters

print(f"binomial(30000, 0.0005)  mean {b.mean():6.3f}   var {b.var():6.3f}")
print(f"poisson(15)              mean {q.mean():6.3f}   var {q.var():6.3f}\n")
ks = range(11, 20)
print("k          ", "  ".join(f"{k:>6}" for k in ks))
print("binomial   ", "  ".join(f"{binom.pmf(k, N, p):6.4f}" for k in ks))
print("poisson    ", "  ".join(f"{poisson.pmf(k, N*p):6.4f}" for k in ks))
```

```
binomial(30000, 0.0005)  mean 14.989   var 14.942
poisson(15)              mean 15.007   var 15.059

k               11      12      13      14      15      16      17      18      19
binomial    0.0663  0.0829  0.0956  0.1025  0.1025  0.0961  0.0848  0.0706  0.0557
poisson     0.0663  0.0829  0.0956  0.1024  0.1024  0.0960  0.0847  0.0706  0.0557
```

Identical to four decimal places. This is why counts of rare events —
radioactive decays, typos per page, deaths naming a particular substance — are
modelled as Poisson essentially by default. It isn't an arbitrary choice; it's
what "many independent chances, each rare" *becomes*.

The formula, for reference:

```
P(n | μ) = μⁿ · e^(−μ) / n!
```

Note the shape `μⁿ · e^(−μ)`. Hold on to it — step 7 depends entirely on it.

---

## 3. Variance = mean, and why that's the whole game

The Poisson has exactly one parameter, μ. That single number fixes both the
mean *and* the variance:

```
mean = μ        variance = μ        sd = √μ
```

```python
# --- 3 ---
import numpy as np
rng = np.random.default_rng(0)

print(" mean   sample variance   sample sd   sqrt(mean)")
for m in [1, 5, 15, 100, 400]:
    s = rng.poisson(m, 500_000)
    print(f"{m:5d} {s.var():16.2f} {s.std():11.2f} {np.sqrt(m):12.2f}")
```

```
 mean   sample variance   sample sd   sqrt(mean)
    1             1.00        1.00         1.00
    5             4.99        2.23         2.24
   15            15.03        3.88         3.87
  100            99.95       10.00        10.00
  400           399.32       19.98        20.00
```

**This is the most important fact in the document.** Under a Poisson model,
spread is not a free parameter you estimate — it is *implied* by the level. If
you expect 100 deaths, you already know that 90 and 110 are unremarkable
(sd = 10) and 140 is not. Nobody had to measure that.

That is an enormous amount of free inferential leverage, and it is the reason
the whole closed-form apparatus exists. It is also exactly what breaks in step
12, which is why step 12 is a problem at all.

The ratio `variance / mean` gets a name, since we'll need it later:

```
φ  =  variance / mean       (the "variance-to-mean ratio", or dispersion)
```

Poisson means φ = 1. φ > 1 is called **overdispersed** — more spread than the
level justifies. φ < 1 is **underdispersed**.

### "But my data isn't φ = 1 — is Poisson dead?"

Almost certainly not, and this is the most common wrong turn. Measure the real
table against a plain Poisson and it looks disastrous:

```
substances with E > 0 : 113
dispersion vs plain Poisson(E)      phi =     2.38   <- looks catastrophic
dispersion vs the NB marginal       phi =     1.06   <- what the model predicts
```

The resolution is that the model is a *stack*, and only one layer needs φ = 1:

```
lambda_s ~ Gamma(alpha, beta)     between substances -- spread lives HERE, by design
n_s | lambda_s ~ Poisson(lambda*E)  within a substance -- φ = 1 required HERE
-----------------------------------------------------------------------------
n_s ~ NegativeBinomial            what you observe   -- φ > 1 ALWAYS
```

The Poisson claim is *conditional*: given a substance's own true rate, its count
is Poisson. Nothing claims the observed table is Poisson — the model insists it
is not, because substances really do have different rates. That excess width is
what step 9 derives as `(1 + E/beta)` and what step 10 fits. **The overdispersion
you can see is the signal, not a violation.**

So the question is never "is the table overdispersed" (it is, always). It is
"is the *conditional* overdispersed" — and step 12 is what happens when it is.

**The catch: you cannot answer that by testing the data.** Two worlds, 203
substances each — A where the model is correct and substances truly differ, B
where every lambda is exactly 1 but counts arrive in clumps of {0,1,2}:

```
A  heterogeneous, Poisson conditional   mean 10.03  var 16.70  var/mean 1.67
B  homogeneous, clumped conditional     mean 10.02  var 16.85  var/mean 1.68

A (some substances really ARE elevated)  flagged  32.9/203   truly elevated  15.0
B (nothing is elevated -- all false)     flagged  34.2/203   truly elevated   0.0
```

Indistinguishable marginals, indistinguishable flag counts, opposite truths. In
B every flag is false. This is why step 12's correction is computed from the
credit scheme's own definition (`E[c^2]/E[c]`, known before seeing any data)
rather than estimated from residuals — whether the conditional holds is a fact
about the generating mechanism, not something the table can be asked.

The practical test is mechanical, not statistical:

> Overdispersion **between** substances is what GPS is for.
> Overdispersion **within** a substance is what breaks it.

Ask whether two deaths naming a substance arrive independently. If one
contaminated batch kills five people, they don't — and no Gamma prior repairs
that, because the Gamma layer models differences *across* substances, not
clustering *inside* one. That is a real unmodelled exposure here: a localized
outbreak is exactly such a cluster, which is part of why `spacetime.py` and
`relrisk.py` exist as separate detectors.

And if the conditional genuinely is off:

| situation | cost of ignoring | fix |
|---|---|---|
| φ > 1 (clumping, contagion, shared supply) | intervals too narrow → false positives | divide n and E by φ (step 12); or NB / quasi-Poisson |
| φ < 1 (bounded or regularized counts) | intervals too wide → missed signals | NB **cannot** fit it (var = μ + μ²/α ≥ μ); needs Conway–Maxwell-Poisson |
| φ ≈ 1.1–1.3 | mild | usually leave it; the model is an approximation regardless |

Note the asymmetry: overdispersion costs false alarms, underdispersion costs
sensitivity. For a surveillance tool that has to defend every flag, erring
underdispersed is the cheaper mistake.

### What is lambda, exactly?

`lambda` is the **true multiplier on E** — one fixed but unknown number per
substance per window. For PCP, `E = 90.57` is what you'd see if PCP held its
2023-24 share; `lambda = 1.5` means its true 2025 share was 1.5x that. The
observed 137 is one noisy draw from `Poisson(lambda * 90.57)`.

The counterfactual makes it concrete: **imagine re-running 2025 a thousand
times** — same population, same supply, same underlying rate. `lambda` stays
pinned wherever it is; the counts scatter around `lambda * E` with variance
`lambda * E`. You observed one of those draws.

What's easy to get backwards: **lambda is not wobbling.** The Gamma over it
expresses *your uncertainty about a fixed number*, not variability in the number
itself. PCP's posterior `Gamma(144.3, 98.1)` says "lambda is one specific value
and I'm 90% sure it lies in 1.27-1.68."

### So what does a violation actually look like?

Variance far above the mean — but there is a trap. **A substance that is
genuinely rising or falling also produces variance >> mean** when you pool its
quarters, because `lambda` really moved. That is the signal, not a defect. The
test has to be run *around a fitted trend*.

Do that on the real quarterly counts, 2023Q1-2025Q4:

```
                      23Q1   23Q2   23Q3   23Q4   24Q1   24Q2   24Q3   24Q4   25Q1   25Q2   25Q3   25Q4
Fentanyl              509    479    494    438    445    388    273    334    311    293    272    225
para-Fluorofentanyl   36     4      15     18     22     20     42     12     5      3      1      3
Methamphetamine       444    439    463    388    432    390    376    396    367    341    331    306
Morphine              2      7      6      8      9      6      20     2      11     12     6      0
```

```
            raw D/df   D/df after removing a log-linear trend
Fentanyl                 7.68              2.19
para-Fluorofentanyl     10.27             10.55
Acetyl fentanyl          6.26              5.33
Morphine                 4.29              4.26
Cocaine                  2.07              2.01
PCP                      2.05              0.86
Methamphetamine          0.57              0.51
```

`D/df = 1` is exactly Poisson. The gap between the columns is the whole story:

- **Fentanyl 7.68 -> 2.19.** Most of the excess was the trend — a clean monotone
  slide from 509 to 225. `lambda` genuinely fell. Not a model failure.
- **PCP 2.05 -> 0.86.** *All* of it was trend. Around its rising line PCP is
  textbook Poisson — the top-ranked substance is also the best-behaved one, so
  its EB05 of 1.27 is earned.
- **para-Fluorofentanyl 10.27 -> 10.55.** Detrending changes nothing, and the
  series shows why: `36, 4, ... 42, 12, 5, 3, 1, 3`. That is not a rate, it is
  **batches** — deaths sharing a contaminated supply arrive together, which is
  precisely the non-independence the conditional Poisson forbids. The textbook
  violation.
- **Methamphetamine 0.57 -> 0.51.** The other direction: *under*dispersed,
  smoother than random.

The cost, since interval width scales as `1/sqrt(effective n)`:

```
para-Fluorofentanyl  n=  12 E= 66.26   phi=1: EB05 0.17 EB95 0.37  |  phi=10.55: EB05 0.31 EB95 0.99   x3.5
Morphine             n=  29 E= 23.53   phi=1: EB05 0.87 EB95 1.50  |  phi= 4.26: EB05 0.65 EB95 1.59   x1.5
PCP                  n= 137 E= 90.57   phi=1: EB05 1.27 EB95 1.68  |  phi= 0.86: EB05 1.29 EB95 1.67   x0.9
Methamphetamine      n=1345 E=1304.90  phi=1: EB05 0.98 EB95 1.08  |  phi= 0.51: EB05 1.00 EB95 1.06   x0.7
```

Reassuring for the current table: the badly-violated substances are *collapsing*
rather than rising, so the overconfidence does not corrupt the top of the
ranking. The substances near the decision boundary are clean or conservative.

**This refines the "you cannot test it" claim above.** That holds for the
*cross-section* — one `(n, E)` pair per substance, which is all
`rank_substances` sees. The *time series* does carry partial information, as the
table shows. It is not full identification (a `lambda` wiggling erratically for
real reasons still mimics clustering), but it is cheap and real.

This is now built: `emerging trends dispersion` runs the test above, and
`emerging trends rank` carries a `poisson` column (`ok` / `clustered` /
`smooth`) from it. It reports without correcting — see
[`docs/findings/dispersion.md`](findings/dispersion.md) for why applying a
per-substance correction needs its own backtest first.

---

## 4. Rates, expected counts, and why `n/E` fails

In `trends.py`, for each substance:

- **`n`** — observed count in the recent window (`n_recent`).
- **`E`** — the *expected* count: the substance's baseline count rescaled to the
  recent window's size. `rank_substances` computes it as
  `n_base * (n_recent_all / n_base_all)` — "if this substance had held its
  baseline share of all OD deaths, this is how many we'd have seen."
- **`λ`** (lambda) — the thing we actually want: the true multiplier. λ = 1
  means "exactly on baseline," λ = 2 means "twice baseline."

The model is `n ~ Poisson(λ·E)`. The obvious estimate of λ is `n/E`. Watch it
fail:

```python
# --- 4 ---
import numpy as np
rng = np.random.default_rng(0)

print("     E    5th pct   median   95th pct   sd of n/E   1/sqrt(E)")
for E in [0.4, 4, 40, 400]:
    r = rng.poisson(E, 500_000) / E          # truth is lambda = 1 in every row
    print(f"{E:6.1f} {np.percentile(r,5):10.2f} {np.percentile(r,50):8.2f}"
          f" {np.percentile(r,95):10.2f} {r.std():11.2f} {1/np.sqrt(E):11.2f}")
```

```
     E    5th pct   median   95th pct   sd of n/E   1/sqrt(E)
   0.4       0.00     0.00       5.00        1.58        1.58
   4.0       0.25     1.00       2.00        0.50        0.50
  40.0       0.75     1.00       1.27        0.16        0.16
 400.0       0.92     1.00       1.08        0.05        0.05
```

**Every row has a true λ of exactly 1.** Nothing is rising anywhere.

Yet at E = 0.4, one run in twenty produces `n/E` of 5.00 or more. At E = 400,
the same 95th percentile is 1.08. The noise in `n/E` scales as `1/√E` — the
estimator is 32× noisier for the small substance than the large one.

So a table ranked by `n/E` is not a ranking of "how much is this rising." It is
very close to a ranking of "how small is this substance," with the genuinely
large risers buried underneath. Nitrous oxide's 10.20 from step 0 is this table's
first row happening in real life.

**The fix has to do two things at once:** pull unreliable estimates toward
something sensible, and do it *more* for small substances than large ones. That
is shrinkage, and Bayes' rule delivers it without being asked.

---

## 5. Bayes' rule

Bayes' rule says how to update a belief when evidence arrives:

```
posterior  ∝  prior × likelihood
```

- **prior** — what you believed before seeing this substance's data.
- **likelihood** — how well each candidate explanation predicts the data you saw.
  Note this is *not* a probability of the explanation; it's the probability of
  the *data*, computed under that explanation.
- **posterior** — the updated belief.

The `∝` ("proportional to") hides a division by a normalizing constant that
makes everything sum to 1. Concretely, with just two candidate explanations:

```python
# --- 5 ---
from scipy.stats import poisson

E, n = 2.0, 5                      # expected 2 deaths, observed 5
hyp   = {"flat  (lambda=1)": 1.0, "rising(lambda=3)": 3.0}
prior = {k: 0.5 for k in hyp}

unnorm = {k: prior[k] * poisson.pmf(n, lam * E) for k, lam in hyp.items()}
Z = sum(unnorm.values())
for k in hyp:
    print(f"{k}   prior {prior[k]:.2f}"
          f"   P(n=5 | this) {poisson.pmf(n, hyp[k]*E):.5f}"
          f"   -> posterior {unnorm[k]/Z:.3f}")
print(f"\nnormalizer Z = P(n=5) = {Z:.5f}")
```

```
flat  (lambda=1)   prior 0.50   P(n=5 | this) 0.03609   -> posterior 0.183
rising(lambda=3)   prior 0.50   P(n=5 | this) 0.16062   -> posterior 0.817

normalizer Z = P(n=5) = 0.09836
```

Start at 50/50. Observing 5 deaths against an expectation of 2 moves you to
82% on "rising." The mechanics are just: multiply each prior by how well that
hypothesis predicted the data, then rescale so the column sums to 1.

**The one conceptual jump for what follows:** λ isn't really a choice between
two values — it's any positive number. So instead of two probabilities we need a
probability *density* over the whole positive line: a curve saying which λ
values are plausible. The posterior stops being a number and becomes a
distribution. `prior × likelihood`, normalized, still does the work.

---

## 6. The Gamma distribution

We need a distribution over λ. Requirements: λ can't be negative (a rate of
−2 is meaningless), and it should be flexible enough to be either vague or
sharp. The **Gamma distribution** is the standard choice, with two parameters:

```
Gamma(α, β)      density ∝ λ^(α−1) · e^(−βλ)

mean = α / β         variance = α / β²        sd = √α / β
```

- **α** ("shape") — controls the *form* of the curve. Below 1 it piles up at
  zero; at 1 it's a plain exponential decay; above 1 it develops a hump.
- **β** ("rate") — a scale knob. Doubling β halves everything.

The pairing that matters here: **α controls how tight the distribution is,
independent of where it's centred.** Fix the mean at α/β = 1 and raise both
together, and the curve concentrates:

```python
# --- 6 ---
import numpy as np
from scipy.stats import gamma as G

print("alpha  beta    mean     sd    5th   50th   95th   shape of the curve")
for a, b in [(0.5,0.5), (1,1), (3,3), (7.3,7.5), (21.5,21.5), (3,1)]:
    d = G(a, scale=1/b)
    shape = "spike at 0" if a < 1 else "exponential decay" if a == 1 else "humped"
    print(f"{a:5.1f} {b:5.1f} {a/b:7.2f} {np.sqrt(a)/b:6.2f}"
          f" {d.ppf(.05):6.2f} {d.ppf(.5):6.2f} {d.ppf(.95):6.2f}   {shape}")
```

```
alpha  beta    mean     sd    5th   50th   95th   shape of the curve
  0.5   0.5    1.00   1.41   0.00   0.45   3.84   spike at 0
  1.0   1.0    1.00   1.00   0.05   0.69   3.00   exponential decay
  3.0   3.0    1.00   0.58   0.27   0.89   2.10   humped
  7.3   7.5    0.97   0.36   0.47   0.93   1.63   humped
 21.5  21.5    1.00   0.22   0.67   0.98   1.38   humped
  3.0   1.0    3.00   1.73   0.82   2.67   6.30   humped
```

Every row but the last is centred at 1, but the 5th-to-95th span narrows from
`0.00–3.84` to `0.67–1.38`. The fourth row is roughly the prior `trends.py`
actually fits on this data (`Gamma(7.267, 7.547)`), and it says: *typical
substances sit near 1, and 90% of them fall between 0.47 and 1.63.*

`scipy` note: scipy parameterizes Gamma by *scale*, which is `1/β`. Hence
`scale=1.0/pr["beta"]` throughout `trends.py`. Mixing these up is the single
most common bug in this area.

**Interpreting α and β, which pays off next:** think of `Gamma(α, β)` as the
belief you'd hold after having already seen **α events over β units of
exposure**. `Gamma(7.3, 7.5)` is "I've effectively watched 7.3 events across 7.5
expected — I'm fairly confident the rate is near 1, but not certain." That
reading is about to become literal.

---

## 7. Conjugacy: the trick that makes it fast

Combine what we have. Prior on λ is Gamma(α, β); the data likelihood is Poisson.
Drop every factor that doesn't involve λ (the `∝` absorbs them):

```
likelihood   P(n | λ)  ∝  λⁿ · e^(−λE)          ← Poisson, from step 2
prior        P(λ)      ∝  λ^(α−1) · e^(−βλ)     ← Gamma,   from step 6

multiply:              ∝  λ^(α+n−1) · e^(−(β+E)λ)
```

Look at the result and compare it to the Gamma template `λ^(α−1)·e^(−βλ)`. It is
*the same form*, with `α → α+n` and `β → β+E`. So:

```
posterior  =  Gamma(α + n, β + E)
```

This property — prior and posterior in the same family — is called
**conjugacy**, and the Gamma is said to be *conjugate* to the Poisson. It is why
`_gps_scores` is six lines with no integration, no MCMC, no numerical anything,
and why a 45-quarter backtest sweep costs about a second.

Don't take my word for the algebra. Compute the posterior by brute force on a
grid and compare:

```python
# --- 7 ---
import numpy as np
from scipy.stats import gamma as G, poisson

a, b = 3.0, 2.0        # prior Gamma(alpha=3, beta=2)
n, E = 5, 2.0          # observed 5, expected 2

grid = np.linspace(1e-6, 20, 200_000)
post = G.pdf(grid, a, scale=1/b) * poisson.pmf(n, grid * E)   # prior x likelihood
post /= np.trapezoid(post, grid)                              # normalize

closed = G.pdf(grid, a + n, scale=1/(b + E))                  # the claim

print(f"brute-force posterior mean   {np.trapezoid(grid*post, grid):.6f}")
print(f"Gamma({a+n:.0f}, {b+E:.0f}) mean            {(a+n)/(b+E):.6f}")
print(f"largest density disagreement {np.abs(post-closed).max():.2e}")
```

```
brute-force posterior mean   2.000000
Gamma(8, 4) mean            2.000000
largest density disagreement 1.55e-15
```

Agreement to machine precision. The closed form isn't an approximation — it's
the same object, written down instead of computed.

**And now the α/β interpretation pays off.** `α + n`, `β + E`: you *literally add
your observed events to α and your observed exposure to β.* The prior was worth
"α pseudo-events over β pseudo-exposure," and the data adds n real events over E
real exposure. The posterior mean is

```
(α + n) / (β + E)
```

which is a weighted blend of the prior mean α/β and the raw estimate n/E, with
the weights set by how much exposure each side brings. Large E → the data
dominates. Small E → the prior dominates.

That is shrinkage. Nobody designed it; it fell out of the algebra.

---

## 8. Shrinkage, seen directly

Here is the payoff. Six substances, run through `Gamma(7.267, 7.547)` — the
prior `trends.py` fits on the real data. The first four have an **identical raw
ratio of 5.0** and differ only in how much evidence backs it:

```python
# --- 8 ---
import numpy as np
from scipy.stats import gamma as G
from scipy.special import digamma

a, b = 7.267, 7.547            # the prior trends.py actually fits

print("     n       E   raw n/E    posterior          EBGM   EB05   EB95")
for n, E in [(2, 0.4), (20, 4), (200, 40), (2000, 400), (0, 2.0), (137, 90.57)]:
    pa, pb = a + n, b + E
    d = G(pa, scale=1/pb)
    print(f"{n:6d} {E:7.1f} {n/E:9.2f}    Gamma({pa:7.1f},{pb:6.1f})"
          f" {np.exp(digamma(pa)-np.log(pb)):6.2f} {d.ppf(.05):6.2f} {d.ppf(.95):6.2f}")
```

```
     n       E   raw n/E    posterior          EBGM   EB05   EB95
     2     0.4      5.00    Gamma(    9.3,   7.9)   1.10   0.62   1.86
    20     4.0      5.00    Gamma(   27.3,  11.5)   2.32   1.67   3.15
   200    40.0      5.00    Gamma(  207.3,  47.5)   4.35   3.87   4.87
  2000   400.0      5.00    Gamma( 2007.3, 407.5)   4.92   4.75   5.11
     0     2.0      0.00    Gamma(    7.3,   9.5)   0.71   0.36   1.28
   137    90.6      1.51    Gamma(  144.3,  98.1)   1.47   1.27   1.68
```

Read the EBGM column top to bottom. Same raw ratio of 5.00 every time, and the
estimate climbs 1.10 → 2.32 → 4.35 → 4.92. With n=2 the prior swamps the data
and the estimate is dragged almost all the way back to 1. With n=2000 the prior
is a rounding error and the estimate sits essentially at the raw 5.0.

**The amount of shrinkage is automatic** — set by the ratio of real exposure E
to prior pseudo-exposure β. There is no threshold to tune, no minimum-count
rule. Small substances get pulled back hard because β=7.5 is large relative to
their E; large substances aren't, because it isn't.

Two more rows worth noting:

- **n=0, E=2** — a substance that vanished. The posterior is `Gamma(7.3, 9.5)`,
  mean 0.77, EB95 = 1.28. Zero observations is *evidence*, and it's handled by
  the same formula with no special case.
- **n=137, E=90.6** — this is PCP, from the real table. EBGM 1.47, EB05 1.27.
  Compare against the live `trends rank` output: `ebgm 1.47`, `eb05 1.27`. You
  have just reproduced the top row of the production ranking by hand.

---

## 9. Where does the prior come from? The negative binomial

One thing is still unexplained. Step 8 used `Gamma(7.267, 7.547)`, which I
pulled from thin air. Where does a prior come from if you don't have one?

**The empirical Bayes answer: estimate it from the 203 substances themselves.**
Each substance is one draw from a population of substances, and the population
tells you what a typical λ looks like. The prior isn't a belief you bring — it's
a summary of the table you already have.

To fit α and β we need `P(n | E, α, β)` — the probability of a count *without
knowing λ*, since λ is exactly what we can't see. Getting it means averaging the
Poisson likelihood over every possible λ, weighted by the prior:

```
P(n | E) = ∫ P(n | λE) · Gamma(λ; α, β) dλ
```

That integral has a closed form, and it is the **negative binomial**:

```
P(n) = Γ(n+α) / (Γ(α)·n!) · (β/(β+E))^α · (E/(β+E))^n
```

Compare to `_nb_logpmf` at `trends.py:293` — it is this, in logs, term for term
(comments added here):

```python
return (gammaln(n + a) - gammaln(a) - gammaln(n + 1)      # log Γ(n+α)/(Γ(α)n!)
        + a * np.log(b / (b + e))                          # log (β/(β+E))^α
        + n * np.log(e / (b + e)))                         # log (E/(β+E))^n
```

(`Γ` is the gamma *function*, a continuous extension of factorial:
`Γ(k) = (k−1)!` for whole numbers. `gammaln` is its log, used because the raw
values overflow. Logs also turn the product over 203 substances into a sum,
which is what the optimizer wants.)

Confirm the whole chain by simulating it — draw λ from the prior, then a count
from Poisson, and check the resulting counts against the negative binomial:

```python
# --- 9 ---
import numpy as np
from scipy.stats import nbinom
rng = np.random.default_rng(0)

a, b, E = 7.267, 7.547, 3.0

lam = rng.gamma(a, 1/b, 500_000)      # step 1: draw a rate from the prior
sim = rng.poisson(lam * E)            # step 2: draw a count given that rate

ks = range(0, 9)
print("k            ", "  ".join(f"{k:>6}" for k in ks))
print("simulated    ", "  ".join(f"{(sim==k).mean():6.4f}" for k in ks))
print("nbinom pmf   ", "  ".join(f"{nbinom.pmf(k, a, b/(b+E)):6.4f}" for k in ks))
print(f"\nmean {sim.mean():.3f}   formula a*E/b            = {a*E/b:.3f}")
print(f"var  {sim.var():.3f}   formula mean * (1 + E/b) = {a*E/b*(1+E/b):.3f}")
print(f"variance / mean = {sim.var()/sim.mean():.3f}  -> wider than Poisson, which is 1.000")
```

```
k                  0       1       2       3       4       5       6       7       8
simulated     0.0878  0.1811  0.2139  0.1864  0.1364  0.0891  0.0511  0.0278  0.0141
nbinom pmf    0.0878  0.1816  0.2135  0.1876  0.1369  0.0878  0.0510  0.0275  0.0140

mean 2.894   formula a*E/b            = 2.889
var  4.051   formula mean * (1 + E/b) = 4.037
variance / mean = 1.400  -> wider than Poisson, which is 1.000
```

**The key structural fact, and remember it for step 12:** the negative binomial
is *wider than Poisson* — φ = 1.40 here, not 1.00. That extra width is
`(1 + E/β)`, and it exists because different substances genuinely have different
λ. The NB's job is to absorb between-substance spread.

That is a virtue when the spread is real. It becomes the failure mechanism in
step 12, when it absorbs spread that isn't.

---

## 10. Empirical Bayes: fitting the prior from the data

Now `_fit_gps_prior` (`trends.py:309`) is readable. It searches for the (α, β)
that make the observed 203 (n, E) pairs as likely as possible under the negative
binomial — plain maximum likelihood, by Nelder–Mead, on the summed log
likelihood.

Does that actually work? Test it the only honest way: generate data from a
*known* prior, hide the prior, and see if the fit finds it.

```python
# --- 10 ---
import numpy as np
from scipy.special import gammaln
from scipy.optimize import minimize
rng = np.random.default_rng(0)

def nb_logpmf(n, e, a, b):                       # == trends.py:293
    e = np.maximum(e, 1e-9)
    return (gammaln(n + a) - gammaln(a) - gammaln(n + 1)
            + a * np.log(b / (b + e)) + n * np.log(e / (b + e)))

TRUE_A, TRUE_B = 7.267, 7.547
print("substances   recovered alpha   recovered beta")
for m in [50, 113, 1000, 20000]:
    E   = np.exp(rng.normal(1.5, 1.6, m)).clip(0.05, 4000)
    lam = rng.gamma(TRUE_A, 1/TRUE_B, m)         # nature's hidden rates
    n   = rng.poisson(lam * E).astype(float)     # all we ever observe is (n, E)
    r = minimize(lambda p: -nb_logpmf(n, E, *np.exp(p)).sum(), [0., 0.],
                 method="Nelder-Mead",
                 options={"maxiter":10000, "maxfev":10000, "xatol":1e-8, "fatol":1e-8})
    a, b = np.exp(r.x)
    print(f"{m:10d} {a:17.2f} {b:16.2f}")
print(f"{'truth':>10} {TRUE_A:17.2f} {TRUE_B:16.2f}")
```

```
substances   recovered alpha   recovered beta
        50              8.56             9.13
       113              8.57             8.91
      1000              6.89             7.11
     20000              7.36             7.68
     truth              7.27             7.55
```

It works, and the row to study is **113** — the number of substances with a
non-zero baseline in the real fit. It recovers 8.57 against a truth of 7.27:
right neighbourhood, roughly 18% off. That is honest resolution at this sample
size, and it is precisely the argument in `_fit_gps_prior`'s docstring for
fitting *one* Gamma rather than DuMouchel's two-component mixture. Two
components need five parameters. If 113 substances leave this much wobble in
two, five is not identifiable — which is what the docstring reports happening.

Notes on the real implementation:

- **Substances with E = 0 are excluded from the fit** but still scored. They
  carry no information about the prior (no baseline exposure), but the posterior
  `Gamma(α+n, β+0)` is perfectly well defined for them. `_nb_logpmf` floors E at
  `1e-9` so this doesn't divide by zero, and the docstring's reading is right:
  n=0 costs nothing, n>0 becomes very surprising, which is what "appeared from
  nothing" should mean.
- **Optimizing in log space** (`np.exp(params)` inside the objective) keeps α
  and β positive without a constrained optimizer.
- Estimating the prior from the same data you then apply it to feels circular.
  It is, mildly — that's the "empirical" in empirical Bayes, and it's a
  well-understood, widely accepted trade. With 113 substances informing two
  parameters, no single substance meaningfully sets its own prior.

---

## 11. EBGM and EB05: what the ranking promises

Every substance now has a posterior `Gamma(α+n, β+E)`. To rank, you need one
number out of each distribution — and which number you pick is a real choice.

**EBGM** — Empirical Bayes Geometric Mean, `exp(E[log λ])`:

```python
"ebgm": np.exp(digamma(post_a) - np.log(post_b))
```

The posterior at these counts is right-skewed — a long tail toward large λ. An
arithmetic mean gets dragged out by that tail; the geometric mean doesn't. It's
the stable centre. (`digamma` is the derivative of `log Γ`, and gives
`E[log λ]` exactly for a Gamma.)

**EB05** — the **5th percentile** of the posterior:

```python
"eb05": gamma_dist.ppf(0.05, post_a, scale=1.0 / post_b)
```

Not a centre at all — a *conservative lower bound*. It reads:

> Given everything observed, there is a 95% posterior probability that this
> substance's true λ is **at least** this large.

So `EB05 > 1` is a specific, falsifiable claim: 95% posterior confidence the
substance is above baseline. Ranking by it means ranking by *how much you can
defend*, not by best guess. A small substance can post a high EBGM but never a
high EB05, because its posterior is too wide for the 5th percentile to climb.

That is the ranking's whole defence against step 4's failure. Demonstrate it —
build a synthetic table where the truth is known, and see which statistic finds
it:

```python
# --- 11 ---
import numpy as np
from scipy.stats import gamma as G
from scipy.special import gammaln
from scipy.optimize import minimize
rng = np.random.default_rng(0)

def nb_logpmf(n, e, a, b):
    e = np.maximum(e, 1e-9)
    return (gammaln(n+a) - gammaln(a) - gammaln(n+1)
            + a*np.log(b/(b+e)) + n*np.log(e/(b+e)))
def fit(n, e):
    r = minimize(lambda p: -nb_logpmf(n, e, *np.exp(p)).sum(), [0., 0.],
                 method="Nelder-Mead", options={"maxiter":10000, "maxfev":10000})
    return np.exp(r.x)

m   = 203
E   = np.exp(rng.normal(1.5, 1.6, m)).clip(0.05, 4000)
lam = rng.gamma(7.267, 1/7.547, m)          # the TRUTH -- never visible in real life
n   = rng.poisson(lam * E).astype(float)

a, b = fit(n, E)
raw  = n / np.maximum(E, 1e-9)
eb05 = G.ppf(0.05, a + n, scale=1/(b + E))
print(f"fitted prior Gamma({a:.2f}, {b:.2f})   [generated from Gamma(7.27, 7.55)]\n")

print("  ranked by RAW n/E                |   ranked by EB05")
print("     n      E  ratio  TRUE lam     |      n       E  EB05  TRUE lam")
for x, y in zip(np.argsort(-raw)[:8], np.argsort(-eb05)[:8]):
    print(f"  {n[x]:5.0f} {E[x]:6.2f} {raw[x]:6.2f} {lam[x]:9.2f}     |"
          f"  {n[y]:5.0f} {E[y]:7.1f} {eb05[y]:5.2f} {lam[y]:9.2f}")

for name, s in [("raw n/E", raw), ("EB05", eb05)]:
    t = np.argsort(-s)[:10]
    print(f"\ntop 10 by {name:8s}  mean TRUE lambda {lam[t].mean():.2f}"
          f"   truly above 1.5: {(lam[t] > 1.5).sum()}/10"
          f"   median E {np.median(E[t]):6.1f}")
```

```
fitted prior Gamma(7.15, 7.72)   [generated from Gamma(7.27, 7.55)]

  ranked by RAW n/E                |   ranked by EB05
     n      E  ratio  TRUE lam     |      n       E  EB05  TRUE lam
      1   0.21   4.68      0.40     |    208   110.4  1.62      1.81
      2   0.78   2.56      1.23     |    149    82.7  1.51      1.82
      9   3.65   2.47      0.76     |     58    32.0  1.32      1.68
      1   0.42   2.41      1.06     |     51    28.9  1.26      1.81
      7   3.07   2.28      1.52     |    113    80.0  1.17      1.24
      3   1.33   2.25      1.35     |     61    41.4  1.12      1.44
      3   1.38   2.17      1.05     |     36    21.8  1.11      1.38
      3   1.39   2.16      1.59     |     23    12.3  1.09      1.48

top 10 by raw n/E   mean TRUE lambda 1.16   truly above 1.5: 2/10   median E    1.4

top 10 by EB05      mean TRUE lambda 1.49   truly above 1.5: 4/10   median E   38.5
```

The `TRUE lam` columns are the verdict, and they are unambiguous.

The raw ratio's **number one substance** — ratio 4.68, an apparent quadrupling —
has a true λ of **0.40**. It is *declining by 60%*, and it tops the table on the
strength of a single death. Its top 10 averages a true λ of 1.16, barely above
baseline, with a median E of 1.4: it has essentially ranked the substances by
smallness, as step 4 predicted.

EB05's top 10 averages a true λ of 1.49, catches twice as many genuinely
elevated substances, and its picks have a median E of 38.5. Its number one has a
true λ of 1.81 — the real thing.

Note also what EB05 is *not* doing: it isn't just "rank by n." Row 5 (n=113,
EB05 1.17) sits below row 4 (n=51, EB05 1.26), because 51-against-28.9 is a
proportionally bigger excess than 113-against-80. It trades off effect size
against evidence, which is exactly the job.

---

## 12. Why weighted counts break it, and the φ fix

Now the question that started this. `GPS_V2_DESIGN.md` §1 proposes replacing
"every mention counts 1" with position-based credit from the cause line
(`trends.py:157`):

```python
SOLO_CREDIT     = 2   # unambiguous: nothing else named
NAMED_CREDIT    = 1   # named, and not (last on a line of 3+)
TRAILING_CREDIT = 0   # last on a line of 3+ -- the adulterant pattern
```

The count becomes a sum of credits: `W = c₁ + c₂ + … + c_N`, where N is Poisson
and each `cᵢ ∈ {0,1,2}`.

**Is W still Poisson?** Step 3 said the entire closed form rests on
`variance = mean`. So check it. A sum like this is a **compound Poisson**, and
its moments are:

```
E[W]   = μ · E[c]
Var(W) = μ · E[c²]

φ = Var(W)/E[W] = E[c²] / E[c]
```

```python
# --- 12a ---
import numpy as np
rng = np.random.default_rng(0)

mu = 500.0
print("credit scheme        phi = E[c^2]/E[c]   mean(W)   var(W)   var/mean")
for label, draw in [
    ("c = 1 (raw count)",  lambda k: np.ones(k)),
    ("c in {0,1}",         lambda k: rng.choice([0, 1], k)),
    ("c in {0,1,2}",       lambda k: rng.choice([0, 1, 2], k)),
    ("c = 2 always",       lambda k: np.full(k, 2.0)),
    ("c ~ Uniform[0,1]",   lambda k: rng.uniform(0, 1, k)),
]:
    c = draw(400_000)
    phi = (c**2).mean() / c.mean()
    W = np.array([draw(k).sum() for k in rng.poisson(mu, 40_000)])
    print(f"{label:20s} {phi:16.2f} {W.mean():9.1f} {W.var():8.1f} {W.var()/W.mean():10.2f}")
```

```
credit scheme        phi = E[c^2]/E[c]   mean(W)   var(W)   var/mean
c = 1 (raw count)                1.00     499.9    504.9       1.01
c in {0,1}                       1.00     250.0    251.8       1.01
c in {0,1,2}                     1.67     500.2    826.0       1.65
c = 2 always                     2.00     999.8   1992.0       1.99
c ~ Uniform[0,1]                 0.67     249.9    167.0       0.67
```

The predicted φ matches the simulated var/mean in every row. Read off the
consequences:

- **`c ∈ {0,1}` stays exactly Poisson.** Deleting points from a Poisson process
  leaves a Poisson process — the *thinning theorem*. Discarding mentions is free.
- **`c = 2 always` gives φ = 2.** And this is the intuition to keep: *doubling a
  Poisson is not a Poisson.* If N ~ Poisson(500), then 2N has mean 1000 but
  variance 2000, and can never take an odd value. Poisson(1000) has variance
  1000. Scaling a count doubles the mean and quadruples the variance.
- **`c ∈ {0,1,2}` lands at φ = 1.67**, between the two.
- **Continuous weights in [0,1] give φ = 0.67** — *under*dispersed, and
  therefore *safer* on this axis than the integer scheme.

That last row matters, because the design note's original rationale was
"integers safe, floats dangerous." **That framing is wrong.** Integrality is
irrelevant; what governs φ is the *scale* of the weights. `{0,1}` (integers) is
perfectly safe and `{0,1,2}` (also integers) is not.

### What actually goes wrong

The failure is indirect and worth tracing, because it isn't "the numbers get
noisier." Recall step 9: the negative binomial's job is to absorb spread. It
cannot tell *within*-substance noise from *between*-substance heterogeneity. So
extra dispersion from weighting gets misread as "substances genuinely differ a
lot," and the fitted prior widens — α falls. A small α means a vague prior,
which means **weak shrinkage** (step 8), which means the small-count noise EB05
exists to suppress comes straight back.

Measure it. Truth: λ = 1 for all 203 substances, so *every* `EB05 > 1` is a
false alarm.

```python
# --- 12b ---
import numpy as np
from scipy.stats import gamma as G
from scipy.special import gammaln
from scipy.optimize import minimize
rng = np.random.default_rng(0)

def nb_logpmf(n, e, a, b):
    e = np.maximum(e, 1e-9)
    return (gammaln(n+a) - gammaln(a) - gammaln(n+1)
            + a*np.log(b/(b+e)) + n*np.log(e/(b+e)))
def fit(n, e):
    r = minimize(lambda p: -nb_logpmf(n, e, *np.exp(p)).sum(), [0., 0.],
                 method="Nelder-Mead", options={"maxiter":10000, "maxfev":10000})
    return np.exp(r.x)

E = np.exp(rng.normal(1.5, 1.6, 203)).clip(0.05, 4000)

def sweep(label, draw, phi=1.0, reps=300):
    """Truth: lambda = 1 for all 203 substances. Every EB05 > 1 is a false alarm."""
    fa, al = [], []
    for _ in range(reps):
        W  = np.array([draw(k).sum() for k in rng.poisson(E)], float) / phi
        Wb = np.array([draw(k).sum() for k in rng.poisson(E)], float) / phi
        a, b = fit(W, Wb)
        fa.append((G.ppf(0.05, a + W, scale=1/(b + Wb)) > 1).sum())
        al.append(a)
    print(f"{label:42s} alpha {np.median(al):6.2f}   false EB05>1 {np.mean(fa):5.1f}/203")

d012 = lambda k: rng.choice([0, 1, 2], k)
sweep("raw counts, c = 1            (phi 1.00)", lambda k: np.ones(k))
sweep("thinned, c in {0,1}          (phi 1.00)", lambda k: rng.choice([0,1], k))
sweep("credits c in {0,1,2}         (phi 1.67)", d012)
sweep("credits {0,1,2}, /= phi      (phi 1.67)", d012, phi=5/3)
```

```
raw counts, c = 1            (phi 1.00)    alpha  15.32   false EB05>1   6.8/203
thinned, c in {0,1}          (phi 1.00)    alpha   7.17   false EB05>1   8.0/203
credits c in {0,1,2}         (phi 1.67)    alpha   2.23   false EB05>1  29.5/203
credits {0,1,2}, /= phi      (phi 1.67)    alpha   8.31   false EB05>1   7.3/203
```

Uncorrected weighting collapses α from 15.3 to 2.2 and **quadruples the false
alarm rate**, 6.8 → 29.5. The mechanism is exactly as traced: dispersion →
diffuse prior → weak shrinkage → noise returns.

### The correction

If W has dispersion φ, then `W/φ` has mean `μ/φ` and variance `μ/φ` — back to
equidispersed. Divide **both** n and E by φ:

- `λ = n/E` is unchanged, so the estimate isn't biased;
- the *effective sample size* drops by a factor of φ, so shrinkage strengthens.

That is the right trade and the right direction: overdispersed data carries less
information, and this is how you tell the model so. The last row confirms it —
α back to 8.31, false alarms back to 7.3, level with raw counts.

Two things make this defensible rather than a fudge:

1. **φ is measured, not tuned.** It's `E[c²]/E[c]` computed directly off the
   observed credit distribution. This is `_credit_dispersion` at `trends.py:191`,
   and on the real cause lines it comes out at **φ = 1.422** — reported in the
   `trends rank --weighted` header.
2. **Feeding fractional n through `_nb_logpmf` is fine.** `gammaln` is defined on
   the reals, so nothing breaks mechanically. What changes is interpretive: the
   output is a *quasi-likelihood* rather than a normalized pmf over the integers.
   That is the standard, accepted treatment for overdispersed count data.

Note the irony worth keeping: the correct fix **requires** the floats the
original note ruled out.

### What it looks like on the real data

```
$ emerging trends rank
GPS prior: Gamma(alpha=7.267, beta=7.547), mean 0.963, sd 0.357, fitted on 113 substances

$ emerging trends rank --weighted
GPS prior: Gamma(alpha=8.554, beta=8.923), mean 0.959, sd 0.328, fitted on 106 substances
dispersion correction: phi=1.422 (1.0 = no overdispersion from the credit weighting)
```

With the correction applied, the weighted fit lands at α = 8.55 — the same
neighbourhood as the unweighted 7.27, not collapsed to 2. That is what a
correction working looks like.

### Caveat: φ is not one number

`_credit_dispersion` computes a single pooled φ, but the per-substance values
vary a lot:

```
Methamphetamine  1.61     Fentanyl  1.37     Alprazolam       1.05
Ethanol          1.57     Heroin    1.33     PCP              1.00
Cocaine          1.44     Morphine  1.28     Acetyl fentanyl  1.00
```

The established substances are often named solo (credit 2) and run 1.4–1.6; PCP
and acetyl fentanyl sit at exactly 1.00 because their credits are only ever 0 or
1. A single global φ over- and under-corrects at the two ends. Fixing that
properly means a per-substance φ_s, which is straightforward to compute but
changes the shrinkage in ways that need their own backtest. Filed as known and
unaddressed, not solved.

---

## Cheat sheet

| Symbol | Meaning | In the code |
|---|---|---|
| `n` | observed count, recent window | `n_recent` |
| `E` | expected count from baseline share | `expected` |
| `λ` | true rate multiplier; 1 = on baseline | never observed |
| `α, β` | prior Gamma shape and rate | `pr["alpha"]`, `pr["beta"]` |
| `φ` | variance-to-mean ratio; 1 = Poisson | `_credit_dispersion` |
| `Γ` | gamma function, continuous factorial | `gammaln` |

| Fact | Why it matters |
|---|---|
| `n ~ Poisson(λE)` | the model; variance = mean, no free spread parameter |
| `sd(n/E) = 1/√E` | raw ratios are `1/√E` noisy → small substances dominate any raw ranking |
| Gamma × Poisson = Gamma | conjugacy; closed-form posterior, no MCMC |
| posterior = `Gamma(α+n, β+E)` | prior is α pseudo-events over β pseudo-exposure; data just adds on |
| posterior mean = `(α+n)/(β+E)` | blend of prior mean and `n/E`, weighted by exposure — this *is* shrinkage |
| marginal of n is negative binomial | lets α, β be fit by MLE from all substances at once |
| NB absorbs overdispersion | a virtue for real heterogeneity; the failure route in step 12 |
| EB05 = 5th posterior percentile | "95% confident λ is at least this" — a defensible floor, not a best guess |
| `φ = E[c²]/E[c]` | weights break `variance = mean`; divide n and E by φ to restore it |

**Where each piece lives**

| Concept | Location |
|---|---|
| negative binomial log-likelihood | `trends.py:293` `_nb_logpmf` |
| MLE fit of the prior | `trends.py:309` `_fit_gps_prior` |
| EBGM / EB05 / EB95 | `trends.py:358` `_gps_scores` |
| n, E, and the φ division | `trends.py:376` `rank_substances` |
| credit scheme and φ | `trends.py:157`–`191` |
| the design argument | `docs/GPS_V2_DESIGN.md` §1 |

**Further reading.** DuMouchel (1999), *Bayesian Data Mining in Large Frequency
Tables*, is the original GPS paper and the source of the EBGM/EB05 conventions
used here. Its two-component mixture is the model `_fit_gps_prior`'s docstring
explains this codebase deliberately does not use.
