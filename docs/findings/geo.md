# Spatial localization

Findings for `emerging/geo.py`. Results in [`results/geo/`](../../results/geo/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

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
