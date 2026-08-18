# Spatial variation in temporal trends

Findings for `emerging/svtt.py` (tests: `tests/test_svtt.py`). Results in [`results/svtt/`](../../results/svtt/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## Where is mortality rising fastest?

`svtt.py` is the one analysis here that uses a population denominator, and the
only one that tests a **slope** rather than a level. Each circle of zipcodes
keeps its own baseline rate and is tested only on its trend, so an area that
has always been bad cannot signal — the objection that made `geo.py`
case-control, handled inside the model instead of by avoiding the denominator.

The MLE of a log-linear Poisson with a population offset reduces to a
one-dimensional monotone root find (population-weighted mean time under trend
`b` must equal case-weighted mean time; the derivative is a variance, so the
root is unique). That is what makes 23,406 circles affordable without fitting a
GLM per circle.

### The window matters more than the method

The first run, over 2015–2025, returned a county-wide trend of **+14.5%/yr**.
That is an artifact. LA's overdose rate is a hump, not a trend:

| 2013 | 2016 | 2019 | 2022 | 2024 | 2025 |
|---|---|---|---|---|---|
| 5.3 | 7.6 | 12.7 | **29.5** | 26.2 | 21.7 |

A log-linear fit across all of that describes neither the rise nor the fall,
and "rising fastest" then conflates *rose earliest during the ramp* with *still
rising now*. The default window is therefore the **decline phase** (2022–2025,
county −9.3%/yr, monotone), and `scan` prints the annual rate series it fitted
so the misspecification cannot recur silently. `--first-year 2015` gives the
whole-epidemic version if you want it.

### The Pomona Valley is going the other way

| Year | Cluster per 100k | Rest of county |
|---|---|---|
| 2022 | 11.5 | 30.0 |
| 2023 | 12.5 | 28.8 |
| 2024 | 23.0 | 26.3 |
| 2025 | 20.4 | 21.7 |

Five zips — 91765 Diamond Bar, 91766/91768 Pomona, 91748 Rowland Heights,
91789 Walnut — **rising 25.2%/yr while the county falls 9.8%/yr**, p = 0.016
adjusted for every circle searched.

**It is a catch-up, not a new epicentre.** The cluster started at a third of
the county rate and converged upward; it is still slightly *below* county
average in level today. A level-based scan, including a population-denominated
one, would miss it entirely. That is the case SVTT exists for.

Three caveats to carry with the number:

- **The rise is one year.** 2023→2024 is the whole step and 2025 falls back. A
  single slope summarises that shape rather than describing it — quote the
  table, not just the +25%.
- **The boundary is not stable.** Dropping 2025 moves the top cluster to a
  broad 58-zip eastern-county region (+7.9% vs −7.9%). Every Pomona zip is in
  both, so the *region* is robust and the *extent* is not. Say "eastern San
  Gabriel / Pomona Valley", not five zipcodes.
- **156 deaths.** Enough for a slope, not enough to subdivide.

### The cluster, profiled — and it is a methamphetamine story

`svtt profile` characterises the most likely cluster against the rest of the
county, because "a place is rising" is not reportable until you know what is in
it. For the Pomona cluster:

- **It is Pomona.** 91768 (+19) and 91766 (+16) supply 35 of the 44-death
  increase; the other three zips contribute 9 between them and are along for
  the ride via the circle geometry.
- **A level shift, not a burst** — 2022–23 runs 2–9 deaths a quarter, 2024–25
  runs 6–17, stepping at 2024Q1 with no reversion.
- **Not a pronouncement artifact** — 3.8% at shared facility addresses against
  18.6% county-wide, at 26/27/44/47 distinct addresses across the four years.
- **No substance drives it** — heroin 1.46×, PCP 1.44×, meth 1.16×, fentanyl
  0.89×, cocaine 0.74× against county shares. More of the same supply.
- **Who** — 55% LATINE against 38% county-wide, median age 48 against 44.

Running SVTT per substance separates the two epidemics:

| | county trend | most likely cluster | p |
|---|---|---|---|
| Fentanyl | **−16.25%/yr** | 91765, 7 zips, +3.8%/yr | **0.63** |
| Methamphetamine | **−7.56%/yr** | 91741, 42 zips (Glendora–Azusa–Claremont–Pomona), **+5.98%/yr** | **0.02** |

**Fentanyl is falling uniformly** — no significant spatial variation in its
trend anywhere in the county. **Methamphetamine is not**: it falls half as fast
and is *rising* across a broad eastern San Gabriel Valley region containing the
all-cause cluster. The county-wide decline is fentanyl-led; the eastern-county
rise is stimulant-involved.

That changes what a response would be. A naloxone-and-test-strip package targets
the part of the epidemic already receding fastest; the part that is growing has
no equivalent pharmacological answer.

### What the denominator does and does not fix

`census.py` supplies ACS resident population, 301 LA zips × 2013–2024, and
every one of the 280 zips holding deaths is covered. But residential population
is the wrong denominator for a place people travel to: 90013 has 14,891
residents and a much larger daytime population, so its per-resident rate is
place-based intensity, not individual risk. A *trend* is more robust to this
than a level, since the bias has to be changing to matter — and the Pomona
cluster is suburban, where that assumption is at its strongest. Downtown is
where it is weakest.

`tests/test_svtt.py` checks the trend MLE against the score equation it solves,
that the LLR is exactly zero when two areas share a slope but differ 10× in
level, and that the null is calibrated under a deliberate 5× spread in baseline
level — the confound a level-based scan would fail.
