# Literature review: rare-event detection and spatial clustering

Reviewed 2026-08-10. Scope: methods for (1) detecting a rare, emerging
substance from counts, and (2) deciding whether its deaths are spatially
localized. The question driving the review was narrow — **is there something
better than what this repository already implements?**

Short answer: the empirical-Bayes core is not superseded, but there is one
clearly better method on each side, and both are free, mature, and consume
data we already have.

---

## 1. Rare-event / signal detection

### 1.1 What we implement, and where it sits

`trends.py` uses a single-component **gamma-Poisson shrinker (GPS)** —
DuMouchel's MGPS with one mixture component instead of two, because the
two-component version is not identifiable at 205 substances (see
`_fit_gps_prior`). Ranking is by **EB05**, the 5th posterior percentile.

This is squarely inside the standard pharmacovigilance disproportionality
family: GPS/EBGM (DuMouchel), BCPNN (Bate et al.), PRR, ROR. Nothing about
using it on medical-examiner narratives rather than spontaneous reports
changes its statistical behaviour — the numerator is a substance, the
denominator is all events, and the shrinkage does the same job.

### 1.2 The field's current verdict on disproportionality: "flawed but valuable"

The 2024–26 literature is sharply critical of disproportionality analysis, but
**the criticism is of DA used alone**, not of the estimator:

- [*Defending disproportionality analyses: flawed but valuable*](https://www.sciencedirect.com/science/article/pii/S0022202X26001296)
- [*Charting and Sidestepping the Pitfalls of Disproportionality Analysis*, Drug Safety 2025](https://link.springer.com/article/10.1007/s40264-025-01604-y)
- [READUS-PV reporting guideline and reflections on it](https://pmc.ncbi.nlm.nih.gov/articles/PMC12014456/)
- [Conducting and interpreting disproportionality analyses, Frontiers 2023](https://www.frontiersin.org/journals/drug-safety-and-regulation/articles/10.3389/fdsfr.2023.1323057/full)

The recurring findings: DA is "insufficient to confirm safety signals and
rarely supports regulatory decisions on their own"; ROR in particular is
"highly sensitive to rare events and small sample sizes, where even a few
additional reports can disproportionately inflate values."

**Implication for us: mostly reassuring.** We already triangulate — as-of
backtest, polysubstance verdicts, geographic test, recording diagnostics. The
`count_ratio` column and the lidocaine adulterant finding are precisely the
kind of corroboration this literature demands. Worth adopting **READUS-PV**
reporting conventions in the manuscript.

### 1.3 The upgrade: tree-based / tree-temporal scan statistics (TreeScan)

This is the strongest finding of the review. TreeScan scans a **hierarchy** of
outcomes, evaluating individual nodes *and groups of related nodes*
simultaneously, with Monte Carlo inference that accounts for the multiplicity
of overlapping nodes and time windows.

- [Tree-based scan statistic vs GPS and BCPNN, Br J Clin Pharmacol 2023](https://bpspubs.onlinelibrary.wiley.com/doi/10.1111/bcp.15810) — TreeScan outperformed both for statin adverse events
- [Scoping review of 44 TreeScan studies, Front Pharmacol 2026](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2026.1905426/full)
- [TreeScan for Drugs, FDA Sentinel Initiative](https://www.sentinelinitiative.org/methods-data-tools/methods/treescan-drugs)
- [Zero-inflated variant (TreeScan-ZIP), Sci Rep 2022](https://www.nature.com/articles/s41598-022-19998-5)
- [Earlier GPS vs TreeScan comparison in health plan data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3834945/)

**We already have the hierarchy.** Every substance in `substance_mentions`
carries an NFLIS category — 23 of them, fully populated:
`Fentanyl and Fentanyl-related` (9 substances), `Benzodiazepines` (19),
`Narcotic Analgesics` (18), `Synthetic Cathinones` (3), and so on. That is a
root → category → substance tree, which is TreeScan's native input.

It fixes three distinct weaknesses in the current implementation:

| Weakness now | What TreeScan does |
|---|---|
| **Aggregation power.** Every analog scored alone; nitazenes are "individually too rare to score" and needed a hand-coded family workaround in `watch` | Scans nodes *and* groups. Three nitazenes at n=1 become a testable node at n=3, and the whole fentanyl-analog branch is tested as a unit |
| **Multiple testing.** 205 substances × 45 quarters, entirely unadjusted | Monte Carlo controls for multiplicity across all nodes and windows |
| **Repeated looks.** The as-of sweep re-tests every quarter against a hand-picked EB05 > 1.5 line | **Recurrence interval** — "the duration of surveillance required for the expected number of clusters at least as unusual as observed to equal 1 by chance" — a principled, interpretable alarm calibration |

**The closest published analogue is nearly exact.** NYC DOHMH ran a
[prospective tree-temporal scan on emerging SARS-CoV-2 variants and *Salmonella* clusters](https://pmc.ncbi.nlm.nih.gov/articles/PMC11984460/):
hierarchical nomenclature, rare emerging entities, prospective weekly looks.
It flagged EG.5.1 at a recurrence interval of 35 years, **seven weeks before
WHO designated it a variant of interest**. Swap lineages for drug analogs and
that is this project.

Caveats, taken from the sources rather than assumed:

- The scoping review states TreeScan is for "exploratory detection, not
  real-time surveillance" and yields "hypothesis signals rather than causal
  estimates" — the same epistemic status as EB05, so nothing is lost.
- The NYC paper's **nomenclature-lag** limitation transfers directly: NFLIS
  categories update on DEA's schedule, so a genuinely novel analog may arrive
  uncategorised and sit at the root until the catalog catches up.
- It is a standalone C++ binary driven by parameter files, not a Python
  library. It breaks the pure-Python pipeline, though it is scriptable.
- Bernoulli and Poisson models want a comparator or an expected-count model;
  the tree-temporal self-controlled variant is the one that fits us best.

### 1.4 Considered and rejected: MaxSPRT

[MaxSPRT / CMaxSPRT](https://ncbi.nlm.nih.gov/pmc/articles/PMC6955154) with
alpha spending is the canonical fix for sequential repeated looks, and there
is a [2026 scoping review of sequential methods in post-marketing surveillance](https://journals.sagepub.com/doi/10.1177/20420986261454206).

Rejected for now: MaxSPRT is designed for a pre-specified hypothesis with a
defined surveillance start date. Applying it to 205 substances scanned
simultaneously is awkward, and TreeScan's recurrence interval buys the same
protection in a form that fits the problem. Revisit only if we commit to
prospective monitoring of a short fixed watchlist.

### 1.5 A gap we under-weighted: right-truncation

Our handling is the crudest available — a hard `CUTOFF` that discards every
quarter whose toxicology is unsettled, which throws away the most
decision-relevant data.

[CDC NCHS publishes provisional overdose counts](https://www.cdc.gov/nchs/nvss/vsrr/drug-overdose-data.Htm)
using regression-based nowcasting against a fitted reporting-delay
distribution. Fitting a delay distribution would let us score 2026 with honest
uncertainty instead of deleting it.

Cautionary counterpoint, and it applies to us:
[*The 2025 Drug Overdose Spike That Wasn't*](https://pmc.ncbi.nlm.nih.gov/articles/PMC13066679/)
documents the CDC nowcast breaking down when the epidemic changed regime —
which, per `regime`, ours demonstrably is (−18% denominator across the current
windows).

### 1.6 Domain context

- [Global patterns of NPS in fatalities: systematic review and meta-analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12822845/)
  — pooled NPS detection rose from 2.5% (2008–2014) to 9.3% (2015–2025). LA's
  narrative-derived rate should be sanity-checked against this.
- [Legislating novel psychoactive substances: 15 years of UK mortality data](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2025.1708335/full)
  — a natural experiment on whether scheduling changes move mortality; relevant
  to the unanswered "policy" half of the `regime` question.
- [Aegis 2025 mid-year NPS update](https://www.aegislabs.com/clinical-update/2025-mid-year-update-on-nps-trends/)
  — commercial lab trend data; reports tianeptine up 36% Q1→Q2 2025 nationally.
  We have **1 lifetime tianeptine mention**, which is either a real regional
  difference or a panel gap.

---

## 2. Spatial clustering

### 2.1 What we implement

`geo.py` runs a **random-labelling permutation test** on a marked point
pattern: mean nearest-neighbour distance among a substance's deaths, compared
against 999 resamples of the same n drawn from all overdose deaths. Reported
as `localization = null NND / observed NND`.

Its named relatives: **Cuzick–Edwards** k-nearest-neighbour test (1990, same
label-permutation null, counts rather than distances), **Diggle–Chetwynd**
D(s) = K_cases − K_controls (1991), **Clark–Evans** nearest-neighbour index
(1954, from which the ratio form is borrowed), and **Kelsall–Diggle / sparr**
kernel relative-risk surfaces.

### 2.2 Our test is not obsolete

The definitive power study —
[1,220,000 simulated datasets under 51 clustering models](https://ij-healthgeographics.biomedcentral.com/articles/10.1186/1476-072X-2-9)
— compared Cuzick–Edwards k-NN, the spatial scan statistic, Tango's MEET,
Besag–Newell's R, Moran's I and others. Findings: the spatial scan statistic is
best for **localized** clusters, Tango's MEET for **global** clustering, and
"with appropriate choice of parameter, Besag-Newell's R and Cuzick-Edwards'
k-NN also perform well." Our variant sits in a family with demonstrated power.
No reason to remove it.

### 2.3 The upgrade: SaTScan's *spatial variation in temporal trends*

This answers the question the project actually asks, which our current test
does not.

Our test is static: *is this substance concentrated?* The
[spatial-variation-in-temporal-trends model](https://www.satscan.org/cgi-bin/satscan/register.pl/SaTScan_Users_Guide.pdf?todo=process_userguide_download)
instead "looks for geographical areas with unusually high or low temporal
trends" — *where is it growing fastest*. That fuses the temporal detector and
the spatial test into one statistic rather than running them separately and
reading both by eye.

Also directly relevant: the **space-time permutation scan statistic** "uses
only case numbers with no need for population-at-risk data" and "adjusts
nonparametrically for both purely spatial and purely temporal variation."
That is exactly our constraint — no usable denominator, which is why deaths
were used as the control pool in the first place.

- [Prospective spatiotemporal cluster detection with SaTScan: tutorial, JMIR 2024](https://publichealth.jmir.org/2024/1/e50653)
- [SaTScan systematic review, Communications Medicine 2024](https://www.nature.com/articles/s43856-024-00699-1)
- [Space-time interactions in fatal opioid overdoses, Riverside County 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12799247/)
- [Are relative risk estimates from scan statistics biased?](https://pmc.ncbi.nlm.nih.gov/articles/PMC4047196/) — yes, somewhat; read cluster RRs with care
- [Maximum linkage space-time permutation scan statistics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4071024/)

### 2.4 Promote PCP specifically to `sparr`

Kernel relative-risk surfaces were originally rejected here on sample size,
citing n ≥ 50–100. **PCP now has n = 367**, comfortably clear.
[`sparr`](https://tilmandavies.github.io/sparr/reference/sparr-package.html) is
actively maintained (docs updated April 2025), and its adaptive-bandwidth
machinery (Davies & Hazelton 2010; Davies, Jones & Hazelton 2016 symmetric
adaptive smoothing) is the right way to turn "PCP is 1.38× tighter" into an
actual map of *where*. See also the
[tutorial on kernel estimation of spatial relative risk](https://arxiv.org/pdf/1707.06888).

### 2.5 Considered and rejected: LGCP / inlabru

[Bayesian inference for case-control point patterns with inlabru](https://arxiv.org/abs/2503.14954)
(March 2025) is the most statistically complete option — log-Gaussian Cox
models, Matérn spatial random effects, covariates, full posteriors, INLA
inference. Related: [multivariate point patterns in disease mapping](https://arxiv.org/pdf/1903.11647).

Rejected for triage use: it is R, requires SPDE mesh construction, and at
n = 9–39 the spatial random field would be prior-dominated — a beautiful
posterior driven mostly by the prior. The right tool if a single substance
becomes a funded investigation with covariates; the wrong tool for a
comparable column across 42 substances.

---

## 3. Verdict

| Priority | Change | Rationale | Cost |
|---|---|---|---|
| **1** | **TreeScan, prospective tree-temporal** | Fixes aggregation power, multiple testing and repeated looks at once. Hierarchy already in our data. Direct precedent on emerging variants. | Free binary, ~1 day |
| **2** | **SaTScan spatial variation in temporal trends** | Answers "where is it growing"; needs no denominator | Free binary, ~1 day |
| **3** | **Delay-distribution nowcast** | Replaces discarding 2026 with using it under uncertainty | ~2 days, in-house |
| 4 | `sparr` surface for PCP only | n = 367 clears the threshold; turns a ratio into a map | R detour |
| 5 | READUS-PV reporting conventions | Manuscript credibility | Writing only |
| — | **Keep** GPS/EBGM | Criticism targets DA *alone*; we triangulate | — |
| — | **Skip** MaxSPRT, LGCP/inlabru | Redundant / prior-dominated at our n | — |

The honest headline: **what is built here is a defensible instance of standard
practice, and its principal deficiency is that it tests one substance at a
time.** Removing that constraint is the single highest-value change available,
and it is most likely to change conclusions exactly where we are currently
blind — the fentanyl-analog and nitazene branches.
