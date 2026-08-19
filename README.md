# Emerging-substance detection (LA County)

Recovers named substances from LACME cause-of-death narratives and ranks which
ones are rising. Scoped to LA County only.

## Why this exists

The modeling panel carries 8 coded substance flags and collapses everything else
into `Others` (4.8% of post-2022 deaths, with no indication of what the "other"
was). The narrative text is far richer — `EFFECTS OF FENTANYL,
METHAMPHETAMINE, AND LIDOCAINE` — and names xylazine, bromazolam, carfentanil,
mitragynine, para-fluorofentanyl and ~200 more. Extraction turns 8 columns into
**219 substances at 99.8% decedent coverage**, which is the input every method
in this directory depends on.

## Layout

```
emerging/
  paths.py  config.py  viz.py  cli.py   shared: where, what window, what colour
  ingest/      cohort  lexicon  extract  census
  core/        tree  spatial  causeline  concentration
  analysis/    trends  treescan  geo  spacetime  svtt  relrisk  nowcast  polysubstance
  validation/  treescan_validate  ground_truth  benchmark
results/<analysis>/    tables and figures for that module, together
docs/findings/<analysis>.md    what it found and what to distrust
tests/test_<analysis>.py
```

One name per question, four ways: `emerging/analysis/svtt.py` ↔
`results/svtt/` ↔ `docs/findings/svtt.md` ↔ `tests/test_svtt.py`.

`config.py` holds the study-design parameters that must agree across analyses —
`CUTOFF` above all, which was previously defined three times. `viz.py` holds the
palette, which was also defined three times and imported from `trends.py` by
five more modules. Neither imports anything from the package, so an analysis can
depend on them without depending on another analysis.

## Pipeline

```bash
emerging pipeline              # every stage, in dependency order
emerging pipeline --dry-run    # print the order and stop
```

Or one at a time — `emerging --help` lists all of them:

```bash
emerging lexicon       build     # alias -> canonical map
emerging extract       extract   # decedent x substance table
emerging extract       validate  # score vs the 8 coded flags
emerging trends        rank      # ranking, default = weighted+role+dual+spatial+treescan-veto
emerging trends        backtest  # does the default detector fire?
emerging trends        alarms    # every plain-EB05 breach, ever  (run before plot)
emerging ground-truth  score     # alarm episodes vs known emergences
emerging benchmark     run       # every detector variant, side by side
emerging benchmark     plot
emerging trends        watch     # national watchlist vs LA
emerging trends        regime    # recording artifacts that move EB05
emerging trends        plot      # five figures
emerging tree          show      # the hierarchy the scan runs over
emerging tree          check     # family members absent from LA
emerging treescan      scan      # tree-temporal scan, current quarter
emerging treescan      backtest  # as-of sweep + head-to-head vs EB05
emerging treescan      plot
emerging treescan-validate compare   # vs the TreeScan C++ binary
emerging polysubstance profile   # cause vs passenger
emerging polysubstance plot
emerging geo           cluster   # is it localized?
emerging geo           plot
emerging census        fetch     # ACS population per zip per year
emerging spacetime     scan      # where, and is it growing there?
emerging spacetime     plot
emerging svtt          scan      # where is mortality rising fastest?
emerging svtt          profile
emerging svtt          plot
emerging nowcast       triangle  # reporting-delay curve from vintages
emerging nowcast       estimate  # is the CUTOFF defensible?
emerging nowcast       plot
emerging relrisk       surface   # relative-risk map (PCP by default)
emerging relrisk       plot
```

`trends alarms` must run before `trends plot` — it writes `alarm_history.csv`
and `eb05_sweep.csv`, which supply the breach annotations and the alarm-timeline
figure, and `plot` warns and skips them rather than recomputing a 45-quarter
sweep. `emerging pipeline` encodes that order so it cannot be got wrong by hand.

The old `python -m emerging.analysis.trends rank` form still works; every module
keeps its own Typer app.

| Module | Role |
|---|---|
| `ingest/lexicon.py` | 4,089 aliases → 3,127 canonicals from the DEA NFLIS catalog + a curated supplement (ethanol, ME misspellings, metabolite/isomer rollups) |
| `ingest/extract.py` | greedy longest-n-gram matcher over the narrative fields, with an audited fuzzy fallback for ME typos |
| `ingest/census.py` | ACS resident population per zip per year — the denominator the rest of the pipeline deliberately did without |
| `ingest/cohort.py` | overdose-death cohort definition — **a vendored copy**, see below |
| `core/tree.py` | the substance hierarchy — NFLIS categories plus a pharmacological family layer; **the domain-knowledge core**, a scan can only find branches this file encodes |
| `core/spatial.py` | projection, geocoded cohort, facility flag, and the case-fraction guard — shared by every spatial analysis |
| `core/causeline.py` | the ME's cause line turned into a per-case ordered substance sequence — position, not presence |
| `core/concentration.py` | spatial concentration of a substance's *recent* deaths as one standardized number, with closed-form null moments instead of a permutation test |
| `analysis/trends.py` | gamma-Poisson empirical-Bayes ranking, as-of backtest, alarm history, watchlist, recording diagnostics, figures |
| `analysis/treescan.py` | tree-temporal scan statistic: every node × every recent window, with a Monte Carlo p-value adjusted for the whole search |
| `analysis/polysubstance.py` | is a flagged substance a cause of death or a passenger — co-occurrence plus position on the ME's cause line |
| `analysis/geo.py` | case-control permutation test for spatial localization against the overdose-death background |
| `analysis/spacetime.py` | space-time scan over circles of zipcodes × recent quarters — separates "over-represented here" from "growing here" |
| `analysis/svtt.py` | spatial variation in temporal trends — where a *slope* differs from the county, not where a *level* is high |
| `analysis/nowcast.py` | reporting-delay curve from extract vintages, with a guard that refuses to apply it when the target extract does not obey it |
| `analysis/relrisk.py` | Kelsall–Diggle log relative-risk surface, CV-selected bandwidth, random-labelling tolerance contours |
| `validation/treescan_validate.py` | exports our tree/counts to TreeScan's input format and diffs the two implementations — needs the binary, which is not in the repo |
| `validation/ground_truth.py` | scores detected alarm episodes against hand-curated emergence intervals by temporal overlap |
| `validation/benchmark.py` | every detector variant on one table, read at a fixed line and at a matched alarm budget — including the tree scan read at its leaves, at its branches, and as a whole hierarchy |
| `config.py` | study-design parameters shared across analyses (`CUTOFF`, `SINCE`, the windows) |
| `viz.py` | the plot palette |
| `paths.py` | where the data lives; `results_dir()` gives each analysis its own output directory |

## Setup

```bash
pip install -e .
```

All commands are run from the repository root. Paths are fixed (`emerging/paths.py`):

| Path | In git? | Notes |
|---|---|---|
| `data/raw/lacme/` | **no** | LACME extract — **contains decedent names** |
| `data/processed/` | **no** | one row per decedent per substance; regenerate in ~1 min |
| `data/raw/nflis/` | yes | DEA NFLIS catalog, public |
| `data/raw/_geo/` | yes | LA County zipcode boundaries, public |
| `data/raw/census/` | yes | ACS population per zip per year, public aggregate |
| `results/<analysis>/` | yes | aggregate only — no case numbers, no coordinates |

**Do not add `.gitignore` exceptions for the first two.** Regenerate the
processed tables with `emerging lexicon build` then `emerging extract extract`.

### The cohort definition is duplicated, and that is a liability

`ingest/cohort.py` (`filter_alcohol_only`, `_norm_zip`) is a **vendored copy**. The
authoritative version lives in the `predict_rosla` repository at
`rosla_nowcast/modeling/zipcode/monthly/signal_check/ccf.py` and is itself kept
in sync with the epi project. It is duplicated here only because this
repository is standalone.

The point of the shared definition is that counts here stay commensurable with
`deaths` in `panel_zip_quarter.parquet`. If upstream changes and this copy does
not, every count silently stops matching and **nothing fails loudly**. Diff
before publishing any number that gets compared to the modeling work.
Copied 2026-08-10 from `predict_rosla` @ e357baf.

## Findings

Each analysis has a findings document — what it concluded, how it was
calibrated, and which artifacts would have manufactured a false result.

| | |
|---|---|
| [Extraction and validation](docs/findings/extract.md) | precision/recall against the ME's coded flags; why coverage is 99.8% |
| [EB05 ranking](docs/findings/trends.md) | does the ranking work, how the expected count is computed, four caveats that change conclusions |
| [Poisson assumption check](docs/findings/dispersion.md) | whether each substance's deaths arrive independently; 9 of 21 do not |
| [Tree-temporal scan](docs/findings/treescan.md) | scanning the hierarchy; validation against the TreeScan binary |
| [Cause or passenger](docs/findings/polysubstance.md) | co-occurrence and position on the cause line |
| [Spatial localization](docs/findings/geo.md) | is it localized, and the provenance of the method |
| [Space-time scan](docs/findings/spacetime.md) | where, and is it growing there |
| [Spatial variation in trends](docs/findings/svtt.md) | where mortality is rising fastest; the Pomona Valley |
| [Reporting delay](docs/findings/nowcast.md) | how complete the recent quarters are; the correction that is withheld |
| [Relative-risk surface](docs/findings/relrisk.md) | PCP; three bugs that each produced a plausible wrong map |
| [Ground-truth overlap](docs/findings/ground_truth.md) | alarm episodes against curated emergence intervals; the presence-vs-growth mismatch |
| [Spatial concentration](docs/findings/concentration.md) | the GPS v2 §3 score; exact null moments, and the count floor below which it says nothing |
| [Detector head-to-head](docs/findings/benchmark.md) | seven variants scored twice; why shrinkage wins and none of the three extensions does |

- [`docs/RESEARCH_LOG.md`](docs/RESEARCH_LOG.md) — current state, game plan, and
  the artifacts that would have manufactured false results
- [`docs/LITERATURE_REVIEW.md`](docs/LITERATURE_REVIEW.md) — rare-event detection
  and spatial clustering methods; what is better than this and what is not

## What this does not do

- **Detect a substance nobody assayed.** If the ME's panel does not test for it
  and the narrative does not name it, no count model recovers it. That bound is
  chemistry and lab scope, not statistics.
- **Spatial *detection*.** Case-control relative-risk surfaces (`sparr`) and
  Kulldorff scans need roughly n≥50–100 cases; only para-fluorofentanyl and
  mitragynine clear that. Spatial remains a characterization tool for
  substances the temporal detector already named — which is exactly what
  `geo.py` does, and it is deliberately downstream of `trends`, never a
  detector in its own right. `core/concentration.py` folds a *cheap* spatial
  term into the ranking itself (GPS v2 §3), and lowers that bar to n≥10 rather
  than removing it: xylazine, at 9 lifetime mentions, scores identically zero
  in all 45 quarters. The benchmark measures what the term buys where it does
  apply, which is about half a quarter of lead time on two substances.
  Xylazine is reachable by *aggregation* rather than by geography — the tree
  scan's `Cutting agents and adulterants` branch catches it at 21%
  attributable, which is the one mechanism in this repository that gets
  anything out of a 9-mention substance (`docs/findings/benchmark.md`).
- **Occupancy-detection modeling.** Separating "absent" from "not tested for"
  needs detection probability to vary across units. One ME office running one
  evolving panel gives no such variation — every zip shares the same lab. This
  needs multi-jurisdiction data to be identifiable.
- **Substance co-occurrence *graph*.** Post-2022 prevalence is methamphetamine
  54% / fentanyl 54% / everything else under 6% — a two-hub star graph where
  centrality trajectories carry little discriminating information, so it is
  not a detector. Pairwise co-occurrence is still highly informative as an
  *interpretive* layer, and that is what `polysubstance.py` uses it for.
- **Identify the policy behind a regime change.** `regime` flags the quarters
  where the denominator, the panel, or the onset pattern shifted. Matching
  those to naloxone saturation, a scheduling change, or an ME protocol
  revision needs sources outside this dataset.
- **Distinguish "absent" from "not assayed".** `watch` reports zeros for
  medetomidine, dexmedetomidine, BTMPS and desomorphine. Every one of those
  names resolves in the lexicon, so the scan looked and found nothing in the
  narratives — but an assay the ME does not run produces an identical zero.
  Only the ME's panel can separate them.
