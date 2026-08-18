# Relative-risk surface

Findings for `emerging/relrisk.py` (tests: `tests/test_relrisk.py`). Results in [`results/relrisk/`](../../results/relrisk/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## The relative-risk surface (PCP)

`geo.py` says PCP is 1.38× more clustered than the overdose background and
anchors on 90013; `spacetime.py` puts a 2.81× cylinder around 90002. Both are
summaries — one number, one circle. `relrisk.py` estimates the surface:

    r(x) = log( f(x) / g(x) )

for case density `f` (deaths naming the substance) and control density `g`
(the other overdose deaths). `r = 0` means "represented here exactly as it is
county-wide". Same compositional question as the rest of the pipeline,
answered continuously. Written in Python rather than driving R's `sparr`.

### Bandwidth is chosen, not assumed

`geo.py` rejected kernel methods partly because at small n the bandwidth drives
the answer — and it does: PCP's peak relative risk runs **3.7× at 1 km to 2.0×
at 6 km**. So the default comes from **likelihood cross-validation on the
case/control labels** (Kelsall–Diggle 1995, what `sparr::LSCV.risk` implements),
which targets the *ratio* rather than either density on its own.

It selects **4.0 km**, and the optimum is flat — 3.5, 4 and 5 km sit within two
log-likelihood units. That flatness is itself a result: PCP's excess is a broad
regional gradient, not a hotspot, and the data will not support finer
structure. `surface` prints the full sensitivity table next to the chosen value
so the reader sees what the choice bought. Silverman's normal-reference rule
(4.19 km) is reported for contrast and deliberately not used — it assumes a
Gaussian cloud, and LA is a 100 km multi-centre sprawl.

### The map and `geo.py` disagree, and both are right

`geo.py` names 90013 (Skid Row) as PCP's top zip. That is the **modal** zip by
raw count — 23 deaths. By **share of local overdose deaths** the picture is
different:

| Zip | | OD deaths | PCP | PCP share |
|---|---|---|---|---|
| 90044 | South LA | 66 | 10 | 15.2% |
| 90744 | Wilmington | 49 | 7 | 14.3% |
| 90806 | Long Beach | 65 | 9 | 13.8% |
| 90013 | Skid Row | 218 | 23 | 10.6% |
| — | county-wide | 7,581 | 367 | **4.84%** |

The surface puts the significant excess in a contiguous **South LA → Harbor
corridor** (90731, 90744, 90732, 90802, 90813, 90044) at about 2.1×, with the
valleys and the Antelope Valley depleted. Skid Row sits inside an elevated
downtown region but is not the peak: it has the most PCP deaths because it has
the most deaths.

12.6% of the county is above the county-wide rate at pointwise p ≤ 0.05. These
are **pointwise** contours over 41,000 cells with no multiplicity correction —
read contiguous regions, not pixels, and use `spacetime.py` for a single
adjusted verdict.

### Three bugs that each produced a plausible wrong map

- **`eps = 1e-12`** floored the log-ratio near −27 wherever a cell had zero
  cases, so a few empty-desert cells set the colour scale to ±20 and hid the
  real ±1 range. Replaced with a half-death pseudo-count centred so that
  `f/g = p_global` gives exactly 0 at any local density. This also regularised
  the upper tail — bandwidth sensitivity fell from 13.4×–2.0× to 3.7×–2.0×, so
  a good part of "the bandwidth drives everything" was really this bug.
- **No low-density mask**, so the ratio was estimated where almost no deaths
  informed it. Now requires ≥ 5 deaths within one bandwidth.
- **`gaussian_filter` returns the input dtype**, so an integer count grid comes
  back truncated and a smoothed density of 0.3 becomes 0. Caught by the
  calibration test firing at 84% instead of 5%; `log_rr` now coerces.
