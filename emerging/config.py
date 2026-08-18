"""
Study-design parameters shared by every analysis in this repository.

These are *data judgements*, not derived quantities. They were previously
defined three times over (`trends.py`, `geo.py`, `polysubstance.py`) with
identical values and separate comments, which meant the study's end-of-data
boundary could be advanced in one module and silently left behind in the other
two. Nothing would have failed; the modules would simply have started
answering questions about different windows.

`polysubstance.py` used to carry a note explaining that it duplicated `CUTOFF`
rather than importing it, to keep the dependency one-way. That concern is real
and this module is the answer to it: `config` imports nothing from the package,
so every analysis can depend on it without any of them depending on each other.

Anything that is meaningful only inside one method stays with that method --
`ALARM_THRESHOLD` is an EB05 reading line and lives in `trends.py`, the
projection and permutation settings are spatial and live in `core/spatial.py`.
Only parameters that must agree across analyses belong here.
"""

from __future__ import annotations

# Hard completeness cutoff: quarters ending after this date are dropped before
# anything is computed. Toxicology through end-2025 is considered settled;
# 2026 is not. This is a *data judgement*, not a derived quantity -- advance it
# by hand when the LACME extract is refreshed and the newer quarters have had
# time to adjudicate.
#
# `nowcast.py` exists to interrogate exactly this value: it measures the
# reporting-delay curve from extract vintages and reports whether the cutoff is
# conservative, correct, or optimistic. Advance it only after re-running
# `nowcast estimate`.
CUTOFF = "2025-12-31"

# Start of the case-control window used by the spatial analyses -- the era the
# ME's coded flags come from. Earlier deaths are in the cohort but not in the
# contrast.
SINCE = "2023-01-01"

# Additional reporting-lag haircut, in quarters, applied on top of the cutoff.
# Defaults to 0 because `CUTOFF` already encodes the settledness judgement --
# these are two ways of expressing the same correction, and stacking them
# discards a year of settled data. Set it only if you drop the cutoff.
LAG_QUARTERS = 0

# The comparison that defines "rising": the most recent RECENT_QUARTERS scored
# against the BASELINE_QUARTERS before them.
RECENT_QUARTERS = 4
BASELINE_QUARTERS = 8
