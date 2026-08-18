"""
The shared plot palette.

Seven hex codes and one colormap, previously defined identically in
`trends.py`, `geo.py` and `polysubstance.py`, and imported *from* `trends.py`
by five more modules. That import is what made rendering the relative-risk
map depend on the empirical-Bayes ranking machinery: `relrisk.py` reached into
a 1,300-line analysis module to learn that ink is #0b0b0b.

This module imports nothing from the package and nothing from matplotlib's
pyplot, so depending on it costs a figure nothing.
"""

from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap

# Validated categorical palette (dataviz skill reference instance; checked with
# scripts/validate_palette.js --mode light: all checks pass, contrast WARN
# relieved by direct labels on every mark).
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#b8b7b2"

# Sequential ramp for share/density fills. Light background -> brand blue ->
# near-black, so it stays legible under PowerNorm compression.
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f2f6fc", "#2a78d6", "#12365f"])
