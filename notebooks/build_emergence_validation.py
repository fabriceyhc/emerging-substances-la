"""
Regenerates notebooks/emergence_validation.ipynb.

Run from anywhere (`python3 notebooks/build_emergence_validation.py`):
Part 1 gets one markdown+code cell pair per substance that has ever fired
EB05, EB05+, TreeScan, NB-Trend or EB05++TS, each seeded with that
substance's `known_emergences.csv` interval where one exists. Re-run this
after a new substance fires (`emerging benchmark method-timeline` /
`emerging trends emergence-table` updated) or after editing
`known_emergences.csv` -- it overwrites the notebook wholesale, including
any `spans=[...]` edits made directly in a Part-1 cell that were never
promoted into `known_emergences.csv` (see the notebook's own top cell).
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def md(*lines):
    return {"cell_type": "markdown", "metadata": {}, "source": [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines else [])}

def code(*lines):
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [],
            "source": [l + "\n" for l in lines[:-1]] + ([lines[-1]] if lines else [])}

cells = []

cells.append(md(
"[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
"(https://colab.research.google.com/github/fabriceyhc/emerging-substances-la/blob/main/notebooks/emergence_validation.ipynb)",
"",
"# Emergence validation: annotated alarm review",
"",
"Reviews every substance that at least one of five detection methods has ever alarmed on, annotated against the same raw death-count line.",
"",
"**Data sources** (regenerate with the commands shown if stale):",
"",
"| file | command | grain |",
"|---|---|---|",
"| `results/trends/rising_substances.csv` | `emerging trends rank` | current EB05++TS score per substance, used only to order Part 1's cells and Part 2's grid |",
"| `results/trends/emergence_by_quarter.csv` | `emerging trends emergence-table` | substance x quarter, raw death counts |",
"| `data/raw/ground_truth/known_emergences.csv` | hand-curated (see `emerging/validation/ground_truth.py`) | the project's actual reference intervals -- seeds each Part 1 cell's `spans=[...]` where a substance already has one |",
"| `results/trends/emergence_by_year.csv` | `emerging trends emergence-table` | substance x year, raw death counts (the requested deliverable table; not plotted here, quarterly is used instead so alarm markers land on the right x-position) |",
"| `results/benchmark/method_timeline.csv` | `emerging benchmark method-timeline` | substance x as-of-quarter x method, score + alarm flag |",
"",
"**The six methods**, an escalating ladder of sophistication (`docs/GPS_V2_DESIGN.md`), not six arbitrary picks:",
"",
"| label | what it is | id in `emerging/validation/benchmark.py` |",
"|---|---|---|",
"| `n/E` | raw ratio, no shrinkage at all -- plotted for context but **excluded** from \"has this substance ever fired\" below: with no shrinkage at all it alarms on nearly every new arrival (872 of the sweep's ~3,700 substance-quarters), so using it as a filter would defeat the point of picking substances worth reviewing | `ratio` |",
"| `EB05` | classic single-gate EB05, raw counts | `eb05` |",
"| `EB05+` | dual-gated EB05 (share up **and** own-count up), raw counts | `eb05-dual` |",
"| `TreeScan` | substance-leaf recurrence interval, solo | `treescan` |",
"| `NB-Trend` | log-linear Poisson trend, Wald z of the fitted slope | `nb-trend` |",
"| `EB05++TS` | **the deployed detector** -- weighted + role-discounted + spatial-fused dual-gated EB05, TreeScan-vetoed at RI 10 | not a single benchmark.py id; read from `trends alarms`'s own cache (`eb05_sweep.csv`), the one reading here expensive enough that `alarms` already pays that cost |",
"",
"Every method's raw score is on an incomparable scale to the others (a ratio, a posterior percentile, a recurrence interval, a Wald z) -- this notebook never plots them against each other, only *when each one alarmed*, as markers on the one line that is comparable: the substance's own raw death count.",
"",
"**Part 1's cells are generated, not hand-written** -- one substance, one cell, each seeded with that substance's `known_emergences.csv` interval if it has one. Editing a cell's `spans=[...]` list and re-running it is the intended way to try out a candidate window visually; once you're satisfied a window is right, the place that opinion actually belongs is a row in `known_emergences.csv` itself (`docs/findings/ground_truth.md` documents the process), not this notebook -- regenerating Part 1 later (when a newly-fired substance needs adding) re-seeds every cell from that file and will overwrite in-notebook edits that were never promoted there.",
))

cells.append(code(
"from pathlib import Path",
"",
"import matplotlib.pyplot as plt",
"import pandas as pd",
"",
"%matplotlib inline",
"",
"# Locally (VS Code, JupyterLab) this notebook is assumed to stay at",
"# notebooks/emergence_validation.ipynb, one level under the repo root, and",
"# nothing below runs. Colab starts with an empty runtime with no checkout of",
"# this repo at all, so clone it there -- the three CSVs this notebook reads",
"# are aggregate counts already committed to the repo (`emerging/paths.py`'s",
"# rule: no case numbers, no coordinates in `results/`), so cloning is safe.",
"try:",
"    import google.colab",
"    IN_COLAB = True",
"except ImportError:",
"    IN_COLAB = False",
"",
"if IN_COLAB:",
"    import os",
"    REPO = Path(\"/content/emerging-substances-la\")",
"    if REPO.exists() and not (REPO / \".git\").exists():",
"        # A previous run's clone was interrupted (e.g. it hung and the cell",
"        # was stopped) and left a partial, non-git directory behind -- retry",
"        # from scratch rather than silently reusing a broken checkout.",
"        !rm -rf {REPO}",
"    if not REPO.exists():",
"        # Fail fast instead of hanging: with no credential helper configured,",
"        # a clone that unexpectedly needs auth would otherwise sit forever on",
"        # a username prompt Colab's `!` shell can never answer.",
"        os.environ[\"GIT_TERMINAL_PROMPT\"] = \"0\"",
"        !git clone -q https://github.com/fabriceyhc/emerging-substances-la.git {REPO}",
"    assert (REPO / \".git\").exists(), \"clone failed -- see the git output above\"",
"    ROOT = REPO",
"else:",
"    ROOT = Path(\"..\")",
"RESULTS = ROOT / \"results\"",
))

cells.append(code(
"timeline = pd.read_csv(RESULTS / \"benchmark\" / \"method_timeline.csv\", parse_dates=[\"as_of\"])",
"quarterly = pd.read_csv(RESULTS / \"trends\" / \"emergence_by_quarter.csv\", index_col=0)",
"quarterly.columns = pd.to_datetime(quarterly.columns)",
"ranking = pd.read_csv(RESULTS / \"trends\" / \"rising_substances.csv\", index_col=0)[\"eb05\"]",
"",
"# \"Fired\" = alarmed under at least one of the four non-ensemble methods or",
"# the deployed detector -- n/E excluded (see the table above). This is a",
"# strict superset of `alarm_history.csv`'s substances (which is EB05++TS",
"# alone): TreeScan and NB-Trend in particular each flag substances the",
"# deployed detector's own dual gate + veto never do, and part of this",
"# review's point is to see where the methods disagree, not just where the",
"# deployed one already agrees with itself.",
"ALARM_METHODS = [\"EB05\", \"EB05+\", \"TreeScan\", \"NB-Trend\", \"EB05++TS\"]",
"ever_alarmed = timeline.loc[timeline[\"method\"].isin(ALARM_METHODS) & timeline[\"alarm\"],",
"                           \"substance\"].unique()",
"",
"# Ordered by current EB05++TS score, descending -- same convention as Part",
"# 2's grid and `geo export`'s one-hot columns, so all three read substances",
"# in the same priority order. A few names here alarmed only under TreeScan",
"# or NB-Trend and never scored highly on EB05++TS itself, so they sit near",
"# the bottom of this list despite being exactly the disagreements worth",
"# reviewing -- that ordering choice, not a bug.",
"FIRED = (ranking.reindex(ever_alarmed)",
"         .sort_values(ascending=False, na_position=\"last\").index.tolist())",
"print(f\"{len(FIRED)} substances have ever fired one of {ALARM_METHODS}:\")",
"print(FIRED)",
))

cells.append(code(
"gt = pd.read_csv(ROOT / \"data\" / \"raw\" / \"ground_truth\" / \"known_emergences.csv\")",
"",
"",
"def _quarter(s):",
"    \"\"\"'2020Q2' -> Timestamp('2020-04-01'); blank/NaN -> None (open-ended).\"\"\"",
"    return None if pd.isna(s) else pd.Period(s, freq=\"Q\").start_time",
"",
"",
"# One row per substance in the CSV today, so each seed is a single-span list",
"# -- `plot_substance` itself takes a *list* of (start, end) pairs, because a",
"# substance can have more than one real episode (Acetyl fentanyl's notes",
"# flag a possible second one the current single-row schema can't hold yet).",
"KNOWN_SPANS = {r[\"substance\"]: [(_quarter(r[\"interval_start\"]), _quarter(r[\"interval_end\"]))]",
"               for _, r in gt[~gt[\"no_emergence\"]].iterrows()}",
"# Documented negatives: an interval is in the CSV only as a parsing",
"# placeholder (see known_emergences.csv's own docstring) and would mislead",
"# if shaded here as though it were a real reference window.",
"NEGATIVE_CASES = set(gt.loc[gt[\"no_emergence\"], \"substance\"])",
))

cells.append(md(
"## Part 1 -- annotated review: substances that fired",
"",
"One cell per substance below -- any substance at least one of EB05, EB05+, TreeScan, NB-Trend or EB05++TS has ever alarmed on, not only the deployed detector's own list (35 today; see the count printed above). The black line is the raw quarterly death count (`emergence_by_quarter.csv`); each colored row of markers above it is one method, placed at the as-of quarters where *that method's own* alarm condition was true (`method_timeline.csv`'s `alarm` column) -- not at a shared y-value tied to the line, so a quiet quarter's markers don't collide with the line itself. `n/E` still plots as context (grey dots) even though it doesn't gate which substances appear here.",
"",
"The grey shaded band(s), where present, are `spans` -- reference interval(s) either seeded from `known_emergences.csv` or left for you to fill in. Edit the list in a cell and re-run just that cell to try a different window; `end=None` means open-ended (shades through to the last available quarter).",
))

cells.append(code(
"METHOD_STYLE = {",
"    \"n/E\":      dict(marker=\"o\", color=\"#b8b7b2\"),",
"    \"EB05\":     dict(marker=\"s\", color=\"#2a78d6\"),",
"    \"EB05+\":    dict(marker=\"^\", color=\"#1baf7a\"),",
"    \"TreeScan\": dict(marker=\"D\", color=\"#eda100\"),",
"    \"NB-Trend\": dict(marker=\"v\", color=\"#7a4fd1\"),",
"    \"EB05++TS\": dict(marker=\"*\", color=\"#eb6834\"),",
"}",
"METHOD_ORDER = list(METHOD_STYLE)",
"",
"",
"def plot_substance(substance, spans=None):",
"    \"\"\"`spans`: a list of (start, end) reference intervals to shade, each",
"    a date-like string/Timestamp; `end=None` means open-ended, shaded",
"    through to the last available quarter.\"\"\"",
"    series = quarterly.loc[substance]",
"    fig, ax = plt.subplots(figsize=(11, 4.5))",
"",
"    for i, (start, end) in enumerate(spans or []):",
"        ax.axvspan(pd.Timestamp(start),",
"                  pd.Timestamp(end) if end is not None else series.index.max(),",
"                  color=\"#b8b7b2\", alpha=0.25, zorder=0,",
"                  label=\"reference interval\" if i == 0 else None)",
"",
"    ax.plot(series.index, series.values, color=\"#0b0b0b\", linewidth=1.6, zorder=2)",
"    ax.fill_between(series.index, series.values, color=\"#0b0b0b\", alpha=0.05, zorder=1)",
"",
"    sub = timeline[timeline[\"substance\"] == substance]",
"    ymax = max(series.max(), 1)",
"    for i, method in enumerate(METHOD_ORDER):",
"        fired_q = sub[(sub[\"method\"] == method) & sub[\"alarm\"]]",
"        if not len(fired_q):",
"            continue",
"        style = METHOD_STYLE[method]",
"        y = ymax * (1.06 + 0.05 * i)",
"        ax.scatter(fired_q[\"as_of\"], [y] * len(fired_q), label=method,",
"                   marker=style[\"marker\"], color=style[\"color\"],",
"                   edgecolors=\"white\", linewidths=0.6, s=45, zorder=3)",
"",
"    ax.set_title(f\"{substance} -- quarterly deaths, six methods' alarm quarters\",",
"                loc=\"left\", fontsize=12)",
"    ax.set_ylim(0, ymax * 1.4)",
"    for side in (\"top\", \"right\"):",
"        ax.spines[side].set_visible(False)",
"    ax.legend(frameon=False, ncol=7, loc=\"upper left\",",
"             bbox_to_anchor=(0, -0.08), fontsize=8.5)",
"    fig.tight_layout()",
"    plt.show()",
"    plt.close(fig)",
))

import pandas as _pd  # build-time only -- to compute the per-substance cells below

_timeline = _pd.read_csv(ROOT / "results" / "benchmark" / "method_timeline.csv", parse_dates=["as_of"])
_ranking = _pd.read_csv(ROOT / "results" / "trends" / "rising_substances.csv", index_col=0)["eb05"]
_gt = _pd.read_csv(ROOT / "data" / "raw" / "ground_truth" / "known_emergences.csv")

_ALARM_METHODS = ["EB05", "EB05+", "TreeScan", "NB-Trend", "EB05++TS"]
_ever_alarmed = _timeline.loc[_timeline["method"].isin(_ALARM_METHODS) & _timeline["alarm"],
                             "substance"].unique()
_FIRED = (_ranking.reindex(_ever_alarmed)
          .sort_values(ascending=False, na_position="last").index.tolist())


def _quarter_str(s):
    return None if _pd.isna(s) else str(_pd.Period(s, freq="Q").start_time.date())


_known = {r["substance"]: (_quarter_str(r["interval_start"]), _quarter_str(r["interval_end"]))
          for _, r in _gt[~_gt["no_emergence"]].iterrows()}
_negative = set(_gt.loc[_gt["no_emergence"], "substance"])

for _name in _FIRED:
    if _name in _known:
        start, end = _known[_name]
        note = f"`known_emergences.csv`: {start} to {end or 'ongoing'}"
        span_line = f'    ("{start}", {repr(end)}),'
    elif _name in _negative:
        note = "documented negative case in `known_emergences.csv` (`no_emergence=True`) -- no reference interval to shade"
        span_line = "    # documented negative case -- see the note above; no interval to seed"
    else:
        note = "not yet in `known_emergences.csv`"
        span_line = '    # ("YYYY-MM-DD", "YYYY-MM-DD"),  # add a candidate window, or leave empty'
    cells.append(md(f"### {_name}", "", f"*{note}*"))
    cells.append(code(
        f'plot_substance("{_name}", spans=[',
        span_line,
        "])",
    ))

cells.append(md(
"## Part 2 -- appendix: all 211 substances, unannotated",
"",
"Completeness/audit view, not a validation read: the ~176 substances not among Part 1's cells never alarmed under any of the five methods and are mostly flat or sparse lines with nothing to disagree about. No alarm markers here -- just the raw quarterly count, so a reviewer can confirm nothing in the long tail was overlooked, without 211 panels each carrying a legend.",
"",
"Ordered by *current* EB05++TS score, descending -- same convention as Part 1, `geo export`'s one-hot columns, and `emergence_by_year.csv`'s rows.",
))

cells.append(code(
"substances = quarterly.index.tolist()  # already current-EB05++TS-ordered",
"ncols = 14",
"nrows = -(-len(substances) // ncols)",
"",
"fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 1.15, nrows * 0.75))",
"for ax, name in zip(axes.flat, substances):",
"    s = quarterly.loc[name]",
"    ax.plot(s.index, s.values, color=\"#2a78d6\", linewidth=0.9)",
"    ax.fill_between(s.index, s.values, color=\"#2a78d6\", alpha=0.12)",
"    ax.set_title(name, fontsize=5.5, pad=1.5)",
"    ax.set_xticks([])",
"    ax.set_yticks([])",
"    for spine in ax.spines.values():",
"        spine.set_visible(False)",
"for ax in axes.flat[len(substances):]:",
"    ax.axis(\"off\")",
"",
"fig.suptitle(\"All 211 substances, raw quarterly deaths -- ordered by current \"",
"            \"EB05++TS score (audit view, no method annotation)\",",
"            fontsize=10, y=1.005)",
"fig.tight_layout()",
"plt.show()",
))

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out_path = ROOT / "notebooks" / "emergence_validation.ipynb"
with open(out_path, "w") as f:
    json.dump(nb, f, indent=1)
print(f"wrote {out_path}")
