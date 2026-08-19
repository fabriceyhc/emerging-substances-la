"""
The single entry point: `emerging <module> <command>`.

Every analysis module builds its own `typer.Typer()`, which worked but meant
there was no way to see the 31 available commands without reading the README,
and no way to express the one ordering constraint the pipeline actually has
(`trends alarms` writes the sweep that `trends plot` annotates from). Both are
handled here.

The old invocations still work -- `emerging treescan scan` is
unchanged, because each module keeps its own `app`. This only adds a front
door.

Import cost: mounting every sub-app imports geopandas, matplotlib and scipy up
front, so `emerging --help` pays a few seconds. That is the price of one
`--help` that lists everything; the alternative is lazy sub-app construction,
which Typer does not support cleanly enough to be worth the indirection here.
"""

from __future__ import annotations

import typer

from emerging.analysis import (geo, nowcast, polysubstance, relrisk,
                                spacetime, svtt, treescan, trends)
from emerging.core import tree
from emerging.ingest import census, extract, lexicon
from emerging.validation import benchmark, ground_truth, treescan_validate

app = typer.Typer(
    add_completion=False,
    help="Emerging-substance detection in LA County overdose deaths.",
    no_args_is_help=True,
)

# Order follows the pipeline: build the lexicon, extract, rank, then the
# scans that take the extraction as given.
for _name, _mod, _help in [
    ("lexicon", lexicon, "alias -> canonical map from the NFLIS catalog"),
    ("extract", extract, "decedent x substance table, and its validation"),
    ("census", census, "ACS resident population per zip per year"),
    ("tree", tree, "the substance hierarchy the scans run over"),
    ("trends", trends, "gamma-Poisson empirical-Bayes ranking (EB05)"),
    ("treescan", treescan, "tree-temporal scan statistic"),
    ("treescan-validate", treescan_validate, "diff against the TreeScan binary"),
    ("polysubstance", polysubstance, "cause of death, or passenger"),
    ("geo", geo, "case-control test for spatial localization"),
    ("spacetime", spacetime, "space-time scan over zipcode circles"),
    ("svtt", svtt, "spatial variation in temporal trends"),
    ("nowcast", nowcast, "reporting-delay curve and completeness"),
    ("relrisk", relrisk, "Kelsall-Diggle relative-risk surface"),
    ("ground-truth", ground_truth, "score alarm episodes against known "
                                   "local emergences by interval overlap"),
    ("benchmark", benchmark, "every detector variant, scored side by side"),
]:
    app.add_typer(_mod.app, name=_name, help=_help)


# The pipeline, in dependency order. `trends alarms` must precede `trends plot`
# -- it writes alarm_history.csv and eb05_sweep.csv, which supply the breach
# annotations and the alarm-timeline figure. `plot` warns and skips them rather
# than recomputing a 45-quarter sweep, so running out of order silently yields
# a thinner figure. Encoding it here means it cannot be got wrong by hand.
PIPELINE: list[tuple[str, str]] = [
    ("lexicon", "build"),
    ("extract", "extract"),
    ("extract", "validate"),
    ("trends", "rank"),
    ("trends", "backtest"),
    ("trends", "alarms"),      # before `trends plot` and `ground-truth score`
    ("ground-truth", "score"),
    ("benchmark", "run"),
    ("benchmark", "plot"),
    ("trends", "watch"),
    ("trends", "regime"),
    ("trends", "plot"),
    ("tree", "show"),
    ("tree", "check"),
    ("treescan", "scan"),
    ("treescan", "backtest"),
    ("treescan", "plot"),
    ("polysubstance", "profile"),
    ("polysubstance", "plot"),
    ("geo", "cluster"),
    ("geo", "plot"),
    ("census", "fetch"),
    ("spacetime", "scan"),
    ("spacetime", "plot"),
    ("svtt", "scan"),
    ("svtt", "profile"),
    ("svtt", "plot"),
    ("nowcast", "triangle"),
    ("nowcast", "estimate"),
    ("nowcast", "plot"),
    ("relrisk", "surface"),
    ("relrisk", "plot"),
]


@app.command("pipeline")
def pipeline(
    dry_run: bool = typer.Option(False, "--dry-run",
                                 help="print the order and exit"),
    stop_on_error: bool = typer.Option(True, "--keep-going/--stop-on-error",
                                       help="halt at the first failure"),
) -> None:
    """Run every stage in dependency order.

    `treescan-validate compare` is deliberately excluded: it needs the TreeScan
    binary, which is not in the repository.
    """
    import subprocess
    import sys

    for group, cmd in PIPELINE:
        line = f"emerging {group} {cmd}"
        if dry_run:
            typer.echo(line)
            continue
        typer.echo(f"\n=== {line} ===")
        rc = subprocess.call([sys.executable, "-m", "emerging.cli", group, cmd])
        if rc != 0:
            typer.echo(f"FAILED ({rc}): {line}", err=True)
            if stop_on_error:
                raise typer.Exit(rc)


if __name__ == "__main__":
    app()
