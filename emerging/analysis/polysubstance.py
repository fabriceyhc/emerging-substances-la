"""
Classify flagged substances as causes of death or as passengers on somebody
else's overdose.

The EB05 ranking answers "is this substance's share of overdose deaths rising".
It cannot answer "would anyone have died of it". Lidocaine makes the gap
concrete: it ranks as a credible riser, and it appears in 17 of 17 recent deaths
*alongside fentanyl* and in 0 alone. It is a cutting agent. The correct headline
is "the fentanyl supply is being cut differently", not "lidocaine is killing
people" — and nothing in `trends.py` can tell those apart.

Two independent discriminators, both computed here:

**Co-occurrence.** For each substance, the share of its deaths that name no
other substance (`alone_pct`), the share naming fentanyl, and the lift of each
partner over that partner's own prevalence. A substance that is never alone and
almost always with one specific partner is a marker for that partner's supply.

**Order within CauseA.** The ME writes one combined cause line —
`EFFECTS OF FENTANYL, METHAMPHETAMINE, AND LIDOCAINE` — and orders it by
clinical significance, not alphabetically. Over the last four settled quarters
fentanyl is named first on 82% of multi-substance lines and last on 9%;
lidocaine is named first on 0% and last on 82%. So the ordinal rank of a
substance within its own cause line is a usable proxy for how the pathologist
weighted it. `lead_pct` is how often it is named first, `last_pct` how often
it is named last, and `mean_rank_pct` its average position on a 0 (always
first) to 1 (always last) scale.

**Why position is normalized rather than a mean ordinal.** Line length varies
with substance type, and a raw mean position conflates "named late" with
"appears on a long line" — a substance always named last scores 2 on a
two-substance line and 5 on a five-substance line. That bias runs the wrong
way for us: adulterants co-occur with more substances, so they sit on longer
lines, so a raw ordinal would inflate for exactly the class we are trying to
identify. Ranking all 26 substances both ways gives Spearman 0.79, and the
disagreements are large — methamphetamine is 20th of 26 on the raw scale and
7th normalized, because its lines are short so position 2 of 2 is *last*.

`last_pct` is reported alongside because the normalized mean, while the right
thing to rank on, is hard to read. "0% first, 82% last" states lidocaine's
case more plainly than 0.93 does. Note that on a two-substance line the
normalized score can only be 0 or 1, so short lines push it to the extremes.

**A discriminator that was tried and does not work:** which *field* the
substance was found in. CauseA vs CauseB vs CauseOther looks like it should
separate primary from incidental, and does not — lidocaine is 100% CauseA,
identical to fentanyl, because the ME puts every substance on the same line.
Position within that line is the signal; presence in it is not.

Verdicts are heuristic and stated as such: `principal driver`, `independent`,
`co-intoxicant`, `adulterant/marker`. They are a reading aid for the ranking
table, not a toxicological finding.

Usage:
    emerging polysubstance profile
    emerging polysubstance profile --quarters 8 --top 30
    emerging polysubstance plot
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer

from emerging.config import CUTOFF, RECENT_QUARTERS
from emerging.ingest.extract import MENTIONS_PATH, TEXT_FIELDS, load_cohort
from emerging.paths import results_dir
from emerging.ingest.lexicon import MAX_NGRAM, load_lexicon, tokenize
from emerging.viz import AQUA, BLUE, INK, INK_2, MUTED, ORANGE

app = typer.Typer(add_completion=False)

RESULTS_DIR = results_dir("polysubstance")

PROFILE_PATH = RESULTS_DIR / "profile.csv"

# Below this many cases in the window the percentages are not worth printing.
MIN_CASES = 4

# Verdict thresholds. Deliberately coarse — these separate obvious cases and
# leave everything else in the middle bucket rather than pretending to resolve
# it. See `_verdict`.
LEAD_PRINCIPAL = 40.0     # named first this often -> the ME's primary agent
ALONE_INDEPENDENT = 15.0  # kills without company this often -> stands alone
ALONE_PASSENGER = 5.0     # essentially never alone ...
RANK_PASSENGER = 0.60     # ... and habitually named late -> passenger
FENT_MARKER = 70.0        # ... and that late company is fentanyl -> adulterant


def _rollup_map(m: pd.DataFrame) -> dict[str, str]:
    """canonical -> rollup, so CauseA scan output matches the mentions table."""
    r = m.drop_duplicates(subset="canonical")
    return dict(zip(r["canonical"], r["rollup"].fillna(r["canonical"])))


def causea_order(
    cohort: pd.DataFrame, roll: dict[str, str], field: str = "CauseA"
) -> dict[str, list[str]]:
    """Case -> substances in the order the ME named them on the cause line.

    Runs the same greedy longest-n-gram-first scan as `extract`, but keeps
    token position instead of discarding it. No fuzzy fallback: this only has
    to *order* substances already known to be in the case, and a typo that
    fails to match here simply leaves that substance unranked rather than
    mis-ranked.
    """
    lex = load_lexicon()
    a2c = lex.alias_to_canonical
    out: dict[str, list[str]] = {}
    for case, text in zip(cohort["CaseNumber"], cohort[field].fillna("")):
        toks = tokenize(text)
        seen: list[str] = []
        i = 0
        while i < len(toks):
            for n in range(min(MAX_NGRAM, len(toks) - i), 0, -1):
                hit = a2c.get(" ".join(toks[i:i + n]))
                if hit is not None:
                    c = roll.get(hit, hit)
                    if c not in seen:
                        seen.append(c)
                    i += n
                    break
            else:
                i += 1
        out[case] = seen
    return out


def case_sets(
    mentions_path: Path = MENTIONS_PATH, cutoff: str = CUTOFF,
    quarters: int = RECENT_QUARTERS,
) -> tuple[pd.Series, list[pd.Timestamp]]:
    """Case -> set of substances, restricted to the last `quarters` settled."""
    m = pd.read_parquet(mentions_path)
    m = m[~m["is_class_term"].fillna(False)].copy()
    m["s"] = m["rollup"].fillna(m["canonical"])
    m = m.drop_duplicates(subset=["CaseNumber", "s"])
    m = m[m["quarter"] + pd.offsets.QuarterEnd(0) <= pd.Timestamp(cutoff)]

    keep_q = sorted(m["quarter"].unique())[-quarters:]
    m = m[m["quarter"].isin(keep_q)]
    return m.groupby("CaseNumber")["s"].apply(set), list(keep_q)


def _verdict(alone: float, lead: float, rank: float, fent: float) -> str:
    """Coarse label from aloneness, cause-line lead rate, and position.

    Order matters. `lead` is measured only on multi-substance cause lines (see
    `profile_table`), so winning it means the ME put this substance ahead of
    real competition — that outranks aloneness, which a substance can achieve
    simply by being the only thing anyone tested for.
    """
    if not np.isnan(lead) and lead >= LEAD_PRINCIPAL:
        return "principal driver"
    if alone >= ALONE_INDEPENDENT:
        return "independent"
    if alone < ALONE_PASSENGER and not np.isnan(rank) and rank >= RANK_PASSENGER:
        # Never alone and habitually named last. If the company it keeps is
        # overwhelmingly fentanyl, it is tracking the fentanyl supply; if not,
        # it is incidental to whatever else was on board.
        return ("adulterant/marker" if fent >= FENT_MARKER
                else "secondary/incidental")
    return "co-intoxicant"


def profile_table(
    sets: pd.Series, order: dict[str, list[str]], min_cases: int = MIN_CASES,
) -> pd.DataFrame:
    """One row per substance: aloneness, cause-line position, top partners."""
    n_all = len(sets)
    prevalence = (pd.Series([s for st in sets for s in st]).value_counts()
                  / n_all)

    rows = []
    for sub, prev in prevalence.items():
        cases = sets.index[sets.apply(lambda st: sub in st)]
        n = len(cases)
        if n < min_cases:
            continue

        alone = sum(1 for c in cases if len(sets[c]) == 1)

        # Ordinal position on the cause line, among the substances the scan
        # actually located there. Cases naming only one substance carry no
        # ordering information and are excluded from *both* statistics: a solo
        # mention is trivially "named first", and counting it inflates
        # lead_pct for exactly the solo-capable substances whose aloneness
        # `alone_pct` already reports. Scoring lead and rank on the same
        # multi-substance denominator is what keeps them consistent —
        # methamphetamine led 49% of lines while averaging position 0.75 until
        # this was fixed.
        leads, lasts, rank_pcts = 0, 0, []
        n_ranked = 0
        for c in cases:
            seq = order.get(c, [])
            if sub not in seq or len(seq) < 2:
                continue
            n_ranked += 1
            idx = seq.index(sub)
            if idx == 0:
                leads += 1
            if idx == len(seq) - 1:
                lasts += 1
            rank_pcts.append(idx / (len(seq) - 1))

        # sort_index() before value_counts() breaks ties by name rather than
        # by hash order, so `top_partners` is reproducible run to run. Without
        # it, two partners at the same count swap places between runs and the
        # CSV shows a spurious diff.
        partners = pd.Series(
            [x for c in cases for x in sets[c] if x != sub]
        ).value_counts().sort_index(kind="stable").sort_values(
            ascending=False, kind="stable")

        top = []
        for p, k in partners.head(3).items():
            lift = (k / n) / prevalence[p] if prevalence.get(p, 0) else np.nan
            top.append(f"{p} {100 * k / n:.0f}% (x{lift:.1f})")

        alone_pct = 100 * alone / n
        fent_pct = 100 * partners.get("Fentanyl", 0) / n
        lead_pct = 100 * leads / n_ranked if n_ranked else np.nan
        last_pct = 100 * lasts / n_ranked if n_ranked else np.nan
        mean_rank = float(np.mean(rank_pcts)) if rank_pcts else np.nan

        rows.append({
            "substance": sub, "n_cases": n,
            "alone_pct": alone_pct,
            "with_fentanyl_pct": fent_pct,
            "lead_pct": lead_pct,
            "last_pct": last_pct,
            "mean_rank_pct": mean_rank,
            "n_ranked": n_ranked,
            "verdict": _verdict(alone_pct, lead_pct, mean_rank, fent_pct),
            "top_partners": "; ".join(top),
        })

    # Substance name breaks ties on n_cases. Without it the row order for tied
    # counts is whatever `rows` happened to be built in, which varies run to
    # run (the co-occurrence sets upstream iterate in string-hash order), so
    # the committed CSV reshuffled on every regeneration and a git diff could
    # not distinguish a real change from noise. `kind="stable"` alone is not
    # enough -- it preserves an order that was itself arbitrary.
    return (pd.DataFrame(rows).set_index("substance")
            .sort_index()
            .sort_values("n_cases", ascending=False, kind="stable"))


@app.command("profile")
def profile(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    quarters: int = typer.Option(RECENT_QUARTERS, "--quarters",
                                 help="Trailing settled quarters to profile"),
    min_cases: int = typer.Option(MIN_CASES, "--min-cases"),
    top: int = typer.Option(30, "--top"),
) -> None:
    """Write the co-occurrence / cause-line-position profile."""
    sets, qs = case_sets(mentions, cutoff, quarters)
    cohort = load_cohort()
    m = pd.read_parquet(mentions)
    order = causea_order(cohort, _rollup_map(m))

    prof = profile_table(sets, order, min_cases)
    out_dir.mkdir(parents=True, exist_ok=True)
    prof.to_csv(out_dir / PROFILE_PATH.name)

    q = lambda t: f"{t.year}Q{t.quarter}"  # noqa: E731
    typer.echo(f"window {q(qs[0])}–{q(qs[-1])}: {len(sets):,} deaths naming "
               f"at least one substance")
    cols = ["n_cases", "alone_pct", "with_fentanyl_pct", "lead_pct",
            "last_pct", "mean_rank_pct", "verdict"]
    typer.echo("\n" + prof.head(top)[cols].to_string(
        float_format=lambda v: f"{v:.2f}"))
    typer.echo(f"\nwrote {out_dir / PROFILE_PATH.name}")


@app.command("plot")
def plot(
    mentions: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
    cutoff: str = typer.Option(CUTOFF, "--cutoff"),
    quarters: int = typer.Option(RECENT_QUARTERS, "--quarters"),
    min_cases: int = typer.Option(MIN_CASES, "--min-cases"),
    top: int = typer.Option(16, "--top"),
) -> None:
    """Two-panel figure: what each substance travels with, and where the ME
    ranks it on the cause line."""
    sets, qs = case_sets(mentions, cutoff, quarters)
    cohort = load_cohort()
    m = pd.read_parquet(mentions)
    order = causea_order(cohort, _rollup_map(m))
    prof = profile_table(sets, order, min_cases).head(top)
    prof = prof.iloc[::-1]  # largest at the top of a horizontal chart
    out_dir.mkdir(parents=True, exist_ok=True)

    # Composition per substance: alone / with fentanyl / with others only.
    # On the fentanyl row itself the middle category is meaningless — every
    # multi-substance fentanyl death trivially "contains fentanyl" — so it
    # collapses into "with others".
    comp = []
    for sub in prof.index:
        cases = sets.index[sets.apply(lambda st: sub in st)]
        n = len(cases)
        alone = sum(1 for c in cases if len(sets[c]) == 1)
        fent = 0 if sub == "Fentanyl" else sum(
            1 for c in cases if len(sets[c]) > 1 and "Fentanyl" in sets[c])
        comp.append((100 * alone / n, 100 * fent / n,
                     100 * (n - alone - fent) / n))
    comp = np.array(comp)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(14.5, 0.42 * len(prof) + 2.4),
        gridspec_kw={"width_ratios": [1.55, 1]})
    y = np.arange(len(prof))

    # --- left: composition of each substance's deaths --------------------
    left = np.zeros(len(prof))
    for vals, color, label in [
        (comp[:, 0], AQUA, "named alone"),
        (comp[:, 1], ORANGE, "with fentanyl"),
        (comp[:, 2], MUTED, "with others, no fentanyl"),
    ]:
        # 2px surface gap between segments so adjacent fills stay separable.
        ax1.barh(y, vals, left=left, height=0.62, color=color,
                 edgecolor="white", linewidth=1.4, zorder=3, label=label)
        left = left + vals
    for i, sub in enumerate(prof.index):
        ax1.annotate(f"n={prof.loc[sub, 'n_cases']}", (101, i), fontsize=8,
                     color=INK_2, va="center")
    ax1.set_yticks(y, prof.index, fontsize=9, color=INK)
    ax1.set_ylim(-0.7, len(prof) - 0.3)  # must match ax2 or the rows shear
    ax1.set_xlim(0, 112)
    ax1.set_xticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100%"],
                   fontsize=8, color=INK_2)
    ax1.set_xlabel("share of that substance's deaths", fontsize=9, color=INK_2)
    ax1.set_title("What it travels with", fontsize=10.5, color=INK, loc="left")
    ax1.legend(frameon=False, fontsize=8.5, ncol=3, loc="upper left",
               bbox_to_anchor=(0, -0.06))

    # --- right: position on the ME's cause line --------------------------
    r = prof["mean_rank_pct"].to_numpy()
    ok = ~np.isnan(r)
    ax2.axvspan(0, 0.33, color=ORANGE, alpha=0.08, linewidth=0, zorder=0)
    ax2.axvspan(RANK_PASSENGER, 1.0, color=MUTED, alpha=0.18, linewidth=0,
                zorder=0)
    ax2.hlines(y[ok], 0, r[ok], color=MUTED, linewidth=1, zorder=2)
    ax2.scatter(r[ok], y[ok], s=64, c=BLUE, edgecolors="white", linewidths=1.2,
                zorder=4)
    for i in np.where(ok)[0]:
        ax2.annotate(
            f"{prof['lead_pct'].iloc[i]:.0f}% first · "
            f"{prof['last_pct'].iloc[i]:.0f}% last",
            (r[i], y[i]), textcoords="offset points", xytext=(9, 0),
            fontsize=7.5, color=INK_2, va="center")
    # Verdict in a fixed right-hand column, so it reads as a label on the row
    # rather than as another value on the position scale.
    vcolor = {"principal driver": ORANGE, "independent": AQUA,
              "adulterant/marker": BLUE, "secondary/incidental": INK_2,
              "co-intoxicant": INK_2}
    for i, sub in enumerate(prof.index):
        v = prof.loc[sub, "verdict"]
        ax2.annotate(v, (1.78, y[i]), fontsize=7.5, color=vcolor[v],
                     va="center", annotation_clip=False)
    ax2.set_yticks(y, [""] * len(prof))
    ax2.set_ylim(-0.7, len(prof) - 0.3)
    ax2.set_xlim(0, 1.78)
    ax2.set_xticks([0, 0.5, 1.0],
                   ["named first", "middle", "named last"],
                   fontsize=8, color=INK_2)
    ax2.set_xlabel("mean position on the ME's combined cause line",
                   fontsize=9, color=INK_2)
    ax2.set_title("How the pathologist ranked it", fontsize=10.5, color=INK,
                  loc="left")

    for ax in (ax1, ax2):
        ax.grid(axis="x", color=MUTED, alpha=0.30, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(MUTED)
        ax.tick_params(length=2, colors=INK_2)

    q = lambda t: f"{t.year}Q{t.quarter}"  # noqa: E731
    fig.suptitle(
        "Cause of death or passenger? Substances by company kept and by "
        "cause-line rank\n"
        f"{q(qs[0])}–{q(qs[-1])}, {len(sets):,} LA County overdose deaths. "
        "A substance that is never alone and always named last is an "
        "adulterant marker, not a killer.",
        fontsize=11, color=INK, x=0.008, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "profile.png", dpi=150, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    typer.echo(f"wrote {out_dir / 'polysubstance.png'}")


if __name__ == "__main__":
    app()
