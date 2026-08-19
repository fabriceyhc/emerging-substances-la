"""
Cause-line substance ordering: the ME writes one combined cause string,
ordered by clinical significance rather than alphabetically. This turns that
string into a per-case ordered sequence, for anything that needs to read
significance from *position* rather than from presence.

Moved out of `polysubstance.py`, which used it for the cause-or-passenger
verdict, once `trends.py` needed the same ordering for a position-weighted
case definition (see `docs/GPS_V2_DESIGN.md` #1) -- two analyses reading
position, one scanner.
"""

from __future__ import annotations

import pandas as pd

from emerging.ingest.lexicon import MAX_NGRAM, load_lexicon, tokenize


def rollup_map(m: pd.DataFrame) -> dict[str, str]:
    """canonical -> rollup, so a CauseA scan matches the mentions table."""
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
