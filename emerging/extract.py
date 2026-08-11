"""
Extract a decedent x substance table from LACME cause-of-death narratives.

The panel's 8 coded flags collapse everything outside
{meth, heroin, cocaine, fentanyl, alcohol, Rx opioids, benzos} into `Others`.
The narrative text does not: it names para-fluorofentanyl, bromazolam,
mitragynine, carfentanil, xylazine and ~100 more. This module recovers them.

Matching. Narratives are normalized and tokenized (see `lexicon.normalize`),
then scanned left-to-right with longest-n-gram-first lookup (n=4..1) against the
lexicon. Matched spans are consumed, so `PARA-FLUOROFENTANYL` wins over
`FENTANYL` on the same text and each token contributes to at most one substance.
Tokens that fail to match and contain a hyphen are re-tried split on the hyphen,
which recovers the ME's combination shorthand (`HEROIN-MORPHINE`,
`FENTANYL-METHAMPHETAMINE`).

Cohort. Uses the same definition as the modeling panel — `filter_alcohol_only`
from the monthly signal-check module (single source of truth, unified with the
epi project) and the EventZip -> ZIPCODE -> DeathZip fallback — so substance
counts are commensurable with `deaths` in `panel_zip_quarter.parquet`.

Three commands:
    extract   -- write the long decedent x substance table
    unmatched -- rank tokens the lexicon never consumed; the safety net for
                 substances NFLIS does not carry or the ME spells oddly
    validate  -- precision/recall of the extraction against the 8 coded flags

Usage:
    python -m rosla_nowcast.emerging.extract extract
    python -m rosla_nowcast.emerging.extract unmatched --top 80
    python -m rosla_nowcast.emerging.extract validate
"""

from __future__ import annotations

import collections
from pathlib import Path

import pandas as pd
import typer

from emerging.lexicon import (
    FUZZY_CUTOFF, MAX_NGRAM, OUT_DIR, Lexicon, fuzzy_resolve, load_lexicon,
    tokenize,
)
# Cohort definition shared with the modeling panel — do not re-implement.
from emerging.cohort import _norm_zip, filter_alcohol_only

app = typer.Typer(add_completion=False)

from emerging.paths import LACME_PATH, RESULTS_DIR  # noqa: E402
MENTIONS_PATH = OUT_DIR / "substance_mentions.parquet"
UNMATCHED_PATH = OUT_DIR / "unmatched_tokens.csv"

# Narrative fields scanned, in priority order. `Text` is the ME's combined
# cause string; the Cause* fields are the structured lines. Every field a
# substance matched in is recorded, so downstream can restrict to e.g. CauseA.
TEXT_FIELDS = ["CauseA", "CauseB", "CauseC", "CauseD", "CauseOther",
               "HowInjuryOccurred", "Text"]

# Coded flag -> the canonical substances that should light it up. Used only by
# `validate`, to score extraction against the ME's own coding.
FLAG_TO_CANONICAL: dict[str, set[str]] = {
    "Methamphetamine": {"Methamphetamine"},
    "Heroin": {"Heroin", "6-Monoacetylmorphine"},
    "Cocaine": {"Cocaine", "Cocaethylene", "Benzoylecgonine"},
    "Fentanyl": {"Fentanyl", "Norfentanyl"},
    "Alcohol": {"Ethanol"},
}
# Flags scored by NFLIS category rather than an enumerated set.
FLAG_TO_CATEGORY: dict[str, set[str]] = {
    "Benzodiazepines": {"Benzodiazepines"},
}
# `Fentanyl` scored a second way: does the ME's flag include analogs?
FENTANYL_ANALOG_CATEGORY = "Fentanyl and Fentanyl-related"


def load_cohort(
    lacme_path: Path = LACME_PATH, drop_chronic_alcohol: bool = True
) -> pd.DataFrame:
    """Load LACME overdose deaths with narrative text, zip, and quarter."""
    df = pd.read_csv(lacme_path, low_memory=False)
    if drop_chronic_alcohol:
        df = filter_alcohol_only(df)

    df["death_date"] = pd.to_datetime(df["DeathDate"], format="mixed", errors="coerce")
    df = df.dropna(subset=["death_date"])

    zc = _norm_zip(df["EventZip"])
    if "ZIPCODE" in df.columns:
        zc = zc.fillna(_norm_zip(df["ZIPCODE"]))
    zc = zc.fillna(_norm_zip(df["DeathZip"]))
    df["zipcode"] = zc

    df["quarter"] = df["death_date"].dt.to_period("Q").dt.to_timestamp()
    df["year"] = df["death_date"].dt.year
    return df.reset_index(drop=True)


class _FuzzyCache:
    """Memoized fuzzy alias lookup, with an audit trail of every substitution.

    Fuzzy matching is the one place the extractor guesses, so every guess it
    makes is recorded and written to `fuzzy_resolutions.csv` for review.
    """

    def __init__(self, lex: Lexicon, enabled: bool = True,
                 cutoff: float = FUZZY_CUTOFF):
        self.lex = lex
        self.enabled = enabled
        self.cutoff = cutoff
        self.aliases = tuple(lex.alias_to_canonical)
        # token -> (canonical, matched_alias) or None
        self._cache: dict[str, tuple[str, str] | None] = {}
        self.hits: collections.Counter = collections.Counter()

    def resolve(self, token: str) -> str | None:
        if not self.enabled:
            return None
        if token not in self._cache:
            alias = fuzzy_resolve(token, self.aliases, self.cutoff)
            self._cache[token] = (
                (self.lex.alias_to_canonical[alias], alias) if alias else None
            )
        entry = self._cache[token]
        if entry is None:
            return None
        canon, alias = entry
        self.hits[(token, alias, canon)] += 1
        return canon

    def audit(self) -> pd.DataFrame:
        return pd.DataFrame(
            [{"token": t, "matched_alias": a, "canonical": c, "n": n}
             for (t, a, c), n in self.hits.most_common()],
            columns=["token", "matched_alias", "canonical", "n"],
        )


def match_tokens(
    tokens: list[str], lex: Lexicon, fuzzy: _FuzzyCache | None = None
) -> tuple[set[str], list[str]]:
    """Greedy longest-first n-gram match over a token list.

    Returns (canonical substances found, tokens left unmatched).
    """
    found: set[str] = set()
    unmatched: list[str] = []
    i, n = 0, len(tokens)
    while i < n:
        hit = False
        for size in range(min(MAX_NGRAM, n - i), 0, -1):
            gram = " ".join(tokens[i:i + size])
            canon = lex.alias_to_canonical.get(gram)
            if canon is not None:
                found.add(canon)
                i += size
                hit = True
                break
        if hit:
            continue
        tok = tokens[i]
        # Combination and modifier shorthand: HEROIN-MORPHINE,
        # METHAMPHETAMINE-INDUCED, COCAINE-ASSOCIATED, ETHANOL-OLEANDRIN.
        # Keep whichever parts resolve. Requiring *all* parts to resolve (the
        # original rule) silently discarded a known substance whenever one
        # sibling token was a modifier or an uncatalogued compound.
        # Shredding a real hyphenated chemical name is not a risk here: whole
        # aliases are tried first and consumed, so PARA-FLUOROFENTANYL never
        # reaches this branch.
        if "-" in tok:
            resolved = [lex.alias_to_canonical.get(p) for p in tok.split("-") if p]
            hits = [r for r in resolved if r is not None]
            if hits:
                found.update(hits)
                i += 1
                continue
        # Last resort: ME typo. Guarded and audited — see _FuzzyCache.
        if fuzzy is not None:
            canon = fuzzy.resolve(tok)
            if canon is not None:
                found.add(canon)
                i += 1
                continue
        unmatched.append(tok)
        i += 1
    return found, unmatched


def extract_mentions(
    df: pd.DataFrame,
    lex: Lexicon,
    fields: list[str] = TEXT_FIELDS,
    fuzzy: _FuzzyCache | None = None,
) -> tuple[pd.DataFrame, collections.Counter]:
    """Scan each decedent's narrative fields; return a long table and the
    unmatched-token counter."""
    fields = [f for f in fields if f in df.columns]
    if not fields:
        raise ValueError(f"none of {TEXT_FIELDS} present in the LACME file")

    rows: list[dict] = []
    unmatched: collections.Counter = collections.Counter()

    for rec in df[["CaseNumber", "death_date", "quarter", "year", "zipcode", *fields]].itertuples(index=False):
        per_case: dict[str, set[str]] = {}
        for field in fields:
            text = getattr(rec, field, None)
            if not isinstance(text, str) or not text.strip():
                continue
            found, um = match_tokens(tokenize(text), lex, fuzzy)
            unmatched.update(um)
            for canon in found:
                per_case.setdefault(canon, set()).add(field)
        for canon, hit_fields in per_case.items():
            rows.append({
                "CaseNumber": rec.CaseNumber,
                "death_date": rec.death_date,
                "quarter": rec.quarter,
                "year": rec.year,
                "zipcode": rec.zipcode,
                "canonical": canon,
                "fields": ",".join(sorted(hit_fields)),
            })

    mentions = pd.DataFrame(rows)
    if mentions.empty:
        return mentions, unmatched

    meta_cols = ["category", "date_added", "date_added_censored",
                 "is_metabolite", "parent", "rollup", "is_class_term"]
    mentions = mentions.join(lex.meta[meta_cols], on="canonical")
    return mentions, unmatched


@app.command("extract")
def extract(
    lacme: Path = typer.Option(LACME_PATH, "--lacme"),
    out: Path = typer.Option(MENTIONS_PATH, "--out"),
    unmatched_out: Path = typer.Option(UNMATCHED_PATH, "--unmatched-out"),
    keep_chronic_alcohol: bool = typer.Option(
        False, "--keep-chronic-alcohol",
        help="Skip filter_alcohol_only (default: apply it, matching the panel)",
    ),
    no_fuzzy: bool = typer.Option(
        False, "--no-fuzzy", help="Disable the typo-tolerant fuzzy fallback"
    ),
    fuzzy_cutoff: float = typer.Option(
        FUZZY_CUTOFF, "--fuzzy-cutoff",
        help="difflib similarity required for a typo substitution",
    ),
) -> None:
    """Write the long decedent x substance table."""
    lex = load_lexicon()
    df = load_cohort(lacme, drop_chronic_alcohol=not keep_chronic_alcohol)
    typer.echo(f"cohort: {len(df):,} deaths, "
               f"{df['death_date'].min().date()} -> {df['death_date'].max().date()}")

    fuzzy = _FuzzyCache(lex, enabled=not no_fuzzy, cutoff=fuzzy_cutoff)
    mentions, unmatched = extract_mentions(df, lex, fuzzy=fuzzy)
    out.parent.mkdir(parents=True, exist_ok=True)
    mentions.to_parquet(out, index=False)

    um = (pd.DataFrame(unmatched.most_common(), columns=["token", "n"])
            .query("n >= 3"))
    um.to_csv(unmatched_out, index=False)

    audit = fuzzy.audit()
    if not audit.empty:
        audit_path = out.parent / "fuzzy_resolutions.csv"
        audit.to_csv(audit_path, index=False)
        typer.echo(f"fuzzy substitutions: {len(audit):,} distinct tokens, "
                   f"{int(audit['n'].sum()):,} mentions -> {audit_path}")

    subs = mentions[~mentions["is_class_term"]]
    per_case = subs.groupby("CaseNumber").size()
    covered = df["CaseNumber"].isin(set(subs["CaseNumber"]))
    typer.echo(f"mentions          : {len(mentions):,}")
    typer.echo(f"distinct canonicals: {subs['canonical'].nunique():,} "
               f"(+{int(mentions['is_class_term'].sum()):,} class-term mentions)")
    typer.echo(f"deaths with >=1 substance: {covered.sum():,} ({covered.mean():.1%})")
    typer.echo(f"substances per death: mean {per_case.mean():.2f}, max {per_case.max()}")
    typer.echo(f"unmatched tokens (n>=3): {len(um):,}")
    typer.echo(f"wrote {out}")
    typer.echo(f"wrote {unmatched_out}")


@app.command("unmatched")
def unmatched_cmd(
    path: Path = typer.Option(UNMATCHED_PATH, "--path"),
    top: int = typer.Option(60, "--top"),
) -> None:
    """Show the most frequent tokens the lexicon never consumed.

    This is the review surface for genuinely novel substances: anything the DEA
    catalog does not carry, or that the ME spells in a way SUPPLEMENT misses,
    lands here. Most rows will be anatomy and disease vocabulary.
    """
    if not path.exists():
        typer.echo(f"{path} not found — run `extract` first.")
        raise typer.Exit(1)
    um = pd.read_csv(path)
    typer.echo(um.head(top).to_string(index=False))


@app.command("validate")
def validate(
    lacme: Path = typer.Option(LACME_PATH, "--lacme"),
    mentions_path: Path = typer.Option(MENTIONS_PATH, "--mentions"),
    out_dir: Path = typer.Option(RESULTS_DIR, "--out-dir"),
) -> None:
    """Score the extraction against the ME's own 8 coded flags.

    The coded flags are the only ground truth available, so agreement with them
    is the credibility check for everything downstream. Disagreement is not
    automatically extractor error — the flags are themselves regex-derived — so
    both directions are reported.
    """
    lex = load_lexicon()
    df = load_cohort(lacme)
    mentions = pd.read_parquet(mentions_path)

    by_case: dict[str, set[str]] = (
        mentions.groupby("CaseNumber")["canonical"].agg(set).to_dict()
    )
    cat_by_canon = lex.meta["category"].to_dict()

    rows = []
    specs: list[tuple[str, set[str] | None, set[str] | None]] = [
        *[(f, cs, None) for f, cs in FLAG_TO_CANONICAL.items()],
        *[(f, None, cats) for f, cats in FLAG_TO_CATEGORY.items()],
        ("Fentanyl (incl. analogs)", {"Fentanyl", "Norfentanyl"},
         {FENTANYL_ANALOG_CATEGORY}),
    ]

    for flag, canon_set, cat_set in specs:
        raw_flag = "Fentanyl" if flag.startswith("Fentanyl (") else flag
        if raw_flag not in df.columns:
            continue
        truth = pd.to_numeric(df[raw_flag], errors="coerce").fillna(0).astype(bool)

        def _pred(case: str) -> bool:
            found = by_case.get(case, set())
            if canon_set and found & canon_set:
                return True
            if cat_set and any(cat_by_canon.get(c) in cat_set for c in found):
                return True
            return False

        pred = df["CaseNumber"].map(_pred).astype(bool)
        tp = int((truth & pred).sum())
        fp = int((~truth & pred).sum())
        fn = int((truth & ~pred).sum())
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        f1 = 2 * prec * rec / (prec + rec) if prec and rec else float("nan")
        rows.append({"flag": flag, "n_flagged": int(truth.sum()),
                     "n_extracted": int(pred.sum()), "tp": tp, "fp": fp, "fn": fn,
                     "precision": prec, "recall": rec, "f1": f1})

    res = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "flag_validation.csv", index=False)
    typer.echo(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    typer.echo(f"\nwrote {out_dir / 'flag_validation.csv'}")

    # What did the ME's catch-all `Others` flag actually contain? This is the
    # headline gain from extraction: naming the 4.8%.
    if "Others" in df.columns:
        others = df.loc[
            pd.to_numeric(df["Others"], errors="coerce").fillna(0).astype(bool),
            "CaseNumber",
        ]
        coded = set().union(*FLAG_TO_CANONICAL.values())
        cnt = collections.Counter(
            c for case in others for c in by_case.get(case, set())
            if c not in coded and not str(c).startswith("_class_")
        )
        top = pd.DataFrame(cnt.most_common(25), columns=["canonical", "n_in_Others"])
        top.to_csv(out_dir / "others_flag_contents.csv", index=False)
        typer.echo(f"\nTop substances inside the `Others` flag (n={len(others):,} deaths):")
        typer.echo(top.to_string(index=False))
        typer.echo(f"\nwrote {out_dir / 'others_flag_contents.csv'}")


if __name__ == "__main__":
    app()
