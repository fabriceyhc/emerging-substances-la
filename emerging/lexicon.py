"""
Substance lexicon: alias -> canonical substance, built from the NFLIS substance
catalog plus a curated supplement.

Why NFLIS. The LACME cause-of-death narrative names substances in free text
(`EFFECTS OF FENTANYL, METHAMPHETAMINE, AND LIDOCAINE`), but the panel only
carries 8 coded flags, with everything else collapsed into `Others` (~4.8% of
post-2022 deaths). To recover the rest we need a substance vocabulary. The DEA
NFLIS substance list gives 3,101 substances with synonyms, a drug category, and
a `Date Added to NFLIS-Drug` — an off-the-shelf lexicon plus a novelty axis.

**Read this before using `date_added` as a novelty score.** It records when a
substance entered the *NFLIS catalog* — i.e. when forensic labs could report it
— not when it hit the street. 923 of 3,101 substances carry `10-1998`, the
catalog's inception dump, so that value is a left-censoring floor, not a date.
para-Fluorofentanyl and carfentanil both sit at 10-1998 yet only appear in LA
County deaths from 2021 and 2024 respectively. `date_added_censored` marks these
rows. Treat NFLIS date as *chemical* novelty and the LACME first-seen quarter
(see `trends.py`) as *local epidemiological* novelty; the interesting quadrant is
where they disagree.

Matching design. Aliases are matched as whole word n-grams (n=4..1, longest
first) against normalized narrative text, never as substrings — `FENTANYL` must
not fire inside `NORFENTANYL`, and `PARA-FLUOROFENTANYL` must win over
`FENTANYL` on the same span. Aliases shorter than `--min-alias-len` are dropped
unless they appear in SHORT_ALLOW, because the NFLIS catalog contains 214
aliases of <=5 characters (`FA`, `DOC`, `DOB`, `AMB`, `CBC`) that would collide
with ordinary narrative words.

Usage:
    python -m rosla_nowcast.emerging.lexicon build
    python -m rosla_nowcast.emerging.lexicon inspect fentanyl
    python -m rosla_nowcast.emerging.lexicon inspect --alias METAMPHETAMINE
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False)

from emerging.paths import NFLIS_PATH, PROCESSED_DIR as OUT_DIR  # noqa: E402
LEXICON_PATH = OUT_DIR / "substance_lexicon.parquet"
SUBSTANCE_META_PATH = OUT_DIR / "substance_meta.parquet"

# NFLIS inception dump: every substance the catalog launched with carries this
# stamp, so it means "known by 1998-10", not "added in 1998-10".
NFLIS_INCEPTION = pd.Timestamp("1998-10-01")

MIN_ALIAS_LEN = 4

# Short aliases worth keeping despite MIN_ALIAS_LEN — unambiguous in an
# overdose narrative and actually used by the ME.
SHORT_ALLOW = {
    "PCP", "MDMA", "MDA", "MDPV", "GHB", "LSD", "THC", "DMT", "PMA", "PMMA",
    "METH", "ETOH", "6-MAM", "U4", "2C-B", "TFMPP", "BZP", "DXM", "5F-ADB",
}

# NFLIS aliases that are ordinary English or clinical words in this corpus and
# would generate false positives. Checked against the 400 most frequent corpus
# tokens; extend as `extract.py unmatched` / `--review` surfaces more.
STOPLIST = {
    "ICE", "OTHER", "WATER", "SPICE", "BATH SALTS", "BLUE", "GOLD", "ENERGY",
    "DOC", "DOB", "DOI", "DOP", "DOM", "AMB", "CBC", "CPP", "DMG", "AXID",
    "FA", "EPE", "APB", "AET", "AMT", "DET", "DNP", "DOET", "DPT", "DMA",
    # Digestive enzyme, never an overdose substance here, but a 0.93 fuzzy
    # neighbour of the very common clinical word PANCREATIC.
    "PANCREATIN",
}

# Substances the ME writes that NFLIS does not catalog, plus the spelling
# variants observed in this corpus. Mapped alias -> canonical. Misspellings are
# real and frequent: METAMPHETAMINE(52), METAHMPHETAMINE(32), FENATNYL(18).
SUPPLEMENT: dict[str, str] = {
    # --- ethanol: NFLIS has no plain "Ethanol" entry ---
    "ETHANOL": "Ethanol",
    "ALCOHOL": "Ethanol",
    "ETOH": "Ethanol",
    "ETHYL ALCOHOL": "Ethanol",
    "ISOPROPYL ALCOHOL": "Isopropanol",
    "ISOPROPANOL": "Isopropanol",
    # --- methamphetamine spelling variants ---
    "METAMPHETAMINE": "Methamphetamine",
    "METAHMPHETAMINE": "Methamphetamine",
    "METHAMPHETMINE": "Methamphetamine",
    "METHAMPHETAMNE": "Methamphetamine",
    "METHAPHETAMINE": "Methamphetamine",
    "METHAMPETAMINE": "Methamphetamine",
    "METHAMPHETAMIINE": "Methamphetamine",
    "MEHTAMPHETAMINE": "Methamphetamine",
    "METHAMPHEAMINE": "Methamphetamine",
    "METHAMPHETAMINES": "Methamphetamine",
    "METHAMPETAMINES": "Methamphetamine",
    "METH": "Methamphetamine",
    # --- fentanyl spelling variants ---
    "FENATNYL": "Fentanyl",
    "FETANYL": "Fentanyl",
    "FENTNAYL": "Fentanyl",
    "FENTYNAL": "Fentanyl",
    "FENANYL": "Fentanyl",
    # --- other common variants / trade names ---
    "AMPHETAMINES": "Amphetamine",
    "COCAINE HYDROCHLORIDE": "Cocaine",
    "ECSTASY": "MDMA",
    "KRATOM": "Mitragynine",
    "NITROUS OXIDE": "Nitrous oxide",
    "CARBON MONOXIDE": "Carbon monoxide",
    "DIFLUOROETHANE": "1,1-Difluoroethane",
    "MARIJUANA": "Cannabis",
    "CANNABIS": "Cannabis",
    "TYLENOL": "Acetaminophen",
    # --- spelled-out chemical names NFLIS only carries with a numeric prefix ---
    "METHYLENEDIOXYMETHAMPHETAMINE": "MDMA",
    "METHYLENEDIOXYAMPHETAMINE": "MDA",
    "METHYLENEDIOXYPYROVALERONE": "MDPV",
    "TETRAHYDROCANNABINOL": "Cannabis",
    "GAMMA-HYDROXYBUTYRATE": "GHB",
    "GAMMA-HYDROXYBUTYRIC ACID": "GHB",
    "HYDROXYBUTYRATE": "GHB",
    "ACETYLFENTANYL": "Acetyl fentanyl",
    "ACRYLFENTANYL": "Acryl fentanyl",
    "PARAFLUOROFENTANYL": "para-Fluorofentanyl",
    "PHENCYCLIDINE": "PCP",
    # --- substances the ME names that NFLIS does not catalog ---
    "NICOTINE": "Nicotine",
    "SALICYLATE": "Salicylate",
    "ASPIRIN": "Salicylate",
    "CAFFEINE": "Caffeine",
    "CHLOROETHANE": "Chloroethane",
    # --- exact entries that pre-empt a wrong fuzzy snap ---
    # Reviewed in fuzzy_resolutions.csv: without these, difflib sent
    # DIETHYLAMIDE (from "lysergic acid diethylamide") to Diethylamine,
    # CLORAZEPAM to Lorazepam rather than Clonazepam, and the MDMA fragment
    # DIOXYMETHAMPHETAMINE to Hydroxymethamphetamine.
    "LYSERGIC ACID DIETHYLAMIDE": "LSD",
    "DIETHYLAMIDE": "LSD",
    "CLORAZEPAM": "Clonazepam",
    "DIOXYMETHAMPHETAMINE": "MDMA",
    # Levomethorphan is vanishingly rare in overdose deaths; a bare
    # "METHORPHAN" token is dextromethorphan in every reviewed instance.
    "METHORPHAN": "Dextromethorphan",

    # --- the adjective form of alcohol ---
    # "ALCOHOLIC INTOXICATION" is the single largest extraction gap (67 deaths).
    # Only the two-word acute phrasings are mapped: a bare "ALCOHOLIC" would
    # also fire on ALCOHOLIC CIRRHOSIS / HEPATITIS / CARDIOMYOPATHY, which are
    # chronic conditions contributing to a death, not an acute substance.
    "ALCOHOLIC INTOXICATION": "Ethanol",
    "ALCOHOL INTOXICATION": "Ethanol",
    "ETHANOL INTOXICATION": "Ethanol",

    # --- brand names the examiner uses ---
    "XANAX": "Alprazolam",
    "OXYCONTIN": "Oxycodone",
    "PERCOCET": "Oxycodone",
    "VICODIN": "Hydrocodone",
    "VALIUM": "Diazepam",
    "KLONOPIN": "Clonazepam",
    "ATIVAN": "Lorazepam",
    "ADDERALL": "Amphetamine",
    "SOMA": "Carisoprodol",
    "SUBOXONE": "Buprenorphine",
    "BENADRYL": "Diphenhydramine",

    # --- poisons and plant toxins outside the NFLIS drug catalog ---
    # NFLIS catalogues controlled and diverted *drugs*; these turn up in
    # overdose casework but are not trafficked, so the catalog omits them.
    "ARSENIC": "Arsenic",
    "ZINC": "Zinc",
    "COLCHICINE": "Colchicine",
    "OLEANDRIN": "Oleandrin",
    "CERBERIN": "Cerberin",
    "NERIIFOLIN": "Neriifolin",
    "ANABASINE": "Anabasine",
    "HELIUM": "Helium",
    "AMYL NITRITE": "Amyl nitrite",
    "LABETALOL": "Labetalol",
    "LABATELOL": "Labetalol",

    # --- real drugs NFLIS omits that were silently going unmatched ---
    # Both were also being mis-snapped by the fuzzy layer onto obscure
    # near-neighbours (benztropine -> clobenztropine, loratadine ->
    # desloratadine). An explicit entry beats a lucky fuzzy hit.
    "BENZTROPINE": "Benztropine",
    "LORATADINE": "Loratadine",

    # --- examiner typos too far from the source to fuzzy-match safely ---
    # All verified by hand against the surrounding narrative. They sit at
    # difflib 0.86-0.88, below the 0.90 cutoff; lowering the cutoff to reach
    # them admits far worse (see FUZZY_BLOCK).
    "CARENTANYL": "Carfentanil",
    "FENTANLY": "Fentanyl",
    "FANTANYL": "Fentanyl",
    "FENTANYK": "Fentanyl",
    "FENTENYL": "Fentanyl",
    "FENTALYL": "Fentanyl",
    "COCAINETOXICITY": "Cocaine",
    "ACOHOLIC INTOXICATION": "Ethanol",
    "ACOHOL INTOXICATION": "Ethanol",
}

# Metabolites and analytic breakdown products. Kept as their own canonical rows
# (so nothing is silently merged) but tagged with a parent, so emergence
# analysis can collapse them — otherwise NORFENTANYL reads as its own rising
# substance when it is just evidence of fentanyl.
METABOLITE_PARENT: dict[str, str] = {
    "Norfentanyl": "Fentanyl",
    "6-Monoacetylmorphine": "Heroin",
    "7-Hydroxymitragynine": "Mitragynine",
    "Cocaethylene": "Cocaine",
    "Benzoylecgonine": "Cocaine",
    "Nordiazepam": "Diazepam",
    "Norchlorcyclizine": "Chlorcyclizine",
    "EDDP": "Methadone",
    "Norbuprenorphine": "Buprenorphine",
    "Nortramadol": "Tramadol",
    "O-Desmethyltramadol": "Tramadol",
    "7-Aminoclonazepam": "Clonazepam",
    "Norketamine": "Ketamine",
}

# Isomer / naming-convention collapses. Distinct from METABOLITE_PARENT: these
# are the *same* substance written two ways, not a parent-product relation.
#
# Fluorofentanyl is the motivating case and a cautionary one. LACME writes both
# `PARA-FLUOROFENTANYL` and a bare `FLUOROFENTANYL`; NFLIS catalogs both, so the
# extractor faithfully splits them. But no single death names both, the two
# series start in the same quarter (2021Q1), and the bare form takes over in
# 2023Q2 when the office changed its writing convention — the "para-" series
# appears to collapse to zero while the combined signal keeps climbing.
# para-Fluorofentanyl is the dominant fluorofentanyl isomer in the US supply, so
# the bare form is folded into it. Left unrolled this manufactures both a
# spurious decline and a spurious second emergence.
ISOMER_ROLLUP: dict[str, str] = {
    "Fluorofentanyl": "para-Fluorofentanyl",
}

# Class terms, not substances. Extracted separately so they never enter the
# emergence ranking as if they were compounds, but retained because
# "BENZODIAZEPINES" with no named member is itself informative.
CLASS_TERMS: dict[str, str] = {
    "OPIATE": "_class_opioid",
    "OPIATES": "_class_opioid",
    "OPIOID": "_class_opioid",
    "OPIOIDS": "_class_opioid",
    "BENZODIAZEPINE": "_class_benzodiazepine",
    "BENZODIAZEPINES": "_class_benzodiazepine",
    "BARBITURATE": "_class_barbiturate",
    "BARBITURATES": "_class_barbiturate",
    "AMPHETAMINE-TYPE": "_class_stimulant",
    "SYNTHETIC CANNABINOID": "_class_synthetic_cannabinoid",
    "SYNTHETIC CANNABINOIDS": "_class_synthetic_cannabinoid",
    "DESIGNER DRUG": "_class_designer",
    "DESIGNER DRUGS": "_class_designer",
}

_DASHES = dict.fromkeys(map(ord, "‐‑‒–—―−"), "-")
# A comma between digits belongs to a chemical name (1,2-Dibromo...); a comma
# anywhere else is a list separator.
_COMMA_IN_NUMBER = re.compile(r"(?<=\d),(?=\d)")
_KEEP = re.compile(r"[^A-Z0-9,\-' ]+")
_WS = re.compile(r"\s+")
# Narrative connectives that separate list items; dropped before n-gramming so
# "FENTANYL AND METHAMPHETAMINE" cannot form a spurious 3-gram.
CONNECTIVES = {"AND", "WITH", "OF", "THE", "PLUS", "OR"}

MAX_NGRAM = 4


def normalize(text: str) -> str:
    """Uppercase, strip punctuation that is not part of a chemical name, and
    collapse whitespace. Applied identically to narratives and to aliases so the
    two sides of the match are always in the same space."""
    if not isinstance(text, str):
        return ""
    t = text.upper().translate(_DASHES)
    t = _COMMA_IN_NUMBER.sub("\x00", t)  # protect in-number commas
    t = _KEEP.sub(" ", t)
    t = t.replace(",", " ").replace("\x00", ",")
    return _WS.sub(" ", t).strip()


def tokenize(text: str) -> list[str]:
    """Normalized text -> token list, connectives removed.

    Leading/trailing hyphens are stripped: the ME's list formatting produces
    tokens like `--METHAMPHETAMINE` that would otherwise never match.
    """
    out = []
    for w in normalize(text).split(" "):
        w = w.strip("-")
        if w and w not in CONNECTIVES:
            out.append(w)
    return out


# --- fuzzy fallback -------------------------------------------------------
# The corpus carries a long tail of ME typos — METHAMPEHTAMINE, METHAMPHTAMINE,
# MTHAMPHETAMINE, NETHAMPHETAMINE, ... — that no hand-written SUPPLEMENT will
# ever exhaust. Unmatched long tokens get one difflib pass against the alias
# set. Guards against snapping a real drug onto a similarly spelled other drug
# (ALPRAZOLAM/FLUALPRAZOLAM scores 0.87): high cutoff plus a length-delta cap.
FUZZY_MIN_LEN = 7
FUZZY_CUTOFF = 0.90
FUZZY_MAX_LEN_DELTA = 3

# Tokens that must never be fuzzy-resolved, whatever the cutoff. Each one was
# observed snapping onto a wrong-but-similar substance when the cutoff was
# tested at 0.85, and each is either a distinct real drug, a chemical class, or
# an ordinary English word:
#   ALCOHOLIC   -> ALCOHOL   would tag Ethanol on alcoholic cirrhosis /
#                            hepatitis / cardiomyopathy — chronic conditions,
#                            not acute substances. 104 mentions; the single
#                            most damaging false positive available.
#   CANNABINOID -> CANNABINOL is a class, not that compound
#   HEMOGLOBIN  -> METHEMOGLOBIN, METHYLENE -> METHYLONE, DECLINE -> DECALIN
# Anything mapped deliberately belongs in SUPPLEMENT, not here.
FUZZY_BLOCK = {
    "ALCOHOLIC", "ALCOHOLISM", "ALCOHOLIM", "ALCOHOHOL", "ETHANOLIS",
    "CANNABINOID", "CANNABINOIDS", "HEMOGLOBIN", "METHYLENE", "DECLINE",
    "PETAMINE", "METHAMPHETA", "ALPRAZOLAMAND",
}


def fuzzy_resolve(
    token: str, aliases: tuple[str, ...], cutoff: float = FUZZY_CUTOFF
) -> str | None:
    """Best fuzzy alias for `token`, or None. Caller memoizes."""
    import difflib

    if len(token) < FUZZY_MIN_LEN or token in FUZZY_BLOCK:
        return None
    pool = [a for a in aliases
            if abs(len(a) - len(token)) <= FUZZY_MAX_LEN_DELTA and " " not in a]
    hit = difflib.get_close_matches(token, pool, n=1, cutoff=cutoff)
    return hit[0] if hit else None


@dataclass
class Lexicon:
    """alias -> canonical map plus per-canonical metadata."""

    alias_to_canonical: dict[str, str]
    meta: pd.DataFrame  # indexed by canonical

    @property
    def max_ngram(self) -> int:
        return max((a.count(" ") + 1 for a in self.alias_to_canonical), default=1)

    def canonical(self, alias: str) -> str | None:
        return self.alias_to_canonical.get(normalize(alias))


def _clean_alias(raw: str) -> str:
    return normalize(str(raw))


def _load_nflis(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"Substance Name", "Synonyms", "NFLIS Detailed Drug Category",
                "Date Added to NFLIS-Drug"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"NFLIS file missing columns: {sorted(missing)}")
    df["date_added"] = pd.to_datetime(
        df["Date Added to NFLIS-Drug"], format="%m-%Y", errors="coerce"
    )
    return df


def build_lexicon(
    nflis_path: Path = NFLIS_PATH, min_alias_len: int = MIN_ALIAS_LEN
) -> Lexicon:
    """Assemble the alias->canonical map from NFLIS + SUPPLEMENT + CLASS_TERMS.

    Collision rule: when two canonicals claim the same alias, the shorter
    canonical name wins (so the alias `FENTANYL` resolves to `Fentanyl`, not to
    some `...fentanyl` analog that happens to list it as a synonym), and
    SUPPLEMENT always overrides NFLIS.
    """
    nflis = _load_nflis(nflis_path)

    alias_map: dict[str, str] = {}
    rows: list[dict] = []

    def _add(alias: str, canonical: str, *, force: bool = False) -> None:
        a = _clean_alias(alias)
        if not a or a in STOPLIST:
            return
        if len(a) < min_alias_len and a not in SHORT_ALLOW:
            return
        if a.count(" ") + 1 > MAX_NGRAM:
            return
        prev = alias_map.get(a)
        if prev is None or force or len(canonical) < len(prev):
            alias_map[a] = canonical

    for _, r in nflis.iterrows():
        canonical = str(r["Substance Name"]).strip()
        if not canonical or canonical.lower() == "nan":
            continue
        rows.append({
            "canonical": canonical,
            "category": r["NFLIS Detailed Drug Category"],
            "date_added": r["date_added"],
            "date_added_censored": bool(r["date_added"] == NFLIS_INCEPTION),
            "source": "nflis",
        })
        _add(canonical, canonical)
        if pd.notna(r["Synonyms"]):
            for syn in str(r["Synonyms"]).split(";"):
                _add(syn, canonical)

    meta = pd.DataFrame(rows).drop_duplicates(subset="canonical")

    # Supplement: overrides NFLIS, and contributes new canonicals for anything
    # the catalog lacks (Ethanol chief among them).
    supp_rows = []
    for alias, canonical in SUPPLEMENT.items():
        _add(alias, canonical, force=True)
        if canonical not in set(meta["canonical"]):
            supp_rows.append({
                "canonical": canonical, "category": "Supplement",
                "date_added": pd.NaT, "date_added_censored": False,
                "source": "supplement",
            })
    for alias, canonical in CLASS_TERMS.items():
        _add(alias, canonical, force=True)
        supp_rows.append({
            "canonical": canonical, "category": "Class term",
            "date_added": pd.NaT, "date_added_censored": False,
            "source": "class_term",
        })
    if supp_rows:
        meta = pd.concat(
            [meta, pd.DataFrame(supp_rows).drop_duplicates(subset="canonical")],
            ignore_index=True,
        ).drop_duplicates(subset="canonical")

    meta["is_class_term"] = meta["canonical"].str.startswith("_class_")
    meta["parent"] = meta["canonical"].map(METABOLITE_PARENT)
    meta["is_metabolite"] = meta["parent"].notna()
    # Analysis grain: collapse metabolites onto their parent and isomer/naming
    # variants onto the dominant form, so one real-world substance is one row.
    meta["rollup"] = (meta["canonical"].map(ISOMER_ROLLUP)
                      .fillna(meta["parent"])
                      .fillna(meta["canonical"]))
    meta = meta.set_index("canonical").sort_index()

    return Lexicon(alias_to_canonical=alias_map, meta=meta)


def load_lexicon(
    lexicon_path: Path = LEXICON_PATH, meta_path: Path = SUBSTANCE_META_PATH
) -> Lexicon:
    if not lexicon_path.exists() or not meta_path.exists():
        raise FileNotFoundError(
            f"{lexicon_path} / {meta_path} not found — run "
            "`python -m rosla_nowcast.emerging.lexicon build` first."
        )
    al = pd.read_parquet(lexicon_path)
    meta = pd.read_parquet(meta_path)
    return Lexicon(
        alias_to_canonical=dict(zip(al["alias"], al["canonical"])),
        meta=meta.set_index("canonical") if "canonical" in meta.columns else meta,
    )


@app.command("build")
def build(
    nflis: Path = typer.Option(NFLIS_PATH, "--nflis"),
    out_dir: Path = typer.Option(OUT_DIR, "--out-dir"),
    min_alias_len: int = typer.Option(MIN_ALIAS_LEN, "--min-alias-len"),
) -> None:
    """Build the lexicon and write it to parquet."""
    lex = build_lexicon(nflis, min_alias_len)
    out_dir.mkdir(parents=True, exist_ok=True)

    al = pd.DataFrame(
        sorted(lex.alias_to_canonical.items()), columns=["alias", "canonical"]
    )
    al.to_parquet(out_dir / LEXICON_PATH.name, index=False)
    lex.meta.reset_index().to_parquet(out_dir / SUBSTANCE_META_PATH.name, index=False)

    m = lex.meta
    typer.echo(f"aliases          : {len(al):,}")
    typer.echo(f"canonicals       : {len(m):,}")
    typer.echo(f"  from NFLIS     : {int((m['source'] == 'nflis').sum()):,}")
    typer.echo(f"  supplement     : {int((m['source'] == 'supplement').sum()):,}")
    typer.echo(f"  class terms    : {int((m['source'] == 'class_term').sum()):,}")
    typer.echo(f"metabolites      : {int(m['is_metabolite'].sum()):,}")
    typer.echo(
        f"date_added left-censored (=={NFLIS_INCEPTION.date()}): "
        f"{int(m['date_added_censored'].sum()):,} "
        f"({m['date_added_censored'].mean():.0%} of NFLIS rows)"
    )
    typer.echo(f"max alias n-gram : {lex.max_ngram}")
    typer.echo(f"wrote {out_dir / LEXICON_PATH.name}")
    typer.echo(f"wrote {out_dir / SUBSTANCE_META_PATH.name}")


@app.command("inspect")
def inspect(
    term: str = typer.Argument(..., help="Substring to look up (case-insensitive)"),
    alias: bool = typer.Option(False, "--alias", help="Exact alias lookup instead"),
) -> None:
    """Look up how a term resolves, and show its NFLIS metadata."""
    lex = load_lexicon()
    if alias:
        canon = lex.canonical(term)
        typer.echo(f"{term!r} -> {canon!r}")
        if canon is None:
            raise typer.Exit(1)
        matches = [canon]
    else:
        matches = [c for c in lex.meta.index if term.lower() in str(c).lower()]
        if not matches:
            typer.echo(f"no canonical matching {term!r}")
            raise typer.Exit(1)

    rev: dict[str, list[str]] = {}
    for a, c in lex.alias_to_canonical.items():
        rev.setdefault(c, []).append(a)

    for c in matches[:25]:
        r = lex.meta.loc[c]
        date = "" if pd.isna(r["date_added"]) else str(r["date_added"])[:7]
        cens = " (censored: NFLIS inception)" if r["date_added_censored"] else ""
        typer.echo(f"\n{c}")
        typer.echo(f"  category   : {r['category']}")
        typer.echo(f"  date_added : {date}{cens}")
        if r["is_metabolite"]:
            typer.echo(f"  metabolite of: {r['parent']}")
        typer.echo(f"  aliases    : {', '.join(sorted(rev.get(c, []))[:12])}")
    if len(matches) > 25:
        typer.echo(f"\n... and {len(matches) - 25} more")


if __name__ == "__main__":
    app()
