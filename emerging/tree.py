"""
The substance hierarchy the tree-temporal scan runs over.

A tree scan can only find a branch that the tree encodes. That makes this file
the domain-knowledge core of `treescan.py`, and the place where a finding is
most likely to be created or destroyed by a judgement call rather than by data.
It is deliberately separated from the statistic so it can be edited, diffed and
argued about on its own.

**Why the NFLIS categories are not enough on their own.** The game plan called
for `root -> NFLIS category -> substance`, because the category field is already
in the data and fully populated. Building it showed that the categories are a
forensic-chemistry filing system, not a pharmacological one, and three of the
big ones are incoherent for our purpose:

  - `Depressants and Tranquilizers` holds PCP (724 mentions, a dissociative)
    next to xylazine (9, an alpha-2 agonist) and zolpidem. A branch test on it
    means nothing: PCP is 60% of the node and would drown any tranq signal.
  - `Narcotic Analgesics` mixes prescription opioids (morphine 1,124,
    oxycodone 455) with the nitazenes (3 mentions total) and U-47700. The
    nitazene family — the exact case that motivated a tree scan — is 0.1% of
    its own parent node, so aggregating to the NFLIS parent *loses* the signal
    rather than pooling it.
  - `Other substances` is a 71-member junk drawer holding lidocaine and
    levamisole (adulterants) alongside atorvastatin and vancomycin.

So the NFLIS categories are kept as nodes — they are the published taxonomy and
a reader can check them — and a **family layer** is added beside them for the
groupings the pharmacology actually implies.

**This is a DAG, not a tree.** Families cross NFLIS category boundaries
(xylazine is filed under `Depressants and Tranquilizers`, medetomidine would
be too, and both are also cutting agents), so a substance can have more than
one parent and nodes may overlap. Nothing in the scan statistic requires a
strict tree: the log-likelihood ratio is computed per node from that node's own
leaf set, and the Monte Carlo maximum absorbs however much the nodes correlate.
The strict-tree requirement belongs to TreeScan's *implementation*, not to the
method. Overlapping nodes cost a little multiplicity and buy the ability to ask
pharmacologically real questions.

Node collection:

    root                    every substance mention
    superclass              5 nodes, NFLIS category rolled up
    NFLIS category          23 nodes, the published taxonomy verbatim
    family                  the pharmacological groupings below
    substance               211 leaves

Nodes whose leaf set duplicates another node's are dropped, keeping the more
specific label — `Heroin` the NFLIS category and `Heroin` the substance are
one test, not two, and counting them twice would inflate the multiplicity
correction for no gain.

Usage:

    python -m emerging.tree show
    python -m emerging.tree show --level family
    python -m emerging.tree check
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd
import typer

from emerging.extract import MENTIONS_PATH

app = typer.Typer(add_completion=False)


# NFLIS category -> superclass. This layer inherits NFLIS's misfilings by
# construction (PCP arrives here as a depressant); the family layer below is
# the corrective, not this one.
SUPERCLASS: dict[str, str] = {
    "Fentanyl and Fentanyl-related": "Opioids",
    "Heroin": "Opioids",
    "Narcotic Analgesics": "Opioids",
    "Narcotics": "Opioids",
    "Cocaine Alkaloids": "Stimulants",
    "Phenethylamines (Amphetamines) (S)": "Stimulants",
    "Phenethylamines (Amphetamines) (H)": "Stimulants",
    "Phenethylamines (Other) (S)": "Stimulants",
    "Synthetic Cathinones": "Stimulants",
    "Stimulants": "Stimulants",
    "Benzodiazepines": "Depressants",
    "Depressants and Tranquilizers": "Depressants",
    "Phenethylamines (2C Series) (H)": "Hallucinogens",
    "Phenethylamines (Other) (H)": "Hallucinogens",
    "Hallucinogens": "Hallucinogens",
    "Psychedelics": "Hallucinogens",
    "Tryptamines": "Hallucinogens",
    "Synthetic Cannabinoids": "Hallucinogens",
    "Analgesics": "Medicines and other",
    "Antidepressants": "Medicines and other",
    "Other": "Medicines and other",
    "Other substances": "Medicines and other",
    "Supplement": "Medicines and other",
}

# Pharmacological families. Members are canonical names as `lexicon.py`
# resolves them; `tree check` reports any that never appear so the list stays
# honest rather than aspirational.
#
# Names that are absent from LA are kept on purpose (medetomidine, brorphine).
# A family node holding only absent members contributes nothing and is dropped,
# but the moment one of them is coded the branch already exists and the scan
# tests it in the same quarter — which is the whole point of a standing tree.
FAMILIES: dict[str, list[str]] = {
    "Prescription opioids": [
        "Morphine", "Oxycodone", "Hydrocodone", "Methadone", "Tramadol",
        "Codeine", "Hydromorphone", "Oxymorphone", "Buprenorphine",
        "Tapentadol", "Meperidine", "Propoxyphene", "Opium",
    ],
    "Fentanyl analogs": [
        "para-Fluorofentanyl", "Acetyl fentanyl", "Carfentanil",
        "Benzylfentanyl", "Furanyl fentanyl", "4-ANPP", "Acryl fentanyl",
        "Cyclopropyl fentanyl", "Butyryl fentanyl", "Valeryl fentanyl",
        "Methoxyacetyl fentanyl", "Fluorofentanyl",
    ],
    "Novel synthetic opioids": [
        "U-47700", "2-Methyl AP-237", "AP-237", "Brorphine", "MT-45",
        "U-48800", "Isotonitazene", "Metonitazene", "Protonitazene",
        "N-Pyrrolidino etonitazene", "Etodesnitazene", "Butonitazene",
    ],
    "Kratom": ["Mitragynine", "7-Hydroxymitragynine"],
    "Prescription benzodiazepines": [
        "Alprazolam", "Diazepam", "Clonazepam", "Lorazepam", "Temazepam",
        "Chlordiazepoxide", "Midazolam", "Oxazepam", "Flurazepam",
        "Estazolam", "Clobazam", "Triazolam",
    ],
    "Designer benzodiazepines": [
        "Etizolam", "Bromazolam", "Flualprazolam", "Flubromazolam",
        "Clonazolam", "Diclazepam", "Delorazepam", "Flubromazepam",
        "Adinazolam", "Nitrazolam", "Norflurazepam", "Deschloroetizolam",
        "Bentazepam", "Phenazepam",
    ],
    "Barbiturates": [
        "Pentobarbital", "Phenobarbital", "Butalbital", "Secobarbital",
        "Amobarbital",
    ],
    "Z-drugs": ["Zolpidem", "Eszopiclone", "Zaleplon", "Zopiclone"],
    "Gabapentinoids": ["Gabapentin", "Pregabalin"],
    "Dissociatives": [
        "PCP", "Ketamine", "Dextromethorphan", "Deoxymethoxetamine",
        "Methoxetamine", "2-Fluorodeschloroketamine",
    ],
    # The "tranq" story: xylazine's national displacement by medetomidine is
    # the single thing this project is most often asked about.
    "Veterinary sedatives": ["Xylazine", "Medetomidine", "Dexmedetomidine"],
    "Local anesthetics": [
        "Lidocaine", "Procaine", "Benzocaine", "Tetracaine", "Bupivacaine",
        "Mepivacaine",
    ],
    # Substances that are in the supply because someone put them there, not
    # because anyone sought them. The local anesthetics are also here: the
    # nodes are allowed to overlap.
    "Cutting agents and adulterants": [
        "Lidocaine", "Levamisole", "Chlorcyclizine", "Phenacetin", "Quinine",
        "Procaine", "Benzocaine", "Tetracaine", "Caffeine", "Mannitol",
        "Boric acid", "Xylazine", "Medetomidine", "Dexmedetomidine",
        "BTMPS (Bis(2,2,6,6-tetramethyl-4-piperidinyl) Sebacate)",
    ],
    "Antipsychotics": [
        "Quetiapine", "Olanzapine", "Clozapine", "Aripiprazole",
        "Risperidone", "Haloperidol", "Chlorpromazine", "Perphenazine",
        "Lurasidone", "Thiothixene", "Prochlorperazine", "Ziprasidone",
        "Paliperidone",
    ],
    "Antihistamines": [
        "Diphenhydramine", "Doxylamine", "Hydroxyzine", "Promethazine",
        "Chlorpheniramine", "Brompheniramine", "Loratadine", "Desloratadine",
        "Chlorcyclizine", "Cyclizine", "Meclizine",
    ],
    "Muscle relaxants": [
        "Cyclobenzaprine", "Carisoprodol", "Baclofen", "Methocarbamol",
        "Tizanidine", "Metaxalone", "Meprobamate", "Orphenadrine",
    ],
    "Inhalants": [
        "1,1-Difluoroethane", "Nitrous oxide", "Chloroethane", "Helium",
        "Amyl nitrite", "Butyl Nitrite", "Chloroform", "Toluene",
    ],
    # Not a drug-supply family. It is here as a control: if this branch fires,
    # the ME's panel changed, not the street supply.
    "Cardiovascular medications": [
        "Amlodipine", "Metoprolol", "Propranolol", "Atenolol", "Diltiazem",
        "Verapamil", "Labetalol", "Nifedipine", "Carvedilol", "Digoxin",
        "Benazepril", "Felodipine", "Flecainide", "Isosorbide Mononitrate",
        "Hydrochlorothiazide", "Atorvastatin", "Rosuvastatin", "Prazosin",
        "Terazosin", "Tamsulosin", "Warfarin", "Clonidine", "Guanfacine",
    ],
}

# Families defined by pattern rather than enumeration, because the analogs are
# coined faster than any list is maintained. Matched against the canonical
# name, case-insensitively.
FAMILY_PATTERNS: dict[str, str] = {
    "Nitazenes": r"nitaz",
    "Fentanyl analogs": r"fentanyl",   # unions with the literal list above
}

# Order of preference when two nodes hold exactly the same leaves. More
# specific wins, so the substance `Heroin` survives and the NFLIS category
# `Heroin` is dropped.
LEVEL_RANK = {"substance": 0, "family": 1, "category": 2, "superclass": 3,
              "root": 4}


@dataclass
class Node:
    """One testable node: a name, a level, and the leaves it contains."""

    name: str
    level: str
    leaves: frozenset[str]
    parents: tuple[str, ...] = field(default=())

    @property
    def n_leaves(self) -> int:
        return len(self.leaves)


def substance_categories(mentions_path=MENTIONS_PATH) -> pd.Series:
    """Canonical substance -> its NFLIS category, one row per substance.

    Cocaine is filed under both `Cocaine Alkaloids` (3,969 mentions) and
    `Narcotics` (2), so ties are broken by frequency and then alphabetically —
    a substance must sit in exactly one category or the category nodes stop
    being a partition.
    """
    m = pd.read_parquet(mentions_path)
    m = m[~m["is_class_term"].fillna(False)]
    sub = m["rollup"].fillna(m["canonical"])
    m = m.assign(_sub=sub).drop_duplicates(subset=["CaseNumber", "_sub"])
    freq = m.groupby(["_sub", "category"]).size().rename("n").reset_index()
    freq = freq.sort_values(["_sub", "n", "category"],
                            ascending=[True, False, True], kind="stable")
    return freq.drop_duplicates("_sub").set_index("_sub")["category"]


def build_nodes(cats: pd.Series) -> list[Node]:
    """Assemble the node collection from a substance -> category mapping."""
    subs = sorted(cats.index)
    sub_set = set(subs)

    nodes: list[Node] = []
    cat_of = cats.to_dict()
    sup_of = {s: SUPERCLASS.get(cat_of[s], "Medicines and other") for s in subs}

    nodes.append(Node("root", "root", frozenset(subs)))

    for sup in sorted(set(sup_of.values())):
        members = frozenset(s for s in subs if sup_of[s] == sup)
        nodes.append(Node(sup, "superclass", members, ("root",)))

    for cat in sorted(set(cat_of.values())):
        members = frozenset(s for s in subs if cat_of[s] == cat)
        nodes.append(Node(cat, "category", members,
                          (SUPERCLASS.get(cat, "Medicines and other"),)))

    fam_members: dict[str, set[str]] = {
        name: {s for s in listed if s in sub_set}
        for name, listed in FAMILIES.items()
    }
    for name, pattern in FAMILY_PATTERNS.items():
        rx = re.compile(pattern, re.I)
        fam_members.setdefault(name, set())
        fam_members[name] |= {s for s in subs if rx.search(s)}
    # `Fentanyl analogs` means the analogs, not fentanyl itself — the pattern
    # would otherwise swallow the parent drug and make the node a proxy for it.
    fam_members["Fentanyl analogs"].discard("Fentanyl")
    for name in sorted(fam_members):
        members = fam_members[name]
        if members:
            nodes.append(Node(name, "family", frozenset(members), ("root",)))

    for s in subs:
        nodes.append(Node(s, "substance", frozenset({s}),
                          (cat_of[s],)))

    return _dedupe(nodes)


def _dedupe(nodes: list[Node]) -> list[Node]:
    """Drop nodes whose leaf set duplicates a more specific node's."""
    best: dict[frozenset[str], Node] = {}
    for nd in nodes:
        prev = best.get(nd.leaves)
        if prev is None or LEVEL_RANK[nd.level] < LEVEL_RANK[prev.level]:
            best[nd.leaves] = nd
    keep = {id(nd) for nd in best.values()}
    return [nd for nd in nodes if id(nd) in keep]


def build(mentions_path=MENTIONS_PATH) -> tuple[list[Node], pd.Series]:
    """Convenience: (nodes, substance -> category)."""
    cats = substance_categories(mentions_path)
    return build_nodes(cats), cats


def unresolved_family_members(cats: pd.Series) -> dict[str, list[str]]:
    """Family members that never appear in the data, by family.

    Not an error — most are substances LA has not seen. Reported so the
    constant above stays a claim about pharmacology rather than a wish list.
    """
    sub_set = set(cats.index)
    return {name: sorted(s for s in listed if s not in sub_set)
            for name, listed in FAMILIES.items()
            if any(s not in sub_set for s in listed)}


@app.command("show")
def show(
    mentions: str = typer.Option(str(MENTIONS_PATH), "--mentions"),
    level: str = typer.Option("", "--level",
                              help="root|superclass|category|family|substance"),
    min_leaves: int = typer.Option(0, "--min-leaves"),
) -> None:
    """Print the node collection the scan will test."""
    nodes, cats = build(mentions)
    typer.echo(f"{len(cats)} substances -> {len(nodes)} testable nodes")
    by_level = pd.Series([nd.level for nd in nodes]).value_counts()
    for lv in ["root", "superclass", "category", "family", "substance"]:
        typer.echo(f"  {lv:12s} {int(by_level.get(lv, 0)):4d}")

    sel = [nd for nd in nodes
           if (not level or nd.level == level) and nd.n_leaves >= min_leaves]
    sel.sort(key=lambda nd: (LEVEL_RANK[nd.level], -nd.n_leaves, nd.name))
    typer.echo("")
    for nd in sel:
        if nd.level == "substance":
            typer.echo(f"  [{nd.level:10s}] {nd.name}")
        else:
            shown = sorted(nd.leaves)[:6]
            more = f", +{nd.n_leaves - 6} more" if nd.n_leaves > 6 else ""
            typer.echo(f"  [{nd.level:10s}] {nd.name}  ({nd.n_leaves}): "
                       f"{', '.join(shown)}{more}")


@app.command("check")
def check(mentions: str = typer.Option(str(MENTIONS_PATH), "--mentions")) -> None:
    """Report family members that never appear, and uncategorised substances."""
    nodes, cats = build(mentions)
    missing = unresolved_family_members(cats)
    typer.echo("Family members not present in LACME (absence, not error):")
    for name in sorted(missing):
        typer.echo(f"  {name}: {', '.join(missing[name])}")

    unmapped = sorted(set(cats.values) - set(SUPERCLASS))
    if unmapped:
        typer.echo("\nNFLIS categories with no superclass mapping "
                   "(fall through to 'Medicines and other'):")
        for c in unmapped:
            typer.echo(f"  {c}")
    else:
        typer.echo("\nEvery NFLIS category is mapped to a superclass.")

    dup = len(SUPERCLASS) + len(FAMILIES) + len(cats) + 1 + 5 - len(nodes)
    typer.echo(f"\n{len(nodes)} nodes after dropping {max(dup, 0)} duplicate "
               f"or empty node(s)")


if __name__ == "__main__":
    app()
