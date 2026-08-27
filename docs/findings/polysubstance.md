# Cause or passenger

Findings for `emerging/polysubstance.py`. Results in [`results/polysubstance/`](../../results/polysubstance/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## Is it a cause of death, or a passenger?

`polysubstance.py`. EB05 measures growth in share; it cannot say whether
anyone died *of* the substance. Two independent discriminators, over the last
4 settled quarters:

- **Company kept.** `alone_pct` (deaths naming nothing else) and
  `with_fentanyl_pct`.
- **Position on the cause line.** The ME writes one combined string —
  `EFFECTS OF FENTANYL, METHAMPHETAMINE, AND LIDOCAINE` — ordered by clinical
  significance, not alphabetically. `lead_pct` is how often the substance is
  named first, `last_pct` how often it is named last, and `mean_rank_pct` its
  average position from 0 (always first) to 1 (always last). All three are
  computed only on multi-substance lines, so a solo mention cannot count as
  "named first".

| | n (2025) | alone | with fentanyl | 1st | last | mean rank | verdict |
|---|---|---|---|---|---|---|---|
| Fentanyl | 1101 | 15% | — | 82% | 9% | 0.13 | principal driver |
| Cocaine | 455 | 26% | 49% | 44% | 31% | 0.43 | principal driver |
| Methamphetamine | 1345 | 37% | 49% | 16% | 66% | 0.75 | independent |
| Alprazolam | 110 | 1% | 75% | 21% | 46% | 0.64 | adulterant/marker |
| PCP | 137 | 0% | 51% | 2% | 81% | 0.90 | secondary/incidental |
| **Lidocaine** | **17** | **0%** | **100%** | **0%** | **82%** | **0.93** | **adulterant/marker** |

Lidocaine is a cutting agent. Every one of its 17 recent deaths names
fentanyl, and the ME writes it **last on the cause line in 14 of 17 and first
in none** — the mirror image of fentanyl. The right headline is "the fentanyl
supply is being cut differently", not "lidocaine is killing people".

### Why position is normalized rather than a mean ordinal

A raw mean position (1…N) conflates "named late" with "appears on a long
line": a substance always named last scores 2 on a two-substance line and 5 on
a five-substance line. The bias runs the wrong way here — adulterants
co-occur with more substances and so sit on longer lines, inflating a raw
ordinal exactly for the class we are trying to identify.

Ranking all 26 substances both ways gives **Spearman 0.79**. Methamphetamine
is 20th of 26 raw but 7th normalized (short lines, so position 2 of 2 is
*last*); para-fluorofentanyl is 6th raw but 15th normalized (long lines, so
2.83 is the middle).

The cost: on a two-substance line the score can only be 0 or 1, so short lines
push toward the extremes; single-substance lines are undefined and excluded;
and it is an unweighted mean of ratios. `last_pct` is reported alongside
because the normalized mean is the right thing to rank on but harder to read.

*A discriminator that was tried and does not work:* which field the substance
was found in. Lidocaine is 100% `CauseA`, identical to fentanyl — the ME puts
everything on one line. Position within that line is the signal; presence in
it is not.

## Dyads: which pairs are enriched, and lift is not a danger score

`polysubstance.py dyads`. A different question from the one above: not "is
this substance a cause or a passenger" but "which pairs co-occur more than
their individual prevalences predict". Candidates are pharmacologically
motivated (drawn from the toxicology literature on additive/synergistic
CNS-respiratory depression, plus the cocaine+ethanol → cocaethylene metabolic
interaction), not mined from the data — the point is to test whether a
flagged combination actually shows up as enriched in LA deaths.

`lift = observed n_both / expected-under-independence`, computed over the
full settled history (2012–2025, N=21,717) for power, with a trailing
4-quarter read alongside so a dyad that has faded or just started
concentrating doesn't hide inside a 14-year average.

| dyad | n_both | lift | recent lift | mechanism |
|---|---|---|---|---|
| Fentanyl + Xylazine | 7 | 1.70 | 1.95 | alpha-2 sedation, naloxone-resistant ("tranq dope") |
| Methadone + Rx benzos | 35 | 1.62 | 2.52 | additive respiratory depression |
| Fentanyl + Designer benzos | 95 | 1.59 | 0.28 | same axis, unscheduled analogs |
| Fentanyl + Rx benzos | 741 | 1.26 | 1.42 | FDA black-box combination |
| Fentanyl + Cocaine | 1,904 | 1.06 | 0.95 | "speedball" |
| Fentanyl + Methamphetamine | 5,199 | 0.97 | 0.97 | "goofball" |
| Fentanyl + Gabapentinoids | 8 | 0.25 | 0.00 | 2019 FDA respiratory-depression warning |
| Fentanyl + Z-drugs | 4 | 0.21 | 0.00 | GABA-A axis, same as benzos |
| Fentanyl + Barbiturates | 2 | 0.09 | 0.00 | GABA-A axis, same as benzos |

**The result argues against a naive reading of "potentiation" as "these will
co-occur more."** Several combinations with real synergistic-depression
literature — gabapentinoids, Z-drugs, barbiturates, diphenhydramine — score
lift well *below* 1: they occur with fentanyl *less* than chance in this
corpus. And Fentanyl+Methamphetamine, the single largest polysubstance
combination in LA (5,199 deaths), sits at lift ≈ 1 — not sought together, just
both so prevalent post-2022 (54%/54% of deaths, see the README's "substance
co-occurrence graph" caveat) that independence alone predicts the overlap.
`lift` measures whether a pairing is concentrating in the local supply or
prescribing pattern; it is silent on how dangerous the pairing is on the
occasions it does occur, which is a pharmacological fact this dataset cannot
speak to. Read the two together: a high-lift dyad (xylazine, designer benzos)
is a supply-contamination signal worth tracking as its own trend; a
literature-flagged but lift<1 dyad (gabapentinoids, Z-drugs) is still
dangerous per-occurrence, it is just not the thing currently reshaping the
local supply.

Year-by-year, designer-benzo share of fentanyl deaths peaked in 2019–2020
(1.65–2.78%) and has not been rising since (0.09–1.18% 2023–2025); xylazine's
share of fentanyl deaths has stayed near zero (0–0.21%) every year on record —
both consistent with the xylazine finding already documented in
`RESEARCH_LOG.md` and `ground_truth.md`: LA lags the national tranq-dope
trend rather than following it.

## Metabolites vs. analogs

`lexicon.py`'s `METABOLITE_PARENT` collapses a metabolite onto its parent for
analysis grain — a metabolite-only mention would otherwise read as its own
(spurious) rising substance. This is deliberately **not** how fentanyl
analogs (para-fluorofentanyl, acetyl fentanyl, carfentanil, …) are handled:
an analog is an independently synthesized, independently psychoactive
compound sold in place of fentanyl, not a breakdown product of it, so each is
kept as its own canonical substance and can be tracked as its own emergence
(para-fluorofentanyl is one of the five ground-truth cases in
`ground_truth.md`). `tree.py`'s `Fentanyl analogs` family node aggregates
them for scan power but explicitly discards `Fentanyl` itself, so the
aggregate can't become a proxy for the parent drug.

The one metabolite rollup that didn't hold up against the corpus's own
evidence: **Nordiazepam** was mapped to `Diazepam`, but nordiazepam is the
shared active metabolite of diazepam *and* chlordiazepoxide, and of its 24
mentions only 1 also names Diazepam versus 4 that also name Chlordiazepoxide
(most name neither). There was no data-backed basis for the specific parent,
so the mapping was removed — Nordiazepam is now tracked as its own canonical
substance rather than silently mis-attributed.

The remaining twelve rollups (Norfentanyl→Fentanyl,
6-Monoacetylmorphine→Heroin, 7-Hydroxymitragynine→Mitragynine,
Cocaethylene/Benzoylecgonine→Cocaine, Norchlorcyclizine→Chlorcyclizine,
EDDP→Methadone, Norbuprenorphine→Buprenorphine,
Nortramadol/O-Desmethyltramadol→Tramadol, 7-Aminoclonazepam→Clonazepam,
Norketamine→Ketamine) each have exactly one possible parent in this corpus's
drug panel, so Nordiazepam's ambiguity doesn't apply. Checked directly: most
of them barely co-occur with their named parent at all — 7-Aminoclonazepam
names Clonazepam in 0 of 6 mentions, Norchlorcyclizine names Chlorcyclizine
in 0 of 7 — which is exactly why the rollup is load-bearing rather than
redundant: without it those deaths would be invisible to the parent's count,
not merely double-counted.
