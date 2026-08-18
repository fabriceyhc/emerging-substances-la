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
