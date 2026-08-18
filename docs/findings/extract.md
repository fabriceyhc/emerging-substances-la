# Extraction and validation

Findings for `emerging/extract.py`. Results in [`results/extract/`](../../results/extract/).
Back to the [README](../../README.md); chronology in [RESEARCH_LOG.md](../RESEARCH_LOG.md).

## Validation

Extraction scored against the ME's own coded flags — the only ground truth
available (`results/extract/flag_validation.csv`):

| Flag | Precision | Recall | F1 |
|---|---|---|---|
| Fentanyl (incl. analogs) | 0.997 | 0.998 | 0.997 |
| Heroin | 0.995 | 0.998 | 0.996 |
| Methamphetamine | 0.993 | 0.989 | 0.991 |
| Cocaine | 0.985 | 0.997 | 0.991 |
| Alcohol | 0.986 | 0.962 | 0.974 |
| Benzodiazepines | 0.990 | 0.936 | 0.962 |

Disagreement is not automatically extractor error — the coded flags are
themselves regex-derived — so both directions are reported.

### Why coverage is 99.8% and not 100%

55 of 22,004 deaths (0.25%) get no named substance. Almost none of it is a
fixable extractor defect:

| Cause | n | Fixable? |
|---|---|---|
| The examiner named a **class, not a drug** — "opioid intoxication", "effects of opioids" | 34 | No. The record does not say which opioid. |
| **Chronic alcohol** phrasing — alcoholic hepatitis / pancreatitis / cardiomyopathy, chronic alcoholism | 14 | No, by design. These are chronic conditions contributing to a death, not acute substance involvement. |
| **Mid-word OCR damage** — `ALCOHOLIC INTOX ICATION`, `100COLCHI CINE` | 4 | Not without intra-word repair, which risks more than it recovers. |
| **Not overdoses at all** — anaphylaxis to ceftriaxone / Zosyn, pulmonary coccidioidomycosis | 3 | No — these are cohort-definition edge cases. |

Getting from 99.4% to 99.8% came from three fixes, all worth knowing about:

1. **The hyphen rule was too strict.** It only split a hyphenated token when
   *every* part resolved, so `METHAMPHETAMINE-INDUCED`, `COCAINE-ASSOCIATED`,
   `ETHANOL-OLEANDRIN` and `ALPRAZOLAM--HELIUM` all discarded a known substance
   because one sibling token was a modifier or an uncatalogued compound. It now
   keeps whichever parts resolve. This is safe because whole aliases are matched
   and consumed first, so `PARA-FLUOROFENTANYL` never reaches the split branch.
2. **`ALCOHOLIC INTOXICATION`** — the adjective form was the single largest gap
   (67 deaths). Only the two-word acute phrasings are mapped; a bare `ALCOHOLIC`
   would also fire on alcoholic cirrhosis and hepatitis.
3. **Brand names and non-catalog poisons** — Xanax, OxyContin, arsenic, zinc,
   colchicine, oleandrin. NFLIS catalogues trafficked *drugs*, so poisons and
   plant toxins that appear in overdose casework are simply absent from it.

**Do not fix the rest by loosening the fuzzy cutoff.** Tested at 0.85, it
recovers a few real typos and admits six false positives, the worst being bare
`ALCOHOLIC` → Ethanol at 104 mentions — which would have tagged alcoholic
cirrhosis as acute alcohol involvement and corrupted the ethanol series. Those
tokens are now in `FUZZY_BLOCK` so the mistake cannot be reintroduced by
changing one number. Verified typos go in `SUPPLEMENT` by name instead.

One of those by-name additions changed a result: `CARENTANYL EFFECTS` (2017-06)
is a real carfentanil death, moving carfentanil's first LA appearance from
2024Q1 to 2017Q2 — an isolated case during the national 2016–17 carfentanil
wave, then a seven-year gap before the 2024 cluster.

Two review surfaces, both worth reading before trusting a result:
`extract unmatched` ranks tokens the lexicon never consumed (where a genuinely
novel substance would first surface), and `fuzzy_resolutions.csv` logs every
typo substitution the matcher guessed at.
