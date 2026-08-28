# Expanded structural-control audit

## Current state

- Saved IEDB assay rows scanned: **1139**
- Deduplicated exact peptides in the frozen registry: **272**
- Peptides predicted for HLA-DRB1*15:01: **96**
- Selected stratum/layer mappings: **16**
- Unique AlphaFold Server jobs: **14**
- Complete exact-sequence AlphaFold jobs currently available: **14**
- Layered structural geometry rows currently available: **1050**

The discovery ranking is copied unchanged from the 2026-08-15 audit. Control selection did not inspect ranking, similarity, RMSD, AlphaFold confidence, or geometry fields.

## Frozen strata

- 21 aa: 5 selected; primary=5, binding sensitivity=0, length sensitivity=0, length+binding sensitivity=0, shortfall=0 (target_met).
- 23 aa: 5 selected; primary=5, binding sensitivity=0, length sensitivity=0, length+binding sensitivity=0, shortfall=0 (target_met).
- 25 aa: 3 selected; primary=3, binding sensitivity=0, length sensitivity=0, length+binding sensitivity=0, shortfall=2 (partial_no_further_relaxation).
- 32 aa: 3 selected; primary=0, binding sensitivity=0, length sensitivity=3, length+binding sensitivity=0, shortfall=2 (partial_no_further_relaxation).

The 32-aa primary stratum remains **not assessable** because the direct IEDB MHC-II range ends at 30 aa and no eligible 31-33 aa primary control exists. Its 25-30 aa controls remain a separate length-sensitivity layer.

## Registry audit counts

- eligible_pre_prediction: **96**
- excluded_invalid_coordinates: **16**
- excluded_mbp_plp_mog: **32**
- excluded_outside_predeclared_submission_lengths: **84**
- excluded_study_candidate: **44**

## Descriptive structural comparisons

A positive background-minus-target delta means the candidate pair has a lower exposed-position RMSD than the median of its equal-weighted controls. A negative value means the controls have the lower median RMSD. These are descriptive comparisons only.

| Discovery rank | Pair | Layer | Unique controls | Target median RMSD (A) | Background median RMSD (A) | Background - target (A) |
|---:|---|---|---:|---:|---:|---:|
| 2 | EBV_TCELL_2268741::HUMAN_MYELIN_117032 | 32-aa length sensitivity only | 3 | 1.208 | 13.974631 | 12.766631 |
| 4 | EBV_TCELL_119155::HUMAN_MYELIN_117032 | 32-aa length sensitivity only | 3 | 10.296 | 15.663122 | 5.367122 |
| 6 | EBV_TCELL_149795::HUMAN_MYELIN_117032 | 32-aa length sensitivity only | 3 | 13.033 | 15.716257 | 2.683257 |
| 9 | EBV_TCELL_119155::HUMAN_MYELIN_112037 | primary exact-bin, length +/-1 | 3 | 6.177 | 9.336278 | 3.159278 |
| 11 | EBV_TCELL_1862913::HUMAN_MYELIN_112037 | primary exact-bin, length +/-1 | 3 | 6.493 | 9.305639 | 2.812639 |
| 12 | EBV_TCELL_2268934::HUMAN_MYELIN_112226 | primary exact-bin, length +/-1 | 5 | 12.381 | 8.799549 | -3.581451 |
| 13 | EBV_TCELL_149795::HUMAN_MYELIN_116079 | primary exact-bin, length +/-1 | 5 | 13.341 | 13.379727 | 0.038727 |
| 14 | EBV_TCELL_1862913::HUMAN_MYELIN_117032 | 32-aa length sensitivity only | 3 | 13.964 | 18.800966 | 4.836965 |

## Interpretation limit

These controls provide descriptive pMHC structural context only. They do not establish peptide presentation, TCR binding, shared-TCR recognition, cross-reactivity, activation, molecular mimicry, or an MS mechanism. No p-value is reported.

## Reproduce

```bash
PYTHONPATH=src python3 src/build_structural_control_expansion.py --response-file processed/structural_control_expansion_2026-08-15/iedb_nextgen_mhcii_raw.json
```
