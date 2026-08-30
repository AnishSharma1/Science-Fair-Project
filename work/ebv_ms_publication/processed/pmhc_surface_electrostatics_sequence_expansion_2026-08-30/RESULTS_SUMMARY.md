# Results Summary

## Main finding

The corrected common-field electrostatic analysis evaluates each frozen candidate within its exact-HLA N3 panel. A primary full-pMHC rank from 1 to 3 is the locked rank-context support threshold; formal support also requires complete register QC.

| Target | HLA | Pair | Full-pMHC rank (eps-in 2) | eps-in 4 / 8 | HLA-subtracted rank | Rank-only context | Formal status |
|---|---|---|---:|---:|---:|---|---|
| HY03_SEQ_01 | HLA-DRB1*03:01 | `IVRQSRGDR` / `LLKDAIGEG` | 19/26 | 18/26 / 15/26 | 18/26 | `electrostatic_context_not_supportive` | `not_evaluable` |
| HY03_SEQ_02 | HLA-DRB1*03:01 | `VTLTSYWRR` / `IAIHHPWIR` | 24/26 | 25/26 / 25/26 | 13/26 | `electrostatic_context_not_supportive` | `not_evaluable` |
| HY08_SEQ_01 | HLA-DRB1*08:01 | `LRALLARSH` / `LEARLSRMH` | 13/26 | 11/26 / 10/26 | 24/26 | `electrostatic_context_not_supportive` | `not_evaluable` |
| HY08_SEQ_02 | HLA-DRB1*08:01 | `VRRRVLVQQ` / `FSRVVHLYR` | 24/26 | 24/26 / 24/26 | 17/26 | `electrostatic_context_not_supportive` | `not_evaluable` |
| HY13_SEQ_01 | HLA-DRB1*13:03 | `WMCMTVRHR` / `IICYNWLHR` | 12/26 | 12/26 / 12/26 | 9/26 | `electrostatic_context_not_supportive` | `not_evaluable` |
| HY15_SEQ_01 | HLA-DRB1*15:01 | `ILIYNGWYA` / `IAIHHPWIR` | 8/26 | 8/26 / 8/26 | 23/26 | `electrostatic_context_not_supportive` | `not_evaluable` |

## Interpretation

- These pairs remain **sequence-supported hypotheses**, not sequence-plus-electrostatics-supported leads.
- The electrostatic result should lower their priority relative to candidates that eventually show agreement across sequence, register, and a validated surface endpoint.
- It does not prove nonrecognition. N3 pairs have unknown TCR-recognition status, and modeled fields can miss induced fit, water networks, dynamics, and an actual TCR footprint.
- Both formal lead gates remain `not_evaluable` because their registers are not robust under the frozen V3 rule.
- The new electrostatic endpoint was not run on the three development controls, so it is not positive-control validated and cannot unlock discovery or specificity claims.

## QC correction

The initial target-surface point analysis was discarded after common-solvent-accessibility QC failed. The final analysis uses 25 position-matched P2/P3/P5/P7/P8 field points that pass in all 60 models per panel. The correction was geometry-only and did not use electrostatic scores.
