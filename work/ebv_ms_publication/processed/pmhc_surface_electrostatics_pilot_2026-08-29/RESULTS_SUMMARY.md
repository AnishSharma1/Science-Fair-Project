# Results Summary

## Main finding

The corrected common-field electrostatic analysis **does not support either BALF5--TALDO1 lead as unusually similar within its frozen exact-HLA N3 panel**.

| HLA | Pair | Full-pMHC rank (eps-in 2) | eps-in 4 / 8 | HLA-subtracted rank | Rank-only context | Formal status |
|---|---|---:|---:|---:|---|---|
| DRB1*13:03 | `YHFVKKHVH` / `LSFDKDAMV` | 21/26 | 21/26 / 21/26 | 25/26 | `electrostatic_context_not_supportive` | `not_evaluable` |
| DRB1*15:01 | `VYHFVKKHV` / `IYNYYKKFS` | 19/26 | 19/26 / 18/26 | 5/26 | `electrostatic_context_not_supportive` | `not_evaluable` |

The unfavorable full-pMHC rank class is stable across solute dielectric values 2, 4, and 8 for both leads. DRB1*15:01 improves to rank 5 after subtracting the modeled HLA field, but it still misses the frozen top-three support rule. DRB1*13:03 does not improve after subtraction.

## Interpretation

- These pairs remain **sequence-supported hypotheses**, not sequence-plus-electrostatics-supported leads.
- The electrostatic result should lower their priority relative to candidates that eventually show agreement across sequence, register, and a validated surface endpoint.
- It does not prove nonrecognition. N3 pairs have unknown TCR-recognition status, and modeled fields can miss induced fit, water networks, dynamics, and an actual TCR footprint.
- Both formal lead gates remain `not_evaluable` because their registers are not robust under the frozen V3 rule.
- The new electrostatic endpoint was not run on the three development controls, so it is not positive-control validated and cannot unlock discovery or specificity claims.

## QC correction

The initial target-surface point analysis was discarded after common-solvent-accessibility QC failed. The final analysis uses 25 position-matched P2/P3/P5/P7/P8 field points that pass in all 60 models per panel. The correction was geometry-only and did not use electrostatic scores.
