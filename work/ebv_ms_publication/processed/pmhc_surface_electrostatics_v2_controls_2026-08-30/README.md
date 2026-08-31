# Control-First Surface Electrostatics V2

This additive package replaces the earlier sparse field descriptor with a dense, model-specific near-surface pMHC map. It evaluates only the frozen Hy.2E11, Ob.1A12, and Hy.1B11 development controls. No discovery candidate was scored.

## Result

The locked control gate is **fail**. Candidate evaluation is **not allowed**. A failed gate permanently retires electrostatics from candidate ranking; a supportive gate would support only a supplementary descriptor.

## Locked Control Results

| Layer | Positive control | Peptide electrostatics | Composite electrostatics | Surface shape |
|---|---|---:|---:|---:|
| pdb | HY2E11_BALF5_MBP | rank 5 (fail) | rank 6 (fail) | rank 4 (fail) |
| pdb | OB1A12_ENGA_MBP | rank 6 (fail) | rank 6 (fail) | rank 5 (fail) |
| af_271828 | HY1B11_UL15_MBP | rank 12 (fail) | rank 10 (fail) | rank 10 (fail) |
| af_314159 | HY1B11_UL15_MBP | rank 3 (fail) | rank 2 (fail) | rank 6 (fail) |
| af_271828 | HY1B11_PMM_MBP | not evaluable | not evaluable | not evaluable |
| af_314159 | HY1B11_PMM_MBP | not evaluable | not evaluable | not evaluable |
| af_271828 | HY2E11_BALF5_MBP | not evaluable | not evaluable | not evaluable |
| af_314159 | HY2E11_BALF5_MBP | rank 18 (fail) | rank 6 (fail) | rank 1 (pass) |
| af_271828 | OB1A12_ENGA_MBP | not evaluable | not evaluable | not evaluable |
| af_314159 | OB1A12_ENGA_MBP | rank 1 (pass) | rank 1 (pass) | rank 10 (fail) |

`not evaluable` means at least one member of the complete comparison panel failed the predeclared 90% pairwise surface-map coverage rule. A positive row can therefore have a raw rank in the detailed tables while its gate result remains not evaluable. Completed rank, sensitivity, or resampling failures are retained as failures and make the overall gate fail.

The full pMHC nonlinear-PB primary calculation used protein dielectric 4, solvent dielectric 78.5, 0.15 M monovalent salt, pH 7.4 PARSE/PROPKA charges, and a 0.5 A maximum APBS grid spacing. Every required physical and sampling sensitivity is reported separately. There is no weighted composite.

Hy.1B11 remains unavailable as an experimental-PDB oracle because only two unique exact-HLA decoys were frozen. Its two AlphaFold panels remain mandatory. N3 and structural decoys are not specificity negatives.

Development-control pMHC surface resemblance only; not evidence of presentation, TCR recognition, activation, specificity, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
