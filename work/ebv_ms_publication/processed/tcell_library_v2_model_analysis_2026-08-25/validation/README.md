# Gold-standard positive-control capture audit

## What counts

The denominator is locked before reading model ranks. A system must have exact peptide identities and HLA arms, recognition of both pMHCs by the same human T-cell clone, and experimental pMHC structures for both arms. Supportive T-cell studies, antibody-only pairs, protein-level associations, canonical tiles, and computational discoveries do not count.

The current denominator contains **one independent system**: Hy.2E11 recognition of DRB5*01:01-BALF5 and DRB1*15:01-MBP, anchored by PDB 1H15 and 1BX2.

## Frozen-method recovery

- Available-set capture@1: **2/2 seeds**.
- Available-set capture@3: **2/2 seeds**.
- Fully evaluable seeds passing the predeclared rule: **1/1**.
- Strict two-seed status: **not_evaluable_incomplete_calibration**.
- Model or score changed to fit the positive: **no**.

Seed 104729 is an available-set result only because controls are missing. Its rank cannot be relabeled as a formal 1-of-26 result. The successful capture shows that the frozen pMHC geometry method recognizes this established example; one independent system is not enough to estimate general sensitivity or validate new biological claims.

Recovery of an experimentally established positive control calibrates known-positive capture; it does not establish presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism for discovery candidates.
