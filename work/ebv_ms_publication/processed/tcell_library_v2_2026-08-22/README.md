# EBV-MS T-cell Library V2 (2026-08-22)

This additive package leaves all V1 analyses unchanged. It freezes a 40-EBV/40-self breadth panel, live IEDB register predictions, 320 AlphaFold Server inputs, and a separate 24-job native-HLA positive-control calibration.

## Current status

- Literature registry: 7 independent biological systems; one strict E1 system.
- Master library: 237 exact IEDB/literature records plus explicitly labeled canonical coverage tiles.
- Frozen panel: 40 EBV peptides across 20 proteins and 40 self peptides across 11 proteins.
- IEDB: 320 completed prediction records with exact `seq_num` and saved raw responses.
- AlphaFold: 320 discovery jobs and 24 calibration jobs are prepared, **not submitted**.
- Geometry/recovery: pending model downloads; no positive has been declared recovered.

Canonical tiling rows for protein/region-level evidence are not experimentally positive epitopes. Antibody-only EBNA1-GlialCAM is documented but excluded from the T-cell denominator. The EBNA1-MBP structural lead remains a computational hypothesis and is separate from positive recovery.

Computational pMHC geometry only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism.
