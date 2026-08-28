# Methods

## Ranking scope

Each HLA allele is a separate 40-by-40 EBV-self screen. Complete pairs are ranked within that HLA by median exposed-position P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove fit, then ensemble IQR, then frozen pair ID. No geometry, rank, percentile, or score is pooled across alleles.

## Control recalibration

The endpoint was retained because it captured the locked Hy.2E11 BALF5-MBP gold-standard system at rank 1 in both available fixed-seed sets. Seed 104759 is the formal reference because it contains the complete positive plus all 25 predeclared score-blind full decoys. Seed 104729 remains an incomplete sensitivity analysis and does not set the control index.

The control separation index is `(decoy median - pair RMSD) / (decoy median - positive median)`. A value of 1 equals the modeled gold-positive median and 0 equals the full-decoy median. Values outside 0 to 1 are retained. This index is a descriptive method reference, not a calibrated probability, false-discovery rate, or biological validation score. Within-HLA rank remains the primary result.

## Missingness

Pairs without complete geometry are retained in `missing_unranked_pairs.csv` and excluded from the relevant HLA rank denominator. Missing results are not imputed.

## Interpretation

Computational pMHC geometry prioritization only; control-reference metrics are not probabilities and do not establish presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism.
