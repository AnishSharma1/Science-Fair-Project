# Methods

## Control-locked method selection

Only the five nonstructural methods from the completed held-out HLA-II benchmark v2 were eligible. Each biological system contributed one conservative score: its worst required positive rank across ligands and seeds. Methods were selected by system capture at 3, then worst system rank, system-weighted reciprocal rank, and a fixed lexical tie-break. Discovery files were read only after this selection. The selected method was TCR-facing BLOSUM62 similarity.

## Same-register ranking

Each pair's two predicted nine-residue HLA-II binding cores were compared directly P1-to-P1 through P9-to-P9. The primary score uses only P2/P3/P5/P7/P8. Every HLA allele is a separate 40-by-40 EBV-self ranking. Higher primary similarity is better; exact score ties retain a shared score rank and tie size, while lexical pair ID supplies a reproducible display order. Geometry and binding percentile do not enter the ranking.

## Validation status

The selected method recovered every positive panel at rank 1 across three independent control systems, but this is below the predeclared minimum of six systems for definitive validation. Results are therefore provisional and do not overwrite earlier structural rankings.

## Interpretation

Descriptive same-register pMHC sequence prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
