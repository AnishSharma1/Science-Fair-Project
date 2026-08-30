# Methods

6 frozen sequence-supported targets were evaluated separately within their HLA alleles. Each target was compared with its previously frozen 5-by-5 exact-HLA N3 panel (25 comparator pairs). N3 denotes unknown TCR recognition and is not a specificity-negative class.

Five AlphaFold models per peptide arm were aligned to a panel reference using the first 85 C-alpha atoms of each HLA chain. PDB2PQR 3.7.1 assigned PARSE charges and radii after PROPKA titration at pH 7.4. APBS 3.4.1 solved the linearized Poisson--Boltzmann equation at 298.15 K, 0.15 M monovalent salt, solvent dielectric 78.5, solvent radius 1.4 A, and solute dielectric 2 (primary), 4, and 8 (sensitivity).

The final comparison shell contains five local samples for each declared TCR-facing position P2/P3/P5/P7/P8. Each point was shifted outward in 0.25 A increments until it cleared the 1.4 A probe surface by at least 0.25 A in every one of the 60 models in its panel. Full-pMHC and HLA-only potentials were sampled at the identical 25 coordinates. The primary metric is the 25th percentile Hodgkin similarity across all 25 EBV-model by self-model combinations; lower-quartile Carbo similarity, sign agreement, and upper-quartile potential RMSE are secondary.

This is a descriptive modeled-pMHC comparison. It does not measure presentation, TCR binding, activation, specificity, cross-reactivity, molecular mimicry, MS causation, probability, or false-discovery rate.
