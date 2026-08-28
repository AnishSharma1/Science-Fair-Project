# Methods

## Eligibility

One biological system equals one paired human alpha-beta TCR or explicitly identified clone. Every strict system has functional recognition evidence for all declared arms, exact HLA alpha/beta allotypes, exact peptide sequences, resolved nine-residue registers, and an experimental structure for each pMHC arm. Hy.1B11 has two required self-microbial comparisons but one system vote. EBNA1-ANO2 remains prospective.

## Controls and leakage prevention

N3 comparators are selected before geometry is read by peptide-length difference, binding-percentile difference, a fixed seeded SHA-256 tie-break, and candidate ID. N3 rows are unknown recognition and are never used as specificity negatives. All future weight selection must hold out the complete biological system across PDB and AlphaFold layers.

## Features

The locked feature family contains exposed-position C-alpha RMSD, C-alpha-to-side-chain-centroid vector RMSD, a five-property physicochemical mismatch, and anchor-position C-alpha RMSD. Each feature becomes a within-panel average-tie rank percentile. Candidate weights are nonnegative quarter increments summing to one. The exposed-CA baseline remains frozen until held-out testing is complete.

## Trust rule

Every required pair must rank at most 3 in an evaluable PDB oracle panel and in both fixed AlphaFold seeds. A completed rank above 3 fails; a missing layer or seed is not evaluable. Only an overall pass permits separate within-HLA discovery reranking. Cross-allele consensus is prohibited.

Computational pMHC geometry prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS mechanism.
