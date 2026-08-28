# Methods

## Ranking endpoint

Both input universes were rescored with normalized BLOSUM62 similarity at predicted HLA-II core positions P2/P3/P5/P7/P8. Higher values rank first. This endpoint was selected before the discovery merge because it placed every declared positive first in all eight panels of the three-system held-out control benchmark. Binding percentile, AlphaFold geometry, original local alignments, and original rank numbers do not enter the combined score.

## Legacy eligibility

An older pair enters only when both P1-P9 cores contain exactly nine residues, the EBV core is a primary HLA-DRB1*15:01 prediction, and the self core is either a primary HLA-DRB1*15:01 prediction or the experimental primary-allele reference. DRB5 calibration rows and unresolved or sensitivity-only registers remain in the eligibility audit but cannot enter this DRB1*15:01 rank.

## Merge and ties

Exact full-peptide duplicates are retained once after proving that both predicted cores agree; V2 is the canonical record and the legacy identifier remains linked. Every remaining unique pair is ranked once. Equal primary scores share `combined_score_rank`; lexical pair ID provides only a deterministic display order in `combined_rank`.

## Limitation

The legacy and V2 libraries were assembled under different candidate-selection rules. The resulting percentile is descriptive of this frozen expanded universe, not a population probability, false-discovery rate, or uniformly sampled biological background.

Descriptive same-register pMHC sequence prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
