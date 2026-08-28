# Combined DRB1*15:01 same-register ranking

This additive package unifies the V2 DRB1*15:01 screen with every eligible pair from the earlier 636-pair register-aware universe. It does not average old and new rank numbers. Every admitted pair is rescored using the control-selected TCR-facing BLOSUM62 method, exact peptide-pair duplicates are collapsed, and one new rank is assigned across the combined universe.

- `combined_ranked_pairs.csv`: complete combined ranking
- `top_25_combined.csv`: compact review table
- `legacy_eligibility_audit_636.csv`: inclusion or exclusion reason for all old pairs
- `exact_duplicate_crosswalk.csv`: deduplicated old/new pair identities
- `ranking_basis.json`: machine-readable ranking rule

An asterisk marks a computationally prioritized pair without exact experimental paired-recognition evidence. The result remains provisional because the selected endpoint has three strict control systems, below the six-system definitive gate.

The legacy and V2 candidate-selection rules were not identical. Therefore, the combined rank is an expanded-library prioritization, and its percentile describes only this frozen 2,043-pair universe rather than a uniformly sampled biological background.

Descriptive same-register pMHC sequence prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
