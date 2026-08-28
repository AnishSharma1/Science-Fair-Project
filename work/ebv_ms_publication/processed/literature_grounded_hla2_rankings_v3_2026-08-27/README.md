# Literature-Grounded HLA-II Ranking V3

This additive package ranks 6,400 V2 pairs separately within four HLA-DR alleles and ranks the frozen 2,043-pair DRB1*15:01 combined universe. Existing rankings and benchmark outputs were not modified.

The primary key is TCR-facing P2/P3/P5/P7/P8 BLOSUM62. Physicochemical mismatch, identity, local modeled surface, and lexical pair ID break ties in that order. Local surface is a separate model-derived annotation summarized conservatively at the 75th percentile across model combinations. Binding percentile is reported separately and never enters the rank.

Evidence tier `M` means missing structure or failure of the strict register-robustness rule; it does not remove a pair from the primary sequence rank. The existing three systems are development controls only. `definitive_validation_gate.json` therefore remains not evaluable and discovery unlock is false.

Descriptive HLA-specific pMHC sequence prioritization with modeled local-surface annotations only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
