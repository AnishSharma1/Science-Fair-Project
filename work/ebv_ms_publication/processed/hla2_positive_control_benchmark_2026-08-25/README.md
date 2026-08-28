# Held-out human HLA-II positive-control benchmark

This additive package expands the strict denominator from one to **3 independent human TCR systems** and **4 required positive comparisons**. The original discovery rankings are unchanged.

## Current status

- Overall trust status: **not_evaluable**.
- Discovery reranking allowed: **false**.
- New AlphaFold jobs: **48**, split as **30 + 18**.
- AlphaFold state: **prepared, not submitted**.
- The two original missing Hy.2E11 jobs were not retried, replaced, or silently excluded.

## Experimental positive geometry

- PAIR_HY2E11_BALF5_MBP: exposed CA 0.908 A; side-chain vector 1.624 A
- PAIR_OB1A12_ENGA_MBP: exposed CA 0.426 A; side-chain vector 1.534 A
- PAIR_HY1B11_UL15_MBP: exposed CA 0.313 A; side-chain vector 0.967 A
- PAIR_HY1B11_PMM_MBP: exposed CA 0.430 A; side-chain vector 0.405 A

## PDB oracle ranks

- PAIR_HY2E11_BALF5_MBP: 7 exact-HLA decoys; available exposed-CA rank 3; complete
- PAIR_OB1A12_ENGA_MBP: 5 exact-HLA decoys; available exposed-CA rank 1; complete
- PAIR_HY1B11_UL15_MBP: 2 exact-HLA decoys; available exposed-CA rank 1; not_evaluable_insufficient_exact_hla_decoys
- PAIR_HY1B11_PMM_MBP: 2 exact-HLA decoys; available exposed-CA rank 3; not_evaluable_insufficient_exact_hla_decoys

Hy.2E11 and Ob.1A12 have evaluable PDB oracle panels under the locked five-decoy minimum. Hy.1B11 remains PDB-not-evaluable because only two exact-HLA structural pair decoys are available for each positive. The new Ob.1A12 and Hy.1B11 AlphaFold panels are pending submission and download. These missing layers intentionally block a formal pass.

Computational pMHC geometry prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS mechanism.
