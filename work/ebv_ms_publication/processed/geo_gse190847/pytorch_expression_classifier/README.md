# PyTorch expression classifier for GSE190847

This analysis uses a seven-gene HLA-II/APC expression panel from RMA log2 microarray B-cell data.

It is intentionally bounded: it tests out-of-sample group signal in peripheral B cells and does not measure EBV infection, pMHC presentation, TCR binding, T-cell activation, or MS causality.

- LOOCV AUC: 0.662
- Balanced accuracy at 0.5: 0.685
- Empirical permutation p(AUC >= observed): 0.089 using 100 permutations
- PyTorch version: 2.8.0
