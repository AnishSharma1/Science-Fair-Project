# Multi-feature sequence-structure concordance rankings

This additive package places five interpretable measurements beside every pair and ranks only pairs with complete same-register sequence and comparable structural data. Each HLA remains separate. The expanded DRB1*15:01 universe is also ranked separately.

The three sequence metrics are converted to within-universe percentiles and averaged. The two structural RMSDs are converted the same way and averaged. The primary rank minimizes the worse family percentile, so a pair must be reasonably strong in both sequence and structure. The balanced family mean is the secondary score and lexical pair ID is the deterministic final display tie-break.

"Whole-protein RMSD" is not available because the AlphaFold models contain pMHC complexes, not the complete EBV and human source proteins. The scientifically relevant broad structural endpoint here is full P1-P9 peptide-core C-alpha RMSD after aligning equivalent HLA-groove atoms. Whole-pMHC RMSD is not used because the shared HLA scaffold would dominate it.

The five features describe complementary proxies of the displayed pMHC surface. They are not a mechanistic model of what a TCR "considers," and no TCR sequence, TCR contact map, binding affinity, or activation assay enters this ranking. The control audit is retrospective and exploratory; it cannot freeze weights or unlock discovery.

Reproduce with:

```bash
PYTHONPATH=src python3 src/build_multifeature_concordance_rankings.py
```

Retrospective exploratory audit on the three existing positive-control systems; it does not freeze weights, satisfy the six-system definitive gate, or unlock discovery.
Descriptive same-register pMHC sequence-structure concordance prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
