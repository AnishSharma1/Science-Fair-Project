# Control-calibrated HLA-specific rankings

- Ranking scope: **within each HLA only**
- Cross-allele consensus: **not used**
- Formal control reference: **seed 104759**
- Gold-standard available-set capture@1: **2/2 seeds**
- Independent gold-standard systems: **1**

Open `RESULTS_SUMMARY.md` for four separate top-10 tables. Full rankings are under `rankings/`; `top_25_by_hla.csv` is a compact combined view that preserves the HLA boundary.

The earlier full-ensemble package is preserved unchanged. This additive version reuses its frozen pair geometries and control results, changes the reporting and calibration layer, and does not rerank across HLA alleles.

Reproduce:

```bash
PYTHONPATH=src python3 src/build_control_calibrated_hla_rankings.py
```

Computational pMHC geometry prioritization only; control-reference metrics are not probabilities and do not establish presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism.
