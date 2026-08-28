# EBV-MS T-cell Library V2 model analysis (2026-08-25)

- Discovery downloads: **319/320**
- Geometry-evaluable discovery jobs: **319**
- Geometry-complete discovery pairs: **6360/6,400**
- Calibration downloads: **22/24**
- Geometry-complete calibration comparisons: **61/72**
- Strict positive recovery: **not_evaluable_incomplete_calibration**
- Gold-standard available-set capture@1: **2/2 seeds**
- Gold-standard independent systems: **1**

The analysis uses all valid AlphaFold model combinations and preserves the one missing discovery job and two missing calibration jobs as technical missingness. Original Downloads and the frozen V2 package are read-only inputs.

Distinct unplanned duplicate runs are excluded from primary ranking and reported in `supplemental/` as technical-repeat sensitivity.

The locked positive-control audit is in `validation/`. It verifies exact biological identity and experimental structures before reading ranks, and it records that the model and score were not changed to fit the positive.

Reproduce:

```bash
PYTHONPATH=src python3 src/run_tcell_library_v2_model_analysis.py
```

Computational pMHC geometry only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism.
