# Multi-allele EBV–MS analysis package (2026-08-21)

- Canonical AF3 jobs complete: **143/150**
- Sample QC rows currently available: **715/750**
- Allele-specific panel predictions: **150/150**
- Pair-universe rows: **1875/1,875**
- Geometry sample combinations: **41,875**
- Frozen controls: **9**
- Fixed-seed robustness jobs prepared: **30/30 maximum**

The package is reproducible from saved inputs and raw IEDB responses. Downloads are referenced read-only; all derived files live here. The 20-job retry JSON is copied into this package when available, but a prepared/uploaded JSON is not a completed AlphaFold run.

Computational pMHC geometry only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS mechanism.

Reproduce:

```bash
PYTHONPATH=src python3 src/run_multiallele_manuscript_analysis.py
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py'
```
