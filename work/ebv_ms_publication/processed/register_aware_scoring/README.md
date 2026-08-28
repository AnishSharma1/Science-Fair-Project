# Register-aware scoring outputs

`register_aware_pair_scores.csv` contains every frozen PASS-universe pair.
Sequence and physicochemical descriptors are recalculated only for existing
local-alignment coordinates that occupy the same P1-P9 register index in both
resolved peptide cores. P1/P4/P6/P9 are labelled HLA-anchor positions and
P2/P3/P5/P7/P8 candidate-exposed positions; these labels do not assert TCR
contacts.

Rows with one or two retained coordinates are `limited_coverage_report_only`.
Only rows with at least three retained coordinates are eligible for the robust
primary sequence/chemistry ranking. The score is a descriptor, not a molecular
mimicry probability or evidence of presentation, shared-TCR binding,
activation, cross-reactivity, or MS mechanism.

`whole_local_alignment_geometry_rmsd_context` is included strictly as context:
it was precomputed across the original local alignment. It is not a
same-register structural score and cannot enter the primary score until the
original per-residue pMHC coordinates are recovered and re-analysed.

A positive computational prioritization result requires an independently
interpretable target, complete frozen strict decoys, robust score coverage,
and persistence under retained register-sensitivity analysis. Otherwise the
result is negative/mixed: the current evidence does not support robust
register-aware prioritization. This is not evidence against all EBV-myelin biology.

Regenerate:

```bash
PYTHONPATH=src python3 src/build_register_aware_score_table.py
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v
```
