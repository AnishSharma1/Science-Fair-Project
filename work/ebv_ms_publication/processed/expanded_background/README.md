# Expanded human-background comparator inputs

- Source-table records scanned: **102**
- Human-background records reviewed: **37**
- Direct IEDB MHC-II and pMHC batch inputs: **27**
- Retained but not modeled for missing verified natural flanks: **8**

The selected records are a pre-score human-background comparator arm. They
were chosen from source/provenance and direct IEDB MHC-II eligibility only,
before binding predictions, structure results, register scores, or candidate
priority values were inspected. Full IEDB peptide sequences are preserved.

This batch can expand strict-decoy feasibility only after all three gates are
met: matching IEDB binding-rank bins, same-register eligibility, and passed
pMHC structure QA. It does not create a biological negative control or
evidence of molecular mimicry.

Reproduce:

```bash
PYTHONPATH=src python3 src/build_expanded_background_inputs.py
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v
```
