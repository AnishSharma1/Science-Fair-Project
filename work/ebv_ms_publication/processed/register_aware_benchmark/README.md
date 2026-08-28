# Register-aware matched-decoy benchmark

## What this artifact does

This is a complete-universe feasibility and pMHC candidate-prioritization
artifact. It does not establish shared-TCR binding, peptide presentation,
T-cell activation, cross-reactivity, or an MS mechanism.

## Frozen computational boundary

- Register hierarchy: `processed/register_sensitivity/experimental_register_overrides.csv`
  is applied before the stored IEDB top-core hypotheses.
- PDB 1BX2 makes `VHFFKNIVT` at positions 5--13 the exact MBP(85--99)
  DRB1*15:01 reference register. PDB 1H15 makes BALF5 `YHFVKKHVH` a
  DRB5*01:01 calibration-only register, not a primary-screen override.
- Exact gH `EKQLFYYIGTMLPN` and candidate-MBP `QRPGFGYGGRASDYKSAHK`
  records are sensitivity-only because no exact DRB1*15:01 experimental
  register was established.
- All candidates without a registry decision retain `IEDB recommended_binding HLA-DRB1*15:01 top-core hypothesis` as a
  computational hypothesis.
- A pair is assessable only when both resolved cores are unique and contained
  in the manifest peptide, neither arm is calibration-only or sensitivity-only,
  and at least one pre-existing local alignment maps to the same P1--P9 position.
- Decoy matching: Both EBV and human peptide lengths within one residue and zero IEDB predicted-binding-rank-bin mismatches; composition distance and peptide pLDDT order only already-eligible decoys.
- Selection boundary: Decoys were selected from assessable background rows using only peptide length, amino-acid composition, peptide pLDDT, and IEDB binding-rank bins; never a pMHC priority score, rank, or structural similarity value.

## Results counts

- pass_pair_universe: 636
- assessable_pairs: 33
- eligible_background_pairs: 23
- annotated_target_records: 349
- ready_targets: 1
- partial_targets: 5
- no_eligible_decoy_targets: 4
- not_assessable_targets: 339

Annotated target records are source/context overlays, not independent
biological positives. In particular, BALF5--MBP overlapping records remain one
DR15-haplotype calibration system, not replicated validation systems.

## Interpretation

The output measures whether the current computational inputs can support the
predeclared strict comparison. It does not calculate a new pMHC effect size or
rerank candidates. A target labelled `ready` has five covariate-matched,
score-blind decoys available; this is a design-readiness property, not support
for molecular mimicry.

## Reproduce

```bash
PYTHONPATH=src python3 src/build_register_aware_benchmark.py
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v
```

## Remaining primary-evidence gates

1. Establish allele-resolved experimental registers for the exact gH and
   candidate-MBP peptides before they enter a same-register benchmark.
2. Obtain direct DRB1*15:01 evidence for BALF5 before promoting it beyond
   DR15-haplotype calibration.
3. Use exact pMHC structures, not the register table alone, to assign any
   peptide-specific solvent exposure or candidate TCR-facing positions.
