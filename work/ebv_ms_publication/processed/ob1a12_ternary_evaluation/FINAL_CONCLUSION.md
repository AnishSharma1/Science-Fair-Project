# Final conclusion: Ob.1A12 ternary ColabFold test

## Decision

**Inconclusive for EBV--MBP cross-reactivity. Do not report a positive docking
result from this ternary ColabFold run.**

This is a calibration failure of the computational assay, not evidence against
the biological hypothesis.

## What passed

- In all fifteen models, the HLA component fit the experimental 1YMM HLA
  scaffold closely (median HLA CA RMSD 0.577--0.589 A).
- The TCR chains themselves were confidently folded (median TCR pLDDT
  85.85--87.30).
- The actual ColabFold global rank-1 scores were strong and nearly matched:
  MBP positive control pLDDT/pTM/ipTM = 90.5/0.808/0.778; EBV hypothesis =
  90.2/0.806/0.775; EBV negative control = 89.3/0.785/0.754.

## What failed

The exact positive control is known experimentally: Ob.1A12 bound to
DRB1*15:01--MBP in PDB 1YMM. Yet its predicted TCR pose did not recover the
experimental orientation (median HLA-aligned TCR alpha RMSD 10.227 A; beta
RMSD 18.040 A). The five ranks also disagree markedly on the positive-control
peptide pose; only rank 1 recovers it closely (0.535 A).

Because the model cannot reproduce the known receptor--pMHC orientation, its
EBV contact counts and global ipTM cannot be used as a selective binding
comparison. In a five-chain complex, ipTM covers multiple interfaces, including
the stable HLA alpha--beta interface; it is not a peptide-specific TCR binding
score. Rank 1 does not rescue the argument: EBV has 38 peptide-contact atoms
versus 45 for the positive control and 39 for the negative control.

## Permitted statement

"The EBV BALF5 and MBP peptide pair remains a structurally and
evidence-prioritized molecular-mimicry hypothesis, but this ColabFold ternary
screen was not calibrated to infer Ob.1A12 cross-recognition."

## Best next step

Do not spend more ColabFold credits on the same unconstrained five-chain
prediction. A new structural claim requires either restrained/template-guided
docking anchored to 1YMM with independent scoring, or an experimental assay
(pMHC tetramer, activation assay, or binding measurement).
