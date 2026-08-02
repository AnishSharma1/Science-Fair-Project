# Pre-specified similarity analysis plan (version 0.1; 2026-08-01)

This plan is written before interpreting the exploratory ranking.

## Primary comparison

- EBV set: unique HLA-DRB1*15:01 EBV T-cell peptides with a positive IEDB
  qualitative outcome.
- Human set: unique myelin-candidate peptides from the human workbook whose
  IEDB source mapping has validated coordinates.
- Pairwise metric: the maximum local sequence identity across alignments with
  at least five aligned residues, using the documented score +2 for a match,
  −1 for a mismatch, and −2 for a gap.
- Summary statistic: the mean of the EBV peptide-level maxima.

## Null model

For each of 10,000 permutations, independently shuffle the residues within
each EBV peptide. This preserves peptide length and amino-acid composition but
destroys residue order. Recompute the same peptide-level maximum and mean.
The one-sided empirical p-value is:

`(1 + number of null means >= observed mean) / (1 + 10,000)`

For this push, the executable null-model check uses a faster, pre-specified
sensitivity metric: the maximum fraction of contiguous 5-mers in each EBV
peptide shared with any human myelin peptide. This is a computational check on
the same composition-preserving null, not a replacement for the alignment
analysis planned for the final manuscript.

This is a sequence-resemblance null, not a biological model of MHC binding,
antigen processing, TCR recognition, or MS risk.

## Guardrails

- The threshold of five aligned residues is fixed before analysis.
- No structural score, simulated value, or neural-network output enters this
  test.
- A small p-value would indicate unusual sequence resemblance under this null,
  not molecular mimicry or disease causation.
- The final manuscript should report sensitivity to alignment scoring and to
  excluding repeated/near-duplicate human epitopes.
