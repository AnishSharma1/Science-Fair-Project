# Initial rebuilt analysis result

The composition-preserving 10,000-permutation sensitivity screen found **zero
exact contiguous 5-mer overlaps** between the 12 unique positive EBV
HLA-DRB1*15:01 peptides and the 53 coordinate-validated human myelin peptides.
The observed mean maximum 5-mer-overlap fraction was 0.0; the empirical
one-sided p-value under the residue-shuffle null was 1.0.

This is useful as a guardrail: the current data do not support a claim of
strong exact local sequence identity. It does **not** rule out biochemical
similarity, MHC-binding similarity, TCR cross-reactivity, or molecular mimicry;
those require different models and experimental evidence.

The result is stored in
`processed/null_model_5mer_results.json`. It is an exploratory checkpoint, not
the final statistical analysis, because the peptide set is small and the human
workbook is an aggregate epitope inventory.
