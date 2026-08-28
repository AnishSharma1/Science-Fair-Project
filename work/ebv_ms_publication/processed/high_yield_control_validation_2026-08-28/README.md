# High-Yield Control Validation

This additive package evaluates 12 frozen discovery candidates against 25
score-blind, exact-HLA N3 comparison pairs each. It reuses the existing V3
model-derived features and recomputes all panel-level percentiles and ranks.

## Result

- Complete panels: 12/12
- Total rows: 312 (12 targets and 300 N3 comparisons)
- Targets ranked in the top three: 8/12
- New untouched strict positive-control systems admitted: 0
- Specificity gate: not evaluable; no explicit N1/N2 registry
- Discovery unlock: false

N3 means recognition is unknown. These rows are fair computational ranking
comparators, not biological negative controls. Sequence-lane target structure
is retained only as register-sensitivity diagnostics and abstains from the V3
primary structural tie-break.

The originally proposed global exclusion of every frozen target arm was not
feasible: fewer than five unique eligible self arms remained. The frozen
resolution excludes the current panel's target arms and all confirmed-control
ligands, while explicitly flagging N3 rows that contain an arm from another
high-yield target. The feasibility audit uses no geometry.

## Primary-source census

The existing Hy.2E11, Ob.1A12, and Hy.1B11 systems remain development
controls. DQ8 A2.13/ET650-4 systems are retained as structure-resolved leads
but fail the frozen distinct-source rule. Exact-allele searches for the four
studied DR alleles produced no newly admissible system with all strict fields.

N3 panels provide descriptive, HLA-specific computational ranking context only; they are not evidence of presentation, TCR recognition, activation, specificity, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
