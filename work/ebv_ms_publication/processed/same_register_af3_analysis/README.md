# Direct same-register AF3 analysis

This analysis compares the frozen 32-pair shortlist only. A pair enters the
table when both pMHCs have complete AF3 samples and each peptide has a unique,
fully contained primary-allele computational P1--P9 core hypothesis.

Unlike the earlier local-alignment diagnostic, this endpoint compares all nine
equivalent register positions directly. HLA-DRA/DRB groove C-alpha atoms are
fitted first; P1--P9, anchor-position (P1/P4/P6/P9), and candidate-exposed-
position (P2/P3/P5/P7/P8) C-alpha RMSDs are then reported across every
available AF3 sample combination.

- Eligible predeclared pairs: **9**
- AF3 cross-sample geometry comparisons: **600**
- Pairs with at least one completed matched background comparator: **0**

The background rule was frozen before structural inspection: human peptide
length within one residue and the same IEDB predicted-binding-rank bin. No
completed matched comparator exists for any eligible pair, so no controlled
structural effect size, enrichment test, or p-value is computed. Where frozen
background cores exist despite failed structure generation, the sequence-only
table reports a descriptive same-register chemistry rank; these small control
sets do not support inferential statistics.

AlphaFold samples and seed repeats are technical sensitivity analyses, not
independent biological replicates. These files do not establish peptide
presentation, a shared TCR surface, T-cell activation, cross-reactivity,
molecular mimicry, or an MS mechanism.
