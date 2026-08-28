# Multi-allele EBV–MS pMHC analysis: manuscript-ready working text

## Methods

A fixed panel of 25 EBV and 25 CNS/self peptides was modeled independently with HLA-DRB1*13:03, HLA-DRB1*03:01, and HLA-DRB1*08:01. The prespecified primary endpoint was transfer of the EBV_TCELL_950–HUMAN_MYELIN_112214 EBNA1–MBP pair; all other 625 within-allele combinations were exploratory. AlphaFold Server outputs were accepted only when all five CIF, confidence-summary, and full-data files were present and exact DRA, DRB, and peptide sequences agreed with the request. The highest-ranked clash-free sample represented each job for descriptive QC; geometry used all valid clash-free sample combinations without post-hoc confidence thresholds.

IEDB recommended-binding predictions were generated separately for each allele. Exact seq_num values, raw responses, percentile ranks, predicted cores, and binding-rank bins were retained. The two 10-residue EBV peptides used source-verified natural flanks, but a structural register was eligible only when the predicted P1–P9 core was unique and fully contained in the modeled peptide.

Structures were superposed within allele using equivalent HLA-DRA and HLA-DRB groove Cα atoms. P1–P9 geometry was summarized for the full core, anchor positions P1/P4/P6/P9, and candidate TCR-facing positions P2/P3/P5/P7/P8. The primary metric was the median exposed-position Cα RMSD across valid sample combinations. No raw metric was pooled across alleles.

Three non-CNS human controls per evaluable allele were chosen before inspecting any AlphaFold geometry. Controls matched MBP length within one residue and the allele-specific binding-rank bin and were ordered by amino-acid composition distance, peptide length, and numeric IEDB identifier. Robustness inputs use fixed seeds 104729 and 104759.

## Current results

The current download set contains 143/150 canonical complete jobs and 1 duplicate sensitivity run. Missing jobs were not imputed. Exact sequence QA was performed for 143 jobs and 715 samples.

- HLA-DRB1*13:03: 50/50 modeled peptides had a unique fully contained predicted core.
- HLA-DRB1*03:01: 49/50 modeled peptides had a unique fully contained predicted core.
- HLA-DRB1*08:01: 50/50 modeled peptides had a unique fully contained predicted core.

Primary anchor geometry from the discovery-seed matrix:

- HLA-DRB1*13:03: 0.864191 Å; primary robustness status is pending fixed-seed controls.
- HLA-DRB1*03:01: 1.175893 Å; primary robustness status is pending fixed-seed controls.
- HLA-DRB1*08:01: 1.622352 Å; primary robustness status is pending fixed-seed controls.

These results describe computational pMHC geometry only. They do not establish natural presentation, TCR binding, T-cell activation, cross-reactivity, molecular mimicry, or an MS disease mechanism. The primary anchor-versus-control result remains incomplete until the fixed-seed robustness jobs are modeled and analyzed.
