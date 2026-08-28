# Methods

## Input freeze and identity

The frozen V2 package defines 320 discovery jobs (80 peptides across four HLA-DRB1 alleles), 6,400 within-allele EBV-self pairs, 24 fixed-seed native-HLA calibration jobs, and 72 calibration comparisons. Downloaded folders were matched by the request name inside `job_request.json`. A complete bundle required five CIF files, five confidence summaries, five full-data JSON files, and one request file. Duplicate copies were selected by the lexicographically smallest content fingerprint before reading model scores.

## Locked gold-standard denominator

Gold-standard eligibility was determined without reading model scores. It required exact peptide identities and HLA arms, recognition of both pMHCs by the same human T-cell clone, and experimental pMHC structures for both arms. The current denominator contains one independent system, Hy.2E11 DRB5*01:01-BALF5 versus DRB1*15:01-MBP. Supportive T-cell studies, antibody-only pairs, protein-level associations, canonical tiles, and computational discoveries were excluded. Recovery was reported at ranks 1 and 3 for each available seed; incomplete control sets were not relabeled as formal 26-pair tests.

## Model QC

Every CIF was required to contain exactly HLA-DRA chain A, the expected HLA-DRB chain B, and the exact peptide chain C. All chain sequences had to match the frozen request. Samples with an AlphaFold clash flag were excluded from geometry but retained in QC. No pLDDT, ipTM, PAE, contact-probability, or geometry threshold was used to select samples. The highest-ranking clash-free exact sample is reported only as a descriptive representative; full-ensemble analyses use every valid sample.

## Register-aware geometry

IEDB `recommended_binding` top cores and exact `seq_num` mappings were frozen before structural analysis. Equivalent HLA-DRA and HLA-DRB groove C-alpha atoms from the first 85 residues of each chain were superposed by the Kabsch algorithm. The primary endpoint was median exposed-position P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove fit. Full-core, anchor-position P1/P4/P6/P9, pseudo-C-beta, and exposed-side-chain-centroid RMSDs are supporting descriptors. Every valid 5-by-5 model combination was calculated; medians and interquartile ranges summarize technical model sensitivity and are not inferential confidence intervals.

## Ranking and calibration

Discovery pairs were ranked separately within each allele by median exposed-position C-alpha RMSD, then ensemble IQR, then frozen pair ID. Cross-allele consensus minimizes the worst allele-specific percentile, then the median percentile, and requires all four allele results; raw geometry is never pooled across alleles. The BALF5-MBP calibration remains frozen as rank <=3 of 26 on both seeds plus positive RMSD below the equal-weight median of 25 full decoys. Missing calibration models are not imputed, so an incomplete seed is reported as not evaluable.

## Unplanned duplicate sensitivity

Complete duplicate folders with identical content were treated as copies. Distinct content fingerprints under the same frozen job identity were excluded from primary selection and analyzed separately as unplanned technical-repeat sensitivity. They cannot increase the primary sample size or serve as biological replication.

## Interpretation boundary

Computational pMHC geometry only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism.
