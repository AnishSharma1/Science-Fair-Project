# External validation rank-recovery benchmark

This analysis overlays pre-annotated literature/context records onto the predeclared pMHC shortlist.

Counts are annotated pair records, not independent positive examples: the BALF5--MBP family contains overlapping records and the newer source annotations do not establish direct EBV--myelin cross-reactive pairs.

It is a prioritization benchmark only. It does not test TCR binding, T-cell activation, affinity, or patient pathogenicity.

- classic_BALF5_MBP_pair_recovery: 10 annotated pair records in a 32-pair universe; mean rank 7.1; top-10 annotated records 8; top-10 empirical p=0.0001
- any_external_overlay: 16 annotated pair records in a 32-pair universe; mean rank 11; top-10 annotated records 10; top-10 empirical p=0.0003
- strict_new_literature_overlay: 5 annotated pair records in a 32-pair universe; mean rank 18.2; top-10 annotated records 2; top-10 empirical p=0.5036
