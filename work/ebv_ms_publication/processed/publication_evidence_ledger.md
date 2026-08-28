# EBV--MBP structural hypothesis: evidence ledger

## Claim that can be made now

The EBV BALF5 peptide `TGGVYHFVKKHVHES` and MBP(85--99)
`ENPVVHFFKNIVTPR` are an **established molecular-mimicry positive-control
system**. The experimentally verified mechanism is allele-pair switching
within the DR15 haplotype: Hy.2E11 recognizes MBP on HLA-DRB1*15:01 (DR2b)
and BALF5 on HLA-DRB5*01:01 (DR2a). This project may use that system to
calibrate a computational screening workflow; it must not claim discovery of
the BALF5--MBP relationship.

## Experimental anchors

| Item | Status | Proper use |
|---|---|---|
| PDB 1YMM: Ob.1A12 TCR + DRB1*15:01 + MBP | Experimental crystal structure | Positive control and TCR orientation reference |
| PDB 2WBJ: Ob.1A12 TCR + DRB1*15:01 + microbial peptide | Experimental crystal structure | Demonstrates the receptor can recognize a non-self ligand; not EBV evidence |
| PDB 1H15: DRB5*01:01 + EBV BALF5(627--641) | Experimental crystal structure | Correct EBV pMHC positive-control arm |
| PDB 1BX2: DRB1*15:01 + MBP(85--99) | Experimental crystal structure | Correct MBP pMHC positive-control arm |
| IEDB assay 1686220, EBV BALF5 627--641 | Curated positive T-cell assay; DRB1*15:01 | Candidate provenance, not a paired-TCR sequence |

## Hy.2E11 receptor provenance

Primary-source chain annotations identify Hy.2E11 as V-alpha 3.1 with junction
`TDSGGSYIPTFGRGTSLIVHP` and V-beta 4/4.3 with junction
`PSGQGTYGYTFGSGTRLTVV`. They do not provide a complete allele-resolved
alpha/beta sequence. Consequently, no Hy.2E11 ternary input is generated;
see `experimental_positive_control/HY2E11_SEQUENCE_AUDIT.md`.

## Computational evidence retained

- All 86 rank-1 pMHC structures passed expected-peptide and chain-layout QA.
- The corrected full screen ranks the BALF5--MBP family first.  For
  `EBV_TCELL_63843` versus `HUMAN_MYELIN_114806`, the six-residue local
  alignment has property similarity 0.8431 and a groove-normalized peptide CA
  RMSD of 0.381 A.
- In the experimental 1H15 versus 1BX2 structures, the homologous
  seven-position peptide core has a 0.838 A CA RMSD after fitting the HLA
  peptide-binding platforms. This is a geometric positive-control metric only.
- Fixed-template compatibility screen: the EBV peptide produces 29
  Ob.1A12--peptide heavy-atom contacts with no TCR clashes; the exact MBP
  positive control has 36 contacts and no TCR clashes.  The EBV scaffold has
  six peptide--HLA overlaps, so it is not a final docked model.

## Dated lead-focused structural-control audit (2026-08-15)

The frozen discovery ranking was retained and the two leads were analyzed in
separate layers; they were not pooled, averaged, or assigned equal evidentiary
weight. See
`lead_focused_robustness_2026-08-15/LEAD_FOCUSED_FINDINGS.md`,
`lead_focused_robustness_2026-08-15/control_rank_and_leave_one_out.csv`,
`lead_focused_robustness_2026-08-15/technical_bootstrap_summary.csv`,
`lead_focused_robustness_2026-08-15/job_pair_stability.csv`, and
`lead_focused_robustness_2026-08-15/pose_cluster_membership.csv`.

- **Rank 1, primary strict-control lead:**
  `EBV_TCELL_950::HUMAN_MYELIN_112214` is the only strict-control lead
  (`consistent_positive`). Its target median was 0.643 A versus an
  equal-weight strict-background median of 7.964 A (background-minus-target
  delta 7.321 A); leave-one-control-out deltas ranged from 5.261 to 10.086 A.
  The exploratory target rank was 1 of 4 (empirical tail fraction 0.25), and
  the 10,000-iteration technical-stability interval was 1.620 to 12.759 A.
  Pose/job consistency: all four target job-pair medians were below the
  equal-weight background median (0.427--0.864 A), and the largest pose
  cluster contained 23 models from 7 jobs.
- **Rank 2, supplemental length sensitivity:**
  `EBV_TCELL_2268741::HUMAN_MYELIN_117032` is
  `length_sensitivity_only__mixed_positive`, not a second primary lead. Its
  deliberate length-sensitivity layer had a target median of 6.396 A, an
  equal-weight background median of 13.975 A, and a delta of 7.578 A;
  leave-one-control-out deltas ranged from 5.303 to 9.944 A. Its exploratory
  tail fraction was 0.25, but target job-pair medians ranged from 0.465 to
  12.328 A and its technical-stability interval crossed zero (-0.517 to
  16.880 A), so it remains sensitivity-only.

These empirical tail fractions are exploratory ranks, not p-values. The
bootstrap intervals quantify technical stability across saved AlphaFold
jobs/models only, not biological replication. This audit is descriptive pMHC
geometry and does not establish peptide presentation, TCR binding, activation,
cross-reactivity, molecular mimicry, or an MS mechanism.

## Claims that must not be made

- A ColabFold pMHC or ternary model does not demonstrate TCR binding.
- It does not establish affinity, activation, molecular mimicry in patients,
  or a causal role in multiple sclerosis.
- Ob.1A12 is an MBP/DRB1*15:01 positive-control receptor. It is not an
  experimentally verified receptor for this EBV peptide.
- The established Hy.2E11 cross-reaction does not present both ligands on
  DRB1*15:01; replacing DRB5*01:01 with DRB1*15:01 changes the biological
  question.

## Ternary ColabFold calibration result

The three-condition, five-rank ternary screen was completed. Although it
recovered the HLA scaffold, it did not recover the experimental Ob.1A12 TCR
orientation for the exact MBP positive control. Therefore its EBV versus
negative-control contact differences are not valid evidence for selective
cross-recognition. See `ob1a12_ternary_evaluation/FINAL_CONCLUSION.md`.

## Decisive next evidence

1. Obtain a primary-source sequence for the Hy.2E11 alpha and beta chains;
   do not infer or assemble them from partial receptor annotations.
2. Use Hy.2E11 with the experimental pMHC pairs as a calibration benchmark.
3. Define performance before screening new candidates: the method must rank
   the known positive above matched negative pMHC pairs.
4. Only then screen additional EBV--autoantigen candidates, reporting hits as
   computational priorities that require activation or binding experiments.
