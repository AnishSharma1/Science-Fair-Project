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
