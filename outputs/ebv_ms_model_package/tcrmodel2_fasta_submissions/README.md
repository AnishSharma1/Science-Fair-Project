# TCRmodel2 submission bundle

Use the **TCR-pMHCII** tab at [TCRmodel2](https://tcrmodel.ibbr.umd.edu/). Upload one FASTA file per job. The five records are deliberately ordered: TCR alpha, TCR beta, peptide, MHC alpha, MHC beta. Headers are descriptive only; the server uses record order.

Although the current server help describes accepting class-II peptides of at least 9 aa, the failed submissions show that this path is enforcing 11 aa. Every pending FASTA therefore uses exactly 11 residues: a 9-aa core with one N- and one C-terminal flank. The already successful 01_CAL_1YMM file is left unchanged.

## Submit first

1. 02_CAL_1ZGL.fasta — exclude 1ZGL.
2. 03b_CAL_2WBJ_CANONICAL_OB.fasta — revised 2WBJ calibration using canonical OB.1A12 variable domains; exclude 2WBJ.
3. 04_HY_MBP_DRB1.fasta — exclude 1YMM, 1ZGL, 2WBJ, 1H15.

Only interpret the Hy.2E11 jobs if the template-excluded calibrators produce credible structures. No resulting score proves binding or cross-reactivity.

Download all five models and result JSON for each submitted job. Record model confidence, TCR-pMHC ipTM, I-pLDDT, template list, and the unusual-docking warning, if any.
