# pMHC modeling readiness gate

The next push can begin structural modeling only for rows in
`processed/pmhc_candidate_manifest.csv` that pass every gate below.

1. Exact HLA-DRB1*15:01 and MHC-II class are recorded.
2. Peptide sequence contains only standard amino acids and retains its IEDB ID.
3. EBV rows retain assay-level evidence; human rows retain a validated source
   accession and coordinates.
4. The peptide is represented as a class-II binding register, not as an
   unexamined 9-mer class-I-style input.
5. The MHC-II alpha/beta chains, peptide, and residue numbering are recorded.
   The two protein chains must be independently identified as DRA-like and
   DRB-like; two DRB-like chains are an invalid class-II layout.
6. Any predicted structure is labeled predicted and is never mixed with an
   experimental PDB without a separate indicator.
7. pLDDT/PAE and model settings are preserved for every run.
8. TCR docking remains a separate, later analysis; no generic TCR score is a
   primary pMHC result.

Legacy PDBs with filename-only mappings remain quarantined until these gates
are satisfied.

The preferred experimental geometry reference is PDB 8TBP, an
HLA-DRB1*15:01 complex with separate alpha, beta, and peptide chains.
