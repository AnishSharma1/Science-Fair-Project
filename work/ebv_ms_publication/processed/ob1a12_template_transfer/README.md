# Ob.1A12 template-transfer screen

This folder tests a narrowly defined structural hypothesis: whether the
experimentally solved MBP-reactive human TCR **Ob.1A12** can be used as an
orientation template for completed DRB1*15:01 pMHC models.

## Experimental anchor

- **PDB 1YMM:** Ob.1A12 TCR + HLA-DRA*01:01/DRB1*15:01 + MBP peptide
  `ENPVVHFFKNIVTPRGGSGGGGG` (the first 14 residues are the project MBP core).
- TCR alpha CDR3: `CATDTTSGTYKYIF`; TCR beta CDR3: `CSARDLTSGANNEQFF`.
- **PDB 2WBJ:** the same Ob.1A12 TCR and DRB1*15:01 with a microbial peptide,
  providing an independent example of a non-self ligand for this receptor.

## Conditions

1. `positive_mbp`: exact MBP core positive control.
2. `test_ebv_mimic`: EBV BALF5 peptide `TGGVYHFVKKHVHES`.
3. `negative_ebv_control`: unrelated DRB1*15:01 EBV peptide
   `VTNILIYNGWYADS`.

## How to interpret the CSV

The script rigidly transfers the crystal TCR after a sequence-based fit of both
HLA chains. It does **not** perform flexible docking, energy minimization,
binding-affinity calculation, or a cross-reactivity prediction.

Only `eligible_as_initial_tcr_geometry = True` (combined HLA CA RMSD <= 6 A)
may be carried forward as an initial docking geometry. A false value means the
completed ColabFold pMHC is structurally incompatible with this particular
template transfer and must be re-modeled using an HLA template before TCR work.
