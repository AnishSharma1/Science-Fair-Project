# Experimental pMHC-II template validation

The official RCSB search identified two directly useful HLA-DR15 templates:

| PDB | Structure | Use in this project |
|---|---|---|
| [8TBP](https://www.rcsb.org/structure/8TBP) | HLA-DRB1*15:01 with a separate 15-residue Smith-antigen peptide chain; X-ray structure | Primary template candidate |
| [5V4M](https://www.rcsb.org/structure/5V4M) | HLA-DR15 with alpha3(135–145); X-ray structure, but the peptide is part of a chimeric beta-chain construct | Secondary geometric/reference template |

8TBP is the preferred template because its metadata explicitly names
HLA-DRB1*15:01 and provides separate alpha, beta, and peptide chains. One
N-terminal peptide residue is not resolved in the coordinate ATOM records even
though the polymer entity sequence contains 15 residues; this must be handled
explicitly when extracting a template.

The local predicted structures were not promoted to experimental evidence. The
local QA found 13 exact peptide-to-manifest matches, 17 clean-layout structures
with unmapped peptides, 10 ambiguous docking complexes, and 10 incomplete
single-chain exports. Their full paths and statuses are recorded in
`processed/pmhc_structure_qa.csv`.
