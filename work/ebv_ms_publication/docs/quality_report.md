# Initial quality-control report (2026-08-01)

## Source inventory

| Table | Rows | Unique peptides | Primary use |
|---|---:|---:|---|
| Human HLA-DRB1*15:01 workbook block | 102 | 102 | Human self-antigen inventory |
| EBV HLA-DRB1*15:01 T-cell API records | 31 | 19 | Primary EBV candidate evidence |
| EBV HLA-DRB1*15:01 MHC-ligand API records | 21 | 21 | Orthogonal binding evidence |
| Human HLA-DRB1*15:01 T-cell API records | 1,139 | 272 | Broad comparator, not yet modeled |

The human workbook block contains 65 myelin-candidate rows and 37 human
background rows. The EBV T-cell table contains 12 unique peptides with a
positive qualitative outcome after collapsing duplicate peptide sequences;
the remaining records are retained as negative controls. All EBV records in
the normalized tables report MHC class II and HLA-DRB1*15:01 explicitly.

## Important interpretation limits

1. The workbook counts references and assays, but its rows are aggregate IEDB
   epitope records. They are not 102 independent biological replicates.
2. The API returns multiple assay rows for the same peptide. Sequence-level
   analyses therefore deduplicate by peptide and retain the assay provenance.
3. The EBV positive set is small and heterogeneous across source proteins and
   publications. It supports candidate prioritization, not a prevalence or
   disease-risk estimate.
4. The exploratory similarity table uses a transparent local-alignment score
   (match +2, mismatch −1, gap −2) and is explicitly not an MHC-binding or
   T-cell-receptor model. Short matches are prevented from dominating the
   shortlist by requiring at least five aligned residues for ranking.
5. No AlphaFold, docking, neural-network, or simulated score has entered these
   tables.

## Immediate validation flags

- Confirm whether each human workbook peptide is truly a myelin protein and
  map the accession/coordinates before treating it as a self-antigen.
- Preserve exact source-protein accessions for all EBV rows; avoid collapsing
  distinct viral proteins solely by display name.
- Add a negative-control design and a pre-specified similarity null model before
  any significance claim.
- Reconcile API records against the IEDB web record pages and publication full
  texts for the final candidate shortlist.
