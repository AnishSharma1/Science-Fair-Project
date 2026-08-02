# Human source-protein registry validation

The 33 distinct UniProt accessions embedded in the human IEDB workbook labels
were queried against UniProtKB on 2026-08-01. All 33 resolved to Homo sapiens
entries; 31 were reviewed Swiss-Prot entries and 2 were unreviewed TrEMBL
entries.

The three intended myelin proteins are represented by the expected accessions:

- MBP: P02686
- PLP1: P60201
- MOG: Q16653

This confirms that the workbook labels point to real human proteins. It does
not by itself validate every peptide coordinate; coordinate-level validation is
recorded separately in `processed/human_epitope_accession_map.csv`, where 14
epitopes remain quarantined because the IEDB mapping lacks start/end positions.

Source: [UniProtKB REST API](https://rest.uniprot.org/).
