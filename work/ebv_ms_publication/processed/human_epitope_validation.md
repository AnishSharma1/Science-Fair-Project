# Human epitope provenance validation

- Workbook epitopes checked: **102**
- Epitopes with at least one Homo sapiens source mapping: **102**
- Human source mappings retained: **349**
- Coordinate-length mismatches: **14**
- Epitopes quarantined for missing coordinates: **14**
- Epitopes with multiple human source mappings: **64**

Every retained mapping has a source organism of Homo sapiens. Most have
a coordinate span equal to the reported peptide length; the mismatches
must be reviewed rather than silently discarded. Multiple mappings
are preserved because the same epitope may be represented by several
sequence accessions or isoform records. The analysis should use the
mapping table rather than silently choosing one accession.

This validates provenance and coordinates; it does not prove that the
peptide is disease-specific, pathogenic, or naturally presented in MS.
