# Data dictionary

| Field | Meaning | Required for primary use |
|---|---|---|
| `iedb_assay_id` | IEDB assay-level identifier | Yes for assay tables |
| `iedb_epitope_id` | IEDB epitope identifier | Yes |
| `peptide` | Uppercase amino-acid sequence | Yes |
| `peptide_length` | Number of residues | Yes |
| `source_antigen_name` | Protein/antigen associated with the sequence | Yes |
| `source_organism` | Organism of the source antigen | Yes |
| `mhc_class` | I or II | Yes |
| `mhc_allele` | Exact HLA allele string | Yes |
| `outcome` | IEDB qualitative assay outcome | Yes for T-cell tables |
| `pubmed_id` / `reference_id` | Literature provenance | At least one |
| `assay_type` | T-cell response or MHC ligand | Yes |
| `source_url` | API or database provenance | Yes |
| `retrieved_date` | Date the source was retrieved | Yes |

The aggregate human workbook additionally includes `n_references`, `n_assays`,
`candidate_class`, and `source_file`. These are descriptive fields and must not
be mistaken for independent assay replicates.
