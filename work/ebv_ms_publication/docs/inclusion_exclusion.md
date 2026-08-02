# Inclusion, exclusion, and quarantine rules

These rules are written before re-running comparisons so that candidate
selection cannot be changed to fit a preferred result.

## Include in the primary MHC-II table

1. Peptide is a linear amino-acid sequence containing only the 20 standard
   residues.
2. HLA restriction is explicitly HLA-DRB1*15:01 (not a broad DRB1 family label).
3. MHC class is explicitly II.
4. EBV records identify the source organism as Human herpesvirus 4 (Epstein–Barr
   virus), or human records identify the source organism as Homo sapiens.
5. The row retains an IEDB assay or epitope identifier and a source reference.

## Exclude from the primary table

- MHC-I records, HLA-A*02:01 records, and any mixed-class analysis.
- The old EBV Top-100 FASTA as evidence of MHC-II presentation; it is mostly
  9–10-mers and is retained as a legacy input only.
- Records with missing allele, class, organism, or assay outcome.
- Nonstandard residues, modified peptides, or sequences that cannot be mapped
  to a source protein without a documented reason.
- Susceptibility genes (for example, HLA or immune-regulatory genes) treated as
  if they were myelin antigens. They can be background proteins only.

## Quarantine pending validation

- All legacy PDBs whose peptide/allele mapping is inferred only from filenames.
- Any AlphaFold pLDDT/PAE or docking score not traceable to an actual run with
  recorded inputs and software settings.
- The existing neural-network classifier until leakage, grouping, and labels
  are independently audited.
