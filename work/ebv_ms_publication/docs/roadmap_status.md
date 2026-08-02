# Publication rebuild status

## Completed in this push

- Locked the primary question and MHC-II/HLA-DRB1*15:01 scope.
- Created an untouched legacy quarantine record and checksum.
- Wrote inclusion, exclusion, and quarantine rules before re-analysis.
- Preserved raw IEDB API responses for EBV and human HLA-DRB1*15:01 records.
- Normalized the human workbook block and assay-level API tables.
- Added a provenance manifest and data dictionary.
- Ran sequence-level exploratory comparisons with explicit scoring rules.
- Ran quality checks for class, allele, sequence alphabet, and provenance.
- Validated IEDB source-antigen mappings for all 102 human epitopes; quarantined
  14 records with missing coordinate metadata.
- Resolved all 33 distinct human source-protein accessions through UniProtKB;
  confirmed MBP, PLP1, and MOG identifiers.
- Rechecked the eight EBV source publications through PubMed.
- Ran a 10,000-permutation exact 5-mer sensitivity null model.
- Added a transparent physicochemical similarity model and 10,000 label
  permutations; the current signal was not positive.
- Built deterministic length-matched EBV negative controls.
- Built the 86-row pMHC candidate manifest and readiness gate for the next
  structural push.
- Audited 50 local MHC-II PDB files against peptide sequence and DRA/DRB chain
  roles; found duplicated-DRB/missing-DRA layouts in eight peptide-matched
  legacy models.
- Retrieved experimental RCSB templates 8TBP and 5V4M and generated explicit
  pMHC modeling inputs.

## Next gates before structural modeling

1. Validate human peptide-to-protein mappings and accessions.
2. Reconcile the EBV shortlist against IEDB record pages and primary papers.
3. Define a null model and statistical analysis plan for peptide similarity.
4. Add exact controls and sensitivity analyses without reusing legacy labels.
5. Only then decide whether any pMHC structural modeling is justified.
