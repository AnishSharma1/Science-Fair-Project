# Structure QA interpretation

- `validated_reference_structure` means that the coordinate file has a short
  peptide chain whose exact sequence matches a validated candidate, one
  DRA-like chain, and one DRB-like chain. It still does not validate the exact
  allele or model accuracy.
- `peptide_match_but_missing_alpha_chain` means the peptide matches, but both
  protein-sized chains look DRB-like relative to the experimental 8TBP
  template. These files must not be used as MHC-II complexes.
- `legacy_ambiguous_complex` includes files with extra/merged chains such as
  docking outputs; these are not pMHC-only structures.
- `incomplete_structure` includes single-chain tetramer exports where chain
  separation was lost.
- `unmapped_peptide` means the structure contains a plausible peptide chain,
  but its sequence is absent from the new evidence-backed candidate manifest.

All predicted structures remain secondary, hypothesis-generating models.
