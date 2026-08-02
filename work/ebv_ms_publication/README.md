# EBV–MS molecular-mimicry rebuild

This is the clean publication workspace for the rebuilt analysis. The primary
scope is **MHC-II presentation by HLA-DRB1*15:01**, with the older MHC-I and
generic-TCR analyses retained only as legacy material until independently
validated.

## Current status

- Raw IEDB API responses are preserved in `raw/`.
- Normalized, row-level tables are generated in `processed/`.
- The original pipeline is not overwritten; its path and checksum are recorded
  in `processed/manifest.json`.
- No simulated values, legacy labels, or structural scores are publication
  inputs at this stage.
- The next structural phase should begin from
  `processed/pmhc_candidate_manifest.csv` and
  `processed/pmhc_candidate_peptides.fasta`, not from the legacy filename-only
  PDB labels.

## Reproducibility

Run `src/prepare_datasets.py` with the bundled Python environment after raw
inputs change. The script only creates derived files under `processed/`.

## Evidence rule

Every candidate must retain an IEDB identifier, peptide sequence, source
antigen, HLA restriction, evidence type, and source URL. A peptide can be used
for a primary result only when its class, allele, organism, and assay evidence
are explicit. Anything else remains exploratory or quarantined.
