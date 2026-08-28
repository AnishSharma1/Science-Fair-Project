# Central Model Library Design

## Goal

Create one readable model-library entry point for the EBV-MS project without moving, renaming, or deleting original ColabFold, TCRmodel2, or AlphaFold 3 files.

## Location and preservation rule

The library will live at `ebv_ms_publication/Model_Library`. It will contain lightweight Finder aliases to source folders plus inventories and README files. The source folders remain authoritative and unmodified:

- legacy ColabFold pMHC artifacts: `ebv_ms_publication/processed`
- TCRmodel2 ternary predictions and DockQ work: `Documents/New project/outputs/ebv_ms_model_package`
- AF3 downloaded results: `iCloud Downloads/folds_2026_08_10_*`
- AF3 submissions and manifests: `ebv_ms_publication/af3_migration_2026-08-10`

## Library layout

```text
Model_Library/
  00_README.md
  01_AF3_pMHC/
    README.md
    original-downloads.alias
    submission-batches.alias
    af3_download_inventory.tsv
  02_ColabFold_pMHC_Legacy/
    README.md
    legacy-colabfold-source.alias
    legacy_colabfold_inventory.tsv
  03_TCR_pMHC_TCRmodel2/
    README.md
    tcrmodel2-source.alias
    tcrmodel2_inventory.tsv
  04_Experimental_References/
    README.md
    reference-structures.alias
  05_Analysis_and_Manifests/
    README.md
    master_model_inventory.tsv
```

## Classification rules

- `AF3 pMHC` means a three-chain HLA-DRA/HLA-DRB1*15:01/peptide prediction. These are not TCR-docking results.
- `ColabFold pMHC legacy` means prior pMHC inputs, QA tables, geometry-triage outputs, and template/scaffold explorations stored with the former workflow.
- `TCR-pMHC TCRmodel2` means five-component ternary models and their calibration/DockQ files. They must not be labelled as ColabFold outputs.
- Experimental PDB references and stated positive/negative controls are listed separately from predicted structures.
- Every inventory reports source path, method/family, role, status, and original name. It does not recompute scientific metrics or alter structural coordinates.

## Acceptance checks

1. The six library sections and their README files exist.
2. Each alias resolves to its intended original source folder.
3. All AF3 download folders and all source model-package folders are represented in inventories.
4. The master inventory has one header row and links to each item-specific inventory.
5. No original source file has been moved, renamed, or deleted.
