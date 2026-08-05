# TCR–pMHC calibration package

This package begins the structure-calibration stage. `source_fasta/` contains the exact chain records downloaded from RCSB on 2026-08-04. Use the chain map below to select biological components and remove only explicitly documented construct artifacts.

## Run order

1. CAL_1YMM
2. CAL_1ZGL
3. CAL_2WBJ
4. Hy.2E11 MBP and BALF5 reproduction runs

Do not interpret Hy.2E11 ternary confidence until the three calibration cases have been scored against their experimental structures.

## 2WBJ handling

2WBJ fuses the peptide to the N-terminus of the deposited TCR-beta construct. For a prediction input, use the peptide `MDFARVHFISALHG` and the beta chain beginning `AVVSQHPS...`; omit `MDFARVHFISALHGSGGGSGGGGG` from the beta chain. This preserves the deposited peptide while preventing the artificial linker from being treated as a TCR sequence.

`FARVHFISALHG` is the literature-reported ENGA epitope core. The leading `MD` is retained only for strict recovery of the deposited 2WBJ peptide construct.

## Evaluation rule

For each method, record chain mapping, preprocessing, version, pMHC geometry, TCR–pMHC interface recovery, and any confidence score. Calculate RMSD only after aligning the same biological components; never include expression tags or fusion linkers.
