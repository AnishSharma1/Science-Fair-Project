# AlphaFold Server batch upload

Upload `alphafold_server_batch.json` as one batch. It contains 90 jobs: each of the 18 biological conditions is duplicated across five fixed, one-seed jobs (90 model seeds total). AlphaFold Server accepts only one seed per job.

If account 1 has already accepted the first 30 jobs, upload `alphafold_server_account_2_remaining_30.json` to account 2 and `alphafold_server_account_3_remaining_30.json` to account 3. These two files contain exactly the remaining 60 jobs shown after `CAL_1YMM` in the master batch.

All protein chains have `useStructureTemplate: false`. This blocks direct template use, including the three calibration PDB entries. It does not make the models fully training-set-blind, so report these as template-excluded calibration results rather than held-out performance.

The batch intentionally excludes influenza and reovirus Hy.2E11 mimics because their allele restriction is not sufficiently resolved for an allele-specific structural input.

After the server finishes, download the entire result archive without renaming files. The analysis order is: pMHC controls, ternary calibrators, Hy.2E11 positives, then decoys. Do not treat a high-confidence decoy as a demonstrated binder or a low-confidence decoy as a demonstrated non-binder.
