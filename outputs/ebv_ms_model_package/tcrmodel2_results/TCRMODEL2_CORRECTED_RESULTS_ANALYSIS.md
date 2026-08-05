# TCRmodel2 corrected-results analysis

## Inputs processed

The five archives in iCloud `Downloads/TCRdocks` were copied and extracted without changing the originals:

| Run | Model confidence | TCR-pMHC ipTM | pMHC templates |
|---|---:|---:|---|
| 1ZGL calibration, reference excluded | 0.857 | 0.816 | 6CQL, 6R0E, 1FYT, 2IAM |
| 2WBJ calibration, reference excluded | 0.771 | 0.693 | 6R0E, 1FYT, 4Y19, 6CQL |
| Hy.2E11 / MBP / DRB1 | 0.845 | 0.826 | 2WBJ, 6R0E, 1FYT, 4Y19 |
| Hy.2E11 / BALF5 / DRB5 replicate 1 | 0.845 | 0.818 | 1ZGL, 6CQL, 6R0E, 1FYT |
| Hy.2E11 / BALF5 / DRB5 replicate 2 | 0.862 | 0.839 | 1ZGL, 6CQL, 6R0E, 1FYT |

The exclusion worked: neither corrected calibrator contains its own direct reference (`1ZGL` or `2WBJ`) in the pMHC-template list.

## Corrected calibration geometry

| Calibration | pMHC RMSD | TCR placement RMSD | TCR-only RMSD | Whole-complex RMSD |
|---|---:|---:|---:|---:|
| 1YMM, earlier run | 10.67 A | 19.85 A | 17.31 A | 17.05 A |
| 1ZGL, corrected exclusion run | 9.93 A | 9.57 A | 6.14 A | 8.30 A |
| 2WBJ, corrected exclusion run | 30.92 A | 51.23 A | 26.38 A | 35.26 A |

The corrected 1ZGL result is a meaningful improvement and is the only reasonably successful template-excluded control so far. 1YMM and 2WBJ do not recover their known structures closely enough to establish general docking accuracy.

## Hy.2E11 hypothesis

The two BALF5 jobs are independent duplicate submissions after the required 11-aa truncation. Their similar top scores (0.845 and 0.862) support internal reproducibility of the BALF5 prediction.

At a 4.5-A heavy-atom contact cutoff, the MBP model and both BALF5 models share 12 Hy.2E11 peptide-contacting residues (Jaccard overlap 0.706). The shared set is:

`D30S, D31I, D32N, D110S, D111G, D112G, D134S, D135Y, E109W, E111S, E112G, E135Y`

Here `D` and `E` are the server output's TCR alpha and beta chains, respectively. This is an internally consistent **candidate shared-recognition interface**, not evidence that the TCR binds both pMHCs.

## Decision and next step

The immediate next calibration should be a corrected 1YMM rerun using the current server-shaped FASTA (standard five headers and MHC binding domains), followed by coordinate comparison. If it improves materially, the Hy interface can be treated as a stronger computational hypothesis. If it remains poor, retain the Hy result only as a lead for experimental validation rather than a docking conclusion.

## Supporting outputs

- `results_analysis/tcrmodel2_calibrator_structural_metrics.tsv`
- `results_analysis/TCRMODEL2_HY_CONTACT_CONSENSUS.md`
- `results_analysis/tcrmodel2_hy_contact_consensus.json`
