# TCRmodel2 calibrator analysis

## Archives processed

| Extracted folder | Source archive | Top-model confidence | TCR-pMHC ipTM |
|---|---|---:|---:|
| `CAL_1YMM` | `final_models.tar.gz` | 0.830 | 0.803 |
| `CAL_1ZGL` | `final_models (2).tar.gz` | 0.865 | 0.822 |
| `CAL_2WBJ` | `final_models (1).tar.gz` | 0.815 | 0.760 |

Each archive was unpacked intact and contains five ranked PDB models, `statistics.json`, TCR template information, and its modeling log.

## Template audit

| Calibration | pMHC templates actually used | Valid template-excluded calibration? |
|---|---|---|
| 1YMM | 2WBJ, 6R0E, 1FYT, 4Y19 | Yes with respect to 1YMM; no exact 1YMM template was used. |
| 1ZGL | **1ZGL**, 6CQL, 6R0E, 1FYT | **No** — exact-reference pMHC template leakage. |
| 2WBJ | **2WBJ**, 6R0E, 1FYT, 4Y19 | **No** — exact-reference pMHC template leakage. |

The apparently strong 1ZGL score therefore cannot be used as independent validation.

## Coordinate comparison of ranked_0 with the experimental PDB

C-alpha atoms were sequence-aligned. pMHC RMSD first fits MHC alpha, MHC beta, and peptide. TCR placement RMSD then evaluates the TCR in that pMHC-aligned frame.

| Calibration | pMHC RMSD | TCR placement RMSD | TCR-only RMSD | Whole-complex RMSD |
|---|---:|---:|---:|---:|
| 1YMM | 10.67 A | 19.85 A | 17.31 A | 17.05 A |
| 1ZGL | 19.13 A | 41.25 A | 17.41 A | 25.93 A |
| 2WBJ | 25.77 A | 62.63 A | 32.86 A | 40.82 A |

## Decision

The server scores are more encouraging than the AF3 scores, but these coordinate comparisons do **not** show reliable recovery of the known complexes. 1ZGL and 2WBJ are additionally invalid as independent calibrators because their exact pMHC references appeared in the template list.

Do not interpret the Hy.2E11 models as evidence of cross-reactivity yet. Before running or using jobs 4–5, resubmit the 1ZGL and 2WBJ controls with their reference PDB IDs absent from the **actual results-page template list**. If the exclusion field does not take effect, a different modeling setup will be necessary.

## Reproducible files

- `tcrmodel2_calibrator_structural_metrics.tsv` — the numerical comparison table
- `compare_tcrmodel2_calibrators.py` — analysis script
- `CAL_*/` — extracted, untouched server outputs
