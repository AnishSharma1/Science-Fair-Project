# AF3 reanalysis: what the downloads do and do not support

## Bottom line

The 90 AlphaFold Server jobs contain a reproducible **Hy.2E11–BALF5 docking hypothesis**, but the controls show that this template-excluded AF3 setup does **not** recover known pMHC or ternary geometry. The hypothesis is therefore not structural evidence for EBV/MBP cross-reactivity and should be tested independently with TCRmodel2 before any biological interpretation.

## Calibration check: high confidence did not mean correct structure

Exact-sequence pMHC controls were aligned to their deposited structures using C-alpha atoms:

| Control | Experimental reference | Mean pMHC RMSD across 5 seed jobs |
|---|---:|---:|
| MBP / DRB1*15:01 | 1YMM | 21.43 A |
| MBP / DRB5*01:01 | 1ZGL | 23.83 A |
| ENGA / DRB1*15:01 | 2WBJ | 25.14 A |

The ternary controls were likewise poor: mean TCR placement RMSD after pMHC alignment was 62.93 A (1YMM), 69.12 A (1ZGL), and 54.60 A (2WBJ). This invalidates use of the AF3 confidence values as calibrated docking accuracy for the Hy.2E11 models.

## The useful internal-consistency signal

Five independently seeded Hy.2E11/BALF5 jobs were produced for both peptide representations:

| Hy.2E11 condition | Mean TCR-pMHC pair ipTM | Stable TCR residues contacting peptide in all 5 seed jobs |
|---|---:|---|
| BALF5 biological 15-mer | 0.801 | A92D, A93S, A94G, A95G, A96S, A97Y, B28Q, B30T, B73L, B96W, B97P, B103Y, B105Y |
| BALF5 deposited 14-mer | 0.825 | Same 13 residues |

The full-15 and deposited-core-14 BALF5 inputs have an identical stable TCR-residue contact set across seeds. This says AF3 is internally consistent about its own BALF5 docking solution, including an alpha-chain cluster at A92-A97 and beta-chain positions B96/B103/B105.

The MBP run shares seven stable contact residues with that BALF5 solution: **A93S, A94G, A95G, A96S, A97Y, B96W, and B103Y**. These are candidate shared-recognition residues for an independent model comparison.

## Why this is still not selectivity evidence

The same contact-overlap test is not specific enough to separate positives from the designed decoys. For example, Hy.2E11/ENGA shares six of the eight MBP stable contact residues (Jaccard 0.667), whereas MBP shares seven of 13 with BALF5 (Jaccard 0.500). The decoys are not experimentally verified negatives in any case.

So the correct conclusion is narrow: **AF3 offers a repeatable candidate interface, not a credible ranking of cognate versus noncognate pMHCs.**

## Next discriminating test

Use the prepared template-excluded TCRmodel2 calibration FASTAs first. If that method reproduces the three experimental controls, compare whether it independently returns the same candidate Hy.2E11 alpha-chain A92-A97 and beta-chain B96/B103/B105 contact cluster for MBP and BALF5. Agreement across the two methods would justify a focused experimental or computational follow-up; disagreement means retire the AF3 contact hypothesis.

## Files behind this assessment

- `af3_calibrator_structural_metrics.tsv` — ternary structural calibration
- `af3_pmhc_control_structural_metrics.tsv` — pMHC-only structural calibration
- `hy_tcr_peptide_contact_consensus.json` — residue-level five-seed contact counts
- `AF3_HY_TCR_PEPTIDE_CONTACT_CONSENSUS.md` — detailed contact consensus
- `analyze_hy_interfaces.py` — reproducible contact-analysis script
