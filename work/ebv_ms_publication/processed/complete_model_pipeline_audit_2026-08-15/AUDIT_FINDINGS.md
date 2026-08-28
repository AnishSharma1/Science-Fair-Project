# Complete pMHC modeling audit

## Direct answer

The Monday AF3 pMHC models were preserved, but the focused same-register analysis did not read their folder. It read the later `json1folds/json2folds` collection instead. This omitted seven otherwise eligible frozen-shortlist pairs from the nine-pair result table. The combined audit now includes the Monday models and removes copied duplicate jobs before scoring.

## Reconciled model inventory

- Saved job directories discovered: **174** (focused_af3_rerun: 60, monday_af3_full: 80, new_background_af3: 10, new_background_af3_seed03_completion: 24).
- Unique candidate/peptide/seed jobs after deduplication: **138** (focused_af3_rerun: 24, monday_af3_full: 80, new_background_af3: 10, new_background_af3_seed03_completion: 24).
- Canonical study jobs with complete parsed pMHC outputs: **133**.
- Frozen shortlist pairs: **32**; primary-allele register-eligible: **16**; structurally evaluated after the audit: **16**.

## What Monday changes

Monday unlocks **7** eligible pairs that were absent from the focused rerun table:

- Rank 17: `EBV_TCELL_119155::HUMAN_MYELIN_116995`
- Rank 25: `EBV_TCELL_119155::HUMAN_MYELIN_117032`
- Rank 20: `EBV_TCELL_119155::HUMAN_MYELIN_67907`
- Rank 21: `EBV_TCELL_119155::HUMAN_MYELIN_112037`
- Rank 9: `EBV_TCELL_2268933::HUMAN_MYELIN_5516`
- Rank 29: `EBV_TCELL_2268934::HUMAN_MYELIN_112226`
- Rank 16: `EBV_TCELL_119155::HUMAN_MYELIN_116976`

## Best current pairings

The audit rank is lexicographic: fraction of AF3 comparisons below 2 A (higher first), median exposed-position RMSD (lower first), then exposed-position property similarity (higher first). It is a transparent prioritization order, not a biological probability.

| Audit rank | Pair | Robustness tier | <2 A fraction | >=5 A fraction | Median exposed RMSD (A) | Exposed property similarity |
|---:|---|---|---:|---:|---:|---:|
| 1 | `EBV_TCELL_950::HUMAN_MYELIN_112214` | tier_A_robust | 0.90 | 0.10 | 0.758 | 0.874 |
| 2 | `EBV_TCELL_2268741::HUMAN_MYELIN_117032` | tier_B_mixed | 0.60 | 0.40 | 1.208 | 0.693 |
| 3 | `EBV_TCELL_119155::HUMAN_MYELIN_116995` | tier_C_unstable_partial_pose | 0.26 | 0.60 | 8.809 | 0.717 |
| 4 | `EBV_TCELL_119155::HUMAN_MYELIN_117032` | tier_C_unstable_partial_pose | 0.12 | 0.76 | 10.296 | 0.663 |
| 5 | `EBV_TCELL_119155::HUMAN_MYELIN_67907` | tier_C_unstable_partial_pose | 0.08 | 0.84 | 10.180 | 0.744 |
| 6 | `EBV_TCELL_149795::HUMAN_MYELIN_117032` | tier_C_unstable_partial_pose | 0.06 | 0.81 | 13.033 | 0.744 |
| 7 | `EBV_TCELL_1862913::HUMAN_MYELIN_112319` | tier_C_unstable_partial_pose | 0.05 | 0.92 | 6.650 | 0.552 |
| 8 | `EBV_TCELL_950::HUMAN_MYELIN_114366` | tier_D_no_consistent_pose | 0.00 | 1.00 | 5.630 | 0.707 |
| 9 | `EBV_TCELL_119155::HUMAN_MYELIN_112037` | tier_D_no_consistent_pose | 0.00 | 0.82 | 6.177 | 0.744 |
| 10 | `EBV_TCELL_2268933::HUMAN_MYELIN_5516` | tier_D_no_consistent_pose | 0.00 | 1.00 | 6.351 | 0.683 |

## Claim boundary and remaining limitation

Length/bin-matched structural background comparisons are now available for **8** of the **16** primary-allele register-eligible shortlist pairs. Background summaries weight each unique comparator candidate once; repeated server jobs and model samples remain technical sensitivity analyses rather than biological replicates. These are descriptive controls with no p-value. The score sheet therefore prioritizes candidates for a computational paper and future testing; it does not establish presentation, shared-TCR recognition, cross-reactivity, molecular mimicry, or an MS mechanism.
