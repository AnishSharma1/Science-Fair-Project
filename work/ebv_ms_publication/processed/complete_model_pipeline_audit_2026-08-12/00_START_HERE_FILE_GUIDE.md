# Complete pMHC Audit: Start Here

This folder contains the reconciled record of the ColabFold shortlist and the
AlphaFold 3 pMHC analysis. The numbered filenames are human-readable copies.
The original technical filenames are preserved because the analysis scripts
may depend on them.

For protein names and exact parent-protein residue intervals, begin with
**`00_PROTEIN_REGION_NAMING_GUIDE.md`** and use the matching files ending in
**`_WITH_PROTEIN_REGIONS.csv`**. The original numeric IDs remain in those
tables as stable audit identifiers.

## The three files most people need

1. **`01_COMPLETE_32_PAIR_SCORECARD_WITH_PROTEIN_REGIONS.csv` — Master candidate scorecard**
   - One row = one predeclared EBV–human peptide pair.
   - Contains all 32 shortlisted pairs, including pairs that could not be
     discovery-ranked because their HLA register or allele context was not
     sufficiently resolved.
   - Use this file to answer: “What are all the pairings, which were eligible,
     and which pairing ranked best?”

2. **`06_AUDIT_RESULTS_SUMMARY.md` — Plain-language audit conclusion**
   - Explains what happened to the Monday models.
   - Lists the seven eligible pairings restored by the audit.
   - Summarizes the leading pairs and the remaining scientific limitation.

3. **`02_ALL_1000_STRUCTURE_COMPARISONS_WITH_PROTEIN_REGIONS.csv` — Full structural evidence**
   - One row = one AF3 model-sample comparison between the EBV pMHC and human
     pMHC belonging to the same candidate pair.
   - Contains all 1,000 direct P1–P9 geometry comparisons used to summarize the
     16 register-eligible pairs.
   - Use this file to verify medians, outliers, and seed/sample sensitivity.

## Supporting audit files

### `03_ALL_150_SAVED_AF3_JOB_FOLDERS.csv` — Saved-model provenance ledger

- One row = one AF3 job folder found on the computer.
- Includes 80 Monday jobs, 60 focused-rerun jobs, and 10 background/other jobs.
- Marks which folders are complete and which folders are copied duplicates of
  the same candidate, peptide, and server seed.
- Use this file to answer: “Where did this model come from, and was it counted
  more than once?”

### `04_UNIQUE_AF3_JOB_QUALITY_SUMMARY.csv` — One-row-per-job quality report

- One row = one unique, mapped, complete AF3 study job after exact copies were
  removed.
- Contains 109 job rows.
- Summarizes confidence, peptide pLDDT, peptide–HLA contact coverage, and pose
  consistency across the five samples produced by that job.
- Use this file to compare overall job quality or identify unstable jobs.

### `05_INDIVIDUAL_AF3_MODEL_QUALITY_METRICS.csv` — Model-sample quality report

- One row = one individual AF3 model sample.
- Contains 545 rows: five samples for each of the 109 unique mapped jobs.
- Retains sample-level ranking score, ipTM, pTM, peptide pLDDT, and contact
  information.
- Use this file when a job-level median is not detailed enough.

## How to read the master scorecard

### Identity columns

- `pair_id`: stable name joining the EBV and human candidate IDs.
- `ebv_candidate_id` / `human_candidate_id`: project identifiers for the two
  peptides.
- `ebv_peptide` / `human_peptide`: full modeled peptide sequences.
- `ebv_p1_p9_core` / `human_p1_p9_core`: the frozen nine-residue register used
  for the equivalent-position comparison.

### Eligibility and provenance columns

- `register_eligible_primary_allele`: whether both members had a resolved
  primary-allele P1–P9 hypothesis suitable for discovery ranking.
- `monday_af3_both_models_available`: whether Monday contained both pMHCs.
- `focused_rerun_both_models_available`: whether the focused rerun contained
  both pMHCs.
- `combined_af3_both_models_available`: whether the reconciled model collection
  contains both pMHCs.
- `audit_status`: why a pair was structurally evaluated or why it was not
  discovery-rankable.

### Same-register sequence columns

- `same_register_property_similarity`: chemistry similarity across all P1–P9
  positions.
- `candidate_exposed_identity_fraction`: exact identity across the descriptive
  P2/P3/P5/P7/P8 position set.
- `candidate_exposed_property_similarity`: amino-acid property similarity
  across that same position set.

“Candidate-exposed” is a computational label. It does not prove that these
residues contact a TCR.

### Structural columns

- `candidate_exposed_ca_rmsd_A_median`: median geometric difference across
  P2/P3/P5/P7/P8 after the two HLA grooves are fitted. Lower is more similar.
- `candidate_exposed_rmsd_lt_2A_fraction`: fraction of comparisons below 2 Å.
  Higher means the similar pose was recovered more consistently.
- `candidate_exposed_rmsd_ge_5A_fraction`: fraction of grossly different poses.
  A large value signals instability even if a few favorable samples exist.
- `af3_cross_sample_geometry_count`: number of model-sample comparisons behind
  the summary statistics.

### Robustness labels

- `tier_A_robust`: at least 80% below 2 Å, no more than 20% at or above 5 Å,
  and a median below 2 Å.
- `tier_B_mixed`: at least 50% below 2 Å and a median below 2 Å, but substantial
  discordant poses remain.
- `tier_C_unstable_partial_pose`: at least one favorable subset exists, but the
  overall ensemble does not consistently recover it.
- `tier_D_no_consistent_pose`: no comparison falls below 2 Å.

These are technical-consistency tiers, not probabilities of biological
cross-reactivity.

## Current interpretation

- **Sole robust lead:** `EBV_TCELL_950::HUMAN_MYELIN_112214`
  (EBNA1–MBP).
- **Mixed secondary result:**
  `EBV_TCELL_2268741::HUMAN_MYELIN_117032`.
- **Monday-only recovered pairs:** all are unstable or show no consistent pose;
  none displaces the lead.
- **Main remaining limitation:** no eligible pair has a completed
  length/bin-matched AF3 structural background comparator.

The outputs support computational pMHC candidate prioritization. They do not
establish natural peptide presentation, shared-TCR recognition, T-cell
activation, cross-reactivity, molecular mimicry, or an MS mechanism.

## Original technical filenames

| Human-readable copy | Original technical filename |
|---|---|
| `01_COMPLETE_32_PAIR_SCORECARD.csv` | `master_pair_score_sheet.csv` |
| `02_ALL_1000_STRUCTURE_COMPARISONS.csv` | `combined_same_register_geometry.csv` |
| `03_ALL_150_SAVED_AF3_JOB_FOLDERS.csv` | `complete_af3_job_inventory.csv` |
| `04_UNIQUE_AF3_JOB_QUALITY_SUMMARY.csv` | `canonical_af3_job_summary.csv` |
| `05_INDIVIDUAL_AF3_MODEL_QUALITY_METRICS.csv` | `canonical_af3_sample_metrics.csv` |
| `06_AUDIT_RESULTS_SUMMARY.md` | `AUDIT_FINDINGS.md` |
