# AlphaFold Server pMHC descriptive analysis

This folder extracts model availability, exact chain/peptide sequence QA, server confidence fields, peptide confidence, peptide--HLA contact proxies, and within-job pose consistency from completed pMHC downloads.

It does **not** perform docking or screening, generate a biological candidate ranking, or infer presentation, shared-TCR binding, cross-reactivity, molecular mimicry, or MS mechanism. The server ranking score selects a representative model only within its own five-sample job.

## Completed model groups

- legacy_candidate_pmhc: **60** completed jobs
- new_human_background_pmhc: **9** completed jobs

`af3_pmhc_job_summary.csv` has one row per completed job; `af3_pmhc_sample_metrics.csv` has one row per server sample; and `af3_pmhc_cohort_summary.csv` reports descriptive cohort medians only. The 18 missing seed-03 background jobs remain in `../expanded_background/alphafold_server_seed_03_download_inventory.csv` and are neither imputed nor retried.
