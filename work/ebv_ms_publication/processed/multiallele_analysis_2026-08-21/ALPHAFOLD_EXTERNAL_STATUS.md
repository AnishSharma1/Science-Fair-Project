# AlphaFold Server external status

- The original 20-job retry was partially completed by AlphaFold Server.
- Thirteen additional retry outputs were downloaded on 2026-08-22; the canonical matrix is now **143/150**.
- The remaining **7** jobs are treated as persistent missingness and are never imputed.
- The 30-job fixed-seed robustness JSON is prepared but intentionally remains unsubmitted until its frozen control manifest is reviewed.
- A prepared JSON is not counted as a completed AlphaFold job; only downloaded folders passing reinventory enter the model matrix.
