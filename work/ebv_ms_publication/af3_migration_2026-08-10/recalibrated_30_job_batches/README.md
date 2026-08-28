# Recalibrated AlphaFold Server queue: 90 jobs in three batches

This supersedes the earlier 20-job batching scheme. Upload these three files:

1. `pmhc_priority_batch_01_one_seed.json` (30 jobs)
2. `pmhc_priority_batch_02_one_seed.json` (30 jobs)
3. `pmhc_priority_batch_03_one_seed.json` (30 jobs)

The first 86 jobs retain the pre-existing pMHC priority order. The final four
jobs in batch 3 are published GlialCAM 370-389 single-phosphoserine variants:

| Job | Peptide modification |
| --- | --- |
| `GLIALCAM_370_389_pS376` | S376, peptide position 7, phosphoserine (`CCD_SEP`) |
| `GLIALCAM_370_389_pS377` | S377, peptide position 8, phosphoserine (`CCD_SEP`) |
| `GLIALCAM_370_389_pS383` | S383, peptide position 14, phosphoserine (`CCD_SEP`) |
| `GLIALCAM_370_389_pS384` | S384, peptide position 15, phosphoserine (`CCD_SEP`) |

All use GlialCAM 370-389: `ATGRTHSSPPRAPSSPGRSR`. These are exploratory
pMHC structures for a published antibody-mimicry/PTM system, not evidence of
HLA-DRB1*15:01 presentation, TCR recognition, or the MS mechanism.

Each job remains one automatically selected AlphaFold Server seed. Preserve the
job-request JSON, mmCIF, confidence JSON, and PAE for each completed job.
