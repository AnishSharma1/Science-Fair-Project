# Corrected pMHC-II structure-generation protocol

## What is ready

The project now has 86 auditable HLA-DRA : HLA-DRB1*15:01 : peptide inputs:

- 12 EBV Tier 1 DRB1*15:01 T-cell candidates
- 21 EBV Tier 2 DRB1*15:01 MHC-ligand candidates
- 53 coordinate-validated human myelin candidates

The mature DRA and DRB1*15:01 sequences were taken from chains A and B of the
experimental 8TBP complex. Every peptide is the complete IEDB sequence; no
class-I 9-mer trimming was applied.

## Files

- `processed/pmhc_colabfold_inputs.fasta`: all candidates in ColabFold's
  colon-delimited multimer format
- `processed/colabfold_inputs/`: one FASTA per candidate
- `processed/pmhc_colabfold_manifest.csv`: candidate-to-input provenance
- `processed/pmhc_colabfold_metadata.json`: sequence and template metadata

## Recommended pilot

Run three EBV Tier 1 candidates first (`EBV_TCELL_2268741`,
`EBV_TCELL_2268933`, and `EBV_TCELL_2268720`) plus one human myelin control.
Use AlphaFold2-Multimer v3, at least five seeds, and record the full output
archive, model ranking, pTM/ipTM, per-residue pLDDT, and PAE. Do not use the
old simplified PDBs as prediction outputs.

## ColabFold settings

1. Use the sequence from the candidate FASTA file as the query sequence.
2. Keep the chain delimiter `:` intact.
3. Select `alphafold2_multimer_v3` (or `auto` when the notebook maps complexes
   to multimer v3).
4. For the initial pilot, use `mmseqs2_uniref_env`, `unpaired_paired`, five
   seeds, and no Amber relaxation.
5. Save the result archive and its JSON score/PAE files. A coordinate file is
   not publication-ready until it passes the project QA script and the peptide
   chain is verified against the manifest.

This workspace contains the corrected inputs and provenance, but no prediction
engine was available locally, so this preparation step intentionally does not
invent coordinates or call a model a finished structure.

## Batch mode for all candidates

Use either `processed/pmhc_colabfold_inputs.fasta` or the two-column file
`processed/pmhc_colabfold_batch.csv` with the supplied cell
`docs/colabfold_batch_cell.py`. The cell uploads one file once, validates all
86 IDs and chain delimiters, runs them through one ColabFold-Multimer
workflow, and copies the complete output to `My Drive/ebv_ms_pmhc_batch` (or a
timestamped sibling if that folder already exists).

The batch keeps the publication-oriented settings: Multimer v3, five model
weights, three recycles, one seed, no Amber relaxation, and full output
retention. It may take hours because it still predicts 86 distinct complexes;
"one batch" avoids repetitive manual setup but does not make 86 predictions a
single model evaluation. If the MMseqs2 service throttles a long submission,
split the same FASTA or CSV into smaller blocks without changing the IDs.
