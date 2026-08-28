# Protein-region naming guide

Use the files ending in **`_WITH_PROTEIN_REGIONS.csv`** when reading, discussing, plotting, or drafting the paper. They retain every original stable ID but add the protein name and exact peptide coordinates near the start of each row.

## Example

`EBV_TCELL_119155` is now labeled:

**EBNA4 (EBNA3B), residues 281-300 — `FIEFVGWLCKKDHTHIREWF`**

Its mapping is a unique exact sequence match to NCBI Protein accession `Q3KST1.1`.

## Main result labels

| Discovery rank | Readable pairing | Stable pair ID |
|---:|---|---|
| 1 | EBNA1 residues 482-496 vs MBP residues 245-263 | `EBV_TCELL_950::HUMAN_MYELIN_112214` |
| 2 | Triplex capsid protein 1 residues 289-302 vs PLP1 residues 218-249 | `EBV_TCELL_2268741::HUMAN_MYELIN_117032` |
| 3 | EBNA4 residues 281-300 vs MBP residues 195-213 | `EBV_TCELL_119155::HUMAN_MYELIN_116995` |
| 4 | EBNA4 residues 281-300 vs PLP1 residues 218-249 | `EBV_TCELL_119155::HUMAN_MYELIN_117032` |
| 5 | EBNA4 residues 281-300 vs MBP residues 278-297 | `EBV_TCELL_119155::HUMAN_MYELIN_67907` |

## Numbering rule

- Coordinates are 1-based and inclusive.
- Each coordinate is tied to the accession listed in the same row.
- EBV coordinates come from an exact peptide match to the recorded NCBI/UniProt protein sequence.
- Human myelin coordinates normally use the canonical UniProt sequence named in the candidate manifest.
- MBP has multiple isoforms, so the same peptide can have different residue numbers in different papers. For example, `LSRFSWGAEGQRPGFGYGG` is MBP residues 245-263 in canonical UniProt `P02686`; shorter MBP isoforms use lower numbers. Never quote an MBP residue interval without its accession or isoform convention.
- `HUMAN_MYELIN_118650` does not occur in canonical MOG `Q16653`; it is labeled MOG residues 145-160 using IEDB's coordinate-validated accession `AAB08089.1`, and this exception is explicitly flagged in the crosswalk.

## Files

- `01_COMPLETE_32_PAIR_SCORECARD_WITH_PROTEIN_REGIONS.csv`: the main 32-pair score sheet with readable pair names.
- `02_ALL_1000_STRUCTURE_COMPARISONS_WITH_PROTEIN_REGIONS.csv`: every AF3 cross-model comparison with readable pair names.
- `03_ALL_150_SAVED_AF3_JOB_FOLDERS_WITH_PROTEIN_REGIONS.csv`: every saved job with a readable candidate/control name.
- `04_UNIQUE_AF3_JOB_QUALITY_SUMMARY_WITH_PROTEIN_REGIONS.csv`: one-row-per-job quality summary with readable names.
- `05_INDIVIDUAL_AF3_MODEL_QUALITY_METRICS_WITH_PROTEIN_REGIONS.csv`: all individual AF3 models with readable names.
- `../protein_region_annotations/candidate_protein_region_annotations.csv`: master ID-to-protein-region crosswalk, including discovery candidates, background comparators, GlialCAM controls, and the excluded non-pMHC decoy.

## Scientific boundary

These labels identify the parent protein segment represented by each peptide. They do not establish T-cell cross-reactivity, molecular mimicry, antigen processing, or disease mechanism.
