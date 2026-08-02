# Locked study scope (version 0.1; 2026-08-01)

## Primary question

Among peptides with defensible HLA-DRB1*15:01 MHC-II evidence, do EBV-derived
peptides share sequence or structural features with human myelin-associated
peptides that are relevant to multiple sclerosis, and can those candidates be
prioritized for follow-up testing?

## Primary analysis population

- HLA: HLA-DRB1*15:01, explicitly treated as a one-chain allele restriction.
- Presentation class: MHC-II only.
- Viral source: Human herpesvirus 4 (Epstein–Barr virus).
- Human comparison: myelin basic protein (MBP), proteolipid protein (PLP), and
  myelin-oligodendrocyte glycoprotein (MOG) peptides, with other human peptides
  retained as background controls.
- Evidence hierarchy: assay-level IEDB records first; aggregate workbook rows
  are used for the human candidate inventory and cross-checks.

## What this study will not claim

- It will not claim that peptide similarity proves molecular mimicry,
  pathogenicity, disease causation, or clinical risk.
- It will not treat HLA-A*02:02 or HLA-DRB1*15:02 as risk alleles. In this
  project they are legacy control labels and are excluded from the primary
  analysis until their exact intended meaning is documented.
- It will not combine MHC-I and MHC-II scores, or compare unrelated generic TCR
  docking runs across MHC classes.
- It will not use simulated scores, random energies, placeholder binding
  values, or unverified structure labels as results.

## Pre-registered primary outcome

The primary output is a ranked, provenance-preserving list of EBV–human peptide
pairs, accompanied by transparent sequence-similarity metrics and sensitivity
analyses. Structural modeling is secondary and hypothesis-generating only.
