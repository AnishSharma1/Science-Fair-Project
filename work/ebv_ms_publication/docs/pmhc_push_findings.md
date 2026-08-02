# pMHC push findings (2026-08-01)

The structural audit found a major legacy-modeling issue.

## What the coordinate files actually contain

- The older `ebv pdbs/` and `myelin pdbs/` regular files often contain two
  DRB1-like chains and no DRA-like alpha chain. They are therefore not valid
  MHC-II alpha/beta/peptide complexes, even when the peptide sequence matches a
  current candidate.
- The `MHC Docking Analysis` files contain ambiguous extra/merged chains and
  are not pMHC-only structures.
- The older `Tetramer Docking Analysis` exports have lost chain separation and
  are incomplete for pMHC analysis.
- The `new ebv pdbs/` and `new myelin pdbs/` files have the expected DRA-like +
  DRB-like + peptide layout, but only four unique current candidate peptides
  have exact matches among those files.

## Promoted local reference structures

Only these candidate structures currently pass peptide-sequence and chain-role
QA:

- EBV_MHC_48976
- HUMAN_MYELIN_118650
- HUMAN_MYELIN_112782
- HUMAN_MYELIN_115622
- HUMAN_MYELIN_5516

These are predicted reference models, not experimental structures. They are
useful for checking file conventions and template geometry, not for claiming
binding affinity or disease relevance.

## Modeling consequence

The main pMHC study cannot reuse the old top-five structures as if they were
valid complexes. We need to construct new HLA-DRA + HLA-DRB1*15:01 + peptide
models, using the validated candidate manifest and the experimental 8TBP
template as the primary geometry reference.
