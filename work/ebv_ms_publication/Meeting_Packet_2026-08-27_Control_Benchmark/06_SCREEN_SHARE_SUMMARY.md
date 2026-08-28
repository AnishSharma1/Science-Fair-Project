# Held-out HLA-II control benchmark

## What worked

- 48/48 new jobs: exact complete five-model bundles
- 240/240 new model samples: exact-chain and clash-free
- 8/8 completed held-out panels: positive ranked in top 3
- Worst completed rank: 9 to 3 versus exposed C-alpha RMSD alone
- Completed-panel failures: 0

## What prevents a validation claim

- 2 incomplete AlphaFold panels
- 2 Hy.1B11 PDB panels with only two exact-HLA decoys
- TCR-facing sequence identity: 8/8 capture, worst rank 2
- Composite: 8/8 capture, worst rank 3

## Verdict

**Encouraging control recovery; formal trust gate remains `not_evaluable`.**

Weights remain unfrozen. Discovery rankings remain unchanged. No cross-allele consensus was created.

Computational pMHC prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS mechanism.
