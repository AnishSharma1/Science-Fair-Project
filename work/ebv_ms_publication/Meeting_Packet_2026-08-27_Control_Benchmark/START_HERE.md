# August 27 control-benchmark call packet

## Open first

1. `01_CALL_CARD.md`
2. `06_SCREEN_SHARE_SUMMARY.md`
3. `02_CONTROL_RESULTS.csv` if exact panel-level ranks are requested

## Run after the call

- `07_DEEP_RESEARCH_REQUEST_BENCHMARK_V2.md`: paste this complete request into a deep-research tool to identify additional controls, specificity negatives, exact-HLA DQ structures, and a prospective benchmark-v2 design.

## Bottom line

The new control-based score recovered the established positive within the top 3 in all **8 of 8 completed held-out panels**, reducing the worst completed-panel rank from **9 to 3** relative to exposed C-alpha RMSD alone.

That is encouraging control recovery, but it is **not a formal benchmark pass**. Two AlphaFold panels remain incomplete and both Hy.1B11 PDB panels have too few exact-HLA decoys. Simple TCR-facing sequence identity also recovered **8 of 8** controls, with a worst rank of 2, so the present benchmark does not yet show that AlphaFold-derived 3D geometry adds independent predictive value.

## Current decision state

- Trust gate: `not_evaluable`
- Completed-panel failures: 0
- Candidate weights frozen: no
- Discovery reranking allowed: no
- Discovery rankings changed: no
- Cross-allele consensus created: no

## Claim boundary

These are computational pMHC control-recovery results. They do not establish antigen presentation, TCR binding, T-cell activation, cross-reactivity, molecular mimicry, or an MS mechanism.
