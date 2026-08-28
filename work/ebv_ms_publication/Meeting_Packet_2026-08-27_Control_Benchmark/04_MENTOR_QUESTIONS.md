# Decisions to ask the mentor

1. **Baseline gate:** Should the next benchmark require the composite to improve on TCR-facing sequence identity, or merely match it while adding mechanistic interpretability?
2. **Feature dominance:** Two outer folds selected 100% physicochemical mismatch; the Hy.1B11 fold selected 50% exposed C-alpha and 50% mismatch. Is that acceptable, or should every frozen score require a nonzero structural contribution?
3. **Specificity:** Which published same-TCR nonrecognized ligands or same-assay/HLA negatives can serve as N1/N2 controls? The present N3 decoys test ranking, not specificity.
4. **DQ structure layer:** Should we expand the DQA1*01:02/DQB1*05:02 structural pool until each Hy.1B11 positive has at least five decoys, or predeclare that layer as unavailable?
5. **Incomplete AlphaFold panels:** Preserve this benchmark as an honest incomplete snapshot, then create a new version with a completely fresh manifest, or continue using only the completed layers?
6. **Discovery unlock:** Confirm that any future pass should regenerate rankings separately within each HLA, with no cross-allele consensus.

## Recommended position

Do not retroactively change the current trust gate after seeing the results. Treat the sequence-baseline requirement and any minimum structural-weight rule as prospective additions to a version 2 benchmark.
