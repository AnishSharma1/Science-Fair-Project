# Register-resolution results

Overall sensitivity status: `declared_register_dependent`.

This package preserves the eight top-three N3 results, then asks how strongly each result depends on the assumed nine-residue HLA-II register. Alternate windows use sequence features only; their structures were not remodeled.

## Counts

- All-window robust: 0/8
- Local +/-1 robust only: 0/8
- Declared-window only: 8/8
- Not supportive at the declared window: 0/8

## Candidate results

| Target | HLA | Declared rank | Worst local rank | Local capture | All-window capture | Status |
|---|---|---:|---:|---:|---:|---|
| HY03_SEQ_01 | HLA-DRB1*03:01 | 1 | 26 | 2/9 | 6/35 | declared_window_only |
| HY03_SEQ_02 | HLA-DRB1*03:01 | 2 | 20 | 2/9 | 4/49 | declared_window_only |
| HY08_SEQ_01 | HLA-DRB1*08:01 | 1 | 23 | 3/9 | 10/49 | declared_window_only |
| HY08_SEQ_02 | HLA-DRB1*08:01 | 1 | 22 | 3/9 | 16/91 | declared_window_only |
| HY13_SEQ_01 | HLA-DRB1*13:03 | 1 | 15 | 2/9 | 6/49 | declared_window_only |
| HY13_SEQ_02 | HLA-DRB1*13:03 | 2 | 23 | 1/9 | 2/49 | declared_window_only |
| HY15_SEQ_01 | HLA-DRB1*15:01 | 1 | 24 | 4/9 | 7/42 | declared_window_only |
| HY15_SEQ_02 | HLA-DRB1*15:01 | 1 | 15 | 4/9 | 10/49 | declared_window_only |

## Interpretation

A declared-window result can remain promising even when shifts fail, but it is then register-dependent and should be tested with nested peptides before biological interpretation. The 25 N3 pairs retain unknown recognition status and do not become specificity negatives.

Window rescoring tests dependence on the assumed P1-P9 register using the frozen sequence ranking and N3 panel. Alternate windows were not structurally modeled or experimentally resolved. Results do not establish presentation, TCR recognition, specificity, cross-reactivity, molecular mimicry, or MS mechanism.
