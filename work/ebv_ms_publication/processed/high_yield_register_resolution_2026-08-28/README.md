# High-yield register resolution (2026-08-28)

This additive package follows the high-yield N3 analysis. It does not change V1-V3 rankings, replace the prior 12-target package, freeze weights, or unlock discovery.

## Question

For each of the eight sequence-supported candidates, does its top-three rank survive plausible or exhaustive changes to the assumed P1-P9 register when compared with the same frozen 25 N3 pairs?

## Method

- Recompute TCR-facing BLOSUM62, physicochemical mismatch, five-position identity, full-core identity, and full-core BLOSUM62 for every fully contained 9-mer combination.
- Rank each alternative against the candidate's unchanged 25-pair N3 panel with the frozen V3 lexicographic sequence order.
- Abstain from structural tie-breaking for every alternate-register row because alternate registers were not modeled.
- Report both the local +/-1-by-+/-1 neighborhood and the exhaustive window set.
- Provide a proposed nested-peptide design, marked proposed and not ordered.

HLA-II ligands commonly occur as nested peptides of variable length, so a parent peptide alone does not prove its binding register. Experimental work has used nested register peptides and binding measurements to separate alternative MHC-II registers. See [Chicz et al., 1992](https://pubmed.ncbi.nlm.nih.gov/1380674/) and [Mohan et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3256971/).

## Files

- `protocol_lock.json`: frozen scope, inputs, ranking logic, and claim gates.
- `frozen_sequence_target_registry.csv`: exact eight candidates and source sequences.
- `all_window_panel_ranks.csv`: exhaustive register-pair sensitivity matrix.
- `local_shift_panel_ranks.csv`: declared +/-1 window matrix.
- `target_register_summary.csv`: candidate-level robustness and worst ranks.
- `experimental_register_priority.csv`: assay order based on local-shift robustness; not a discovery rerank.
- `experimental_peptide_panel.csv`: proposed parent and nested 9-mer sequences; nothing has been ordered.
- `register_resolution_gate.json`: machine-readable status and permanent lock flags.
- `RESULTS_SUMMARY.md`: concise result table.
- `SHA256SUMS.csv`: deterministic artifact checksums.

## Claim boundary

Window rescoring tests dependence on the assumed P1-P9 register using the frozen sequence ranking and N3 panel. Alternate windows were not structurally modeled or experimentally resolved. Results do not establish presentation, TCR recognition, specificity, cross-reactivity, molecular mimicry, or MS mechanism.
