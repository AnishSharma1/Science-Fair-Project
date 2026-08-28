# Pre-meeting rigor design

## Objective

Produce reproducible, meeting-ready artifacts that make the current pMHC prioritization project easier to audit and improve without converting computational outputs into biological claims.

## Approved scope

1. Record HLA-DRB1*15:01 MHC-II binding/core predictions for every candidate peptide through the official IEDB API, retaining the raw response, requested method, and candidate-to-response mapping.
2. Enumerate every possible 9-mer window for every peptide independently of the predictor, so register sensitivity can be inspected even when multiple binding cores remain plausible.
3. Construct a decoy-readiness table for each literature-annotated screen pair using only prespecified nuisance variables: paired peptide length, amino-acid composition, model confidence, and predicted-binding rank bin. Never use the pMHC priority score to select a decoy.
4. Deduplicate the external-validation overlay into evidence clusters, labeling the BALF5--MBP family as one calibration system and newer source annotations as contextual overlays rather than direct cross-reactivity positives.
5. Update the meeting brief with the actual predicted cores and results.

## Architecture

`src/premeeting_rigor.py` contains deterministic, dependency-free operations: core-window enumeration, IEDB TSV parsing, amino-acid composition distance, binding-rank binning, and decoy candidate ordering. `src/build_premeeting_rigor_artifacts.py` performs file and network I/O, writes all derived artifacts below `processed/`, and is rerunnable from the raw candidate manifest. The source response is preserved so the artifacts remain auditable if IEDB updates its models.

## Evidence boundaries

- An IEDB-predicted 9-mer core is a **computational register hypothesis**.
- A predicted binding rank is a matching covariate, not experimental HLA presentation evidence.
- A decoy table is an analysis design input, not a validated negative set.
- Cluster-aware overlay summaries test literature-annotation recovery only; they do not test shared TCR recognition, activation, affinity, or MS causation.

## Validation

- Unit tests must demonstrate inclusive 9-mer enumeration, robust parsing of a known IEDB response row, and decoy ordering that does not depend on the pMHC priority heuristic.
- The generator must fail if the API response cannot be mapped one-to-one to the manifest candidates.
- The generated summary must retain the IEDB method, endpoint, retrieval time, and raw response path.
