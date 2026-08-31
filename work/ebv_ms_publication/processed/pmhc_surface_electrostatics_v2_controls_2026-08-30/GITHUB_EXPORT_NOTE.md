# GitHub Results Export

This directory is the publishable results and provenance layer of the verified
control-first surface-electrostatics V2 package. The complete 6 GB archive is
retained in the authoritative iCloud research project at:

`processed/pmhc_surface_electrostatics_v2_controls_2026-08-30/`

GitHub includes the protocol and environment locks, registries, checksums,
rank and sensitivity tables, score matrix, gates, verification record, figure,
and executable source/tests. It omits the following reproducible heavy
intermediates:

- `aligned_models/`
- `model_calculation_records/`
- `raw_calculations/`
- `sampled_surface_vectors/`
- `io.mc`

`SHA256SUMS.csv` describes the complete authoritative archive, including those
omitted intermediates. The locked control gate is `fail`, candidate evaluation
is not allowed, and electrostatics is retired from candidate ranking.
