# pMHC Surface Electrostatics Analysis

This additive package evaluates 6 frozen sequence-supported candidates against their previously frozen 25-pair, exact-HLA N3 panels. It uses PDB2PQR 3.7.1 with PROPKA at pH 7.4 and PARSE parameters, followed by APBS 3.4.1 linearized Poisson--Boltzmann calculations under a shared panel grid.

## Result

Overall gate: `not_evaluable`.

The target ranks are reported in `target_electrostatic_summary.csv`. When the frozen V3 register-robustness flag is false, electrostatic ranks remain sensitivity evidence and the formal gate abstains. N3 pairs have unknown recognition status and were not treated as specificity negatives.

## Reproducibility

- `protocol_lock.json` freezes the method and claim boundary.
- `environment_manifest.json` records tool versions, archive checksum, and the container image digest.
- `protonation_provenance.csv` and `target_histidine_propka.csv` retain charge-state provenance.
- `apbs_provenance.csv` links every APBS input, log, transient grid checksum, and sampled vector.
- Full OpenDX files were deleted after checksum and deterministic point sampling to avoid a multi-gigabyte package; every input needed to regenerate them is retained.
- No existing V1--V3 package or discovery ranking was modified.
