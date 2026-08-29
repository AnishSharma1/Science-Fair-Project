# Eight-candidate computational evidence dossier

This additive package audits the eight sequence-supported high-yield HLA-II candidates. It does not rerank the discovery universe or create a composite score.

## Result state

- Stage-one gate: `complete`
- High priority: `0`
- Medium priority: `2`
- Hold: `6`
- Stage two: `not_evaluable_pending_experimental_binding_and_register`

## Evidence layers

- IEDB assay provenance is classified by exact sequence, exact HLA, MHC class, and host.
- NetMHCIIpan 4.3 EL/BA and MixMHC2pred 2.1 remain separate predictors.
- HLA Ligand Atlas release 2020.12 hits distinguish exact, nested, monoallelic, and multiallelic evidence.
- All eight submitter-provided PXD068488 DR2a/DR2b processed peptide tables were searched; raw spectra were not reprocessed.
- gnomAD r4 canonical-transcript missense variants were checked with a locked common-frequency threshold of 1%.
- Human rarity database scope: `full_cached_public_source`.
- EBV reciprocal rarity scope: `full_cached_public_source`.
- Immunopeptidome absence is missing evidence, never a negative.

## Files

- `protocol_lock.json` and `source_manifest.json`: frozen rules, versions, sources, and checksums.
- `raw_response_manifest.csv`: checksum and exact target linkage for every cached raw file.
- `peptide_arm_evidence.csv`: sixteen peptide-arm evidence records.
- `candidate_evidence_matrix.csv`: eight candidate-level records.
- `iedb_assay_provenance.csv`: exact assay classifications with raw-response linkage.
- `predictor_register_comparison.csv`: independent predictor scores and core agreement.
- `immunopeptidome_hits.csv`: observed exact/nested ligand records.
- `proteome_rarity_summary.csv` and `proteome_nearest_neighbors.csv`: unconditioned sequence rarity diagnostics.
- `presentation_conditioned_rarity.csv`: separate exact-HLA candidate-library rarity diagnostic.
- `conservation_results.csv`: sequence conservation plus mapped common overlapping missense variants.
- `stage1_assay_recommendations.csv`, `stage1_assay_gate.json`, and `stage2_tcell_gate.json`: experimental funnel.
- `candidate_dossiers/`: human-readable evidence sheets.
- `SHA256SUMS.csv`: deterministic checksums.

## Claim boundary

This dossier prioritizes peptide-HLA binding and register experiments. It does not establish natural presentation, TCR recognition, specificity, cross-reactivity, molecular mimicry, probability, false-discovery rate, or an MS mechanism.
