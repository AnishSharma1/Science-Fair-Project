# Pre-meeting rigor results

**Status:** Analysis-design evidence for the August 10 meeting. These outputs are computational hypotheses and bookkeeping checks, not evidence of peptide presentation, a shared TCR surface, T-cell activation, cross-reactivity, or disease mechanism.

## What is now complete

| Deliverable | Verified result | Meeting implication |
|---|---|---|
| HLA-II binding/register coverage | All 86 manifest candidates have an IEDB `recommended_binding` result for HLA-DRB1*15:01: 76 direct inputs, 7 long peptides tiled with overlapping 30-mers, and 3 short peptides evaluated with source-coordinate-verified natural flanks | The screen no longer ignores length-ineligible candidates, but IEDB calls remain provisional register hypotheses |
| Manifest-contained top cores | 84/86 top cores are fully contained in their manifest peptide; 2/86 extend beyond a 10-mer into verified native flank | Do not force a P1--P9 comparison for the two flank-dependent cases |
| Register-position diagnostic | Of 32 shortlisted pairs, 30 are assessable from a unique manifest-contained IEDB top core; 28/30 retain **zero** of their existing local alignments at the same P1--P9 position, 1 retains 2, and 1 retains 5 | The current ranking is not yet a register-aware result; this is a clear reason to redesign rather than overclaim |
| Top BALF5--MBP-like screen hits | Ranks 1--3 each have 6 original local alignments but 0 same-register alignments under their IEDB top-core hypotheses | Keep them as structural/modeling candidates, not register-equivalent mimicry candidates |
| Strict matched-decoy feasibility | For 16 literature/context-labeled pairs, 14 have no eligible decoy and 2 have only partial sets (3 decoys total). No pair has the prespecified five decoys | The present 32-pair shortlist is too narrow for a five-decoy benchmark under strict length and binding-bin matching |
| Validation independence audit | BALF5--MBP collapses to one calibration system; Drosu and Wang annotations are source/context overlays, not direct cross-reactive pair validation | The external-overlay language has been corrected before it reaches the meeting |

## Register diagnostic: the result to lead with

The diagnostic retains the original 32-pair ordering and does **not** rerank it. For each already locally aligned residue pair, it asks whether the two residues would occupy the same P1--P9 position under the IEDB top-core call. It also writes all 3,714 pairs of manifest-contained 9-mer windows, so an apparently favorable alternative cannot be selected after inspecting the answer.

- Rank 1: BALF5-like EBV `TGGVYHFVKKHVHES` vs MBP-region `VVHFFKNIVTPRT` has 6 original local alignments but **0/6** same-register alignments under top cores `VYHFVKKHV` (start 4) and `VHFFKNIVT` (start 2).
- Ranks 2 and 3 show the same top-core result: **0/6**.
- Rank 4, `EBV_TCELL_2268683::HUMAN_MYELIN_115641`, retains **5/5** under its top-core hypothesis. Its panel annotation is only source/context overlay, so it is a review candidate, not a validated molecular-mimicry positive.
- The result is sensitive to possible windows for many pairs. That sensitivity is recorded, not optimized: it must be resolved with a predeclared hierarchy approved by Yicong/Olivia.

## Decoy decision that cannot be deferred

The only strict prespecified decoy rule used here was:

- both peptide lengths within one residue; and
- zero mismatch in the IEDB percentile-rank bin (strong <=2, intermediate <=10, weak >10).

Even before adding composition and model-confidence constraints to the eligibility filter, the current screening universe cannot deliver five decoys per labeled target. The decision for the meeting is therefore binary:

1. Expand the **pre-scoring** HLA-DRB1*15:01 candidate universe and rerun the screen; or
2. Change a matching tolerance in advance, state the trade-off, and keep it fixed for every target.

Do not replace the shortfall with random decoys.

## Files to open in the meeting

1. `processed/register_sensitivity/register_prediction_summary.csv` — one audited prediction record per manifest candidate.
2. `processed/register_sensitivity/register_aware_shortlist_diagnostic.csv` — the 32-pair top-core diagnostic.
3. `processed/register_sensitivity/register_window_pair_sensitivity.csv` — all retained 9-mer window combinations.
4. `processed/matched_decoys/decoy_feasibility_summary.csv` — the decoy shortfall.
5. `processed/validation_hygiene/validation_evidence_clusters.csv` — independence and claim boundaries.
6. `processed/register_sensitivity/positive_control_allele_context.csv` — BALF5/MBP DRB5 versus DRB1 computational context.

## Exact decisions to obtain

1. Is the IEDB top-core call, with full-window sensitivity retained, an acceptable preliminary register hierarchy? If not, what predictor/experimental consensus rule replaces it?
2. May the BALF5--MBP pair be used only as DR15-haplotype structural calibration because BALF5 and MBP have different presenting alleles in the established system?
3. Should the next iteration expand the pre-scoring candidate universe or relax a specific decoy tolerance? Which one and why?
4. Which register positions can be summarized as pocket-facing versus candidate TCR-facing for this allele context, if any?
5. Does the rank-4 same-register result merit experimental/literature review, or should it remain solely a contextual overlay?

## Reproducibility

- IEDB raw response: `processed/register_sensitivity/iedb_mhcii_drb1501_raw.tsv`
- Flank provenance: `raw/iedb_natural_flank_extensions.csv`
- Generation scripts: `src/build_premeeting_rigor_artifacts.py`, `src/build_register_aware_pair_diagnostic.py`, and `src/build_positive_control_allele_context.py`
- Tests: `tests/test_premeeting_rigor.py` and `tests/test_register_aware_diagnostic.py`
