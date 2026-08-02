# Publication strategy

## Working title

Calibrated computational prioritization of EBV--CNS molecular mimicry candidates in the HLA-DR15 antigen-presentation pathway

## Central thesis

This paper should present a reproducible computational prioritization workflow,
not a new discovery of the BALF5--MBP mimicry pair and not a standalone proof
of patient TCR cross-reactivity.

The strongest defensible claim is:

> We built an auditable HLA-DR15-focused workflow that uses experimentally
> established EBV--MBP pMHC mimicry as a positive-control anchor, validates
> modeled pMHC inputs, ranks EBV--CNS candidate pairs by local pMHC surface
> similarity, and explicitly rejects uncalibrated ternary predictions as
> TCR-recognition evidence.

## What makes the paper publishable

1. It corrects the biological mechanism: BALF5--MBP mimicry is an allele-pair
   switching system across DR2a and DR2b, not a same-allele DRB1*15:01 result.
2. It is reproducible and auditable: every included peptide keeps provenance,
   sequence, source antigen, HLA restriction, and QA status.
3. It has a real experimental anchor: PDB 1H15 and 1BX2 provide the known
   pMHC positive-control comparison.
4. It reports negative calibration honestly: the Ob.1A12 ternary ColabFold run
   failed its positive-control orientation benchmark.
5. It connects to current mechanism: recent literature supports EBV-driven
   changes in DR15 antigen presentation and DRB1*15:01-restricted EBV CD4 T-cell
   responses.

## Paper structure

### Abstract

State the need for transparent prioritization of EBV--CNS molecular mimicry
candidates in HLA-DR15. Report the experimental positive-control calibration,
the 86-model pMHC QA set, the top BALF5-family screen result, and the failure
of uncalibrated ternary ColabFold to support TCR-recognition claims.

### Introduction

Frame EBV and HLA-DR15 as convergent MS risk factors. Introduce class-II antigen
presentation and CD4 T cells. Then narrow to the specific methodological gap:
computational mimicry studies often overstate peptide similarity or docking
scores without calibration against known positive and negative controls.

### Methods

Describe dataset curation, inclusion/exclusion, pMHC model QA, experimental
positive-control alignment, candidate ranking, null/negative controls, GEO
context analysis, and ternary-model rejection criteria.

### Results

1. Experimental BALF5--MBP structures define the positive-control geometry.
2. The curated 86-model pMHC set passes peptide and chain-layout QA.
3. Corrected local pMHC screening prioritizes BALF5-family EBV--myelin pairs.
4. The Ob.1A12 ternary workflow fails calibration and is not used as binding
   evidence.
5. Public expression context is negative or contextual, not direct validation.

### Discussion

Emphasize that the workflow prioritizes hypotheses rather than proving
pathogenic cross-reactivity. The paper's contribution is methodological rigor:
it shows what can be inferred from pMHC modeling, what cannot, and what
experimental data would be needed next.

## Figure plan

| Figure | Purpose | Current artifact |
|---|---|---|
| Figure 1 | Experimental positive-control pMHC equivalence | `processed/experimental_positive_control/figure_1_experimental_positive_control_300dpi.png` |
| Figure 2 | Full pMHC screen ranking with top candidates and controls | `processed/publication_figures/figure_2_pmhc_screen_shortlist_300dpi.png` |
| Figure 3 | Claim/evidence ladder showing permitted inference levels | `processed/publication_figures/figure_3_claim_ladder_300dpi.png` |
| Figure 4 or supplement | PyTorch transcriptomic-context classifier for GSE190847 | `processed/publication_figures/figure_4_gse190847_pytorch_expression_300dpi.png` |
| Figure 5 | External validation overlay against literature-positive candidates | `processed/publication_figures/figure_5_external_validation_overlay_300dpi.png` |
| Supplement 1 | pMHC structure QA summary | `processed/colabfold_pmhc_peptide_qa.csv` |
| Supplement 2 | Ternary calibration failure | `processed/ob1a12_ternary_evaluation/FINAL_CONCLUSION.md` |
| Supplement 3 | GEO contextual negative result | `processed/geo_gse190847/` |

## Submission path

The first realistic target should be a computational biology, immunology
methods, or student-accessible preprint-style venue. The manuscript should not
be framed as a clinical MS mechanism paper unless a mentor helps add direct
experimental validation.

## Decisive remaining work

1. Decide whether the PyTorch expression classifier stays as Figure 4 or moves
   to the supplement.
2. Write the manuscript from the skeleton without adding unsupported claims.
3. Add a reproducibility appendix listing exact input files and scripts.
4. Ask a mentor to review the claim matrix before submission.
5. If Hy.2E11 full receptor sequences become available, run a new calibrated
   receptor-level benchmark as a separate revision, not as an inferred result.
