# Yicong Liu (Robbin) meeting brief

**Meeting:** Monday, August 10, 2026, 9:30 AM Central / 4:30 PM Sweden  
**Format:** Zoom invitation in the email thread  
**Likely participants:** Yicong Liu; Olivia Thomas may join if available  
**Purpose:** Agree on a biologically defensible next validation step for the EBV--HLA-DR15 pMHC prioritization workflow.

## The 90-second opening

The project is a reproducible, HLA class-II pMHC candidate-prioritization workflow. It uses a known EBV BALF5--MBP molecular-mimicry system as a structural calibration anchor, performs chain and peptide QA on 86 modeled pMHC structures, and ranks 32 predeclared EBV--myelin pairs by local sequence, physicochemical, and HLA-fitted peptide-backbone features.

The claim boundary is deliberately narrow: the current score identifies local similarity in modeled pMHCs. It does **not** establish equivalent HLA-binding registers, a common TCR-accessible surface, TCR binding, T-cell activation, or patient mechanism. A five-chain TCR-pMHC ColabFold experiment failed recovery of the Ob.1A12 positive-control pose and is excluded from evidence.

Before this meeting, I added IEDB HLA-DRB1*15:01 top-core hypotheses for all 86 candidates, a 32-pair same-register diagnostic, a strict matched-decoy feasibility screen, and an independence audit. The main result is cautionary: the top three BALF5--MBP-like screen pairs retain zero same-register local alignments under their top-core hypotheses, and no labeled target has five strictly matched decoys in the current universe. I would like advice on the correct redesign rather than a retrospective justification.

## Current evidence and its correct interpretation

| Finding | What it supports | What it does not support |
|---|---|---|
| Experimental BALF5--MBP pMHC structures have a 0.838 A seven-position core CA RMSD after HLA-groove fitting | The structural metric recovers a known, literature-established mimicry anchor | New evidence of Hy.2E11 affinity, activation, or cross-reactivity |
| 86 rank-1 pMHC models passed structural QA | The screen's modeled inputs have the expected chain/peptide layout | Biological presentation or TCR recognition |
| 32 EBV--myelin pairs were ranked before literature overlay | A reproducible prioritization set | A discovery-validation cohort |
| Classic BALF5--MBP records recover strongly in the overlay benchmark | Positive-control calibration/literature-overlay recovery | Independent external validation, because records overlap the established system |
| Strict new-literature overlay: 2 of 5 positives in top 10; empirical p = 0.5036 | A transparent, non-confirmatory independent overlay result | Generalizable validation or evidence of enrichment |
| All 86 candidates now have auditable IEDB DRB1*15:01 core/binding hypotheses | A reproducible starting point for register review | Experimental presentation or a final register call |
| 28 of 30 assessable shortlist pairs have zero same-register local alignments under IEDB top cores | The current local score is not register-aware | That an alternative window should be selected post hoc |
| No labeled pair has five strict matched decoys in the 32-pair universe | A real design constraint exposed before analysis | Permission to use unmatched random controls |

## Decisions requested from Yicong and Olivia

1. **Register assignment:** Is the IEDB top-core call, retaining all window sensitivity, an acceptable preliminary hierarchy? If not, what experimental/predictor consensus rule should replace it?
2. **Allele logic:** Is BALF5--MBP best presented strictly as a DR15-haplotype molecular-mimicry calibration anchor, given BALF5 is presented by DRB5*01:01 and MBP by DRB1*15:01? What comparisons are legitimate across those two class-II molecules?
3. **Comparison rule:** Which positions must match by register before a local structural score can be called biologically comparable? How should likely HLA-pocket-facing and TCR-exposed positions be treated differently?
4. **Decoy design:** Should the next step expand the pre-scoring candidate universe or relax one specific matching tolerance? The present 32-pair universe gives no five-decoy set under strict length and binding-bin matching.
5. **Evaluation endpoint:** Is rank enrichment against matched decoys the primary endpoint? Which secondary endpoints would be worth reporting?
6. **Claim boundary:** If a register-aware score separates positive controls from matched decoys, what exact wording would be acceptable? If it does not, should the manuscript remain a pMHC prioritization study only?

## Screens to have open

1. `processed/experimental_positive_control/figure_1_experimental_positive_control_300dpi.png`
2. `processed/publication_figures/figure_2_pmhc_screen_shortlist_300dpi.png`
3. `processed/publication_figures/figure_3_claim_ladder_300dpi.png`
4. `docs/meeting_prep/2026-08-10-premeeting-rigor-results.md`
5. `processed/register_sensitivity/register_aware_shortlist_diagnostic.csv`
6. `processed/matched_decoys/decoy_feasibility_summary.csv`
7. `processed/validation_hygiene/validation_evidence_clusters.csv`
8. `processed/register_sensitivity/positive_control_allele_context.csv`

## Call flow

| Time | Topic | Desired output |
|---|---|---|
| 0:00--0:02 | Opening and narrow claim boundary | Shared understanding of what the present data do and do not show |
| 0:02--0:06 | BALF5--MBP calibration and allele context | Confirm correct use of the positive control |
| 0:06--0:12 | Register diagnostic | Mentor-approved way to identify comparable P1--P9 positions |
| 0:12--0:20 | Decoy shortfall and protocol | Exact matching variables, candidate universe, and primary endpoint |
| 0:20--0:25 | Olivia's HLA/mimicry read | Biological constraints and manuscript framing |
| Last 2 min | Read back decisions | Concrete next analysis and allowed claim |

## Pre-call checklist

- [ ] Read this brief aloud once; keep the opening under 90 seconds.
- [ ] Open all eight artifacts above in separate tabs/windows.
- [ ] Join Zoom by 9:20 AM Central with headphones and screen sharing tested.
