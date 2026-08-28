# Yicong call card: EBV--MS pMHC register redesign

## Opening (about 75 seconds)

> I built this as a reproducible HLA class-II pMHC candidate-prioritization workflow, not as a claim of TCR cross-reactivity. We used the established BALF5--MBP system as a structural calibration anchor, then ranked 32 predeclared EBV--myelin pairs from local sequence, physicochemical, and groove-fitted peptide-backbone features.
>
> The key pre-meeting result is a limitation that changes the next step. I added IEDB DRB1*15:01 core hypotheses for all 86 candidates and checked whether residues the original score aligned occupy the same P1--P9 position. Of 30 assessable pairs, 28 retain no same-register local alignments; the three BALF5--MBP-like shortlist hits retain zero of six. So I am not calling those register-equivalent or TCR-facing mimics.
>
> I also tested strict matched-decoy feasibility. In the current 32-pair universe, no labeled target has five decoys matched on peptide length and binding-rank bin. I want to agree on the register hierarchy, the proper DR15 allele framing for BALF5--MBP, and whether to expand the candidate universe or predeclare a single tolerance change before rerunning anything.

## Three results to show

1. **Structural calibration:** BALF5--MBP experimental pMHCs have a 0.838 A seven-position core CA RMSD after HLA-groove fitting. This calibrates the metric; it does not validate a new TCR interaction.
2. **Register diagnostic:** all 86 candidates have auditable IEDB DRB1*15:01 hypotheses. Of 32 pairs, 30 are assessable; 28/30 retain 0 same-register local alignments under top-core hypotheses. Rank 4 retains 5/5, but is only a source/context overlay and needs review.
3. **Decoy feasibility:** 14/16 labeled pairs have zero strictly matched decoys; two have partial sets totaling three decoys. No target has the predeclared five.

## Ask for decisions, not reassurance

1. What register assignment hierarchy should govern the next analysis: experimental structure, exact published annotation, predictor consensus, and then sensitivity for remaining ambiguity?
2. Are BALF5/DRB5*01:01 and MBP/DRB1*15:01 legitimate only as a DR15-haplotype structural calibration, rather than same-allele validation?
3. Which P1--P9 positions should be summarized as HLA-pocket-facing versus candidate TCR-exposed in this class-II context, and which comparisons are biologically defensible?
4. Should the decoy redesign expand the pre-scoring candidate universe, or relax exactly one matching tolerance? Which primary endpoint should be locked before scoring?

## If asked directly about TCR evidence

> We do not have it. The five-chain Ob.1A12 ColabFold calibration did not recover the experimental TCR pose, so no predicted TCR contacts, affinity, activation, or cross-reactivity claims are being carried forward.

## Screens/files to open

- `docs/meeting_prep/2026-08-10-premeeting-rigor-results.md`
- `processed/register_sensitivity/register_aware_shortlist_diagnostic.csv`
- `processed/register_sensitivity/register_window_pair_sensitivity.csv`
- `processed/matched_decoys/decoy_feasibility_summary.csv`
- `processed/register_sensitivity/positive_control_allele_context.csv`
- `docs/meeting_prep/2026-08-10-revised-claim-language.md`

## Safe closing sentence

> My proposed next deliverable is a preregistered-in-spirit, register-aware pMHC prioritization benchmark with mentor-approved registers and matched decoys; it would still be a candidate-prioritization result, not a TCR or disease-mechanism result.
