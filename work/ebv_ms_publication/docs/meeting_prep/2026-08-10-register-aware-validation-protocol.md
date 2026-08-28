# Register-aware matched-decoy validation protocol

**Status:** Meeting-ready proposal; do not treat provisional register windows as biological annotations until reviewed by Yicong Liu and Olivia Thomas.

## Objective

Test whether the pMHC prioritization score retains useful discrimination after the comparison is restricted to plausible equivalent HLA-II P1--P9 registers and evaluated against biologically matched decoys.

## Pre-meeting implementation status

- IEDB `recommended_binding` HLA-DRB1*15:01 hypotheses were obtained for all 86 manifest peptides. The raw response, 30-mer handling for long peptides, and source-verified native flanks for 10-mers are retained in `processed/register_sensitivity/` and `raw/iedb_natural_flank_extensions.csv`.
- A descriptive same-register diagnostic was run on the existing 32-pair shortlist. It does not alter the ranking and retains every manifest-contained 9-mer window combination for sensitivity.
- The diagnostic shows that 28/30 assessable pairs have zero locally aligned residues at the same top-core P1--P9 position. This is a redesign signal, not a null biological result.
- The strict candidate universe supplies no complete five-decoy set. Expand the universe or change a prespecified tolerance only after mentor approval.

## Non-negotiable claim boundary

This protocol can evaluate pMHC-level prioritization only. Neither a passing nor failing result measures TCR affinity, T-cell activation, cross-reactivity, or patient pathogenicity.

## Analysis arms

| Arm | Inputs | Role | Interpretation |
|---|---|---|---|
| Structural calibration | Experimental BALF5--MBP pMHC structures | Verify core-surface metric on the established system | Calibration only; not independent validation |
| Register-aware screen | Predeclared 32 EBV--myelin pairs | Recalculate scores only for positions declared comparable by register | Candidate prioritization |
| Matched decoys | One or more decoys per evaluated pair | Test whether the score exceeds confounds that drive easy negatives | Primary discrimination benchmark |
| Independent literature systems | Distinct, known cross-reactive pMHC systems with suitable data | Test transfer beyond BALF5--MBP | External validation only if systems are genuinely independent |

## Register assignment hierarchy

For each peptide, use the first available source in this order and retain the source and confidence in the output table.

1. **Experimental pMHC structure:** derive the register from the resolved peptide position in the HLA groove.
2. **Published allele-specific experimental annotation:** use only when the exact peptide and class-II context match.
3. **Allele-specific prediction:** generate all plausible 9-mer cores and retain both the top call and near-tied alternatives. Record the predictor name, endpoint, retrieval date, input sequence, score, and version whenever the service exposes one.
4. **Motif-only provisional call:** permitted for the discussion pilot, never as a publication result without Yicong/Olivia approval.

If multiple registers remain plausible, perform a sensitivity analysis across all retained calls. Do not choose the register that gives the best similarity score after seeing the ranking.

The current IEDB top core is therefore an entry in this hierarchy, not the final register annotation.

## Register-aware scoring rule

1. Compare only residue pairs occupying the same declared register position.
2. Keep HLA-pocket-facing and likely TCR-exposed positions as separate summaries; do not collapse them without a predeclared rationale.
3. Preserve the current sequence, physicochemical, and HLA-fitted peptide-backbone components, but label the result `register_aware_pmhc_prioritization_score`.
4. Report the original score and register-aware score side-by-side. A rank change is an expected diagnostic, not a failure.
5. Exclude any pair with no defensible, comparable register rather than imputing equivalence.

## Matched-decoy construction

For each evaluated EBV or EBV--myelin pair, select decoys from the same predeclared candidate universe without inspecting the final similarity score.

Required matching variables:

- peptide length: exact where possible; otherwise within one residue;
- amino-acid composition: minimize composition distance rather than using an unrelated random sequence;
- predicted HLA-DRB1*15:01 binding: same ordinal affinity bin, defined before scoring;
- model confidence: same predeclared confidence bin for peptide residues;
- evidence tier and organism/antigen class: preserve the same input quality context as the evaluated peptide.

Recommended minimum: five decoys per evaluated positive/control when the candidate universe allows it. If fewer can be matched, report the shortfall and do not replace them with unmatched random peptides.

## Primary and secondary endpoints

**Primary endpoint:** predeclared rank/enrichment of the known control or literature-supported pair relative to its matched-decoy set.

**Secondary endpoints:**

- paired score difference between each control and its decoys;
- empirical permutation p-value using within-matching-stratum label swaps;
- sensitivity to alternative plausible registers;
- component-level contributions (sequence, physicochemical, backbone geometry) without interpreting them as causal mechanisms.

Do not use raw ColabFold ipTM, contact count, or a five-chain TCR prediction as a selective-binding endpoint; the Ob.1A12 calibration did not pass.

## Pre-specified decision rule for mentor approval

Proceed to a full register-aware benchmark only after the meeting confirms:

1. the register assignment source hierarchy;
2. BALF5--MBP's correct DR15-haplotype framing;
3. the matched-decoy variables and acceptable tolerance bins;
4. the primary endpoint and number of decoys;
5. the strongest claim permitted for a positive, null, or mixed result.

Until then, the pilot is a method-design artifact and not an analysis result.
