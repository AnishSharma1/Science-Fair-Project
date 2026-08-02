# Manuscript skeleton

## Title

Calibrated computational prioritization of EBV--CNS molecular mimicry candidates in the HLA-DR15 antigen-presentation pathway

## Abstract draft

Epstein--Barr virus (EBV) infection and the HLA-DR15 haplotype are major
converging risk factors for multiple sclerosis (MS), but computational
molecular-mimicry studies can overstate peptide or docking similarities when
not calibrated against experimental controls. Here, we rebuilt an auditable
HLA class-II workflow for prioritizing EBV--CNS candidate mimicry pairs. We
used the established Hy.2E11-associated BALF5--MBP system as an experimental
positive-control anchor, correctly treating it as allele-pair switching between
DRB5*01:01 and DRB1*15:01 rather than a same-allele DRB1*15:01 result.
Experimental pMHC structures 1H15 and 1BX2 showed close seven-position peptide
core equivalence after HLA-groove superposition. We then applied predefined QA
to 86 rank-1 modeled pMHC structures and ranked EBV--myelin candidates by local
peptide geometry and physicochemical similarity. The corrected screen
prioritized BALF5-family candidates, but an Ob.1A12 five-chain ternary
ColabFold workflow failed its known positive-control orientation benchmark and
was therefore not used as evidence of TCR recognition. This study provides a
transparent prioritization framework and evidence ledger for EBV--CNS pMHC
mimicry hypotheses while defining the experimental data required for
cross-reactivity claims.

## Claim boundary

This manuscript can claim a calibrated computational prioritization workflow.
It cannot claim new discovery of BALF5--MBP mimicry, TCR binding, T-cell
activation, patient pathogenicity, or causal proof of MS mechanism.

## Results text scaffold

### Result 1: Experimental positive-control calibration

The BALF5--MBP system was first recoded as an experimental positive-control
anchor. The correct biological relationship involves BALF5 presented by
HLA-DRB5*01:01 and MBP presented by HLA-DRB1*15:01. This distinction is central
because replacing DRB5*01:01 with DRB1*15:01 changes the mechanism being
tested.

### Result 2: pMHC model quality gate

Before structural comparison, modeled pMHC files were required to contain the
expected peptide and a DRA-like/DRB-like/peptide chain layout. This gate
removed earlier legacy structures with duplicated DRB-like chains or ambiguous
chain layouts from publication evidence.

### Result 3: Candidate prioritization screen

The corrected screen compared Tier-1 EBV T-cell peptides against myelin/CNS
candidates using local peptide geometry after HLA-groove normalization and a
physicochemical similarity score. BALF5-family comparisons ranked highly, but
the output is interpreted as candidate prioritization only.

### Result 4: Ternary model calibration failure

The Ob.1A12 ternary ColabFold screen was tested against the known MBP positive
control. Because the predicted positive-control TCR orientation did not recover
the experimental 1YMM pose, EBV ternary scores were rejected as evidence of
selective recognition.

### Result 5: External literature validation overlay

The ranked pMHC shortlist was overlaid with independently reported
literature-positive candidates. Classic BALF5--MBP positive-control pairs were
strongly enriched near the top of the shortlist, supporting the workflow as a
structural-mimicry prioritization tool. In contrast, newer DRB1*15:01 EBV
glycoprotein positives from the antigen-presentation literature were not
broadly enriched. This difference should be interpreted as a boundary of the
metric: it recovers the classic pMHC surface-mimicry signal better than it
captures all EBV--MS mechanisms.

### Result 6: Contextual expression evidence and PyTorch linear probe

The GSE190847 B-cell expression analysis provided contextual evidence only. It
did not directly measure EBV infection, pMHC presentation, or TCR recognition,
and its HLA-II/APC module was not significant after correction.

A bounded PyTorch linear-probe analysis was added to test whether the same
pre-registered seven-gene HLA-II/APC panel carried out-of-sample signal for
untreated PPMS versus healthy controls. This classifier is useful only as
transcriptomic context. It should be reported with leave-one-out performance,
permutation testing, and the explicit caveat that it does not validate EBV
infection, antigen presentation, or molecular mimicry.

## Discussion scaffold

The main contribution is not a new biological proof but a cleaner evidentiary
standard. A publishable computational mimicry paper must separate four layers:
epidemiologic/clinical context, experimental pMHC anchors, computational
candidate ranking, and direct receptor or activation validation. This project
supports the first three layers and explicitly identifies the fourth as future
work.

## Methods checklist

- Dataset provenance and IEDB inclusion/exclusion rules.
- Protein/accession validation.
- pMHC ColabFold input generation.
- Peptide and chain-layout QA.
- Experimental pMHC superposition metric.
- Candidate ranking score and caveats.
- Negative controls and null models.
- External validation overlay and rank-recovery permutation test.
- Ternary calibration benchmark and rejection rule.
- GEO contextual analysis.
- PyTorch expression classifier: feature panel, leave-one-out protocol,
  permutation test, and limitations.
- Reproducibility and file manifest.
