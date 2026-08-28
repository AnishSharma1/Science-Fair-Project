# Same-register AF3 publication-readiness decision

## Direct conclusion

The dated lead-focused audit supports **one strict-control computational pMHC
lead**:

- EBV: `EBV_TCELL_950`, EBNA-1 peptide `AEGLRALLARSHVER`
- Human: `HUMAN_MYELIN_112214`, MBP peptide `LSRFSWGAEGQRPGFGYGG`
- Predicted P1--P9 cores: `LRALLARSH` versus `WGAEGQRPG`
- HLA context recorded for both inputs: HLA-DRB1*15:01

Rank 1 is the only strict-control lead. Rank 2
(`EBV_TCELL_2268741::HUMAN_MYELIN_117032`) is a separate,
length-sensitivity/job-dependence result and is not a co-primary lead. The
paper can be framed as a register-aware computational pMHC prioritization
study that nominates a testable EBNA-1--MBP structural-resemblance hypothesis.
It does not establish peptide presentation, TCR binding, activation,
cross-reactivity, molecular mimicry, or an MS mechanism.

## Dated lead-focused structural-control result (2026-08-15)

The frozen discovery order was preserved. The audit analyzed the two leads in
separate layers and did not pool or average them.

### Rank 1: primary strict-control panel

`EBV_TCELL_950::HUMAN_MYELIN_112214` was classified
`consistent_positive` using only its three strict primary controls.

- target median: 0.643 A
- equal-weight background median: 7.964 A
- background-minus-target delta: 7.321 A
- leave-one-control-out delta range: 5.261 to 10.086 A
- exploratory target rank: 1 of 4; empirical tail fraction: 0.25
- 10,000-iteration technical-stability interval: 1.620 to 12.759 A

The target job-pair medians were 0.427--0.864 A and were below the
equal-weight background median in all four target job-pairs. The largest pose
cluster contained 23 models from 7 jobs. These are pose/job consistency
annotations for saved AlphaFold technical samples, not biological replicates.

### Rank 2: supplemental length sensitivity/job dependence

`EBV_TCELL_2268741::HUMAN_MYELIN_117032` used only the deliberate
`length_sensitivity_exact_bin_pm7` layer and was classified
`length_sensitivity_only__mixed_positive`.

- target median: 6.396 A
- equal-weight background median: 13.975 A
- background-minus-target delta: 7.578 A
- leave-one-control-out delta range: 5.303 to 9.944 A
- exploratory target rank: 1 of 4; empirical tail fraction: 0.25
- 10,000-iteration technical-stability interval: -0.517 to 16.880 A

Its target job-pair medians ranged from 0.465 to 12.328 A and its technical
interval crosses zero. It therefore remains sensitivity-only and must not be
combined with, or given the evidentiary weight of, rank 1.

The empirical tail fractions are exploratory ranks, not p-values. The
bootstrap intervals quantify technical stability across saved jobs/models only
and are not biological replication. The full audit is documented in
`../lead_focused_robustness_2026-08-15/LEAD_FOCUSED_FINDINGS.md`,
`control_rank_and_leave_one_out.csv`, `technical_bootstrap_summary.csv`,
`job_pair_stability.csv`, and `pose_cluster_membership.csv` in that directory.

## Why the result is not molecular-mimicry proof

1. **Computational pMHC geometry.** AlphaFold pMHC models are descriptive
   structure hypotheses, not experimental evidence of peptide presentation.
2. **Computational registers.** Both core assignments are unique IEDB
   predictions, not exact experimentally resolved registers for this pair.
3. **No receptor evidence.** pMHC geometry does not establish shared-TCR
   binding, T-cell activation, or cross-reactivity.
4. **Technical, not biological, replication.** The saved AlphaFold jobs and
   models are technical samples, not patients, experiments, or independent
   biological systems.
5. **Limited frozen controls.** Three controls set the empirical resolution;
   an empirical tail fraction of 0.25 is the smallest possible fraction here.
6. **The myelin evidence is aggregate.** The EBV peptide has a
   DRB1*15:01 T-cell-assay record, while the MBP peptide is a Tier-3 aggregate
   DRB1*15:01 epitope record rather than direct cross-reactivity evidence.

The classic published BALF5--MBP system reached a higher evidentiary bar: the
same patient-derived TCR recognized both pMHCs, and experimental crystal
structures supported surface equivalence. That system is an external
calibration precedent, not evidence for the new EBNA-1--MBP pair.

## Figure placement

- **Main results:** rank-1 primary-control panel in the main results.
- **Supplement:** rank-2 sensitivity/job-dependence and PyTorch classifier in
  the supplement.

## Defensible paper claim

> A predeclared, register-aware computational analysis retained one
> HLA-DRB1*15:01 EBNA-1--MBP pMHC pair as a strict-control lead, with a
> 7.321-A background-minus-target median difference and an exploratory
> empirical tail fraction of 0.25 across three frozen controls. A separate
> rank-2 length-sensitivity analysis was mixed across saved jobs. These
> results nominate a testable pMHC structural-resemblance hypothesis but do
> not establish TCR cross-reactivity or molecular mimicry.

## Publication decision

- **Suitable now:** computational methods/resource paper, hypothesis paper, or
  preprint centered on same-register prioritization, transparent primary and
  sensitivity-control results, and the single strict-control EBNA-1--MBP lead.
- **Not suitable now:** a paper claiming discovery or proof of EBV--myelin
  molecular mimicry, TCR cross-reactivity, or an MS mechanism.
- **Submission-readiness limitation:** report the three-control empirical
  resolution, rank-2 length mismatch/job dependence, and technical-sample
  boundary directly; do not replace them with unmatched controls or favorable
  post-hoc windows.

## Literature calibration

- Lang et al., *Nature Immunology* (2002), [A functional and structural basis
  for TCR cross-reactivity in multiple sclerosis](https://www.nature.com/articles/ni835).
- Abramson et al., *Nature* (2024), [Accurate structure prediction of
  biomolecular interactions with AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w).
- Terwilliger et al., *Nature Methods* (2024), [AlphaFold predictions are
  valuable hypotheses and accelerate but do not replace experimental structure
  determination](https://www.nature.com/articles/s41592-023-02087-4).
