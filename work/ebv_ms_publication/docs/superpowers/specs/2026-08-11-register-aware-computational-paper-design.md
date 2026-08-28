# Register-Aware Computational Paper Design

## Objective

Produce a reproducible computational pMHC candidate-prioritization paper for
EBV and myelin peptides restricted to HLA-DRB1*15:01. The paper must be
scientifically useful whether the final benchmark identifies robust candidate
pairs or shows that the current evidence cannot support them.

## Permitted conclusion

The strongest positive conclusion is: register-aware pMHC features prioritize
specific EBV-myelin pairs relative to score-blind matched decoys under the
defined computational model. These are candidates for experimental testing.

The project will not claim shared-TCR binding, T-cell activation,
cross-reactivity, molecular mimicry as an established biological fact, or an
EBV-driven MS mechanism.

## Evidence hierarchy

1. Exact same-allele experimental register evidence overrides a predictor.
   MBP(85-99) `ENPVVHFFKNIVTPR` uses its PDB 1BX2 DRB1*15:01 reference
   register `VHFFKNIVT`.
2. Different-allele structural evidence is calibration only. BALF5 PDB 1H15
   remains DRB5*01:01/DR2a calibration and cannot validate a DRB1*15:01
   same-allele result.
3. Exact peptides without an established DRB1*15:01 register remain
   sensitivity-only. This currently includes gH `EKQLFYYIGTMLPN` and MBP
   `QRPGFGYGGRASDYKSAHK`.
4. Other candidates use the retained IEDB top-core prediction only as a
   computational hypothesis, with its source and sensitivity status exposed in
   every result table.

## Analysis design

### 1. Freeze inputs and eligibility

Use the source-traceable candidate manifest, geometry matrix, IEDB prediction
record, and experimental-register override registry. A primary-analysis pair
must have two primary-allele, manifest-contained cores and at least one
previously aligned residue at the same P1-P9 position. Exclude calibration-only
and sensitivity-only records from primary score and decoy eligibility while
keeping them in auditable output tables.

### 2. Register-aware feature scoring

Recalculate each component from the pre-existing local alignment after
filtering to same-register positions. Keep three outputs separate:

- all same-register aligned positions;
- anchor-focused positions P1/P4/P6/P9;
- non-anchor, candidate-exposed positions P2/P3/P5/P7/P8.

The anchor/exposure labels are descriptive structural annotations, not TCR
contact calls. Do not select favorable alternative cores after looking at the
score; report retained register sensitivity separately.

### 3. Matched-decoy benchmark

For each externally annotated/control target, select five decoys from the full
pre-score PASS universe only when both peptide lengths are within one residue
and both IEDB binding-rank bins match. Use composition distance and peptide
pLDDT only to order otherwise eligible decoys. Never use the priority score,
rank, RMSD, or register-aware score to select them.

Primary endpoint: target rank/enrichment relative to its strict decoy set.
Secondary endpoints: paired score difference, within-stratum permutation test,
and score sensitivity across retained ambiguous registers.

### 4. Result decision gate

**Positive computational result:** one or more non-control EBV-myelin pairs
remain same-register eligible, have complete strict decoy sets, and show a
predeclared robust score/enrichment pattern that persists under register
sensitivity analysis.

**Negative or mixed computational result:** no such pair survives the complete
eligibility and decoy criteria. Report this as a method result demonstrating
why unregistered local similarity is insufficient for molecular-mimicry
prioritization.

Either outcome is useful; neither establishes biological cross-reactivity.

## Deliverables

1. Frozen input/override registry and data dictionary.
2. Register-aware score table with component and position-class columns.
3. Strict-decoy and permutation result tables with explicit shortfalls.
4. Sensitivity appendix for alternative plausible registers.
5. One figure showing the evidence hierarchy and one figure showing primary
   versus sensitivity-only outcomes.
6. Manuscript results text for both positive and negative decision branches.

## Prospective wet-lab direction

Wet-lab work is explicitly outside this paper. The future validation section
will propose peptide-HLA binding/elution assays, register-resolving
competition/mutagenesis, and clone-defined TCR assays with HLA-restriction
controls for the final computational candidates.
