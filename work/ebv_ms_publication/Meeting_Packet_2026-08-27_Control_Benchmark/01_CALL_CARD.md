# Call card: held-out HLA-II positive-control benchmark

## Opening, about 60 seconds

> I replaced the single exposed control with three independent human HLA-II TCR systems and kept each biological system completely out while selecting weights. Across the eight panels that were formally complete, the known positive ranked in the top 3 every time. The worst held-out rank improved from 9 with exposed C-alpha geometry alone to 3 with the candidate composite.
>
> I am treating that as encouraging control recovery, not as a fully validated structural model. The formal gate is still not evaluable because two AlphaFold panels preserve missing jobs and the Hy.1B11 PDB panels have only two exact-HLA decoys each. Also, TCR-facing sequence identity alone recovered all eight completed controls with a worst rank of 2, so I have not shown that predicted 3D geometry adds value beyond sequence and chemistry. I left all discovery rankings unchanged.

## Three results to show

1. **Technical integrity:** all 48 new jobs matched the frozen requests exactly; all 240 new model samples passed exact-chain and clash checks.
2. **Held-out recovery:** the candidate composite achieved capture-at-3 in 8/8 completed panels. Completed ranks were 1, 1, 2, 2, 3, 3, 1, and 1.
3. **Baseline warning:** TCR-facing sequence identity also achieved 8/8, with a worst rank of 2. The structural pipeline has recovered the controls, but has not yet demonstrated independent structural signal.

## Exact status by biological system

### Hy.2E11, BALF5-MBP

- PDB oracle: complete, held-out composite rank 1 of 8
- AlphaFold seed 104759: complete, rank 1 of 26
- AlphaFold seed 104729: incomplete, available exposed-geometry rank 1 of 17
- System status: `not_evaluable`

### Ob.1A12, EngA-MBP

- PDB oracle: complete, held-out composite rank 2 of 6
- AlphaFold seed 104759: complete, rank 2 of 26
- AlphaFold seed 104729: incomplete, available exposed-geometry rank 2 of 21
- System status: `not_evaluable`

### Hy.1B11, UL15-MBP and PMM-MBP

- UL15 AlphaFold: rank 1 of 26 in both seeds
- PMM AlphaFold: rank 3 of 26 in both seeds
- Both PDB panels: not evaluable because each has only two exact-HLA decoys
- System status: `not_evaluable`

## The interpretation to defend

> The controls show that the provisional scoring procedure can recover established HLA-II cross-recognition systems in every completed held-out panel. They do not yet validate the full structural model, because required panels are incomplete and a simple sequence baseline performs at least as well.

## Decisions to request

1. Should the next preregistered benchmark require the composite to beat TCR-facing sequence identity, not only exposed C-alpha RMSD?
2. Is the strong weight on physicochemical mismatch biologically acceptable, or evidence that the benchmark is mostly recovering sequence resemblance?
3. Which same-TCR N1 or same-assay/HLA N2 negatives should be added for a real specificity test?
4. Should we expand the exact DQ structural decoy library, or predeclare an AlphaFold-only DQ validation layer?
5. For a future clean benchmark, should the incomplete AlphaFold panels be rebuilt from a new frozen manifest rather than patched into this one?

## Safe closing sentence

> My proposed next step is to keep the current discovery rankings frozen, strengthen the benchmark against the sequence-only baseline and specificity negatives, and only then regenerate rankings separately for each HLA.
