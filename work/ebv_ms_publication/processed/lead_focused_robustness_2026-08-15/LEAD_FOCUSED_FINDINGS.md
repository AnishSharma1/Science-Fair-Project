# Lead-focused structural robustness findings

## Scope and methods

This audit keeps the frozen discovery ranking and analyzes the two leads separately. Rank 1 uses only its three strict primary controls. Rank 2 uses only `length_sensitivity_exact_bin_pm7` and is supplemental/sensitivity-only. The two layers were not pooled, averaged together, or given equal evidentiary weight.

Every saved path for the ten involved entities was rediscovered from its AlphaFold request and required exactly five CIFs, five full-data files, and five confidence files. Request identity used candidate ID, the exact ordered three-chain sequences, and model seed. Jobs were deduplicated only when that request identity and all five ordered CIF SHA-256 hashes matched. The manifest retains 15 distinct jobs and records 4 identical saved-path copies.

For every retained model, the exact frozen P1-P9 core had one unambiguous placement. CIFs were parsed with the existing project parser. Exposed-position RMSD used P2/P3/P5/P7/P8 after fitting the HLA groove. Confidence was retained only as an annotation and did not select models or influence the 2.0-A complete-linkage pose clusters.

The hierarchical bootstrap used 10,000 iterations per lead and seed 20260815. It resampled unique jobs and then the five technical models within jobs, retaining all three controls at equal top-level weight.

## Results

### Rank 1: strict primary-control lead

- Pair: `EBV_TCELL_950::HUMAN_MYELIN_112214`
- Classification: `consistent_positive`
- Target median: 0.643 A
- Equal-weight background median: 7.964 A
- Background-minus-target delta: 7.321 A
- Leave-one-control-out delta range: 5.261 to 10.086 A
- Exploratory target rank: 1 of 4; empirical tail fraction 0.25
- Technical-stability delta interval: 1.620 to 12.759 A

### Rank 2: supplemental length sensitivity

- Pair: `EBV_TCELL_2268741::HUMAN_MYELIN_117032`
- Classification: `length_sensitivity_only__mixed_positive`
- Target median: 6.396 A
- Equal-weight background median: 13.975 A
- Background-minus-target delta: 7.578 A
- Leave-one-control-out delta range: 5.303 to 9.944 A
- Exploratory target rank: 1 of 4; empirical tail fraction 0.25
- Technical-stability delta interval: -0.517 to 16.880 A

The empirical tail fractions are exploratory ranks, not p-values. With three controls, 0.25 is the smallest possible fraction. The bootstrap intervals quantify technical stability across saved AlphaFold jobs/models only; they are not biological replication or p-values.

## Limitations and claim boundary

This is descriptive computational pMHC geometry. AlphaFold jobs and models are technical samples, not biological replicates. The small frozen control sets limit empirical resolution. Rank 2 has a deliberate peptide-length mismatch and remains sensitivity-only. These outputs do not establish peptide presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or an MS mechanism.

## Exact reproduce command

From the project root, with the output folder absent (the generator refuses overwrite):

```bash
PYTHONPATH=src python3 -m lead_focused_robustness --project-root "$PWD"
```

Project audited: `/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication`
