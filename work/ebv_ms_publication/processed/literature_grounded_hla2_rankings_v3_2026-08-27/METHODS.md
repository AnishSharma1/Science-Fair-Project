# Methods

Equivalent HLA-II grooves were aligned by Kabsch fitting of the first 85 C-alpha atoms from each alpha and beta chain. Each peptide P1-P9 surface fingerprint contains side-chain solvent accessibility, centroid position, centroid orientation, the exposed-residue distance matrix, scaled charge/hydropathy/donor/acceptor/aromaticity mismatch, and exposed-backbone geometry. Glycine uses its C-alpha atom as the side-chain fallback. Solvent accessibility uses a deterministic 32-point Shrake-Rupley approximation with a 1.4 A probe.

Every structural feature is summarized by its model-combination 75th percentile. Raw q25, median, q75, IQR, extrema, and per-feature within-HLA percentiles are retained. Register sensitivity enumerates all contained nine-residue window pairs for sequence and evaluates declared +/- 1 windows on the fixed modeled coordinates. These are sensitivity calculations, not alternative-register structure predictions.

Primary ranking is lexicographic: BLOSUM62 descending, physicochemical mismatch ascending, identity descending, local-surface percentile ascending, pair ID ascending. Binding predictions, full-core RMSD, and anchor RMSD are diagnostics only.

Descriptive HLA-specific pMHC sequence prioritization with modeled local-surface annotations only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
