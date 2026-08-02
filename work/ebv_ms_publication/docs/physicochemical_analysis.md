# Physicochemical similarity checkpoint

The rebuilt screen adds a transparent descriptor model to sequence identity.
For each local alignment of at least five residues, it averages four
properties: normalized Kyte–Doolittle hydropathy, charge compatibility,
aromatic compatibility, and coarse size compatibility.

The positive EBV group had 12 unique peptides; the negative EBV group had 7.
Against 53 coordinate-validated myelin peptides, the positive-minus-negative
difference in peptide-level maximum similarity was **−0.0399**, with a
10,000-permutation label-test p-value of **0.697**. This does not support a
positive similarity signal in the current small dataset.

The result is exploratory. The descriptor model is not trained, and it is not
a substitute for an MHC-binding predictor, structural model, or TCR assay.
The exact-control limitation is important: no negative EBV peptides of length
10 were available, so the two length-10 positives were matched to the nearest
available negative length (14 residues).
