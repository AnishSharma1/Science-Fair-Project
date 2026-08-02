# Matched negative-control design

The current control design is deliberately explicit:

- EBV positives and negatives are deduplicated at the peptide level.
- Each positive is paired to a negative EBV peptide with the smallest absolute
  length difference; ties are resolved alphabetically.
- Reuse is allowed when the negative set is smaller than the positive set, and
  the control table records the pairing rather than hiding it.
- Human myelin peptides are compared separately with coordinate-validated human
  background peptides; they are not treated as biological replicates.
- Label permutations preserve the observed positive/negative group sizes.

This is adequate for a transparent exploratory checkpoint but not yet an
adequately powered inferential study. The next data expansion should seek
additional allele-matched negative EBV T-cell assays, especially at lengths
10–14, before making a publication-level claim.
