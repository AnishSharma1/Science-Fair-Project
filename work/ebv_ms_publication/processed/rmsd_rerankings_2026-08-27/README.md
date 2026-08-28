# RMSD sensitivity rankings

This additive package ranks every pair with a comparable structural measurement by median exposed-position P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove alignment. Lower RMSD ranks first. V2 alleles remain separate; the expanded DRB1*15:01 universe also receives its own structural rank.

Missing RMSDs are never imputed. Legacy whole-local-alignment RMSDs are not used. Eligible legacy RMSDs are recomputed from saved AF3 structures using the same register positions and atom-level endpoint.

The frozen RMSD endpoint failed the completed three-system control benchmark, capturing only 2 of 8 positive panels within the top three. Therefore these outputs are sensitivity analyses and do not replace the control-supported same-register BLOSUM62 rankings.

Exploratory structural sensitivity ranking; the frozen exposed-C-alpha RMSD endpoint did not pass the three-system control benchmark and is not the control-supported primary ranking.
Descriptive same-register pMHC sequence prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
