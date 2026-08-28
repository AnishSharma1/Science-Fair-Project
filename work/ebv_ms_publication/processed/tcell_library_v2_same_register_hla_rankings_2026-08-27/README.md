# T-cell Library V2 same-register rankings

This additive package ranks all V2 EBV-self pairs separately within each HLA using the control-selected P2/P3/P5/P7/P8 BLOSUM62 method. Existing structural packages remain unchanged.

- Full table: `all_hla_ranked_pairs.csv`
- Separate HLA tables: `rankings/`
- Compact view: `top_25_by_hla.csv`
- Exact top 10 for every HLA: `top_10_exact_epitopes_by_hla.csv`
- Exact selected candidate: `selected_epitope_pair.csv` and `selected_epitope_pair.json`
- Control comparison: `control_method_comparison.csv`
- Validation status: `ranking_gate.json`

The ranking is provisional because the method has three strict control systems rather than the six required for definitive validation. Score ties are explicit; `hla_rank` is a deterministic display order and `hla_score_rank` is the tied scientific rank. An asterisk marks a computationally prioritized exact peptide pair whose paired recognition has not been experimentally confirmed; it does not erase evidence attached to either individual peptide.

Reproduce with:

```bash
PYTHONPATH=src python3 src/build_same_register_hla_rankings_v2.py
```

Descriptive same-register pMHC sequence prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
