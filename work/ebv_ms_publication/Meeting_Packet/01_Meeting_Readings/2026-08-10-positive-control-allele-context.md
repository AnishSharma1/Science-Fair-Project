# BALF5--MBP positive-control allele context

## Why this matters

The established BALF5--MBP mimicry anchor should not be silently treated as a same-allele HLA-DRB1*15:01 benchmark. In the experimental calibration context used by this project, BALF5 is associated with HLA-DRB5*01:01 (DR2a) and MBP with HLA-DRB1*15:01 (DR2b). These are linked on the DR15 haplotype but have distinct peptide repertoires. A primary study of these molecules reports complementary HLA-DR15 peptide presentation rather than interchangeable binding behavior ([Scholz et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5566978/)).

## Reproducible IEDB computational context check

All three rows were queried with IEDB MHC-II `recommended_binding`, retaining the returned top core, percentile rank, and IC50 hypothesis. IEDB documents this endpoint and its supported peptide lengths; it is a predictor, not experimental presentation evidence ([IEDB MHC-II API documentation](https://tools.iedb.org/main/tools-api/)).

| Peptide / allele context | IEDB top core | Percentile rank | IC50 (nM) | Correct use |
|---|---:|---:|---:|---|
| BALF5 `TGGVYHFVKKHVHES` / DRB5*01:01 | `YHFVKKHVH` | 11.0 | 126.27 | Computational DR2a context for the structural calibration system |
| BALF5 `TGGVYHFVKKHVHES` / DRB1*15:01 | `VYHFVKKHV` | 20.0 | 314.19 | Contrast only; not evidence that the experimental BALF5 complex is DRB1-presented |
| MBP `ENPVVHFFKNIVTPR` / DRB1*15:01 | `VHFFKNIVT` | 0.08 | 11.91 | Computational DR2b context for the MBP side |

The output is saved in `processed/register_sensitivity/positive_control_allele_context.csv`; the per-query raw API records are in `processed/register_sensitivity/iedb_mhcii_positive_control_raw.txt`.

## What this does and does not justify

It justifies a precise meeting question: **what kind of cross-allele, DR15-haplotype structural comparison is legitimate for calibration?** It does not justify calling the two peptides register-equivalent, proving a shared TCR footprint, or extending the result to the DRB1-only screen.

Suggested manuscript-safe sentence:

> Experimental BALF5--MBP structures were retained as a DR15-haplotype structural calibration anchor. Because the established complexes use DRB5*01:01 and DRB1*15:01 contexts, respectively, they were not treated as same-allele validation of the DRB1*15:01 computational screen.
