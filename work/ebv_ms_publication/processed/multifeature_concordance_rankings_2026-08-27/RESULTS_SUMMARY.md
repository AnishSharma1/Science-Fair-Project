# Multi-feature sequence-structure concordance rankings

Five displayed-pMHC proxies are shown and combined: TCR-facing physicochemical mismatch, TCR-facing BLOSUM62, TCR-facing same-register identity, exposed-register RMSD, and full P1-P9 core RMSD. Lower percentiles and lower concordance ranks are better.

The primary objective minimizes the worse of the sequence-family and structure-family percentiles. This rewards pairs that are jointly good in both families instead of allowing one excellent family to hide one poor family.

Retrospective control audit: **5/8 panels captured at rank 3**. This audit is exploratory and cannot freeze weights or unlock discovery.

## HLA-DRB1*15:01

Ranked: **1560**.

| Rank | Epitope pair | Sequence family | Structure family | Worst family | Exposed RMSD (A) | Full-core RMSD (A) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | BZLF1 198-210 / MAG 612-626* | 0.007 | 0.009 | 0.009 | 0.420 | 0.400 |
| 2 | EBNA3C 325-339 / MBP 84-97* | 0.014 | 0.024 | 0.024 | 0.539 | 0.500 |
| 3 | EBNA1 482-496 / ANO2 79-93* | 0.032 | 0.010 | 0.032 | 0.468 | 0.371 |
| 4 | BALF4/gB 103-116 / CRYAB 1-15* | 0.017 | 0.043 | 0.043 | 0.652 | 0.617 |
| 5 | EBNA3C 325-339 / MBP 84-98* | 0.014 | 0.043 | 0.043 | 0.686 | 0.594 |
| 6 | BNRF1 133-147 / ANO2 104-118* | 0.062 | 0.060 | 0.062 | 0.733 | 0.725 |
| 7 | EBNA3C 141-155 / MAG 612-626* | 0.042 | 0.073 | 0.073 | 0.810 | 0.798 |
| 8 | BHRF1 61-75 / TALDO1 216-230* | 0.062 | 0.089 | 0.089 | 0.865 | 0.927 |
| 9 | BALF5 627-641 / TALDO1 216-230* | 0.100 | 0.018 | 0.100 | 0.492 | 0.471 |
| 10 | EBNA1 482-496 / CLDN11 193-207* | 0.109 | 0.013 | 0.109 | 0.433 | 0.443 |

## HLA-DRB1*13:03

Ranked: **1600**.

| Rank | Epitope pair | Sequence family | Structure family | Worst family | Exposed RMSD (A) | Full-core RMSD (A) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | BARF1 1-15 / MAG 408-422* | 0.029 | 0.010 | 0.029 | 0.447 | 0.434 |
| 2 | BXLF2/gH 127-140 / CNP 136-150* | 0.056 | 0.045 | 0.056 | 0.518 | 0.719 |
| 3 | BNRF1 37-51 / PLP1 143-154* | 0.040 | 0.068 | 0.068 | 0.701 | 0.754 |
| 4 | BHRF1 122-133 / MAG 408-422* | 0.080 | 0.032 | 0.080 | 0.551 | 0.547 |
| 5 | BZLF1 198-210 / CNTN2 343-357* | 0.050 | 0.080 | 0.080 | 0.689 | 0.883 |
| 6 | BMRF1 126-140 / TALDO1 216-230* | 0.056 | 0.081 | 0.081 | 0.849 | 0.733 |
| 7 | BNRF1 37-51 / PLP1 105-118* | 0.019 | 0.093 | 0.093 | 0.920 | 0.776 |
| 8 | EBNA3C 141-155 / CNTN2 343-357* | 0.106 | 0.080 | 0.106 | 0.694 | 0.871 |
| 9 | BaRF1 185-199 / TALDO1 216-230* | 0.110 | 0.073 | 0.110 | 0.730 | 0.754 |
| 10 | EBNA2 280-290 / ANO2 129-143* | 0.121 | 0.095 | 0.121 | 0.850 | 0.846 |

## HLA-DRB1*03:01

Ranked: **1600**.

| Rank | Epitope pair | Sequence family | Structure family | Worst family | Exposed RMSD (A) | Full-core RMSD (A) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | BNRF1 37-51 / CLDN11 193-207* | 0.003 | 0.018 | 0.018 | 0.452 | 0.577 |
| 2 | BNRF1 37-51 / PLP1 105-118* | 0.015 | 0.021 | 0.021 | 0.518 | 0.565 |
| 3 | BALF2 854-866 / PLP1 179-192* | 0.014 | 0.043 | 0.043 | 0.660 | 0.806 |
| 4 | EBNA2 139-153 / MOG 97-109* | 0.046 | 0.060 | 0.060 | 0.692 | 1.232 |
| 5 | BXLF2/gH 127-140 / MBP 84-100* | 0.005 | 0.074 | 0.074 | 1.046 | 1.218 |
| 6 | BXLF2/gH 127-140 / MBP 86-100* | 0.005 | 0.074 | 0.074 | 1.045 | 1.227 |
| 7 | EBNA1 475-489 / MOG 97-109* | 0.070 | 0.074 | 0.074 | 1.000 | 1.466 |
| 8 | BALF4/gB 103-116 / CRYAB 1-15* | 0.015 | 0.075 | 0.075 | 1.064 | 1.187 |
| 9 | EBNA2 280-290 / CLDN11 193-207* | 0.081 | 0.055 | 0.081 | 0.814 | 0.906 |
| 10 | BALF4/gB 576-590 / ANO2 79-93* | 0.098 | 0.053 | 0.098 | 0.811 | 0.854 |

## HLA-DRB1*08:01

Ranked: **1600**.

| Rank | Epitope pair | Sequence family | Structure family | Worst family | Exposed RMSD (A) | Full-core RMSD (A) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | BARF1 1-15 / MAG 408-422* | 0.029 | 0.003 | 0.029 | 0.316 | 0.413 |
| 2 | BARF1 1-15 / ANO2 79-93* | 0.040 | 0.019 | 0.040 | 0.444 | 0.621 |
| 3 | BZLF1 198-210 / CNTN2 343-357* | 0.047 | 0.047 | 0.047 | 0.645 | 0.777 |
| 4 | BXLF2/gH 127-140 / CNP 136-150* | 0.058 | 0.027 | 0.058 | 0.490 | 0.696 |
| 5 | BARF1 1-15 / CLDN11 65-79* | 0.014 | 0.074 | 0.074 | 0.866 | 0.797 |
| 6 | BNRF1 81-95 / MAG 205-219* | 0.033 | 0.092 | 0.092 | 0.847 | 0.954 |
| 7 | EBNA3C 141-155 / CNTN2 343-357* | 0.106 | 0.086 | 0.106 | 0.830 | 0.913 |
| 8 | BFRF1 125-139 / ANO2 104-118* | 0.111 | 0.035 | 0.111 | 0.569 | 0.723 |
| 9 | EBNA1 482-496 / CLDN11 193-207* | 0.117 | 0.094 | 0.117 | 0.917 | 0.904 |
| 10 | EBNA2 446-459 / MAG 408-422* | 0.121 | 0.084 | 0.121 | 0.924 | 0.819 |

## Combined DRB1*15:01 universe

Ranked: **1913**; missing comparable structure: **130**.

| Rank | Epitope pair | Sequence family | Structure family | Worst family | Sequence rank | RMSD rank |
|---:|---|---:|---:|---:|---:|---:|
| 1 | BZLF1 198-210 / MAG 612-626* | 0.008 | 0.008 | 0.008 | 11 | 22 |
| 2 | EBNA3C 325-339 / MBP 84-97* | 0.015 | 0.023 | 0.023 | 13 | 53 |
| 3 | EBNA1 482-496 / ANO2 79-93* | 0.032 | 0.009 | 0.032 | 1 | 32 |
| 4 | EBNA3C 325-339 / MBP 84-98* | 0.015 | 0.042 | 0.042 | 15 | 95 |
| 5 | BALF4/gB 103-116 / CRYAB 1-15* | 0.017 | 0.042 | 0.042 | 9 | 87 |
| 6 | BNRF1 133-147 / ANO2 104-118* | 0.061 | 0.058 | 0.061 | 42 | 114 |
| 7 | EBNA3C 141-155 / MAG 612-626* | 0.042 | 0.071 | 0.071 | 49 | 140 |
| 8 | EBNA1 482-496 / MBP 245-263* | 0.049 | 0.072 | 0.072 | 162 | 123 |
| 9 | EBNA1 482-496 / PLP1 181-200* | 0.084 | 0.011 | 0.084 | 206 | 27 |
| 10 | BHRF1 61-75 / TALDO1 216-230* | 0.061 | 0.085 | 0.085 | 34 | 156 |
| 11 | BALF5 627-641 / TALDO1 216-230* | 0.097 | 0.017 | 0.097 | 19 | 40 |
| 12 | EBNA1 482-496 / CLDN11 193-207* | 0.110 | 0.012 | 0.110 | 222 | 23 |
| 13 | LMP2 72-86 / ANO2 154-168* | 0.071 | 0.112 | 0.112 | 67 | 213 |
| 14 | BZLF1 61-75 / MBP 84-97* | 0.117 | 0.088 | 0.117 | 122 | 183 |
| 15 | BZLF1 61-75 / MBP 84-98* | 0.117 | 0.098 | 0.117 | 124 | 198 |
| 16 | BFRF3 289-302 / CLDN11 193-207* | 0.117 | 0.046 | 0.117 | 228 | 92 |
| 17 | BARF1 25-39 / CRYAB 1-15* | 0.078 | 0.120 | 0.120 | 25 | 233 |
| 18 | EBNA2 280-290 / TALDO1 216-230* | 0.081 | 0.123 | 0.123 | 59 | 243 |
| 19 | LMP2 73-87 / PLP1 105-118* | 0.128 | 0.084 | 0.128 | 371 | 169 |
| 20 | BALF2 854-868 / ANO2 79-93* | 0.131 | 0.046 | 0.131 | 147 | 46 |
| 21 | BALF4/gB 103-116 / MOG 221-235* | 0.068 | 0.134 | 0.134 | 45 | 257 |
| 22 | LMP1 43-53 / CLDN11 193-207* | 0.140 | 0.041 | 0.140 | 300 | 61 |
| 23 | LMP2 73-87 / PLP1 140-152* | 0.140 | 0.109 | 0.140 | 450 | 216 |
| 24 | BZLF1 198-210 / MAG 205-219* | 0.038 | 0.141 | 0.141 | 24 | 267 |
| 25 | BHRF1 61-75 / PLP1 143-154* | 0.145 | 0.030 | 0.145 | 243 | 60 |

V2 missing comparable structure across all HLAs: **40**.

> Full-core RMSD means P1-P9 peptide-core C-alpha RMSD after HLA-groove alignment. It is not whole-source-protein RMSD and not whole-pMHC RMSD.

> Descriptive same-register pMHC sequence-structure concordance prioritization only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate.
