# V2 full-ensemble results summary

- Discovery jobs: **319/320** downloaded; **319** geometry-evaluable.
- Discovery pairs: **6360/6,400** geometry-complete.
- Calibration jobs: **22/24** downloaded; **22** geometry-evaluable.
- Calibration comparisons: **61/72** geometry-complete.
- Strict positive-recovery status: **not_evaluable_incomplete_calibration**.

## Gold-standard positive-control capture

The denominator was locked before ranking and contains one independent experimentally established system: Hy.2E11 recognition of DRB5*01:01-BALF5 and DRB1*15:01-MBP, with experimental pMHC structures 1H15 and 1BX2.

- Available-set capture@1: **2/2 seeds**.
- Fully evaluable seeds passing the predeclared rule: **1/1**.
- Strict two-seed status: **not_evaluable_incomplete_calibration**.
- Model or score changed to fit the positive: **no**.

This confirms capture of the known control in both available seed sets. It does not estimate broad sensitivity because the gold-standard denominator contains only one independent system, and one seed remains incomplete.

## Top cross-allele consensus pairs

Consensus minimizes the worst within-allele percentile, then the median, and requires all four alleles.

| Rank | EBV | Self | EBV protein | Self protein | Median percentile | Worst percentile |
|---:|---|---|---|---|---:|---:|
| 1 | ILCFVMAARQRLQDI | KTTICGKGLSATVT | EBNA3C | PLP1 | 0.0544 | 0.0776 |
| 2 | EKQLFYYIGTMLPN | KTTICGKGLSATVT | BXLF2_gH | PLP1 | 0.0694 | 0.0866 |
| 3 | IMNILRIYYSPSIM | KTTICGKGLSATVT | BFRF3 | PLP1 | 0.0410 | 0.0882 |
| 4 | ILCFVMAARQRLQDI | AYHYRKRGVHLAQGF | EBNA3C | ANO2 | 0.0447 | 0.1090 |
| 5 | QHYREVAAAKSSE | KTTICGKGLSATVT | BZLF1 | PLP1 | 0.0794 | 0.1194 |
| 6 | PPSIDPADLDESWD | KTTICGKGLSATVT | EBNA2 | PLP1 | 0.0538 | 0.1366 |
| 7 | QHYREVAAAKSSE | GKWLGHPDKFVG | BZLF1 | PLP1 | 0.0863 | 0.1388 |
| 8 | DNEIFLTKKMTEVCQ | KTTICGKGLSATVT | BALF4_gB | PLP1 | 0.0913 | 0.1572 |
| 9 | EKQLFYYIGTMLPN | GKWLGHPDKFVG | BXLF2_gH | PLP1 | 0.0935 | 0.1588 |
| 10 | TVFYNIPPMPL | GKWLGHPDKFVG | EBNA2 | PLP1 | 0.1361 | 0.1614 |

> Computational pMHC geometry only; not evidence of presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or MS disease mechanism.
