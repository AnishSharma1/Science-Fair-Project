# AF3 Hy.2E11 TCR-peptide contact consensus

Heavy-atom contacts use a 4.5 Å cutoff in the best-ranked model from each of five independent seed jobs. A stable residue/contact occurs in all five seed jobs. This is a hypothesis-generating consistency analysis; the template-excluded calibrators did not recover known ternary geometry, so these contacts cannot establish binding, specificity, or cross-reactivity.

## HY_MBP_DRB1

- Models: 5
- Mean TCR-peptide residue-pair contacts/model: 21.2
- Stable TCR residues: A93S, A94G, A95G, A96S, A97Y, B96W, B98S, B103Y
- Stable peptide residues: P4V, P7F, P9K, P10N, P12V
- Stable TCR–peptide pairs: A96S-P9K, A97Y-P7F, A97Y-P9K

## HY_BALF5_FULL15

- Models: 5
- Mean TCR-peptide residue-pair contacts/model: 26.2
- Stable TCR residues: A92D, A93S, A94G, A95G, A96S, A97Y, B28Q, B30T, B73L, B96W, B97P, B103Y, B105Y
- Stable peptide residues: P5Y, P6H, P7F, P8V, P9K, P10K, P11H, P12V
- Stable TCR–peptide pairs: A92D-P9K, A93S-P6H, A94G-P6H, A94G-P7F, A95G-P5Y, A95G-P6H, A95G-P7F, A96S-P7F, A96S-P9K, A97Y-P7F, A97Y-P8V, A97Y-P9K, A97Y-P10K, B28Q-P12V, B30T-P10K, B30T-P11H, B73L-P12V, B96W-P9K, B96W-P10K, B96W-P11H, B97P-P11H, B103Y-P9K, B103Y-P11H, B105Y-P9K

## HY_BALF5_CORE14

- Models: 5
- Mean TCR-peptide residue-pair contacts/model: 24.6
- Stable TCR residues: A92D, A93S, A94G, A95G, A96S, A97Y, B28Q, B30T, B73L, B96W, B97P, B103Y, B105Y
- Stable peptide residues: P4Y, P5H, P6F, P7V, P8K, P9K, P10H, P11V
- Stable TCR–peptide pairs: A92D-P8K, A93S-P5H, A94G-P5H, A94G-P6F, A95G-P4Y, A95G-P5H, A95G-P6F, A96S-P6F, A96S-P8K, A97Y-P6F, A97Y-P7V, A97Y-P8K, A97Y-P9K, B28Q-P11V, B30T-P9K, B30T-P10H, B73L-P11V, B96W-P8K, B96W-P9K, B97P-P10H, B103Y-P8K, B105Y-P8K

## DECOY_01_HY_MBP_DRB5

- Models: 5
- Mean TCR-peptide residue-pair contacts/model: 23.2
- Stable TCR residues: A29I, A30N, A93S, A94G, A95G, A97Y, B28Q, B30T, B96W
- Stable peptide residues: P4F, P5K, P6N, P7I, P8V
- Stable TCR–peptide pairs: A94G-P5K, A95G-P4F, A95G-P5K, A97Y-P6N, A97Y-P7I, B30T-P8V

## DECOY_02_HY_ENGA_DRB1

- Models: 5
- Mean TCR-peptide residue-pair contacts/model: 25.2
- Stable TCR residues: A93S, A94G, A95G, A97Y, B96W, B98S, B99G
- Stable peptide residues: P3F, P5R, P6V, P7H, P8F, P10S, P11A
- Stable TCR–peptide pairs: A97Y-P8F, B96W-P8F

## Hy.2E11 shared-contact comparison

- HY_MBP_DRB1__HY_BALF5_FULL15: stable-TCR-residue Jaccard = 0.5; shared = A93S, A94G, A95G, A96S, A97Y, B103Y, B96W
- HY_MBP_DRB1__HY_BALF5_CORE14: stable-TCR-residue Jaccard = 0.5; shared = A93S, A94G, A95G, A96S, A97Y, B103Y, B96W
- HY_BALF5_FULL15__HY_BALF5_CORE14: stable-TCR-residue Jaccard = 1.0; shared = A92D, A93S, A94G, A95G, A96S, A97Y, B103Y, B105Y, B28Q, B30T, B73L, B96W, B97P
