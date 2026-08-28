# Targeted literature-control scan and inclusion boundary

## Result

No additional external system was added to the validation dataset before the meeting. The current evidence is insufficient to label any new record an independent EBV--myelin cross-reactivity positive without exact peptide, presenting allele, TCR, assay, and independence evidence.

This is deliberate: importing loosely related EBV, MBP, or HLA-DR15 records would inflate apparent validation while weakening the scientific claim.

## What was checked

| Source | What it contributes | What it cannot contribute here |
|---|---|---|
| [IEDB MHC-II API documentation](https://tools.iedb.org/main/tools-api/) | Supports the auditable `recommended_binding` workflow, including peptide core and rank output | Experimental presentation, a shared register, or cross-reactivity |
| [Scholz et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5566978/) | Supports the biological reason to keep DRB1*15:01 and DRB5*01:01 contexts explicit: their peptide repertoires are complementary | A direct BALF5--MBP TCR assay or a new external validation pair |
| Existing BALF5--MBP experimental-structure materials | One literature-established structural calibration system | Multiple independent validation examples |
| Existing Drosu and Wang panel annotations | Source/context overlays for prioritization bookkeeping | Direct EBV--myelin cross-reactive pair validation, unless exact functional evidence is added and audited |

## Inclusion rule for any future external positive control

Add a new system only when all of the following are documented before screening results are inspected:

1. Exact EBV and human peptide sequences.
2. Exact class-II allele for each peptide and the reason they may be compared.
3. Direct shared-TCR, binding, activation, or functional cross-reactivity evidence—not merely disease association or separate antigen reactivity.
4. Enough structural or experimentally anchored register information to make a P1--P9 comparison defensible.
5. Independence from the BALF5--MBP calibration family and from any duplicate peptide/protein record already counted.

## Meeting question

Ask Yicong and Olivia to name one or two systems that satisfy this rule, if they know of them. Do not pre-commit an external validation result until the system passes this audit.
