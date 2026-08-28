# EBV–MS meeting candidate sheet

## One-sentence purpose

Use this sheet to distinguish a **known BALF5–MBP molecular-mimicry calibration system** from the project’s **computational candidates**. The score is a pMHC similarity screen, not evidence that a shared TCR binds a pair.

## What the columns mean

- **Local match**: number of residues included in the original local comparison.
- **Property similarity**: similarity of the aligned residues’ chemical properties in the computational screen.
- **Backbone RMSD**: geometric difference between peptide segments after fitting the HLA grooves; lower is more geometrically similar.
- **Register status**: whether equivalent P1–P9 positions in the HLA-II groove have been established. For the candidate pairs below, this remains unresolved.

## The positive-control anchor

| System | Why it matters | Structural observation | Meeting-safe interpretation |
|---|---|---|---|
| **EBV BALF5(627–641)** `TGGVYHFVKKHVHES` + **MBP(85–99)** `ENPVVHFFKNIVTPR` | Literature-established Hy.2E11 system across the DR15 haplotype | Experimental pMHC structures; known exposed segments BALF5 `YHFVKKH` and MBP `VHFFKNI`; backbone RMSD 0.838 Å after HLA-groove fitting | Calibration anchor, not a new discovery and not same-allele evidence. Ask which positions and exposed residues are genuinely comparable across DRB5*01:01 and DRB1*15:01. |

## Ranked screen records — interpreted correctly

| Original rank | Pair | Computational observation | How to discuss it |
|---|---|---|---|
| 1 | BALF5 + short MBP-85-region record `VVHFFKNIVTPRT` | 6-residue local match; property similarity 0.8431; backbone RMSD 0.332 Å | Same established BALF5–MBP control family. Not an independent discovery. Register equivalence unresolved. |
| 2 | BALF5 + MBP(85–99) `ENPVVHFFKNIVTPR` | 6-residue local match; property similarity 0.8431; backbone RMSD 0.370 Å | Overlapping version of the same control family. Collapse with rank 1 before any independent-pair count. |
| 3 | BALF5 + overlapping MBP-85-region record `ENPVVHFFKNIVTP` | 6-residue local match; property similarity 0.8431; backbone RMSD 0.381 Å | A third representation of the same control family, useful only as sensitivity/duplicate handling. |
| 4 | **EBV glycoprotein H** `EKQLFYYIGTMLPN` + **MBP** `QRPGFGYGGRASDYKSAHK` | 5-residue local match; property similarity 0.8428; backbone RMSD 1.260 Å | Highest-ranked non-MBP-85 duplicate candidate. Its P1–P9 registers and comparable exposed positions are not yet established. This is the main pair to ask Yicong about. |

## Useful contrast pairs

| Pair | Why it is useful |
|---|---|
| BALF5 + longer overlapping MBP record `ENPVVHFFKNIVTPRTP` (rank 5) | Same local chemistry but a much worse backbone RMSD (9.480 Å). Demonstrates why local sequence/chemistry alone is insufficient. |
| Glycoprotein H + different MBP record `LSRFSWGAEGQRPGFGYGG` (rank 31) | Same EBV glycoprotein-H query as rank 4 but poor geometric agreement (24.463 Å). Potential within-query negative only after HLA-binding matching. |

## The exact question for Yicong

> After collapsing the repeated BALF5–MBP records, glycoprotein H–MBP is our highest non-control candidate. Which P1–P9 register assignments are plausible for each peptide, and do any of the originally matched residues occupy equivalent, plausibly TCR-facing positions rather than different HLA anchor positions?

## Claim boundary

This sheet reports computational pMHC similarity candidates. It does not show shared-TCR binding, T-cell activation, cross-reactivity, or an MS mechanism.
