# GEO biological-context validation plan

## Purpose

GEO analysis is a **context** layer for the EBV--MBP molecular-mimicry
hypothesis. It cannot show pMHC presentation, TCR binding, or cross-reactivity.

## Selected datasets

| Dataset | Biological question | Use | Do not claim |
|---|---|---|---|
| GSE311678 | Which B-cell states are present in MS and how do anti-CD20 treatment and mucosal B-cell trafficking relate to them? | Assess HLA-II antigen-presentation program in identified B-cell states | EBV infection, BALF5 expression, or peptide presentation |
| GSE317492 | What are the transcriptional states of EBV-infected germinal-center B cells during primary infection? | Ask whether lytic EBV state markers co-occur with the host HLA-II antigen-presentation program | That BALF5 transcript proves the BALF5 peptide is presented by DRB1*15:01 |

## Pre-registered gene modules

### Host antigen-presentation module

`HLA-DRA`, `HLA-DRB1`, `CD74`, `CIITA`, `HLA-DMA`, `HLA-DMB`, `CTSS`.

### EBV lytic-context module

`BALF5` (target source gene), `BZLF1`, `BRLF1`, and any dataset-supplied lytic
state annotation.  BALF5 is lytic; lack of BALF5 in a primarily latent B-cell
dataset is expected and does not refute the molecular-mimicry hypothesis.

## Analysis sequence

1. Use the authors' cell-type annotations where supplied; otherwise identify
   B cells before testing either gene module.
2. Calculate per-cell module scores and pseudobulk per donor/state. Never treat
   individual cells as independent biological replicates.
3. In GSE311678, report the antigen-presentation module by B-cell state and
   treatment/disease group only if donor metadata support that comparison.
4. In GSE317492, compare the host antigen-presentation module between annotated
   EBV states; separately report whether viral reads include BALF5.
5. Label this result "biological context" in the manuscript. It is independent
   of, and does not rescue, the uncalibrated ternary docking result.

## Stop rules

- Do not merge the cohorts or calculate a cross-dataset disease effect.
- Do not test a large unregistered list of EBV genes and select the best one.
- Do not report BALF5 absence as evidence of absent antigen presentation.
