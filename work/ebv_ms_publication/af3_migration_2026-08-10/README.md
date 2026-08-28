# AlphaFold Server migration: audited 86 pMHC jobs

## Contents

- `pmhc_86_one_seed_alphafoldserver.json`: AlphaFold Server batch with 86 independent pMHC jobs.
- `pmhc_86_one_seed_job_index.csv`: maps each server job name to its exact chain order and peptide.
- `pmhc_priority_batch_01_one_seed.json` through `pmhc_priority_batch_05_one_seed.json`:
  five upload-ready priority batches (20, 20, 20, 20, and 6 jobs).
- `pmhc_86_priority_order.csv`: the order, batch assignment, and ranking rationale.

## Input policy

Each job contains exactly three protein chains, in this order:

1. mature HLA-DRA (178 aa)
2. mature HLA-DRB1*15:01 (189 aa)
3. the full, untrimmed IEDB peptide (10--40 aa)

`modelSeeds: []` requests AlphaFold Server's one automatically selected seed.
It is a one-seed screening pass, not an ensemble estimate.

## Submission

Upload the JSON through AlphaFold Server's **Upload JSON** control. It contains
86 jobs, below the server's 100-job import cap. Submit only within the service's
published quota for the account being used.

## Batch order

The batches are ranked using only pre-existing project evidence:

1. the BALF5--MBP calibration anchor;
2. candidates in the existing shortlist with nonzero **top-core** same-register
   alignment (not an optimized alternative-window result);
3. remaining existing shortlist candidates, ordered by their already recorded
   shortlist rank;
4. the remaining audited manifest candidates, with Tier 1 evidence ahead of
   lower-tier candidates.

This is an execution order, not a new biological score or a claim that a
candidate binds a TCR.

## Download and QC every completed job

Preserve the downloaded result zip, especially the original `*_job_request.json`,
the mmCIF coordinate file, full confidence JSON, and PAE file. Record whether:

- chain A/B/C are present and chain C exactly matches the input peptide;
- the peptide occupies the HLA-II groove rather than being detached or tangled;
- peptide-residue pLDDT and peptide-to-HLA PAE support a usable local geometry;
- any odd model is rerun with a new seed before pairwise interpretation.

pTM/ipTM are global confidence measures. They are QC fields here, not peptide-HLA
affinity, TCR binding, cross-reactivity, or disease-mechanism evidence.

## Before pairwise re-ranking

Run the same-register analysis on the AF3-derived structures only after the
BLAST-expanded candidate set has been filtered for HLA-DRB1*15:01 presentation
and each compared residue is mapped to a defensible P1--P9 register. Keep
pocket-facing and likely TCR-facing positions separate.
