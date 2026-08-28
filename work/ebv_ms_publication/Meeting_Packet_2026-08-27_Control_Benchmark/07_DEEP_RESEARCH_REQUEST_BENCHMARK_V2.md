# Deep research request: how to validate and advance the EBV-MS HLA-II pMHC prioritization model

## Role

Act as a rigorous computational-immunology research analyst. Use current primary literature, experimentally curated structural databases, and immunological assay databases to determine how this project should proceed. The goal is not to maximize favorable results. The goal is to design a benchmark that can credibly distinguish useful structural signal from simple peptide-sequence similarity and prevent data leakage or circular validation.

Search through the date this request is run. Cite direct links to primary papers and database records. Clearly distinguish verified facts, reasonable inferences, and unresolved items.

## Project objective

This project prioritizes EBV and human-self peptide pairs whose HLA class-II-bound conformations may be similar at positions potentially exposed to a T-cell receptor. It is a computational pMHC prioritization project.

It does **not** currently establish:

- antigen processing or presentation;
- recognition or binding by the same TCR;
- T-cell activation;
- cross-reactivity or molecular mimicry;
- causation or mechanism in multiple sclerosis;
- a probability, false-discovery rate, or clinical risk score.

The eventual discovery output must be ranked separately within each exact HLA alpha/beta allotype. Do not recommend a cross-allele consensus ranking.

## Fixed current benchmark

The current benchmark contains three independent human alpha-beta TCR systems and four required positive comparisons:

1. **Hy.2E11: EBV BALF5-MBP**
   - BALF5 structure: PDB 1H15, HLA-DRA*01:01/DRB5*01:01
   - MBP structure: PDB 1BX2, HLA-DRA*01:01/DRB1*15:01
   - Primary source: DOI 10.1038/ni835
2. **Ob.1A12: bacterial EngA-MBP**
   - EngA structure: PDB 2WBJ
   - MBP structure: PDB 1YMM
   - HLA-DRA*01:01/DRB1*15:01
   - Primary source: DOI 10.1016/j.immuni.2009.01.009
3. **Hy.1B11: HSV UL15-MBP and Pseudomonas PMM-MBP**
   - UL15 structure: PDB 4MAY
   - PMM structure: PDB 4GRL
   - MBP structure: PDB 3PL6
   - HLA-DQA1*01:02/DQB1*05:02
   - Primary source: DOI 10.1038/ncomms3623

Hy.1B11 contributes two required positive comparisons but only one independent biological-system vote.

EBNA1-ANO2 is prospective only. Do not admit it as a strict control unless the same-clone identity, exact peptides, exact nine-residue registers, exact HLA alpha/beta allotypes, functional evidence for both arms, and structures are resolved.

## Fixed computational results

Treat the following as supplied project data. Do not reinterpret missing panels as completed results.

- 48 of 48 newly prepared AlphaFold jobs were recovered as exact five-model bundles.
- All 240 new model samples passed exact-chain and clash checks.
- Including reused frozen models, 350 valid model samples were analyzed.
- Eight required panels were complete and four required evaluations were incomplete.
- The nested leave-one-biological-system-out candidate composite ranked the positive within the top 3 in 8/8 completed panels.
- Composite held-out ranks: 1, 1, 2, 2, 3, 3, 1, and 1; worst rank 3.
- Frozen exposed-position C-alpha RMSD captured 3/8 completed panels; worst rank 9.
- TCR-facing sequence identity captured 8/8 completed panels; worst rank 2.
- Full nine-residue core identity captured 7/8 completed panels; worst rank 4.
- For held-out Hy.2E11 and Ob.1A12, weight selection chose 100% TCR-facing physicochemical mismatch.
- For held-out Hy.1B11, it chose 50% exposed C-alpha RMSD and 50% physicochemical mismatch.
- Side-chain-vector and anchor C-alpha features received zero weight in every outer fold.
- The formal trust gate is `not_evaluable`, not `pass` and not `fail`.
- No completed composite panel failed the top-3 rule.
- Two AlphaFold seed panels remain incomplete because frozen jobs are missing.
- Both Hy.1B11 PDB panels are not evaluable because each has only two eligible exact-HLA decoys rather than the required five.
- Candidate weights remain unfrozen, discovery rankings remain unchanged, and no cross-allele consensus was created.

The central concern is that control recovery may be driven mainly by sequence or physicochemical resemblance rather than independent three-dimensional information.

## Strict definitions

Use these definitions consistently:

- **Independent biological system:** one paired human alpha-beta TCR or explicitly identified human T-cell clone.
- **Required positive comparison:** two exact pMHC ligands functionally recognized by that same TCR or clone.
- **N1 negative:** same TCR explicitly tested and shown not to recognize the ligand under a relevant assay.
- **N2 negative:** negative from the same assay context and exact HLA, even if not tested with precisely the same TCR configuration.
- **N3 comparator:** exact-HLA, score-blind matched peptide with unknown TCR recognition. N3 supports ranking comparisons only and cannot support specificity claims.
- **Strict structural control:** both exact peptide sequences, exact HLA alpha/beta allotypes, experimentally resolved nine-residue registers, paired-TCR/clone identity, functional evidence for both arms, and structural coordinates are verified.

Do not convert shorthand such as "DR15," "DR2," or "DQ1" into exact allotypes without direct evidence. Do not treat a predicted binding register as experimentally resolved. Do not count multiple ligands for one TCR as multiple independent systems.

## Research workstreams

### 1. Find additional independent positive-control systems

Search TCR3d, RCSB PDB, IEDB, VDJdb, McPAS-TCR, ATLAS, and primary literature for human alpha-beta TCRs or explicitly identified human clones that recognize at least two distinct HLA class-II ligands.

For every candidate system, verify:

- exact paired TCR identity or exact clone identity;
- exact peptide sequences and source proteins/organisms;
- exact HLA alpha and beta allotypes for each ligand;
- exact experimentally resolved P1-P9 register for each ligand;
- functional or binding evidence for both arms;
- PDB identifiers and chain assignments;
- whether ligands represent distinct biological sources rather than engineered variants;
- independence from the current three control systems.

Classify each system as `strictly_eligible`, `promising_but_incomplete`, or `ineligible`, with explicit reasons. Do not hide excluded systems.

### 2. Find specificity negatives

Search the papers, supplements, assay tables, and related studies for N1 and N2 negatives associated with Hy.2E11, Ob.1A12, Hy.1B11, or newly identified systems.

Prioritize:

- peptides tested against the same TCR with reported nonrecognition;
- altered peptide ligands with interpretable negative results;
- exact-HLA peptides tested in the same assay;
- negative tetramer, binding, cytokine, proliferation, cytotoxicity, or activation results with clear experimental context.

For every negative, report the assay, peptide, HLA, concentration or relevant condition, outcome, source location, and whether it qualifies as N1, N2, or neither. Explain assay comparability limitations.

### 3. Expand or assess the HLA-DQ PDB oracle layer

Search for structures containing the exact HLA-DQA1*01:02/DQB1*05:02 heterodimer with experimentally resolved peptides and usable groove coordinates.

For each structure, verify exact allotypes from sequence or authoritative annotation, peptide sequence, experimental register, resolution, receptor-bound status, technical duplication, and PDB chains. Identify a score-blind deduplicated structural pool suitable for noncognate pair decoys.

Determine whether at least five defensible decoys can be produced for each Hy.1B11 positive without reusing technical duplicates or uncertain-allotype structures. If not, recommend one of these options with justification:

- retain the PDB layer as `not_evaluable`;
- predeclare a lower decoy threshold for a future benchmark;
- use a larger exact-HLA AlphaFold-only validation layer;
- recruit a different DQ control system with a better structural pool.

Do not lower the current benchmark threshold retroactively.

### 4. Evaluate whether the model adds value beyond sequence

Review methodological literature on TCR-pMHC similarity, molecular-mimicry prediction, peptide structural comparison, HLA-II register-aware scoring, and leakage-controlled TCR benchmarking.

Assess which baselines should be mandatory:

- TCR-facing exact sequence identity;
- full-core identity;
- BLOSUM or substitution-matrix similarity;
- physicochemical distance without structures;
- binding-percentile similarity;
- peptide length and register agreement;
- random ranking;
- simple experimental-PDB geometry when available.

Recommend how to test incremental value from AlphaFold-derived 3D features through ablation, nested validation, permutation testing, uncertainty intervals, and paired comparisons across systems. Address the small number of independent TCR systems.

Answer explicitly:

- Must a structural score strictly outperform the best sequence-only baseline before it can be called structurally validated?
- Would matching the sequence baseline be acceptable if the structural score improves interpretability or external generalization?
- Should a frozen composite require a nonzero structural-feature weight?
- How should model selection avoid exploiting control labels or decoy composition?
- How many independent systems are realistically needed before weights can be frozen?

### 5. Design a prospective benchmark version 2

Propose a fully prospective, preregistered-in-spirit benchmark that is not altered after reading results.

The design should address:

- a fresh complete AlphaFold manifest rather than silently patching missing frozen jobs;
- at least two fixed AlphaFold seeds or an evidence-based alternative;
- five score-blind comparator ligands on each arm and 25 pair decoys per panel;
- exact-HLA matching and register verification;
- system-level leave-one-out or another leakage-resistant split;
- multi-ligand systems receiving one biological vote;
- separation of N1/N2 specificity testing from N3 ranking;
- feature ablations and sequence-only baselines;
- deterministic tie-breaking and checksums;
- explicit handling of missing results;
- predeclared go/no-go criteria for discovery reranking;
- separate final rankings for each HLA only.

Compare at least three possible validation rules, such as:

1. capture-at-3 in every required panel;
2. capture-at-3 plus strict improvement over the strongest sequence-only baseline;
3. capture-at-3 plus noninferiority to sequence, a required nonzero structural component, and successful N1/N2 specificity discrimination.

Recommend one rule and explain its statistical and biological tradeoffs.

### 6. Determine the scientifically defensible interpretation of the current results

Answer these questions directly:

1. Can the current result be described as "control-validated," and if so, exactly what is validated?
2. Is "encouraging held-out control recovery" more accurate than "validated structural model"?
3. Does the dominance of sequence and physicochemical features undermine the structural claim?
4. What claims are appropriate for a mentor call, abstract, methods section, poster, or manuscript?
5. What evidence would be required before applying the score to EBV-self discovery candidates?

Provide recommended language and prohibited language.

## Evidence standards

Prioritize sources in this order:

1. primary functional studies and their supplements;
2. RCSB PDB records and coordinate files;
3. TCR3d and experimentally curated immune-receptor databases;
4. primary benchmarking or methods papers;
5. reviews only for discovery and context.

For every important factual claim, provide a direct citation. Include DOI, PMID, PDB ID, database accession, and exact table/figure/supplement location when possible.

Do not infer missing peptide residues, HLA allotypes, registers, TCR identities, assay outcomes, or negative results. Label unavailable information as `not_reported` or `unresolved`. Distinguish an absence of reported recognition from experimentally demonstrated nonrecognition.

## Required deliverables

Return the final research report in this order:

### A. One-page executive decision

- current scientific verdict;
- whether discovery reranking should remain locked;
- three most important reasons;
- recommended benchmark-v2 route;
- top five next actions.

### B. Current-result interpretation table

Columns:

`claim | supported_now | exact_evidence | limitation | recommended_wording`

### C. Additional positive-control registry

Columns:

`system_id | TCR_or_clone | species | ligand_1 | ligand_2 | exact_HLA_1 | exact_HLA_2 | peptide_1 | peptide_2 | core_1 | core_2 | PDB_1 | PDB_2 | functional_evidence | register_evidence | eligibility | exclusion_reason | primary_citation`

### D. Negative-control registry

Columns:

`TCR_system | peptide | exact_HLA | assay | tested_condition | outcome | negative_tier | comparability | source_location | citation`

### E. Exact-HLA DQ structural-decoy registry

Columns:

`PDB_ID | alpha_allotype | beta_allotype | peptide | core | register_source | resolution | receptor_bound | duplicate_group | usable_decoy | exclusion_reason | evidence_link`

### F. Baseline and ablation recommendation

Specify mandatory baselines, evaluation metrics, statistical comparisons, and the criterion for claiming incremental structural value.

### G. Benchmark-v2 protocol

Provide a numbered, implementation-ready protocol covering curation, freezing, model generation, quality control, scoring, nested validation, missingness, trust-gate logic, and final HLA-specific rankings.

### H. Decision matrix

Compare at least three benchmark-v2 options by rigor, feasibility, cost, time, leakage risk, interpretability, and ability to support specificity versus ranking claims.

### I. Prioritized action plan

Separate actions into:

- before the mentor call;
- next 48 hours;
- next two weeks;
- longer-term experimental validation.

For each action, state the expected output and the decision it enables.

### J. Unresolved questions and search log

List unresolved identities, failed searches, ambiguous allotypes/registers, unavailable supplements, and evidence gaps. Include the exact databases and search strategies used so the research can be reproduced.

## Final quality-control checklist

Before completing the report, verify that:

- every admitted positive system uses the same paired TCR or explicitly identified clone;
- every strict peptide and HLA allotype is exact rather than inferred from shorthand;
- every strict register is experimentally resolved;
- one TCR system receives only one independent vote;
- N3 comparators are never presented as specificity negatives;
- held-out systems are excluded from weight selection and decoy optimization;
- sequence-only baselines are compared directly with the composite;
- missing panels cannot become passes;
- no cross-allele discovery consensus is recommended;
- all disease-mechanism and TCR-cross-reactivity claims remain outside the computational evidence.

Conclude with a clear recommendation: **proceed**, **proceed only after specified benchmark changes**, or **do not proceed with discovery reranking**, and state the evidence threshold that would change that decision.
