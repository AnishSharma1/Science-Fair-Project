# Register-Aware Computational Scoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a reproducible, claim-bounded computational result package that recomputes sequence and physicochemical features only at matched HLA-DRB1*15:01 register positions, evaluates them against pre-existing strict decoys, and reports the present structural-geometry limitation honestly.

**Architecture:** A small pure-Python scoring module will turn an existing local-alignment coordinate map and two resolved P1-P9 cores into per-position, anchor-class, exposed-class, and total sequence/chemistry features. A builder will join those scores to the frozen benchmark universe and retain the existing aggregate geometry value only as labelled context. A separate decoy evaluator will calculate within-target ranks and paired contrasts without reselecting decoys or pooling correlated evidence records.

**Tech Stack:** Python 3 standard library (`csv`, `pathlib`, `statistics`, `unittest`), existing processed CSVs, and existing `src/biochemical_similarity.py` property definitions. No new packages, structure predictions, web calls, or wet-lab data.

## Global Constraints

- Treat `processed/register_aware_benchmark/benchmark_pair_universe.csv`, `target_feasibility.csv`, and `matched_decoy_sets.csv` as frozen inputs to scoring; do not modify decoy selection while calculating scores.
- Primary scoring admits only `register_assessment == "assessable_register_hypothesis"`. Calibration-only and sensitivity-only rows remain exported but never receive a primary rank or decoy comparison.
- Preserve the experimental-register hierarchy already implemented: same-allele 1BX2 MBP evidence can override prediction; DRB5*01:01 1H15 BALF5 is calibration only; unresolved gH and candidate-MBP calls are sensitivity only.
- The available geometry matrix supplies a whole-original-local-alignment RMSD, not per-residue coordinates. Never call that number a same-register geometry score and never include it in a primary composite score until the original model coordinates are recovered and recomputed.
- Use `P1/P4/P6/P9` only as HLA-anchor-position labels and `P2/P3/P5/P7/P8` only as non-anchor/candidate-exposed-position labels. Neither label supports a TCR-contact assertion.
- A row with one or two matched-register residues is reportable as limited-coverage evidence but must not enter the robust primary ranking set. The robust set is predeclared as at least three matched-register aligned residues; current frozen coverage is 12 robust rows and 21 limited-coverage rows.
- Do not claim TCR binding, activation, cross-reactivity, established molecular mimicry, or an MS mechanism. The allowed conclusion remains computational pMHC candidate prioritization relative to matched decoys.
- The copied bundle is not a Git repository. Leave generated artifacts and test output; do not initialize, commit, or publish a repository.

---

### Task 1: Implement deterministic register-filtered feature primitives

**Files:**
- Create: `src/register_aware_scoring.py`
- Create: `tests/test_register_aware_scoring.py`

**Interfaces:**
- Consumes: `original_local_alignment_coordinates`, peptide sequences, and one-based P1-P9 core starts from a benchmark row.
- Produces: `score_same_register_alignment(row) -> dict[str, object]`.
- Required score fields: `same_register_alignment_count`, `same_register_positions`, `same_register_identity_count`, `same_register_identity_fraction`, mean hydropathy/charge/aromatic/size similarities, `same_register_property_similarity`, and the corresponding `all`, `anchor`, and `candidate_exposed` class counts/scores.

- [ ] **Step 1: Write failing unit tests for coordinate/register filtering**

```python
row = {
    "ebv_peptide": "ABCDEFGHIJKLMNO",
    "human_peptide": "ABQDEFGHIJKLMNO",
    "ebv_top_core_start_1_based": "4",
    "human_top_core_start_1_based": "4",
    "original_local_alignment_coordinates": "3:3;4:4;5:5;6:6;7:7;8:8;9:9;10:10;11:11",
}
score = score_same_register_alignment(row)
assert score["all_same_register_positions"] == "P1;P2;P3;P4;P5;P6;P7;P8"
assert score["anchor_same_register_positions"] == "P1;P4;P6"
assert score["candidate_exposed_same_register_positions"] == "P2;P3;P5;P7;P8"
```

Add negative tests for malformed coordinates, positions outside either peptide, and a coordinate whose EBV/Human P1-P9 indices differ. These must fail clearly rather than silently creating a score.

- [ ] **Step 2: Run the new tests and confirm they fail**

Run: `PYTHONPATH=src python3 -m unittest tests.test_register_aware_scoring -v`

Expected: FAIL because the scorer does not yet exist.

- [ ] **Step 3: Implement pure scoring helpers**

Implement and test these helpers in `src/register_aware_scoring.py`:

```python
parse_local_alignment_positions(text) -> list[tuple[int, int]]
register_position(peptide_index_1_based, core_start_1_based) -> int | None
property_similarity(ebv_residue, human_residue) -> dict[str, float]
score_same_register_alignment(row) -> dict[str, object]
```

Reuse the existing Kyte-Doolittle, charge, aromaticity, and size conventions from `src/biochemical_similarity.py`; document the exact similarity formula next to the new helper. Retain every qualifying coordinate in ascending EBV position order. For each class, return an empty score rather than a zero when no positions qualify, so “no measurement” cannot be confused with poor similarity.

Set `score_coverage_status` exactly as follows:

| Condition | Status |
| --- | --- |
| not assessable register hypothesis | `excluded_nonprimary_register_status` |
| zero qualifying coordinates | `no_same_register_alignment` |
| 1–2 qualifying coordinates | `limited_coverage_report_only` |
| at least 3 qualifying coordinates | `robust_primary_ranking_eligible` |

- [ ] **Step 4: Add deterministic property and class tests**

Test identical residues yield identity fraction and each available property similarity of `1.0`; test a known charged mismatch; test P1/P4/P6/P9 versus P2/P3/P5/P7/P8 partitioning; and test that the all-position result equals the weighted mean of nonempty classes. Verify that a calibration-only row cannot be labelled robust merely because its coordinates match.

- [ ] **Step 5: Run the module tests**

Run: `PYTHONPATH=src python3 -m unittest tests.test_register_aware_scoring -v`

Expected: PASS.

### Task 2: Build the auditable primary score table and geometry limitation table

**Files:**
- Create: `src/build_register_aware_score_table.py`
- Modify: `tests/test_register_aware_scoring.py`
- Create: `processed/register_aware_scoring/register_aware_pair_scores.csv`
- Create: `processed/register_aware_scoring/score_coverage_summary.csv`
- Create: `processed/register_aware_scoring/register_aware_geometry_limitations.csv`
- Create: `processed/register_aware_scoring/README.md`

**Interfaces:**
- Consumes: frozen pair universe plus `processed/colabfold_tier1_ebv_myelin_geometry_matrix.csv`.
- Produces: one score row for every universe pair, keyed by `pair_id`, plus an explicit geometry interpretation record.

- [ ] **Step 1: Write failing builder tests**

Create tiny temporary CSV fixtures. Assert that every universe row produces exactly one score row; all source/status columns are copied; and only assessable rows call the feature scorer. Verify that an unmatched `pair_id` raises an error rather than dropping geometry data.

- [ ] **Step 2: Implement an explicit geometry join**

Join the existing geometry column `local_peptide_ca_rmsd_after_hla_fit` by `pair_id` without transforming it. Add these fields to every score row:

```text
whole_local_alignment_geometry_rmsd_context
geometry_context_status
geometry_primary_score_eligible
geometry_context_note
```

For present data, set:

```text
geometry_context_status = whole_original_local_alignment_only
geometry_primary_score_eligible = False
geometry_context_note = Precomputed RMSD spans the original local alignment; no per-residue model coordinates are available to recompute a same-register-only RMSD.
```

Write `register_aware_geometry_limitations.csv` as a one-row machine-readable record of that limitation, the expected coordinate inputs, and the exact recovery condition: recover both original pMHC PDBs for every scored pair, HLA-fit them with the prior procedure, and recompute RMSD using only retained matched-register residue pairs.

- [ ] **Step 3: Implement coverage summaries**

Group score rows by `score_coverage_status`, `register_assessment`, and `same_register_alignment_count`. Include the frozen expected coverage check: 33 assessable rows total, with 11 rows at one matched residue, 10 at two, 9 at three, and 3 at four. Fail the build if these counts drift unexpectedly, because that signals an upstream eligibility or coordinate-map change.

- [ ] **Step 4: Regenerate artifacts**

Run: `PYTHONPATH=src python3 src/build_register_aware_score_table.py`

Expected: `register_aware_pair_scores.csv` contains all 636 PASS-universe rows; 12 are `robust_primary_ranking_eligible`, 21 are `limited_coverage_report_only`, and calibration/sensitivity rows remain present but excluded.

- [ ] **Step 5: Document the table contract**

In the generated README, define all score fields, state that raw scores are not a biological mimicry probability, identify position labels as descriptive, and explain why the geometry RMSD is context only. Include the exact regeneration and test commands.

### Task 3: Evaluate frozen strict decoys without invalid pooled statistics

**Files:**
- Create: `src/run_register_aware_decoy_evaluation.py`
- Create: `tests/test_register_aware_decoy_evaluation.py`
- Create: `processed/register_aware_scoring/decoy_score_comparison.csv`
- Create: `processed/register_aware_scoring/decoy_evaluation_summary.csv`
- Modify: `processed/register_aware_scoring/README.md`

**Interfaces:**
- Consumes: the score table, frozen `matched_decoy_sets.csv`, and `target_feasibility.csv`.
- Produces: target-versus-decoy comparisons only when target and all five preselected decoys have robust primary-ranking scores.

- [ ] **Step 1: Predeclare the score used for decoy comparison**

Use `same_register_property_similarity` as the primary sequence/chemistry endpoint. Do **not** create a composite with geometry, pLDDT, binding rank, annotation evidence, review priority, or decoy-selection variables. Export identity fraction, anchor score, and candidate-exposed score as secondary descriptive endpoints.

- [ ] **Step 2: Write failing tests for strict-set integrity**

Use a fixture with one target and five decoys. Assert that the evaluator:

- requires exactly the recorded five `strict_pass` decoys;
- rejects a target or any decoy with limited coverage;
- calculates `target_minus_decoy_median` and a deterministic within-set rank;
- leaves a target as `not_evaluable` when its strict set is incomplete rather than backfilling with a relaxed decoy.

Also test that two rows belonging to the same annotated evidence system are never treated as independent replicates in a pooled p-value.

- [ ] **Step 3: Implement row-level results, not an inflated cohort p-value**

For every target in `target_feasibility.csv`, emit one row containing eligibility, reason, decoy IDs, target score, decoy median/min/max, target rank among its six-member set, and secondary endpoint values. If fewer than three independent evaluable evidence systems exist, write `global_inference_status = insufficient_independent_systems` and no global p-value.

The evaluator may report an exact rank-based empirical tail fraction within an individual six-member set, labelled `descriptive_within_set_rank_fraction`; it must not call this a generalizable p-value.

- [ ] **Step 4: Run the evaluator and inspect its decision gate**

Run: `PYTHONPATH=src python3 src/run_register_aware_decoy_evaluation.py`

Expected with the current frozen inputs: an auditable explanation of which strict sets lack robust coverage or complete decoys. No selection relaxation, no global significance assertion, and no positive molecular-mimicry conclusion should appear.

- [ ] **Step 5: Update README interpretation language**

Add a short result-decision rule:

- A positive computational prioritization result needs an independently interpretable target, complete frozen strict decoys, robust score coverage, and persistence under the retained register-sensitivity analysis.
- Otherwise the result is negative/mixed: the current candidate/decoy evidence does not support robust register-aware prioritization. This is a valid methodological finding, not evidence against all EBV–myelin biology.

### Task 4: Generate paper-ready provenance and sensitivity reporting

**Files:**
- Create: `src/build_register_aware_reporting_tables.py`
- Create: `processed/register_aware_scoring/evidence_hierarchy_table.csv`
- Create: `processed/register_aware_scoring/register_sensitivity_appendix.csv`
- Create: `processed/register_aware_scoring/paper_result_branch.md`
- Create: `tests/test_register_aware_reporting_tables.py`

**Interfaces:**
- Consumes: the score table, decoy results, experimental override registry, and IEDB prediction summary.
- Produces: manuscript-ready factual tables and a decision-branch text draft, with no new analysis decision made in prose.

- [ ] **Step 1: Write failing provenance tests**

Assert that:

- 1BX2 MBP is presented as same-allele experimental-register support;
- 1H15 BALF5 is present but marked DRB5*01:01 calibration-only;
- gH and candidate-MBP records appear in the sensitivity appendix, not the primary result set;
- every row retains `register_source`, `register_assessment`, and the score/geometry eligibility statuses.

- [ ] **Step 2: Implement tables and decision-branch prose**

`evidence_hierarchy_table.csv` should compactly map every evidence category to eligible use, source, and claim boundary. `register_sensitivity_appendix.csv` should show excluded alternatives with the precise reason, rather than silently dropping them.

`paper_result_branch.md` must select exactly one generated branch from observed table statuses:

```text
positive_computational_prioritization
negative_or_mixed_method_result
```

It must state the preconditions/failures mechanically and must only use the permitted pMHC candidate-prioritization wording. Include a one-sentence prospective wet-lab paragraph clearly outside the present paper.

- [ ] **Step 3: Run complete verification and regenerate all artifacts**

Run:

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v
PYTHONPATH=src python3 src/build_register_aware_score_table.py
PYTHONPATH=src python3 src/run_register_aware_decoy_evaluation.py
PYTHONPATH=src python3 src/build_register_aware_reporting_tables.py
```

Expected: all tests pass; all tables have stable row counts and reproducible ordering; the final branch correctly reflects current data without biological overclaiming.

## Execution Status

Completed inline on 2026-08-11. All four tasks were implemented, artifacts
were regenerated from frozen inputs, and the complete test suite passed
(33 tests). The observed strict-decoy result is negative/mixed: no annotated
target has a complete, robustly scored five-decoy set under the predeclared
rules.

## Plan Self-Review

- **Coverage:** turns the approved paper design into frozen-input scoring, robust-coverage gating, strict-decoy evaluation, sensitivity reporting, and paper-ready outputs.
- **Geometry honesty:** explicitly prevents use of the existing aggregate local-alignment RMSD as a register-filtered structural score. Per-residue geometry becomes a later recover-and-recompute task, not an untracked assumption.
- **Statistical honesty:** preserves row-level decoy evidence but prevents a pooled p-value from correlated annotations or incomplete/limited-coverage strict sets.
- **Claim boundary:** no task can yield a TCR, cross-reactivity, established mimicry, or MS-mechanism claim; all positive language is constrained to computational pMHC prioritization.
- **Current expected outcome:** the plan is designed to produce a defensible negative/mixed result if robust decoy comparisons cannot be completed. That is informative and publishable as a methods/candidate-prioritization result only if positioned and validated appropriately.
