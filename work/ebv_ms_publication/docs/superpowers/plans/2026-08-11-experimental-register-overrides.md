# Experimental Register Override Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add source-traceable experimental-register and exclusion decisions to the benchmark without changing stored IEDB predictor outputs.

**Architecture:** A CSV registry will hold one reviewed decision per relevant candidate. The benchmark builder will resolve each candidate through that registry before falling back to its IEDB top-core hypothesis; an experimental same-allele register may be assessed, while DRB5 calibration and unresolved sensitivity candidates are made explicitly ineligible. The generated report will show this hierarchy and retain the existing pMHC-only claim boundary.

**Tech Stack:** Python 3 standard library (`csv`, `pathlib`, `unittest`); existing benchmark files; no new dependencies.

## Global Constraints

- Never overwrite `processed/register_sensitivity/register_prediction_summary.csv`; it remains the raw IEDB predictor record.
- `HUMAN_MYELIN_13572` uses the exact 1BX2 DRB1*15:01 reference core `VHFFKNIVT`, start 5.
- `EBV_TCELL_63843` uses `YHFVKKHVH`, start 5, only as a DRB5*01:01 calibration record; it cannot enter a DRB1*15:01 same-allele benchmark.
- `EBV_TCELL_2268683` and `HUMAN_MYELIN_115641` are sensitivity-only because exact same-allele experimental registers were not established.
- Do not make TCR binding, activation, cross-reactivity, molecular-mimicry, or MS-mechanism claims.
- No Git repository exists in this copied bundle; leave reproducibility commands and passing test results rather than committing.

---

### Task 1: Register-decision registry and resolver

**Files:**
- Create: `processed/register_sensitivity/experimental_register_overrides.csv`
- Modify: `src/build_register_aware_benchmark.py`
- Modify: `tests/test_register_aware_benchmark.py`

**Interfaces:**
- Consumes: a registry keyed by `candidate_id`, a manifest peptide, and the stored IEDB prediction row.
- Produces: `resolve_candidate_register(candidate_id, peptide, prediction, overrides) -> dict[str, object]` with `register_status`, `register_source`, `core_start_1_based`, and `core_peptide`.

- [ ] **Step 1: Write failing resolver tests**

```python
mbp = resolve_candidate_register("HUMAN_MYELIN_13572", "ENPVVHFFKNIVTPR", prediction, overrides)
assert mbp["register_status"] == "experimental_primary_allele_reference"
assert mbp["core_start_1_based"] == 5

bal = resolve_candidate_register("EBV_TCELL_63843", "TGGVYHFVKKHVHES", prediction, overrides)
assert bal["register_status"] == "calibration_only_nonprimary_allele"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: FAIL because `resolve_candidate_register` and the registry loader do not exist.

- [ ] **Step 3: Create the registry and resolver**

Create rows for 1BX2 MBP, 1H15 BALF5 calibration, and the gH/MBP sensitivity-only exclusions. Validate that a primary-allele core is a 9-mer contained at the recorded one-based start in the manifest peptide. Fall back to the stored IEDB call only when no override exists.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: PASS, including source-priority and nonprimary-calibration tests.

### Task 2: Apply registry decisions to the universe and regenerate outputs

**Files:**
- Modify: `src/build_register_aware_benchmark.py`
- Modify: `tests/test_register_aware_benchmark.py`
- Modify: `processed/register_aware_benchmark/benchmark_pair_universe.csv`
- Modify: `processed/register_aware_benchmark/target_feasibility.csv`
- Modify: `processed/register_aware_benchmark/matched_decoy_sets.csv`
- Modify: `processed/register_aware_benchmark/README.md`
- Modify: `docs/meeting_prep/2026-08-10-revised-claim-language.md`

**Interfaces:**
- Consumes: `resolve_candidate_register` records for each geometry pair.
- Produces: pair rows with per-arm source/status fields and explicit `calibration_only_nonprimary_allele` or `sensitivity_only_unresolved_register` exclusion statuses.

- [ ] **Step 1: Write failing universe test**

```python
rows = build_pair_universe(geometry, predictions, manifest, overrides_by_candidate=overrides)
assert rows[0]["register_assessment"] == "calibration_only_nonprimary_allele"
assert rows[0]["decoy_background_eligible"] is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: FAIL because `build_pair_universe` does not accept override records.

- [ ] **Step 3: Apply decisions before positional comparison**

```python
if "sensitivity_only_unresolved" in {ebv["register_status"], human["register_status"]}:
    assessment = "sensitivity_only_unresolved_register"
elif "calibration_only_nonprimary_allele" in {ebv["register_status"], human["register_status"]}:
    assessment = "calibration_only_nonprimary_allele"
else:
    assessment = compare_same_register_positions(...)
```

Write both arms’ source/status/core fields. A pair can become decoy-background eligible only after the primary-allele register comparison succeeds.

- [ ] **Step 4: Regenerate and inspect the outputs**

Run: `PYTHONPATH=src python3 src/build_register_aware_benchmark.py`

Expected: BALF5 pairs show calibration-only status; gH and candidate-MBP pairs show sensitivity-only status; MBP 13572 identifies PDB 1BX2 as its register source.

- [ ] **Step 5: Full verification**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v`

Expected: all tests pass and the report includes the override hierarchy without strengthening biological claims.

## Plan Self-Review

- Coverage: separates experimental same-allele, cross-allotype calibration, unresolved sensitivity, and predictor-only records without deleting provenance.
- Scope: excludes new structure prediction, new external database collection, TCR modeling, and changes to the score itself.
- Consistency: only same-allele experimental overrides can support a primary-screen register assessment; all other reviewed records remain visible but excluded from decoy eligibility.
