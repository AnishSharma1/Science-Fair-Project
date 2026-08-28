# Register-Aware pMHC Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a reproducible, register-aware pMHC candidate-prioritization benchmark that uses the full already-modeled Tier-1 EBV × myelin pair universe for matched decoys, while preserving conservative claim boundaries.

**Architecture:** The benchmark will reuse the existing IEDB-derived core hypotheses and local-alignment coordinates. A small, pure helper module will determine whether an already-aligned residue pair occupies the same P1–P9 index, and a builder will assemble a pre-score decoy background from every PASS row in the geometry matrix—not merely the 32 previously ranked rows. A report generator will emit auditable tables for targets, selected decoys, exclusions, and readiness without calculating a new biological effect claim.

**Tech Stack:** Python 3 standard library (`csv`, `collections`, `pathlib`, `unittest`); existing project CSV artifacts; no new package dependency.

## Global Constraints

- Scope remains HLA-DRB1*15:01, MHC-II, EBV versus validated myelin-source peptides.
- Use the complete pre-score geometry matrix only where `status == "PASS"`; never choose decoys using `review_priority_heuristic`, rank, local RMSD, or a register-aware score.
- The register source must be retained for every pair. IEDB `recommended_binding` calls remain computational hypotheses, not experimental presentation evidence.
- A pair is assessable only if both selected cores are fully contained, unique, and the original local alignment contains at least one same-register residue; otherwise record a reason for exclusion.
- Strict decoy matching is fixed at both EBV and human peptide lengths within one residue and zero IEDB binding-rank-bin mismatches. Composition distance and peptide pLDDT only order already-eligible decoys.
- Require five decoys for a benchmark-ready target. Preserve partial and zero-decoy targets explicitly; do not substitute random or unmatched controls.
- Output language must say “pMHC candidate prioritization” and must not claim shared TCR binding, activation, cross-reactivity, or an MS mechanism.
- This copied project bundle has no `.git` directory; record validation commands in the output report rather than attempting a commit.

---

## File Structure

- `src/register_aware_benchmark.py` — pure eligibility, covariate, and decoy-selection functions.
- `src/build_register_aware_benchmark.py` — reads existing artifacts and writes the complete-universe benchmark tables.
- `tests/test_register_aware_benchmark.py` — unit tests for register eligibility and strict selection boundaries.
- `processed/register_aware_benchmark/benchmark_pair_universe.csv` — every pre-score PASS pair plus an inclusion/exclusion decision.
- `processed/register_aware_benchmark/matched_decoy_sets.csv` — only predeclared, strict-matched decoys, selected deterministically.
- `processed/register_aware_benchmark/target_feasibility.csv` — one row per external/context target and its five-decoy readiness status.
- `processed/register_aware_benchmark/README.md` — methods, results counts, reproducibility command, and claim boundary.

### Task 1: Create the pure register-aware benchmark helpers

**Files:**
- Create: `src/register_aware_benchmark.py`
- Create: `tests/test_register_aware_benchmark.py`

**Interfaces:**
- Consumes: parsed local alignments as `list[tuple[int, str, int, str]]`; core start positions as `int`; record dictionaries containing the fields used by `premeeting_rigor.ordered_decoys`.
- Produces: `is_assessable_same_register_pair(alignment, ebv_core_start, human_core_start) -> bool` and `strict_eligible_decoys(target, candidates, limit=5) -> tuple[list[dict[str, object]], int]`.

- [ ] **Step 1: Write the failing tests**

```python
from register_aware_benchmark import is_assessable_same_register_pair, strict_eligible_decoys

def test_assessable_pair_requires_at_least_one_same_register_aligned_residue():
    alignment = [(4, "Y", 3, "H"), (7, "F", 6, "F")]
    assert is_assessable_same_register_pair(alignment, 4, 3) is True
    assert is_assessable_same_register_pair(alignment, 4, 5) is False

def test_strict_eligible_decoys_rejects_binding_bin_mismatches():
    target = {"pair_id": "target", "ebv_peptide": "ABCDEFGHIJK", "human_peptide": "LMNOPQRSTUV", "ebv_plddt": 80, "human_plddt": 80, "ebv_binding_rank": 1, "human_binding_rank": 1}
    mismatched = {**target, "pair_id": "bad", "ebv_binding_rank": 11}
    selected, available = strict_eligible_decoys(target, [mismatched])
    assert selected == []
    assert available == 0
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: FAIL because `register_aware_benchmark` does not exist.

- [ ] **Step 3: Write the minimal implementation**

```python
from premeeting_rigor import eligible_decoys
from register_aware_diagnostic import same_register_alignment_count

def is_assessable_same_register_pair(alignment, ebv_core_start, human_core_start):
    return bool(alignment) and same_register_alignment_count(
        alignment, ebv_core_start, human_core_start
    ) > 0

def strict_eligible_decoys(target, candidates, limit=5):
    return eligible_decoys(target, candidates, limit=limit)
```

- [ ] **Step 4: Run the focused test to verify it passes**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: PASS with two passing tests.

- [ ] **Step 5: Record the change**

Add the new files to the working copy. Do not run `git commit`: this bundle is not a Git repository.

### Task 2: Build the full pre-score pair universe with transparent exclusions

**Files:**
- Create: `src/build_register_aware_benchmark.py`
- Modify: `tests/test_register_aware_benchmark.py`
- Create: `processed/register_aware_benchmark/benchmark_pair_universe.csv`

**Interfaces:**
- Consumes: `processed/colabfold_tier1_ebv_myelin_geometry_matrix.csv`, `processed/pmhc_candidate_manifest.csv`, `processed/register_sensitivity/register_prediction_summary.csv`, and `processed/external_validation_panel.csv`.
- Produces: `build_pair_universe(geometry_rows, prediction_by_candidate, manifest_by_candidate, validation_groups_by_candidate=None) -> list[dict[str, object]]`, with columns `pair_id`, `pair_validation`, `register_assessment`, `same_register_alignment_count`, `ebv_binding_rank`, `human_binding_rank`, `ebv_plddt`, `human_plddt`, and `decoy_background_eligible`.

- [ ] **Step 1: Write the failing test**

```python
def test_build_pair_universe_excludes_zero_same_register_pairs_from_background():
    rows = build_pair_universe(
        geometry_rows=[{"status": "PASS", "ebv_candidate_id": "E", "human_candidate_id": "H", "aligned_positions_ebv_to_human": "4Y:3H", "ebv_peptide_mean_plddt": "80", "human_peptide_mean_plddt": "81"}],
        prediction_by_candidate={"E": {"predicted_core_start_positions_1_based": "4", "predicted_core_fully_contained_in_manifest_peptide": "True", "predicted_percentile_rank": "1"}, "H": {"predicted_core_start_positions_1_based": "5", "predicted_core_fully_contained_in_manifest_peptide": "True", "predicted_percentile_rank": "1"}},
        manifest_by_candidate={"E": {"peptide": "ABCDEFGHIJK"}, "H": {"peptide": "LMNOPQRSTUV"}},
    )
    assert rows[0]["register_assessment"] == "no_same_register_local_alignment"
    assert rows[0]["decoy_background_eligible"] is False
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: FAIL because `build_pair_universe` is not defined.

- [ ] **Step 3: Implement deterministic universe construction**

```python
if geometry_row["status"] != "PASS":
    continue
alignment = parse_local_alignment_positions(geometry_row["aligned_positions_ebv_to_human"])
if not unique_contained_core(ebv_prediction) or not unique_contained_core(human_prediction):
    assessment = "unresolved_or_flank_dependent_core"
elif not is_assessable_same_register_pair(alignment, ebv_start, human_start):
    assessment = "no_same_register_local_alignment"
else:
    assessment = "assessable_iedb_top_core_hypothesis"
```

Assign `decoy_background_eligible = True` only for `assessable_iedb_top_core_hypothesis` background rows. Derive `pair_validation` with the same explicit source/context rules currently used in `src/run_external_validation_benchmark.py`; do not infer a new validation label from the score.

- [ ] **Step 4: Build the universe and inspect its invariants**

Run: `PYTHONPATH=src python3 src/build_register_aware_benchmark.py`

Expected: creates `processed/register_aware_benchmark/benchmark_pair_universe.csv`; every row has a register assessment; no non-PASS geometry row is included; no priority-score field appears in the output.

- [ ] **Step 5: Run all unit tests**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v`

Expected: all existing 14 tests plus the new focused tests pass.

### Task 3: Select strict decoys from the expanded universe and calculate readiness

**Files:**
- Modify: `src/build_register_aware_benchmark.py`
- Modify: `tests/test_register_aware_benchmark.py`
- Create: `processed/register_aware_benchmark/matched_decoy_sets.csv`
- Create: `processed/register_aware_benchmark/target_feasibility.csv`

**Interfaces:**
- Consumes: `build_pair_universe()` output.
- Produces: `build_decoy_benchmark(pair_rows, target_decoy_count=5) -> tuple[list[dict[str, object]], list[dict[str, object]]]`. Every external/context target produces one feasibility row; a target that lacks an assessable top-core hypothesis has `readiness_status = "not assessable"` and zero selected decoys.

- [ ] **Step 1: Write the failing test**

```python
def test_decoy_selection_uses_only_assessable_background_rows():
    target = assessable_target_row()
    valid = assessable_background_row("valid")
    invalid = {**assessable_background_row("invalid"), "decoy_background_eligible": False}
    decoys, feasibility = build_decoy_benchmark([target, valid, invalid], target_decoy_count=1)
    assert [row["decoy_pair_id"] for row in decoys] == ["valid"]
    assert feasibility[0]["eligible_decoy_count"] == 1
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: FAIL because `build_decoy_benchmark` is not defined.

- [ ] **Step 3: Implement strict, score-blind selection**

```python
background = [row for row in pair_rows if row["pair_validation"] == "background" and row["decoy_background_eligible"]]
targets = [row for row in pair_rows if row["pair_validation"] != "background"]
for target in targets:
    if target["register_assessment"] != "assessable_iedb_top_core_hypothesis":
        write_not_assessable_feasibility_row(target)
        continue
    selected, available = strict_eligible_decoys(target, background, limit=target_decoy_count)
```

Write `readiness_status` as `ready` only if `available >= 5`, `partial` if `1 <= available < 5`, `no eligible decoy` if `available == 0` after an assessable register call, and `not assessable` for an unresolved/flank-dependent core or zero same-register alignment. Store the exact matching rule and the selection boundary on every output row.

- [ ] **Step 4: Regenerate the tables**

Run: `PYTHONPATH=src python3 src/build_register_aware_benchmark.py`

Expected: decoy rows contain no score/rank field; feasibility counts are calculated from the full PASS, assessable background rather than the old 32-pair shortlist.

- [ ] **Step 5: Run all unit tests**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v`

Expected: all tests pass.

### Task 4: Produce a decision-ready reproducibility report

**Files:**
- Modify: `src/build_register_aware_benchmark.py`
- Create: `processed/register_aware_benchmark/README.md`
- Modify: `docs/meeting_prep/2026-08-10-revised-claim-language.md`

**Interfaces:**
- Consumes: benchmark-universe, decoy-set, and feasibility outputs from Tasks 2–3.
- Produces: a short Markdown report that states counts, the exact frozen analysis boundary, unresolved mentor decisions, and the correct claim limitation.

- [ ] **Step 1: Write the failing report-content test**

```python
def test_report_does_not_overstate_tcr_or_disease_evidence():
    report = render_benchmark_readme(summary={"ready_targets": 0, "partial_targets": 1})
    assert "does not establish shared-TCR binding" in report
    assert "MS mechanism" in report
    assert "ready_targets: 0" in report
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_register_aware_benchmark.py' -v`

Expected: FAIL because `render_benchmark_readme` is not defined.

- [ ] **Step 3: Implement the report and update meeting language**

```python
def render_benchmark_readme(summary):
    return (
        "# Register-aware matched-decoy benchmark\\n\\n"
        f"- ready_targets: {summary['ready_targets']}\\n"
        f"- partial_targets: {summary['partial_targets']}\\n\\n"
        "This is pMHC candidate prioritization only; it does not establish "
        "shared-TCR binding, T-cell activation, cross-reactivity, or an MS mechanism.\\n"
    )
```

Document that IEDB core calls are provisional, include the reproducibility command, and identify only these remaining mentor decisions: the final register hierarchy, BALF5–MBP’s DR15-haplotype framing, and which positions may be described as pocket-facing versus candidate TCR-exposed.

- [ ] **Step 4: Regenerate and read the report**

Run: `PYTHONPATH=src python3 src/build_register_aware_benchmark.py && sed -n '1,220p' processed/register_aware_benchmark/README.md`

Expected: report has no TCR/disease claim, names the complete-universe strict matching rule, and gives an auditable readiness count.

- [ ] **Step 5: Final verification**

Run: `PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v`

Expected: full suite passes. Record the output in the final handoff because a Git commit is unavailable in this bundle.

## Plan Self-Review

- Coverage: implements the predeclared register rule, full-universe expansion, strict score-blind matched-decoy selection, transparent infeasibility reporting, and meeting-safe claim language.
- Scope: intentionally omits a new structural run, a new predictor, TCR modeling, experimental conclusions, and any relaxation of the strict matching rule.
- Consistency: all selection functions use the existing `premeeting_rigor.eligible_decoys` constraints; all new output derives from pre-score PASS geometry plus stored IEDB hypotheses.
