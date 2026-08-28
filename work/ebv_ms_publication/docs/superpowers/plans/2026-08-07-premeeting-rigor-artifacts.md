# Pre-meeting Rigor Artifacts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create auditable HLA-II register, decoy-readiness, and cluster-aware validation artifacts before the Yicong Liu meeting.

**Architecture:** A pure-Python library owns deterministic scientific bookkeeping. A separate executable obtains official IEDB predictions, maps them to the existing candidate manifest, and writes derivative CSV/Markdown files without altering the original ranking or structural outputs.

**Tech Stack:** Python 3 standard library; official IEDB MHC-II REST API; CSV and TSV files.

## Global Constraints

- Preserve every pre-existing `processed/` input unchanged.
- All output paths must be under `processed/register_sensitivity/`, `processed/matched_decoys/`, or `processed/validation_hygiene/`.
- Never use `review_priority_heuristic` to select or order decoys.
- Mark prediction-derived fields as computational hypotheses.
- Do not make TCR-recognition, cross-reactivity, activation, or disease-causation claims.

---

### Task 1: Pure rigor utilities

**Files:**

- Create: `tests/test_premeeting_rigor.py`
- Create: `src/premeeting_rigor.py`

**Interfaces:**

- Produces `enumerate_core_windows(peptide: str) -> list[dict[str, object]]`.
- Produces `parse_iedb_mhcii_tsv(text: str) -> list[dict[str, str]]`.
- Produces `binding_rank_bin(rank: float) -> str`.
- Produces `composition_distance(left: str, right: str) -> float`.
- Produces `ordered_decoys(target: dict[str, object], candidates: list[dict[str, object]], limit: int) -> list[dict[str, object]]`.

- [ ] Write failing tests for exhaustive window enumeration, IEDB TSV parsing, and decoy ordering without priority-score access.
- [ ] Run `python3 -m unittest tests/test_premeeting_rigor.py -v` and verify failure before creating the library.
- [ ] Implement the minimum pure functions.
- [ ] Rerun the test file until it passes.

### Task 2: Reproducible artifact generator

**Files:**

- Create: `src/build_premeeting_rigor_artifacts.py`
- Creates: `processed/register_sensitivity/iedb_mhcii_drb1501_raw.tsv`
- Creates: `processed/register_sensitivity/register_prediction_summary.csv`
- Creates: `processed/register_sensitivity/register_window_catalog.csv`
- Creates: `processed/matched_decoys/decoy_readiness.csv`
- Creates: `processed/validation_hygiene/validation_evidence_clusters.csv`
- Creates: `processed/validation_hygiene/README.md`

- [ ] Write a failing test for mapped prediction rows and rank-bin matching.
- [ ] Run it before writing the executable.
- [ ] Implement official IEDB retrieval, raw-response preservation, manifest mapping, decoy generation, and cluster summaries.
- [ ] Run the generator and full unittest suite.
- [ ] Verify all rows map to source candidates and no decoy output contains the priority heuristic.

### Task 3: Cluster-aware validation wording and meeting integration

**Files:**

- Modify: `src/run_external_validation_benchmark.py`
- Modify: `docs/meeting_prep/2026-08-10-yicong-liu-brief.md`
- Modify: `docs/meeting_prep/2026-08-10-register-aware-validation-protocol.md`

- [ ] Revise wording so BALF5--MBP is one calibration system and newer source overlays are never called direct positive pairs.
- [ ] Link new outputs in the meeting packet.
- [ ] Run the benchmark, generator, and full test suite.
