# Central Model Library Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a centralized, method-separated EBV-MS model library that preserves every original model file.

**Architecture:** The library uses lightweight aliases to authoritative source folders rather than copied structural files. Markdown READMEs give a human-readable map, and TSV inventories provide an inspectable, machine-readable catalog. Original source paths remain unchanged.

**Tech Stack:** macOS Finder aliases, shell inventory commands, Markdown, TSV.

## Global Constraints

- Never move, rename, delete, or modify source model folders.
- Keep AF3 pMHC, legacy ColabFold pMHC, and TCRmodel2 ternary outputs in distinct sections.
- Record external source paths exactly and retain source basename in all inventories.
- Do not interpret structural confidence or geometry while organizing files.
- Validate aliases and every generated inventory after creation.

---

### Task 1: Inventory sources and create library skeleton

**Files:**
- Create: `Model_Library/00_README.md`
- Create: `Model_Library/01_AF3_pMHC/README.md`
- Create: `Model_Library/02_ColabFold_pMHC_Legacy/README.md`
- Create: `Model_Library/03_TCR_pMHC_TCRmodel2/README.md`
- Create: `Model_Library/04_Experimental_References/README.md`
- Create: `Model_Library/05_Analysis_and_Manifests/README.md`

- [ ] **Step 1: Confirm all authoritative source directories exist**

Run:

```bash
test -d "$AF3_DOWNLOAD_ROOT" && test -d "$COLABFOLD_ROOT" && test -d "$TCRMODEL2_ROOT"
```

Expected: success only when all three source roots are present.

- [ ] **Step 2: Create the six section directories and READMEs**

Use exact names from the design document. State method, scope, authoritative source, and the no-move preservation rule in each README.

- [ ] **Step 3: Validate skeleton**

Run:

```bash
find "$LIBRARY_ROOT" -maxdepth 2 -name README.md -type f | wc -l
```

Expected: six README files.

### Task 2: Create method-specific aliases and inventories

**Files:**
- Create: `Model_Library/01_AF3_pMHC/af3_download_inventory.tsv`
- Create: `Model_Library/02_ColabFold_pMHC_Legacy/legacy_colabfold_inventory.tsv`
- Create: `Model_Library/03_TCR_pMHC_TCRmodel2/tcrmodel2_inventory.tsv`

- [ ] **Step 1: Create Finder aliases to each authoritative source folder**

Create four AF3 download aliases, one AF3 submission-batch alias, one legacy ColabFold alias, one TCRmodel2 model-package alias, and one experimental-reference alias. The alias filename must express the method and source role.

- [ ] **Step 2: Generate the AF3 inventory**

Use one TSV row per `folds_2026_08_10_*` job directory, with columns `method`, `model_family`, `source_download_folder`, `job_name`, `source_path`, and `expected_models`. Set `method=AlphaFold3` and `model_family=pMHC`; write `expected_models=5` only when five model CIF files exist, otherwise write the observed count.

- [ ] **Step 3: Generate the legacy ColabFold inventory**

Use one TSV row per recognized top-level legacy artifact or artifact directory. Use `method=ColabFold`, `model_family=pMHC`, and classify entries as `input`, `QA`, `geometry_triage`, `template_transfer`, or `scaffold_exploration` from their folder/name.

- [ ] **Step 4: Generate the TCRmodel2 inventory**

Use one TSV row per `tcrmodel2_results` run and one row per reference structure. Use `method=TCRmodel2`, `model_family=TCR-pMHC`, report the five ranked PDB-model count, and classify calibration/reference outputs separately from hypothesis runs.

- [ ] **Step 5: Validate aliases and inventory headers**

Open each alias target programmatically and check every TSV contains a header plus at least one data row. Do not follow aliases with a destructive command.

### Task 3: Create master catalog and final preservation check

**Files:**
- Create: `Model_Library/05_Analysis_and_Manifests/master_model_inventory.tsv`

- [ ] **Step 1: Build one master row per method-specific inventory**

Columns: `library_section`, `method`, `model_family`, `inventory_file`, `source_root`, `notes`. Include AF3 pMHC, legacy ColabFold pMHC, TCRmodel2 TCR-pMHC, and experimental-reference rows.

- [ ] **Step 2: Write the top-level README navigation**

Document section names, the distinction between pMHC-only and TCR–pMHC predictions, and the rule that the library is a catalog—not a new source of scientific claims.

- [ ] **Step 3: Verify no source path was modified**

Compare the pre-build and post-build directory counts for all source roots and confirm that generated content lives only under `Model_Library` and the planning/spec documentation paths.

- [ ] **Step 4: Final spot-check**

Print the `Model_Library` tree to depth two, inspect the first three rows of each TSV, and resolve every created alias. Expected result: six sections, working aliases, nonempty inventories, and unchanged source-folder counts.

## Self-Review

- Spec coverage: Tasks 1–3 create the agreed library, aliases, all inventories, master manifest, navigation, and preservation verification.
- Placeholder scan: no incomplete actions remain; all generated artifacts and source families are named explicitly.
- Consistency: all model-specific inventories feed the master catalog, and all sections use the names from the approved design.
