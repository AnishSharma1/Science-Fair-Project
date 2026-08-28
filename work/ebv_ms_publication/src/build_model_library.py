#!/usr/bin/env python3
"""Build a non-destructive catalog for EBV-MS structural-model artifacts."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT = Path("/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads/ebv_ms_publication")
DOWNLOADS = Path("/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads")
AF3_ROOTS = sorted(path for path in DOWNLOADS.glob("folds_2026_08_10_*") if path.is_dir())
COLABFOLD_ROOT = PROJECT / "processed"
AF3_SUBMISSIONS = PROJECT / "af3_migration_2026-08-10"
TCR_PACKAGE = Path("/Users/anishsharma/Documents/New project/outputs/ebv_ms_model_package")
LIBRARY = PROJECT / "Model_Library"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, content: str) -> None:
    path.write_text(content.rstrip() + "\n", encoding="utf-8")


def add_link(section: Path, name: str, target: Path) -> None:
    link = section / name
    if link.is_symlink() or link.exists():
        if link.resolve() != target.resolve():
            raise RuntimeError(f"Refusing to replace unexpected existing item: {link}")
        return
    link.symlink_to(target)


def build_readmes(sections: dict[str, Path]) -> None:
    write_text(
        LIBRARY / "00_README.md",
        """# EBV-MS Model Library

This is a centralized catalog. The model files remain in their original locations; each folder link points to an authoritative source.

## Sections

- `01_AF3_pMHC`: new AlphaFold 3 three-chain pMHC downloads and submission batches.
- `02_ColabFold_pMHC_Legacy`: legacy ColabFold pMHC inputs, QA, and exploratory structure artifacts.
- `03_TCR_pMHC_TCRmodel2`: five-component TCR–pMHC TCRmodel2 predictions and calibration runs.
- `04_Experimental_References`: experimental PDB reference structures used for calibration.
- `05_Analysis_and_Manifests`: the cross-method catalog.

The pMHC-only folders do not contain TCR-docking results. This catalog organizes files only; it does not add scientific conclusions.""",
    )
    write_text(
        sections["af3"] / "README.md",
        f"""# AlphaFold 3 pMHC

Authoritative raw downloads are linked here from `{DOWNLOADS}`. The job inventory contains one row per downloaded job and reports the observed number of model CIF files. Submission batches are a separate link so requests are not confused with results.

All models are three-chain pMHC predictions (HLA-DRA, HLA-DRB1*15:01, peptide), not TCR–pMHC docking predictions. Originals are preserved.""",
    )
    write_text(
        sections["colabfold"] / "README.md",
        f"""# Legacy ColabFold pMHC

This section links to the prior ColabFold pMHC workflow in `{COLABFOLD_ROOT}`. Its inventory separates inputs, QA, geometry triage, template transfer, and scaffold exploration so older exploratory work is not mistaken for new AF3 output.

Originals are preserved.""",
    )
    write_text(
        sections["tcr"] / "README.md",
        f"""# TCR–pMHC TCRmodel2

This section links to `{TCR_PACKAGE}`. These are TCRmodel2 ternary predictions and calibration outputs, not ColabFold models. Each run is listed with its observed ranked-PDB count and its calibration or hypothesis role.

Originals are preserved.""",
    )
    write_text(
        sections["references"] / "README.md",
        """# Experimental References

This section links to the experimental PDB reference structures used in the TCRmodel2 calibration work. They are references, not prediction outputs. Originals are preserved.""",
    )
    write_text(
        sections["manifests"] / "README.md",
        """# Analysis and Manifests

`master_model_inventory.tsv` maps each model family to its method-specific inventory and authoritative source root. Use the method-specific inventories for individual files/runs; use this table to navigate the full library.""",
    )


def af3_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for root in AF3_ROOTS:
        for job in sorted(path for path in root.iterdir() if path.is_dir()):
            model_count = len(list(job.glob("*_model_*.cif")))
            rows.append(
                {
                    "method": "AlphaFold3",
                    "model_family": "pMHC",
                    "source_download_folder": root.name,
                    "job_name": job.name,
                    "source_path": str(job),
                    "expected_models": str(model_count),
                }
            )
    return rows


def legacy_rows() -> list[dict[str, str]]:
    artifacts = [
        ("colabfold_inputs", "input"),
        ("pmhc_colabfold_inputs.fasta", "input"),
        ("pmhc_colabfold_batch.csv", "input"),
        ("pmhc_colabfold_manifest.csv", "input"),
        ("pmhc_colabfold_metadata.json", "input"),
        ("colabfold_pmhc_peptide_qa.csv", "QA"),
        ("colabfold_pmhc_peptide_residue_qa.csv", "QA"),
        ("pmhc_structure_qa.csv", "QA"),
        ("pmhc_structure_qa_summary.csv", "QA"),
        ("colabfold_pmhc_pair_triage.csv", "geometry_triage"),
        ("colabfold_tier1_ebv_myelin_geometry_matrix.csv", "geometry_triage"),
        ("ob1a12_template_transfer", "template_transfer"),
        ("ob1a12_ebv_scaffold_test", "scaffold_exploration"),
        ("ob1a12_ternary_colabfold_inputs.fasta", "input"),
        ("ob1a12_ternary_evaluation", "geometry_triage"),
    ]
    rows: list[dict[str, str]] = []
    for relative, role in artifacts:
        source = COLABFOLD_ROOT / relative
        if source.exists():
            rows.append(
                {
                    "method": "ColabFold",
                    "model_family": "pMHC",
                    "artifact_role": role,
                    "artifact_name": source.name,
                    "source_path": str(source),
                }
            )
    return rows


def tcr_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    results = TCR_PACKAGE / "tcrmodel2_results"
    for run in sorted(path for path in results.iterdir() if path.is_dir()):
        rows.append(
            {
                "method": "TCRmodel2",
                "model_family": "TCR-pMHC",
                "artifact_role": "calibration" if run.name.startswith("CAL_") else "hypothesis_run",
                "run_or_reference": run.name,
                "source_path": str(run),
                "ranked_pdb_count": str(len(list(run.glob("ranked_*.pdb")))),
            }
        )
    references = TCR_PACKAGE / "reference_structures"
    for pdb in sorted(references.glob("*.cif")):
        rows.append(
            {
                "method": "experimental_reference",
                "model_family": "TCR-pMHC",
                "artifact_role": "experimental_reference",
                "run_or_reference": pdb.stem,
                "source_path": str(pdb),
                "ranked_pdb_count": "",
            }
        )
    return rows


def main() -> None:
    if len(AF3_ROOTS) != 4:
        raise RuntimeError(f"Expected four AF3 download folders, found {len(AF3_ROOTS)}")
    for source in [COLABFOLD_ROOT, AF3_SUBMISSIONS, TCR_PACKAGE]:
        if not source.is_dir():
            raise RuntimeError(f"Missing authoritative source directory: {source}")

    sections = {
        "af3": LIBRARY / "01_AF3_pMHC",
        "colabfold": LIBRARY / "02_ColabFold_pMHC_Legacy",
        "tcr": LIBRARY / "03_TCR_pMHC_TCRmodel2",
        "references": LIBRARY / "04_Experimental_References",
        "manifests": LIBRARY / "05_Analysis_and_Manifests",
    }
    for section in [LIBRARY, *sections.values()]:
        section.mkdir(parents=True, exist_ok=True)
    build_readmes(sections)

    for index, source in enumerate(AF3_ROOTS, start=1):
        add_link(sections["af3"], f"AF3_download_{index:02d}", source)
    add_link(sections["af3"], "AF3_submission_batches", AF3_SUBMISSIONS)
    add_link(sections["colabfold"], "Legacy_ColabFold_pMHC_source", COLABFOLD_ROOT)
    add_link(sections["tcr"], "TCRmodel2_model_package", TCR_PACKAGE)
    add_link(sections["references"], "Experimental_reference_structures", TCR_PACKAGE / "reference_structures")

    write_tsv(
        sections["af3"] / "af3_download_inventory.tsv",
        ["method", "model_family", "source_download_folder", "job_name", "source_path", "expected_models"],
        af3_rows(),
    )
    write_tsv(
        sections["colabfold"] / "legacy_colabfold_inventory.tsv",
        ["method", "model_family", "artifact_role", "artifact_name", "source_path"],
        legacy_rows(),
    )
    write_tsv(
        sections["tcr"] / "tcrmodel2_inventory.tsv",
        ["method", "model_family", "artifact_role", "run_or_reference", "source_path", "ranked_pdb_count"],
        tcr_rows(),
    )
    write_tsv(
        sections["manifests"] / "master_model_inventory.tsv",
        ["library_section", "method", "model_family", "inventory_file", "source_root", "notes"],
        [
            {"library_section": "01_AF3_pMHC", "method": "AlphaFold3", "model_family": "pMHC", "inventory_file": "../01_AF3_pMHC/af3_download_inventory.tsv", "source_root": str(DOWNLOADS), "notes": "Four download folders; submissions linked separately."},
            {"library_section": "02_ColabFold_pMHC_Legacy", "method": "ColabFold", "model_family": "pMHC", "inventory_file": "../02_ColabFold_pMHC_Legacy/legacy_colabfold_inventory.tsv", "source_root": str(COLABFOLD_ROOT), "notes": "Legacy inputs, QA, triage, and explorations."},
            {"library_section": "03_TCR_pMHC_TCRmodel2", "method": "TCRmodel2", "model_family": "TCR-pMHC", "inventory_file": "../03_TCR_pMHC_TCRmodel2/tcrmodel2_inventory.tsv", "source_root": str(TCR_PACKAGE), "notes": "Ternary predictions and calibration runs."},
            {"library_section": "04_Experimental_References", "method": "experimental_reference", "model_family": "TCR-pMHC", "inventory_file": "../03_TCR_pMHC_TCRmodel2/tcrmodel2_inventory.tsv", "source_root": str(TCR_PACKAGE / "reference_structures"), "notes": "Experimental PDB references, not predictions."},
        ],
    )


if __name__ == "__main__":
    main()
