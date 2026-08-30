"""Prepare, calculate, and analyze the additive pMHC electrostatics pilot.

The workflow is deliberately limited to two frozen BALF5--TALDO1 leads and
their already frozen, score-blind HLA-specific N3 panels. It never modifies
V1--V3 rankings and never converts computational resemblance into a
specificity, cross-reactivity, molecular-mimicry, or MS-mechanism claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from analyze_af3_pmhc_downloads import parse_mmcif
from pmhc_surface_electrostatics import (
    APBSParameters,
    DIELECTRIC_VALUES,
    GridSpec,
    align_model_to_reference,
    assign_electrostatic_context,
    build_apbs_input,
    build_common_accessible_field_patch,
    build_shared_grid,
    candidate_exposed_histidines,
    carbo_similarity,
    dielectric_robustness,
    hodgkin_similarity,
    parse_open_dx,
    potential_rmse,
    rank_panel,
    sign_agreement_fraction,
    summarize_electrostatic_ensemble,
    surface_patch_points,
    trilinear_sample,
    write_model_pdb,
)


ROOT = Path(__file__).resolve().parents[1]
V3_DIR = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
N3_DIR = ROOT / "processed/high_yield_control_validation_2026-08-28"
V2_RESULTS_DIR = ROOT / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26"
DEFAULT_OUT = ROOT / "processed/pmhc_surface_electrostatics_pilot_2026-08-29"
TARGET_IDS = ("HY13_SEQ_02", "HY15_SEQ_02")
SEED = 271828
PDB2PQR = Path(
    "/Users/anishsharma/.cache/ebv_ms_tools/pmhc_electrostatics/"
    "pdb2pqr-3.7.1-py313-venv/bin/pdb2pqr"
)
APBS_INSTALL = Path(
    "/Users/anishsharma/.cache/ebv_ms_tools/pmhc_electrostatics/"
    "apbs-3.4.1-linux/APBS-3.4.1.Linux"
)
APBS_ARCHIVE = Path(
    "/Users/anishsharma/.cache/ebv_ms_tools/pmhc_electrostatics/"
    "APBS-3.4.1.Linux.zip"
)
APBS_IMAGE = "ubuntu:22.04"
APBS_IMAGE_DIGEST = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
CLAIM_BOUNDARY = (
    "Descriptive, model-derived HLA-specific pMHC local-field resemblance only; "
    "not evidence of presentation, TCR recognition, activation, specificity, "
    "cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate."
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fields or sorted({key for row in rows for key in row}))
    if not names:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _safe_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def _relative(path: Path, root: Path) -> str:
    return str(path.relative_to(root))


def validate_panel_rows(panel_id: str, rows: Sequence[Mapping[str, Any]]) -> None:
    if len(rows) != 26:
        raise ValueError(f"{panel_id} must contain exactly 26 rows")
    if sum(str(row.get("row_role")) == "target" for row in rows) != 1:
        raise ValueError(f"{panel_id} must contain exactly one target row")
    if sum(str(row.get("row_role")) == "n3" for row in rows) != 25:
        raise ValueError(f"{panel_id} must contain exactly 25 N3 rows")
    if len({str(row["pair_id"]) for row in rows}) != 26:
        raise ValueError(f"{panel_id} contains duplicate pair IDs")
    alleles = {str(row["allele"]) for row in rows}
    if len(alleles) != 1:
        raise ValueError(f"{panel_id} is not exact-HLA separated")


def parse_propka_histidines(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    pattern = re.compile(r"^\s+HIS\s+(\d+)\s+C\s+(-?\d+(?:\.\d+)?)\s+6\.50\s*$")
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            rows.append(
                {
                    "sequence_position_1_based": int(match.group(1)),
                    "predicted_pka": float(match.group(2)),
                }
            )
    unique = {(row["sequence_position_1_based"], row["predicted_pka"]): row for row in rows}
    return [unique[key] for key in sorted(unique)]


def build_lead_gate(
    *,
    target_id: str,
    rank: int | None,
    register_qc: bool,
    model_qc: bool,
    dielectric_robust: bool,
) -> dict[str, Any]:
    status = assign_electrostatic_context(rank, register_qc, model_qc)
    rank_only_context = (
        "not_evaluable" if rank is None
        else "electrostatic_context_supportive" if int(rank) <= 3
        else "electrostatic_context_not_supportive"
    )
    return {
        "target_id": target_id,
        "status": status,
        "rank_only_context": rank_only_context,
        "primary_full_pmhc_hodgkin_rank": rank,
        "register_qc_complete": bool(register_qc),
        "model_and_calculation_qc_complete": bool(model_qc),
        "rank_class_robust_across_solute_dielectrics_2_4_8": bool(dielectric_robust),
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "cross_reactivity_claim_allowed": False,
        "molecular_mimicry_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }


def _checksums(output_dir: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(item for item in output_dir.rglob("*") if item.is_file()):
        if path.name == "SHA256SUMS.csv" or path.suffix == ".dx":
            continue
        rows.append({"relative_path": _relative(path, output_dir), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    write_csv(output_dir / "SHA256SUMS.csv", rows, ("relative_path", "sha256", "bytes"))
    return rows


def _tool_version(command: Sequence[str]) -> str:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    return (completed.stdout or completed.stderr).splitlines()[0].strip() if completed.returncode == 0 else "unavailable"


def _panel_rows() -> dict[str, list[dict[str, str]]]:
    rows = read_csv(N3_DIR / "panel_feature_matrix.csv")
    panels = {target_id: [row for row in rows if row["target_id"] == target_id] for target_id in TARGET_IDS}
    for target_id, panel in panels.items():
        validate_panel_rows(target_id, panel)
    return panels


def _arm_registry(panel_id: str, rows: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    arms: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        for arm_class, prefix in (("ebv", "ebv"), ("self", "self")):
            key = (arm_class, row[f"{prefix}_candidate_id"])
            arm = {
                "panel_id": panel_id,
                "allele": row["allele"],
                "arm_class": arm_class,
                "candidate_id": row[f"{prefix}_candidate_id"],
                "protein": row[f"{prefix}_protein"],
                "sequence": row[f"{prefix}_sequence"],
                "core": row[f"{prefix}_core_p1_p9"],
                "core_start_1_based": int(row[f"{prefix}_declared_core_start_1_based"]),
                "model_count": int(row["left_model_count" if arm_class == "ebv" else "right_model_count"]),
            }
            if key in arms and arms[key] != arm:
                raise ValueError(f"inconsistent arm record for {panel_id} {key}")
            arms[key] = arm
    result = [arms[key] for key in sorted(arms)]
    if len([row for row in result if row["arm_class"] == "ebv"]) != 6:
        raise ValueError(f"{panel_id} does not contain six unique EBV arms")
    if len([row for row in result if row["arm_class"] == "self"]) != 6:
        raise ValueError(f"{panel_id} does not contain six unique self arms")
    return result


def _heavy_bounds(model: Mapping[str, Sequence[Mapping[str, Any]]]) -> tuple[np.ndarray, np.ndarray]:
    atoms = np.asarray(
        [
            atom["xyz"]
            for chain in model.values()
            for residue in chain
            for atom in residue["atoms"]
            if str(atom.get("element", "")).upper() != "H"
        ],
        dtype=float,
    )
    return atoms.min(axis=0), atoms.max(axis=0)


def prepare(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite existing package: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    panels = _panel_rows()
    frozen_targets = [row for row in read_csv(N3_DIR / "frozen_target_registry.csv") if row["target_id"] in TARGET_IDS]
    if {row["target_id"] for row in frozen_targets} != set(TARGET_IDS):
        raise ValueError("the two frozen lead records were not recovered exactly")
    write_csv(output_dir / "frozen_target_registry.csv", frozen_targets)
    frozen_pairs = [row for target_id in TARGET_IDS for row in panels[target_id]]
    write_csv(output_dir / "frozen_panel_pairs.csv", frozen_pairs)

    upstream = {
        "v3_checksums_sha256": sha256_file(V3_DIR / "SHA256SUMS.csv"),
        "n3_package_checksums_sha256": sha256_file(N3_DIR / "SHA256SUMS.csv"),
    }
    protocol = {
        "protocol_version": "pmhc_surface_electrostatics_pilot_v1",
        "frozen_at": "2026-08-29",
        "target_ids": list(TARGET_IDS),
        "target_count": 2,
        "panel_definition": "one frozen target plus 25 previously score-blind exact-HLA N3 pairs",
        "n3_recognition_status": "unknown_not_specificity_negative",
        "seed": SEED,
        "chain_roles": {"A": "mhc_alpha", "B": "mhc_beta", "C": "peptide"},
        "groove_alignment": "first 85 CA atoms from chains A and B; deterministic Kabsch fit",
        "candidate_exposed_positions": [2, 3, 5, 7, 8],
        "patch_definition": "solvent-accessible heavy-atom surface points from target EBV model 0 side chains at P2/P3/P5/P7/P8; fixed within panel",
        "patch_samples_per_atom": 48,
        "primary_metric": "full_pmhc_hodgkin_similarity_q25_across_25_model_combinations",
        "secondary_metrics": ["Carbo_similarity", "sign_agreement", "potential_RMSE_q75", "HLA_subtracted_potential"],
        "primary_rank_direction": "higher_Hodgkin_is_better",
        "apbs": {
            "version": "3.4.1",
            "equation": "linearized_Poisson_Boltzmann",
            "solute_dielectrics": list(DIELECTRIC_VALUES),
            "primary_solute_dielectric": 2.0,
            "solvent_dielectric": 78.5,
            "temperature_K": 298.15,
            "monovalent_salt_M": 0.15,
            "solvent_radius_A": 1.4,
            "maximum_fine_grid_spacing_A": 0.5,
            "shared_grid_per_panel": True,
        },
        "pdb2pqr": {"version": "3.7.1", "force_field": "PARSE", "pH": 7.4, "titration": "PROPKA_3.5.1"},
        "gate": "rank_1_to_3_and_register_QC_and_complete_model_calculation_QC",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
        "upstream": upstream,
    }
    write_json(output_dir / "protocol_lock.json", protocol)

    environment = {
        "pdb2pqr_executable": str(PDB2PQR),
        "pdb2pqr_version": _tool_version([str(PDB2PQR), "--version"]),
        "pdb2pqr_python": "3.13",
        "apbs_archive": str(APBS_ARCHIVE),
        "apbs_archive_sha256": sha256_file(APBS_ARCHIVE),
        "apbs_install": str(APBS_INSTALL),
        "apbs_container_image": APBS_IMAGE,
        "apbs_container_image_digest": APBS_IMAGE_DIGEST,
        "apbs_container_platform": "linux/amd64",
        "native_macos_apbs_status": "unavailable_upstream_intel_python_and_FETK_ARPACK_link_incompatibility",
        "container_apbs_version_verified": True,
    }
    write_json(output_dir / "environment_manifest.json", environment)

    audit_rows = read_csv(V3_DIR / "model_surface_audit.csv")
    audit_index: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in audit_rows:
        audit_index.setdefault((row["allele"], row["candidate_id"]), []).append(row)
    model_rows: list[dict[str, Any]] = []
    arm_rows: list[dict[str, Any]] = []
    alignment_rows: list[dict[str, Any]] = []
    for panel_id in TARGET_IDS:
        panel = panels[panel_id]
        target = next(row for row in panel if row["row_role"] == "target")
        arms = _arm_registry(panel_id, panel)
        arm_rows.extend(arms)
        reference_key = (target["allele"], target["ebv_candidate_id"])
        references = sorted(audit_index.get(reference_key, []), key=lambda row: int(row["sample_index"]))
        if len(references) != 5:
            raise ValueError(f"{panel_id} target EBV arm lacks five models")
        reference = parse_mmcif(Path(references[0]["model_path"]))
        aligned_models: list[dict[str, list[dict[str, Any]]]] = []
        panel_model_records: list[dict[str, Any]] = []
        for arm in arms:
            source_models = sorted(
                audit_index.get((arm["allele"], arm["candidate_id"]), []),
                key=lambda row: int(row["sample_index"]),
            )
            if len(source_models) != 5 or any(row["surface_status"] != "parsed_complete" for row in source_models):
                raise ValueError(f"{panel_id} {arm['candidate_id']} lacks a complete five-model ensemble")
            for source in source_models:
                observed_sha = sha256_file(Path(source["model_path"]))
                if observed_sha != source["model_sha256"]:
                    raise ValueError(f"source model checksum changed: {source['model_path']}")
                model = parse_mmcif(Path(source["model_path"]))
                if "".join(str(residue["aa"]) for residue in model["C"]) != arm["sequence"]:
                    raise ValueError(f"model peptide mismatch for {panel_id} {arm['candidate_id']}")
                aligned, fit_rmsd = align_model_to_reference(model, reference)
                model_key = _safe_token(f"{panel_id}__{arm['arm_class']}__{arm['candidate_id']}__s{source['sample_index']}")
                full_pdb = output_dir / "prepared_structures" / panel_id / f"{model_key}__full.pdb"
                hla_pdb = output_dir / "prepared_structures" / panel_id / f"{model_key}__hla_only.pdb"
                write_model_pdb(aligned, full_pdb, include_peptide=True)
                write_model_pdb(aligned, hla_pdb, include_peptide=False)
                record = {
                    **arm,
                    "model_key": model_key,
                    "sample_index": int(source["sample_index"]),
                    "source_model_path": source["model_path"],
                    "source_model_sha256": source["model_sha256"],
                    "alignment_reference_candidate_id": target["ebv_candidate_id"],
                    "alignment_reference_sample_index": 0,
                    "groove_fit_rmsd_A": fit_rmsd,
                    "full_pdb": _relative(full_pdb, output_dir),
                    "full_pdb_sha256": sha256_file(full_pdb),
                    "hla_only_pdb": _relative(hla_pdb, output_dir),
                    "hla_only_pdb_sha256": sha256_file(hla_pdb),
                }
                model_rows.append(record)
                panel_model_records.append(record)
                alignment_rows.append(
                    {
                        "panel_id": panel_id,
                        "model_key": model_key,
                        "groove_fit_rmsd_A": fit_rmsd,
                        "status": "pass" if fit_rmsd <= 1.0 else "fail",
                    }
                )
                aligned_models.append(aligned)
        patch = surface_patch_points(
            reference,
            core_start_1_based=int(target["ebv_declared_core_start_1_based"]),
            samples_per_atom=48,
        )
        patch_rows = [
            {"panel_id": panel_id, "point_index": index, "x_A": point[0], "y_A": point[1], "z_A": point[2]}
            for index, point in enumerate(patch)
        ]
        write_csv(output_dir / "patch_definitions" / f"{panel_id}__patch_points.csv", patch_rows)
        grid = build_shared_grid([_heavy_bounds(model) for model in aligned_models], patch)
        write_json(output_dir / "grid_definitions" / f"{panel_id}__grid.json", asdict(grid))
        if any(row["status"] != "pass" for row in alignment_rows if row["panel_id"] == panel_id):
            raise ValueError(f"{panel_id} contains a failed groove alignment")
        if len(panel_model_records) != 60:
            raise ValueError(f"{panel_id} must prepare exactly 60 arm models")

    write_csv(output_dir / "panel_arm_registry.csv", arm_rows)
    write_csv(output_dir / "model_registry.csv", model_rows)
    write_csv(output_dir / "alignment_qc.csv", alignment_rows)
    status = {
        "status": "prepared",
        "target_count": 2,
        "panel_count": 2,
        "pair_count": len(frozen_pairs),
        "arm_count": len(arm_rows),
        "model_count": len(model_rows),
        "geometry_selection_order": "comparators_frozen_upstream_before_this_geometry_stage",
        "calculation_complete": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "prepare_status.json", status)
    _checksums(output_dir)
    return status


def _run_logged(command: Sequence[str], stdout_path: Path, stderr_path: Path) -> int:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, text=True, check=False)
    return completed.returncode


def _protonate_model(output_dir: Path, model: Mapping[str, str], form: str) -> dict[str, Any]:
    source = output_dir / model[f"{form}_pdb"]
    directory = output_dir / "raw_calculations/protonation" / model["panel_id"]
    stem = f"{model['model_key']}__{form}"
    pqr = directory / f"{stem}.pqr"
    stdout = directory / f"{stem}.stdout.txt"
    stderr = directory / f"{stem}.stderr.txt"
    command = [
        str(PDB2PQR), "--ff=PARSE", "--keep-chain", "--titration-state-method=propka",
        "--with-ph=7.4", str(source), str(pqr),
    ]
    exit_code = _run_logged(command, stdout, stderr)
    log = pqr.with_suffix(".log")
    return {
        "panel_id": model["panel_id"],
        "model_key": model["model_key"],
        "form": form,
        "status": "complete" if exit_code == 0 and pqr.exists() else "failed",
        "exit_code": exit_code,
        "pqr": _relative(pqr, output_dir) if pqr.exists() else "",
        "pqr_sha256": sha256_file(pqr) if pqr.exists() else "",
        "propka_log": _relative(log, output_dir) if log.exists() else "",
        "stdout": _relative(stdout, output_dir),
        "stderr": _relative(stderr, output_dir),
        "command": " ".join(command),
    }


def _load_patch(output_dir: Path, panel_id: str) -> np.ndarray:
    rows = read_csv(output_dir / "patch_definitions" / f"{panel_id}__patch_points.csv")
    rows.sort(key=lambda row: int(row["point_index"]))
    return np.asarray([[float(row["x_A"]), float(row["y_A"]), float(row["z_A"])] for row in rows])


def _pqr_atoms(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, list[np.ndarray]]]:
    coordinates: list[list[float]] = []
    radii: list[float] = []
    ca_by_chain: dict[str, list[np.ndarray]] = {"A": [], "B": [], "C": []}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        fields = line.split()
        if len(fields) < 11:
            raise ValueError(f"malformed PQR atom row in {path}")
        point = [float(fields[6]), float(fields[7]), float(fields[8])]
        coordinates.append(point)
        radii.append(float(fields[10]))
        if fields[2] == "CA" and fields[4] in ca_by_chain:
            ca_by_chain[fields[4]].append(np.asarray(point, dtype=float))
    if not coordinates:
        raise ValueError(f"no PQR atoms found in {path}")
    return np.asarray(coordinates, dtype=float), np.asarray(radii, dtype=float), ca_by_chain


def _archive_invalid_initial_analysis(output_dir: Path) -> Path:
    archive = output_dir / "discarded_initial_surface_patch_analysis"
    if archive.exists():
        raise FileExistsError(f"discard archive already exists: {archive}")
    archive.mkdir(parents=True)
    candidates = [
        "raw_calculations/apbs", "sampled_potentials.csv", "apbs_provenance.csv",
        "calculation_status.json", "model_combination_electrostatic_metrics.csv",
        "panel_electrostatic_ranks.csv", "target_electrostatic_summary.csv",
        "electrostatic_context_gate.json", "specificity_gate.json",
        "development_control_calibration_gate.json", "analysis_status.json",
        "figure_status.json", "figures", "lead_dossiers", "README.md",
    ]
    moved = []
    for relative in candidates:
        source = output_dir / relative
        if not source.exists():
            continue
        destination = archive / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))
        moved.append(relative)
    write_json(
        archive / "DISCARD_REASON.json",
        {
            "status": "invalid_discarded_before_final_interpretation",
            "reason": "the initial target-surface sampling points were not solvent-accessible across all comparator models",
            "score_used_to_define_correction": False,
            "replacement_rule": "geometry-only common solvent-accessible positional field shell",
            "moved_artifacts": moved,
        },
    )
    return archive


def _freeze_common_field_patches(output_dir: Path) -> dict[str, Any]:
    models = read_csv(output_dir / "model_registry.csv")
    targets = {row["target_id"]: row for row in read_csv(output_dir / "frozen_target_registry.csv")}
    protonation_rows = read_csv(output_dir / "protonation_provenance.csv")
    protonation = {(row["model_key"], row["form"]): row for row in protonation_rows}
    qc_rows: list[dict[str, Any]] = []
    patch_status: dict[str, Any] = {}
    for panel_id in TARGET_IDS:
        target = targets[panel_id]
        panel_models = [row for row in models if row["panel_id"] == panel_id]
        atomsets: list[tuple[np.ndarray, np.ndarray]] = []
        pqr_by_model: dict[str, tuple[np.ndarray, np.ndarray, dict[str, list[np.ndarray]]]] = {}
        for model in panel_models:
            pqr_row = protonation.get((model["model_key"], "full"))
            if not pqr_row or pqr_row["status"] != "complete":
                raise ValueError(f"missing full protonation result for {model['model_key']}")
            parsed = _pqr_atoms(output_dir / pqr_row["pqr"])
            pqr_by_model[model["model_key"]] = parsed
            atomsets.append((parsed[0], parsed[1]))
        reference_model = next(
            row for row in panel_models
            if row["candidate_id"] == target["ebv_candidate_id"] and int(row["sample_index"]) == 0
        )
        _, _, ca = pqr_by_model[reference_model["model_key"]]
        if len(ca["A"]) < 85 or len(ca["B"]) < 85:
            raise ValueError(f"{panel_id} reference PQR lacks groove CA coordinates")
        start = int(reference_model["core_start_1_based"]) - 1
        core_ca = np.vstack(ca["C"][start : start + 9])
        groove_ca = np.vstack(ca["A"][:85] + ca["B"][:85])
        patch, metadata = build_common_accessible_field_patch(
            core_ca,
            groove_ca,
            atomsets,
            minimum_height_A=2.0,
            maximum_height_A=25.0,
            height_step_A=0.25,
            probe_radius_A=1.4,
            minimum_clearance_A=0.25,
        )
        patch_rows = [
            {
                "panel_id": panel_id,
                "point_index": index,
                "x_A": point[0],
                "y_A": point[1],
                "z_A": point[2],
                **metadata[index],
            }
            for index, point in enumerate(patch)
        ]
        write_csv(output_dir / "patch_definitions" / f"{panel_id}__patch_points.csv", patch_rows)
        grid = build_shared_grid(
            [(coordinates.min(axis=0), coordinates.max(axis=0)) for coordinates, _ in atomsets],
            patch,
        )
        write_json(output_dir / "grid_definitions" / f"{panel_id}__grid.json", asdict(grid))
        for model in panel_models:
            coordinates, radii, _ = pqr_by_model[model["model_key"]]
            clearance = np.min(
                np.linalg.norm(patch[:, None, :] - coordinates[None, :, :], axis=2)
                - (radii[None, :] + 1.4),
                axis=1,
            )
            qc_rows.append(
                {
                    "panel_id": panel_id,
                    "model_key": model["model_key"],
                    "point_count": len(patch),
                    "solvent_accessible_point_count": int(np.sum(clearance >= 0.25 - 1e-8)),
                    "solvent_accessible_fraction": float(np.mean(clearance >= 0.25 - 1e-8)),
                    "minimum_clearance_A": float(clearance.min()),
                    "status": "pass" if np.all(clearance >= 0.25 - 1e-8) else "fail",
                }
            )
        heights = [float(row["height_A"]) for row in metadata]
        patch_status[panel_id] = {
            "point_count": len(patch),
            "height_min_A": min(heights),
            "height_median_A": float(np.median(heights)),
            "height_max_A": max(heights),
            "all_models_all_points_accessible": all(row["status"] == "pass" for row in qc_rows if row["panel_id"] == panel_id),
        }
    write_csv(output_dir / "patch_accessibility_qc.csv", qc_rows)
    if len(qc_rows) != 120 or any(row["status"] != "pass" for row in qc_rows):
        raise ValueError("common field patch accessibility QC failed")
    protocol_path = output_dir / "protocol_lock.json"
    protocol = json.loads(protocol_path.read_text())
    protocol["patch_definition"] = (
        "25 fixed position-matched field points: five local points for each of P2/P3/P5/P7/P8, "
        "each shifted along the HLA-to-peptide outward normal by the smallest 0.25 A increment "
        "that gives at least 0.25 A clearance beyond the 1.4 A probe surface in every panel model"
    )
    protocol["patch_accessibility_requirement"] = "all_points_accessible_in_all_60_models_per_panel"
    protocol["protocol_amendment"] = {
        "date": "2026-08-29",
        "reason": "initial target-surface points failed common solvent-accessibility QC",
        "uses_electrostatic_scores": False,
        "initial_analysis_status": "discarded_invalid",
    }
    write_json(protocol_path, protocol)
    status = {
        "status": "complete",
        "panel_status": patch_status,
        "model_qc_row_count": len(qc_rows),
        "all_models_all_points_accessible": True,
        "geometry_only_no_electrostatic_scores_used": True,
    }
    write_json(output_dir / "patch_accessibility_status.json", status)
    return status


def refreeze_common_field_patch(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    _archive_invalid_initial_analysis(output_dir)
    status = _freeze_common_field_patches(output_dir)
    _checksums(output_dir)
    return status


def _apbs_model_dielectric(
    output_dir: Path,
    model: Mapping[str, str],
    protonation: Mapping[tuple[str, str], Mapping[str, str]],
    dielectric: float,
    grid: GridSpec,
    patch: np.ndarray,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    directory = output_dir / "raw_calculations/apbs" / model["panel_id"] / model["model_key"] / f"eps{int(dielectric)}"
    directory.mkdir(parents=True, exist_ok=True)
    potentials: dict[str, np.ndarray] = {}
    provenance: list[dict[str, Any]] = []
    for form in ("full", "hla_only"):
        pqr_row = protonation[(model["model_key"], form)]
        if pqr_row["status"] != "complete":
            provenance.append({"model_key": model["model_key"], "form": form, "solute_dielectric": dielectric, "status": "not_run_protonation_failed"})
            continue
        pqr_host = output_dir / pqr_row["pqr"]
        prefix_host = directory / form
        input_host = directory / f"{form}.in"
        params = APBSParameters(solute_dielectric=float(dielectric))
        pqr_container = Path("/work") / pqr_host.relative_to(output_dir)
        prefix_container = Path("/work") / prefix_host.relative_to(output_dir)
        input_host.write_text(build_apbs_input(pqr_container, prefix_container, grid, params), encoding="ascii")
        stdout = directory / f"{form}.stdout.txt"
        stderr = directory / f"{form}.stderr.txt"
        command = [
            "docker", "run", "--rm", "--platform", "linux/amd64",
            "-v", f"{APBS_INSTALL}:/apbs:ro", "-v", f"{output_dir}:/work", "-w", "/work",
            "-e", "LD_LIBRARY_PATH=/apbs/lib", APBS_IMAGE, "/apbs/bin/apbs",
            str(Path("/work") / input_host.relative_to(output_dir)),
        ]
        exit_code = _run_logged(command, stdout, stderr)
        dx = Path(str(prefix_host) + ".dx")
        status = "complete" if exit_code == 0 and dx.exists() else "failed"
        record = {
            "panel_id": model["panel_id"],
            "model_key": model["model_key"],
            "candidate_id": model["candidate_id"],
            "arm_class": model["arm_class"],
            "sample_index": model["sample_index"],
            "form": form,
            "solute_dielectric": dielectric,
            "status": status,
            "exit_code": exit_code,
            "input": _relative(input_host, output_dir),
            "input_sha256": sha256_file(input_host),
            "stdout": _relative(stdout, output_dir),
            "stderr": _relative(stderr, output_dir),
            "dx_sha256": sha256_file(dx) if dx.exists() else "",
            "dx_bytes": dx.stat().st_size if dx.exists() else "",
            "dx_retention": "transient_deleted_after_checksum_and_sampling" if dx.exists() else "not_created",
            "container_image_digest": APBS_IMAGE_DIGEST,
        }
        provenance.append(record)
        if status == "complete":
            potentials[form] = trilinear_sample(parse_open_dx(dx), patch)
            dx.unlink()
    samples: list[dict[str, Any]] = []
    if set(potentials) == {"full", "hla_only"}:
        normalized = potentials["full"] - potentials["hla_only"]
        for index in range(len(patch)):
            samples.append(
                {
                    "panel_id": model["panel_id"],
                    "model_key": model["model_key"],
                    "candidate_id": model["candidate_id"],
                    "arm_class": model["arm_class"],
                    "sample_index": model["sample_index"],
                    "solute_dielectric": dielectric,
                    "point_index": index,
                    "full_potential_kT_per_e": float(potentials["full"][index]),
                    "hla_only_potential_kT_per_e": float(potentials["hla_only"][index]),
                    "hla_normalized_potential_kT_per_e": float(normalized[index]),
                }
            )
    return samples, provenance


def calculate(output_dir: Path = DEFAULT_OUT, workers: int = 4) -> dict[str, Any]:
    if not (output_dir / "prepare_status.json").exists():
        raise FileNotFoundError("prepare stage has not completed")
    if not PDB2PQR.exists() or not (APBS_INSTALL / "bin/apbs").exists():
        raise FileNotFoundError("version-pinned PDB2PQR/APBS tools are unavailable")
    models = read_csv(output_dir / "model_registry.csv")
    protonation_path = output_dir / "protonation_provenance.csv"
    protonation_reused = False
    if protonation_path.exists():
        protonation_rows = read_csv(protonation_path)
        protonation_reused = (
            len(protonation_rows) == 2 * len(models)
            and all(row["status"] == "complete" for row in protonation_rows)
            and all(
                (output_dir / row["pqr"]).exists()
                and sha256_file(output_dir / row["pqr"]) == row["pqr_sha256"]
                for row in protonation_rows
            )
        )
    else:
        protonation_rows = []
    if not protonation_reused:
        protonation_rows = []
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = [pool.submit(_protonate_model, output_dir, model, form) for model in models for form in ("full", "hla_only")]
            for future in as_completed(futures):
                protonation_rows.append(future.result())
        protonation_rows.sort(key=lambda row: (row["panel_id"], row["model_key"], row["form"]))
        write_csv(protonation_path, protonation_rows)
    protonation = {(row["model_key"], row["form"]): row for row in protonation_rows}

    histidine_rows: list[dict[str, Any]] = []
    target_by_panel = {row["target_id"]: row for row in read_csv(output_dir / "frozen_target_registry.csv")}
    for model in models:
        if model["candidate_id"] != target_by_panel[model["panel_id"]]["ebv_candidate_id"]:
            continue
        pqr_row = protonation[(model["model_key"], "full")]
        for row in parse_propka_histidines(output_dir / pqr_row["propka_log"] if pqr_row["propka_log"] else Path("/nonexistent")):
            core_position = row["sequence_position_1_based"] - int(model["core_start_1_based"]) + 1
            histidine_rows.append(
                {
                    "panel_id": model["panel_id"],
                    "allele": model["allele"],
                    "model_key": model["model_key"],
                    "sample_index": model["sample_index"],
                    **row,
                    "core_position": core_position if 1 <= core_position <= 9 else "outside_declared_core",
                    "candidate_exposed_position": core_position in (2, 3, 5, 7, 8),
                    "protonated_at_pH_7_4_by_pka_rule": float(row["predicted_pka"]) > 7.4,
                }
            )
    write_csv(
        output_dir / "target_histidine_propka.csv",
        histidine_rows,
        ("panel_id", "allele", "model_key", "sample_index", "sequence_position_1_based", "core_position", "candidate_exposed_position", "predicted_pka", "protonated_at_pH_7_4_by_pka_rule"),
    )

    patch_status_path = output_dir / "patch_accessibility_status.json"
    if not patch_status_path.exists() or not json.loads(patch_status_path.read_text()).get("all_models_all_points_accessible"):
        _freeze_common_field_patches(output_dir)

    grid_by_panel = {
        panel_id: GridSpec(**json.loads((output_dir / "grid_definitions" / f"{panel_id}__grid.json").read_text()))
        for panel_id in TARGET_IDS
    }
    patch_by_panel = {panel_id: _load_patch(output_dir, panel_id) for panel_id in TARGET_IDS}
    sample_rows: list[dict[str, Any]] = []
    apbs_rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = [
            pool.submit(
                _apbs_model_dielectric, output_dir, model, protonation, dielectric,
                grid_by_panel[model["panel_id"]], patch_by_panel[model["panel_id"]],
            )
            for model in models
            for dielectric in DIELECTRIC_VALUES
        ]
        for future in as_completed(futures):
            samples, provenance = future.result()
            sample_rows.extend(samples)
            apbs_rows.extend(provenance)
    sample_rows.sort(key=lambda row: (row["panel_id"], row["model_key"], float(row["solute_dielectric"]), int(row["point_index"])))
    apbs_rows.sort(key=lambda row: (row.get("panel_id", ""), row["model_key"], float(row["solute_dielectric"]), row["form"]))
    write_csv(output_dir / "sampled_potentials.csv", sample_rows)
    write_csv(output_dir / "apbs_provenance.csv", apbs_rows)
    expected_pqr = 2 * len(models)
    expected_apbs = 2 * len(models) * len(DIELECTRIC_VALUES)
    expected_vectors = len(models) * len(DIELECTRIC_VALUES)
    complete_vectors = len({(row["model_key"], row["solute_dielectric"]) for row in sample_rows})
    status_value = "complete" if (
        len(protonation_rows) == expected_pqr
        and all(row["status"] == "complete" for row in protonation_rows)
        and len(apbs_rows) == expected_apbs
        and all(row["status"] == "complete" for row in apbs_rows)
        and complete_vectors == expected_vectors
    ) else "not_evaluable_incomplete_calculations"
    status = {
        "status": status_value,
        "protonation_run_count": len(protonation_rows),
        "expected_protonation_run_count": expected_pqr,
        "protonation_results_reused": protonation_reused,
        "apbs_run_count": len(apbs_rows),
        "expected_apbs_run_count": expected_apbs,
        "complete_sampled_vector_count": complete_vectors,
        "expected_sampled_vector_count": expected_vectors,
        "transient_dx_files_retained": False,
        "raw_input_log_and_sample_linkage_retained": True,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "calculation_status.json", status)
    _checksums(output_dir)
    return status


def _vector_index(rows: Sequence[Mapping[str, str]]) -> dict[tuple[str, str, int, float, str], np.ndarray]:
    grouped: dict[tuple[str, str, int, float], list[Mapping[str, str]]] = {}
    for row in rows:
        key = (row["panel_id"], row["candidate_id"], int(row["sample_index"]), float(row["solute_dielectric"]))
        grouped.setdefault(key, []).append(row)
    result = {}
    for key, values in grouped.items():
        values.sort(key=lambda row: int(row["point_index"]))
        result[(*key, "full")] = np.asarray([float(row["full_potential_kT_per_e"]) for row in values])
        result[(*key, "hla_normalized")] = np.asarray([float(row["hla_normalized_potential_kT_per_e"]) for row in values])
    return result


def _write_figures(output_dir: Path, target_rows: Sequence[Mapping[str, Any]], pair_summaries: Sequence[Mapping[str, Any]]) -> str:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return "not_evaluable_matplotlib_unavailable"
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for target in target_rows:
        panel_id = target["target_id"]
        values = [
            row for row in pair_summaries
            if row["panel_id"] == panel_id and row["potential_mode"] == "full" and float(row["solute_dielectric"]) == 2.0
        ]
        values.sort(key=lambda row: int(row["electrostatic_rank"]))
        colors = ["#c23b3b" if row["row_role"] == "target" else "#5b778e" for row in values]
        fig, ax = plt.subplots(figsize=(8.0, 4.5))
        ax.bar(range(1, 27), [float(row["hodgkin_similarity_q25"]) for row in values], color=colors)
        ax.set_xlabel("Panel rank")
        ax.set_ylabel("Conservative Hodgkin similarity (q25)")
        ax.set_title(f"{target['allele']} local electrostatic context")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(figures / f"{panel_id}__primary_rank.png", dpi=180)
        plt.close(fig)
    return "complete"


def analyze(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    calculation = json.loads((output_dir / "calculation_status.json").read_text())
    if calculation["status"] != "complete":
        raise RuntimeError("calculation stage is incomplete; analysis cannot manufacture missing values")
    pairs = read_csv(output_dir / "frozen_panel_pairs.csv")
    samples = read_csv(output_dir / "sampled_potentials.csv")
    vectors = _vector_index(samples)
    combination_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for pair in pairs:
        panel_id = pair["target_id"]
        for dielectric in DIELECTRIC_VALUES:
            for mode in ("full", "hla_normalized"):
                ensemble = []
                for left_index in range(5):
                    for right_index in range(5):
                        left = vectors[(panel_id, pair["ebv_candidate_id"], left_index, float(dielectric), mode)]
                        right = vectors[(panel_id, pair["self_candidate_id"], right_index, float(dielectric), mode)]
                        metrics = {
                            "hodgkin_similarity": hodgkin_similarity(left, right),
                            "carbo_similarity": carbo_similarity(left, right),
                            "sign_agreement_fraction": sign_agreement_fraction(left, right),
                            "potential_rmse": potential_rmse(left, right),
                        }
                        ensemble.append(metrics)
                        combination_rows.append(
                            {
                                "panel_id": panel_id,
                                "allele": pair["allele"],
                                "pair_id": pair["pair_id"],
                                "row_role": pair["row_role"],
                                "solute_dielectric": dielectric,
                                "potential_mode": mode,
                                "ebv_sample_index": left_index,
                                "self_sample_index": right_index,
                                **metrics,
                            }
                        )
                summary_rows.append(
                    {
                        "panel_id": panel_id,
                        "allele": pair["allele"],
                        "pair_id": pair["pair_id"],
                        "row_role": pair["row_role"],
                        "ebv_candidate_id": pair["ebv_candidate_id"],
                        "self_candidate_id": pair["self_candidate_id"],
                        "solute_dielectric": dielectric,
                        "potential_mode": mode,
                        **summarize_electrostatic_ensemble(ensemble),
                        "n3_is_not_specificity_negative": pair["row_role"] == "n3",
                        "claim_boundary": CLAIM_BOUNDARY,
                    }
                )
    ranked_rows: list[dict[str, Any]] = []
    for panel_id in TARGET_IDS:
        for dielectric in DIELECTRIC_VALUES:
            for mode in ("full", "hla_normalized"):
                group = [
                    row for row in summary_rows
                    if row["panel_id"] == panel_id and float(row["solute_dielectric"]) == float(dielectric) and row["potential_mode"] == mode
                ]
                ranked_rows.extend(rank_panel(group))
    ranked_rows.sort(key=lambda row: (row["panel_id"], float(row["solute_dielectric"]), row["potential_mode"], int(row["electrostatic_rank"])))
    write_csv(output_dir / "model_combination_electrostatic_metrics.csv", combination_rows)
    write_csv(output_dir / "panel_electrostatic_ranks.csv", ranked_rows)

    targets = read_csv(output_dir / "frozen_target_registry.csv")
    target_summaries: list[dict[str, Any]] = []
    lead_gates = []
    for target in targets:
        panel_id = target["target_id"]
        rows = [row for row in ranked_rows if row["panel_id"] == panel_id and row["row_role"] == "target"]
        rank_lookup = {(row["potential_mode"], float(row["solute_dielectric"])): int(row["electrostatic_rank"]) for row in rows}
        class_by_dielectric = {
            dielectric: "supportive" if rank_lookup[("full", dielectric)] <= 3 else "not_supportive"
            for dielectric in DIELECTRIC_VALUES
        }
        robust = dielectric_robustness(class_by_dielectric)
        register_qc = _truth(target.get("register_robust"))
        model_qc = all(key in rank_lookup for key in (("full", 2.0), ("full", 4.0), ("full", 8.0), ("hla_normalized", 2.0)))
        primary_rank = rank_lookup.get(("full", 2.0))
        gate = build_lead_gate(
            target_id=panel_id,
            rank=primary_rank,
            register_qc=register_qc,
            model_qc=model_qc,
            dielectric_robust=robust,
        )
        lead_gates.append(gate)
        primary_row = next(row for row in rows if row["potential_mode"] == "full" and float(row["solute_dielectric"]) == 2.0)
        primary_panel = [
            row for row in ranked_rows
            if row["panel_id"] == panel_id and row["potential_mode"] == "full" and float(row["solute_dielectric"]) == 2.0
        ]
        best_n3 = max(
            (row for row in primary_panel if row["row_role"] == "n3"),
            key=lambda row: (float(row["hodgkin_similarity_q25"]), str(row["pair_id"])),
        )
        margin = float(primary_row["hodgkin_similarity_q25"]) - float(best_n3["hodgkin_similarity_q25"])
        target_summaries.append(
            {
                "target_id": panel_id,
                "allele": target["allele"],
                "ebv_core": target["ebv_core"],
                "self_core": target["self_core"],
                "primary_full_pmhc_rank": primary_rank,
                "primary_full_pmhc_hodgkin_q25": primary_row["hodgkin_similarity_q25"],
                "best_n3_pair_id": best_n3["pair_id"],
                "best_n3_hodgkin_q25": best_n3["hodgkin_similarity_q25"],
                "target_margin_over_best_n3": margin,
                "hla_normalized_rank_eps2": rank_lookup.get(("hla_normalized", 2.0)),
                "full_rank_eps4": rank_lookup.get(("full", 4.0)),
                "full_rank_eps8": rank_lookup.get(("full", 8.0)),
                "dielectric_rank_class_robust": robust,
                "register_robust": register_qc,
                "final_status": gate["status"],
                "rank_only_context": gate["rank_only_context"],
                "interpretation": "rank is descriptive sensitivity context; final gate abstains when register robustness is unresolved",
                "claim_boundary": CLAIM_BOUNDARY,
            }
        )
    write_csv(output_dir / "target_electrostatic_summary.csv", target_summaries)
    gate_status = (
        "electrostatic_context_supportive"
        if all(row["status"] == "electrostatic_context_supportive" for row in lead_gates)
        else "not_evaluable"
        if any(row["status"] == "not_evaluable" for row in lead_gates)
        else "electrostatic_context_not_supportive"
    )
    ranking_gate = {
        "status": gate_status,
        "lead_results": lead_gates,
        "panel_count": 2,
        "n3_specificity_role": "unknown_recognition_comparators_only",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "electrostatic_context_gate.json", ranking_gate)
    write_json(
        output_dir / "specificity_gate.json",
        {
            "status": "not_evaluable_no_N1_or_N2_assay_negatives",
            "N3_rows_used_for_specificity": 0,
            "specificity_claim_allowed": False,
            "claim_boundary": CLAIM_BOUNDARY,
        },
    )
    prior_gate_path = V2_RESULTS_DIR / "pilot_attribution_gate.json"
    prior_gate = json.loads(prior_gate_path.read_text()) if prior_gate_path.exists() else {}
    write_json(
        output_dir / "development_control_calibration_gate.json",
        {
            "status": "not_evaluable_new_electrostatic_endpoint_not_run_on_controls",
            "control_role": "existing_three_system_controls_are_development_only",
            "prior_v2_pilot_status": prior_gate.get("status", "not_available"),
            "independent_validation_claim_allowed": False,
            "reason": "This bounded two-lead pilot does not constitute an untouched positive-control validation of electrostatics.",
        },
    )
    figure_status = _write_figures(output_dir, targets, ranked_rows)
    try:
        matplotlib_version = importlib.metadata.version("matplotlib")
    except importlib.metadata.PackageNotFoundError:
        matplotlib_version = "not_installed"
    write_json(output_dir / "figure_status.json", {"status": figure_status, "matplotlib_version": matplotlib_version})
    environment_path = output_dir / "environment_manifest.json"
    environment = json.loads(environment_path.read_text())
    environment["matplotlib_version"] = matplotlib_version
    write_json(environment_path, environment)
    for target, result in zip(targets, target_summaries):
        text = f"""# {target['allele']}: {target['ebv_core']} vs {target['self_core']}

- Full-pMHC electrostatic rank at solute dielectric 2: **{result['primary_full_pmhc_rank']} of 26**.
- Margin over the best N3 pair: **{float(result['target_margin_over_best_n3']):.4f} Hodgkin units**.
- HLA-subtracted sensitivity rank: **{result['hla_normalized_rank_eps2']} of 26**.
- Solute-dielectric sensitivity ranks (2/4/8): **{result['primary_full_pmhc_rank']} / {result['full_rank_eps4']} / {result['full_rank_eps8']}**.
- Formal context status: **{result['final_status']}**.
- Register robust under the frozen V3 rule: **{result['register_robust']}**.

The numerical rank describes resemblance within this frozen exact-HLA computational panel. It does not show presentation, TCR recognition, cross-reactivity, molecular mimicry, or an MS mechanism.
"""
        dossier = output_dir / "lead_dossiers" / f"{target['target_id']}.md"
        dossier.parent.mkdir(parents=True, exist_ok=True)
        dossier.write_text(text, encoding="utf-8")
    readme = f"""# pMHC Surface Electrostatics Pilot

This additive package evaluates two frozen BALF5--TALDO1 leads against their previously frozen 25-pair, exact-HLA N3 panels. It uses PDB2PQR 3.7.1 with PROPKA at pH 7.4 and PARSE parameters, followed by APBS 3.4.1 linearized Poisson--Boltzmann calculations under a shared panel grid.

## Result

Overall gate: `{gate_status}`.

The target ranks are reported in `target_electrostatic_summary.csv`. Because the frozen V3 register-robustness flag is false for both leads, electrostatic ranks remain sensitivity evidence and the formal gate abstains. N3 pairs have unknown recognition status and were not treated as specificity negatives.

## Reproducibility

- `protocol_lock.json` freezes the method and claim boundary.
- `environment_manifest.json` records tool versions, archive checksum, and the container image digest.
- `protonation_provenance.csv` and `target_histidine_propka.csv` retain charge-state provenance.
- `apbs_provenance.csv` links every APBS input, log, transient grid checksum, and sampled vector.
- Full OpenDX files were deleted after checksum and deterministic point sampling to avoid a multi-gigabyte package; every input needed to regenerate them is retained.
- No existing V1--V3 package or discovery ranking was modified.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    by_target = {row["target_id"]: row for row in target_summaries}
    dr13 = by_target["HY13_SEQ_02"]
    dr15 = by_target["HY15_SEQ_02"]
    results_summary = f"""# Results Summary

## Main finding

The corrected common-field electrostatic analysis **does not support either BALF5--TALDO1 lead as unusually similar within its frozen exact-HLA N3 panel**.

| HLA | Pair | Full-pMHC rank (eps-in 2) | eps-in 4 / 8 | HLA-subtracted rank | Rank-only context | Formal status |
|---|---|---:|---:|---:|---|---|
| DRB1*13:03 | `{dr13['ebv_core']}` / `{dr13['self_core']}` | {dr13['primary_full_pmhc_rank']}/26 | {dr13['full_rank_eps4']}/26 / {dr13['full_rank_eps8']}/26 | {dr13['hla_normalized_rank_eps2']}/26 | `{dr13['rank_only_context']}` | `{dr13['final_status']}` |
| DRB1*15:01 | `{dr15['ebv_core']}` / `{dr15['self_core']}` | {dr15['primary_full_pmhc_rank']}/26 | {dr15['full_rank_eps4']}/26 / {dr15['full_rank_eps8']}/26 | {dr15['hla_normalized_rank_eps2']}/26 | `{dr15['rank_only_context']}` | `{dr15['final_status']}` |

The unfavorable full-pMHC rank class is stable across solute dielectric values 2, 4, and 8 for both leads. DRB1*15:01 improves to rank {dr15['hla_normalized_rank_eps2']} after subtracting the modeled HLA field, but it still misses the frozen top-three support rule. DRB1*13:03 does not improve after subtraction.

## Interpretation

- These pairs remain **sequence-supported hypotheses**, not sequence-plus-electrostatics-supported leads.
- The electrostatic result should lower their priority relative to candidates that eventually show agreement across sequence, register, and a validated surface endpoint.
- It does not prove nonrecognition. N3 pairs have unknown TCR-recognition status, and modeled fields can miss induced fit, water networks, dynamics, and an actual TCR footprint.
- Both formal lead gates remain `not_evaluable` because their registers are not robust under the frozen V3 rule.
- The new electrostatic endpoint was not run on the three development controls, so it is not positive-control validated and cannot unlock discovery or specificity claims.

## QC correction

The initial target-surface point analysis was discarded after common-solvent-accessibility QC failed. The final analysis uses 25 position-matched P2/P3/P5/P7/P8 field points that pass in all 60 models per panel. The correction was geometry-only and did not use electrostatic scores.
"""
    (output_dir / "RESULTS_SUMMARY.md").write_text(results_summary, encoding="utf-8")
    methods = """# Methods

Two frozen BALF5--TALDO1 targets were evaluated separately for HLA-DRB1*13:03 and HLA-DRB1*15:01. Each target was compared with its previously frozen 5-by-5 exact-HLA N3 panel (25 comparator pairs). N3 denotes unknown TCR recognition and is not a specificity-negative class.

Five AlphaFold models per peptide arm were aligned to a panel reference using the first 85 C-alpha atoms of each HLA chain. PDB2PQR 3.7.1 assigned PARSE charges and radii after PROPKA titration at pH 7.4. APBS 3.4.1 solved the linearized Poisson--Boltzmann equation at 298.15 K, 0.15 M monovalent salt, solvent dielectric 78.5, solvent radius 1.4 A, and solute dielectric 2 (primary), 4, and 8 (sensitivity).

The final comparison shell contains five local samples for each declared TCR-facing position P2/P3/P5/P7/P8. Each point was shifted outward in 0.25 A increments until it cleared the 1.4 A probe surface by at least 0.25 A in every one of the 60 models in its panel. Full-pMHC and HLA-only potentials were sampled at the identical 25 coordinates. The primary metric is the 25th percentile Hodgkin similarity across all 25 EBV-model by self-model combinations; lower-quartile Carbo similarity, sign agreement, and upper-quartile potential RMSE are secondary.

This is a descriptive modeled-pMHC comparison. It does not measure presentation, TCR binding, activation, specificity, cross-reactivity, molecular mimicry, MS causation, probability, or false-discovery rate.
"""
    (output_dir / "METHODS.md").write_text(methods, encoding="utf-8")
    result = {
        "status": gate_status,
        "target_count": 2,
        "panel_pair_count": 52,
        "model_combination_metric_rows": len(combination_rows),
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "analysis_status.json", result)
    _checksums(output_dir)
    return result


def run_all(output_dir: Path = DEFAULT_OUT, workers: int = 4) -> dict[str, Any]:
    prepare(output_dir)
    calculate(output_dir, workers=workers)
    return analyze(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "refreeze-patch", "calculate", "analyze", "all"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.stage == "prepare":
        result = prepare(args.output_dir)
    elif args.stage == "refreeze-patch":
        result = refreeze_common_field_patch(args.output_dir)
    elif args.stage == "calculate":
        result = calculate(args.output_dir, workers=args.workers)
    elif args.stage == "analyze":
        result = analyze(args.output_dir)
    else:
        result = run_all(args.output_dir, workers=args.workers)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
