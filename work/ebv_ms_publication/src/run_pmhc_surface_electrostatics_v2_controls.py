"""Control-first dense pMHC surface electrostatics V2 workflow.

Stages are deliberately separate. ``prepare`` freezes and materializes only
positive-control inputs; ``calculate`` runs version-pinned PDB2PQR/APBS;
``analyze`` applies the locked panel ranks; and ``verify`` audits the package.
No candidate or discovery package is an allowed input.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import math
import os
import random
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from hla2_positive_control_benchmark import parse_mmcif_atoms, residue_sequence
from hla2_positive_control_benchmark_v2 import blosum62_similarity
from pmhc_surface_electrostatics import (
    carbo_similarity,
    hodgkin_similarity,
    parse_open_dx,
    potential_rmse,
    sign_agreement_fraction,
    trilinear_sample,
    write_model_pdb,
)
from pmhc_surface_electrostatics_v2 import (
    align_local_groove,
    build_apbs_surface_input,
    build_control_gate,
    canonical_groove_frame,
    dense_lateral_grid,
    hierarchical_pair_summary,
    label_surface_regions,
    sample_outer_surface,
    standardize_pmhc_chains,
    validate_control_only_paths,
)


ROOT = Path(__file__).resolve().parents[1]
V1 = ROOT / "processed/hla2_positive_control_benchmark_2026-08-25"
V2_PROTOCOL = ROOT / "processed/hla2_positive_control_benchmark_v2_pilot_2026-08-26"
V2_RESULTS = ROOT / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26"
DEFAULT_OUT = ROOT / "processed/pmhc_surface_electrostatics_v2_controls_2026-08-30"
PDB2PQR = Path.home() / ".cache/ebv_ms_tools/pmhc_electrostatics/pdb2pqr-3.7.1-py313-venv/bin/pdb2pqr"
APBS_INSTALL = Path.home() / ".cache/ebv_ms_tools/pmhc_electrostatics/apbs-3.4.1-linux/APBS-3.4.1.Linux"
NATIVE_APBS = Path.home() / ".cache/ebv_ms_tools/apbs-3.4.1-osx-arm64/bin/apbs"
APBS_IMAGE = "ubuntu:22.04"
APBS_IMAGE_DIGEST = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
BOOTSTRAP_SEED = 271828
CLAIM_BOUNDARY = (
    "Development-control pMHC surface resemblance only; not evidence of presentation, "
    "TCR recognition, activation, specificity, cross-reactivity, molecular mimicry, "
    "MS mechanism, probability, or false-discovery rate."
)
INPUT_PATHS = (V1, V2_PROTOCOL, V2_RESULTS)
PRIMARY_VARIANT = "npbe_eps4_grid0.50"
PHYSICAL_VARIANTS = (
    (PRIMARY_VARIANT, 4.0, False, 0.50, False),
    ("npbe_eps2_grid0.50", 2.0, False, 0.50, False),
    ("npbe_eps8_grid0.50", 8.0, False, 0.50, False),
    ("lpbe_eps4_grid0.50", 4.0, True, 0.50, False),
    ("npbe_eps4_grid0.35", 4.0, False, 0.35, False),
    ("npbe_eps4_grid0.50_fixed_charge", 4.0, False, 0.50, True),
)
MAP_DENSITIES = {"coarse": 0.50, "fine": 0.20}
SURFACE_OFFSETS_A = (0.0, 0.5, 1.0)
REFERENCE_BY_FAMILY = {"DR": "1H15", "DQ": "4MAY"}
REQUIRED_PAIRS = (
    "PAIR_HY2E11_BALF5_MBP",
    "PAIR_OB1A12_ENGA_MBP",
    "PAIR_HY1B11_UL15_MBP",
    "PAIR_HY1B11_PMM_MBP",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fields or sorted({key for row in rows for key in row}))
    if not names:
        raise ValueError(f"field names are required for {path}")
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


def _relative(path: Path, root: Path = ROOT) -> str:
    return str(path.resolve().relative_to(root.resolve()))


def _safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def _sample_index(path: Path) -> int:
    match = re.search(r"_model_(\d+)\.cif$", path.name)
    if not match:
        raise ValueError(f"cannot read model index from {path.name}")
    return int(match.group(1))


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _input_snapshot() -> list[dict[str, Any]]:
    rows = []
    for directory in INPUT_PATHS:
        for path in sorted(item for item in directory.rglob("*") if item.is_file()):
            rows.append({
                "input_package": directory.name,
                "relative_path": _relative(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            })
    return rows


def _write_checksums(output_dir: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(item for item in output_dir.rglob("*") if item.is_file()):
        if path.name == "SHA256SUMS.csv" or path.suffix == ".dx":
            continue
        rows.append({"relative_path": _relative(path, output_dir), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    write_csv(output_dir / "SHA256SUMS.csv", rows, ("relative_path", "sha256", "bytes"))
    return rows


def _family(beta_allele: str) -> str:
    if "DRB" in beta_allele:
        return "DR"
    if "DQB" in beta_allele:
        return "DQ"
    raise ValueError(f"unsupported HLA-II beta allele: {beta_allele}")


def _heavy_coordinates(model: Mapping[str, Sequence[Mapping[str, Any]]], chain: str) -> np.ndarray:
    points = [
        atom["xyz"]
        for residue in model[chain]
        for atom in residue["atoms"]
        if str(atom.get("element", "")).upper() != "H"
    ]
    return np.asarray(points, dtype=float)


def _ca_coordinates(residues: Sequence[Mapping[str, Any]]) -> np.ndarray:
    points = []
    for residue in residues:
        atom = next((atom for atom in residue["atoms"] if atom["name"] == "CA"), None)
        if atom is None:
            raise ValueError("missing CA")
        points.append(atom["xyz"])
    return np.asarray(points, dtype=float)


def _protocol_lock(snapshot_sha256: str) -> dict[str, Any]:
    return {
        "protocol_version": "pmhc_surface_electrostatics_v2_controls",
        "frozen_at": "2026-08-30",
        "input_snapshot_sha256": snapshot_sha256,
        "control_systems": ["Hy.2E11", "Ob.1A12", "Hy.1B11"],
        "required_positive_pairs": list(REQUIRED_PAIRS),
        "experimental_pdb_count": 10,
        "positive_experimental_pdb_count": 7,
        "alphafold_source_model_count": 360,
        "alphafold_seeds": [271828, 314159],
        "alignment": "sequence-equivalent HLA N/CA/C/O atoms within 12 A of reference P1-P9 core",
        "canonical_frame": "P1-to-P9 longitudinal; groove-to-peptide outward; beta-ward transverse",
        "lateral_extent": {"longitudinal": "P1-4A through P9+4A", "transverse_A": [-14.0, 14.0]},
        "map_density_spacing_A": MAP_DENSITIES,
        "minimum_composite_points": 500,
        "minimum_peptide_points": 100,
        "minimum_pairwise_coverage": 0.90,
        "surface_crossing": "last APBS smol 0.5 crossing along canonical outward axis",
        "surface_offsets_A": list(SURFACE_OFFSETS_A),
        "region_rule": "nearest heavy atom; peptide, hla_alpha, or hla_beta",
        "peptide_pair_region_rule": "union of peptide-labelled coordinates in the two compared maps",
        "helix_pair_region_rule": "intersection of non-peptide-labelled coordinates in the two maps",
        "electrostatic_metric": "Hodgkin similarity; higher is better",
        "shape_metrics": ["surface_height_RMSE_A", "surface_normal_mean_dot"],
        "physical_variants": [
            {"variant": name, "solute_dielectric": eps, "linear_pb": linear, "maximum_grid_spacing_A": spacing, "fixed_charge": fixed}
            for name, eps, linear, spacing, fixed in PHYSICAL_VARIANTS
        ],
        "primary_variant": PRIMARY_VARIANT,
        "primary_map_density": "fine",
        "primary_surface_offset_A": 0.5,
        "ensemble_rule": "lower of median left-arm marginal medians and median right-arm marginal medians",
        "resampling": {"draws": 1000, "seed": BOOTSTRAP_SEED, "threshold_top3_fraction": 0.80},
        "gate": "all required PDB-evaluable and both-seed AF peptide/composite/shape ranks top3 across sensitivities",
        "candidate_files_permitted": False,
        "candidate_evaluation_runs_in_this_package": False,
        "weighted_composite_created": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }


def _load_pdb_registry() -> list[dict[str, str]]:
    rows = read_csv(V1 / "benchmark/pdb_structural_ligand_registry.csv")
    selected = [row for row in rows if _truth(row.get("selected_for_oracle_pool", True))]
    if len(selected) != 10 or len({row["pdb_id"] for row in selected}) != 10:
        raise ValueError("expected ten unique frozen experimental structures")
    return selected


def _reference_models(pdb_rows: Sequence[Mapping[str, str]]) -> tuple[dict[str, dict[str, list[dict[str, Any]]]], dict[str, int]]:
    by_pdb = {row["pdb_id"].upper(): row for row in pdb_rows}
    models: dict[str, dict[str, list[dict[str, Any]]]] = {}
    starts: dict[str, int] = {}
    for family, pdb_id in REFERENCE_BY_FAMILY.items():
        row = by_pdb[pdb_id]
        parsed = parse_mmcif_atoms(V1 / "sources/pdb" / f"{pdb_id}.cif")
        model, _ = _standardize_registry_model(parsed, row)
        models[family] = model
        starts[family] = int(row["core_start_1_based"])
    return models, starts


def _standardize_registry_model(
    parsed: Mapping[str, Sequence[Mapping[str, Any]]],
    row: Mapping[str, str],
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Handle both separate peptide chains and curated peptide-beta fusions."""
    alpha_chain = row["mhc_alpha_chain"]
    beta_chain = row["mhc_beta_chain"]
    peptide_chain = row["peptide_chain"]
    if peptide_chain != beta_chain:
        return standardize_pmhc_chains(
            parsed,
            alpha_chain=alpha_chain,
            beta_chain=beta_chain,
            peptide_chain=peptide_chain,
        )
    sequence = residue_sequence(parsed[beta_chain])
    peptide = row["peptide_sequence"]
    start = sequence.find(peptide)
    if start != 0:
        raise ValueError(f"linked peptide is not the curated beta-chain prefix in {row['pdb_id']}")
    split = len(peptide)
    standardized = {
        "A": [dict(residue) for residue in parsed[alpha_chain]],
        "B": [dict(residue) for residue in parsed[beta_chain][split:]],
        "C": [dict(residue) for residue in parsed[beta_chain][:split]],
    }
    excluded = sorted(chain for chain in parsed if chain not in {alpha_chain, beta_chain})
    return standardized, {
        "source_alpha_chain": alpha_chain,
        "source_beta_chain": beta_chain,
        "source_peptide_chain": peptide_chain,
        "peptide_beta_fusion_split_after_residue": split,
        "excluded_chain_ids": excluded,
        "tcr_or_other_protein_chains_removed": bool(excluded),
    }


def _write_aligned_model(
    output_dir: Path,
    model_key: str,
    standardized: Mapping[str, Sequence[Mapping[str, Any]]],
    reference: Mapping[str, Sequence[Mapping[str, Any]]],
    reference_core_start: int,
) -> tuple[Path, dict[str, Any]]:
    aligned, qc = align_local_groove(
        standardized,
        reference,
        reference_core_start_1_based=reference_core_start,
        cutoff_A=12.0,
    )
    path = output_dir / "aligned_models" / f"{_safe(model_key)}.pdb"
    write_model_pdb(aligned, path, include_peptide=True)
    return path, qc


def prepare(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    validate_control_only_paths(INPUT_PATHS, ROOT)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty V2 package: {output_dir}")
    snapshot = _input_snapshot()
    write_csv(output_dir / "input_immutability_snapshot.csv", snapshot)
    snapshot_sha = sha256_file(output_dir / "input_immutability_snapshot.csv")
    protocol = _protocol_lock(snapshot_sha)
    write_json(output_dir / "protocol_lock.json", protocol)
    pdb_rows = _load_pdb_registry()
    references, reference_starts = _reference_models(pdb_rows)

    registry_dir = output_dir / "registries"
    registry_dir.mkdir(parents=True, exist_ok=True)
    for source, name in (
        (V1 / "benchmark/pdb_structural_ligand_registry.csv", "experimental_structure_registry.csv"),
        (V2_PROTOCOL / "registry/positive_pair_registry.csv", "positive_pair_registry.csv"),
        (V2_PROTOCOL / "controls/comparison_universe.csv", "alphafold_pair_registry.csv"),
        (V2_PROTOCOL / "benchmark/pdb_oracle_frozen_pairings.csv", "pdb_pair_registry.csv"),
    ):
        shutil.copy2(source, registry_dir / name)

    model_rows: list[dict[str, Any]] = []
    alignment_rows: list[dict[str, Any]] = []
    for row in pdb_rows:
        source = V1 / "sources/pdb" / f"{row['pdb_id']}.cif"
        parsed = parse_mmcif_atoms(source)
        standardized, removal = _standardize_registry_model(parsed, row)
        if row["core_sequence"] not in residue_sequence(standardized["C"]):
            raise ValueError(f"curated core missing from {row['pdb_id']}")
        family = _family(row["mhc_beta_allele"])
        model_key = f"pdb|{row['ligand_id']}"
        aligned_path, qc = _write_aligned_model(output_dir, model_key, standardized, references[family], reference_starts[family])
        model_rows.append({
            "model_key": model_key,
            "layer": "experimental_pdb",
            "ligand_id": row["ligand_id"],
            "pmhc_id": row["ligand_id"],
            "panel_seed": "",
            "sample_index": 0,
            "hla_family": family,
            "mhc_alpha_allele": row["mhc_alpha_allele"],
            "mhc_beta_allele": row["mhc_beta_allele"],
            "peptide_sequence": row["peptide_sequence"],
            "core_sequence": row["core_sequence"],
            "core_start_1_based": int(row["core_start_1_based"]),
            "source_path": _relative(source),
            "source_sha256": sha256_file(source),
            "aligned_pdb": _relative(aligned_path, output_dir),
            "aligned_pdb_sha256": sha256_file(aligned_path),
            "excluded_chain_ids": ";".join(removal["excluded_chain_ids"]),
            "tcr_or_other_chains_removed": removal["tcr_or_other_protein_chains_removed"],
        })
        alignment_rows.append({"model_key": model_key, **qc})

    jobs = read_csv(V2_RESULTS / "inventory/job_inventory.csv")
    complete_jobs = [row for row in jobs if row["download_status"] == "complete_exact"]
    if len(complete_jobs) != 72:
        raise ValueError(f"expected 72 complete AlphaFold jobs, found {len(complete_jobs)}")
    for job in complete_jobs:
        directory = Path(job["canonical_directory"])
        cifs = sorted(directory.glob("*_model_*.cif"), key=_sample_index)
        if len(cifs) != 5:
            raise ValueError(f"{job['job_name']} does not contain five model CIFs")
        family = _family(job["mhc_beta_allele"])
        for cif in cifs:
            sample_index = _sample_index(cif)
            parsed = parse_mmcif_atoms(cif)
            standardized, removal = standardize_pmhc_chains(parsed, alpha_chain="A", beta_chain="B", peptide_chain="C")
            observed = residue_sequence(standardized["C"])
            if observed != job["peptide_sequence"]:
                raise ValueError(f"peptide mismatch in {cif}")
            model_key = f"af|{job['panel_seed']}|{job['ligand_id']}|m{sample_index}"
            aligned_path, qc = _write_aligned_model(output_dir, model_key, standardized, references[family], reference_starts[family])
            model_rows.append({
                "model_key": model_key,
                "layer": "alphafold",
                "ligand_id": job["ligand_id"],
                "pmhc_id": job["ligand_id"],
                "panel_seed": int(job["panel_seed"]),
                "sample_index": sample_index,
                "hla_family": family,
                "mhc_alpha_allele": job["mhc_alpha_allele"],
                "mhc_beta_allele": job["mhc_beta_allele"],
                "peptide_sequence": job["peptide_sequence"],
                "core_sequence": job["core_sequence"],
                "core_start_1_based": int(job["core_start_1_based"]),
                "source_path": str(cif),
                "source_sha256": sha256_file(cif),
                "aligned_pdb": _relative(aligned_path, output_dir),
                "aligned_pdb_sha256": sha256_file(aligned_path),
                "excluded_chain_ids": ";".join(removal["excluded_chain_ids"]),
                "tcr_or_other_chains_removed": removal["tcr_or_other_protein_chains_removed"],
            })
            alignment_rows.append({"model_key": model_key, **qc})
    if len(model_rows) != 370 or sum(row["layer"] == "alphafold" for row in model_rows) != 360:
        raise ValueError("expected 10 PDB plus 360 AlphaFold models")
    if any(float(row["fit_rmsd_A"]) > 3.0 for row in alignment_rows):
        raise ValueError("local groove alignment exceeded the frozen 3 A QC ceiling")
    write_csv(output_dir / "model_registry.csv", model_rows)
    write_csv(output_dir / "alignment_qc.csv", alignment_rows)
    environment = {
        "python": subprocess.run(["python3", "--version"], capture_output=True, text=True).stdout.strip(),
        "numpy": np.__version__,
        "pdb2pqr": str(PDB2PQR),
        "pdb2pqr_sha256": sha256_file(PDB2PQR),
        "apbs": str(APBS_INSTALL / "bin/apbs"),
        "apbs_sha256": sha256_file(APBS_INSTALL / "bin/apbs"),
        "container_image": APBS_IMAGE,
        "container_digest": APBS_IMAGE_DIGEST,
    }
    write_json(output_dir / "environment_lock.json", environment)
    status = {
        "status": "prepared",
        "experimental_structure_count": 10,
        "experimental_positive_structure_count": 7,
        "alphafold_job_count": 72,
        "alphafold_source_model_count": 360,
        "candidate_or_discovery_inputs_read": False,
        "geometry_scores_read_during_freeze": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "prepare_status.json", status)
    _write_checksums(output_dir)
    return status


def _run_logged(command: Sequence[str], stdout_path: Path, stderr_path: Path) -> int:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, text=True, check=False)
    return result.returncode


def _protonate(output_dir: Path, model: Mapping[str, str]) -> dict[str, Any]:
    source = output_dir / model["aligned_pdb"]
    directory = output_dir / "raw_calculations/protonation"
    pqr = directory / f"{_safe(model['model_key'])}.pqr"
    stdout = directory / f"{_safe(model['model_key'])}.stdout.txt"
    stderr = directory / f"{_safe(model['model_key'])}.stderr.txt"
    command = [
        str(PDB2PQR), "--ff=PARSE", "--keep-chain", "--titration-state-method=propka",
        "--with-ph=7.4", str(source), str(pqr),
    ]
    exit_code = _run_logged(command, stdout, stderr)
    return {
        "model_key": model["model_key"],
        "status": "complete" if exit_code == 0 and pqr.exists() else "failed",
        "exit_code": exit_code,
        "pqr": _relative(pqr, output_dir) if pqr.exists() else "",
        "pqr_sha256": sha256_file(pqr) if pqr.exists() else "",
        "propka_log": _relative(pqr.with_suffix(".log"), output_dir) if pqr.with_suffix(".log").exists() else "",
        "stdout": _relative(stdout, output_dir),
        "stderr": _relative(stderr, output_dir),
    }


def _pqr_rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text(encoding="ascii", errors="strict").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        fields = line.split()
        if len(fields) < 11:
            raise ValueError(f"malformed PQR row in {path}")
        rows.append({
            "record": fields[0], "serial": int(fields[1]), "atom": fields[2], "resname": fields[3],
            "chain": fields[4], "resid": int(fields[5]), "x": float(fields[6]), "y": float(fields[7]),
            "z": float(fields[8]), "charge": float(fields[9]), "radius": float(fields[10]),
        })
    if not rows:
        raise ValueError(f"no atoms in {path}")
    return rows


def _pqr_atom_key(row: Mapping[str, Any]) -> tuple[str, int, str]:
    return str(row["chain"]), int(row["resid"]), str(row["atom"])


def _charge_signature(rows: Sequence[Mapping[str, Any]], chains: set[str]) -> str:
    values = [
        f"{row['chain']}:{row['resid']}:{row['atom']}:{float(row['charge']):.4f}:{float(row['radius']):.4f}"
        for row in rows if row["chain"] in chains
    ]
    return hashlib.sha256("|".join(values).encode("ascii")).hexdigest()


def _write_pqr(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for row in rows:
        atom = str(row["atom"])
        lines.append(
            f"{row['record']:<6}{int(row['serial']):5d} {atom:>4s} {row['resname']:>3s} {row['chain']:1s}"
            f"{int(row['resid']):4d}    {float(row['x']):8.3f}{float(row['y']):8.3f}{float(row['z']):8.3f}"
            f" {float(row['charge']):7.4f} {float(row['radius']):6.4f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def _canonicalize_pqr(source: Path, destination: Path, core_start_1_based: int) -> None:
    core, _regions, alpha, beta, _all = _pqr_geometry(source, core_start_1_based)
    frame = canonical_groove_frame(core, alpha, beta)
    rotation = np.column_stack((frame.longitudinal, frame.transverse, frame.outward))
    transformed = []
    for original in _pqr_rows(source):
        row = dict(original)
        point = np.asarray([row["x"], row["y"], row["z"]], dtype=float)
        canonical = (point - frame.origin) @ rotation
        row["x"], row["y"], row["z"] = (float(value) for value in canonical)
        transformed.append(row)
    _write_pqr(destination, transformed)


def _fixed_charge_pqrs(
    output_dir: Path,
    models: Sequence[Mapping[str, str]],
    protonation: Mapping[str, Mapping[str, str]],
) -> list[dict[str, Any]]:
    parsed = {model["model_key"]: _pqr_rows(output_dir / protonation[model["model_key"]]["pqr"]) for model in models}
    hla_groups: dict[tuple[str, str, str, str], list[Mapping[str, str]]] = {}
    ligand_groups: dict[tuple[str, str, str], list[Mapping[str, str]]] = {}
    for model in models:
        layer_seed = model["panel_seed"] if model["layer"] == "alphafold" else "pdb"
        hla_groups.setdefault((model["layer"], layer_seed, model["mhc_alpha_allele"], model["mhc_beta_allele"]), []).append(model)
        ligand_groups.setdefault((model["layer"], layer_seed, model["ligand_id"]), []).append(model)
    hla_donor: dict[tuple[str, str, str, str], Mapping[str, str]] = {
        key: sorted(group, key=lambda row: row["model_key"])[0] for key, group in hla_groups.items()
    }
    peptide_donor: dict[tuple[str, str, str], Mapping[str, str]] = {}
    for key, group in ligand_groups.items():
        signatures: dict[str, list[Mapping[str, str]]] = {}
        for model in group:
            signature = _charge_signature(parsed[model["model_key"]], {"C"})
            signatures.setdefault(signature, []).append(model)
        winning = min(signatures.items(), key=lambda item: (-len(item[1]), item[0]))[1]
        peptide_donor[key] = sorted(winning, key=lambda row: row["model_key"])[0]
    provenance = []
    for model in models:
        layer_seed = model["panel_seed"] if model["layer"] == "alphafold" else "pdb"
        hkey = (model["layer"], layer_seed, model["mhc_alpha_allele"], model["mhc_beta_allele"])
        pkey = (model["layer"], layer_seed, model["ligand_id"])
        hdonor = hla_donor[hkey]
        pdonor = peptide_donor[pkey]
        hcharges = {_pqr_atom_key(row): (row["charge"], row["radius"]) for row in parsed[hdonor["model_key"]] if row["chain"] in {"A", "B"}}
        pcharges = {_pqr_atom_key(row): (row["charge"], row["radius"]) for row in parsed[pdonor["model_key"]] if row["chain"] == "C"}
        fixed = []
        replaced = 0
        for original in parsed[model["model_key"]]:
            row = dict(original)
            donor = hcharges if row["chain"] in {"A", "B"} else pcharges
            if _pqr_atom_key(row) in donor:
                row["charge"], row["radius"] = donor[_pqr_atom_key(row)]
                replaced += 1
            fixed.append(row)
        path = output_dir / "raw_calculations/fixed_charge" / f"{_safe(model['model_key'])}.pqr"
        _write_pqr(path, fixed)
        canonical_pqr = output_dir / "raw_calculations/canonical_pqr" / f"{_safe(model['model_key'])}.pqr"
        canonical_fixed_pqr = output_dir / "raw_calculations/canonical_fixed_charge" / f"{_safe(model['model_key'])}.pqr"
        _canonicalize_pqr(output_dir / protonation[model["model_key"]]["pqr"], canonical_pqr, int(model["core_start_1_based"]))
        _canonicalize_pqr(path, canonical_fixed_pqr, int(model["core_start_1_based"]))
        provenance.append({
            "model_key": model["model_key"], "status": "complete", "fixed_pqr": _relative(path, output_dir),
            "fixed_pqr_sha256": sha256_file(path), "hla_charge_donor_model_key": hdonor["model_key"],
            "peptide_modal_charge_donor_model_key": pdonor["model_key"], "replaced_atom_count": replaced,
            "total_atom_count": len(fixed), "canonical_pqr": _relative(canonical_pqr, output_dir),
            "canonical_pqr_sha256": sha256_file(canonical_pqr),
            "canonical_fixed_pqr": _relative(canonical_fixed_pqr, output_dir),
            "canonical_fixed_pqr_sha256": sha256_file(canonical_fixed_pqr),
        })
    return provenance


def _pqr_geometry(path: Path, core_start_1_based: int) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    rows = _pqr_rows(path)
    atoms_by_region = {
        "hla_alpha": np.asarray([[row["x"], row["y"], row["z"]] for row in rows if row["chain"] == "A" and not row["atom"].startswith("H")]),
        "hla_beta": np.asarray([[row["x"], row["y"], row["z"]] for row in rows if row["chain"] == "B" and not row["atom"].startswith("H")]),
        "peptide": np.asarray([[row["x"], row["y"], row["z"]] for row in rows if row["chain"] == "C" and not row["atom"].startswith("H")]),
    }
    ca = {chain: {} for chain in ("A", "B", "C")}
    for row in rows:
        if row["chain"] in ca and row["atom"] == "CA":
            ca[row["chain"]][row["resid"]] = np.asarray([row["x"], row["y"], row["z"]])
    start = int(core_start_1_based)
    core = np.vstack([ca["C"][index] for index in range(start, start + 9)])
    groove_centroid = core.mean(axis=0)
    alpha = np.vstack([point for point in ca["A"].values() if np.min(np.linalg.norm(core - point, axis=1)) <= 12.0])
    beta = np.vstack([point for point in ca["B"].values() if np.min(np.linalg.norm(core - point, axis=1)) <= 12.0])
    all_coordinates = np.asarray([[row["x"], row["y"], row["z"]] for row in rows])
    return core, atoms_by_region, alpha, beta, all_coordinates


def _grid_definition(
    all_atoms: np.ndarray,
    base_points: np.ndarray,
    outward: np.ndarray,
    maximum_spacing_A: float,
) -> dict[str, tuple[Any, ...]]:
    surface_extent = np.vstack((base_points - 10.0 * outward, base_points + 27.0 * outward))
    fine_min = surface_extent.min(axis=0) - 3.0
    fine_max = surface_extent.max(axis=0) + 3.0
    fine_length = fine_max - fine_min
    dime = tuple(int(math.ceil(float(length) / maximum_spacing_A / 32.0) * 32 + 1) for length in fine_length)
    actual_length = np.asarray([(dimension - 1) * maximum_spacing_A for dimension in dime])
    fine_center = (fine_min + fine_max) / 2.0
    actual_fine_min = fine_center - actual_length / 2.0
    actual_fine_max = fine_center + actual_length / 2.0
    coarse_min = np.minimum(all_atoms.min(axis=0), actual_fine_min) - 20.0
    coarse_max = np.maximum(all_atoms.max(axis=0), actual_fine_max) + 20.0
    return {
        "dime": dime,
        "fine_lengths_A": tuple(float(value) for value in actual_length),
        "fine_center_A": tuple(float(value) for value in fine_center),
        "coarse_lengths_A": tuple(float(value) for value in (coarse_max - coarse_min)),
        "coarse_center_A": tuple(float(value) for value in ((coarse_min + coarse_max) / 2.0)),
    }


def _apbs_command(output_dir: Path, input_path: Path, container_name: str | None) -> list[str]:
    if container_name is None:
        return [str(NATIVE_APBS), str(input_path)]
    container_path = Path("/work") / input_path.relative_to(output_dir)
    return ["docker", "exec", "-e", "LD_LIBRARY_PATH=/apbs/lib", container_name, "/apbs/bin/apbs", str(container_path)]


def _gzip_file(path: Path) -> Path:
    target = path.with_suffix(path.suffix + ".gz")
    with path.open("rb") as source, target.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as destination:
            shutil.copyfileobj(source, destination)
    path.unlink()
    return target


def _calculate_model(
    output_dir: Path,
    model: Mapping[str, str],
    protonation: Mapping[str, str],
    fixed_charge: Mapping[str, str],
    container_name: str | None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    record_path = output_dir / "model_calculation_records" / f"{_safe(model['model_key'])}.json"
    if record_path.exists():
        cached = json.loads(record_path.read_text(encoding="utf-8"))
        expected_vectors = [
            output_dir / "sampled_surface_vectors" / f"{_safe(model['model_key'])}__{variant}.npz"
            for variant, _dielectric, _linear, _spacing, _fixed in PHYSICAL_VARIANTS
        ]
        if cached.get("model_status", {}).get("status") == "complete" and all(path.exists() for path in expected_vectors):
            return cached["model_status"], cached["provenance"]
    primary_pqr = output_dir / fixed_charge["canonical_pqr"]
    core, atoms_by_region, alpha, beta, all_atoms = _pqr_geometry(primary_pqr, int(model["core_start_1_based"]))
    frame = canonical_groove_frame(core, alpha, beta)
    lateral = {name: dense_lateral_grid(core, frame, spacing_A=spacing) for name, spacing in MAP_DENSITIES.items()}
    directory = output_dir / "raw_calculations/apbs" / _safe(model["model_key"])
    grids: dict[str, tuple[Any, Any]] = {}
    surface_cache: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    provenance: list[dict[str, Any]] = []
    for variant, dielectric, linear, spacing, fixed in PHYSICAL_VARIANTS:
        pqr = output_dir / (fixed_charge["canonical_fixed_pqr"] if fixed else fixed_charge["canonical_pqr"])
        grid = _grid_definition(all_atoms, lateral["fine"][0], frame.outward, spacing)
        variant_dir = directory / variant
        variant_dir.mkdir(parents=True, exist_ok=True)
        potential_prefix = variant_dir / "potential"
        accessibility_prefix = variant_dir / "accessibility"
        input_path = variant_dir / "apbs.in"
        write_accessibility = variant in {PRIMARY_VARIANT, "npbe_eps4_grid0.35"}
        tool_pqr = Path("/work") / pqr.relative_to(output_dir) if container_name else pqr
        tool_potential = Path("/work") / potential_prefix.relative_to(output_dir) if container_name else potential_prefix
        tool_accessibility = Path("/work") / accessibility_prefix.relative_to(output_dir) if container_name else accessibility_prefix
        input_path.write_text(build_apbs_surface_input(
            tool_pqr,
            tool_potential,
            tool_accessibility,
            dime=grid["dime"], lengths_A=grid["fine_lengths_A"], center_A=grid["fine_center_A"],
            coarse_lengths_A=grid["coarse_lengths_A"], coarse_center_A=grid["coarse_center_A"],
            solute_dielectric=dielectric, linear=linear, write_accessibility=write_accessibility,
        ), encoding="ascii")
        stdout = variant_dir / "stdout.txt"
        stderr = variant_dir / "stderr.txt"
        exit_code = _run_logged(_apbs_command(output_dir, input_path, container_name), stdout, stderr)
        potential_dx = Path(str(potential_prefix) + ".dx")
        accessibility_dx = Path(str(accessibility_prefix) + ".dx")
        status = "complete" if exit_code == 0 and potential_dx.exists() and (not write_accessibility or accessibility_dx.exists()) else "failed"
        record = {
            "model_key": model["model_key"], "variant": variant, "status": status, "exit_code": exit_code,
            "input": _relative(input_path, output_dir), "input_sha256": sha256_file(input_path),
            "stdout": _relative(stdout, output_dir), "stderr": _relative(stderr, output_dir),
            "potential_dx_sha256": sha256_file(potential_dx) if potential_dx.exists() else "",
            "accessibility_dx_sha256": sha256_file(accessibility_dx) if accessibility_dx.exists() else "",
            "execution_backend": "docker_linux_amd64" if container_name else "native_osx_arm64_conda_forge",
            "apbs_executable_sha256": sha256_file(APBS_INSTALL / "bin/apbs") if container_name else sha256_file(NATIVE_APBS),
            "container_digest": APBS_IMAGE_DIGEST if container_name else "",
        }
        provenance.append(record)
        if status != "complete":
            continue
        potential_grid = parse_open_dx(potential_dx)
        if write_accessibility:
            accessibility_grid = parse_open_dx(accessibility_dx)
            grids["fine_grid" if spacing < 0.5 else "standard_grid"] = (accessibility_grid, grid)
            compressed = _gzip_file(accessibility_dx)
            record["retained_accessibility_map"] = _relative(compressed, output_dir)
            record["retained_accessibility_map_sha256"] = sha256_file(compressed)
        else:
            source_key = "fine_grid" if spacing < 0.5 else "standard_grid"
            if source_key not in grids:
                raise ValueError("accessibility-producing APBS variant must run before reuse")
            accessibility_grid = grids[source_key][0]
        source_key = "fine_grid" if spacing < 0.5 else "standard_grid"
        arrays: dict[str, np.ndarray] = {}
        for density, (base_points, _) in lateral.items():
            cache_key = (source_key, density)
            if cache_key not in surface_cache:
                sampled = sample_outer_surface(
                    accessibility_grid, potential_grid, base_points, frame.outward,
                    search_min_A=-8.0, search_max_A=24.0, search_step_A=0.25, offset_A=0.0,
                )
                coverage = np.asarray([bool(row["covered"]) for row in sampled], dtype=bool)
                heights = np.full(len(sampled), np.nan)
                normals = np.full((len(sampled), 3), np.nan)
                surface_points = np.full((len(sampled), 3), np.nan)
                for index, row in enumerate(sampled):
                    if not row["covered"]:
                        continue
                    heights[index] = row["surface_height_A"]
                    normals[index] = row["normal"]
                    surface_points[index] = row["surface_point"]
                labels = np.asarray(["missing"] * len(sampled), dtype="U12")
                if np.any(coverage):
                    labels[coverage] = label_surface_regions(surface_points[coverage], atoms_by_region)
                surface_cache[cache_key] = {
                    "coverage": coverage, "height": heights, "normal": normals,
                    "surface_point": surface_points, "label": labels,
                }
            geometry = surface_cache[cache_key]
            for offset in SURFACE_OFFSETS_A:
                coverage = geometry["coverage"]
                potentials = np.full(len(coverage), np.nan)
                if np.any(coverage):
                    sample_points = geometry["surface_point"][coverage] + float(offset) * geometry["normal"][coverage]
                    potentials[coverage] = trilinear_sample(potential_grid, sample_points)
                token = f"{variant}__{density}__offset{offset:.1f}"
                arrays[f"{token}__coverage"] = coverage
                arrays[f"{token}__height"] = geometry["height"]
                arrays[f"{token}__normal"] = geometry["normal"]
                arrays[f"{token}__potential"] = potentials
                arrays[f"{token}__label"] = geometry["label"]
        potential_dx.unlink()
        npz_path = output_dir / "sampled_surface_vectors" / f"{_safe(model['model_key'])}__{variant}.npz"
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(npz_path, **arrays)
        record["sampled_vector_file"] = _relative(npz_path, output_dir)
        record["sampled_vector_sha256"] = sha256_file(npz_path)
        record["potential_dx_retention"] = "deleted_after_checksum_and_sampling"
    complete = sum(row["status"] == "complete" for row in provenance)
    model_status = {
        "model_key": model["model_key"], "status": "complete" if complete == len(PHYSICAL_VARIANTS) else "failed",
        "complete_variant_count": complete, "expected_variant_count": len(PHYSICAL_VARIANTS),
    }
    write_json(record_path, {"model_status": model_status, "provenance": provenance})
    return model_status, provenance


def _start_apbs_container(output_dir: Path) -> str:
    name = f"pmhc_v2_apbs_{os.getpid()}"
    command = [
        "docker", "run", "-d", "--rm", "--platform", "linux/amd64", "--name", name,
        "-v", f"{APBS_INSTALL}:/apbs:ro", "-v", f"{output_dir}:/work", "-w", "/work",
        "-e", "LD_LIBRARY_PATH=/apbs/lib", APBS_IMAGE, "sleep", "infinity",
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"could not start APBS container: {result.stderr.strip()}")
    return name


def calculate(output_dir: Path = DEFAULT_OUT, workers: int = 4) -> dict[str, Any]:
    if not (output_dir / "prepare_status.json").exists():
        raise FileNotFoundError("prepare stage has not completed")
    backend = os.environ.get("PMHC_V2_APBS_BACKEND", "native" if NATIVE_APBS.exists() else "docker")
    if backend not in {"native", "docker"}:
        raise ValueError("PMHC_V2_APBS_BACKEND must be native or docker")
    apbs_executable = NATIVE_APBS if backend == "native" else APBS_INSTALL / "bin/apbs"
    if not PDB2PQR.exists() or not apbs_executable.exists():
        raise FileNotFoundError("version-pinned PDB2PQR/APBS tools are unavailable")
    write_json(output_dir / "execution_environment_addendum.json", {
        "actual_apbs_version": "3.4.1",
        "execution_backend": "native_osx_arm64_conda_forge" if backend == "native" else "docker_linux_amd64",
        "apbs_executable": str(apbs_executable),
        "apbs_executable_sha256": sha256_file(apbs_executable),
        "container_digest": APBS_IMAGE_DIGEST if backend == "docker" else None,
        "physical_protocol_changed": False,
        "reason": "Docker Desktop VM storage attachment failure" if backend == "native" else "locked container backend",
    })
    models = read_csv(output_dir / "model_registry.csv")
    protonation_path = output_dir / "protonation_provenance.csv"
    if protonation_path.exists():
        protonation_rows = read_csv(protonation_path)
    else:
        protonation_rows = []
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = [pool.submit(_protonate, output_dir, model) for model in models]
            for future in as_completed(futures):
                protonation_rows.append(future.result())
        protonation_rows.sort(key=lambda row: row["model_key"])
        write_csv(protonation_path, protonation_rows)
    if len(protonation_rows) != len(models) or any(row["status"] != "complete" for row in protonation_rows):
        status = {"status": "not_evaluable_protonation_incomplete", "complete": sum(row["status"] == "complete" for row in protonation_rows), "expected": len(models)}
        write_json(output_dir / "calculation_status.json", status)
        return status
    protonation = {row["model_key"]: row for row in protonation_rows}
    fixed_path = output_dir / "fixed_charge_provenance.csv"
    fixed_rows = read_csv(fixed_path) if fixed_path.exists() else _fixed_charge_pqrs(output_dir, models, protonation)
    if not fixed_path.exists():
        write_csv(fixed_path, fixed_rows)
    fixed = {row["model_key"]: row for row in fixed_rows}
    model_status_rows: list[dict[str, Any]] = []
    apbs_rows: list[dict[str, Any]] = []
    container = _start_apbs_container(output_dir) if backend == "docker" else None
    try:
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = [pool.submit(_calculate_model, output_dir, model, protonation[model["model_key"]], fixed[model["model_key"]], container) for model in models]
            for future in as_completed(futures):
                model_status, provenance = future.result()
                model_status_rows.append(model_status)
                apbs_rows.extend(provenance)
    finally:
        if container is not None:
            subprocess.run(["docker", "stop", container], capture_output=True, text=True, check=False)
    model_status_rows.sort(key=lambda row: row["model_key"])
    apbs_rows.sort(key=lambda row: (row["model_key"], row["variant"]))
    write_csv(output_dir / "model_calculation_status.csv", model_status_rows)
    write_csv(output_dir / "apbs_provenance.csv", apbs_rows)
    complete = sum(row["status"] == "complete" for row in model_status_rows)
    status = {
        "status": "complete" if complete == len(models) else "not_evaluable_incomplete_calculations",
        "complete_model_count": complete, "expected_model_count": len(models),
        "pdb2pqr_run_count": len(protonation_rows), "apbs_run_count": len(apbs_rows),
        "expected_apbs_run_count": len(models) * len(PHYSICAL_VARIANTS),
        "execution_backend": "native_osx_arm64_conda_forge" if backend == "native" else "docker_linux_amd64",
        "transient_potential_grids_retained": False, "accessibility_maps_retained_compressed": True,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "calculation_status.json", status)
    _write_checksums(output_dir)
    return status


def _surface_arrays(
    output_dir: Path,
    model_key: str,
    variant: str,
    density: str,
    offset: float,
) -> dict[str, np.ndarray]:
    path = output_dir / "sampled_surface_vectors" / f"{_safe(model_key)}__{variant}.npz"
    if not path.exists():
        raise FileNotFoundError(path)
    token = f"{variant}__{density}__offset{offset:.1f}"
    with np.load(path) as data:
        return {name: data[f"{token}__{name}"].copy() for name in ("coverage", "height", "normal", "potential", "label")}


def _normalize_lateral_arrays(
    arrays: Mapping[str, np.ndarray],
    *,
    spacing_A: float,
    transverse_half_width_A: float = 14.0,
    longitudinal_margin_A: float = 4.0,
    normalized_p1_p9_span_A: float = 28.0,
) -> dict[str, np.ndarray]:
    """Resample a model-specific P1-P9 lattice onto one canonical index grid."""
    transverse_count = int(round((2.0 * transverse_half_width_A) / spacing_A)) + 1
    point_count = len(arrays["coverage"])
    if point_count % transverse_count:
        raise ValueError("surface lattice cannot be reshaped with the locked transverse extent")
    longitudinal_count = point_count // transverse_count
    source_x = -longitudinal_margin_A + np.arange(longitudinal_count, dtype=float) * spacing_A
    inferred_span = float(source_x[-1] - longitudinal_margin_A)
    if inferred_span <= 0:
        raise ValueError("invalid inferred P1-P9 span")
    canonical_x = np.where(
        source_x < 0.0,
        source_x,
        np.where(
            source_x <= inferred_span,
            source_x * normalized_p1_p9_span_A / inferred_span,
            normalized_p1_p9_span_A + source_x - inferred_span,
        ),
    )
    target_x = np.arange(
        -longitudinal_margin_A,
        normalized_p1_p9_span_A + longitudinal_margin_A + spacing_A / 2.0,
        spacing_A,
    )
    upper = np.clip(np.searchsorted(canonical_x, target_x, side="left"), 0, longitudinal_count - 1)
    lower = np.maximum(upper - 1, 0)
    choose_upper = np.abs(canonical_x[upper] - target_x) < np.abs(target_x - canonical_x[lower])
    nearest = np.where(choose_upper, upper, lower)
    denominator = canonical_x[upper] - canonical_x[lower]
    weight = np.divide(
        target_x - canonical_x[lower], denominator,
        out=np.zeros_like(target_x), where=denominator != 0,
    )

    result: dict[str, np.ndarray] = {}
    for name in ("coverage", "label"):
        source = np.asarray(arrays[name]).reshape(longitudinal_count, transverse_count)
        result[name] = source[nearest].reshape(-1)
    for name in ("height", "potential"):
        source = np.asarray(arrays[name], dtype=float).reshape(longitudinal_count, transverse_count)
        left = source[lower]
        right = source[upper]
        interpolated = left * (1.0 - weight[:, None]) + right * weight[:, None]
        interpolated = np.where(np.isfinite(interpolated), interpolated, np.where(np.isfinite(left), left, right))
        result[name] = interpolated.reshape(-1)
    source_normal = np.asarray(arrays["normal"], dtype=float).reshape(longitudinal_count, transverse_count, 3)
    normal = source_normal[lower] * (1.0 - weight[:, None, None]) + source_normal[upper] * weight[:, None, None]
    norms = np.linalg.norm(normal, axis=2, keepdims=True)
    normal = np.divide(normal, norms, out=np.zeros_like(normal), where=norms > 0)
    result["normal"] = normal.reshape(-1, 3)
    finite = (
        np.isfinite(result["height"])
        & np.isfinite(result["potential"])
        & np.all(np.isfinite(result["normal"]), axis=1)
    )
    result["coverage"] = np.asarray(result["coverage"], dtype=bool) & finite
    return result


def _score_surface_pair(left: Mapping[str, np.ndarray], right: Mapping[str, np.ndarray]) -> dict[str, Any]:
    if len(left["coverage"]) != len(right["coverage"]):
        return {"status": "not_evaluable_lattice_size_mismatch"}
    intersection = left["coverage"] & right["coverage"]
    denominator = min(int(np.sum(left["coverage"])), int(np.sum(right["coverage"])))
    pairwise_coverage = float(np.sum(intersection) / denominator) if denominator else 0.0
    peptide = intersection & ((left["label"] == "peptide") | (right["label"] == "peptide"))
    helix = intersection & (left["label"] != "peptide") & (right["label"] != "peptide")
    composite_count = int(np.sum(intersection))
    peptide_count = int(np.sum(peptide))
    if composite_count < 500 or peptide_count < 100 or pairwise_coverage < 0.90:
        return {
            "status": "not_evaluable_surface_qc", "composite_point_count": composite_count,
            "peptide_point_count": peptide_count, "helix_point_count": int(np.sum(helix)),
            "pairwise_map_coverage": pairwise_coverage,
        }
    left_potential = left["potential"]
    right_potential = right["potential"]
    height_delta = left["height"][intersection] - right["height"][intersection]
    normals = np.sum(left["normal"][intersection] * right["normal"][intersection], axis=1)
    result = {
        "status": "complete", "composite_point_count": composite_count, "peptide_point_count": peptide_count,
        "helix_point_count": int(np.sum(helix)), "pairwise_map_coverage": pairwise_coverage,
        "peptide_hodgkin": hodgkin_similarity(left_potential[peptide], right_potential[peptide]),
        "composite_hodgkin": hodgkin_similarity(left_potential[intersection], right_potential[intersection]),
        "helix_hodgkin": hodgkin_similarity(left_potential[helix], right_potential[helix]),
        "surface_height_rmse_A": float(np.sqrt(np.mean(height_delta * height_delta))),
        "surface_normal_mean_dot": float(np.mean(normals)),
        "peptide_carbo": carbo_similarity(left_potential[peptide], right_potential[peptide]),
        "peptide_sign_agreement": sign_agreement_fraction(left_potential[peptide], right_potential[peptide]),
        "peptide_potential_rmse": potential_rmse(left_potential[peptide], right_potential[peptide]),
    }
    return result


def _ensemble_metric(matrix: np.ndarray, *, lower_is_better: bool = False) -> dict[str, Any]:
    transformed = -np.asarray(matrix, dtype=float) if lower_is_better else np.asarray(matrix, dtype=float)
    summary = hierarchical_pair_summary(transformed)
    if lower_is_better:
        summary["conservative_score"] = -float(summary["conservative_score"])
        summary["left_marginal_medians"] = [-float(value) for value in summary["left_marginal_medians"]]
        summary["right_marginal_medians"] = [-float(value) for value in summary["right_marginal_medians"]]
        summary["left_marginal_summary"] = -float(summary["left_marginal_summary"])
        summary["right_marginal_summary"] = -float(summary["right_marginal_summary"])
    return summary


def _rank_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result = [dict(row) for row in rows]
    endpoints = {
        "peptide": lambda row: (-float(row["peptide_hodgkin"]), str(row["pair_id"])),
        "composite": lambda row: (-float(row["composite_hodgkin"]), str(row["pair_id"])),
        "helix": lambda row: (-float(row["helix_hodgkin"]), str(row["pair_id"])),
        "height": lambda row: (float(row["surface_height_rmse_A"]), str(row["pair_id"])),
        "normal": lambda row: (-float(row["surface_normal_mean_dot"]), str(row["pair_id"])),
        "shape": lambda row: (float(row["surface_height_rmse_A"]), -float(row["surface_normal_mean_dot"]), str(row["pair_id"])),
    }
    for endpoint, key in endpoints.items():
        complete = sorted((row for row in result if row["status"] == "complete"), key=key)
        ranks = {row["pair_id"]: index for index, row in enumerate(complete, start=1)}
        for row in result:
            row[f"{endpoint}_rank"] = ranks.get(row["pair_id"], "")
    return result


def _model_lookup(models: Sequence[Mapping[str, str]]) -> dict[tuple[str, str, str], list[str]]:
    lookup: dict[tuple[str, str, str], list[tuple[int, str]]] = {}
    for row in models:
        seed = row["panel_seed"] if row["layer"] == "alphafold" else "pdb"
        key = (row["layer"], seed, row["pmhc_id"])
        lookup.setdefault(key, []).append((int(row["sample_index"]), row["model_key"]))
    return {key: [model_key for _, model_key in sorted(values)] for key, values in lookup.items()}


def _score_all_panels(output_dir: Path) -> tuple[list[dict[str, Any]], dict[tuple[Any, ...], dict[str, np.ndarray]]]:
    models = read_csv(output_dir / "model_registry.csv")
    lookup = _model_lookup(models)
    af_pairs = read_csv(output_dir / "registries/alphafold_pair_registry.csv")
    pdb_pairs = read_csv(output_dir / "registries/pdb_pair_registry.csv")
    summary_rows: list[dict[str, Any]] = []
    matrices: dict[tuple[Any, ...], dict[str, np.ndarray]] = {}
    for variant, _eps, _linear, _spacing, _fixed in PHYSICAL_VARIANTS:
        for density in MAP_DENSITIES:
            for offset in SURFACE_OFFSETS_A:
                cache: dict[str, dict[str, np.ndarray]] = {}
                def vector(model_key: str) -> dict[str, np.ndarray]:
                    if model_key not in cache:
                        cache[model_key] = _normalize_lateral_arrays(
                            _surface_arrays(output_dir, model_key, variant, density, offset),
                            spacing_A=MAP_DENSITIES[density],
                        )
                    return cache[model_key]
                for pair in af_pairs:
                    seed = pair["panel_seed"]
                    left_models = lookup.get(("alphafold", seed, pair["left_pmhc_id"]), [])
                    right_models = lookup.get(("alphafold", seed, pair["right_pmhc_id"]), [])
                    if len(left_models) != 5 or len(right_models) != 5:
                        summary_rows.append({**pair, "layer": "alphafold", "variant": variant, "density": density, "offset_A": offset, "status": "not_evaluable_missing_five_model_arm"})
                        continue
                    cell_rows = [_score_surface_pair(vector(left), vector(right)) for left in left_models for right in right_models]
                    if any(row["status"] != "complete" for row in cell_rows):
                        incomplete = [row for row in cell_rows if row["status"] != "complete"]
                        worst = min(incomplete, key=lambda row: (row.get("peptide_point_count", 0), row.get("pairwise_map_coverage", 0)))
                        summary_rows.append({**pair, "layer": "alphafold", "variant": variant, "density": density, "offset_A": offset, **worst})
                        continue
                    metric_matrices = {
                        metric: np.asarray([float(row[metric]) for row in cell_rows]).reshape(5, 5)
                        for metric in ("peptide_hodgkin", "composite_hodgkin", "helix_hodgkin", "surface_height_rmse_A", "surface_normal_mean_dot")
                    }
                    key = ("alphafold", pair["positive_pair_id"], seed, pair["pair_id"], variant, density, float(offset))
                    matrices[key] = metric_matrices
                    row = {**pair, "layer": "alphafold", "variant": variant, "density": density, "offset_A": offset, "status": "complete"}
                    for metric, matrix in metric_matrices.items():
                        aggregate = _ensemble_metric(matrix, lower_is_better=metric == "surface_height_rmse_A")
                        row[metric] = aggregate["conservative_score"]
                        row[f"{metric}_left_marginals"] = ";".join(f"{value:.8g}" for value in aggregate["left_marginal_medians"])
                        row[f"{metric}_right_marginals"] = ";".join(f"{value:.8g}" for value in aggregate["right_marginal_medians"])
                    row["composite_point_count"] = min(int(cell["composite_point_count"]) for cell in cell_rows)
                    row["peptide_point_count"] = min(int(cell["peptide_point_count"]) for cell in cell_rows)
                    row["pairwise_map_coverage"] = min(float(cell["pairwise_map_coverage"]) for cell in cell_rows)
                    summary_rows.append(row)
                for pair in pdb_pairs:
                    left_models = lookup.get(("experimental_pdb", "pdb", pair["left_ligand_id"]), [])
                    right_models = lookup.get(("experimental_pdb", "pdb", pair["right_ligand_id"]), [])
                    if len(left_models) != 1 or len(right_models) != 1:
                        summary_rows.append({**pair, "layer": "experimental_pdb", "variant": variant, "density": density, "offset_A": offset, "status": "not_evaluable_missing_structure"})
                        continue
                    score = _score_surface_pair(vector(left_models[0]), vector(right_models[0]))
                    summary_rows.append({**pair, "layer": "experimental_pdb", "variant": variant, "density": density, "offset_A": offset, **score})
    ranked: list[dict[str, Any]] = []
    groups: dict[tuple[str, str, str, str, float], list[dict[str, Any]]] = {}
    for row in summary_rows:
        panel_seed = row.get("panel_seed", "pdb") if row["layer"] == "alphafold" else "pdb"
        groups.setdefault((row["layer"], row["positive_pair_id"], str(panel_seed), row["variant"], row["density"], float(row["offset_A"])), []).append(row)
    for group in groups.values():
        ranked.extend(_rank_rows(group))
    return ranked, matrices


def _write_matrix_archive(
    output_dir: Path,
    matrices: Mapping[tuple[Any, ...], Mapping[str, np.ndarray]],
) -> None:
    metrics = (
        "peptide_hodgkin",
        "composite_hodgkin",
        "helix_hodgkin",
        "surface_height_rmse_A",
        "surface_normal_mean_dot",
    )
    keys = sorted(matrices, key=lambda key: tuple(str(value) for value in key))
    tensor = np.empty((len(keys), len(metrics), 5, 5), dtype=np.float64)
    manifest = []
    for index, key in enumerate(keys):
        layer, positive_pair_id, panel_seed, pair_id, variant, density, offset = key
        for metric_index, metric in enumerate(metrics):
            tensor[index, metric_index] = np.asarray(matrices[key][metric], dtype=np.float64)
        manifest.append({
            "matrix_index": index,
            "layer": layer,
            "positive_pair_id": positive_pair_id,
            "panel_seed": panel_seed,
            "pair_id": pair_id,
            "variant": variant,
            "density": density,
            "offset_A": offset,
            "metric_axis_order": ";".join(metrics),
            "left_model_axis": "sample_index_0_to_4",
            "right_model_axis": "sample_index_0_to_4",
        })
    with (output_dir / "region_score_matrices.npy").open("wb") as handle:
        np.save(handle, tensor, allow_pickle=False)
    write_csv(output_dir / "region_score_matrix_manifest.csv", manifest)


def _resampling_stability(
    ranked_rows: Sequence[Mapping[str, Any]],
    matrices: Mapping[tuple[Any, ...], Mapping[str, np.ndarray]],
) -> list[dict[str, Any]]:
    primary = [
        row for row in ranked_rows
        if row["layer"] == "alphafold" and row["variant"] == PRIMARY_VARIANT
        and row["density"] == "fine" and float(row["offset_A"]) == 0.5 and row["status"] == "complete"
    ]
    groups: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for row in primary:
        groups.setdefault((row["positive_pair_id"], row["panel_seed"]), []).append(row)
    results = []
    for (positive_pair_id, seed), panel in sorted(groups.items()):
        if len(panel) != 26:
            continue
        expected_keys = [
            ("alphafold", positive_pair_id, seed, row["pair_id"], PRIMARY_VARIANT, "fine", 0.5)
            for row in panel
        ]
        if any(key not in matrices for key in expected_keys):
            continue
        identity = f"{positive_pair_id}|{seed}|{BOOTSTRAP_SEED}"
        rng = np.random.default_rng(int(hashlib.sha256(identity.encode()).hexdigest()[:16], 16))
        counts = {endpoint: 0 for endpoint in ("peptide", "composite", "shape")}
        for _ in range(1000):
            left_indices = rng.integers(0, 5, size=5)
            right_indices = rng.integers(0, 5, size=5)
            scored = []
            for row in panel:
                key = ("alphafold", positive_pair_id, seed, row["pair_id"], PRIMARY_VARIANT, "fine", 0.5)
                values = matrices[key]
                entry = {"pair_id": row["pair_id"], "comparison_role": row["comparison_role"]}
                for metric in ("peptide_hodgkin", "composite_hodgkin", "surface_height_rmse_A", "surface_normal_mean_dot"):
                    matrix = values[metric][np.ix_(left_indices, right_indices)]
                    entry[metric] = _ensemble_metric(matrix, lower_is_better=metric == "surface_height_rmse_A")["conservative_score"]
                scored.append(entry)
            ranked = _rank_rows([{**entry, "status": "complete", "helix_hodgkin": 0.0} for entry in scored])
            positive = next(row for row in ranked if row["comparison_role"] == "positive")
            counts["peptide"] += int(int(positive["peptide_rank"]) <= 3)
            counts["composite"] += int(int(positive["composite_rank"]) <= 3)
            counts["shape"] += int(int(positive["shape_rank"]) <= 3)
        for endpoint, count in counts.items():
            results.append({
                "positive_pair_id": positive_pair_id, "panel_seed": seed, "endpoint": endpoint,
                "draw_count": 1000, "top3_count": count, "top3_fraction": count / 1000.0,
                "seed": BOOTSTRAP_SEED,
            })
    return results


def _gate_requirements(
    ranked_rows: Sequence[Mapping[str, Any]],
    resampling: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    resample = {(row["positive_pair_id"], row["panel_seed"], row["endpoint"]): row for row in resampling}
    requirements = []
    primary_rows = [
        row for row in ranked_rows if row["variant"] == PRIMARY_VARIANT and row["density"] == "fine"
        and float(row["offset_A"]) == 0.5 and row.get("comparison_role", row.get("pair_role")) == "positive"
    ]
    for row in primary_rows:
        pair_id = row["positive_pair_id"]
        if row["layer"] == "experimental_pdb" and pair_id not in {"PAIR_HY2E11_BALF5_MBP", "PAIR_OB1A12_ENGA_MBP"}:
            continue
        layer = f"af_{row['panel_seed']}" if row["layer"] == "alphafold" else "pdb"
        all_panel_rows = [
            candidate for candidate in ranked_rows
            if candidate["layer"] == row["layer"] and candidate["positive_pair_id"] == pair_id
            and (row["layer"] != "alphafold" or candidate["panel_seed"] == row["panel_seed"])
        ]
        primary_panel = [
            candidate for candidate in all_panel_rows
            if candidate["variant"] == PRIMARY_VARIANT and candidate["density"] == "fine"
            and float(candidate["offset_A"]) == 0.5
        ]
        condition_groups: dict[tuple[str, str, float], list[Mapping[str, Any]]] = {}
        for candidate in all_panel_rows:
            condition_groups.setdefault(
                (candidate["variant"], candidate["density"], float(candidate["offset_A"])), []
            ).append(candidate)
        panel_complete = bool(primary_panel) and all(candidate["status"] == "complete" for candidate in primary_panel)
        for endpoint in ("peptide", "composite", "shape"):
            complete = panel_complete and row["status"] == "complete" and row.get(f"{endpoint}_rank", "") != ""
            sensitivity = bool(condition_groups) and all(
                all(candidate["status"] == "complete" for candidate in condition_rows)
                and any(
                    candidate.get("comparison_role", candidate.get("pair_role")) == "positive"
                    and candidate.get(f"{endpoint}_rank", "") != ""
                    and int(candidate[f"{endpoint}_rank"]) <= 3
                    for candidate in condition_rows
                )
                for condition_rows in condition_groups.values()
            )
            stability = 1.0 if layer == "pdb" else float(resample.get((pair_id, row["panel_seed"], endpoint), {}).get("top3_fraction", 0.0))
            requirements.append({
                "layer": layer, "pair_id": pair_id, "endpoint": endpoint,
                "status": "complete" if complete else "missing", "rank": int(row[f"{endpoint}_rank"]) if complete else None,
                "sensitivity_top3": sensitivity, "resampling_top3_fraction": stability,
            })
    return requirements


def _baseline_rows(output_dir: Path, ranked_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    source_path = V2_RESULTS / "benchmark/method_rank_long.csv"
    existing = read_csv(source_path)
    rows = []
    for row in existing:
        rows.append({
            "layer": "alphafold", "positive_pair_id": row["positive_pair_id"], "panel_seed": row["panel_seed"],
            "method": row["method"], "positive_rank": row["positive_rank"], "source": _relative(source_path),
        })
    for row in ranked_rows:
        if row["variant"] != PRIMARY_VARIANT or row["density"] != "fine" or float(row["offset_A"]) != 0.5:
            continue
        if row.get("comparison_role", row.get("pair_role")) != "positive":
            continue
        for endpoint in ("peptide", "composite", "helix", "height", "normal", "shape"):
            rows.append({
                "layer": row["layer"], "positive_pair_id": row["positive_pair_id"],
                "panel_seed": row.get("panel_seed", "pdb") if row["layer"] == "alphafold" else "pdb",
                "method": f"surface_v2_{endpoint}", "positive_rank": row.get(f"{endpoint}_rank", ""),
                "source": "this_package_primary_surface_analysis",
            })
    return rows


def _write_figures(output_dir: Path, requirements: Sequence[Mapping[str, Any]]) -> str:
    complete = [row for row in requirements if row.get("rank") not in (None, "")]
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return _write_svg_rank_figure(output_dir, complete)
    labels = [f"{row['layer']}\n{row['pair_id'].replace('PAIR_', '')}\n{row['endpoint']}" for row in complete]
    values = [int(row["rank"]) for row in complete]
    figure, axis = plt.subplots(figsize=(max(10, len(values) * 0.45), 5))
    colors = ["#2f7d4a" if value <= 3 else "#b43c32" for value in values]
    axis.bar(range(len(values)), values, color=colors)
    axis.axhline(3, color="#222222", linewidth=1, linestyle="--")
    axis.set_ylabel("Positive-control rank")
    axis.set_xticks(range(len(values)), labels, rotation=90, fontsize=7)
    axis.set_title("Control-first surface electrostatics V2: primary ranks")
    figure.tight_layout()
    path = output_dir / "figures/primary_control_ranks.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return "complete"


def _write_svg_rank_figure(output_dir: Path, rows: Sequence[Mapping[str, Any]]) -> str:
    """Write a deterministic dependency-free fallback for the primary rank figure."""
    width = max(1200, 95 * len(rows) + 180)
    height = 760
    chart_left, chart_top, chart_bottom = 90, 70, 430
    chart_height = chart_bottom - chart_top
    maximum_rank = max(26, max((int(row["rank"]) for row in rows), default=1))
    bar_width = 52
    step = (width - chart_left - 40) / max(1, len(rows))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="90" y="34" font-family="Arial, sans-serif" font-size="22" font-weight="bold">Control-first surface electrostatics V2: primary ranks</text>',
        f'<line x1="{chart_left}" y1="{chart_bottom}" x2="{width - 30}" y2="{chart_bottom}" stroke="#222"/>',
        f'<line x1="{chart_left}" y1="{chart_top}" x2="{chart_left}" y2="{chart_bottom}" stroke="#222"/>',
    ]
    for tick in (0, 3, 10, 20, maximum_rank):
        y = chart_bottom - (tick / maximum_rank) * chart_height
        parts.extend((
            f'<line x1="{chart_left - 5}" y1="{y:.2f}" x2="{width - 30}" y2="{y:.2f}" stroke="#dddddd"/>',
            f'<text x="{chart_left - 12}" y="{y + 5:.2f}" text-anchor="end" font-family="Arial, sans-serif" font-size="13">{tick}</text>',
        ))
    top3_y = chart_bottom - (3 / maximum_rank) * chart_height
    parts.append(f'<line x1="{chart_left}" y1="{top3_y:.2f}" x2="{width - 30}" y2="{top3_y:.2f}" stroke="#111" stroke-width="2" stroke-dasharray="7 5"/>')
    for index, row in enumerate(rows):
        rank = int(row["rank"])
        x = chart_left + step * index + (step - bar_width) / 2
        bar_height = (rank / maximum_rank) * chart_height
        y = chart_bottom - bar_height
        color = "#2f7d4a" if rank <= 3 else "#b43c32"
        label = f"{row['layer']} | {str(row['pair_id']).replace('PAIR_', '')} | {row['endpoint']}"
        parts.extend((
            f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_width}" height="{bar_height:.2f}" fill="{color}"/>',
            f'<text x="{x + bar_width / 2:.2f}" y="{max(chart_top + 14, y - 7):.2f}" text-anchor="middle" font-family="Arial, sans-serif" font-size="13">{rank}</text>',
            f'<text transform="translate({x + bar_width / 2:.2f},{chart_bottom + 15}) rotate(60)" text-anchor="start" font-family="Arial, sans-serif" font-size="11">{html.escape(label)}</text>',
        ))
    parts.extend((
        '<rect x="90" y="700" width="16" height="16" fill="#2f7d4a"/><text x="114" y="713" font-family="Arial, sans-serif" font-size="13">raw rank within top 3</text>',
        '<rect x="310" y="700" width="16" height="16" fill="#b43c32"/><text x="334" y="713" font-family="Arial, sans-serif" font-size="13">raw rank above 3</text>',
        '<text x="90" y="742" font-family="Arial, sans-serif" font-size="12">Gate results also require complete panels, sensitivity robustness, and AlphaFold resampling stability.</text>',
        '</svg>',
    ))
    path = output_dir / "figures/primary_control_ranks.svg"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return "complete_svg_fallback"


def _requirement_result(value: Mapping[str, Any] | None) -> str:
    if value is None or value.get("status") != "complete":
        return "not evaluable"
    rank = int(value["rank"])
    sensitivity = bool(value.get("sensitivity_top3", False))
    stable = not str(value.get("layer", "")).startswith("af_") or float(
        value.get("resampling_top3_fraction", 0.0)
    ) >= 0.80
    outcome = "pass" if rank <= 3 and sensitivity and stable else "fail"
    return f"rank {rank} ({outcome})"


def _readme_rank_table(requirements: Sequence[Mapping[str, Any]]) -> str:
    indexed = {
        (str(row["layer"]), str(row["pair_id"]), str(row["endpoint"])): row
        for row in requirements
    }
    rows = []
    ordered = (
        ("pdb", "PAIR_HY2E11_BALF5_MBP"),
        ("pdb", "PAIR_OB1A12_ENGA_MBP"),
        ("af_271828", "PAIR_HY1B11_UL15_MBP"),
        ("af_314159", "PAIR_HY1B11_UL15_MBP"),
        ("af_271828", "PAIR_HY1B11_PMM_MBP"),
        ("af_314159", "PAIR_HY1B11_PMM_MBP"),
        ("af_271828", "PAIR_HY2E11_BALF5_MBP"),
        ("af_314159", "PAIR_HY2E11_BALF5_MBP"),
        ("af_271828", "PAIR_OB1A12_ENGA_MBP"),
        ("af_314159", "PAIR_OB1A12_ENGA_MBP"),
    )
    for layer, pair_id in ordered:
        cells = [
            _requirement_result(indexed.get((layer, pair_id, endpoint)))
            for endpoint in ("peptide", "composite", "shape")
        ]
        rows.append(f"| {layer} | {pair_id.replace('PAIR_', '')} | {' | '.join(cells)} |")
    return "\n".join((
        "| Layer | Positive control | Peptide electrostatics | Composite electrostatics | Surface shape |",
        "|---|---|---:|---:|---:|",
        *rows,
    ))


def _finder_duplicate_paths(output_dir: Path) -> list[str]:
    duplicate = re.compile(r" \d+(?:\.[^.]+)?$")
    return [
        _relative(path, output_dir)
        for path in sorted(output_dir.rglob("*"))
        if path.is_file() and duplicate.search(path.name)
    ]


def _gate_invariants(gate: Mapping[str, Any], requirements: Sequence[Mapping[str, Any]]) -> bool:
    failed = int(gate.get("failed_result_count", -1))
    missing = int(gate.get("missing_result_count", -1))
    required = int(gate.get("required_result_count", -1))
    status = str(gate.get("status", ""))
    counted_failures = sum(
        row.get("status") == "complete" and (
            row.get("rank") in (None, "") or int(row["rank"]) > 3
            or not _truth(row.get("sensitivity_top3", False))
            or (
                str(row.get("layer", "")).startswith("af_")
                and float(row.get("resampling_top3_fraction") or 0.0) < 0.80
            )
        )
        for row in requirements
    )
    counted_missing = sum(row.get("status") != "complete" for row in requirements)
    expected_status = "fail" if counted_failures else "not_evaluable" if counted_missing or not requirements else "supportive"
    return all((
        required == len(requirements),
        failed == counted_failures,
        missing == counted_missing,
        status == expected_status,
        bool(gate.get("candidate_evaluation_allowed", False)) == (status == "supportive"),
        bool(gate.get("electrostatics_retired_from_candidate_ranking", False)) == (status == "fail"),
        gate.get("weights_frozen") is False,
        gate.get("discovery_unlock_allowed") is False,
        gate.get("specificity_claim_allowed") is False,
    ))


def analyze(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    calculation = json.loads((output_dir / "calculation_status.json").read_text())
    if calculation.get("status") != "complete":
        gate = build_control_gate([{"layer": "all", "pair_id": "all", "endpoint": "all", "status": "missing", "rank": None, "sensitivity_top3": False, "resampling_top3_fraction": None}])
        write_json(output_dir / "control_gate.json", gate)
        return gate
    ranked_rows, matrices = _score_all_panels(output_dir)
    _write_matrix_archive(output_dir, matrices)
    write_csv(output_dir / "hierarchical_ensemble_summary.csv", ranked_rows)
    write_csv(output_dir / "sensitivity_ranks.csv", [
        row for row in ranked_rows
        if row.get("comparison_role", row.get("pair_role")) == "positive"
    ])
    coverage_rows = [{
        key: row.get(key, "") for key in (
            "layer", "positive_pair_id", "panel_seed", "pair_id", "variant", "density", "offset_A", "status",
            "composite_point_count", "peptide_point_count", "helix_point_count", "pairwise_map_coverage",
        )
    } for row in ranked_rows]
    write_csv(output_dir / "surface_map_coverage.csv", coverage_rows)
    resampling = _resampling_stability(ranked_rows, matrices)
    write_csv(
        output_dir / "resampling_results.csv",
        resampling,
        ("positive_pair_id", "panel_seed", "endpoint", "draw_count", "top3_count", "top3_fraction", "seed"),
    )
    requirements = _gate_requirements(ranked_rows, resampling)
    write_csv(output_dir / "gate_requirement_results.csv", requirements)
    gate = build_control_gate(requirements)
    gate["claim_boundary"] = CLAIM_BOUNDARY
    gate["hy1b11_pdb_oracle_status"] = "not_evaluable_availability"
    gate["candidate_files_read"] = False
    write_json(output_dir / "control_gate.json", gate)
    write_json(output_dir / "specificity_gate.json", {
        "status": "not_evaluable_no_verified_n1_n2_negatives", "n3_counted_as_specificity_negative": False,
        "specificity_claim_allowed": False, "claim_boundary": CLAIM_BOUNDARY,
    })
    write_csv(output_dir / "baseline_comparisons.csv", _baseline_rows(output_dir, ranked_rows))
    write_csv(output_dir / "oracle_availability.csv", [
        {"positive_pair_id": "PAIR_HY2E11_BALF5_MBP", "unique_decoy_count": 7, "status": "evaluable_required"},
        {"positive_pair_id": "PAIR_OB1A12_ENGA_MBP", "unique_decoy_count": 5, "status": "evaluable_required"},
        {"positive_pair_id": "PAIR_HY1B11_UL15_MBP", "unique_decoy_count": 2, "status": "not_evaluable_availability"},
        {"positive_pair_id": "PAIR_HY1B11_PMM_MBP", "unique_decoy_count": 2, "status": "not_evaluable_availability"},
    ])
    figure_status = _write_figures(output_dir, requirements)
    write_json(output_dir / "figure_status.json", {"status": figure_status})
    rank_table = _readme_rank_table(requirements)
    readme = f"""# Control-First Surface Electrostatics V2

This additive package replaces the earlier sparse field descriptor with a dense, model-specific near-surface pMHC map. It evaluates only the frozen Hy.2E11, Ob.1A12, and Hy.1B11 development controls. No discovery candidate was scored.

## Result

The locked control gate is **{gate['status']}**. Candidate evaluation is **{'allowed for a separate future package' if gate['candidate_evaluation_allowed'] else 'not allowed'}**. A failed gate permanently retires electrostatics from candidate ranking; a supportive gate would support only a supplementary descriptor.

## Locked Control Results

{rank_table}

`not evaluable` means at least one member of the complete comparison panel failed the predeclared 90% pairwise surface-map coverage rule. A positive row can therefore have a raw rank in the detailed tables while its gate result remains not evaluable. Completed rank, sensitivity, or resampling failures are retained as failures and make the overall gate fail.

The full pMHC nonlinear-PB primary calculation used protein dielectric 4, solvent dielectric 78.5, 0.15 M monovalent salt, pH 7.4 PARSE/PROPKA charges, and a 0.5 A maximum APBS grid spacing. Every required physical and sampling sensitivity is reported separately. There is no weighted composite.

Hy.1B11 remains unavailable as an experimental-PDB oracle because only two unique exact-HLA decoys were frozen. Its two AlphaFold panels remain mandatory. N3 and structural decoys are not specificity negatives.

{CLAIM_BOUNDARY}
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    status = {
        "status": "complete", "control_gate_status": gate["status"],
        "candidate_evaluation_allowed": gate["candidate_evaluation_allowed"],
        "weights_frozen": False, "discovery_unlock_allowed": False,
        "candidate_or_discovery_files_read": False, "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(output_dir / "analysis_status.json", status)
    _write_checksums(output_dir)
    return status


def verify(output_dir: Path = DEFAULT_OUT) -> dict[str, Any]:
    validate_control_only_paths(INPUT_PATHS, ROOT)
    models = read_csv(output_dir / "model_registry.csv")
    snapshot = read_csv(output_dir / "input_immutability_snapshot.csv")
    snapshot_pass = all((ROOT / row["relative_path"]).exists() and sha256_file(ROOT / row["relative_path"]) == row["sha256"] for row in snapshot)
    checksums = read_csv(output_dir / "SHA256SUMS.csv")
    checksum_pass = all((output_dir / row["relative_path"]).exists() and sha256_file(output_dir / row["relative_path"]) == row["sha256"] for row in checksums)
    calculation = json.loads((output_dir / "calculation_status.json").read_text())
    gate = json.loads((output_dir / "control_gate.json").read_text())
    requirements = read_csv(output_dir / "gate_requirement_results.csv")
    resampling = read_csv(output_dir / "resampling_results.csv")
    manifest = read_csv(output_dir / "region_score_matrix_manifest.csv")
    matrices = np.load(output_dir / "region_score_matrices.npy", allow_pickle=False)
    backend = json.loads((output_dir / "execution_environment_addendum.json").read_text())
    concordance = json.loads((output_dir / "backend_concordance.json").read_text())
    determinism_path = output_dir / "rebuild_determinism.json"
    determinism = json.loads(determinism_path.read_text()) if determinism_path.exists() else {}
    assertions = {
        "input_snapshot_unchanged": snapshot_pass,
        "package_checksums_valid": checksum_pass,
        "model_count_370": len(models) == 370,
        "alphafold_model_count_360": sum(row["layer"] == "alphafold" for row in models) == 360,
        "experimental_structure_count_10": sum(row["layer"] == "experimental_pdb" for row in models) == 10,
        "candidate_paths_absent": not any(any(token in str(path).lower() for token in ("candidate_evidence", "high_yield_candidate", "discovery")) for path in output_dir.rglob("*")),
        "control_gate_present": (output_dir / "control_gate.json").exists(),
        "specificity_gate_present": (output_dir / "specificity_gate.json").exists(),
        "calculation_complete_370_models": calculation.get("status") == "complete" and int(calculation.get("complete_model_count", -1)) == 370,
        "calculation_complete_2220_apbs_runs": int(calculation.get("apbs_run_count", -1)) == 2220 and int(calculation.get("expected_apbs_run_count", -2)) == 2220,
        "matrix_archive_matches_manifest": matrices.ndim == 4 and matrices.shape[0] == len(manifest) and [int(row["matrix_index"]) for row in manifest] == list(range(len(manifest))),
        "thirty_locked_gate_requirements": len(requirements) == 30,
        "gate_invariants_valid": _gate_invariants(gate, requirements),
        "resampling_locked": all(int(row["draw_count"]) == 1000 and int(row["seed"]) == BOOTSTRAP_SEED for row in resampling),
        "backend_protocol_unchanged": backend.get("actual_apbs_version") == "3.4.1" and backend.get("physical_protocol_changed") is False,
        "backend_concordance_pass": concordance.get("status") == "pass_bitwise_identical_sampled_vectors" and concordance.get("all_array_keys_identical") is True and float(concordance.get("maximum_absolute_float_difference", -1.0)) == 0.0,
        "deterministic_rebuild_pass": determinism.get("status") == "pass" and bool(determinism.get("all_files_identical")),
        "no_finder_suffix_duplicates": not _finder_duplicate_paths(output_dir),
    }
    status = "pass" if all(assertions.values()) else "fail"
    result = {"status": status, "assertions": assertions, "claim_boundary": CLAIM_BOUNDARY}
    write_json(output_dir / "verification_status.json", result)
    _write_checksums(output_dir)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "calculate", "analyze", "verify", "all"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.stage == "prepare":
        result = prepare(args.output_dir)
    elif args.stage == "calculate":
        result = calculate(args.output_dir, workers=args.workers)
    elif args.stage == "analyze":
        result = analyze(args.output_dir)
    elif args.stage == "verify":
        result = verify(args.output_dir)
    else:
        prepare(args.output_dir)
        calculate(args.output_dir, workers=args.workers)
        analyze(args.output_dir)
        result = verify(args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
