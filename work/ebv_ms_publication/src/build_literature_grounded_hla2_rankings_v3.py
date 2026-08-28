"""Build additive, literature-grounded HLA-II discovery rankings (V3).

V3 deliberately keeps sequence and modeled structure as separate evidence
families. TCR-facing BLOSUM62 is the primary key; physicochemical mismatch,
identity, and a local peptide-surface fingerprint only break exact upstream
ties. Binding predictions are reported but never enter the mimicry rank.

The output is descriptive pMHC prioritization. It is not evidence of antigen
presentation, TCR binding, activation, cross-reactivity, molecular mimicry,
multiple-sclerosis mechanism, probability, or false-discovery rate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import multiprocessing
from functools import lru_cache
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree

from analyze_af3_pmhc_downloads import ca_coordinates, kabsch, parse_mmcif
from build_same_register_hla_rankings_v2 import (
    ALLELES,
    ALLELE_SLUGS,
    read_csv,
    sequence_metrics,
    sha256_file,
    write_csv,
    write_json,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V2_SEQUENCE = ROOT / "processed/tcell_library_v2_same_register_hla_rankings_2026-08-27/all_hla_ranked_pairs.csv"
DEFAULT_COMBINED_SEQUENCE = ROOT / "processed/drb1501_combined_same_register_rankings_2026-08-27/combined_ranked_pairs.csv"
DEFAULT_V2_SAMPLES = ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25/qc/model_sample_qc.csv"
DEFAULT_V2_REGISTERS = ROOT / "processed/tcell_library_v2_2026-08-22/allele_register_predictions_320.csv"
DEFAULT_LEGACY_SAMPLES = ROOT / "processed/complete_model_pipeline_audit_2026-08-15/canonical_af3_sample_metrics.csv"
DEFAULT_CONTROL_FEATURES = ROOT / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26/benchmark/af3_pair_feature_matrix.csv"
DEFAULT_OUT = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
DEFAULT_CACHE = ROOT / "processed/.cache/literature_grounded_hla2_v3_surface"

TCR_FACING_INDICES = (1, 2, 4, 6, 7)
ANCHOR_INDICES = (0, 3, 5, 8)
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
SURFACE_FEATURES = (
    "sidechain_exposure_mismatch_fraction",
    "exposure_weighted_centroid_rmsd_A",
    "exposure_weighted_orientation_rmsd_A",
    "exposed_distance_matrix_rmsd_A",
    "exposure_weighted_chemistry_mismatch",
    "exposure_weighted_backbone_rmsd_A",
)
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
CLAIM_BOUNDARY = (
    "Descriptive HLA-specific pMHC sequence prioritization with modeled local-surface "
    "annotations only; not evidence of presentation, TCR binding, activation, "
    "cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate."
)


@dataclass(frozen=True)
class SurfaceModel:
    sequence: str
    groove_ca: np.ndarray
    peptide_ca: np.ndarray
    sidechain_centroids: np.ndarray
    sidechain_orientations: np.ndarray
    sidechain_sasa: np.ndarray


_GROOVE_ALIGNMENT_CACHE: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
_SURFACE_SLICE_CACHE: dict[tuple[int, int], SurfaceModel] = {}


def _atom_coordinate(residue: Mapping[str, Any], atom_name: str) -> np.ndarray:
    for atom in residue["atoms"]:
        if str(atom["name"]) == atom_name:
            return np.asarray(atom["xyz"], dtype=float)
    raise ValueError(f"missing {atom_name} atom")


def sidechain_centroid(residue: Mapping[str, Any]) -> np.ndarray:
    """Return the heavy-atom side-chain centroid, with glycine falling back to CA."""
    coordinates = [
        np.asarray(atom["xyz"], dtype=float)
        for atom in residue["atoms"]
        if str(atom["name"]) not in BACKBONE_ATOMS and str(atom["element"]) != "H"
    ]
    if not coordinates:
        return _atom_coordinate(residue, "CA")
    return np.mean(np.asarray(coordinates), axis=0)


def fibonacci_sphere(count: int = 64) -> np.ndarray:
    if count < 8:
        raise ValueError("surface sphere requires at least eight points")
    indices = np.arange(count, dtype=float)
    z = 1.0 - 2.0 * (indices + 0.5) / count
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    angle = math.pi * (3.0 - math.sqrt(5.0)) * indices
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def _sidechain_sasa(
    peptide: Sequence[Mapping[str, Any]],
    all_chains: Sequence[Sequence[Mapping[str, Any]]],
    *,
    probe_A: float = 1.4,
    sphere_points: int = 32,
) -> np.ndarray:
    """Approximate Shrake-Rupley side-chain SASA against the complete pMHC."""
    atoms: list[tuple[np.ndarray, float, int, str]] = []
    peptide_atom_ids: dict[int, list[int]] = defaultdict(list)
    peptide_ids = {id(residue): index for index, residue in enumerate(peptide)}
    for chain in all_chains:
        for residue in chain:
            for atom in residue["atoms"]:
                element = str(atom["element"]).upper()
                if element == "H":
                    continue
                atom_index = len(atoms)
                atoms.append(
                    (
                        np.asarray(atom["xyz"], dtype=float),
                        VDW_RADII.get(element, 1.70) + probe_A,
                        id(residue),
                        str(atom["name"]),
                    )
                )
                if id(residue) in peptide_ids:
                    peptide_atom_ids[peptide_ids[id(residue)]].append(atom_index)
    coordinates = np.asarray([item[0] for item in atoms], dtype=float)
    radii = np.asarray([item[1] for item in atoms], dtype=float)
    tree = cKDTree(coordinates)
    sphere = fibonacci_sphere(sphere_points)
    output = []
    for residue_index, residue in enumerate(peptide):
        eligible = [
            atom_index
            for atom_index in peptide_atom_ids[residue_index]
            if atoms[atom_index][3] not in BACKBONE_ATOMS
        ]
        if not eligible:  # Glycine: use CA exposure as a reproducible fallback.
            eligible = [
                atom_index
                for atom_index in peptide_atom_ids[residue_index]
                if atoms[atom_index][3] == "CA"
            ]
        sasa = 0.0
        for atom_index in eligible:
            center, radius = coordinates[atom_index], radii[atom_index]
            samples = center + sphere * radius
            accessible = 0
            neighborhoods = tree.query_ball_point(samples, float(radii.max()))
            for sample, neighbors in zip(samples, neighborhoods):
                candidates = np.asarray(
                    [neighbor for neighbor in neighbors if neighbor != atom_index], dtype=int
                )
                if not len(candidates):
                    accessible += 1
                    continue
                distances = np.linalg.norm(coordinates[candidates] - sample, axis=1)
                accessible += int(not np.any(distances < radii[candidates]))
            sasa += 4.0 * math.pi * radius * radius * accessible / sphere_points
        output.append(sasa)
    return np.asarray(output, dtype=float)


def surface_model_from_mmcif(path: Path) -> SurfaceModel:
    model = parse_mmcif(path)
    if set(model) != {"A", "B", "C"}:
        raise ValueError(f"expected exact A/B/C pMHC chain layout in {path}")
    peptide = model["C"]
    peptide_ca = ca_coordinates(peptide)
    centroids = np.vstack([sidechain_centroid(residue) for residue in peptide])
    return SurfaceModel(
        sequence="".join(str(residue["aa"]) for residue in peptide),
        groove_ca=np.vstack([ca_coordinates(model[chain][:85]) for chain in ("A", "B")]),
        peptide_ca=peptide_ca,
        sidechain_centroids=centroids,
        sidechain_orientations=centroids - peptide_ca,
        sidechain_sasa=_sidechain_sasa(peptide, [model["A"], model["B"], peptide]),
    )


def _cached_surface_model(path: Path, cache_root: Path) -> tuple[SurfaceModel, str]:
    """Load an exact-input surface cache or build it from the source mmCIF."""
    source_sha256 = sha256_file(path)
    cache_path = cache_root / f"{source_sha256}.npz"
    if cache_path.exists():
        with np.load(cache_path, allow_pickle=False) as cached:
            return SurfaceModel(
                sequence=str(cached["sequence"].item()),
                groove_ca=cached["groove_ca"],
                peptide_ca=cached["peptide_ca"],
                sidechain_centroids=cached["sidechain_centroids"],
                sidechain_orientations=cached["sidechain_orientations"],
                sidechain_sasa=cached["sidechain_sasa"],
            ), source_sha256
    model = surface_model_from_mmcif(path)
    cache_root.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        sequence=np.asarray(model.sequence),
        groove_ca=model.groove_ca,
        peptide_ca=model.peptide_ca,
        sidechain_centroids=model.sidechain_centroids,
        sidechain_orientations=model.sidechain_orientations,
        sidechain_sasa=model.sidechain_sasa,
    )
    return model, source_sha256


def slice_surface_model(model: SurfaceModel, start_1_based: int) -> SurfaceModel:
    start = start_1_based - 1
    stop = start + 9
    if start < 0 or stop > len(model.sequence):
        raise ValueError("register is not fully contained in modeled peptide")
    return SurfaceModel(
        sequence=model.sequence[start:stop],
        groove_ca=model.groove_ca,
        peptide_ca=model.peptide_ca[start:stop],
        sidechain_centroids=model.sidechain_centroids[start:stop],
        sidechain_orientations=model.sidechain_orientations[start:stop],
        sidechain_sasa=model.sidechain_sasa[start:stop],
    )


def _cached_surface_slice(model: SurfaceModel, start_1_based: int) -> SurfaceModel:
    key = (id(model), start_1_based)
    if key not in _SURFACE_SLICE_CACHE:
        _SURFACE_SLICE_CACHE[key] = slice_surface_model(model, start_1_based)
    return _SURFACE_SLICE_CACHE[key]


def _weighted_rmsd(left: np.ndarray, right: np.ndarray, weights: np.ndarray) -> float:
    squared = np.sum((np.asarray(left) - np.asarray(right)) ** 2, axis=1)
    return float(np.sqrt(np.sum(weights * squared)))


def _exposure_weights(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    combined = np.maximum(0.0, np.asarray(left)) + np.maximum(0.0, np.asarray(right))
    if float(combined.sum()) <= 0.0:
        return np.full(len(combined), 1.0 / len(combined))
    return combined / combined.sum()


def _unit_vectors(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return np.divide(values, norms, out=np.zeros_like(values), where=norms > 1e-12)


@lru_cache(maxsize=20)
def _chemistry_descriptor(aa: str) -> np.ndarray:
    # Scaled formal charge, Kyte-Doolittle hydropathy, donors, acceptors, aromaticity.
    charge = {**{aa: 0.0 for aa in "ACFGILMNPQSTVWY"}, "D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
    hydro = dict(zip("ACDEFGHIKLMNPQRSTVWY", (1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3)))
    donor = {aa: float(aa in "HKNRQSTWY") for aa in hydro}
    acceptor = {aa: float(aa in "DEHNQSTY") for aa in hydro}
    aromatic = {aa: float(aa in "FWY") for aa in hydro}
    if aa not in hydro:
        raise ValueError(f"unsupported amino acid {aa}")
    return np.asarray(((charge[aa] + 1) / 2, (hydro[aa] + 4.5) / 9, donor[aa], acceptor[aa], aromatic[aa]))


def surface_pair_metrics(left: SurfaceModel, right: SurfaceModel) -> dict[str, float]:
    """Compare two exact-register local surfaces after equivalent-groove fitting."""
    if len(left.sequence) != 9 or len(right.sequence) != 9:
        raise ValueError("surface comparison requires exact nine-residue cores")
    if left.groove_ca.shape != right.groove_ca.shape:
        raise ValueError("equivalent HLA-II groove arrays are required")
    alignment_key = (id(left.groove_ca), id(right.groove_ca))
    if alignment_key not in _GROOVE_ALIGNMENT_CACHE:
        _GROOVE_ALIGNMENT_CACHE[alignment_key] = kabsch(right.groove_ca, left.groove_ca)
    rotation, translation = _GROOVE_ALIGNMENT_CACHE[alignment_key]
    right_ca = right.peptide_ca @ rotation + translation
    right_centroids = right.sidechain_centroids @ rotation + translation
    right_orientations = right.sidechain_orientations @ rotation
    weights = _exposure_weights(left.sidechain_sasa, right.sidechain_sasa)
    left_relative = left.sidechain_sasa / max(1e-12, float(left.sidechain_sasa.sum()))
    right_relative = right.sidechain_sasa / max(1e-12, float(right.sidechain_sasa.sum()))
    left_distance = np.linalg.norm(left.sidechain_centroids[:, None, :] - left.sidechain_centroids[None, :, :], axis=2)
    right_distance = np.linalg.norm(right_centroids[:, None, :] - right_centroids[None, :, :], axis=2)
    pair_weights = np.outer(weights, weights)
    np.fill_diagonal(pair_weights, 0.0)
    pair_weights /= max(1e-12, float(pair_weights.sum()))
    chemistry = np.asarray([
        np.mean(np.abs(_chemistry_descriptor(a) - _chemistry_descriptor(b)))
        for a, b in zip(left.sequence, right.sequence)
    ])
    fitted_groove = right.groove_ca @ rotation + translation
    return {
        "hla_groove_ca_rmsd_A": float(np.sqrt(np.mean(np.sum((left.groove_ca - fitted_groove) ** 2, axis=1)))),
        "sidechain_exposure_mismatch_fraction": float(np.sum(np.abs(left_relative - right_relative)) / 2.0),
        "exposure_weighted_centroid_rmsd_A": _weighted_rmsd(left.sidechain_centroids, right_centroids, weights),
        "exposure_weighted_orientation_rmsd_A": _weighted_rmsd(_unit_vectors(left.sidechain_orientations), _unit_vectors(right_orientations), weights),
        "exposed_distance_matrix_rmsd_A": float(np.sqrt(np.sum(pair_weights * (left_distance - right_distance) ** 2))),
        "exposure_weighted_chemistry_mismatch": float(np.sum(weights * chemistry)),
        "exposure_weighted_backbone_rmsd_A": _weighted_rmsd(left.peptide_ca, right_ca, weights),
        "full_core_ca_rmsd_A": float(np.sqrt(np.mean(np.sum((left.peptide_ca - right_ca) ** 2, axis=1)))),
        "anchor_ca_rmsd_A": float(np.sqrt(np.mean(np.sum((left.peptide_ca[list(ANCHOR_INDICES)] - right_ca[list(ANCHOR_INDICES)]) ** 2, axis=1)))),
        "left_exposure_weights": ";".join(f"{value:.6f}" for value in left_relative),
        "right_exposure_weights": ";".join(f"{value:.6f}" for value in right_relative),
    }


def enumerate_register_windows(sequence: str) -> list[tuple[int, str]]:
    if len(sequence) < 9:
        return []
    return [(start + 1, sequence[start : start + 9]) for start in range(len(sequence) - 8)]


def allowed_register_starts(sequence_length: int, declared_start_1_based: int) -> tuple[int, ...]:
    maximum = sequence_length - 8
    if not 1 <= declared_start_1_based <= maximum:
        raise ValueError("declared register is not fully contained")
    return tuple(
        start
        for start in range(declared_start_1_based - 1, declared_start_1_based + 2)
        if 1 <= start <= maximum
    )


def assign_evidence_tier(
    sequence_percentile: float | None,
    structure_percentile: float | None,
    register_robust: bool,
) -> str:
    if sequence_percentile is None or structure_percentile is None or not register_robust:
        return "M"
    if sequence_percentile <= 0.10 and structure_percentile <= 0.10:
        return "A"
    if sequence_percentile <= 0.25 and structure_percentile <= 0.25:
        return "B"
    if sequence_percentile <= 0.25:
        return "C"
    if structure_percentile <= 0.25:
        return "D"
    return "E"


def _number(row: Mapping[str, Any], field: str, default: float = math.inf) -> float:
    value = row.get(field, "")
    if value is None or str(value).strip() == "":
        return default
    number = float(value)
    return number if math.isfinite(number) else default


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def rank_v3_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Apply the locked lexicographic V3 ranking separately within each HLA."""
    by_hla: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen: set[str] = set()
    for source in rows:
        row = dict(source)
        pair_id = str(row["pair_id"])
        if pair_id in seen:
            raise ValueError(f"duplicate pair ID {pair_id}")
        seen.add(pair_id)
        allele = str(row["allele"])
        if row.get("ebv_allele", allele) != allele or row.get("self_allele", allele) != allele:
            raise ValueError(f"pair {pair_id} is not exact-HLA")
        by_hla[allele].append(row)
    output: list[dict[str, Any]] = []
    for allele in sorted(by_hla):
        ordered = sorted(
            by_hla[allele],
            key=lambda row: (
                -_number(row, "tcr_facing_blosum62_similarity", -math.inf),
                _number(row, "tcr_face_physicochemical_mismatch"),
                -_number(row, "tcr_facing_sequence_identity", -math.inf),
                _number(
                    row,
                    "ranking_local_surface_percentile"
                    if "ranking_local_surface_percentile" in row
                    else "local_surface_percentile",
                ),
                str(row["pair_id"]),
            ),
        )
        for index, row in enumerate(ordered, start=1):
            row["primary_rank"] = index
            row["rank_scope"] = "within_hla_only"
            row["primary_method"] = "tcr_facing_blosum62_lexicographic_v3"
            row["deterministic_tie_break"] = "physicochemical_then_identity_then_local_surface_then_pair_id"
            output.append(row)
    return output


def _average_tie_percentiles(
    rows: Sequence[Mapping[str, Any]], field: str, *, higher_is_better: bool = False
) -> dict[str, float]:
    eligible = [row for row in rows if math.isfinite(_number(row, field))]
    ordered = sorted(
        eligible,
        key=lambda row: (
            -_number(row, field) if higher_is_better else _number(row, field),
            str(row["pair_id"]),
        ),
    )
    denominator = max(1, len(ordered) - 1)
    output: dict[str, float] = {}
    position = 0
    while position < len(ordered):
        end = position + 1
        value = _number(ordered[position], field)
        while end < len(ordered) and _number(ordered[end], field) == value:
            end += 1
        percentile = ((position + end - 1) / 2.0) / denominator
        for row in ordered[position:end]:
            output[str(row["pair_id"])] = round(percentile, 12)
        position = end
    return output


def _empirical_percentile(value: float, reference: Sequence[float], *, higher: bool = False) -> float:
    ordered = sorted(float(item) for item in reference if math.isfinite(float(item)))
    if not ordered:
        return math.nan
    if higher:
        better = sum(item > value for item in ordered)
        equal = sum(item == value for item in ordered)
    else:
        better = sum(item < value for item in ordered)
        equal = sum(item == value for item in ordered)
    return ((better + (equal - 1) / 2.0) / max(1, len(ordered) - 1))


def _quantile_summary(values: Sequence[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    if not len(array):
        raise ValueError("cannot summarize an empty ensemble")
    q25, median, q75 = np.quantile(array, [0.25, 0.5, 0.75])
    return {
        "min": float(array.min()),
        "q25": float(q25),
        "median": float(median),
        "q75": float(q75),
        "max": float(array.max()),
        "iqr": float(q75 - q25),
    }


def _pair_surface_summary(
    left_models: Sequence[SurfaceModel],
    right_models: Sequence[SurfaceModel],
    left_start: int,
    right_start: int,
) -> dict[str, Any]:
    if not left_models or not right_models:
        return {"surface_status": "missing_model", "model_combination_count": 0}
    left_cores = [_cached_surface_slice(model, left_start) for model in left_models]
    right_cores = [_cached_surface_slice(model, right_start) for model in right_models]
    combinations = [surface_pair_metrics(left, right) for left in left_cores for right in right_cores]
    output: dict[str, Any] = {
        "surface_status": "complete",
        "left_model_count": len(left_cores),
        "right_model_count": len(right_cores),
        "model_combination_count": len(combinations),
    }
    for feature in ("hla_groove_ca_rmsd_A", *SURFACE_FEATURES, "full_core_ca_rmsd_A", "anchor_ca_rmsd_A"):
        for statistic, value in _quantile_summary([float(row[feature]) for row in combinations]).items():
            output[f"{feature}_{statistic}"] = round(value, 8)
    left_exposure = np.mean(
        [_cached_surface_slice(model, left_start).sidechain_sasa for model in left_models], axis=0
    )
    right_exposure = np.mean(
        [_cached_surface_slice(model, right_start).sidechain_sasa for model in right_models], axis=0
    )
    output["left_mean_sidechain_sasa_P1_P9_A2"] = ";".join(f"{value:.4f}" for value in left_exposure)
    output["right_mean_sidechain_sasa_P1_P9_A2"] = ";".join(f"{value:.4f}" for value in right_exposure)
    return output


def _candidate_register_start(sequence: str, core: str) -> tuple[int | None, str]:
    starts = [start for start, window in enumerate_register_windows(sequence) if window == core]
    if len(starts) == 1:
        return starts[0], "unique_exact_core_in_peptide"
    if not starts:
        return None, "declared_core_not_found_in_peptide"
    return None, "declared_core_occurs_multiple_times"


def _load_v2_models(
    sample_path: Path,
    cache_root: Path,
) -> tuple[dict[tuple[str, str], list[SurfaceModel]], list[dict[str, Any]]]:
    groups: dict[tuple[str, str], list[SurfaceModel]] = defaultdict(list)
    audit = []
    for row in read_csv(sample_path):
        if row.get("cohort") != "discovery" or row.get("sample_status") != "pass_exact_clash_free":
            continue
        path = Path(row["cif_path"])
        model, model_sha256 = _cached_surface_model(path, cache_root)
        if model.sequence != row["observed_peptide"]:
            raise ValueError(f"V2 sequence mismatch for {path}")
        key = (row["allele"], row["entity_id"])
        groups[key].append(model)
        audit.append({
            "library": "tcell_library_v2",
            "allele": row["allele"],
            "candidate_id": row["entity_id"],
            "sample_index": row["sample_index"],
            "model_path": str(path),
            "model_sha256": model_sha256,
            "sequence": model.sequence,
            "surface_status": "parsed_complete",
        })
    return groups, audit


def _load_legacy_models(
    sample_path: Path,
    cache_root: Path,
) -> tuple[dict[str, list[SurfaceModel]], list[dict[str, Any]]]:
    groups: dict[str, list[SurfaceModel]] = defaultdict(list)
    audit = []
    for row in read_csv(sample_path):
        if row.get("sequence_layout_status") != "pass_exact_three_chain_peptide_match":
            continue
        path = Path(row["model_path"])
        model, model_sha256 = _cached_surface_model(path, cache_root)
        if model.sequence != row["observed_peptide"]:
            raise ValueError(f"legacy sequence mismatch for {path}")
        groups[row["candidate_id"]].append(model)
        audit.append({
            "library": "legacy_drb1501",
            "allele": "HLA-DRB1*15:01",
            "candidate_id": row["candidate_id"],
            "sample_index": row["sample_index"],
            "model_path": str(path),
            "model_sha256": model_sha256,
            "sequence": model.sequence,
            "surface_status": "parsed_complete",
        })
    return groups, audit


def _normalize_v2_row(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        **dict(row),
        "ebv_core_p1_p9": row["ebv_predicted_core"],
        "self_core_p1_p9": row["self_predicted_core"],
        "pair_coordinate_label": row.get("pair_coordinate_label", ""),
        "source_membership": "v2",
    }


def _declared_registers(
    row: Mapping[str, Any],
    v2_registers: Mapping[tuple[str, str], Mapping[str, Any]],
) -> tuple[int | None, int | None, bool, str]:
    allele = str(row["allele"])
    source = str(row.get("source_membership", "v2"))
    if source in {"v2", "v2_only", "v2_and_legacy_exact_duplicate"}:
        left = v2_registers.get((allele, str(row["ebv_candidate_id"])))
        right = v2_registers.get((allele, str(row["self_candidate_id"])))
        if not left or not right:
            return None, None, False, "missing_v2_register_record"
        statuses = {left.get("register_resolution"), right.get("register_resolution")}
        unique = statuses == {"resolved_unique_fully_contained"}
        return int(left["core_start"]), int(right["core_start"]), unique, (
            "iedb_resolved_unique_fully_contained" if unique else "v2_register_not_uniquely_resolved"
        )
    left_start, left_status = _candidate_register_start(str(row["ebv_sequence"]), str(row["ebv_core_p1_p9"]))
    right_start, right_status = _candidate_register_start(str(row["self_sequence"]), str(row["self_core_p1_p9"]))
    # Legacy registers were used for sensitivity coverage but lack a frozen V3 IEDB
    # unique-resolution record, so they cannot meet the strict robustness definition.
    return left_start, right_start, False, f"legacy_not_iedb_unique:{left_status};{right_status}"


def _models_for_pair(
    row: Mapping[str, Any],
    v2_models: Mapping[tuple[str, str], Sequence[SurfaceModel]],
    legacy_models: Mapping[str, Sequence[SurfaceModel]],
) -> tuple[Sequence[SurfaceModel], Sequence[SurfaceModel], str]:
    source = str(row.get("source_membership", "v2"))
    if source in {"v2", "v2_only", "v2_and_legacy_exact_duplicate"}:
        allele = str(row["allele"])
        return (
            v2_models.get((allele, str(row["ebv_candidate_id"])), ()),
            v2_models.get((allele, str(row["self_candidate_id"])), ()),
            "tcell_library_v2_models",
        )
    return (
        legacy_models.get(str(row["ebv_candidate_id"]), ()),
        legacy_models.get(str(row["self_candidate_id"]), ()),
        "legacy_drb1501_models",
    )


def _base_tier(sequence_percentile: float, structure_percentile: float, unstable: bool) -> str:
    if unstable and sequence_percentile <= 0.25:
        return "C"
    return assign_evidence_tier(sequence_percentile, structure_percentile, True)


def _annotate_declared_surface(
    rows: Sequence[Mapping[str, Any]],
    v2_registers: Mapping[tuple[str, str], Mapping[str, Any]],
    v2_models: Mapping[tuple[str, str], Sequence[SurfaceModel]],
    legacy_models: Mapping[str, Sequence[SurfaceModel]],
) -> list[dict[str, Any]]:
    annotated = []
    for source in rows:
        row = dict(source)
        left_start, right_start, registry_unique, register_status = _declared_registers(row, v2_registers)
        row.update({
            "ebv_declared_core_start_1_based": left_start if left_start is not None else "",
            "self_declared_core_start_1_based": right_start if right_start is not None else "",
            "declared_register_registry_unique": registry_unique,
            "declared_register_status": register_status,
        })
        left_models, right_models, model_source = _models_for_pair(row, v2_models, legacy_models)
        row["surface_model_source"] = model_source
        if left_start is None or right_start is None:
            row.update({"surface_status": "register_unavailable", "model_combination_count": 0})
        else:
            row.update(_pair_surface_summary(left_models, right_models, left_start, right_start))
        annotated.append(row)

    for allele in sorted({str(row["allele"]) for row in annotated}):
        group = [row for row in annotated if row["allele"] == allele]
        complete = [row for row in group if row["surface_status"] == "complete"]
        feature_maps = {
            feature: _average_tie_percentiles(complete, f"{feature}_q75")
            for feature in SURFACE_FEATURES
        }
        sequence_map = _average_tie_percentiles(
            group, "tcr_facing_blosum62_similarity", higher_is_better=True
        )
        for row in group:
            pair_id = str(row["pair_id"])
            row["sequence_family_percentile"] = sequence_map[pair_id]
            if row["surface_status"] != "complete":
                row["local_surface_percentile"] = ""
                row["surface_ensemble_uncertainty"] = ""
                continue
            values = []
            for feature in SURFACE_FEATURES:
                percentile = feature_maps[feature][pair_id]
                row[f"{feature}_q75_percentile"] = percentile
                values.append(percentile)
            row["local_surface_percentile"] = round(float(np.mean(values)), 12)
            normalized_iqrs = [
                _number(row, f"{feature}_iqr") / max(1e-8, _number(row, f"{feature}_q75"))
                for feature in SURFACE_FEATURES
            ]
            row["surface_ensemble_uncertainty"] = round(float(np.mean(normalized_iqrs)), 12)
        complete_uncertainty = [row for row in complete if str(row.get("surface_ensemble_uncertainty", ""))]
        uncertainty_map = _average_tie_percentiles(
            complete_uncertainty, "surface_ensemble_uncertainty"
        )
        for row in group:
            pair_id = str(row["pair_id"])
            if pair_id in uncertainty_map:
                row["surface_ensemble_uncertainty_percentile"] = uncertainty_map[pair_id]
                row["surface_ensemble_stable"] = uncertainty_map[pair_id] <= 0.75
            else:
                row["surface_ensemble_uncertainty_percentile"] = ""
                row["surface_ensemble_stable"] = False
    return annotated


def _register_sensitivity(
    rows: Sequence[dict[str, Any]],
    v2_models: Mapping[tuple[str, str], Sequence[SurfaceModel]],
    legacy_models: Mapping[str, Sequence[SurfaceModel]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    by_hla = defaultdict(list)
    for row in rows:
        by_hla[str(row["allele"])].append(row)
    references: dict[str, dict[str, list[float]]] = {}
    for allele, group in by_hla.items():
        complete = [row for row in group if row["surface_status"] == "complete"]
        references[allele] = {
            "blosum": [_number(row, "tcr_facing_blosum62_similarity") for row in group],
            **{
                feature: [_number(row, f"{feature}_q75") for row in complete]
                for feature in SURFACE_FEATURES
            },
        }

    summary_rows: list[dict[str, Any]] = []
    allowed_rows: list[dict[str, Any]] = []
    for row in rows:
        allele = str(row["allele"])
        left_sequence, right_sequence = str(row["ebv_sequence"]), str(row["self_sequence"])
        all_sequence_scores = [
            sequence_metrics(left_core, right_core)["tcr_facing_blosum62_similarity"]
            for _left_start, left_core in enumerate_register_windows(left_sequence)
            for _right_start, right_core in enumerate_register_windows(right_sequence)
        ]
        declared_sequence_percentile = float(row["sequence_family_percentile"])
        structural_available = row["surface_status"] == "complete"
        declared_structure_percentile = (
            float(row["local_surface_percentile"]) if structural_available else math.nan
        )
        declared_unstable = not bool(row.get("surface_ensemble_stable", False))
        declared_tier = (
            _base_tier(declared_sequence_percentile, declared_structure_percentile, declared_unstable)
            if structural_available
            else "M"
        )
        tier_values = []
        left_start = int(row["ebv_declared_core_start_1_based"]) if str(row["ebv_declared_core_start_1_based"]) else None
        right_start = int(row["self_declared_core_start_1_based"]) if str(row["self_declared_core_start_1_based"]) else None
        left_models, right_models, _model_source = _models_for_pair(row, v2_models, legacy_models)
        if left_start is not None and right_start is not None and structural_available:
            for alt_left in allowed_register_starts(len(left_sequence), left_start):
                for alt_right in allowed_register_starts(len(right_sequence), right_start):
                    seq = sequence_metrics(
                        left_sequence[alt_left - 1 : alt_left + 8],
                        right_sequence[alt_right - 1 : alt_right + 8],
                    )
                    surface = _pair_surface_summary(left_models, right_models, alt_left, alt_right)
                    feature_percentiles = [
                        _empirical_percentile(
                            float(surface[f"{feature}_q75"]), references[allele][feature]
                        )
                        for feature in SURFACE_FEATURES
                    ]
                    sequence_percentile = _empirical_percentile(
                        float(seq["tcr_facing_blosum62_similarity"]),
                        references[allele]["blosum"],
                        higher=True,
                    )
                    structure_percentile = float(np.mean(feature_percentiles))
                    tier = _base_tier(sequence_percentile, structure_percentile, False)
                    tier_values.append(tier)
                    allowed_rows.append({
                        "allele": allele,
                        "pair_id": row["pair_id"],
                        "ebv_window_start_1_based": alt_left,
                        "self_window_start_1_based": alt_right,
                        "is_declared_window_pair": alt_left == left_start and alt_right == right_start,
                        "ebv_window_core": left_sequence[alt_left - 1 : alt_left + 8],
                        "self_window_core": right_sequence[alt_right - 1 : alt_right + 8],
                        "tcr_facing_blosum62_similarity": round(float(seq["tcr_facing_blosum62_similarity"]), 12),
                        "sequence_family_percentile": round(sequence_percentile, 12),
                        "local_surface_percentile": round(structure_percentile, 12),
                        "evidence_tier_without_register_gate": tier,
                        "interpretation": "Sensitivity on fixed modeled coordinates; alternate windows are not alternate-register structure predictions.",
                    })
        stable_tier = bool(tier_values) and set(tier_values) == {declared_tier}
        register_robust = bool(row["declared_register_registry_unique"]) and stable_tier
        row["declared_evidence_tier_before_register_gate"] = declared_tier
        row["register_tier_stable_across_allowed_windows"] = stable_tier
        row["register_robust"] = register_robust
        row["evidence_tier"] = (
            assign_evidence_tier(
                declared_sequence_percentile,
                declared_structure_percentile if structural_available else None,
                register_robust,
            )
            if not (structural_available and declared_unstable and declared_sequence_percentile <= 0.25 and register_robust)
            else "C"
        )
        row["ranking_local_surface_percentile"] = (
            row["local_surface_percentile"] if structural_available and register_robust else ""
        )
        row["structural_ranking_status"] = (
            "ranked_local_surface_register_robust"
            if structural_available and register_robust
            else (
                "surface_annotated_not_ranked_register_uncertain"
                if structural_available
                else "not_ranked_missing_surface_model"
            )
        )
        row["missingness_reason"] = "" if structural_available else str(row["surface_status"])
        summary_rows.append({
            "allele": allele,
            "pair_id": row["pair_id"],
            "declared_register_status": row["declared_register_status"],
            "declared_register_registry_unique": row["declared_register_registry_unique"],
            "all_fully_contained_window_pair_count": len(all_sequence_scores),
            "all_window_tcr_facing_blosum62_min": round(min(all_sequence_scores), 12),
            "all_window_tcr_facing_blosum62_max": round(max(all_sequence_scores), 12),
            "allowed_window_pair_count": len(tier_values),
            "allowed_window_evidence_tiers": ";".join(tier_values),
            "declared_evidence_tier_before_register_gate": declared_tier,
            "register_tier_stable_across_allowed_windows": stable_tier,
            "register_robust": register_robust,
            "sensitivity_scope": "all_windows_sequence_and_declared_plus_or_minus_one_surface_on_fixed_models",
        })
    return summary_rows, allowed_rows


def build_v3_universe(
    rows: Sequence[Mapping[str, Any]],
    v2_registers: Mapping[tuple[str, str], Mapping[str, Any]],
    v2_models: Mapping[tuple[str, str], Sequence[SurfaceModel]],
    legacy_models: Mapping[str, Sequence[SurfaceModel]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    annotated = _annotate_declared_surface(rows, v2_registers, v2_models, legacy_models)
    sensitivity, allowed = _register_sensitivity(annotated, v2_models, legacy_models)
    ranked = rank_v3_rows(annotated)
    for row in ranked:
        count = sum(item["allele"] == row["allele"] for item in ranked)
        row["primary_percentile"] = round((int(row["primary_rank"]) - 1) / max(1, count - 1), 12)
        row["binding_percentile_role"] = "reported_presentation_plausibility_only_not_used_in_rank"
        row["computational_pair_marker"] = "*"
        row["pair_evidence_status"] = "computational_pair_no_exact_paired_recognition_evidence"
        row["claim_boundary"] = CLAIM_BOUNDARY
    return ranked, sensitivity, allowed


def _write_union_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table {path}")
    fields: list[str] = []
    seen = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    write_csv(path, rows, fields)


def _development_control_audit(
    feature_path: Path,
    method_rank_path: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    features = read_csv(feature_path)
    output: list[dict[str, Any]] = []
    panel_fields = ("system_id", "positive_pair_id", "panel_seed")
    panels: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in features:
        panels[tuple(str(row[field]) for field in panel_fields)].append(dict(row))
    v3_ranks = {}
    for panel, rows in sorted(panels.items()):
        for structure_field in ("exposed_ca_rmsd_A_q75", "exposed_sidechain_vector_rmsd_A_q75"):
            mapping = _average_tie_percentiles(rows, structure_field)
            for row in rows:
                row.setdefault("control_surface_proxy_percentiles", []).append(mapping[row["pair_id"]])
        for row in rows:
            row["local_surface_percentile"] = float(np.mean(row["control_surface_proxy_percentiles"]))
        ranked = rank_v3_rows([{**row, "allele": "CONTROL_EXACT_HLA"} for row in rows])
        positive = next(row for row in ranked if row["pair_role"] == "positive")
        v3_ranks[panel] = int(positive["primary_rank"])
        output.append({
            "system_id": panel[0],
            "positive_pair_id": panel[1],
            "panel_seed": panel[2],
            "method": "v3_blosum_lexicographic_with_available_surface_proxy",
            "positive_rank": positive["primary_rank"],
            "capture_at_3": int(positive["primary_rank"]) <= 3,
            "comparison_count": len(rows),
            "control_role": "development_only_not_independent_v3_validation",
        })
    existing = read_csv(method_rank_path)
    for row in existing:
        output.append({
            "system_id": row["system_id"],
            "positive_pair_id": row["positive_pair_id"],
            "panel_seed": row["panel_seed"],
            "method": row["method"],
            "positive_rank": int(row["positive_rank"]),
            "capture_at_3": str(row["capture_at_3"]).lower() == "true",
            "comparison_count": int(row["comparison_count"]),
            "control_role": "development_only_not_independent_v3_validation",
        })
    gate = {
        "benchmark_version": "literature_grounded_hla2_ranking_v3_development_audit",
        "status": "development_controls_only_not_independent_validation",
        "independent_system_count": 0,
        "development_system_count": len({panel[0] for panel in panels}),
        "development_panel_count": len(panels),
        "v3_development_capture_at_3_count": sum(rank <= 3 for rank in v3_ranks.values()),
        "v3_development_rank_1_count": sum(rank == 1 for rank in v3_ranks.values()),
        "v3_development_panel_ranks": [
            {"system_id": key[0], "positive_pair_id": key[1], "panel_seed": key[2], "rank": rank}
            for key, rank in sorted(v3_ranks.items())
        ],
        "surface_proxy_limitation": "Existing control matrix lacks V3 solvent-accessibility fingerprints; the V3 primary sequence keys are exact, while the final tie-break uses the available exposed-CA and side-chain-vector q75 proxy.",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return sorted(output, key=lambda row: (row["method"], row["system_id"], row["positive_pair_id"], int(row["panel_seed"]))), gate


def _top_by_hla(rows: Sequence[Mapping[str, Any]], count: int = 10) -> list[dict[str, Any]]:
    by_hla: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        by_hla[str(row["allele"])].append(row)
    return [dict(row) for allele in sorted(by_hla) for row in by_hla[allele][:count]]


def _candidate_selection(rows: Sequence[Mapping[str, Any]], per_hla: int = 5) -> list[dict[str, Any]]:
    output = []
    for allele in sorted({str(row["allele"]) for row in rows}):
        group = [
            row for row in rows
            if row["allele"] == allele
            and row.get("evidence_tier") in {"A", "B", "C"}
            and _truth(row.get("register_robust"))
            and row.get("surface_status") == "complete"
        ]
        group.sort(
            key=lambda row: (
                {"A": 0, "B": 1, "C": 2}[str(row["evidence_tier"])],
                int(row["primary_rank"]),
                max(
                    _number(row, "ebv_binding_percentile_rank"),
                    _number(row, "self_binding_percentile_rank"),
                ),
                str(row["pair_id"]),
            )
        )
        for row in group[:per_hla]:
            selected = dict(row)
            selected["candidate_selection_basis"] = (
                "sequence_strong_structurally_annotated_register_robust; binding percentiles "
                "reported separately for experimental-planning context"
            )
            selected["recommended_followup_scope"] = "experimental shortlist; reserve electrostatics or MD for 10-20 final candidates"
            output.append(selected)
    return output


def _input_record(path: Path, role: str) -> dict[str, Any]:
    return {"role": role, "path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}


def _checksum_table(out: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.csv":
            rows.append({"relative_path": str(path.relative_to(out)), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    return rows


def _cache_worker(arguments: tuple[str, str]) -> tuple[str, int]:
    path, cache_root = arguments
    model, digest = _cached_surface_model(Path(path), Path(cache_root))
    return digest, len(model.sequence)


def prewarm_surface_cache(
    *,
    v2_samples_path: Path = DEFAULT_V2_SAMPLES,
    legacy_samples_path: Path = DEFAULT_LEGACY_SAMPLES,
    cache_root: Path = DEFAULT_CACHE,
    workers: int = 4,
) -> dict[str, int]:
    paths = [
        row["cif_path"]
        for row in read_csv(v2_samples_path)
        if row.get("cohort") == "discovery" and row.get("sample_status") == "pass_exact_clash_free"
    ]
    paths.extend(
        row["model_path"]
        for row in read_csv(legacy_samples_path)
        if row.get("sequence_layout_status") == "pass_exact_three_chain_peptide_match"
    )
    with multiprocessing.Pool(processes=workers) as pool:
        results = list(
            pool.imap_unordered(
                _cache_worker,
                [(path, str(cache_root)) for path in paths],
                chunksize=4,
            )
        )
    return {
        "requested_model_count": len(paths),
        "cached_model_count": len(results),
        "unique_model_sha256_count": len({digest for digest, _length in results}),
    }


def _apply_register_aware_structural_abstention(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    updated = []
    for source in rows:
        row = dict(source)
        complete = row.get("surface_status") == "complete"
        robust = _truth(row.get("register_robust"))
        row["ranking_local_surface_percentile"] = (
            row.get("local_surface_percentile", "") if complete and robust else ""
        )
        row["structural_ranking_status"] = (
            "ranked_local_surface_register_robust"
            if complete and robust
            else (
                "surface_annotated_not_ranked_register_uncertain"
                if complete
                else "not_ranked_missing_surface_model"
            )
        )
        updated.append(row)
    return rank_v3_rows(updated)


def finalize_existing_package(out: Path = DEFAULT_OUT) -> dict[str, Any]:
    """Apply the locked M-tier structural abstention to existing computed features."""
    v2 = _apply_register_aware_structural_abstention(read_csv(out / "v3_all_hla_ranked_pairs.csv"))
    combined = _apply_register_aware_structural_abstention(
        read_csv(out / "combined_drb1501_v3_ranked_pairs.csv")
    )
    for universe in (v2, combined):
        counts = defaultdict(int)
        for row in universe:
            counts[str(row["allele"])] += 1
        for row in universe:
            row["primary_percentile"] = round(
                (int(row["primary_rank"]) - 1) / max(1, counts[str(row["allele"])] - 1), 12
            )
    _write_union_csv(out / "v3_all_hla_ranked_pairs.csv", v2)
    for allele in ALLELES:
        _write_union_csv(
            out / "rankings" / f"{ALLELE_SLUGS[allele]}_ranked_pairs.csv",
            [row for row in v2 if row["allele"] == allele],
        )
    _write_union_csv(out / "v3_top_10_by_hla.csv", _top_by_hla(v2))
    _write_union_csv(out / "combined_drb1501_v3_ranked_pairs.csv", combined)
    _write_union_csv(out / "combined_drb1501_top_10_v3.csv", combined[:10])
    _write_union_csv(out / "candidate_selection_table.csv", _candidate_selection(v2))
    _write_union_csv(
        out / "missing_or_register_uncertain_v2.csv",
        [row for row in v2 if row.get("evidence_tier") == "M"],
    )
    manifest = json.loads((out / "analysis_manifest.json").read_text(encoding="utf-8"))
    manifest["candidate_selection_count"] = len(_candidate_selection(v2))
    manifest["register_uncertain_surface_tie_break_abstention"] = True
    write_json(out / "analysis_manifest.json", manifest)
    _write_union_csv(out / "SHA256SUMS.csv", _checksum_table(out))
    return {
        "v2_pair_count": len(v2),
        "combined_pair_count": len(combined),
        "register_robust_structurally_ranked_v2_count": sum(
            row["structural_ranking_status"] == "ranked_local_surface_register_robust"
            for row in v2
        ),
        "checksum_artifact_count": len(read_csv(out / "SHA256SUMS.csv")),
    }


def run(
    *,
    v2_sequence_path: Path = DEFAULT_V2_SEQUENCE,
    combined_sequence_path: Path = DEFAULT_COMBINED_SEQUENCE,
    v2_samples_path: Path = DEFAULT_V2_SAMPLES,
    v2_registers_path: Path = DEFAULT_V2_REGISTERS,
    legacy_samples_path: Path = DEFAULT_LEGACY_SAMPLES,
    control_features_path: Path = DEFAULT_CONTROL_FEATURES,
    out: Path = DEFAULT_OUT,
    cache_root: Path = DEFAULT_CACHE,
) -> dict[str, Any]:
    """Build the additive V3 package without modifying any input package."""
    out.mkdir(parents=True, exist_ok=True)
    v2_source = [_normalize_v2_row(row) for row in read_csv(v2_sequence_path)]
    combined_source = [dict(row) for row in read_csv(combined_sequence_path)]
    if len(v2_source) != 6400 or len(combined_source) != 2043:
        raise ValueError("V3 requires the frozen 6,400-pair V2 and 2,043-pair combined universes")
    register_rows = read_csv(v2_registers_path)
    v2_registers = {(row["allele"], row["candidate_id"]): row for row in register_rows}
    if len(v2_registers) != 320:
        raise ValueError("expected 320 frozen V2 allele-register records")

    v2_models, v2_model_audit = _load_v2_models(v2_samples_path, cache_root)
    legacy_models, legacy_model_audit = _load_legacy_models(legacy_samples_path, cache_root)
    v2_ranked, v2_sensitivity, v2_allowed = build_v3_universe(
        v2_source, v2_registers, v2_models, legacy_models
    )
    combined_ranked, combined_sensitivity, combined_allowed = build_v3_universe(
        combined_source, v2_registers, v2_models, legacy_models
    )

    _write_union_csv(out / "v3_all_hla_ranked_pairs.csv", v2_ranked)
    for allele in ALLELES:
        _write_union_csv(
            out / "rankings" / f"{ALLELE_SLUGS[allele]}_ranked_pairs.csv",
            [row for row in v2_ranked if row["allele"] == allele],
        )
    _write_union_csv(out / "v3_top_10_by_hla.csv", _top_by_hla(v2_ranked))
    _write_union_csv(out / "combined_drb1501_v3_ranked_pairs.csv", combined_ranked)
    _write_union_csv(out / "combined_drb1501_top_10_v3.csv", combined_ranked[:10])
    _write_union_csv(out / "candidate_selection_table.csv", _candidate_selection(v2_ranked))
    _write_union_csv(out / "register_sensitivity_summary_v2.csv", v2_sensitivity)
    _write_union_csv(out / "register_sensitivity_allowed_windows_v2.csv", v2_allowed)
    _write_union_csv(out / "register_sensitivity_summary_combined_drb1501.csv", combined_sensitivity)
    _write_union_csv(out / "register_sensitivity_allowed_windows_combined_drb1501.csv", combined_allowed)
    _write_union_csv(out / "model_surface_audit.csv", v2_model_audit + legacy_model_audit)
    missing = [row for row in v2_ranked if row["evidence_tier"] == "M"]
    _write_union_csv(out / "missing_or_register_uncertain_v2.csv", missing)

    control_methods = control_features_path.parent / "method_rank_long.csv"
    control_rows, control_gate = _development_control_audit(control_features_path, control_methods)
    _write_union_csv(out / "development_control_method_comparison.csv", control_rows)
    write_json(out / "development_control_audit.json", control_gate)
    write_json(out / "definitive_validation_gate.json", {
        "status": "not_evaluable_insufficient_untouched_strict_systems",
        "required_independent_system_count": 6,
        "target_independent_system_count": 8,
        "required_hla_family_count": 2,
        "current_independent_v3_validation_system_count": 0,
        "current_three_system_role": "development_controls_only",
        "required_rule": "all panels capture_at_3; majority systems improve over BLOSUM62; none worsen; strictly better system-weighted reciprocal rank",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    })
    write_json(out / "specificity_gate.json", {
        "status": "not_evaluable_no_new_v3_n1_n2_registry",
        "ranking_comparators_are_not_specificity_negatives": True,
        "n3_excluded_from_specificity": True,
        "claim_boundary": "Ranking performance is not specificity evidence.",
    })
    write_json(out / "pdb_oracle_availability.json", {
        "hy1b11_status": "not_evaluable_availability",
        "reason": "Fewer than five frozen unique exact-HLA experimental structural decoys; proposed 5JN5 and 5JN6 are unrelated structures.",
        "pdb_5jn5_identity": "human phosphoglucomutase-1 D263Y",
        "pdb_5jn6_identity": "RPA3313 solution NMR structure",
        "availability_is_not_a_pass": True,
    })

    report_paths = [
        ROOT / "../Downloads/A Blueprint for Definitive Validation_ Advancing a Structural HLA-II pMHC Model Beyond Sequence-Based Baselines.pdf",
        ROOT / "../Downloads/Elicit - EBV–MS HLA-II pMHC Validation and Benchmark-v2 Report.pdf",
    ]
    # Use the exact user-supplied iCloud locations when present.
    report_paths = [
        Path("/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads/A Blueprint for Definitive Validation_ Advancing a Structural HLA-II pMHC Model Beyond Sequence-Based Baselines.pdf"),
        Path("/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/Downloads/Elicit - EBV–MS HLA-II pMHC Validation and Benchmark-v2 Report.pdf"),
    ]
    protocol = {
        "protocol_version": "LITERATURE_GROUNDED_HLA2_RANKING_V3_2026-08-27",
        "frozen_before_v3_geometry_read": True,
        "ranking_keys_in_order": [
            "tcr_facing_blosum62_similarity_desc",
            "tcr_face_physicochemical_mismatch_asc",
            "tcr_facing_sequence_identity_desc",
            "local_surface_percentile_asc",
            "pair_id_lexical_asc",
        ],
        "surface_features": list(SURFACE_FEATURES),
        "surface_ensemble_summary": "model_combination_q75",
        "solvent_exposure_method": "deterministic_32_point_shrake_rupley_probe_1.4_A",
        "register_sensitivity": "all fully contained P1-P9 windows for sequence; declared plus/minus one for surface using fixed coordinates",
        "binding_percentile_used_in_rank": False,
        "current_controls_role": "development_only",
        "future_validation_minimum_systems": 6,
        "future_validation_target_systems": 8,
        "cross_allele_consensus_created": False,
        "protein_language_model_used": False,
        "tcr_docking_used": False,
        "research_reports_are_inputs_not_instructions": [
            _input_record(path, "research_input_only") for path in report_paths if path.exists()
        ],
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "platform": platform.platform(),
        },
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "protocol_lock.json", protocol)

    input_records = [
        _input_record(v2_sequence_path, "frozen_v2_sequence_universe"),
        _input_record(combined_sequence_path, "frozen_combined_drb1501_universe"),
        _input_record(v2_samples_path, "v2_model_inventory"),
        _input_record(v2_registers_path, "v2_register_registry"),
        _input_record(legacy_samples_path, "legacy_model_inventory"),
        _input_record(control_features_path, "development_control_features"),
        _input_record(control_methods, "development_control_method_ranks"),
    ]
    manifest = {
        "package": "literature_grounded_hla2_rankings_v3",
        "status": "complete_additive_descriptive_ranking",
        "v2_pair_count": len(v2_ranked),
        "v2_surface_complete_count": sum(row["surface_status"] == "complete" for row in v2_ranked),
        "v2_evidence_tier_counts": {tier: sum(row["evidence_tier"] == tier for row in v2_ranked) for tier in "ABCDE" + "M"},
        "combined_drb1501_pair_count": len(combined_ranked),
        "combined_surface_complete_count": sum(row["surface_status"] == "complete" for row in combined_ranked),
        "candidate_selection_count": len(_candidate_selection(v2_ranked)),
        "input_records": input_records,
        "existing_packages_modified": False,
        "read_through_surface_cache": str(cache_root),
        "cache_is_acceleration_only_not_an_analysis_output": True,
        "discovery_unlock_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "analysis_manifest.json", manifest)

    readme = f"""# Literature-Grounded HLA-II Ranking V3

This additive package ranks {len(v2_ranked):,} V2 pairs separately within four HLA-DR alleles and ranks the frozen {len(combined_ranked):,}-pair DRB1*15:01 combined universe. Existing rankings and benchmark outputs were not modified.

The primary key is TCR-facing P2/P3/P5/P7/P8 BLOSUM62. Physicochemical mismatch, identity, local modeled surface, and lexical pair ID break ties in that order. Local surface is a separate model-derived annotation summarized conservatively at the 75th percentile across model combinations. Binding percentile is reported separately and never enters the rank.

Evidence tier `M` means missing structure or failure of the strict register-robustness rule; it does not remove a pair from the primary sequence rank. The existing three systems are development controls only. `definitive_validation_gate.json` therefore remains not evaluable and discovery unlock is false.

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")
    methods = f"""# Methods

Equivalent HLA-II grooves were aligned by Kabsch fitting of the first 85 C-alpha atoms from each alpha and beta chain. Each peptide P1-P9 surface fingerprint contains side-chain solvent accessibility, centroid position, centroid orientation, the exposed-residue distance matrix, scaled charge/hydropathy/donor/acceptor/aromaticity mismatch, and exposed-backbone geometry. Glycine uses its C-alpha atom as the side-chain fallback. Solvent accessibility uses a deterministic 32-point Shrake-Rupley approximation with a 1.4 A probe.

Every structural feature is summarized by its model-combination 75th percentile. Raw q25, median, q75, IQR, extrema, and per-feature within-HLA percentiles are retained. Register sensitivity enumerates all contained nine-residue window pairs for sequence and evaluates declared +/- 1 windows on the fixed modeled coordinates. These are sensitivity calculations, not alternative-register structure predictions.

Primary ranking is lexicographic: BLOSUM62 descending, physicochemical mismatch ascending, identity descending, local-surface percentile ascending, pair ID ascending. Binding predictions, full-core RMSD, and anchor RMSD are diagnostics only.

{CLAIM_BOUNDARY}
"""
    (out / "METHODS.md").write_text(methods, encoding="utf-8")
    results = f"""# Results Summary

- V2 pairs ranked: {len(v2_ranked):,}
- V2 pairs with complete local-surface ensembles: {sum(row['surface_status'] == 'complete' for row in v2_ranked):,}
- V2 register-robust pairs: {sum(bool(row['register_robust']) for row in v2_ranked):,}
- Combined DRB1*15:01 pairs ranked: {len(combined_ranked):,}
- Combined pairs with complete local-surface ensembles: {sum(row['surface_status'] == 'complete' for row in combined_ranked):,}
- Development-control V3 capture at 3: {control_gate['v3_development_capture_at_3_count']}/{control_gate['development_panel_count']} panels
- Independent V3 validation systems: 0; definitive gate: not evaluable

The rankings are ready for descriptive review, but the future untouched six-system benchmark is still required before V3 can be called independently validated.
"""
    (out / "RESULTS_SUMMARY.md").write_text(results, encoding="utf-8")
    _write_union_csv(out / "SHA256SUMS.csv", _checksum_table(out))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--prewarm-cache", action="store_true")
    parser.add_argument("--finalize-existing", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    arguments = parser.parse_args()
    if arguments.prewarm_cache and arguments.finalize_existing:
        parser.error("choose only one of --prewarm-cache and --finalize-existing")
    if arguments.prewarm_cache:
        result = prewarm_surface_cache(workers=arguments.workers)
    elif arguments.finalize_existing:
        result = finalize_existing_package(arguments.out)
    else:
        result = run(out=arguments.out)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
