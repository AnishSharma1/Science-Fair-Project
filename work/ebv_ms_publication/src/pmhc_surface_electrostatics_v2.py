"""Pure geometry, surface sampling, ensemble, and gate logic for electrostatics V2.

This module has no project-path constants and never reads discovery data.  The
command-line workflow supplies only frozen positive-control registries.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from pmhc_surface_electrostatics import OpenDXGrid, trilinear_sample


FORBIDDEN_PATH_TOKENS = (
    "candidate_evidence",
    "high_yield_candidate",
    "literature_grounded_hla2_rankings",
    "discovery",
)


@dataclass(frozen=True)
class GrooveFrame:
    origin: np.ndarray
    longitudinal: np.ndarray
    transverse: np.ndarray
    outward: np.ndarray


BACKBONE_ATOMS = ("N", "CA", "C", "O")


def standardize_pmhc_chains(
    model: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    alpha_chain: str,
    beta_chain: str,
    peptide_chain: str,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Keep only curated pMHC chains and rename them to generic A/B/C roles."""
    selected = (alpha_chain, beta_chain, peptide_chain)
    if len(set(selected)) != 3 or any(chain not in model for chain in selected):
        raise ValueError("curated alpha, beta, and peptide chains must be distinct and present")
    standardized = {
        role: [
            {
                **{key: value for key, value in residue.items() if key != "atoms"},
                "atoms": [{**atom, "xyz": tuple(float(x) for x in atom["xyz"])} for atom in residue["atoms"]],
            }
            for residue in model[source]
        ]
        for role, source in zip(("A", "B", "C"), selected)
    }
    excluded = sorted(chain for chain in model if chain not in selected)
    return standardized, {
        "source_alpha_chain": alpha_chain,
        "source_beta_chain": beta_chain,
        "source_peptide_chain": peptide_chain,
        "excluded_chain_ids": excluded,
        "tcr_or_other_protein_chains_removed": bool(excluded),
    }


def _sequence(residues: Sequence[Mapping[str, Any]]) -> str:
    return "".join(str(residue["aa"]) for residue in residues)


def _alignment_map(reference: str, moving: str) -> dict[int, int]:
    m, n = len(reference), len(moving)
    scores = np.zeros((m + 1, n + 1), dtype=int)
    scores[:, 0] = np.arange(m + 1) * -2
    scores[0, :] = np.arange(n + 1) * -2
    trace = np.zeros((m + 1, n + 1), dtype=np.int8)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            options = (
                scores[i - 1, j - 1] + (2 if reference[i - 1] == moving[j - 1] else -1),
                scores[i - 1, j] - 2,
                scores[i, j - 1] - 2,
            )
            direction = int(np.argmax(options))
            trace[i, j] = direction
            scores[i, j] = options[direction]
    result: dict[int, int] = {}
    i, j = m, n
    while i or j:
        direction = int(trace[i, j]) if i and j else (1 if i else 2)
        if direction == 0:
            result[i - 1] = j - 1
            i -= 1
            j -= 1
        elif direction == 1:
            i -= 1
        else:
            j -= 1
    return result


def _atom(residue: Mapping[str, Any], name: str) -> np.ndarray | None:
    found = next((item for item in residue["atoms"] if str(item["name"]) == name), None)
    return None if found is None else np.asarray(found["xyz"], dtype=float)


def _core_ca(model: Mapping[str, Sequence[Mapping[str, Any]]], start_1_based: int) -> np.ndarray:
    start = int(start_1_based) - 1
    peptide = model.get("C", ())
    if start < 0 or start + 9 > len(peptide):
        raise ValueError("reference P1-P9 core is not contained in peptide chain")
    coordinates = [_atom(residue, "CA") for residue in peptide[start : start + 9]]
    if any(point is None for point in coordinates):
        raise ValueError("reference P1-P9 core is missing a CA atom")
    return np.vstack(coordinates)


def align_local_groove(
    moving: Mapping[str, Sequence[Mapping[str, Any]]],
    reference: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    reference_core_start_1_based: int,
    cutoff_A: float = 12.0,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Fit sequence-equivalent HLA backbone atoms lying near the reference core."""
    reference_core = _core_ca(reference, reference_core_start_1_based)
    target_coordinates: list[np.ndarray] = []
    moving_coordinates: list[np.ndarray] = []
    selected_residue_count = 0
    for chain in ("A", "B"):
        if chain not in moving or chain not in reference:
            raise ValueError(f"missing HLA chain {chain}")
        mapping = _alignment_map(_sequence(reference[chain]), _sequence(moving[chain]))
        for reference_index, reference_residue in enumerate(reference[chain]):
            heavy = [
                np.asarray(atom["xyz"], dtype=float)
                for atom in reference_residue["atoms"]
                if str(atom.get("element", "")).upper() != "H"
            ]
            if not heavy or min(float(np.linalg.norm(point - core_point)) for point in heavy for core_point in reference_core) > cutoff_A:
                continue
            if reference_index not in mapping:
                continue
            moving_residue = moving[chain][mapping[reference_index]]
            matched_this_residue = 0
            for atom_name in BACKBONE_ATOMS:
                target = _atom(reference_residue, atom_name)
                source = _atom(moving_residue, atom_name)
                if target is not None and source is not None:
                    target_coordinates.append(target)
                    moving_coordinates.append(source)
                    matched_this_residue += 1
            selected_residue_count += int(matched_this_residue > 0)
    if len(target_coordinates) < 12:
        raise ValueError("fewer than 12 equivalent local groove backbone atoms were found")
    target = np.vstack(target_coordinates)
    source = np.vstack(moving_coordinates)
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    u, _, vt = np.linalg.svd((source - source_center).T @ (target - target_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    translation = target_center - source_center @ rotation
    fitted = source @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    aligned: dict[str, list[dict[str, Any]]] = {}
    for chain, residues in moving.items():
        aligned[chain] = []
        for residue in residues:
            copied = {key: value for key, value in residue.items() if key != "atoms"}
            copied["atoms"] = [
                {
                    **atom,
                    "xyz": tuple(float(value) for value in (np.asarray(atom["xyz"]) @ rotation + translation)),
                }
                for atom in residue["atoms"]
            ]
            aligned[chain].append(copied)
    return aligned, {
        "matched_groove_residue_count": selected_residue_count,
        "matched_backbone_atom_count": len(target_coordinates),
        "reference_core_distance_cutoff_A": float(cutoff_A),
        "fit_rmsd_A": rmsd,
    }


def build_apbs_surface_input(
    pqr_path: Path,
    potential_prefix: Path,
    accessibility_prefix: Path,
    *,
    dime: tuple[int, int, int],
    lengths_A: tuple[float, float, float],
    center_A: tuple[float, float, float],
    coarse_lengths_A: tuple[float, float, float] | None = None,
    coarse_center_A: tuple[float, float, float] | None = None,
    solute_dielectric: float = 4.0,
    solvent_dielectric: float = 78.5,
    temperature_K: float = 298.15,
    salt_M: float = 0.15,
    linear: bool = False,
    write_accessibility: bool = True,
) -> str:
    """Render a deterministic single-grid APBS calculation for V2."""
    equation = "lpbe" if linear else "npbe"
    accessibility = f"    write smol dx {accessibility_prefix}\n" if write_accessibility else ""
    coarse_lengths = coarse_lengths_A or lengths_A
    coarse_center = coarse_center_A or center_A
    return f"""read
    mol pqr {pqr_path}
end
elec
    mg-auto
    dime {dime[0]} {dime[1]} {dime[2]}
    cglen {coarse_lengths[0]:.6f} {coarse_lengths[1]:.6f} {coarse_lengths[2]:.6f}
    fglen {lengths_A[0]:.6f} {lengths_A[1]:.6f} {lengths_A[2]:.6f}
    cgcent {coarse_center[0]:.6f} {coarse_center[1]:.6f} {coarse_center[2]:.6f}
    fgcent {center_A[0]:.6f} {center_A[1]:.6f} {center_A[2]:.6f}
    mol 1
    {equation}
    bcfl sdh
    pdie {solute_dielectric:g}
    sdie {solvent_dielectric:g}
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp {temperature_K:g}
    ion charge 1 conc {salt_M:g} radius 2.0
    ion charge -1 conc {salt_M:g} radius 1.8
    calcenergy no
    calcforce no
    write pot dx {potential_prefix}
{accessibility}end
quit
"""


def _unit(vector: np.ndarray, label: str) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm <= 1e-10:
        raise ValueError(f"cannot define {label} axis from a zero vector")
    return vector / norm


def canonical_groove_frame(
    core_ca: np.ndarray,
    alpha_groove_ca: np.ndarray,
    beta_groove_ca: np.ndarray,
) -> GrooveFrame:
    """Create the preregistered P1-to-P9, beta-ward, outward groove frame."""
    core = np.asarray(core_ca, dtype=float)
    alpha = np.asarray(alpha_groove_ca, dtype=float)
    beta = np.asarray(beta_groove_ca, dtype=float)
    if core.shape != (9, 3):
        raise ValueError("core_ca must be a 9-by-3 array")
    if alpha.ndim != 2 or beta.ndim != 2 or alpha.shape[1:] != (3,) or beta.shape[1:] != (3,):
        raise ValueError("groove coordinates must be N-by-3 arrays")
    longitudinal = _unit(core[-1] - core[0], "longitudinal")
    groove_centroid = np.vstack((alpha, beta)).mean(axis=0)
    outward_raw = core.mean(axis=0) - groove_centroid
    outward_raw -= longitudinal * float(np.dot(outward_raw, longitudinal))
    outward = _unit(outward_raw, "outward")
    beta_ward = beta.mean(axis=0) - alpha.mean(axis=0)
    beta_ward -= longitudinal * float(np.dot(beta_ward, longitudinal))
    beta_ward -= outward * float(np.dot(beta_ward, outward))
    transverse = _unit(beta_ward, "transverse")
    return GrooveFrame(
        origin=core[0].copy(),
        longitudinal=longitudinal,
        transverse=transverse,
        outward=outward,
    )


def dense_lateral_grid(
    core_ca: np.ndarray,
    frame: GrooveFrame,
    *,
    longitudinal_margin_A: float = 4.0,
    transverse_half_width_A: float = 14.0,
    spacing_A: float = 0.75,
) -> tuple[np.ndarray, list[dict[str, float | int]]]:
    """Generate the deterministic two-dimensional lattice below the outer surface."""
    if spacing_A <= 0:
        raise ValueError("spacing_A must be positive")
    core = np.asarray(core_ca, dtype=float)
    length = float(np.dot(core[-1] - core[0], frame.longitudinal))
    longitudinal = np.arange(
        -longitudinal_margin_A,
        length + longitudinal_margin_A + spacing_A / 2.0,
        spacing_A,
    )
    transverse = np.arange(
        -transverse_half_width_A,
        transverse_half_width_A + spacing_A / 2.0,
        spacing_A,
    )
    points: list[np.ndarray] = []
    metadata: list[dict[str, float | int]] = []
    for i, along in enumerate(longitudinal):
        for j, across in enumerate(transverse):
            points.append(frame.origin + along * frame.longitudinal + across * frame.transverse)
            metadata.append({
                "longitudinal_index": i,
                "transverse_index": j,
                "longitudinal_A": float(along),
                "transverse_A": float(across),
            })
    return np.asarray(points, dtype=float), metadata


def _grid_spacing(grid: OpenDXGrid) -> float:
    return float(min(np.linalg.norm(np.asarray(delta, dtype=float)) for delta in grid.deltas))


def _accessibility_normal(grid: OpenDXGrid, point: np.ndarray) -> np.ndarray:
    step = max(_grid_spacing(grid) * 0.5, 0.05)
    gradient = []
    for axis in np.eye(3):
        pair = np.vstack((point + axis * step, point - axis * step))
        sampled = trilinear_sample(grid, pair)
        gradient.append(float(sampled[0] - sampled[1]) / (2.0 * step))
    return _unit(np.asarray(gradient), "surface normal")


def sample_outer_surface(
    accessibility: OpenDXGrid,
    potential: OpenDXGrid,
    lateral_points: np.ndarray,
    outward_axis: np.ndarray,
    *,
    search_min_A: float = -8.0,
    search_max_A: float = 24.0,
    search_step_A: float = 0.25,
    offset_A: float = 0.5,
) -> list[dict[str, Any]]:
    """Find the last solute-to-solvent crossing and sample just outside it."""
    outward = _unit(np.asarray(outward_axis, dtype=float), "outward")
    heights = np.arange(search_min_A, search_max_A + search_step_A / 2.0, search_step_A)
    bases = np.asarray(lateral_points, dtype=float)
    line = bases[:, None, :] + heights[None, :, None] * outward[None, None, :]
    try:
        values = trilinear_sample(accessibility, line.reshape(-1, 3)).reshape(len(bases), len(heights))
    except ValueError:
        return [{"covered": False, "missing_reason": "surface_search_outside_grid"} for _ in bases]
    transitions = (values[:, :-1] <= 0.5) & (values[:, 1:] > 0.5)
    covered = np.any(transitions, axis=1)
    reverse_index = np.argmax(transitions[:, ::-1], axis=1)
    indices = transitions.shape[1] - 1 - reverse_index
    row_indices = np.arange(len(bases))
    low = values[row_indices, indices]
    high = values[row_indices, indices + 1]
    delta = high - low
    fractions = np.full(len(delta), 0.5, dtype=float)
    np.divide(0.5 - low, delta, out=fractions, where=np.abs(delta) > 1e-12)
    crossing_heights = heights[indices] + fractions * search_step_A
    surface_points = bases + crossing_heights[:, None] * outward
    valid_indices = np.flatnonzero(covered)
    results: list[dict[str, Any]] = [
        {"covered": False, "missing_reason": "no_outer_surface_crossing"} for _ in bases
    ]
    if not len(valid_indices):
        return results
    valid_points = surface_points[valid_indices]
    step = max(_grid_spacing(accessibility) * 0.5, 0.05)
    gradient = np.empty((len(valid_points), 3), dtype=float)
    try:
        for axis_index, axis in enumerate(np.eye(3)):
            plus = trilinear_sample(accessibility, valid_points + axis * step)
            minus = trilinear_sample(accessibility, valid_points - axis * step)
            gradient[:, axis_index] = (plus - minus) / (2.0 * step)
        norms = np.linalg.norm(gradient, axis=1)
        normal_ok = norms > 1e-10
        normals = np.zeros_like(gradient)
        normals[normal_ok] = gradient[normal_ok] / norms[normal_ok, None]
        flip = (normals @ outward) < 0.0
        normals[flip] *= -1.0
        sample_points = valid_points + float(offset_A) * normals
        sampled_potentials = trilinear_sample(potential, sample_points)
    except ValueError:
        return [{"covered": False, "missing_reason": "normal_or_offset_outside_grid"} for _ in bases]
    for local_index, global_index in enumerate(valid_indices):
        if not normal_ok[local_index]:
            results[int(global_index)] = {"covered": False, "missing_reason": "undefined_surface_normal"}
            continue
        results[int(global_index)] = {
            "covered": True,
            "surface_height_A": float(crossing_heights[global_index]),
            "surface_point": valid_points[local_index],
            "normal": normals[local_index],
            "sample_point": sample_points[local_index],
            "potential_kT_per_e": float(sampled_potentials[local_index]),
        }
    return results


def label_surface_regions(
    surface_points: np.ndarray,
    atoms_by_region: Mapping[str, np.ndarray],
) -> list[str]:
    """Assign peptide/alpha/beta labels by the nearest non-hydrogen atom."""
    ordered_regions = ("peptide", "hla_alpha", "hla_beta")
    missing = [region for region in ordered_regions if region not in atoms_by_region]
    if missing:
        raise ValueError(f"missing atom regions: {missing}")
    labels: list[str] = []
    for point in np.asarray(surface_points, dtype=float):
        distances = []
        for region in ordered_regions:
            atoms = np.asarray(atoms_by_region[region], dtype=float)
            if atoms.ndim != 2 or atoms.shape[1:] != (3,) or not len(atoms):
                raise ValueError(f"region {region} requires a nonempty N-by-3 atom array")
            distances.append(float(np.min(np.linalg.norm(atoms - point, axis=1))))
        labels.append(ordered_regions[int(np.argmin(distances))])
    return labels


def hierarchical_pair_summary(score_matrix: np.ndarray) -> dict[str, Any]:
    """Summarize a 5x5 ensemble without treating 25 cells as replicates."""
    matrix = np.asarray(score_matrix, dtype=float)
    if matrix.shape != (5, 5) or not np.all(np.isfinite(matrix)):
        raise ValueError("a complete finite 5-by-5 score matrix is required")
    left = np.median(matrix, axis=1)
    right = np.median(matrix, axis=0)
    left_summary = float(np.median(left))
    right_summary = float(np.median(right))
    return {
        "left_marginal_medians": left.tolist(),
        "right_marginal_medians": right.tolist(),
        "left_marginal_summary": left_summary,
        "right_marginal_summary": right_summary,
        "conservative_score": min(left_summary, right_summary),
        "model_combination_count": 25,
        "independent_replicate_count": 0,
    }


def build_control_gate(requirements: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Apply the fail/missing/supportive decision order to frozen requirements."""
    rows = [dict(row) for row in requirements]
    missing = [row for row in rows if str(row.get("status")) != "complete"]
    failures = []
    for row in rows:
        if str(row.get("status")) != "complete":
            continue
        if row.get("rank") is None or int(row["rank"]) > 3:
            failures.append({**row, "failure_reason": "rank_above_3"})
            continue
        if not bool(row.get("sensitivity_top3", False)):
            failures.append({**row, "failure_reason": "sensitivity_rank_class_not_robust"})
            continue
        fraction = row.get("resampling_top3_fraction")
        if str(row.get("layer", "")).startswith("af_") and (
            fraction is None or float(fraction) < 0.80
        ):
            failures.append({**row, "failure_reason": "resampling_stability_below_0_80"})
    if failures:
        status = "fail"
    elif missing or not rows:
        status = "not_evaluable"
    else:
        status = "supportive"
    return {
        "benchmark_version": "PMHC_SURFACE_ELECTROSTATICS_V2_CONTROLS",
        "status": status,
        "required_result_count": len(rows),
        "failed_result_count": len(failures),
        "missing_result_count": len(missing),
        "failures": failures,
        "missing": missing,
        "candidate_evaluation_allowed": status == "supportive",
        "electrostatics_retired_from_candidate_ranking": status == "fail",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
    }


def validate_control_only_paths(paths: Sequence[Path], project_root: Path) -> None:
    """Reject inputs whose resolved project-relative name identifies discovery data."""
    root = Path(project_root).resolve()
    for path in paths:
        resolved = Path(path).resolve()
        try:
            relative = resolved.relative_to(root)
        except ValueError as error:
            raise ValueError(f"control input is outside project root: {resolved}") from error
        lowered = str(relative).lower()
        if any(token in lowered for token in FORBIDDEN_PATH_TOKENS):
            raise ValueError(f"candidate/discovery input is forbidden: {relative}")
