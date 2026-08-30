"""Pure electrostatic comparison logic for the additive pMHC surface pilot.

The functions in this module do not invoke PDB2PQR or APBS and do not read
discovery labels. They define the frozen numerical behavior used by the
prepare, calculate, and analyze workflow.
"""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


CANDIDATE_EXPOSED_POSITIONS = (2, 3, 5, 7, 8)
DIELECTRIC_VALUES = (2.0, 4.0, 8.0)

AA1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
VDW_RADII_A = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


@dataclass(frozen=True)
class APBSParameters:
    solute_dielectric: float = 2.0
    solvent_dielectric: float = 78.5
    temperature_K: float = 298.15
    positive_ion_concentration_M: float = 0.15
    negative_ion_concentration_M: float = 0.15
    positive_ion_radius_A: float = 2.0
    negative_ion_radius_A: float = 1.8
    solvent_radius_A: float = 1.4


@dataclass(frozen=True)
class GridSpec:
    dime: tuple[int, int, int]
    coarse_length_A: tuple[float, float, float]
    coarse_center_A: tuple[float, float, float]
    fine_length_A: tuple[float, float, float]
    fine_center_A: tuple[float, float, float]
    coarse_min_A: tuple[float, float, float]
    coarse_max_A: tuple[float, float, float]
    fine_spacing_A: tuple[float, float, float]


@dataclass(frozen=True)
class OpenDXGrid:
    origin: np.ndarray
    deltas: np.ndarray
    values: np.ndarray


def _groove_ca(model: Mapping[str, Sequence[Mapping[str, Any]]], residue_count: int) -> np.ndarray:
    coordinates = []
    for chain in ("A", "B"):
        residues = model.get(chain, ())
        if len(residues) < residue_count:
            raise ValueError(f"chain {chain} has fewer than {residue_count} groove residues")
        for residue in residues[:residue_count]:
            atom = next((item for item in residue["atoms"] if item["name"] == "CA"), None)
            if atom is None:
                raise ValueError(f"chain {chain} groove residue is missing CA")
            coordinates.append(atom["xyz"])
    return np.asarray(coordinates, dtype=float)


def align_model_to_reference(
    model: Mapping[str, Sequence[Mapping[str, Any]]],
    reference: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    groove_residue_count: int = 85,
) -> tuple[dict[str, list[dict[str, Any]]], float]:
    """Align a three-chain pMHC model to a reference HLA groove."""
    moving = _groove_ca(model, groove_residue_count)
    target = _groove_ca(reference, groove_residue_count)
    moving_center = moving.mean(axis=0)
    target_center = target.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - moving_center).T @ (target - target_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    translation = target_center - moving_center @ rotation
    fitted = moving @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))
    aligned: dict[str, list[dict[str, Any]]] = {}
    for chain in ("A", "B", "C"):
        if chain not in model:
            continue
        aligned[chain] = []
        for residue in model[chain]:
            copied = {key: value for key, value in residue.items() if key != "atoms"}
            copied["atoms"] = [
                {
                    **atom,
                    "xyz": tuple(float(value) for value in (np.asarray(atom["xyz"], dtype=float) @ rotation + translation)),
                }
                for atom in residue["atoms"]
            ]
            aligned[chain].append(copied)
    return aligned, rmsd


def _fibonacci_sphere(count: int) -> np.ndarray:
    if count < 1:
        raise ValueError("sphere sample count must be positive")
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    points = []
    for index in range(count):
        y = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        angle = golden_angle * index
        points.append((math.cos(angle) * radius, y, math.sin(angle) * radius))
    return np.asarray(points, dtype=float)


def surface_patch_points(
    model: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    core_start_1_based: int,
    positions: Sequence[int] = CANDIDATE_EXPOSED_POSITIONS,
    samples_per_atom: int = 48,
    solvent_radius_A: float = 1.4,
) -> np.ndarray:
    """Return solvent-accessible points around declared peptide-facing side chains."""
    peptide = model.get("C", ())
    start = int(core_start_1_based) - 1
    if start < 0 or start + 9 > len(peptide):
        raise ValueError("declared P1-P9 core is not fully contained in peptide chain C")
    all_atoms: list[tuple[np.ndarray, str, str, int]] = []
    for chain in ("A", "B", "C"):
        for residue_index, residue in enumerate(model.get(chain, ())):
            for atom in residue["atoms"]:
                if str(atom.get("element", "")).upper() == "H":
                    continue
                all_atoms.append(
                    (np.asarray(atom["xyz"], dtype=float), str(atom.get("element", "C")).upper(), chain, residue_index)
                )
    coordinates = np.vstack([item[0] for item in all_atoms])
    radii = np.asarray([VDW_RADII_A.get(item[1], 1.70) + solvent_radius_A for item in all_atoms])
    sphere = _fibonacci_sphere(samples_per_atom)
    patch: list[np.ndarray] = []
    requested = {start + int(position) - 1 for position in positions}
    for residue_index in sorted(requested):
        residue = peptide[residue_index]
        source_atoms = [
            atom for atom in residue["atoms"]
            if str(atom.get("element", "")).upper() != "H" and str(atom["name"]) not in BACKBONE_ATOMS
        ]
        if not source_atoms:
            source_atoms = [atom for atom in residue["atoms"] if str(atom["name"]) == "CA"]
        for atom in source_atoms:
            center = np.asarray(atom["xyz"], dtype=float)
            radius = VDW_RADII_A.get(str(atom.get("element", "C")).upper(), 1.70) + solvent_radius_A
            candidates = center + sphere * radius
            for point in candidates:
                distances = np.linalg.norm(coordinates - point, axis=1)
                source_matches = np.linalg.norm(coordinates - center, axis=1) < 1e-8
                occluded = np.any((distances < radii - 1e-6) & ~source_matches)
                if not occluded:
                    patch.append(point)
    if not patch:
        raise ValueError("no solvent-accessible peptide patch points were generated")
    array = np.asarray(patch, dtype=float)
    order = np.lexsort((array[:, 2], array[:, 1], array[:, 0]))
    return array[order]


def build_common_accessible_field_patch(
    core_ca_coordinates: np.ndarray,
    groove_ca_coordinates: np.ndarray,
    molecular_atoms: Sequence[tuple[np.ndarray, np.ndarray]],
    *,
    positions: Sequence[int] = CANDIDATE_EXPOSED_POSITIONS,
    lateral_offset_A: float = 1.5,
    minimum_height_A: float = 2.0,
    maximum_height_A: float = 20.0,
    height_step_A: float = 0.25,
    probe_radius_A: float = 1.4,
    minimum_clearance_A: float = 0.25,
) -> tuple[np.ndarray, list[dict[str, float | int | str]]]:
    """Build a fixed position-matched field shell outside every panel model.

    Five local points are anchored to each requested core position. Each is
    shifted along the HLA-to-peptide outward normal by the smallest frozen
    height that is solvent-accessible in every supplied molecular model.
    """
    core = np.asarray(core_ca_coordinates, dtype=float)
    groove = np.asarray(groove_ca_coordinates, dtype=float)
    if core.shape != (9, 3) or groove.ndim != 2 or groove.shape[1] != 3 or not len(groove):
        raise ValueError("core must be 9-by-3 and groove must be a nonempty N-by-3 array")
    if not molecular_atoms:
        raise ValueError("at least one molecular atom set is required")
    normalized_atoms: list[tuple[np.ndarray, np.ndarray]] = []
    for coordinates, radii in molecular_atoms:
        coordinates = np.asarray(coordinates, dtype=float)
        radii = np.asarray(radii, dtype=float).reshape(-1)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3 or len(coordinates) != len(radii):
            raise ValueError("each atom set requires matching N-by-3 coordinates and N radii")
        normalized_atoms.append((coordinates, radii))
    longitudinal = core[-1] - core[0]
    longitudinal /= np.linalg.norm(longitudinal)
    outward = core.mean(axis=0) - groove.mean(axis=0)
    outward -= longitudinal * float(np.dot(outward, longitudinal))
    if np.linalg.norm(outward) <= 1e-8:
        raise ValueError("cannot resolve an outward HLA-to-peptide normal")
    outward /= np.linalg.norm(outward)
    transverse = np.cross(outward, longitudinal)
    transverse /= np.linalg.norm(transverse)
    offsets = (
        ("center", np.zeros(3)),
        ("longitudinal_minus", -longitudinal * lateral_offset_A),
        ("longitudinal_plus", longitudinal * lateral_offset_A),
        ("transverse_minus", -transverse * lateral_offset_A),
        ("transverse_plus", transverse * lateral_offset_A),
    )
    heights = np.arange(minimum_height_A, maximum_height_A + height_step_A / 2.0, height_step_A)
    points: list[np.ndarray] = []
    metadata: list[dict[str, float | int | str]] = []
    for position in positions:
        if int(position) < 1 or int(position) > 9:
            raise ValueError("core positions must fall within P1-P9")
        for offset_label, lateral in offsets:
            base = core[int(position) - 1] + lateral
            selected: tuple[np.ndarray, float, float] | None = None
            for height in heights:
                point = base + outward * float(height)
                clearance = min(
                    float(np.min(np.linalg.norm(coordinates - point, axis=1) - (radii + probe_radius_A)))
                    for coordinates, radii in normalized_atoms
                )
                if clearance >= minimum_clearance_A - 1e-10:
                    selected = (point, float(height), clearance)
                    break
            if selected is None:
                raise ValueError(
                    f"no common solvent-accessible point found for P{position} {offset_label} "
                    f"within {maximum_height_A:.2f} A"
                )
            point, height, clearance = selected
            points.append(point)
            metadata.append(
                {
                    "core_position": int(position),
                    "local_offset": offset_label,
                    "height_A": height,
                    "minimum_clearance_A": clearance,
                    "outward_normal_x": float(outward[0]),
                    "outward_normal_y": float(outward[1]),
                    "outward_normal_z": float(outward[2]),
                }
            )
    return np.asarray(points, dtype=float), metadata


def write_model_pdb(
    model: Mapping[str, Sequence[Mapping[str, Any]]],
    path: Path,
    *,
    include_peptide: bool = True,
) -> None:
    """Write aligned AlphaFold coordinates as a deterministic standard PDB."""
    chains = ("A", "B", "C") if include_peptide else ("A", "B")
    lines: list[str] = []
    serial = 1
    for chain in chains:
        for residue_index, residue in enumerate(model.get(chain, ()), start=1):
            residue_name = AA1_TO_3.get(str(residue.get("aa", "")).upper())
            if residue_name is None:
                raise ValueError(f"unsupported amino acid {residue.get('aa')!r}")
            bfactors = list(residue.get("bfactors", ()))
            for atom_index, atom in enumerate(residue["atoms"]):
                x, y, z = (float(value) for value in atom["xyz"])
                atom_name = str(atom["name"])
                element = str(atom.get("element", atom_name[:1])).upper()
                bfactor = float(bfactors[atom_index]) if atom_index < len(bfactors) else 0.0
                lines.append(
                    f"ATOM  {serial:5d} {atom_name:>4s} {residue_name:>3s} {chain}{residue_index:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{bfactor:6.2f}          {element:>2s}"
                )
                serial += 1
        lines.append("TER")
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def candidate_exposed_histidines(core: str) -> tuple[int, ...]:
    """Return one-based candidate-exposed histidine positions in a P1-P9 core."""
    core = str(core).upper()
    if len(core) != 9:
        raise ValueError("an HLA-II core must contain exactly nine residues")
    return tuple(position for position in CANDIDATE_EXPOSED_POSITIONS if core[position - 1] == "H")


def _paired_values(left: Sequence[float], right: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    left_array = np.asarray(left, dtype=float).reshape(-1)
    right_array = np.asarray(right, dtype=float).reshape(-1)
    if left_array.shape != right_array.shape or not left_array.size:
        raise ValueError("potential arrays must be nonempty and have identical shapes")
    if not np.all(np.isfinite(left_array)) or not np.all(np.isfinite(right_array)):
        raise ValueError("potential arrays must contain only finite values")
    return left_array, right_array


def hodgkin_similarity(left: Sequence[float], right: Sequence[float]) -> float:
    """Return the Hodgkin electrostatic similarity index in [-1, 1]."""
    left_array, right_array = _paired_values(left, right)
    denominator = float(np.dot(left_array, left_array) + np.dot(right_array, right_array))
    if denominator <= 1e-30:
        return 1.0
    return float(2.0 * np.dot(left_array, right_array) / denominator)


def carbo_similarity(left: Sequence[float], right: Sequence[float]) -> float:
    """Return the Carbo/cosine similarity of two potential vectors."""
    left_array, right_array = _paired_values(left, right)
    left_norm = float(np.linalg.norm(left_array))
    right_norm = float(np.linalg.norm(right_array))
    if left_norm <= 1e-15 and right_norm <= 1e-15:
        return 1.0
    if left_norm <= 1e-15 or right_norm <= 1e-15:
        return 0.0
    return float(np.dot(left_array, right_array) / (left_norm * right_norm))


def sign_agreement_fraction(left: Sequence[float], right: Sequence[float]) -> float:
    left_array, right_array = _paired_values(left, right)
    return float(np.mean(np.sign(left_array) == np.sign(right_array)))


def potential_rmse(left: Sequence[float], right: Sequence[float]) -> float:
    left_array, right_array = _paired_values(left, right)
    return float(np.sqrt(np.mean((left_array - right_array) ** 2)))


def summarize_electrostatic_ensemble(rows: Sequence[Mapping[str, Any]]) -> dict[str, float | int]:
    """Summarize a complete five-by-five pair ensemble conservatively."""
    if len(rows) != 25:
        raise ValueError("electrostatic ensemble requires exactly 25 model combinations")
    output: dict[str, float | int] = {"model_combination_count": 25}
    for field in ("hodgkin_similarity", "carbo_similarity", "sign_agreement_fraction"):
        values = np.asarray([float(row[field]) for row in rows], dtype=float)
        output[f"{field}_min"] = float(values.min())
        output[f"{field}_q25"] = float(np.quantile(values, 0.25))
        output[f"{field}_median"] = float(np.quantile(values, 0.5))
        output[f"{field}_q75"] = float(np.quantile(values, 0.75))
        output[f"{field}_max"] = float(values.max())
        output[f"{field}_iqr"] = float(np.quantile(values, 0.75) - np.quantile(values, 0.25))
    values = np.asarray([float(row["potential_rmse"]) for row in rows], dtype=float)
    output.update(
        {
            "potential_rmse_min": float(values.min()),
            "potential_rmse_q25": float(np.quantile(values, 0.25)),
            "potential_rmse_median": float(np.quantile(values, 0.5)),
            "potential_rmse_q75": float(np.quantile(values, 0.75)),
            "potential_rmse_max": float(values.max()),
            "potential_rmse_iqr": float(np.quantile(values, 0.75) - np.quantile(values, 0.25)),
        }
    )
    return output


def _multigrid_dimension(length_A: float, maximum_spacing_A: float) -> int:
    intervals = max(32, int(math.ceil(float(length_A) / float(maximum_spacing_A))))
    intervals = int(math.ceil(intervals / 32.0) * 32)
    return intervals + 1


def build_shared_grid(
    molecular_bounds: Sequence[tuple[np.ndarray, np.ndarray]],
    patch_points: np.ndarray,
    *,
    padding_A: float = 12.0,
    maximum_spacing_A: float = 0.5,
) -> GridSpec:
    """Build one deterministic APBS coarse/fine grid shared by a whole panel."""
    if not molecular_bounds:
        raise ValueError("at least one molecular bound is required")
    patch = np.asarray(patch_points, dtype=float)
    if patch.ndim != 2 or patch.shape[1] != 3 or not len(patch):
        raise ValueError("patch points must be a nonempty N-by-3 array")
    mins = np.vstack([np.asarray(bounds[0], dtype=float) for bounds in molecular_bounds])
    maxs = np.vstack([np.asarray(bounds[1], dtype=float) for bounds in molecular_bounds])
    fine_min = patch.min(axis=0) - padding_A
    fine_max = patch.max(axis=0) + padding_A
    coarse_min = np.minimum(mins.min(axis=0) - padding_A, fine_min)
    coarse_max = np.maximum(maxs.max(axis=0) + padding_A, fine_max)
    fine_length = fine_max - fine_min
    dime = tuple(_multigrid_dimension(value, maximum_spacing_A) for value in fine_length)
    spacing = tuple(float(value / (dimension - 1)) for value, dimension in zip(fine_length, dime))
    return GridSpec(
        dime=dime,
        coarse_length_A=tuple(float(value) for value in coarse_max - coarse_min),
        coarse_center_A=tuple(float(value) for value in (coarse_max + coarse_min) / 2.0),
        fine_length_A=tuple(float(value) for value in fine_length),
        fine_center_A=tuple(float(value) for value in (fine_max + fine_min) / 2.0),
        coarse_min_A=tuple(float(value) for value in coarse_min),
        coarse_max_A=tuple(float(value) for value in coarse_max),
        fine_spacing_A=spacing,
    )


def _triplet(values: Sequence[float | int]) -> str:
    return " ".join(f"{float(value):.8g}" for value in values)


def build_apbs_input(
    pqr_path: Path,
    potential_prefix: Path,
    grid: GridSpec,
    parameters: APBSParameters,
) -> str:
    """Render the frozen linearized Poisson-Boltzmann APBS input."""
    return f"""read
    mol pqr {pqr_path}
end
elec name pmhc_potential
    mg-auto
    dime {' '.join(str(value) for value in grid.dime)}
    cglen {_triplet(grid.coarse_length_A)}
    fglen {_triplet(grid.fine_length_A)}
    cgcent {_triplet(grid.coarse_center_A)}
    fgcent {_triplet(grid.fine_center_A)}
    mol 1
    lpbe
    bcfl sdh
    pdie {parameters.solute_dielectric:.8g}
    sdie {parameters.solvent_dielectric:.8g}
    srfm smol
    chgm spl2
    sdens 10.0
    srad {parameters.solvent_radius_A:.8g}
    swin 0.3
    temp {parameters.temperature_K:.8g}
    ion charge 1 conc {parameters.positive_ion_concentration_M:.8g} radius {parameters.positive_ion_radius_A:.8g}
    ion charge -1 conc {parameters.negative_ion_concentration_M:.8g} radius {parameters.negative_ion_radius_A:.8g}
    calcenergy no
    calcforce no
    write pot dx {potential_prefix}
end
quit
"""


def parse_open_dx(path: Path) -> OpenDXGrid:
    """Parse an axis-aligned APBS OpenDX scalar grid."""
    lines = path.read_text(encoding="ascii", errors="strict").splitlines()
    counts: tuple[int, int, int] | None = None
    origin: np.ndarray | None = None
    deltas: list[np.ndarray] = []
    item_count: int | None = None
    data_start: int | None = None
    for index, line in enumerate(lines):
        fields = line.split()
        if line.startswith("object 1 class gridpositions counts"):
            counts = tuple(int(value) for value in fields[-3:])
        elif fields[:1] == ["origin"]:
            origin = np.asarray([float(value) for value in fields[1:4]], dtype=float)
        elif fields[:1] == ["delta"] and len(deltas) < 3:
            deltas.append(np.asarray([float(value) for value in fields[1:4]], dtype=float))
        elif "data follows" in line and "class array" in line:
            match = re.search(r"items\s+(\d+)\s+data follows", line)
            if not match:
                raise ValueError("OpenDX array header lacks an item count")
            item_count = int(match.group(1))
            data_start = index + 1
            break
    if counts is None or origin is None or len(deltas) != 3 or item_count is None or data_start is None:
        raise ValueError("incomplete OpenDX grid header")
    if item_count != int(np.prod(counts)):
        raise ValueError("OpenDX item count does not match grid dimensions")
    data: list[float] = []
    for line in lines[data_start:]:
        if line.startswith(("attribute", "object", "component")):
            break
        data.extend(float(value) for value in line.split())
        if len(data) >= item_count:
            break
    if len(data) != item_count:
        raise ValueError("OpenDX scalar array is incomplete")
    delta_array = np.asarray(deltas, dtype=float)
    off_axis = delta_array - np.diag(np.diag(delta_array))
    if np.max(np.abs(off_axis)) > 1e-12 or np.any(np.diag(delta_array) <= 0):
        raise ValueError("only positive axis-aligned OpenDX grids are supported")
    return OpenDXGrid(
        origin=origin,
        deltas=np.diag(delta_array),
        values=np.asarray(data, dtype=float).reshape(counts),
    )


def trilinear_sample(grid: OpenDXGrid, points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("sample points must be an N-by-3 array")
    fractional = (points - grid.origin) / grid.deltas
    maximum = np.asarray(grid.values.shape, dtype=float) - 1.0
    if np.any(fractional < -1e-8) or np.any(fractional > maximum + 1e-8):
        raise ValueError("sample point lies outside the OpenDX grid")
    fractional = np.clip(fractional, 0.0, maximum)
    lower = np.floor(fractional).astype(int)
    upper = np.minimum(lower + 1, np.asarray(grid.values.shape) - 1)
    weight = fractional - lower
    sampled = np.zeros(len(points), dtype=float)
    for xbit in (0, 1):
        for ybit in (0, 1):
            for zbit in (0, 1):
                indices = np.column_stack(
                    (
                        np.where(xbit, upper[:, 0], lower[:, 0]),
                        np.where(ybit, upper[:, 1], lower[:, 1]),
                        np.where(zbit, upper[:, 2], lower[:, 2]),
                    )
                )
                factors = (
                    np.where(xbit, weight[:, 0], 1.0 - weight[:, 0])
                    * np.where(ybit, weight[:, 1], 1.0 - weight[:, 1])
                    * np.where(zbit, weight[:, 2], 1.0 - weight[:, 2])
                )
                sampled += grid.values[indices[:, 0], indices[:, 1], indices[:, 2]] * factors
    return sampled


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _number(value: Any, default: float = math.inf) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _seeded_hash(seed: int, target_id: str, candidate_id: str) -> str:
    return hashlib.sha256(f"{seed}|{target_id}|{candidate_id}".encode("utf-8")).hexdigest()


def select_supported_comparators(
    target: Mapping[str, Any],
    candidates: Sequence[Mapping[str, Any]],
    *,
    excluded_candidate_ids: set[str],
    count: int = 5,
    seed: int = 271828,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Select exact-HLA, predictor-supported comparators without geometry access."""
    target_length = len(str(target["sequence"]))
    target_binding = np.mean(
        [
            _number(target.get("netmhciipan_el_percentile")),
            _number(target.get("mixmhc2pred_percentile")),
        ]
    )
    provenance: list[dict[str, Any]] = []
    eligible: list[dict[str, Any]] = []
    for source in candidates:
        row = dict(source)
        identifier = str(row.get("candidate_id", ""))
        reason = "eligible"
        if row.get("allele") != target.get("allele"):
            reason = "wrong_hla"
        elif identifier in excluded_candidate_ids or identifier == target.get("candidate_id"):
            reason = "excluded_target_or_control"
        elif not _truth(row.get("binding_consensus")):
            reason = "binding_not_supported"
        elif not _truth(row.get("register_consensus")) or not _truth(row.get("declared_core_match")):
            reason = "register_not_supported"
        elif int(_number(row.get("model_count"), -1)) != 5 or row.get("surface_status") != "complete":
            reason = "incomplete_model_ensemble"
        elif len(str(row.get("core", ""))) != 9 or str(row.get("core", "")) not in str(row.get("sequence", "")):
            reason = "invalid_core_or_sequence"
        binding = np.mean(
            [
                _number(row.get("netmhciipan_el_percentile")),
                _number(row.get("mixmhc2pred_percentile")),
            ]
        )
        annotated = {
            **row,
            "eligibility_reason": reason,
            "selection_length_difference": abs(len(str(row.get("sequence", ""))) - target_length),
            "selection_binding_percentile_difference": abs(float(binding) - float(target_binding)),
            "selection_seeded_hash": _seeded_hash(seed, str(target.get("candidate_id", "")), identifier),
            "selection_uses_geometry": False,
        }
        provenance.append(annotated)
        if reason == "eligible":
            eligible.append(annotated)
    eligible.sort(
        key=lambda row: (
            int(row["selection_length_difference"]),
            float(row["selection_binding_percentile_difference"]),
            str(row["selection_seeded_hash"]),
            str(row["candidate_id"]),
        )
    )
    selected: list[dict[str, Any]] = []
    seen_accession_core: set[tuple[str, str]] = set()
    for row in eligible:
        key = (str(row.get("accession") or row.get("protein", "")), str(row["core"]))
        if key in seen_accession_core:
            row["eligibility_reason"] = "duplicate_accession_core"
            continue
        seen_accession_core.add(key)
        if len(selected) < count:
            row["selection_order"] = len(selected) + 1
            selected.append(row)
        else:
            row["eligibility_reason"] = "eligible_not_selected"
    selected_ids = {str(row["candidate_id"]) for row in selected}
    for row in provenance:
        if str(row.get("candidate_id")) in selected_ids:
            row["selected"] = True
        else:
            row["selected"] = False
    if len(selected) != count:
        return [], provenance
    return selected, provenance


def rank_panel(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    if any(not math.isfinite(_number(row.get("hodgkin_similarity_q25"))) for row in rows):
        raise ValueError("a complete panel requires finite q25 Hodgkin similarities")
    ordered = sorted(
        (dict(row) for row in rows),
        key=lambda row: (-float(row["hodgkin_similarity_q25"]), str(row["pair_id"])),
    )
    for rank, row in enumerate(ordered, start=1):
        row["electrostatic_rank"] = rank
        row["electrostatic_percentile"] = 0.0 if len(ordered) == 1 else (rank - 1) / (len(ordered) - 1)
    return ordered


def assign_electrostatic_context(rank: int | None, register_qc: bool, model_qc: bool) -> str:
    if rank is None or not register_qc or not model_qc:
        return "not_evaluable"
    return "electrostatic_context_supportive" if int(rank) <= 3 else "electrostatic_context_not_supportive"


def dielectric_robustness(status_by_dielectric: Mapping[float, str]) -> bool:
    expected = {2.0, 4.0, 8.0}
    if set(float(value) for value in status_by_dielectric) != expected:
        return False
    statuses = set(status_by_dielectric.values())
    return len(statuses) == 1 and "not_evaluable" not in statuses
