"""Pure contracts for the held-out human HLA-II positive-control benchmark.

The module deliberately has no dependency on discovery rankings. Control and
comparator selection accepts only sequence, HLA, provenance, and binding fields;
geometry is read only after a panel has been frozen.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import shlex
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np


ANCHOR_INDICES = np.asarray([0, 3, 5, 8])
EXPOSED_INDICES = np.asarray([1, 2, 4, 6, 7])
CLAIM_BOUNDARY = (
    "Computational pMHC geometry prioritization only; not evidence of presentation, "
    "TCR binding, activation, cross-reactivity, molecular mimicry, or MS mechanism."
)
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
FORMAL_CHARGE = {aa: 0.0 for aa in HYDROPATHY}
FORMAL_CHARGE.update({"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0})
H_BOND_DONOR = {aa: float(aa in "RKNQHSTWY") for aa in HYDROPATHY}
H_BOND_ACCEPTOR = {aa: float(aa in "DENQHSTY") for aa in HYDROPATHY}
AROMATIC = {aa: float(aa in "FWY") for aa in HYDROPATHY}


@dataclass(frozen=True)
class PmhcGeometry:
    ligand_id: str
    core_sequence: str
    groove_ca: np.ndarray
    core_ca: np.ndarray
    sidechain_vectors: np.ndarray


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "true", "yes"}


def validate_control_registry(
    systems: Sequence[Mapping[str, Any]],
    ligands: Sequence[Mapping[str, Any]],
    positive_pairs: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    system_ids = [str(row.get("system_id", "")) for row in systems]
    ligand_ids = [str(row.get("ligand_id", "")) for row in ligands]
    pair_ids = [str(row.get("pair_id", "")) for row in positive_pairs]
    for label, values in (("system", system_ids), ("ligand", ligand_ids), ("pair", pair_ids)):
        if not all(values) or len(values) != len(set(values)):
            raise ValueError(f"{label} IDs must be non-empty and unique")
    known_systems = set(system_ids)
    known_ligands = set(ligand_ids)
    for ligand in ligands:
        if str(ligand.get("system_id", "")) not in known_systems:
            raise ValueError(f"unknown system for ligand {ligand['ligand_id']}")
        sequence = str(ligand.get("sequence", ""))
        core = str(ligand.get("core", ""))
        if len(core) != 9 or core not in sequence:
            raise ValueError(f"exact nine-residue core is invalid for {ligand['ligand_id']}")
        if not all(str(ligand.get(field, "")) for field in (
            "mhc_alpha_allele", "mhc_beta_allele", "pdb_id"
        )):
            raise ValueError(f"incomplete pMHC identity for {ligand['ligand_id']}")
    strict_ids = {
        str(row["system_id"]) for row in systems
        if str(row.get("eligibility")) == "strict"
    }
    pairs_by_system = defaultdict(int)
    for pair in positive_pairs:
        system_id = str(pair.get("system_id", ""))
        if system_id not in strict_ids:
            raise ValueError(f"positive pair is not assigned to a strict system: {pair['pair_id']}")
        left = str(pair.get("left_ligand_id", ""))
        right = str(pair.get("right_ligand_id", ""))
        if left not in known_ligands or right not in known_ligands:
            raise ValueError(f"unknown ligand in positive pair {pair['pair_id']}")
        pairs_by_system[system_id] += 1
    if set(pairs_by_system) != strict_ids:
        raise ValueError("every strict system requires at least one positive pair")
    if any(int(row.get("independent_system_weight", 0)) != 1 for row in systems if str(row.get("eligibility")) == "strict"):
        raise ValueError("every strict biological system must receive exactly one vote")
    return {
        "strict_independent_system_count": len(strict_ids),
        "strict_positive_pair_count": len(positive_pairs),
        "strict_system_ids": sorted(strict_ids),
        "prospective_system_ids": sorted(
            str(row["system_id"]) for row in systems
            if str(row.get("eligibility")) == "prospective"
        ),
        "positive_pairs_by_system": dict(sorted(pairs_by_system.items())),
    }


def validate_comparator_registry(
    rows: Sequence[Mapping[str, Any]], *, expected_pair_ids: Sequence[str]
) -> Dict[str, Any]:
    expected = {str(pair_id) for pair_id in expected_pair_ids}
    observed = {str(row.get("positive_pair_id", "")) for row in rows}
    if not expected or observed != expected:
        raise ValueError(
            f"comparator pair IDs do not match the expected set: expected={sorted(expected)} "
            f"observed={sorted(observed)}"
        )
    grouped = defaultdict(list)
    for row in rows:
        pair_id = str(row.get("positive_pair_id", ""))
        arm = str(row.get("comparator_arm", ""))
        if arm not in {"microbial", "self"}:
            raise ValueError(f"invalid comparator arm for {pair_id}: {arm}")
        candidate_id = str(row.get("candidate_id", ""))
        sequence = str(row.get("sequence", ""))
        core = str(row.get("predicted_core", ""))
        if not candidate_id:
            raise ValueError(f"empty comparator candidate ID for {pair_id}/{arm}")
        starts = [
            index + 1 for index in range(max(0, len(sequence) - 8))
            if sequence[index:index + 9] == core
        ]
        if len(core) != 9 or len(starts) != 1:
            raise ValueError(f"exact unique nine-residue core is invalid for {candidate_id}")
        try:
            core_start = int(row.get("core_start_1_based", 0))
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid core start for {candidate_id}") from error
        if core_start != starts[0]:
            raise ValueError(f"core start does not match the peptide sequence for {candidate_id}")
        if str(row.get("register_resolution", "")) != "resolved_unique_fully_contained":
            raise ValueError(f"unresolved comparator register for {candidate_id}")
        if not str(row.get("seq_num", "")) or not str(row.get("raw_response_file", "")):
            raise ValueError(f"missing raw-response linkage for {candidate_id}")
        if str(row.get("negative_tier", "")) != "N3":
            raise ValueError(f"comparator {candidate_id} must be labeled N3")
        if str(row.get("recognition_status", "")) != "unknown_not_specificity_negative":
            raise ValueError(f"comparator {candidate_id} cannot be used as a specificity negative")
        if not _truth(row.get("selection_is_score_blind")):
            raise ValueError(f"comparator {candidate_id} was not selected score-blind")
        if any("geometry" in str(field).lower() for field in row):
            raise ValueError(f"geometry field leaked into comparator selection for {candidate_id}")
        grouped[(pair_id, arm)].append(candidate_id)
    for pair_id in sorted(expected):
        for arm in ("microbial", "self"):
            candidate_ids = grouped[(pair_id, arm)]
            if len(candidate_ids) != 5 or len(set(candidate_ids)) != 5:
                raise ValueError(f"exactly five unique {arm} comparators are required for {pair_id}")
    return {
        "comparison_pair_count": len(expected),
        "comparator_row_count": len(rows),
        "unique_comparator_count": len({str(row["candidate_id"]) for row in rows}),
    }


def _descriptor(aa: str) -> np.ndarray:
    if aa not in HYDROPATHY:
        raise ValueError(f"unsupported amino acid: {aa}")
    return np.asarray([
        (FORMAL_CHARGE[aa] + 1.0) / 2.0,
        (HYDROPATHY[aa] + 4.5) / 9.0,
        H_BOND_DONOR[aa],
        H_BOND_ACCEPTOR[aa],
        AROMATIC[aa],
    ])


def physicochemical_mismatch(left_core: str, right_core: str) -> float:
    if len(left_core) != 9 or len(right_core) != 9:
        raise ValueError("physicochemical comparison requires two exact nine-residue cores")
    values = [
        float(np.mean(np.abs(_descriptor(left_core[index]) - _descriptor(right_core[index]))))
        for index in EXPOSED_INDICES
    ]
    return float(np.mean(values))


def kabsch(source: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3:
        raise ValueError("Kabsch arrays must have equal N x 3 shapes")
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    u, _, vt = np.linalg.svd((source - source_center).T @ (target - target_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = u @ vt
    return rotation, target_center - source_center @ rotation


def _rmsd(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def pair_features(left: PmhcGeometry, right: PmhcGeometry) -> Dict[str, float]:
    if left.groove_ca.shape != right.groove_ca.shape:
        raise ValueError("equivalent MHC-II groove arrays are required")
    for value in (left.core_ca, right.core_ca, left.sidechain_vectors, right.sidechain_vectors):
        if np.asarray(value).shape != (9, 3):
            raise ValueError("core coordinate and vector arrays must be 9 x 3")
    rotation, translation = kabsch(right.groove_ca, left.groove_ca)
    fitted_ca = right.core_ca @ rotation + translation
    fitted_vectors = right.sidechain_vectors @ rotation
    return {
        "exposed_ca_rmsd_A": _rmsd(left.core_ca[EXPOSED_INDICES], fitted_ca[EXPOSED_INDICES]),
        "exposed_sidechain_vector_rmsd_A": _rmsd(
            left.sidechain_vectors[EXPOSED_INDICES], fitted_vectors[EXPOSED_INDICES]
        ),
        "tcr_face_physicochemical_mismatch": physicochemical_mismatch(
            left.core_sequence, right.core_sequence
        ),
        "anchor_ca_rmsd_A": _rmsd(left.core_ca[ANCHOR_INDICES], fitted_ca[ANCHOR_INDICES]),
    }


def rank_feature_percentiles(
    rows: Sequence[Mapping[str, Any]], features: Sequence[str]
) -> List[Dict[str, Any]]:
    output = [dict(row) for row in rows]
    count = len(output)
    for feature in features:
        ordered = sorted(range(count), key=lambda index: (float(output[index][feature]), index))
        position = 0
        while position < count:
            end = position + 1
            value = float(output[ordered[position]][feature])
            while end < count and float(output[ordered[end]][feature]) == value:
                end += 1
            average_rank_zero_based = (position + end - 1) / 2.0
            percentile = average_rank_zero_based / max(1, count - 1)
            for order_index in ordered[position:end]:
                output[order_index][f"{feature}_percentile"] = round(percentile, 12)
            position = end
    return output


def generate_weight_grid(features: Sequence[str]) -> List[Dict[str, float]]:
    features = tuple(features)
    if not features:
        raise ValueError("at least one feature is required")
    units = 4
    combinations: List[Tuple[int, ...]] = []

    def visit(prefix: Tuple[int, ...], remaining_features: int, remaining_units: int) -> None:
        if remaining_features == 1:
            combinations.append(prefix + (remaining_units,))
            return
        for value in range(remaining_units + 1):
            visit(prefix + (value,), remaining_features - 1, remaining_units - value)

    visit(tuple(), len(features), units)
    return [
        {feature: value / units for feature, value in zip(features, combination)}
        for combination in combinations
    ]


def _group_key(row: Mapping[str, Any]) -> Tuple[str, str, str, str]:
    return (
        str(row["system_id"]),
        str(row.get("positive_pair_id", "")),
        str(row["layer"]),
        str(row["panel_seed"]),
    )


def _evaluate_weights(
    rows: Sequence[Mapping[str, Any]], features: Sequence[str], weights: Mapping[str, float]
) -> List[Dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[_group_key(row)].append(dict(row))
    results = []
    for (system_id, positive_pair_id, layer, panel_seed), group in sorted(grouped.items()):
        ranked = rank_feature_percentiles(group, features)
        for row in ranked:
            row["composite_score"] = sum(
                float(weights[feature]) * float(row[f"{feature}_percentile"])
                for feature in features
            )
        ordered = sorted(ranked, key=lambda row: (row["composite_score"], str(row["pair_id"])))
        positives = [index + 1 for index, row in enumerate(ordered) if row["pair_role"] == "positive"]
        if not positives:
            continue
        results.append({
            "system_id": system_id,
            "positive_pair_id": positive_pair_id,
            "layer": layer,
            "panel_seed": panel_seed,
            "positive_rank": max(positives),
            "positive_pair_count": len(positives),
            "comparison_count": len(ordered),
        })
    return results


def select_weights(
    rows: Sequence[Mapping[str, Any]], features: Sequence[str]
) -> Dict[str, float]:
    best_weights = None
    best_objective = None
    for weights in generate_weight_grid(features):
        evaluations = _evaluate_weights(rows, features, weights)
        if not evaluations:
            objective = (math.inf, math.inf, math.inf, math.inf, math.inf, tuple(weights.values()))
        else:
            capture_count = sum(int(row["positive_rank"] <= 3) for row in evaluations)
            worst_rank = max(int(row["positive_rank"]) for row in evaluations)
            mean_rr = sum(1.0 / int(row["positive_rank"]) for row in evaluations) / len(evaluations)
            nonzero = sum(value > 0 for value in weights.values())
            exposed_weight = weights.get("exposed_ca_rmsd_A", 0.0)
            objective = (
                -capture_count, worst_rank, -mean_rr, nonzero, -exposed_weight,
                tuple(weights[feature] for feature in features),
            )
        if best_objective is None or objective < best_objective:
            best_objective = objective
            best_weights = dict(weights)
    if best_weights is None:
        raise ValueError("could not select a weight vector")
    return best_weights


def leave_one_system_out(
    rows: Sequence[Mapping[str, Any]], features: Sequence[str]
) -> List[Dict[str, Any]]:
    systems = sorted({str(row["system_id"]) for row in rows})
    output = []
    for held_out in systems:
        training = [row for row in rows if str(row["system_id"]) != held_out]
        testing = [row for row in rows if str(row["system_id"]) == held_out]
        training_ids = sorted({str(row["system_id"]) for row in training})
        weights = select_weights(training, features)
        for result in _evaluate_weights(testing, features, weights):
            output.append({
                "held_out_system_id": held_out,
                "training_system_ids": ";".join(training_ids),
                **result,
                **{f"weight_{feature}": weights[feature] for feature in features},
                "evaluation_status": "complete",
                "capture_at_3": int(result["positive_rank"]) <= 3,
            })
    return output


def select_score_blind_comparators(
    target: Mapping[str, Any], candidates: Sequence[Mapping[str, Any]], *, count: int, seed: int
) -> List[Dict[str, Any]]:
    if count < 1:
        raise ValueError("comparator count must be positive")
    target_sequence = str(target["sequence"])
    target_rank = float(target.get("binding_percentile", 0.0))
    selected = []
    for candidate in candidates:
        candidate_id = str(candidate["candidate_id"])
        sequence = str(candidate["sequence"])
        digest = hashlib.sha256(f"{seed}|{candidate_id}".encode("utf-8")).hexdigest()
        selected.append({
            "candidate_id": candidate_id,
            "sequence": sequence,
            "binding_percentile": candidate.get("binding_percentile", ""),
            "selection_length_difference": abs(len(sequence) - len(target_sequence)),
            "selection_binding_percentile_difference": abs(
                float(candidate.get("binding_percentile", 0.0)) - target_rank
            ),
            "selection_seeded_hash": digest,
            "selection_is_score_blind": True,
        })
    selected.sort(key=lambda row: (
        row["selection_length_difference"],
        row["selection_binding_percentile_difference"],
        row["selection_seeded_hash"],
        row["candidate_id"],
    ))
    if len(selected) < count:
        raise ValueError("insufficient exact-HLA comparators")
    return selected[:count]


def build_pdb_oracle_pairings(
    positive_pair: Mapping[str, Any], structural_ligands: Sequence[Mapping[str, Any]]
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Enumerate score-blind exact-HLA structural pairings for one positive."""
    ligand_ids = [str(row.get("ligand_id", "")) for row in structural_ligands]
    if not all(ligand_ids) or len(ligand_ids) != len(set(ligand_ids)):
        raise ValueError("structural ligand IDs must be non-empty and unique")
    left_hla = (
        str(positive_pair["left_mhc_alpha_allele"]),
        str(positive_pair["left_mhc_beta_allele"]),
    )
    right_hla = (
        str(positive_pair["right_mhc_alpha_allele"]),
        str(positive_pair["right_mhc_beta_allele"]),
    )
    by_hla = defaultdict(list)
    for row in structural_ligands:
        by_hla[(str(row["mhc_alpha_allele"]), str(row["mhc_beta_allele"]))].append(
            str(row["ligand_id"])
        )
    left_positive = str(positive_pair["left_ligand_id"])
    right_positive = str(positive_pair["right_ligand_id"])
    if left_positive not in by_hla[left_hla] or right_positive not in by_hla[right_hla]:
        raise ValueError(f"positive ligands are absent from the exact-HLA structural pool for {positive_pair['pair_id']}")
    if left_hla == right_hla:
        pool = sorted(by_hla[left_hla])
        candidates = [
            (pool[left], pool[right])
            for left in range(len(pool))
            for right in range(left + 1, len(pool))
        ]
        positive_key = tuple(sorted((left_positive, right_positive)))
        candidates = [pair for pair in candidates if tuple(sorted(pair)) != positive_key]
    else:
        candidates = [
            (left, right)
            for left in sorted(by_hla[left_hla])
            for right in sorted(by_hla[right_hla])
            if (left, right) != (left_positive, right_positive)
        ]
    rows = [{
        "system_id": positive_pair["system_id"],
        "positive_pair_id": positive_pair["pair_id"],
        "pair_id": f"{positive_pair['pair_id']}|pdb|positive",
        "pair_role": "positive",
        "left_ligand_id": left_positive,
        "right_ligand_id": right_positive,
        "selection_is_score_blind": True,
    }]
    rows.extend({
        "system_id": positive_pair["system_id"],
        "positive_pair_id": positive_pair["pair_id"],
        "pair_id": f"{positive_pair['pair_id']}|pdb|{left}|{right}",
        "pair_role": "pdb_exact_hla_decoy",
        "left_ligand_id": left,
        "right_ligand_id": right,
        "selection_is_score_blind": True,
    } for left, right in candidates)
    decoy_count = len(candidates)
    status = "complete" if decoy_count >= 5 else "not_evaluable_insufficient_exact_hla_decoys"
    return rows, {
        "positive_pair_id": positive_pair["pair_id"],
        "decoy_count": decoy_count,
        "comparison_count": len(rows),
        "evaluation_status": status,
    }


def build_af3_job_batches(
    ligands: Sequence[Mapping[str, Any]],
    hla_sequences: Mapping[Tuple[str, str], Mapping[str, str]],
    *,
    panel_seeds: Sequence[int],
    batch_size: int = 30,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[List[Dict[str, Any]]]]:
    if batch_size < 1 or batch_size > 30:
        raise ValueError("AlphaFold Server batches must contain 1 to 30 jobs")
    jobs = []
    manifest = []
    for ligand in sorted(ligands, key=lambda row: str(row["ligand_id"])):
        alpha = str(ligand["mhc_alpha_allele"])
        beta = str(ligand["mhc_beta_allele"])
        hla = hla_sequences.get((alpha, beta))
        if hla is None:
            raise ValueError(f"missing HLA-II sequences for {alpha}/{beta}")
        chains = (
            str(hla["mhc_alpha_sequence"]),
            str(hla["mhc_beta_sequence"]),
            str(ligand["sequence"]),
        )
        if not all(chains):
            raise ValueError(f"empty AlphaFold chain for {ligand['ligand_id']}")
        for seed in panel_seeds:
            token = "".join(character.lower() if character.isalnum() else "_" for character in str(ligand["ligand_id"]))
            token = "_".join(part for part in token.split("_") if part)
            name = f"ebvms_hla2_control_{token}_s{int(seed)}"
            jobs.append({
                "name": name,
                "modelSeeds": [int(seed)],
                "sequences": [
                    {"proteinChain": {"sequence": sequence, "count": 1}}
                    for sequence in chains
                ],
                "dialect": "alphafoldserver",
                "version": 1,
            })
            manifest.append({
                "job_name": name,
                "ligand_id": ligand["ligand_id"],
                "system_id": ligand.get("system_id", ""),
                "ligand_role": ligand.get("ligand_role", ""),
                "mhc_alpha_allele": alpha,
                "mhc_beta_allele": beta,
                "peptide_sequence": ligand["sequence"],
                "core_sequence": ligand.get(
                    "core_sequence", ligand.get("core", ligand.get("predicted_core", ""))
                ),
                "core_start_1_based": ligand.get("core_start_1_based", ""),
                "register_resolution": ligand.get("register_resolution", ""),
                "register_source": ligand.get("register_source", ""),
                "seq_num": ligand.get("seq_num", ""),
                "raw_response_file": ligand.get("raw_response_file", ""),
                "panel_seed": int(seed),
                "chain_roles": "mhc_alpha;mhc_beta;peptide",
                "status": "prepared_not_submitted",
                "claim_boundary": CLAIM_BOUNDARY,
            })
    if len({job["name"] for job in jobs}) != len(jobs):
        raise ValueError("AlphaFold job names must be unique")
    batches = [jobs[index:index + batch_size] for index in range(0, len(jobs), batch_size)]
    return jobs, manifest, batches


def validate_af3_job_package(
    batches: Sequence[Sequence[Mapping[str, Any]]],
    manifest: Sequence[Mapping[str, Any]],
    hla_sequences: Mapping[Tuple[str, str], Mapping[str, str]],
) -> Dict[str, Any]:
    if not batches or any(not batch or len(batch) > 30 for batch in batches):
        raise ValueError("AlphaFold batches must each contain 1 to 30 jobs")
    jobs = [dict(job) for batch in batches for job in batch]
    job_names = [str(job.get("name", "")) for job in jobs]
    manifest_by_name = {str(row.get("job_name", "")): row for row in manifest}
    if not all(job_names) or len(job_names) != len(set(job_names)):
        raise ValueError("AlphaFold job names must be non-empty and unique")
    if len(manifest_by_name) != len(manifest) or set(job_names) != set(manifest_by_name):
        raise ValueError("AlphaFold jobs and manifest rows do not match one-to-one")
    required_job_keys = {"name", "modelSeeds", "sequences", "dialect", "version"}
    for job in jobs:
        name = str(job["name"])
        row = manifest_by_name[name]
        if set(job) != required_job_keys or job["dialect"] != "alphafoldserver" or job["version"] != 1:
            raise ValueError(f"invalid AlphaFold Server schema for {name}")
        if str(row.get("status", "")) != "prepared_not_submitted":
            raise ValueError(f"invalid submission status for {name}")
        if str(row.get("chain_roles", "")) != "mhc_alpha;mhc_beta;peptide":
            raise ValueError(f"invalid chain roles for {name}")
        alpha = str(row["mhc_alpha_allele"])
        beta = str(row["mhc_beta_allele"])
        hla = hla_sequences.get((alpha, beta))
        if hla is None:
            raise ValueError(f"unknown exact HLA alpha/beta pair for {name}")
        sequences = job.get("sequences", [])
        if len(sequences) != 3 or any(
            set(item) != {"proteinChain"}
            or set(item["proteinChain"]) != {"sequence", "count"}
            or item["proteinChain"]["count"] != 1
            for item in sequences
        ):
            raise ValueError(f"invalid protein-chain schema for {name}")
        observed = [item["proteinChain"]["sequence"] for item in sequences]
        expected = [
            str(hla["mhc_alpha_sequence"]),
            str(hla["mhc_beta_sequence"]),
            str(row["peptide_sequence"]),
        ]
        if observed != expected:
            raise ValueError(f"chain order or exact sequence mismatch for {name}")
        if job["modelSeeds"] != [int(row["panel_seed"])]:
            raise ValueError(f"model seed mismatch for {name}")
        peptide = str(row["peptide_sequence"])
        core = str(row.get("core_sequence", ""))
        try:
            start = int(row.get("core_start_1_based", 0)) - 1
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid register start for {name}") from error
        if len(core) != 9 or start < 0 or peptide[start:start + 9] != core:
            raise ValueError(f"exact nine-residue register mismatch for {name}")
        if not str(row.get("register_resolution", "")) or not str(row.get("register_source", "")):
            raise ValueError(f"missing register provenance for {name}")
        if str(row.get("ligand_role", "")) == "N3_comparator" and (
            not str(row.get("seq_num", "")) or not str(row.get("raw_response_file", ""))
        ):
            raise ValueError(f"missing comparator raw-response linkage for {name}")
    return {
        "job_count": len(jobs),
        "manifest_row_count": len(manifest),
        "batch_sizes": [len(batch) for batch in batches],
        "unique_job_name_count": len(set(job_names)),
        "chain_order": ["mhc_alpha", "mhc_beta", "peptide"],
        "submission_state": "prepared_not_submitted",
        "jobs_sha256": hashlib.sha256(stable_json(jobs).encode("utf-8")).hexdigest(),
    }


def build_trust_gate(
    results: Sequence[Mapping[str, Any]], *, required_system_ids: Sequence[str]
) -> Dict[str, Any]:
    by_system = defaultdict(list)
    for row in results:
        by_system[str(row["system_id"])].append(row)
    system_rows = []
    for system_id in required_system_ids:
        rows = by_system.get(str(system_id), [])
        failed = any(
            str(row.get("evaluation_status")) == "complete"
            and int(row.get("positive_rank", 10**9)) > 3
            for row in rows
        )
        incomplete = not rows or any(str(row.get("evaluation_status")) != "complete" for row in rows)
        status = "fail" if failed else ("not_evaluable" if incomplete else "pass")
        system_rows.append({"system_id": str(system_id), "trust_status": status})
    if any(row["trust_status"] == "fail" for row in system_rows):
        overall = "fail"
    elif any(row["trust_status"] == "not_evaluable" for row in system_rows):
        overall = "not_evaluable"
    else:
        overall = "pass"
    return {
        "benchmark_version": "EBV_MS_HLA2_HELD_OUT_CONTROLS_V1",
        "overall_trust_status": overall,
        "required_system_count": len(required_system_ids),
        "passed_system_count": sum(row["trust_status"] == "pass" for row in system_rows),
        "failed_system_count": sum(row["trust_status"] == "fail" for row in system_rows),
        "not_evaluable_system_count": sum(row["trust_status"] == "not_evaluable" for row in system_rows),
        "system_statuses": system_rows,
        "discovery_reranking_allowed": overall == "pass",
        "incomplete_results_block_formal_pass": True,
        "claim_boundary": CLAIM_BOUNDARY,
    }


def parse_mmcif_atoms(path: Path) -> Dict[str, List[Dict[str, Any]]]:
    """Parse the protein atom loop in AF3 and wwPDB mmCIF files."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = next((index for index, line in enumerate(lines) if line.strip() == "_atom_site.group_PDB"), None)
    if start is None:
        raise ValueError(f"No atom-site loop in {path}")
    headers = []
    index = start
    while index < len(lines) and lines[index].strip().startswith("_atom_site."):
        headers.append(lines[index].strip().replace("_atom_site.", "", 1))
        index += 1
    lookup = {name: position for position, name in enumerate(headers)}
    required = {
        "group_PDB", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id",
        "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv",
    }
    if not required.issubset(lookup):
        raise ValueError(f"Missing atom-site fields in {path}: {sorted(required - set(lookup))}")
    chains = defaultdict(dict)
    chain_order = defaultdict(list)
    while index < len(lines):
        stripped = lines[index].strip()
        index += 1
        if not stripped or stripped.startswith("#"):
            if chain_order:
                break
            continue
        if not (stripped.startswith("ATOM ") or stripped.startswith("HETATM ")):
            if chain_order:
                break
            continue
        fields = shlex.split(stripped)
        if len(fields) < len(headers) or fields[lookup["group_PDB"]] != "ATOM":
            continue
        residue_name = fields[lookup["label_comp_id"]]
        if residue_name not in AA3:
            continue
        chain_id = fields[lookup["label_asym_id"]]
        residue_id = fields[lookup["label_seq_id"]]
        if residue_id not in chains[chain_id]:
            chains[chain_id][residue_id] = {"aa": AA3[residue_name], "atoms": [], "bfactors": []}
            chain_order[chain_id].append(residue_id)
        residue = chains[chain_id][residue_id]
        try:
            xyz = tuple(float(fields[lookup[name]]) for name in ("Cartn_x", "Cartn_y", "Cartn_z"))
            bfactor = float(fields[lookup["B_iso_or_equiv"]])
        except ValueError:
            continue
        residue["atoms"].append({
            "name": fields[lookup["label_atom_id"]],
            "element": fields[lookup["type_symbol"]],
            "xyz": xyz,
        })
        residue["bfactors"].append(bfactor)
    return {
        chain_id: [chains[chain_id][residue_id] for residue_id in chain_order[chain_id]]
        for chain_id in chain_order
    }


def residue_sequence(residues: Sequence[Mapping[str, Any]]) -> str:
    return "".join(str(residue["aa"]) for residue in residues)


def _atom_coordinate(residue: Mapping[str, Any], atom_name: str) -> np.ndarray:
    for atom in residue["atoms"]:
        if atom["name"] == atom_name:
            return np.asarray(atom["xyz"], dtype=float)
    raise ValueError(f"missing {atom_name} atom")


def _sidechain_vector(residue: Mapping[str, Any]) -> np.ndarray:
    ca = _atom_coordinate(residue, "CA")
    coordinates = [
        np.asarray(atom["xyz"], dtype=float)
        for atom in residue["atoms"]
        if atom["name"] not in BACKBONE_ATOMS and atom["element"] != "H"
    ]
    if not coordinates:
        return np.zeros(3)
    return np.mean(coordinates, axis=0) - ca


def _global_alignment_index_map(reference: str, observed: str) -> Dict[int, int]:
    """Map reference positions to observed positions through a global alignment."""
    m, n = len(reference), len(observed)
    scores = np.zeros((m + 1, n + 1), dtype=int)
    scores[:, 0] = np.arange(m + 1) * -2
    scores[0, :] = np.arange(n + 1) * -2
    trace = np.zeros((m + 1, n + 1), dtype=np.int8)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            options = (
                scores[i - 1, j - 1] + (2 if reference[i - 1] == observed[j - 1] else -1),
                scores[i - 1, j] - 2,
                scores[i, j - 1] - 2,
            )
            direction = int(np.argmax(options))
            trace[i, j] = direction
            scores[i, j] = options[direction]
    mapping = {}
    i, j = m, n
    while i or j:
        direction = int(trace[i, j]) if i and j else (1 if i else 2)
        if direction == 0:
            mapping[i - 1] = j - 1
            i -= 1
            j -= 1
        elif direction == 1:
            i -= 1
        else:
            j -= 1
    return mapping


def aligned_chain_ca(
    residues: Sequence[Mapping[str, Any]],
    reference_sequence: str,
    *,
    reference_start: int,
    count: int = 85,
) -> np.ndarray:
    observed = residue_sequence(residues)
    mapping = _global_alignment_index_map(reference_sequence, observed)
    reference_positions = range(reference_start, reference_start + count)
    missing = [position for position in reference_positions if position not in mapping]
    if missing:
        raise ValueError(f"missing aligned MHC-II groove reference positions: {missing}")
    return np.vstack([
        _atom_coordinate(residues[mapping[position]], "CA")
        for position in reference_positions
    ])


def geometry_from_mmcif(
    path: Path,
    *,
    ligand_id: str,
    peptide_sequence: str,
    core_sequence: str,
    mhc_alpha_chain: str,
    mhc_beta_chain: str,
    peptide_chain: str,
    mhc_alpha_reference_sequence: str = "",
    mhc_beta_reference_sequence: str = "",
    mhc_alpha_reference_start: int = 0,
    mhc_beta_reference_start: int = 0,
) -> PmhcGeometry:
    model = parse_mmcif_atoms(path)
    for chain in (mhc_alpha_chain, mhc_beta_chain, peptide_chain):
        if chain not in model:
            raise ValueError(f"missing curated chain {chain} in {path.name}")
    peptide_residues = model[peptide_chain]
    observed = residue_sequence(peptide_residues)
    peptide_start = observed.find(peptide_sequence)
    if peptide_start < 0:
        core_start_observed = observed.find(core_sequence)
        if core_start_observed < 0:
            raise ValueError(f"exact peptide/core not found in {path.name}")
        core_residues = peptide_residues[core_start_observed:core_start_observed + 9]
    else:
        core_start = peptide_sequence.index(core_sequence)
        start = peptide_start + core_start
        core_residues = peptide_residues[start:start + 9]
    if len(core_residues) != 9 or residue_sequence(core_residues) != core_sequence:
        raise ValueError(f"incomplete exact core in {path.name}")
    if mhc_alpha_reference_sequence and mhc_beta_reference_sequence:
        alpha_ca = aligned_chain_ca(
            model[mhc_alpha_chain], mhc_alpha_reference_sequence,
            reference_start=mhc_alpha_reference_start,
        )
        beta_ca = aligned_chain_ca(
            model[mhc_beta_chain], mhc_beta_reference_sequence,
            reference_start=mhc_beta_reference_start,
        )
        groove_ca = np.vstack([alpha_ca, beta_ca])
    else:
        alpha = model[mhc_alpha_chain][:85]
        beta = model[mhc_beta_chain][:85]
        if len(alpha) != 85 or len(beta) != 85:
            raise ValueError(f"MHC-II groove is incomplete in {path.name}")
        groove_ca = np.vstack([_atom_coordinate(residue, "CA") for residue in [*alpha, *beta]])
    core_ca = np.vstack([_atom_coordinate(residue, "CA") for residue in core_residues])
    vectors = np.vstack([_sidechain_vector(residue) for residue in core_residues])
    return PmhcGeometry(ligand_id, core_sequence, groove_ca, core_ca, vectors)


def stable_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def rows_sha256(rows: Iterable[Mapping[str, Any]]) -> str:
    normalized = "\n".join(stable_json(dict(row)) for row in rows)
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()
