"""Versioned contracts for the two-stage HLA-II benchmark v2.

The module is deliberately additive. It does not import, read, or write discovery
rankings, and it never converts pilot results into frozen weights.
"""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from typing import Any, Mapping, Sequence


PILOT_SEEDS = (271828, 314159)
TCR_FACING_POSITIONS = ("P2", "P3", "P5", "P7", "P8")
TCR_FACING_INDICES = (1, 2, 4, 6, 7)
ANCHOR_POSITIONS = ("P1", "P4", "P6", "P9")
STRUCTURAL_FEATURES = (
    "exposed_ca_rmsd_A",
    "exposed_sidechain_vector_rmsd_A",
    "anchor_ca_rmsd_A",
)
FULL_COMPOSITE_FEATURES = (
    "exposed_ca_rmsd_A",
    "exposed_sidechain_vector_rmsd_A",
    "tcr_face_physicochemical_mismatch",
    "anchor_ca_rmsd_A",
)
NONSTRUCTURAL_BASELINES = (
    "tcr_face_physicochemical_mismatch",
    "tcr_facing_sequence_identity",
    "full_core_sequence_identity",
    "tcr_facing_blosum62_similarity",
    "full_core_blosum62_similarity",
)
DIAGNOSTIC_ONLY_FEATURES = (
    "binding_percentile_similarity",
    "peptide_length_register_agreement",
)
CLAIM_BOUNDARY_V2 = (
    "Descriptive computational pMHC prioritization only; not evidence of presentation, "
    "TCR binding, activation, cross-reactivity, molecular mimicry, MS mechanism, "
    "probability, or false-discovery rate."
)


def _as_bool(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "true", "yes"}


_BLOSUM62_ALPHABET = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_BLOSUM62 = {
    (left, right): value
    for left, row in zip(_BLOSUM62_ALPHABET, _BLOSUM62_ROWS.strip().splitlines())
    for right, value in zip(_BLOSUM62_ALPHABET, (int(item) for item in row.split()))
}


def _stable_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_sha256(value: str, label: str) -> None:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")


def build_protocol_lock(
    *,
    strict_system_ids: Sequence[str],
    positive_pair_ids: Sequence[str],
    registry_sha256: str,
    comparator_sha256: str,
    oracle_pairings_sha256: str = "0" * 64,
    software_versions: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Return the immutable, pre-geometry pilot protocol record."""
    systems = sorted(str(value) for value in strict_system_ids)
    pairs = sorted(str(value) for value in positive_pair_ids)
    if len(systems) != 3 or len(set(systems)) != 3:
        raise ValueError("the attribution pilot requires exactly three independent systems")
    if not pairs or len(pairs) != len(set(pairs)):
        raise ValueError("positive pair IDs must be non-empty and unique")
    _validate_sha256(registry_sha256, "registry_sha256")
    _validate_sha256(comparator_sha256, "comparator_sha256")
    _validate_sha256(oracle_pairings_sha256, "oracle_pairings_sha256")
    payload: dict[str, Any] = {
        "benchmark_version": "EBV_MS_HLA2_BENCHMARK_V2_PILOT",
        "benchmark_stage": "three_system_attribution_pilot",
        "status": "prepared_not_submitted",
        "claim_boundary": CLAIM_BOUNDARY_V2,
        "strict_system_ids": systems,
        "positive_pair_ids": pairs,
        "independent_system_count": 3,
        "one_vote_per_system": True,
        "panel_seeds": list(PILOT_SEEDS),
        "reuses_v1_jobs": False,
        "comparators_per_arm": 5,
        "pair_decoys_per_panel": 25,
        "batch_job_limit": 30,
        "tcr_facing_positions": list(TCR_FACING_POSITIONS),
        "anchor_positions": list(ANCHOR_POSITIONS),
        "full_composite_features": list(FULL_COMPOSITE_FEATURES),
        "structural_features": list(STRUCTURAL_FEATURES),
        "nonstructural_baselines": list(NONSTRUCTURAL_BASELINES),
        "diagnostic_only_features": list(DIAGNOSTIC_ONLY_FEATURES),
        "binding_percentile_enters_composite": False,
        "discovery_files_accessible_during_calibration": False,
        "pilot_can_freeze_weights": False,
        "pilot_can_unlock_discovery": False,
        "registry_sha256": registry_sha256,
        "comparator_sha256": comparator_sha256,
        "oracle_pairings_sha256": oracle_pairings_sha256,
        "software_versions": dict(software_versions or {}),
        "tie_break": "lexical_pair_id_after_fixed_numeric_tolerance",
        "random_ranking_rule": "deterministic_sha256_panel_identity_permutation",
        "missing_results_rule": "not_evaluable_never_pass",
        "pdb_oracle_rule": "mandatory_only_if_five_exact_hla_decoys_frozen_before_scoring",
        "pilot_gate_logic": (
            "all_panels_complete_and_capture_at_3;majority_systems_better_than_training_selected_"
            "nonstructural_baseline;none_worse;credited_improvements_require_structural_weight_"
            "at_least_0.25_and_are_removed_by_structural_ablation"
        ),
        "definitive_gate_logic": (
            "minimum_6_strict_systems;target_8;minimum_2_hla2_families;all_panels_complete_"
            "and_capture_at_3;majority_better;none_worse;strictly_better_system_weighted_mrr"
        ),
    }
    payload["protocol_sha256"] = _stable_sha256(payload)
    return payload


def blosum62_similarity(
    left_core: str, right_core: str, *, positions: Sequence[int] | None = None
) -> float:
    if len(left_core) != 9 or len(right_core) != 9:
        raise ValueError("BLOSUM62 comparison requires two exact nine-residue cores")
    indices = tuple(range(9)) if positions is None else tuple(int(value) for value in positions)
    if not indices or any(index < 0 or index >= 9 for index in indices):
        raise ValueError("BLOSUM62 positions are out of range")
    try:
        scores = [_BLOSUM62[(left_core[index], right_core[index])] for index in indices]
        self_scores = [
            max(
                _BLOSUM62[(left_core[index], left_core[index])],
                _BLOSUM62[(right_core[index], right_core[index])],
            )
            for index in indices
        ]
    except KeyError as error:
        raise ValueError(f"unsupported amino acid in BLOSUM62 comparison: {error.args[0]}") from error
    return sum(scores) / sum(self_scores)


def aggregate_system_results(
    panel_rows: Sequence[Mapping[str, Any]],
    *,
    required_seeds: Sequence[int] = PILOT_SEEDS,
    required_pairs_by_system: Mapping[str, Sequence[str]] | None = None,
) -> list[dict[str, Any]]:
    """Collapse panel/seed outcomes to one conservative vote per TCR system."""
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in panel_rows:
        grouped[str(row["system_id"])].append(row)
    output = []
    expected_seeds = {int(seed) for seed in required_seeds}
    for system_id, rows in sorted(grouped.items()):
        observed_pairs = {str(row["positive_pair_id"]) for row in rows}
        required_pairs = set(
            str(value) for value in (required_pairs_by_system or {}).get(system_id, observed_pairs)
        )
        complete = bool(required_pairs)
        for pair_id in required_pairs:
            pair_rows = [row for row in rows if str(row["positive_pair_id"]) == pair_id]
            observed_seeds = {
                int(row["panel_seed"])
                for row in pair_rows
                if str(row.get("evaluation_status", "")) == "complete"
            }
            complete = complete and observed_seeds == expected_seeds
        if not complete:
            output.append({
                "system_id": system_id,
                "evaluation_status": "not_evaluable_missing_required_panels",
                "independent_system_vote": 1,
            })
            continue
        required_rows = [
            row for row in rows
            if str(row["positive_pair_id"]) in required_pairs
            and int(row["panel_seed"]) in expected_seeds
        ]
        composite_score = max(int(row["composite_rank"]) for row in required_rows)
        baseline_score = max(int(row["baseline_rank"]) for row in required_rows)
        ablated_score = max(int(row.get("ablated_rank", row["baseline_rank"])) for row in required_rows)
        improved_rows = [
            row for row in required_rows
            if int(row["composite_rank"]) < int(row["baseline_rank"])
        ]
        structural_weights = [float(row.get("structural_weight", 0.0)) for row in improved_rows]
        output.append({
            "system_id": system_id,
            "evaluation_status": "complete",
            "required_positive_pair_count": len(required_pairs),
            "required_panel_count": len(required_rows),
            "system_score": composite_score,
            "baseline_system_score": baseline_score,
            "ablated_system_score": ablated_score,
            "capture_at_3": composite_score <= 3,
            "system_outcome_vs_baseline": (
                "better" if composite_score < baseline_score
                else "worse" if composite_score > baseline_score
                else "equal"
            ),
            "minimum_structural_weight_on_improvements": min(structural_weights, default=0.0),
            "structural_ablation_removes_improvement": all(
                int(row["composite_rank"]) >= int(row["baseline_rank"])
                or int(row.get("ablated_rank", row["baseline_rank"])) >= int(row["baseline_rank"])
                for row in required_rows
            ),
            "independent_system_vote": 1,
        })
    return output


def _gate_counts(system_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    better = [row for row in system_rows if int(row["system_score"]) < int(row["baseline_system_score"])]
    worse = [row for row in system_rows if int(row["system_score"]) > int(row["baseline_system_score"])]
    equal = [row for row in system_rows if int(row["system_score"]) == int(row["baseline_system_score"])]
    structural_support = all(
        float(row.get("minimum_structural_weight_on_improvements", 0.0)) >= 0.25
        and bool(row.get("structural_ablation_removes_improvement"))
        for row in better
    )
    return {
        "better_system_count": len(better),
        "equal_system_count": len(equal),
        "worse_system_count": len(worse),
        "majority_better": len(better) > len(system_rows) / 2,
        "no_system_worse": not worse,
        "all_systems_capture_at_3": all(int(row["system_score"]) <= 3 for row in system_rows),
        "structural_support_on_every_improvement": bool(better) and structural_support,
        "composite_system_weighted_mrr": sum(1 / int(row["system_score"]) for row in system_rows) / len(system_rows),
        "baseline_system_weighted_mrr": sum(
            1 / int(row["baseline_system_score"]) for row in system_rows
        ) / len(system_rows),
    }


def build_pilot_attribution_gate(
    system_rows: Sequence[Mapping[str, Any]], *, required_system_ids: Sequence[str],
    oracle_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    required = sorted(str(value) for value in required_system_ids)
    by_id = {str(row["system_id"]): row for row in system_rows}
    complete = (
        len(required) == 3
        and set(by_id) == set(required)
        and all(str(by_id[system_id].get("evaluation_status")) == "complete" for system_id in required)
    )
    gate: dict[str, Any] = {
        "benchmark_version": "EBV_MS_HLA2_BENCHMARK_V2_PILOT",
        "required_system_ids": required,
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "cross_allele_consensus_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY_V2,
    }
    if not complete:
        gate.update({
            "pilot_attribution_status": "not_evaluable",
            "blocking_reason": "missing_or_incomplete_required_system",
        })
        return gate
    counts = _gate_counts([by_id[system_id] for system_id in required])
    supportive = (
        counts["all_systems_capture_at_3"]
        and counts["majority_better"]
        and counts["no_system_worse"]
        and counts["structural_support_on_every_improvement"]
    )
    gate.update(counts)
    mandatory_oracles = [row for row in oracle_rows if _as_bool(row.get("mandatory_if_scored"))]
    gate["mandatory_pdb_oracle_count"] = len(mandatory_oracles)
    gate["unavailable_pdb_oracle_count"] = sum(
        str(row.get("oracle_status")) == "not_evaluable_availability" for row in oracle_rows
    )
    if any(str(row.get("oracle_status")) == "fail" for row in mandatory_oracles):
        gate["pilot_attribution_status"] = "fail"
        gate["blocking_reason"] = "mandatory_pdb_oracle_rank_above_3"
    elif any(str(row.get("oracle_status")) != "pass" for row in mandatory_oracles):
        gate["pilot_attribution_status"] = "not_evaluable"
        gate["blocking_reason"] = "mandatory_pdb_oracle_pending_or_incomplete"
    else:
        gate["pilot_attribution_status"] = "supportive" if supportive else "fail"
    return gate


def build_oracle_availability(
    pair_rows: Sequence[Mapping[str, Any]], *, minimum_decoys: int = 5
) -> list[dict[str, Any]]:
    if minimum_decoys < 1:
        raise ValueError("minimum_decoys must be positive")
    output = []
    for row in pair_rows:
        count = int(row.get("eligible_decoy_count", 0))
        result = {
            **dict(row),
            "minimum_required_decoys": minimum_decoys,
            "mandatory_if_scored": count >= minimum_decoys,
        }
        if count < minimum_decoys:
            result["oracle_status"] = "not_evaluable_availability"
        elif row.get("positive_rank") in (None, ""):
            result["oracle_status"] = "required_pending_results"
        else:
            result["oracle_status"] = "pass" if int(row["positive_rank"]) <= 3 else "fail"
        output.append(result)
    return output


def validate_specificity_registry(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    required_fields = (
        "negative_id", "peptide", "exact_hla", "assay", "tested_condition", "outcome",
        "source_location",
    )
    admitted = []
    excluded_n3 = 0
    for row in rows:
        tier = str(row.get("negative_tier", ""))
        if tier == "N3":
            excluded_n3 += 1
            continue
        if tier not in {"N1", "N2"}:
            continue
        missing = [field for field in required_fields if not str(row.get(field, "")).strip()]
        if missing:
            raise ValueError(
                f"specificity negative {row.get('negative_id', '<unknown>')} is missing {missing}"
            )
        admitted.append(dict(row))
    return {
        "specificity_status": "prepared" if admitted else "not_evaluable_no_verified_negatives",
        "admitted_negative_count": len(admitted),
        "admitted_n1_count": sum(str(row["negative_tier"]) == "N1" for row in admitted),
        "admitted_n2_count": sum(str(row["negative_tier"]) == "N2" for row in admitted),
        "excluded_n3_count": excluded_n3,
        "ranking_gate_independent": True,
        "specificity_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY_V2,
    }


def build_definitive_ranking_gate(
    system_rows: Sequence[Mapping[str, Any]],
    system_registry: Sequence[Mapping[str, Any]],
    *,
    minimum_systems: int = 6,
    oracle_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    strict = [row for row in system_registry if str(row.get("eligibility")) == "strict"]
    strict_ids = sorted(str(row["system_id"]) for row in strict)
    families = sorted({str(row.get("hla_family", "")) for row in strict if row.get("hla_family")})
    base: dict[str, Any] = {
        "benchmark_version": "EBV_MS_HLA2_BENCHMARK_V2_DEFINITIVE",
        "minimum_independent_systems": minimum_systems,
        "target_independent_systems": 8,
        "strict_independent_system_count": len(strict_ids),
        "hla_families": families,
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "cross_allele_consensus_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY_V2,
    }
    if len(strict_ids) < minimum_systems:
        base["definitive_status"] = "blocked_registry_size"
        return base
    if len(families) < 2:
        base["definitive_status"] = "blocked_hla_family_diversity"
        return base
    if any(not _as_bool(row.get("distinct_biological_sources_verified")) for row in strict):
        base["definitive_status"] = "blocked_distinct_source_requirement"
        return base
    by_id = {str(row["system_id"]): row for row in system_rows}
    if set(strict_ids) != set(by_id) or any(
        str(by_id[system_id].get("evaluation_status")) != "complete" for system_id in strict_ids
    ):
        base["definitive_status"] = "not_evaluable"
        return base
    selected_rows = [by_id[system_id] for system_id in strict_ids]
    counts = _gate_counts(selected_rows)
    passed = (
        counts["all_systems_capture_at_3"]
        and counts["majority_better"]
        and counts["no_system_worse"]
        and counts["structural_support_on_every_improvement"]
        and counts["composite_system_weighted_mrr"] > counts["baseline_system_weighted_mrr"]
    )
    mandatory_oracles = [row for row in oracle_rows if _as_bool(row.get("mandatory_if_scored"))]
    oracle_failed = any(str(row.get("oracle_status")) == "fail" for row in mandatory_oracles)
    oracle_pending = any(str(row.get("oracle_status")) != "pass" for row in mandatory_oracles)
    base.update(counts)
    base["mandatory_pdb_oracle_count"] = len(mandatory_oracles)
    if oracle_failed:
        base["definitive_status"] = "fail"
    elif oracle_pending:
        base["definitive_status"] = "not_evaluable"
    else:
        base["definitive_status"] = "pass" if passed else "fail"
    final_pass = base["definitive_status"] == "pass"
    base["weights_frozen"] = final_pass
    base["discovery_unlock_allowed"] = final_pass
    return base
