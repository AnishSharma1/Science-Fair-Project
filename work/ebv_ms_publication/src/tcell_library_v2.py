"""Contracts for the additive EBV-MS T-cell library and V2 benchmark.

The functions in this module are deliberately pure: panel and control selection
cannot inspect AlphaFold-derived fields, and all expansion steps are
deterministic.  External IEDB retrieval and file writing live in the runner.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import math
import re
from typing import Any, Iterable, Sequence

import numpy as np


LIBRARY_VERSION = "EBV_MS_TCELL_LIBRARY_2026-08-22"
PANEL_VERSION = "EBV_MS_TCELL_V2_80_2026-08-22"
ALLELES = (
    "HLA-DRB1*15:01",
    "HLA-DRB1*13:03",
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
)
CALIBRATION_SEEDS = (104729, 104759)
EVIDENCE_TIERS = {
    "E1_exact_pmhc_positive",
    "E2_human_tcell_protein_or_region",
    "E3_supportive_tcell",
    "context_only",
}
ALLOWED_KINGDOMS = {"EBV", "human_self"}
SOURCE_CERTAINTY_ORDER = {
    "exact_primary_source": 0,
    "exact_iedb_positive": 1,
    "region_mapped_primary_source": 2,
    "canonical_tiling_for_protein_level_evidence": 3,
    "canonical_exploratory": 4,
}
ANCHOR_POSITIONS = (1, 4, 6, 9)
EXPOSED_POSITIONS = (2, 3, 5, 7, 8)


def _truth(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def validate_registry(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Validate evidence-system identity and strict T-cell denominator rules."""
    ids = [str(row.get("biological_system_id", "")).strip() for row in rows]
    if not all(ids):
        raise ValueError("every registry row requires biological_system_id")
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate biological_system_id")
    for row in rows:
        tier = str(row.get("evidence_tier", ""))
        if tier not in EVIDENCE_TIERS:
            raise ValueError(f"invalid evidence tier: {tier}")
        if not str(row.get("primary_source", "")).strip():
            raise ValueError(f"missing primary source for {row['biological_system_id']}")
        modality = str(row.get("receptor_modality", "")).lower()
        in_denominator = _truth(row.get("tcell_positive_denominator", False))
        if "antibody" in modality and in_denominator:
            raise ValueError("antibody-only evidence cannot enter the T-cell denominator")
        if in_denominator and (tier != "E1_exact_pmhc_positive" or "t_cell" not in modality):
            raise ValueError("strict denominator is restricted to E1 human T-cell systems")
    strict_ids = sorted(
        str(row["biological_system_id"])
        for row in rows
        if _truth(row.get("tcell_positive_denominator", False))
    )
    return {
        "system_count": len(rows),
        "strict_e1_count": len(strict_ids),
        "tcell_denominator_ids": strict_ids,
        "tier_counts": dict(sorted(Counter(str(row["evidence_tier"]) for row in rows).items())),
    }


def _candidate_key(row: dict[str, Any]) -> tuple[Any, ...]:
    """Score-blind candidate key; AlphaFold/geometry fields are never read."""
    certainty = SOURCE_CERTAINTY_ORDER.get(str(row.get("source_certainty", "")), 99)
    return (
        0 if _truth(row.get("required_for_confirmed_system", False)) else 1,
        int(row.get("evidence_priority", 99)),
        0 if _truth(row.get("native_hla_evidence", False)) else 1,
        certainty,
        str(row.get("accession", "")),
        int(row.get("start") or 10**9),
        int(row.get("end") or 10**9),
        str(row.get("sequence", "")),
        str(row.get("candidate_id", "")),
    )


def _freeze_side(
    rows: Sequence[dict[str, Any]], *, desired: int, min_proteins: int, max_per_protein: int
) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row["protein_symbol"])].append(dict(row))
    if len(groups) < min_proteins:
        raise ValueError(f"protein diversity shortfall: {len(groups)} < {min_proteins}")
    for protein in groups:
        groups[protein].sort(key=_candidate_key)

    selected: list[dict[str, Any]] = []
    selected_ids: set[str] = set()
    counts: Counter[str] = Counter()

    required = sorted(
        (row for row in rows if _truth(row.get("required_for_confirmed_system", False))),
        key=_candidate_key,
    )
    for row in required:
        candidate_id = str(row["candidate_id"])
        protein = str(row["protein_symbol"])
        if candidate_id not in selected_ids and counts[protein] < max_per_protein:
            selected.append(dict(row))
            selected_ids.add(candidate_id)
            counts[protein] += 1

    # Breadth first: select one peptide from every protein before any second.
    for protein in sorted(groups):
        if counts[protein]:
            continue
        row = next((r for r in groups[protein] if str(r["candidate_id"]) not in selected_ids), None)
        if row is not None:
            selected.append(dict(row))
            selected_ids.add(str(row["candidate_id"]))
            counts[protein] += 1

    remaining = sorted(
        (row for row in rows if str(row["candidate_id"]) not in selected_ids),
        key=_candidate_key,
    )
    for row in remaining:
        if len(selected) == desired:
            break
        protein = str(row["protein_symbol"])
        if counts[protein] >= max_per_protein:
            continue
        selected.append(dict(row))
        selected_ids.add(str(row["candidate_id"]))
        counts[protein] += 1
    if len(selected) != desired:
        raise ValueError(f"could not freeze {desired} peptides under diversity/cap rules")
    return selected


def freeze_panel(candidates: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Freeze the exact 40-EBV/40-self V2 panel without structural information."""
    ids = [str(row.get("candidate_id", "")) for row in candidates]
    if not all(ids) or len(ids) != len(set(ids)):
        raise ValueError("candidate IDs must be non-empty and unique")
    for row in candidates:
        if row.get("kingdom") not in ALLOWED_KINGDOMS:
            raise ValueError(f"invalid kingdom for {row.get('candidate_id')}")
        sequence = str(row.get("sequence", "")).upper()
        if not sequence or re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
            raise ValueError(f"invalid amino-acid sequence for {row.get('candidate_id')}")
        if len(sequence) < 11 and not _truth(row.get("natural_flanks_verified", False)):
            raise ValueError("peptides shorter than 11 aa require verified natural flanks")

    # Alternative lengths/coordinate records that share the same proposed
    # 9-mer are retained in the master library but only one can be a panel
    # primary.  The score-blind evidence key selects that primary.
    collapsed: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in candidates:
        core = str(row.get("proposed_core", ""))
        collapse_key = (
            str(row["kingdom"]),
            str(row.get("accession", "")),
            core if core else f"NO_CORE::{row['candidate_id']}",
        )
        old = collapsed.get(collapse_key)
        if old is None or _candidate_key(row) < _candidate_key(old):
            collapsed[collapse_key] = row
    panel_candidates = list(collapsed.values())

    ebv = _freeze_side(
        [row for row in panel_candidates if row["kingdom"] == "EBV"],
        desired=40,
        min_proteins=18,
        max_per_protein=3,
    )
    human = _freeze_side(
        [row for row in panel_candidates if row["kingdom"] == "human_self"],
        desired=40,
        min_proteins=11,
        max_per_protein=4,
    )
    panel = ebv + human
    for index, row in enumerate(panel, start=1):
        row["library_version"] = LIBRARY_VERSION
        row["panel_version"] = PANEL_VERSION
        row["panel_index"] = index
        row["seq_num"] = index
        row["sequence_length"] = len(str(row["sequence"]))
        row["selection_key"] = "|".join(map(str, _candidate_key(row)))
        row["selection_is_score_blind"] = True
    return panel


def build_prediction_manifest(
    panel: Sequence[dict[str, Any]], alleles: Sequence[str] = ALLELES
) -> list[dict[str, Any]]:
    if len(panel) != 80:
        raise ValueError("prediction manifest requires the frozen 80-peptide panel")
    rows: list[dict[str, Any]] = []
    for allele in alleles:
        for seq_num, peptide in enumerate(panel, start=1):
            rows.append({
                "panel_version": PANEL_VERSION,
                "allele": allele,
                "candidate_id": peptide["candidate_id"],
                "seq_num": seq_num,
                "sequence": peptide["sequence"],
                "prediction_status": "pending_iedb_recommended_binding",
                "raw_response_file": "",
                "percentile_rank": "",
                "predicted_core": "",
                "core_start": "",
                "binding_rank_bin": "",
                "register_resolution": "pending",
            })
    return rows


def build_discovery_pairs(
    panel: Sequence[dict[str, Any]], alleles: Sequence[str] = ALLELES
) -> list[dict[str, Any]]:
    ebv = [row for row in panel if row["kingdom"] == "EBV"]
    human = [row for row in panel if row["kingdom"] == "human_self"]
    if len(ebv) != 40 or len(human) != 40:
        raise ValueError("discovery universe requires exactly 40 EBV and 40 self peptides")
    pairs: list[dict[str, Any]] = []
    for allele in alleles:
        for viral in ebv:
            for self_row in human:
                pairs.append({
                    "panel_version": PANEL_VERSION,
                    "allele": allele,
                    "ebv_allele": allele,
                    "self_allele": allele,
                    "pair_id": f"{allele}|{viral['candidate_id']}|{self_row['candidate_id']}",
                    "ebv_candidate_id": viral["candidate_id"],
                    "self_candidate_id": self_row["candidate_id"],
                    "geometry_status": "pending_models_and_resolved_registers",
                    "analysis_class": "exploratory_within_allele",
                })
    return pairs


def _safe_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def build_af3_jobs(
    panel: Sequence[dict[str, Any]], hla_sequences: dict[str, dict[str, str]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[list[dict[str, Any]]]]:
    if len(panel) != 80:
        raise ValueError("AlphaFold job builder requires exactly 80 peptides")
    jobs: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    for peptide in panel:
        for allele in ALLELES:
            if allele not in hla_sequences:
                raise ValueError(f"missing HLA sequence for {allele}")
            hla = hla_sequences[allele]
            chains = (hla["dra_sequence"], hla["drb_sequence"], str(peptide["sequence"]))
            if not all(chains):
                raise ValueError(f"empty chain sequence for {allele}")
            name = f"ebvms_v2_{_safe_token(allele)}_{_safe_token(str(peptide['candidate_id']))}"
            job = {
                "name": name,
                "modelSeeds": [],
                "sequences": [
                    {"proteinChain": {"sequence": sequence, "count": 1}}
                    for sequence in chains
                ],
                "dialect": "alphafoldserver",
                "version": 1,
            }
            jobs.append(job)
            inventory.append({
                "panel_version": PANEL_VERSION,
                "job_name": name,
                "allele": allele,
                "candidate_id": peptide["candidate_id"],
                "peptide_sequence": peptide["sequence"],
                "server_seed": "assigned_by_server",
                "model_status": "prepared_not_submitted",
                "expected_sample_count": 5,
                "observed_sample_count": 0,
                "canonical_model_key": "(panel_version,allele,candidate_id,server_seed,sample_index)",
            })
    if len(jobs) != 320 or len({job["name"] for job in jobs}) != 320:
        raise ValueError("expected exactly 320 unique AlphaFold jobs")
    batches = [jobs[index : index + 30] for index in range(0, len(jobs), 30)]
    if [len(batch) for batch in batches] != [30] * 10 + [20]:
        raise AssertionError("unexpected 320-job batch partition")
    return jobs, inventory, batches


def _composition_distance(left: str, right: str) -> float:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    left_counts = Counter(left)
    right_counts = Counter(right)
    return sum(abs(left_counts[aa] / len(left) - right_counts[aa] / len(right)) for aa in alphabet)


def select_native_controls(
    target: dict[str, Any], pool: Sequence[dict[str, Any]], *, count: int = 5
) -> list[dict[str, Any]]:
    """Choose fixed score-blind controls by composition, length, then ID."""
    target_sequence = str(target["sequence"])
    target_bin = str(target["binding_rank_bin"])
    eligible = []
    for row in pool:
        if _truth(row.get("excluded_source", False)):
            continue
        sequence = str(row.get("sequence", ""))
        if not sequence or abs(len(sequence) - len(target_sequence)) > 1:
            continue
        if str(row.get("binding_rank_bin", "")) != target_bin:
            continue
        selected = dict(row)
        selected["composition_distance"] = round(_composition_distance(target_sequence, sequence), 12)
        selected["length_difference"] = abs(len(sequence) - len(target_sequence))
        selected["selection_key"] = (
            selected["composition_distance"],
            selected["length_difference"],
            str(row["candidate_id"]),
        )
        eligible.append(selected)
    eligible.sort(key=lambda row: row["selection_key"])
    if len(eligible) < count:
        raise ValueError("insufficient controls under frozen length and binding-bin rule")
    return eligible[:count]


def build_native_calibration_jobs(
    entities: Sequence[dict[str, Any]], seeds: Sequence[int] = CALIBRATION_SEEDS
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if tuple(seeds) != CALIBRATION_SEEDS:
        raise ValueError(f"native calibration seeds are frozen as {CALIBRATION_SEEDS}")
    if len(entities) != 12:
        raise ValueError("native calibration requires exactly 12 entities")
    per_arm = Counter(str(row["arm"]) for row in entities)
    if per_arm != Counter({"viral": 6, "self": 6}):
        raise ValueError("native calibration requires six viral and six self entities")
    keys = [(str(row["allele"]), str(row["entity_id"])) for row in entities]
    if len(keys) != len(set(keys)):
        raise ValueError("native calibration entity IDs must be unique within allele")
    jobs: list[dict[str, Any]] = []
    manifest: list[dict[str, Any]] = []
    for entity in sorted(entities, key=lambda row: (str(row["arm"]), str(row["entity_id"]))):
        allele = str(entity["allele"])
        if entity["arm"] == "viral" and allele != "HLA-DRB5*01:01":
            raise ValueError("viral calibration arm must use native HLA-DRB5*01:01")
        if entity["arm"] == "self" and allele != "HLA-DRB1*15:01":
            raise ValueError("self calibration arm must use native HLA-DRB1*15:01")
        for seed in seeds:
            name = f"ebvms_native_{_safe_token(str(entity['entity_id']))}_s{seed}"
            jobs.append({
                "name": name,
                "modelSeeds": [int(seed)],
                "sequences": [
                    {"proteinChain": {"sequence": str(sequence), "count": 1}}
                    for sequence in (
                        entity["dra_sequence"],
                        entity["drb_sequence"],
                        entity["sequence"],
                    )
                ],
                "dialect": "alphafoldserver",
                "version": 1,
            })
            manifest.append({
                "job_name": name,
                "entity_id": entity["entity_id"],
                "entity_role": entity["entity_role"],
                "arm": entity["arm"],
                "allele": allele,
                "peptide_sequence": entity["sequence"],
                "server_seed": int(seed),
                "status": "prepared_not_submitted",
                "cross_hla_comparison_allowed": True,
                "claim_boundary": "declared DRB5/DRB1 positive-control calibration only",
            })
    return jobs, manifest


def build_calibration_comparison_universe(
    manifest: Sequence[dict[str, Any]], seeds: Sequence[int] = CALIBRATION_SEEDS
) -> list[dict[str, Any]]:
    """Enumerate the frozen 26-pair primary and 10-pair sensitivity sets per seed."""
    output: list[dict[str, Any]] = []
    for seed in seeds:
        rows = [row for row in manifest if int(row["server_seed"]) == int(seed)]
        viral = [row for row in rows if row["arm"] == "viral"]
        self_rows = [row for row in rows if row["arm"] == "self"]
        viral_positive = [row for row in viral if str(row["entity_role"]).startswith("E1_positive")]
        self_positive = [row for row in self_rows if str(row["entity_role"]).startswith("E1_positive")]
        viral_controls = [row for row in viral if not str(row["entity_role"]).startswith("E1_positive")]
        self_controls = [row for row in self_rows if not str(row["entity_role"]).startswith("E1_positive")]
        if not (len(viral_positive) == len(self_positive) == 1 and len(viral_controls) == len(self_controls) == 5):
            raise ValueError("calibration comparison universe requires one positive and five controls per arm per seed")

        def add(left: dict[str, Any], right: dict[str, Any], analysis_set: str, pair_role: str) -> None:
            output.append({
                "seed": int(seed),
                "analysis_set": analysis_set,
                "pair_role": pair_role,
                "viral_entity_id": left["entity_id"],
                "self_entity_id": right["entity_id"],
                "pair_id": f"s{seed}|{left['entity_id']}|{right['entity_id']}",
                "geometry_status": "pending_two_models_and_equivalent_cross_hla_groove_fit",
                "ranking_endpoint": "median TCR-facing P2,P3,P5,P7,P8 C-alpha RMSD",
            })

        add(viral_positive[0], self_positive[0], "primary_rank_of_26", "E1_positive")
        for left in sorted(viral_controls, key=lambda row: str(row["entity_id"])):
            for right in sorted(self_controls, key=lambda row: str(row["entity_id"])):
                add(left, right, "primary_rank_of_26", "full_decoy")
        for right in sorted(self_controls, key=lambda row: str(row["entity_id"])):
            add(viral_positive[0], right, "single_arm_sensitivity", "positive_viral_vs_self_control")
        for left in sorted(viral_controls, key=lambda row: str(row["entity_id"])):
            add(left, self_positive[0], "single_arm_sensitivity", "viral_control_vs_positive_self")
    return output


def classify_positive_recovery(seed_summaries: Sequence[dict[str, Any]]) -> str:
    if not seed_summaries:
        return "not_evaluable"
    by_seed = {int(row["seed"]): row for row in seed_summaries}
    if set(by_seed) != set(CALIBRATION_SEEDS):
        return "not_evaluable"
    recovered = all(
        int(row["rank_of_26"]) <= 3
        and float(row["positive_rmsd"]) < float(row["equal_weight_control_median"])
        for row in by_seed.values()
    )
    return "recovered" if recovered else "failed_calibration"


def _rmsd(left: np.ndarray, right: np.ndarray) -> float:
    return float(math.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def groove_fitted_rmsd(
    left_groove: np.ndarray,
    left_core: np.ndarray,
    right_groove: np.ndarray,
    right_core: np.ndarray,
) -> dict[str, float]:
    """Fit right to left using equivalent groove atoms, then compare P1-P9."""
    arrays = [np.asarray(value, dtype=float) for value in (left_groove, left_core, right_groove, right_core)]
    left_groove, left_core, right_groove, right_core = arrays
    if left_groove.shape != right_groove.shape or left_groove.ndim != 2 or left_groove.shape[1] != 3:
        raise ValueError("groove coordinates must be equivalent N x 3 arrays")
    if left_core.shape != (9, 3) or right_core.shape != (9, 3):
        raise ValueError("P1-P9 coordinates must be 9 x 3 arrays")
    left_center = left_groove.mean(axis=0)
    right_center = right_groove.mean(axis=0)
    covariance = (right_groove - right_center).T @ (left_groove - left_center)
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    aligned_core = (right_core - right_center) @ rotation + left_center
    anchor_idx = np.array([position - 1 for position in ANCHOR_POSITIONS])
    exposed_idx = np.array([position - 1 for position in EXPOSED_POSITIONS])
    return {
        "full_core_rmsd": _rmsd(left_core, aligned_core),
        "anchor_rmsd": _rmsd(left_core[anchor_idx], aligned_core[anchor_idx]),
        "exposed_rmsd": _rmsd(left_core[exposed_idx], aligned_core[exposed_idx]),
    }


def rows_sha256(rows: Iterable[dict[str, Any]]) -> str:
    """Stable in-memory checksum helper used by integration verification."""
    normalized = "\n".join(
        "\t".join(f"{key}={row[key]}" for key in sorted(row))
        for row in rows
    )
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()
