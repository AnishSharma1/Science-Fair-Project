"""Deterministic N3 ranking context for the frozen high-yield HLA-II shortlist.

N3 comparators have unknown TCR-recognition status. They provide a fair local
ranking background but are not biological negatives and cannot establish
specificity, cross-reactivity, or molecular mimicry.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections import Counter
from typing import Any, Mapping, Sequence


SEED = 271828
CLAIM_BOUNDARY = (
    "N3 panels provide descriptive, HLA-specific computational ranking context only; "
    "they are not evidence of presentation, TCR recognition, activation, specificity, "
    "cross-reactivity, molecular mimicry, MS mechanism, probability, or false-discovery rate."
)
SURFACE_FEATURES = (
    "sidechain_exposure_mismatch_fraction_q75",
    "exposure_weighted_centroid_rmsd_A_q75",
    "exposure_weighted_orientation_rmsd_A_q75",
    "exposed_distance_matrix_rmsd_A_q75",
    "exposure_weighted_chemistry_mismatch_q75",
    "exposure_weighted_backbone_rmsd_A_q75",
)

FROZEN_TARGETS = (
    {
        "target_id": "HY03_SEQ_01",
        "allele": "HLA-DRB1*03:01",
        "lane": "sequence",
        "ebv_core": "IVRQSRGDR",
        "self_core": "LLKDAIGEG",
        "pair_id": "HLA-DRB1*03:01|EBV_IEDB_73663fe4c7bb|SELF_IEDB_842d700c6a24",
    },
    {
        "target_id": "HY03_SEQ_02",
        "allele": "HLA-DRB1*03:01",
        "lane": "sequence",
        "ebv_core": "VTLTSYWRR",
        "self_core": "IAIHHPWIR",
        "pair_id": "HLA-DRB1*03:01|EBV_IEDB_156b80230e03|SELF_CANON_CRYAB_0001_0015",
    },
    {
        "target_id": "HY03_REG_01",
        "allele": "HLA-DRB1*03:01",
        "lane": "register",
        "ebv_core": "IWMCMTVRH",
        "self_core": "LSFDKDAMV",
        "pair_id": "HLA-DRB1*03:01|EBV_IEDB_2fbac9820d7d|SELF_CANON_TALDO1_0108_0122",
    },
    {
        "target_id": "HY08_SEQ_01",
        "allele": "HLA-DRB1*08:01",
        "lane": "sequence",
        "ebv_core": "LRALLARSH",
        "self_core": "LEARLSRMH",
        "pair_id": "HLA-DRB1*08:01|EBV_IEDB_ce143ad82ead|SELF_CANON_ANO2_0079_0093",
    },
    {
        "target_id": "HY08_SEQ_02",
        "allele": "HLA-DRB1*08:01",
        "lane": "sequence",
        "ebv_core": "VRRRVLVQQ",
        "self_core": "FSRVVHLYR",
        "pair_id": "HLA-DRB1*08:01|EBV_IEDB_91869932d81f|SELF_IEDB_96777e34f994",
    },
    {
        "target_id": "HY08_REG_01",
        "allele": "HLA-DRB1*08:01",
        "lane": "register",
        "ebv_core": "LRALLARSH",
        "self_core": "YSKAFTLTI",
        "pair_id": "HLA-DRB1*08:01|EBV_IEDB_ce143ad82ead|SELF_CANON_CNP_0272_0286",
    },
    {
        "target_id": "HY13_SEQ_01",
        "allele": "HLA-DRB1*13:03",
        "lane": "sequence",
        "ebv_core": "WMCMTVRHR",
        "self_core": "IICYNWLHR",
        "pair_id": "HLA-DRB1*13:03|EBV_IEDB_2fbac9820d7d|SELF_IEDB_18eebcfe8ed0",
    },
    {
        "target_id": "HY13_SEQ_02",
        "allele": "HLA-DRB1*13:03",
        "lane": "sequence",
        "ebv_core": "YHFVKKHVH",
        "self_core": "LSFDKDAMV",
        "pair_id": "HLA-DRB1*13:03|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0108_0122",
    },
    {
        "target_id": "HY13_REG_01",
        "allele": "HLA-DRB1*13:03",
        "lane": "register",
        "ebv_core": "LRALLARSH",
        "self_core": "YSKAFTLTI",
        "pair_id": "HLA-DRB1*13:03|EBV_IEDB_ce143ad82ead|SELF_CANON_CNP_0272_0286",
    },
    {
        "target_id": "HY15_SEQ_01",
        "allele": "HLA-DRB1*15:01",
        "lane": "sequence",
        "ebv_core": "ILIYNGWYA",
        "self_core": "IAIHHPWIR",
        "pair_id": "HLA-DRB1*15:01|EBV_IEDB_bce490b86fd9|SELF_CANON_CRYAB_0001_0015",
    },
    {
        "target_id": "HY15_SEQ_02",
        "allele": "HLA-DRB1*15:01",
        "lane": "sequence",
        "ebv_core": "VYHFVKKHV",
        "self_core": "IYNYYKKFS",
        "pair_id": "HLA-DRB1*15:01|EBV_IEDB_35bb9c18fac4|SELF_CANON_TALDO1_0216_0230",
    },
    {
        "target_id": "HY15_REG_01",
        "allele": "HLA-DRB1*15:01",
        "lane": "register",
        "ebv_core": "VTLTSYWRR",
        "self_core": "VLLLESHCA",
        "pair_id": "HLA-DRB1*15:01|EBV_IEDB_156b80230e03|SELF_CANON_MAG_0408_0422",
    },
)


def _number(value: Any, default: float = math.inf) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _slug(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(value).upper())


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _average_tie_percentiles(
    rows: Sequence[Mapping[str, Any]], field: str, *, higher_is_better: bool = False
) -> dict[str, float]:
    eligible = [row for row in rows if math.isfinite(_number(row.get(field)))]
    ordered = sorted(
        eligible,
        key=lambda row: (
            -_number(row.get(field)) if higher_is_better else _number(row.get(field)),
            str(row["pair_id"]),
        ),
    )
    output: dict[str, float] = {}
    index = 0
    while index < len(ordered):
        end = index + 1
        value = _number(ordered[index].get(field))
        while end < len(ordered) and _number(ordered[end].get(field)) == value:
            end += 1
        average_rank = ((index + 1) + end) / 2.0
        percentile = 0.0 if len(ordered) == 1 else (average_rank - 1.0) / (len(ordered) - 1.0)
        for row in ordered[index:end]:
            output[str(row["pair_id"])] = round(percentile, 12)
        index = end
    return output


def _ordinal_ranks(
    rows: Sequence[Mapping[str, Any]], field: str, *, higher_is_better: bool = False
) -> dict[str, int]:
    ordered = sorted(
        rows,
        key=lambda row: (
            -_number(row.get(field), -math.inf)
            if higher_is_better
            else _number(row.get(field)),
            str(row["pair_id"]),
        ),
    )
    return {str(row["pair_id"]): index for index, row in enumerate(ordered, start=1)}


def validate_frozen_targets(targets: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    pair_ids = [str(row["pair_id"]) for row in targets]
    target_ids = [str(row["target_id"]) for row in targets]
    if len(pair_ids) != len(set(pair_ids)) or len(target_ids) != len(set(target_ids)):
        raise ValueError("duplicate frozen target")
    counts = Counter(str(row["allele"]) for row in targets)
    if len(targets) != 12 or set(counts.values()) != {3}:
        raise ValueError("frozen registry must contain exactly three targets per HLA")
    return {"target_count": len(targets), "targets_per_hla": dict(sorted(counts.items()))}


def select_comparator_arms(
    target: Mapping[str, Any],
    candidates: Sequence[Mapping[str, Any]],
    *,
    allele: str,
    arm_class: str,
    excluded_candidate_ids: set[str],
    count: int = 5,
    seed: int = SEED,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Select comparator arms without reading geometry or pair scores."""
    if count < 1:
        raise ValueError("comparator count must be positive")
    target_length = len(str(target["sequence"]))
    target_binding = _number(target.get("binding_percentile"))
    provenance: list[dict[str, Any]] = []
    eligible: list[dict[str, Any]] = []
    for source in candidates:
        row = dict(source)
        candidate_id = str(row["candidate_id"])
        reason = "eligible"
        if str(row.get("allele")) != allele:
            reason = "wrong_hla"
        elif candidate_id in excluded_candidate_ids:
            reason = "excluded_frozen_or_control_arm"
        elif _number(row.get("binding_percentile")) > 20:
            reason = "binding_percentile_above_20"
        elif int(_number(row.get("model_count"), 0)) != 5:
            reason = "incomplete_model_ensemble"
        elif str(row.get("surface_status")) != "complete":
            reason = "incomplete_surface_features"
        digest = hashlib.sha256(f"{seed}|{arm_class}|{candidate_id}".encode("utf-8")).hexdigest()
        record = {
            **row,
            "arm_class": arm_class,
            "eligibility_reason": reason,
            "selection_length_difference": abs(len(str(row.get("sequence", ""))) - target_length),
            "selection_binding_percentile_difference": abs(
                _number(row.get("binding_percentile")) - target_binding
            ),
            "selection_seeded_hash": digest,
            "selection_is_score_blind": True,
            "selected": False,
        }
        provenance.append(record)
        if reason == "eligible":
            eligible.append(record)
    eligible.sort(
        key=lambda row: (
            int(row["selection_length_difference"]),
            float(row["selection_binding_percentile_difference"]),
            str(row["selection_seeded_hash"]),
            str(row["candidate_id"]),
        )
    )
    selected: list[dict[str, Any]] = []
    seen_identity: set[tuple[str, str]] = set()
    for row in eligible:
        identity = (
            _slug(row.get("accession") or row.get("protein")),
            _slug(row.get("core")),
        )
        if identity in seen_identity:
            row["eligibility_reason"] = "duplicate_accession_core"
            continue
        seen_identity.add(identity)
        row["selected"] = True
        row["selection_order"] = len(selected) + 1
        selected.append(dict(row))
        if len(selected) == count:
            break
    selected_ids = {str(row["candidate_id"]) for row in selected}
    for row in provenance:
        if str(row["candidate_id"]) in selected_ids:
            chosen = next(item for item in selected if item["candidate_id"] == row["candidate_id"])
            row["selected"] = True
            row["selection_order"] = chosen["selection_order"]
        elif row["eligibility_reason"] == "eligible":
            identity = (
                _slug(row.get("accession") or row.get("protein")),
                _slug(row.get("core")),
            )
            if sum(
                1
                for item in eligible
                if (
                    _slug(item.get("accession") or item.get("protein")),
                    _slug(item.get("core")),
                )
                == identity
                and str(item["candidate_id"]) in selected_ids
            ):
                row["eligibility_reason"] = "duplicate_accession_core"
            else:
                row["eligibility_reason"] = "eligible_not_selected_after_score_blind_ordering"
    if len(selected) < count:
        raise ValueError(f"insufficient eligible {arm_class} comparators")
    return selected, provenance


def build_n3_panel(
    target_pair: Mapping[str, Any],
    ebv_arms: Sequence[Mapping[str, Any]],
    self_arms: Sequence[Mapping[str, Any]],
    pair_lookup: Mapping[tuple[str, str], Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if len(ebv_arms) != 5 or len(self_arms) != 5:
        raise ValueError("a complete panel requires five comparator arms per side")
    panel = [{**dict(target_pair), "row_role": "target", "n3_specificity_role": "excluded"}]
    for left in ebv_arms:
        for right in self_arms:
            key = (str(left["candidate_id"]), str(right["candidate_id"]))
            if key not in pair_lookup:
                raise ValueError(f"missing exact-HLA V3 pair for {key[0]} and {key[1]}")
            panel.append({
                **dict(pair_lookup[key]),
                "row_role": "n3",
                "n3_specificity_role": "unknown_recognition_not_a_specificity_negative",
            })
    pair_ids = [str(row["pair_id"]) for row in panel]
    if len(panel) != 26 or len(pair_ids) != len(set(pair_ids)):
        raise ValueError("panel must contain one target and 25 unique N3 pairs")
    return panel


def rank_panel_rows(rows: Sequence[Mapping[str, Any]], *, seed: int = SEED) -> list[dict[str, Any]]:
    if not rows:
        return []
    annotated = [dict(row) for row in rows]
    if all(all(math.isfinite(_number(row.get(field))) for field in SURFACE_FEATURES) for row in annotated):
        feature_maps = [_average_tie_percentiles(annotated, field) for field in SURFACE_FEATURES]
        for row in annotated:
            pair_id = str(row["pair_id"])
            row["panel_local_surface_percentile"] = round(
                sum(mapping[pair_id] for mapping in feature_maps) / len(feature_maps), 12
            )
    else:
        surface_map = _average_tie_percentiles(annotated, "local_surface_score")
        for row in annotated:
            row["panel_local_surface_percentile"] = surface_map[str(row["pair_id"])]
    ordered = sorted(
        annotated,
        key=lambda row: (
            -_number(row.get("tcr_facing_blosum62_similarity"), -math.inf),
            _number(row.get("tcr_face_physicochemical_mismatch")),
            -_number(row.get("tcr_facing_sequence_identity"), -math.inf),
            _number(row.get("panel_local_surface_percentile"))
            if _truth(row.get("register_robust"))
            else math.inf,
            str(row["pair_id"]),
        ),
    )
    primary = {str(row["pair_id"]): rank for rank, row in enumerate(ordered, start=1)}
    diagnostic_specs = {
        "panel_local_surface_rank": ("panel_local_surface_percentile", False),
        "panel_exposed_backbone_rank": ("exposure_weighted_backbone_rmsd_A_q75", False),
        "panel_full_core_rmsd_rank": ("full_core_ca_rmsd_A_q75", False),
        "panel_anchor_rmsd_rank": ("anchor_ca_rmsd_A_q75", False),
        "panel_physicochemical_rank": ("tcr_face_physicochemical_mismatch", False),
        "panel_identity_rank": ("tcr_facing_sequence_identity", True),
    }
    diagnostic = {
        output: _ordinal_ranks(annotated, field, higher_is_better=higher)
        for output, (field, higher) in diagnostic_specs.items()
    }
    random_order = sorted(
        annotated,
        key=lambda row: (
            hashlib.sha256(f"{seed}|random|{row['pair_id']}".encode("utf-8")).hexdigest(),
            str(row["pair_id"]),
        ),
    )
    random_rank = {str(row["pair_id"]): rank for rank, row in enumerate(random_order, start=1)}
    output = []
    for row in annotated:
        pair_id = str(row["pair_id"])
        row["panel_primary_rank"] = primary[pair_id]
        row["panel_primary_percentile"] = round((primary[pair_id] - 1) / max(1, len(rows) - 1), 12)
        row["panel_primary_structure_abstained"] = not _truth(row.get("register_robust"))
        for field, mapping in diagnostic.items():
            row[field] = mapping[pair_id]
        row["panel_random_rank"] = random_rank[pair_id]
        row["claim_boundary"] = CLAIM_BOUNDARY
        output.append(row)
    return output


def build_ranking_context_gate(panel_summaries: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    targets = []
    for source in panel_summaries:
        row = dict(source)
        if str(row.get("panel_status")) != "complete":
            status = "not_evaluable"
        else:
            status = (
                "rank_context_supportive"
                if int(row["target_primary_rank"]) <= 3
                else "rank_context_not_supportive"
            )
        targets.append({**row, "ranking_context_status": status})
    return {
        "status": "complete" if targets and all(row["ranking_context_status"] != "not_evaluable" for row in targets) else "not_evaluable",
        "target_count": len(targets),
        "supportive_target_count": sum(row["ranking_context_status"] == "rank_context_supportive" for row in targets),
        "targets": targets,
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "n3_is_unknown_recognition_not_negative_control": True,
        "claim_boundary": CLAIM_BOUNDARY,
    }


def build_specificity_gate(explicit_n1_n2_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    return {
        "status": "not_evaluable_no_explicit_n1_n2" if not explicit_n1_n2_rows else "registry_available_not_analyzed",
        "explicit_n1_n2_count": len(explicit_n1_n2_rows),
        "n3_excluded_from_specificity": True,
        "specificity_claim_allowed": False,
        "discovery_unlock_allowed": False,
        "claim_boundary": "Ranking performance and N3 comparisons are not specificity evidence.",
    }
