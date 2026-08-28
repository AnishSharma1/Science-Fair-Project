"""Register-sensitivity analysis for sequence-supported high-yield HLA-II pairs.

Alternate windows are sequence-only sensitivity analyses. They are not alternate
AlphaFold models, experimentally resolved registers, or evidence of TCR binding.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from build_same_register_hla_rankings_v2 import sequence_metrics
from high_yield_control_validation import CLAIM_BOUNDARY, rank_panel_rows


REGISTER_CLAIM_BOUNDARY = (
    "Window rescoring tests dependence on the assumed P1-P9 register using the frozen "
    "sequence ranking and N3 panel. Alternate windows were not structurally modeled or "
    "experimentally resolved. Results do not establish presentation, TCR recognition, "
    "specificity, cross-reactivity, molecular mimicry, or MS mechanism."
)


def enumerate_windows(sequence: str, *, width: int = 9) -> list[dict[str, Any]]:
    normalized = str(sequence).strip().upper()
    if width < 1:
        raise ValueError("window width must be positive")
    if len(normalized) < width:
        raise ValueError("sequence is shorter than window width")
    return [
        {"start_1_based": start + 1, "core": normalized[start : start + width]}
        for start in range(len(normalized) - width + 1)
    ]


def _integer(row: Mapping[str, Any], field: str) -> int:
    try:
        return int(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"missing or invalid {field}") from error


def _validate_declared_core(
    sequence: str,
    start_1_based: int,
    declared_core: str,
    *,
    side: str,
) -> None:
    observed = sequence[start_1_based - 1 : start_1_based + 8]
    if len(observed) != 9 or observed != declared_core:
        raise ValueError(
            f"{side} declared core does not match the exact sequence at its one-based start"
        )


def evaluate_target_windows(
    target: Mapping[str, Any],
    frozen_panel_rows: Sequence[Mapping[str, Any]],
    *,
    local_shift: int = 1,
) -> list[dict[str, Any]]:
    """Rescore all window pairs against one unchanged N3 panel."""
    if local_shift < 0:
        raise ValueError("local shift must be nonnegative")
    target_rows = [row for row in frozen_panel_rows if str(row.get("row_role")) == "target"]
    n3_rows = [row for row in frozen_panel_rows if str(row.get("row_role")) == "n3"]
    if len(target_rows) != 1 or not n3_rows:
        raise ValueError("panel must contain one target and at least one N3 row")
    pair_id = str(target["pair_id"])
    if str(target_rows[0].get("pair_id")) != pair_id:
        raise ValueError("frozen panel target does not match the target registry")

    ebv_sequence = str(target["ebv_sequence"]).strip().upper()
    self_sequence = str(target["self_sequence"]).strip().upper()
    ebv_declared_start = _integer(target, "ebv_declared_core_start_1_based")
    self_declared_start = _integer(target, "self_declared_core_start_1_based")
    ebv_declared_core = str(target["ebv_core_p1_p9"]).strip().upper()
    self_declared_core = str(target["self_core_p1_p9"]).strip().upper()
    _validate_declared_core(
        ebv_sequence, ebv_declared_start, ebv_declared_core, side="EBV"
    )
    _validate_declared_core(
        self_sequence, self_declared_start, self_declared_core, side="self"
    )

    output: list[dict[str, Any]] = []
    for ebv_window in enumerate_windows(ebv_sequence):
        for self_window in enumerate_windows(self_sequence):
            is_declared = (
                ebv_window["start_1_based"] == ebv_declared_start
                and self_window["start_1_based"] == self_declared_start
            )
            is_local = (
                abs(ebv_window["start_1_based"] - ebv_declared_start) <= local_shift
                and abs(self_window["start_1_based"] - self_declared_start) <= local_shift
            )
            metrics = sequence_metrics(ebv_window["core"], self_window["core"])
            variant = {
                **dict(target_rows[0]),
                **metrics,
                "ebv_core": ebv_window["core"],
                "ebv_core_p1_p9": ebv_window["core"],
                "self_core": self_window["core"],
                "self_core_p1_p9": self_window["core"],
                "register_robust": False,
            }
            ranked = rank_panel_rows([variant, *(dict(row) for row in n3_rows)])
            ranked_target = next(row for row in ranked if str(row["pair_id"]) == pair_id)
            output.append(
                {
                    "target_id": str(target["target_id"]),
                    "allele": str(target["allele"]),
                    "pair_id": pair_id,
                    "ebv_window_start_1_based": ebv_window["start_1_based"],
                    "self_window_start_1_based": self_window["start_1_based"],
                    "ebv_window_core": ebv_window["core"],
                    "self_window_core": self_window["core"],
                    "is_declared_window_pair": is_declared,
                    "is_local_shift_window_pair": is_local,
                    "ebv_shift_from_declared": (
                        ebv_window["start_1_based"] - ebv_declared_start
                    ),
                    "self_shift_from_declared": (
                        self_window["start_1_based"] - self_declared_start
                    ),
                    **{key: round(float(value), 12) for key, value in metrics.items()},
                    "panel_primary_rank": int(ranked_target["panel_primary_rank"]),
                    "panel_primary_percentile": ranked_target["panel_primary_percentile"],
                    "capture_at_3": int(ranked_target["panel_primary_rank"]) <= 3,
                    "frozen_n3_pair_count": len(n3_rows),
                    "structure_abstained_for_alternate_register": True,
                    "sensitivity_scope": "all_fully_contained_p1_p9_window_pairs",
                    "claim_boundary": REGISTER_CLAIM_BOUNDARY,
                }
            )
    return output


def summarize_target_windows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if not rows:
        raise ValueError("target window rows are required")
    target_ids = {str(row["target_id"]) for row in rows}
    if len(target_ids) != 1:
        raise ValueError("window summary must contain one target")
    declared = [row for row in rows if bool(row["is_declared_window_pair"])]
    local = [row for row in rows if bool(row["is_local_shift_window_pair"])]
    if len(declared) != 1 or not local:
        raise ValueError("window summary requires one declared row and local-shift rows")

    declared_capture = bool(declared[0]["capture_at_3"])
    local_capture = all(bool(row["capture_at_3"]) for row in local)
    all_capture = all(bool(row["capture_at_3"]) for row in rows)
    if all_capture:
        status = "all_window_robust"
    elif local_capture:
        status = "local_shift_robust_only"
    elif declared_capture:
        status = "declared_window_only"
    else:
        status = "rank_context_not_supportive"

    best = min(rows, key=lambda row: (int(row["panel_primary_rank"]), str(row["ebv_window_core"]), str(row["self_window_core"])))
    worst = max(rows, key=lambda row: (int(row["panel_primary_rank"]), str(row["ebv_window_core"]), str(row["self_window_core"])))
    return {
        "target_id": str(declared[0]["target_id"]),
        "allele": str(declared[0].get("allele", "")),
        "pair_id": str(declared[0].get("pair_id", "")),
        "declared_ebv_core": str(declared[0].get("ebv_window_core", "")),
        "declared_self_core": str(declared[0].get("self_window_core", "")),
        "declared_window_rank": int(declared[0]["panel_primary_rank"]),
        "declared_window_capture_at_3": declared_capture,
        "local_shift_window_pair_count": len(local),
        "local_shift_capture_at_3_count": sum(bool(row["capture_at_3"]) for row in local),
        "local_shift_capture_fraction": round(
            sum(bool(row["capture_at_3"]) for row in local) / len(local), 12
        ),
        "worst_local_shift_rank": max(int(row["panel_primary_rank"]) for row in local),
        "all_window_pair_count": len(rows),
        "all_window_capture_at_3_count": sum(bool(row["capture_at_3"]) for row in rows),
        "all_window_capture_fraction": round(
            sum(bool(row["capture_at_3"]) for row in rows) / len(rows), 12
        ),
        "best_all_window_rank": int(best["panel_primary_rank"]),
        "best_all_window_ebv_core": str(best["ebv_window_core"]),
        "best_all_window_self_core": str(best["self_window_core"]),
        "worst_all_window_rank": int(worst["panel_primary_rank"]),
        "worst_all_window_ebv_core": str(worst["ebv_window_core"]),
        "worst_all_window_self_core": str(worst["self_window_core"]),
        "register_resolution_status": status,
        "experimentally_resolved_register": False,
        "claim_boundary": REGISTER_CLAIM_BOUNDARY,
    }


def build_register_resolution_gate(
    summaries: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    statuses = [str(row["register_resolution_status"]) for row in summaries]
    if not statuses or "not_evaluable" in statuses:
        status = "not_evaluable"
    elif all(item == "all_window_robust" for item in statuses):
        status = "supportive_all_windows"
    elif all(
        item in {"all_window_robust", "local_shift_robust_only"}
        for item in statuses
    ):
        status = "supportive_local_shifts"
    elif all(item != "rank_context_not_supportive" for item in statuses):
        status = "declared_register_dependent"
    else:
        status = "rank_context_not_supportive"
    return {
        "gate_name": "high_yield_register_resolution_sensitivity",
        "status": status,
        "target_count": len(statuses),
        "all_window_robust_count": statuses.count("all_window_robust"),
        "local_shift_robust_only_count": statuses.count("local_shift_robust_only"),
        "declared_window_only_count": statuses.count("declared_window_only"),
        "rank_context_not_supportive_count": statuses.count("rank_context_not_supportive"),
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "independent_validation_claim_allowed": False,
        "interpretation": REGISTER_CLAIM_BOUNDARY,
        "upstream_n3_claim_boundary": CLAIM_BOUNDARY,
    }


def prioritize_register_confirmation(
    summaries: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Order assay work by register sensitivity, not by a new discovery score."""
    ordered = sorted(
        (dict(row) for row in summaries),
        key=lambda row: (
            -float(row["local_shift_capture_fraction"]),
            int(row["declared_window_rank"]),
            int(row["worst_local_shift_rank"]),
            str(row["target_id"]),
        ),
    )
    for rank, row in enumerate(ordered, start=1):
        row["experimental_priority_rank"] = rank
        row["priority_basis"] = (
            "higher local-shift top-three fraction, then declared rank, then worst local rank"
        )
        row["priority_is_not_a_discovery_rerank"] = True
    return ordered


def build_experimental_peptide_panel(
    targets: Sequence[Mapping[str, Any]],
    *,
    local_shift: int = 1,
) -> list[dict[str, Any]]:
    """Create a proposed nested-peptide design; no synthesis/order is performed."""
    output: list[dict[str, Any]] = []
    for target in targets:
        for side in ("ebv", "self"):
            sequence = str(target[f"{side}_sequence"]).strip().upper()
            declared_start = _integer(target, f"{side}_declared_core_start_1_based")
            output.append(
                {
                    "panel_peptide_id": f"{target['target_id']}|{side.upper()}|PARENT",
                    "target_id": str(target["target_id"]),
                    "allele": str(target["allele"]),
                    "arm": side,
                    "peptide_role": "parent_assay_peptide",
                    "shift_from_declared": "",
                    "core_start_1_based": "",
                    "sequence": sequence,
                    "length": len(sequence),
                    "assay_purpose": "reference binding and response peptide; does not distinguish register alone",
                    "proposed_not_ordered": True,
                    "claim_boundary": REGISTER_CLAIM_BOUNDARY,
                }
            )
            windows = {
                int(row["start_1_based"]): str(row["core"])
                for row in enumerate_windows(sequence)
            }
            for shift in range(-local_shift, local_shift + 1):
                start = declared_start + shift
                if start not in windows:
                    continue
                output.append(
                    {
                        "panel_peptide_id": (
                            f"{target['target_id']}|{side.upper()}|CORE|{start:02d}"
                        ),
                        "target_id": str(target["target_id"]),
                        "allele": str(target["allele"]),
                        "arm": side,
                        "peptide_role": "register_discrimination_core",
                        "shift_from_declared": shift,
                        "core_start_1_based": start,
                        "sequence": windows[start],
                        "length": 9,
                        "assay_purpose": "nested core for direct register-sensitivity comparison",
                        "proposed_not_ordered": True,
                        "claim_boundary": REGISTER_CLAIM_BOUNDARY,
                    }
                )
    identifiers = [str(row["panel_peptide_id"]) for row in output]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("experimental peptide panel IDs must be unique")
    return output
