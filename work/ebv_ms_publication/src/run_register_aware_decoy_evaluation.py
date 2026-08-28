"""Evaluate frozen strict decoy sets without reselection or pooled inference."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
BENCHMARK = PROC / "register_aware_benchmark"
OUT = PROC / "register_aware_scoring"
PRIMARY_SCORE = "same_register_property_similarity"
ROBUST_STATUS = "robust_primary_ranking_eligible"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _secondary(score_row: dict[str, Any] | None, field: str) -> object:
    return "" if score_row is None else score_row.get(field, "")


def _base_row(feasibility: dict[str, Any], decoy_rows: list[dict[str, Any]]) -> dict[str, object]:
    return {
        "target_pair_id": feasibility["target_pair_id"],
        "target_validation_label": feasibility["target_validation_label"],
        "target_register_assessment": feasibility["target_register_assessment"],
        "frozen_selected_decoy_count": feasibility.get("selected_decoy_count", ""),
        "observed_decoy_row_count": len(decoy_rows),
        "strict_decoy_pair_ids": ";".join(
            row["decoy_pair_id"]
            for row in sorted(decoy_rows, key=lambda item: int(item["decoy_ordinal"]))
        ),
        "primary_endpoint": PRIMARY_SCORE,
        "target_score": "",
        "decoy_score_median": "",
        "decoy_score_min": "",
        "decoy_score_max": "",
        "target_minus_decoy_median": "",
        "target_rank_among_six": "",
        "descriptive_within_set_rank_fraction": "",
        "target_identity_fraction": "",
        "target_anchor_property_similarity": "",
        "target_candidate_exposed_property_similarity": "",
    }


def evaluate_targets(
    score_rows: list[dict[str, Any]],
    decoy_rows: list[dict[str, Any]],
    feasibility_rows: list[dict[str, Any]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Return one strict-set decision per target and a no-pooling summary."""
    score_by_pair = {str(row["pair_id"]): row for row in score_rows}
    if len(score_by_pair) != len(score_rows):
        raise ValueError("score table has duplicate pair IDs")
    decoys_by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in decoy_rows:
        decoys_by_target[str(row["target_pair_id"])].append(row)

    comparisons: list[dict[str, object]] = []
    evaluable_labels: set[str] = set()
    for feasibility in feasibility_rows:
        target_id = str(feasibility["target_pair_id"])
        selected = decoys_by_target.get(target_id, [])
        result = _base_row(feasibility, selected)
        target_score_row = score_by_pair.get(target_id)
        if feasibility["target_register_assessment"] != "assessable_register_hypothesis":
            result["evaluation_status"] = "not_evaluable_register_status"
        elif len(selected) != 5 or feasibility.get("selected_decoy_count") != "5":
            result["evaluation_status"] = "not_evaluable_incomplete_strict_decoy_set"
        elif target_score_row is None or target_score_row.get("score_coverage_status") != ROBUST_STATUS:
            result["evaluation_status"] = "not_evaluable_target_score_coverage"
        else:
            decoy_scores = [score_by_pair.get(str(row["decoy_pair_id"])) for row in selected]
            if any(row is None or row.get("score_coverage_status") != ROBUST_STATUS for row in decoy_scores):
                result["evaluation_status"] = "not_evaluable_decoy_score_coverage"
            else:
                target_value = float(target_score_row[PRIMARY_SCORE])
                values = [float(row[PRIMARY_SCORE]) for row in decoy_scores if row is not None]
                ranked = sorted(
                    [(target_value, target_id), *[(value, row["pair_id"]) for value, row in zip(values, decoy_scores)]],
                    key=lambda item: (-item[0], item[1]),
                )
                rank = next(index for index, (_, identifier) in enumerate(ranked, start=1) if identifier == target_id)
                result.update({
                    "evaluation_status": "evaluable_descriptive_only",
                    "target_score": target_value,
                    "decoy_score_median": median(values),
                    "decoy_score_min": min(values),
                    "decoy_score_max": max(values),
                    "target_minus_decoy_median": round(target_value - median(values), 6),
                    "target_rank_among_six": rank,
                    "descriptive_within_set_rank_fraction": rank / 6,
                    "target_identity_fraction": _secondary(target_score_row, "same_register_identity_fraction"),
                    "target_anchor_property_similarity": _secondary(target_score_row, "anchor_same_register_property_similarity"),
                    "target_candidate_exposed_property_similarity": _secondary(
                        target_score_row, "candidate_exposed_same_register_property_similarity"
                    ),
                })
                evaluable_labels.add(str(feasibility["target_validation_label"]))
        comparisons.append(result)

    summary = {
        "primary_endpoint": PRIMARY_SCORE,
        "target_count": len(comparisons),
        "evaluable_target_count": sum(
            row["evaluation_status"] == "evaluable_descriptive_only" for row in comparisons
        ),
        "independent_evaluable_system_count": len(evaluable_labels),
        "global_inference_status": (
            "not_computed_no_pooled_p_value" if len(evaluable_labels) >= 3
            else "insufficient_independent_systems"
        ),
        "global_p_value": "",
        "interpretation": (
            "Within-set ranks are descriptive only. No pooled p-value is calculated "
            "because targets may share evidence systems and strict decoy feasibility is limited."
        ),
    }
    return comparisons, summary


def main() -> None:
    comparisons, summary = evaluate_targets(
        read_csv(OUT / "register_aware_pair_scores.csv"),
        read_csv(BENCHMARK / "matched_decoy_sets.csv"),
        read_csv(BENCHMARK / "target_feasibility.csv"),
    )
    write_csv(OUT / "decoy_score_comparison.csv", comparisons)
    write_csv(OUT / "decoy_evaluation_summary.csv", [summary])
    print(
        f"wrote {len(comparisons)} target comparisons; "
        f"evaluable={summary['evaluable_target_count']}"
    )


if __name__ == "__main__":
    main()
