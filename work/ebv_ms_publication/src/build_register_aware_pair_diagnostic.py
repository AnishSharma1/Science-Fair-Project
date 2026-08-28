"""Build a transparent register-position diagnostic for the 32-pair shortlist.

The output does not replace the original ranking.  It shows how much of each
already-reported local alignment lands in equivalent class-II P1--P9 positions
under the IEDB top-core hypothesis and across every manifest-contained 9-mer
window.  Neither output is evidence of presentation or TCR recognition.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from register_aware_diagnostic import (
    parse_local_alignment_positions,
    same_register_alignment_count,
    window_pair_sensitivity,
)


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "register_sensitivity"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def unique_top_core_start(prediction: dict[str, str]) -> int | None:
    """Return one manifest-contained top-core position, otherwise ``None``."""
    if prediction["predicted_core_fully_contained_in_manifest_peptide"] != "True":
        return None
    positions = [item for item in prediction["predicted_core_start_positions_1_based"].split(";") if item]
    if len(positions) != 1:
        return None
    return int(positions[0])


def top_core_status(ebv_start: int | None, human_start: int | None) -> str:
    if ebv_start is not None and human_start is not None:
        return "assessable_top_core_hypothesis"
    return "not_assessable_top_core_requires_flank_or_register_resolution"


def main() -> None:
    shortlist = read_csv(PROC / "fullscreen_tier1_ebv_myelin_shortlist.csv")
    geometry = {
        (row["ebv_candidate_id"], row["human_candidate_id"]): row
        for row in read_csv(PROC / "colabfold_tier1_ebv_myelin_geometry_matrix.csv")
    }
    manifest = {
        row["candidate_id"]: row
        for row in read_csv(PROC / "pmhc_candidate_manifest.csv")
    }
    predictions = {
        row["candidate_id"]: row
        for row in read_csv(OUT / "register_prediction_summary.csv")
    }

    pair_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    for rank, shortlist_row in enumerate(shortlist, start=1):
        ebv_id = shortlist_row["ebv_candidate_id"]
        human_id = shortlist_row["human_candidate_id"]
        pair_id = f"{ebv_id}::{human_id}"
        ebv_peptide = manifest[ebv_id]["peptide"]
        human_peptide = manifest[human_id]["peptide"]
        geometry_row = geometry[(ebv_id, human_id)]
        alignment = parse_local_alignment_positions(geometry_row["aligned_positions_ebv_to_human"])
        ebv_prediction = predictions[ebv_id]
        human_prediction = predictions[human_id]
        ebv_top_start = unique_top_core_start(ebv_prediction)
        human_top_start = unique_top_core_start(human_prediction)
        status = top_core_status(ebv_top_start, human_top_start)
        top_count: int | str = ""
        top_fraction: float | str = ""
        if status == "assessable_top_core_hypothesis":
            top_count = same_register_alignment_count(alignment, ebv_top_start, human_top_start)
            top_fraction = round(top_count / len(alignment), 6) if alignment else 0.0

        sensitivity = window_pair_sensitivity(ebv_peptide, human_peptide, alignment)
        counts = [row["same_register_alignment_count"] for row in sensitivity]
        max_count = max(counts)
        maximizing = [
            f"EBV:{row['ebv_core_start_1_based']};Human:{row['human_core_start_1_based']}"
            for row in sensitivity
            if row["same_register_alignment_count"] == max_count
        ]
        maximizing_display = (
            "none; no same-register local alignment under any manifest-contained core combination"
            if max_count == 0
            else " | ".join(maximizing)
        )
        count_distribution = ";".join(
            f"{count}:{frequency}"
            for count, frequency in sorted(Counter(counts).items())
        )
        for row in sensitivity:
            sensitivity_rows.append({
                "pair_id": pair_id,
                "original_shortlist_rank": rank,
                "ebv_candidate_id": ebv_id,
                "human_candidate_id": human_id,
                "ebv_core_start_1_based": row["ebv_core_start_1_based"],
                "ebv_core_peptide": ebv_peptide[row["ebv_core_start_1_based"] - 1:row["ebv_core_start_1_based"] + 8],
                "human_core_start_1_based": row["human_core_start_1_based"],
                "human_core_peptide": human_peptide[row["human_core_start_1_based"] - 1:row["human_core_start_1_based"] + 8],
                "same_register_local_alignment_count": row["same_register_alignment_count"],
                "interpretation": "All manifest-contained 9-mer window combinations retained for sensitivity; not a predictor-selected register or a re-ranking score.",
            })
        pair_rows.append({
            "pair_id": pair_id,
            "original_shortlist_rank": rank,
            "original_locally_aligned_residues": len(alignment),
            "ebv_candidate_id": ebv_id,
            "human_candidate_id": human_id,
            "ebv_peptide": ebv_peptide,
            "human_peptide": human_peptide,
            "ebv_iedb_top_core": ebv_prediction["predicted_core_peptide"],
            "ebv_iedb_top_core_start_1_based": ebv_prediction["predicted_core_start_positions_1_based"],
            "human_iedb_top_core": human_prediction["predicted_core_peptide"],
            "human_iedb_top_core_start_1_based": human_prediction["predicted_core_start_positions_1_based"],
            "top_core_register_status": status,
            "top_core_same_register_local_alignment_count": top_count,
            "top_core_same_register_alignment_fraction": top_fraction,
            "all_manifest_window_pair_count": len(sensitivity),
            "all_manifest_same_register_alignment_min": min(counts),
            "all_manifest_same_register_alignment_max": max_count,
            "all_manifest_same_register_alignment_count_distribution": count_distribution,
            "maximizing_window_start_pairs": maximizing_display,
            "interpretation": "Descriptive register-position diagnostic. Retain the full sensitivity table; do not use a favorable window combination to strengthen a claim post hoc.",
        })

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(
        OUT / "register_aware_shortlist_diagnostic.csv",
        pair_rows,
        list(pair_rows[0]),
    )
    write_csv(
        OUT / "register_window_pair_sensitivity.csv",
        sensitivity_rows,
        list(sensitivity_rows[0]),
    )
    print(f"Wrote {OUT / 'register_aware_shortlist_diagnostic.csv'} ({len(pair_rows)} pairs)")
    print(f"Wrote {OUT / 'register_window_pair_sensitivity.csv'} ({len(sensitivity_rows)} window pairs)")


if __name__ == "__main__":
    main()
