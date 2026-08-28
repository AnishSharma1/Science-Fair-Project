"""Retrieve auditable HLA-DRB1*15:01 prediction hypotheses for comparators."""

from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from typing import Any

from build_premeeting_rigor_artifacts import fetch_iedb_predictions, predicted_core_positions
from premeeting_rigor import binding_rank_bin, map_prediction_rows, parse_iedb_mhcii_tsv


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "expanded_background"
MANIFEST = OUT / "background_pmhc_candidate_manifest.csv"
SUBMISSIONS = OUT / "background_iedb_submission_manifest.csv"
RAW = OUT / "iedb_mhcii_drb1501_raw.tsv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty prediction table")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def build_prediction_rows(
    candidates: list[dict[str, str]],
    mapped_rows: list[dict[str, str]],
    raw_response_path: str,
    retrieved_utc: str,
) -> list[dict[str, object]]:
    """Create one source-traceable top-core hypothesis per direct comparator."""
    by_candidate: dict[str, list[dict[str, str]]] = {}
    for row in mapped_rows:
        by_candidate.setdefault(row["candidate_id"], []).append(row)
    output: list[dict[str, object]] = []
    for candidate in candidates:
        predictions = by_candidate.get(candidate["candidate_id"], [])
        if not predictions:
            raise ValueError(f"missing IEDB prediction for {candidate['candidate_id']}")
        prediction = min(predictions, key=lambda row: float(row["rank"]))
        core = prediction["core_peptide"]
        starts = predicted_core_positions(candidate["peptide"], core)
        if not starts:
            raise ValueError(f"IEDB core is not contained in comparator peptide: {candidate['candidate_id']}")
        rank = float(prediction["rank"])
        output.append({
            "candidate_id": candidate["candidate_id"],
            "arm": candidate["arm"],
            "evidence_tier": candidate["evidence_tier"],
            "peptide": candidate["peptide"],
            "peptide_length": candidate["peptide_length"],
            "hla": candidate["hla"],
            "prediction_method_requested": "recommended_binding",
            "prediction_endpoint": "https://tools-cluster-interface.iedb.org/tools_api/mhcii/",
            "prediction_retrieval_utc": retrieved_utc,
            "raw_response_path": raw_response_path,
            "interpretation": "Computational binding/register hypothesis only; not experimental presentation evidence.",
            "prediction_status": "predicted",
            "submission_strategy": prediction["submission_strategy"],
            "submission_segment_count": len(predictions),
            "iedb_seq_num": prediction["seq_num"],
            "predicted_core_peptide": core,
            "predicted_core_start_positions_1_based": starts,
            "prediction_input_peptide": prediction["peptide"],
            "predicted_core_fully_contained_in_manifest_peptide": True,
            "predicted_ic50_nM": prediction.get("ic50", ""),
            "predicted_percentile_rank": prediction["rank"],
            "binding_rank_bin": binding_rank_bin(rank),
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--response-file", type=Path, help="Use a saved IEDB response instead of a live request.")
    args = parser.parse_args()
    candidates = read_csv(MANIFEST)
    submissions = read_csv(SUBMISSIONS)
    if len(candidates) != len(submissions):
        raise ValueError("direct comparator manifest and IEDB submission rows must be one-to-one")
    raw_text = (
        args.response_file.read_text(encoding="utf-8")
        if args.response_file else fetch_iedb_predictions(submissions)
    )
    response_rows = parse_iedb_mhcii_tsv(raw_text)
    mapped = map_prediction_rows(submissions, response_rows)
    RAW.write_text(raw_text, encoding="utf-8")
    retrieved = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    predictions = build_prediction_rows(candidates, mapped, str(RAW), retrieved)
    write_csv(OUT / "background_register_prediction_summary.csv", predictions)
    print(f"wrote {len(predictions)} IEDB prediction records to {OUT}")


if __name__ == "__main__":
    main()
