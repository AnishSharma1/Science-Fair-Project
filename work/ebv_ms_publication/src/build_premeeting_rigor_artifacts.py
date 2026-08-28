"""Build auditable HLA-II register and matched-decoy preparation artifacts.

The resulting files organize prediction-derived hypotheses and analysis-design
covariates. They do not provide evidence of peptide presentation, TCR binding,
cross-reactivity, activation, or multiple-sclerosis mechanism.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path

from premeeting_rigor import (
    binding_rank_bin,
    eligible_decoys,
    enumerate_core_windows,
    iedb_submission_segments,
    map_prediction_rows,
    natural_flank_submission_segment,
    ordered_decoys,
    parse_iedb_mhcii_tsv,
)


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
REGISTER_OUT = PROC / "register_sensitivity"
DECOY_OUT = PROC / "matched_decoys"
HYGIENE_OUT = PROC / "validation_hygiene"
IEDB_ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
IEDB_METHOD = "recommended_binding"
ALLELE = "HLA-DRB1*15:01"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def fetch_iedb_predictions(candidates: list[dict[str, str]]) -> str:
    fasta = "\n".join(
        f">{row['submission_id']}\n{row['peptide']}" for row in candidates
    )
    body = urllib.parse.urlencode({
        "method": IEDB_METHOD,
        "sequence_text": fasta,
        "allele": ALLELE,
        "length": "asis",
    }).encode("utf-8")
    request = urllib.request.Request(
        IEDB_ENDPOINT,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8")


def predicted_core_positions(peptide: str, core: str) -> str:
    starts = [
        str(index + 1)
        for index in range(max(0, len(peptide) - len(core) + 1))
        if peptide[index:index + len(core)] == core
    ]
    return ";".join(starts)


def make_prediction_rows(manifest: list[dict[str, str]], mapped: list[dict[str, str]]) -> list[dict[str, object]]:
    prediction_by_candidate: dict[str, list[dict[str, str]]] = {}
    for row in mapped:
        prediction_by_candidate.setdefault(row["candidate_id"], []).append(row)
    rows = []
    for candidate in manifest:
        base = {
            "candidate_id": candidate["candidate_id"],
            "arm": candidate["arm"],
            "evidence_tier": candidate["evidence_tier"],
            "peptide": candidate["peptide"],
            "peptide_length": candidate["peptide_length"],
            "hla": candidate["hla"],
            "prediction_method_requested": IEDB_METHOD,
            "prediction_endpoint": IEDB_ENDPOINT,
            "prediction_retrieval_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "raw_response_path": str(REGISTER_OUT / "iedb_mhcii_drb1501_raw.tsv"),
            "interpretation": "Computational binding/register hypothesis only; not experimental presentation evidence.",
        }
        predictions = prediction_by_candidate.get(candidate["candidate_id"], [])
        if not predictions:
            length = int(candidate["peptide_length"])
            if length >= 11:
                raise ValueError(f"Missing IEDB response rows for eligible candidate {candidate['candidate_id']}")
            rows.append({
                **base,
                "prediction_status": "not_submitted_length_lt_11",
                "submission_strategy": "requires verified natural flanking residues before MHC-II prediction",
                "submission_segment_count": 0,
                "iedb_seq_num": "",
                "predicted_core_peptide": "",
                "predicted_core_start_positions_1_based": "",
                "prediction_input_peptide": "",
                "predicted_core_fully_contained_in_manifest_peptide": False,
                "predicted_ic50_nM": "",
                "predicted_percentile_rank": "",
                "binding_rank_bin": "not_available",
            })
            continue
        prediction = min(predictions, key=lambda row: float(row["rank"]))
        rank = float(prediction["rank"])
        rows.append({
            **base,
            "prediction_status": "predicted",
            "submission_strategy": prediction["submission_strategy"],
            "submission_segment_count": len(predictions),
            "iedb_seq_num": prediction["seq_num"],
            "predicted_core_peptide": prediction["core_peptide"],
            "predicted_core_start_positions_1_based": predicted_core_positions(candidate["peptide"], prediction["core_peptide"]),
            "prediction_input_peptide": prediction["peptide"],
            "predicted_core_fully_contained_in_manifest_peptide": bool(predicted_core_positions(candidate["peptide"], prediction["core_peptide"])),
            "predicted_ic50_nM": prediction.get("ic50", ""),
            "predicted_percentile_rank": prediction["rank"],
            "binding_rank_bin": binding_rank_bin(rank),
        })
    return rows


def make_window_rows(predictions: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for prediction in predictions:
        peptide = str(prediction["peptide"])
        for window in enumerate_core_windows(peptide):
            rows.append({
                "candidate_id": prediction["candidate_id"],
                "peptide": peptide,
                "window_start_1_based": window["start"],
                "window_end_1_based": window["end"],
                "core_peptide": window["core_peptide"],
                "matches_iedb_predicted_core": str(window["core_peptide"]) == prediction["predicted_core_peptide"],
                "register_status": "candidate window; requires sensitivity analysis and biological review",
            })
    return rows


def make_pair_records(
    shortlist: list[dict[str, str]], annotations: list[dict[str, str]], predictions: list[dict[str, object]], manifest: list[dict[str, str]]
) -> list[dict[str, object]]:
    prediction_by_candidate = {str(row["candidate_id"]): row for row in predictions}
    manifest_by_candidate = {row["candidate_id"]: row for row in manifest}
    annotation_by_pair = {
        (row["ebv_candidate_id"], row["human_candidate_id"]): row for row in annotations
    }
    records = []
    for row in shortlist:
        ebv_id, human_id = row["ebv_candidate_id"], row["human_candidate_id"]
        ebv_prediction, human_prediction = prediction_by_candidate[ebv_id], prediction_by_candidate[human_id]
        if ebv_prediction["prediction_status"] != "predicted" or human_prediction["prediction_status"] != "predicted":
            continue
        annotation = annotation_by_pair.get((ebv_id, human_id), {})
        records.append({
            "pair_id": f"{ebv_id}::{human_id}",
            "ebv_candidate_id": ebv_id,
            "human_candidate_id": human_id,
            "ebv_peptide": manifest_by_candidate[ebv_id]["peptide"],
            "human_peptide": manifest_by_candidate[human_id]["peptide"],
            "ebv_plddt": float(row["ebv_peptide_mean_plddt"]),
            "human_plddt": float(row["human_peptide_mean_plddt"]),
            "ebv_binding_rank": float(ebv_prediction["predicted_percentile_rank"]),
            "human_binding_rank": float(human_prediction["predicted_percentile_rank"]),
            "pair_validation": annotation.get("pair_validation", "background"),
        })
    return records


def make_decoy_rows(records: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    background = [record for record in records if record["pair_validation"] == "background"]
    targets = [record for record in records if record["pair_validation"] != "background"]
    rows = []
    feasibility = []
    for target in targets:
        eligible, available = eligible_decoys(target, background, limit=5)
        feasibility.append({
            "target_pair_id": target["pair_id"],
            "target_validation_label": target["pair_validation"],
            "eligible_decoy_count": available,
            "selected_decoy_count": len(eligible),
            "target_decoy_count": 5,
            "readiness_status": "ready" if available >= 5 else "partial" if available else "no eligible decoy",
            "matching_rule": "Both peptide lengths within one residue and zero predicted-binding-rank-bin mismatches.",
            "interpretation": "Feasibility screen only; mentor approval remains required before analysis.",
        })
        for ordinal, decoy in enumerate(eligible, start=1):
            candidate = next(record for record in background if record["pair_id"] == decoy["pair_id"])
            rows.append({
                "target_pair_id": target["pair_id"],
                "target_validation_label": target["pair_validation"],
                "decoy_ordinal": ordinal,
                "decoy_pair_id": decoy["pair_id"],
                "decoy_ebv_candidate_id": candidate["ebv_candidate_id"],
                "decoy_human_candidate_id": candidate["human_candidate_id"],
                **decoy,
                "selection_boundary": "Ordered only by length, amino-acid composition, pLDDT, and predicted-binding-rank bins; never by pMHC priority score.",
                "interpretation": "Candidate decoy for mentor-reviewed analysis design; not a validated biological negative.",
            })
    return rows, feasibility


def make_cluster_rows(panel: list[dict[str, str]], annotations: list[dict[str, str]]) -> list[dict[str, object]]:
    panel_by_group: dict[str, list[dict[str, str]]] = {}
    for row in panel:
        panel_by_group.setdefault(row["validation_group"], []).append(row)
    rows = []
    for group, entries in sorted(panel_by_group.items()):
        candidate_ids = {entry["candidate_id"] for entry in entries}
        shortlist_pair_count = sum(
            1
            for row in annotations
            if row["ebv_candidate_id"] in candidate_ids or row["human_candidate_id"] in candidate_ids
        )
        if group == "classic_BALF5_MBP_structural_positive":
            role = "single literature-established calibration system"
            independence = "not independent; overlapping BALF5--MBP family records collapse to one system"
        else:
            role = "source/context overlay"
            independence = "not a direct cross-reactivity positive-pair validation set"
        rows.append({
            "validation_group": group,
            "panel_entry_count": len(entries),
            "unique_candidate_count": len(candidate_ids),
            "shortlist_pair_records_touched": shortlist_pair_count,
            "source_labels": " | ".join(sorted({entry["source_label"] for entry in entries})),
            "axes": " | ".join(sorted({entry["axis"] for entry in entries})),
            "evidence_role": role,
            "independence_interpretation": independence,
            "claim_boundary": "Literature-context/calibration bookkeeping only; not evidence of shared TCR recognition, activation, affinity, or disease mechanism.",
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--response-file",
        type=Path,
        help="Use a previously saved IEDB response instead of contacting the API.",
    )
    args = parser.parse_args()

    manifest = read_csv(PROC / "pmhc_candidate_manifest.csv")
    flanks = {row["candidate_id"]: row for row in read_csv(ROOT / "raw" / "iedb_natural_flank_extensions.csv")}
    submission_segments = []
    for candidate in manifest:
        segments = iedb_submission_segments(candidate)
        if not segments:
            flank = flanks.get(candidate["candidate_id"])
            if flank is None:
                raise ValueError(f"Missing verified natural flanks for short peptide {candidate['candidate_id']}")
            segments = [natural_flank_submission_segment(candidate, flank)]
        submission_segments.extend(segments)
    if args.response_file:
        raw_text = args.response_file.read_text(encoding="utf-8")
    else:
        raw_text = fetch_iedb_predictions(submission_segments)
    response_rows = parse_iedb_mhcii_tsv(raw_text)
    mapped = map_prediction_rows(submission_segments, response_rows)

    REGISTER_OUT.mkdir(parents=True, exist_ok=True)
    raw_path = REGISTER_OUT / "iedb_mhcii_drb1501_raw.tsv"
    raw_path.write_text(raw_text, encoding="utf-8")
    predictions = make_prediction_rows(manifest, mapped)
    prediction_fields = list(predictions[0])
    write_csv(REGISTER_OUT / "register_prediction_summary.csv", predictions, prediction_fields)
    window_rows = make_window_rows(predictions)
    write_csv(REGISTER_OUT / "register_window_catalog.csv", window_rows, list(window_rows[0]))

    annotation_path = PROC / "external_validation_benchmark" / "external_validation_pair_annotations.csv"
    if not annotation_path.exists():
        raise FileNotFoundError("Run src/run_external_validation_benchmark.py before building decoy readiness artifacts")
    shortlist = read_csv(PROC / "fullscreen_tier1_ebv_myelin_shortlist.csv")
    annotations = read_csv(annotation_path)
    records = make_pair_records(shortlist, annotations, predictions, manifest)
    decoy_rows, feasibility_rows = make_decoy_rows(records)
    if not feasibility_rows:
        raise ValueError("No validation-labeled pairs with prediction data were available for decoy preparation")
    DECOY_OUT.mkdir(parents=True, exist_ok=True)
    if decoy_rows:
        write_csv(DECOY_OUT / "decoy_readiness.csv", decoy_rows, list(decoy_rows[0]))
    else:
        write_csv(
            DECOY_OUT / "decoy_readiness.csv",
            [],
            ["target_pair_id", "target_validation_label", "decoy_ordinal", "decoy_pair_id"],
        )
    write_csv(DECOY_OUT / "decoy_feasibility_summary.csv", feasibility_rows, list(feasibility_rows[0]))

    panel = read_csv(PROC / "external_validation_panel.csv")
    cluster_rows = make_cluster_rows(panel, annotations)
    write_csv(HYGIENE_OUT / "validation_evidence_clusters.csv", cluster_rows, list(cluster_rows[0]))
    with (HYGIENE_OUT / "README.md").open("w", encoding="utf-8") as handle:
        handle.write("# Cluster-aware validation evidence inventory\n\n")
        handle.write("This inventory collapses overlapping literature annotations into evidence systems. It is not a cross-reactivity validation analysis.\n\n")
        for row in cluster_rows:
            handle.write(f"- `{row['validation_group']}`: {row['evidence_role']}; {row['independence_interpretation']}.\n")

    print(f"Wrote {REGISTER_OUT}")
    print(f"Wrote {DECOY_OUT / 'decoy_readiness.csv'}")
    print(f"Wrote {HYGIENE_OUT}")


if __name__ == "__main__":
    main()
