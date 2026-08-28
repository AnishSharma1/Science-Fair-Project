"""Inventory AlphaFold Server downloads without moving, retrying, or scoring them."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "expanded_background"
MANIFEST = OUT / "alphafold_server_seed_manifest.csv"
RESULTS_FOLDER = ROOT / "Alphafold3_pMHCs" / "json3folds"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def artifact_counts(job_directory: Path) -> tuple[int, int, int, int]:
    return (
        len(list(job_directory.glob("*_model_*.cif"))),
        len(list(job_directory.glob("*_summary_confidences_*.json"))),
        len(list(job_directory.glob("*_full_data_*.json"))),
        len(list(job_directory.glob("*_job_request.json"))),
    )


def inventory_expected_seed_jobs(
    results_folder: Path, manifest_rows: list[dict[str, str]]
) -> list[dict[str, object]]:
    """List every expected seed-3 job, including unavailable results as no-retry rows."""
    directories = {
        directory.name.lower(): directory
        for directory in results_folder.iterdir()
        if directory.is_dir()
    }
    rows: list[dict[str, object]] = []
    for manifest in manifest_rows:
        expected_name = manifest["job_name"].lower()
        job_directory = directories.get(expected_name)
        if job_directory is None:
            cifs = summaries = full_data = requests = 0
            status = "not_downloaded_or_failed_no_retry"
            interpretation = "No model result; exclude from structural summaries without retrying."
            directory_name = ""
        else:
            cifs, summaries, full_data, requests = artifact_counts(job_directory)
            directory_name = job_directory.name
            if (cifs, summaries, full_data, requests) == (5, 5, 5, 1):
                status = "complete_five_sample_result"
                interpretation = "Eligible for later descriptive structural extraction."
            else:
                status = "partial_download_exclude_from_summary"
                interpretation = "Incomplete result; exclude from structural summaries without retrying."
        rows.append({
            "batch_folder": results_folder.name,
            "batch_file": manifest["batch_file"],
            "seed_label": manifest["seed_label"],
            "model_seed": manifest["model_seed"],
            "candidate_id": manifest["candidate_id"],
            "expected_job_name": manifest["job_name"],
            "downloaded_job_directory": directory_name,
            "model_cif_count": cifs,
            "summary_confidences_count": summaries,
            "full_data_count": full_data,
            "job_request_count": requests,
            "completeness_status": status,
            "availability_interpretation": interpretation,
            "claim_boundary": "Availability is a technical run-status field, not biological evidence.",
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    manifest_rows = [row for row in read_csv(MANIFEST) if row["seed_label"] == "03"]
    inventory = inventory_expected_seed_jobs(RESULTS_FOLDER, manifest_rows)
    write_csv(OUT / "alphafold_server_seed_03_download_inventory.csv", inventory)


if __name__ == "__main__":
    main()
