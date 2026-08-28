"""Build fixed-seed AlphaFold Server import files from the frozen pMHC batch."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "expanded_background"
INPUT = OUT / "background_pmhc_colabfold_batch.csv"
SEED_03_INVENTORY = OUT / "alphafold_server_seed_03_download_inventory.csv"
SEED_03_INCOMPLETE = OUT / "alphafold_server_seed_03_incomplete_jobs.json"
SEEDS = (("01", 104729), ("02", 104759), ("03", 104761))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_server_jobs(
    batch_rows: list[dict[str, str]], *, seed: int, seed_label: str
) -> list[dict[str, object]]:
    """Create AlphaFold Server dialect jobs, one DRA/DRB/peptide complex each."""
    jobs: list[dict[str, object]] = []
    for row in batch_rows:
        chains = row["sequence"].split(":")
        if len(chains) != 3 or not all(chains):
            raise ValueError(f"{row.get('id', '<unknown>')} must contain exactly three chains")
        identifier = row["id"]
        jobs.append({
            "name": f"ebvms_bg_{identifier}_s{seed_label}",
            "modelSeeds": [seed],
            "sequences": [
                {"proteinChain": {"sequence": chain, "count": 1}}
                for chain in chains
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        })
    return jobs


def write_json(path: Path, jobs: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(jobs, handle, indent=2)
        handle.write("\n")


def select_incomplete_jobs(
    jobs: list[dict[str, object]], inventory_rows: list[dict[str, str]]
) -> list[dict[str, object]]:
    """Exclude Seed 03 jobs with a verified complete five-sample result."""
    inventory_by_name = {row["expected_job_name"]: row for row in inventory_rows}
    job_names = {str(job["name"]) for job in jobs}
    if set(inventory_by_name) != job_names:
        raise ValueError("Seed 03 inventory must cover every frozen Seed 03 job exactly")
    return [
        job for job in jobs
        if inventory_by_name[str(job["name"])]["completeness_status"]
        != "complete_five_sample_result"
    ]


def main() -> None:
    rows = read_csv(INPUT)
    if len(rows) > 100:
        raise ValueError("AlphaFold Server imports allow at most 100 jobs per file")
    manifest: list[dict[str, object]] = []
    for seed_label, seed in SEEDS:
        filename = f"alphafold_server_seed_{seed_label}_jobs.json"
        jobs = build_server_jobs(rows, seed=seed, seed_label=seed_label)
        write_json(OUT / filename, jobs)
        manifest.extend({
            "batch_file": filename,
            "seed_label": seed_label,
            "model_seed": seed,
            "candidate_id": row["id"],
            "job_name": job["name"],
            "chain_order": "HLA-DRA;HLA-DRB1*15:01;human-background peptide",
        } for row, job in zip(rows, jobs))
        if seed_label == "03" and SEED_03_INVENTORY.exists():
            incomplete_jobs = select_incomplete_jobs(jobs, read_csv(SEED_03_INVENTORY))
            write_json(SEED_03_INCOMPLETE, incomplete_jobs)
    with (OUT / "alphafold_server_seed_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)


if __name__ == "__main__":
    main()
