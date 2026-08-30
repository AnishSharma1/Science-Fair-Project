"""Descriptively extract pMHC QA metrics from completed AlphaFold Server jobs.

This workflow is limited to model availability, sequence/layout QA, confidence,
peptide--HLA contact proxies, and within-job pose consistency. It does not
infer presentation, TCR binding, cross-reactivity, or molecular mimicry.
"""

from __future__ import annotations

import csv
import json
import re
import argparse
from collections import OrderedDict
from pathlib import Path
from statistics import median
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AF3_ROOT = ROOT / "Alphafold3_pMHCs"
OUT = ROOT / "processed" / "alphafold_server_pmhc_descriptive_analysis"
LEGACY_MANIFEST = ROOT / "processed" / "pmhc_candidate_manifest.csv"
BACKGROUND_MANIFEST = ROOT / "processed" / "expanded_background" / "background_pmhc_candidate_manifest.csv"
BACKGROUND_PREDICTIONS = ROOT / "processed" / "expanded_background" / "background_register_prediction_summary.csv"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def request_details(request: Any) -> dict[str, object]:
    """Normalize the array-shaped server request and recover its third chain."""
    job = request[0] if isinstance(request, list) else request
    chains = [entry["proteinChain"]["sequence"] for entry in job.get("sequences", []) if "proteinChain" in entry]
    if len(chains) != 3:
        raise ValueError("Expected exactly three protein chains in the AlphaFold Server request")
    seed = (job.get("modelSeeds") or [""])[0]
    return {"request_name": str(job["name"]), "server_seed": int(seed) if str(seed).isdigit() else "", "requested_dra": chains[0], "requested_drb": chains[1], "requested_peptide": chains[2]}


def parse_mmcif(path: Path) -> dict[str, list[dict[str, object]]]:
    """Read the atom loop needed for chain sequence, coordinates, and pLDDT."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        start = lines.index("_atom_site.group_PDB")
    except ValueError as error:
        raise ValueError(f"No atom-site loop in {path}") from error
    headers, index = [], start
    while index < len(lines) and lines[index].startswith("_atom_site."):
        headers.append(lines[index].removeprefix("_atom_site."))
        index += 1
    lookup = {name: position for position, name in enumerate(headers)}
    required = {"group_PDB", "type_symbol", "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv"}
    if missing := required - set(lookup):
        raise ValueError(f"Missing atom-site fields in {path}: {sorted(missing)}")
    chains: dict[str, OrderedDict[str, dict[str, object]]] = {}
    while index < len(lines):
        line = lines[index]
        if line == "#" or (line and not line.startswith("ATOM") and not line.startswith("HETATM")):
            break
        index += 1
        fields = line.split()
        if len(fields) < len(headers) or fields[lookup["group_PDB"]] != "ATOM":
            continue
        residue_name = fields[lookup["label_comp_id"]]
        if residue_name not in AA3:
            continue
        chain = fields[lookup["label_asym_id"]]
        residue = chains.setdefault(chain, OrderedDict()).setdefault(fields[lookup["label_seq_id"]], {"aa": AA3[residue_name], "atoms": [], "bfactors": []})
        try:
            xyz = tuple(float(fields[lookup[key]]) for key in ("Cartn_x", "Cartn_y", "Cartn_z"))
            bfactor = float(fields[lookup["B_iso_or_equiv"]])
        except ValueError:
            continue
        residue["atoms"].append({"name": fields[lookup["label_atom_id"]], "element": fields[lookup["type_symbol"]], "xyz": xyz})
        residue["bfactors"].append(bfactor)
    return {chain: list(residues.values()) for chain, residues in chains.items()}


def sequence(residues: list[dict[str, object]]) -> str:
    return "".join(str(residue["aa"]) for residue in residues)


def residue_plddt(residue: dict[str, object]) -> float:
    values = [float(value) for value in residue["bfactors"]]
    return sum(values) / len(values)


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return float(np.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))))


def peptide_hla_metrics(peptide: list[dict[str, object]], hla: list[list[dict[str, object]]]) -> dict[str, object]:
    hla_atoms = [atom for chain in hla for residue in chain for atom in residue["atoms"] if atom["element"] != "H"]
    counts = []
    for residue in peptide:
        peptide_atoms = [atom for atom in residue["atoms"] if atom["element"] != "H"]
        counts.append(sum(any(distance(atom["xyz"], hla_atom["xyz"]) <= 4.0 for hla_atom in hla_atoms) for atom in peptide_atoms))
    return {"peptide_mean_plddt": round(sum(residue_plddt(residue) for residue in peptide) / len(peptide), 2), "peptide_min_plddt": round(min(residue_plddt(residue) for residue in peptide), 2), "peptide_residues_with_any_hla_contact": sum(count > 0 for count in counts), "peptide_mean_hla_contacting_heavy_atoms_per_residue": round(sum(counts) / len(counts), 2), "peptide_lowest_hla_contact_positions_1based": ";".join(str(i + 1) for i, count in enumerate(counts) if count == min(counts))}


def ca_coordinates(residues: list[dict[str, object]]) -> np.ndarray:
    coordinates = []
    for residue in residues:
        ca = next((atom["xyz"] for atom in residue["atoms"] if atom["name"] == "CA"), None)
        if ca is None:
            raise ValueError("Missing CA atom")
        coordinates.append(ca)
    return np.asarray(coordinates, dtype=float)


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    source_center, target_center = source.mean(axis=0), target.mean(axis=0)
    u, _, vt = np.linalg.svd((source - source_center).T @ (target - target_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    return rotation, target_center - source_center @ rotation


def peptide_rmsd_after_hla_fit(reference: dict[str, list[dict[str, object]]], other: dict[str, list[dict[str, object]]]) -> float:
    reference_hla = np.vstack([ca_coordinates(reference[chain][:85]) for chain in ("A", "B")])
    other_hla = np.vstack([ca_coordinates(other[chain][:85]) for chain in ("A", "B")])
    rotation, translation = kabsch(other_hla, reference_hla)
    fitted = ca_coordinates(other["C"]) @ rotation + translation
    target = ca_coordinates(reference["C"])
    return float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))


def candidate_id_from_request_name(name: str) -> str:
    upper = name.upper()
    return upper.removeprefix("EBVMS_BG_").rsplit("_S", 1)[0] if upper.startswith("EBVMS_BG_") else upper


def study_candidate_id(request: Any, candidate_metadata: dict[str, dict[str, str]]) -> str | None:
    """Return a project candidate ID without assuming an arbitrary saved job is pMHC."""
    job = request[0] if isinstance(request, list) else request
    candidate_id = candidate_id_from_request_name(str(job.get("name", "")))
    return candidate_id if candidate_id in candidate_metadata else None


def sample_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.json$", path.name)
    if match is None:
        raise ValueError(f"No sample index in {path.name}")
    return int(match.group(1))


def analyze_complete_job(job_directory: Path, source_folder: str, candidate_metadata: dict[str, dict[str, str]]) -> dict[str, list[dict[str, object]] | dict[str, object]]:
    request_path = next(job_directory.glob("*_job_request.json"))
    details = request_details(json.loads(request_path.read_text(encoding="utf-8")))
    candidate_id = candidate_id_from_request_name(str(details["request_name"]))
    cohort = candidate_metadata.get(candidate_id, {}).get("af3_cohort", "excluded_unmapped_decoy_or_other")
    sample_rows, parsed_models = [], {}
    for summary_path in sorted(job_directory.glob("*_summary_confidences_*.json"), key=sample_index):
        index = sample_index(summary_path)
        cif_matches = list(job_directory.glob(f"*_model_{index}.cif"))
        if len(cif_matches) != 1:
            continue
        model = parse_mmcif(cif_matches[0])
        observed = sequence(model.get("C", []))
        status = "pass_exact_three_chain_peptide_match" if set(model) == {"A", "B", "C"} and observed == details["requested_peptide"] else "fail_peptide_sequence_mismatch"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        metrics = peptide_hla_metrics(model["C"], [model["A"], model["B"]]) if status.startswith("pass") else {}
        sample_rows.append({"source_folder": source_folder, "job_directory": job_directory.name, "request_name": details["request_name"], "candidate_id": candidate_id, "af3_cohort": cohort, "server_seed": details["server_seed"], "sample_index": index, "sequence_layout_status": status, "requested_peptide": details["requested_peptide"], "observed_peptide": observed, "ranking_score": summary.get("ranking_score", ""), "iptm": summary.get("iptm", ""), "ptm": summary.get("ptm", ""), "fraction_disordered": summary.get("fraction_disordered", ""), "has_clash": summary.get("has_clash", ""), **metrics})
        if status.startswith("pass"):
            parsed_models[index] = model
    if not sample_rows:
        raise ValueError(f"No CIF/summary pairs found in {job_directory}")
    selected = max(sample_rows, key=lambda row: float(row["ranking_score"]))
    selected_index = int(selected["sample_index"])
    rmsds = [peptide_rmsd_after_hla_fit(parsed_models[selected_index], model) for index, model in parsed_models.items() if index != selected_index] if selected_index in parsed_models else []
    statuses = {row["sequence_layout_status"] for row in sample_rows}
    job_row = {"source_folder": source_folder, "job_directory": job_directory.name, "request_name": details["request_name"], "candidate_id": candidate_id, "af3_cohort": cohort, "server_seed": details["server_seed"], "requested_peptide": details["requested_peptide"], "n_cif_summary_sample_pairs": len(sample_rows), "sequence_layout_status": "pass_exact_three_chain_peptide_match" if statuses == {"pass_exact_three_chain_peptide_match"} else "fail_peptide_sequence_mismatch", "selected_model_index": selected_index, "selected_model_ranking_score": selected["ranking_score"], "selected_model_iptm": selected["iptm"], "selected_model_ptm": selected["ptm"], "selected_model_has_clash": selected["has_clash"], "selected_model_peptide_mean_plddt": selected.get("peptide_mean_plddt", ""), "selected_model_peptide_min_plddt": selected.get("peptide_min_plddt", ""), "selected_model_peptide_residues_with_any_hla_contact": selected.get("peptide_residues_with_any_hla_contact", ""), "selected_model_peptide_mean_hla_contacting_heavy_atoms_per_residue": selected.get("peptide_mean_hla_contacting_heavy_atoms_per_residue", ""), "sample_peptide_mean_plddt_median": round(median(float(row["peptide_mean_plddt"]) for row in sample_rows if row.get("peptide_mean_plddt", "") != ""), 2) if parsed_models else "", "selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_median_A": round(median(rmsds), 3) if rmsds else "", "selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_max_A": round(max(rmsds), 3) if rmsds else "", "interpretation": "Descriptive pMHC model QA and within-job pose consistency only; not a biological ranking or binding conclusion."}
    return {"samples": sample_rows, "job": job_row}


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def summarize_cohorts(job_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Summarize completed model cohorts without comparing peptide biology."""
    summaries = []
    for cohort in sorted({str(row["af3_cohort"]) for row in job_rows}):
        rows = [row for row in job_rows if row["af3_cohort"] == cohort]
        summaries.append({
            "af3_cohort": cohort,
            "completed_jobs": len(rows),
            "selected_model_iptm_median": round(median(float(row["selected_model_iptm"]) for row in rows), 3),
            "selected_model_peptide_mean_plddt_median": round(median(float(row["selected_model_peptide_mean_plddt"]) for row in rows), 2),
            "within_job_peptide_pose_rmsd_median_A": round(median(float(row["selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_median_A"]) for row in rows), 3),
            "interpretation": "Descriptive cohort summary only; no between-peptide biological ranking.",
        })
    return summaries


def background_descriptive_rows(job_rows: list[dict[str, object]], metadata: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    """Join background prediction context without turning it into a decoy score."""
    rows = []
    for job in sorted((row for row in job_rows if row["af3_cohort"] == "new_human_background_pmhc"), key=lambda row: str(row["candidate_id"])):
        source = metadata[str(job["candidate_id"])]
        rows.append({
            "candidate_id": job["candidate_id"],
            "peptide": source.get("peptide", ""),
            "source_antigen": source.get("source_antigen", ""),
            "iedb_predicted_binding_rank_bin": source.get("binding_rank_bin", ""),
            "iedb_predicted_core": source.get("predicted_core_peptide", ""),
            "selected_model_iptm": job["selected_model_iptm"],
            "selected_model_peptide_mean_plddt": job["selected_model_peptide_mean_plddt"],
            "within_job_peptide_pose_rmsd_median_A": job["selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_median_A"],
            "interpretation": "Descriptive pMHC context only; not a decoy decision or biological ranking.",
        })
    return rows


def build_candidate_metadata() -> dict[str, dict[str, str]]:
    metadata = {row["candidate_id"]: {**row, "af3_cohort": "legacy_candidate_pmhc"} for row in read_csv(LEGACY_MANIFEST)}
    predictions = {row["candidate_id"]: row for row in read_csv(BACKGROUND_PREDICTIONS)}
    for row in read_csv(BACKGROUND_MANIFEST):
        metadata[row["candidate_id"]] = {**row, **predictions.get(row["candidate_id"], {}), "af3_cohort": "new_human_background_pmhc"}
    return metadata


def write_readme(job_rows: list[dict[str, object]]) -> None:
    cohorts: dict[str, int] = {}
    for row in job_rows:
        cohorts[str(row["af3_cohort"])] = cohorts.get(str(row["af3_cohort"]), 0) + 1
    lines = ["# AlphaFold Server pMHC descriptive analysis", "", "This folder extracts model availability, exact chain/peptide sequence QA, server confidence fields, peptide confidence, peptide--HLA contact proxies, and within-job pose consistency from completed pMHC downloads.", "", "It does **not** perform docking or screening, generate a biological candidate ranking, or infer presentation, shared-TCR binding, cross-reactivity, molecular mimicry, or MS mechanism. The server ranking score selects a representative model only within its own five-sample job.", "", "## Completed model groups", ""]
    lines.extend(f"- {cohort}: **{count}** completed jobs" for cohort, count in sorted(cohorts.items()))
    lines.extend(["", "`af3_pmhc_job_summary.csv` has one row per completed job; `af3_pmhc_sample_metrics.csv` has one row per server sample; and `af3_pmhc_cohort_summary.csv` reports descriptive cohort medians only. The 18 missing seed-03 background jobs remain in `../expanded_background/alphafold_server_seed_03_download_inventory.csv` and are neither imputed nor retried.", ""])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def merge_shards() -> None:
    """Combine completed source-folder shards without re-reading model files."""
    job_paths = sorted(OUT.glob("af3_pmhc_job_summary_json*_*[0-9].csv"))
    sample_paths = sorted(OUT.glob("af3_pmhc_sample_metrics_json*_*[0-9].csv"))
    jobs = [row for path in job_paths for row in read_csv(path)]
    samples = [row for path in sample_paths for row in read_csv(path)]
    if not jobs or not samples:
        raise ValueError("No analysis shards found to merge")
    write_csv(OUT / "af3_pmhc_job_summary.csv", jobs)
    write_csv(OUT / "af3_pmhc_sample_metrics.csv", samples)
    write_csv(OUT / "af3_pmhc_cohort_summary.csv", summarize_cohorts(jobs))
    write_csv(OUT / "new_background_af3_descriptive_summary.csv", background_descriptive_rows(jobs, build_candidate_metadata()))
    write_readme(jobs)
    print(f"Merged {len(jobs)} completed-job rows and {len(samples)} sample rows")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-folder", action="append", choices=("json1folds", "json2folds", "json3folds"))
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--merge-shards", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.merge_shards:
        merge_shards()
        return
    metadata = build_candidate_metadata()
    all_samples, all_jobs = [], []
    excluded = []
    source_folders = args.source_folder or ("json1folds", "json2folds", "json3folds")
    for source_folder in source_folders:
        job_directories = sorted(path for path in (AF3_ROOT / source_folder).iterdir() if path.is_dir() and list(path.glob("*_job_request.json")))
        if args.limit:
            job_directories = job_directories[args.start:args.start + args.limit]
        else:
            job_directories = job_directories[args.start:]
        for job_directory in job_directories:
            request_path = next(job_directory.glob("*_job_request.json"))
            request = json.loads(request_path.read_text(encoding="utf-8"))
            candidate_id = study_candidate_id(request, metadata)
            if candidate_id is None:
                raw_job = request[0] if isinstance(request, list) else request
                excluded.append({"source_folder": source_folder, "job_directory": job_directory.name, "request_name": raw_job.get("name", ""), "exclusion_reason": "Unmapped saved job; not part of either AF3 pMHC study cohort and not analyzed."})
                continue
            result = analyze_complete_job(job_directory, source_folder, metadata)
            all_samples.extend(result["samples"])
            all_jobs.append(result["job"])
    if not all_jobs:
        raise ValueError("No completed jobs selected")
    suffix = f"_{source_folders[0]}_{args.start:03d}" if len(source_folders) == 1 and args.limit else ""
    write_csv(OUT / f"af3_pmhc_sample_metrics{suffix}.csv", all_samples)
    write_csv(OUT / f"af3_pmhc_job_summary{suffix}.csv", all_jobs)
    if excluded:
        write_csv(OUT / "af3_pmhc_excluded_saved_jobs.csv", excluded)
    print(f"Wrote {len(all_jobs)} completed-job rows and {len(all_samples)} sample rows to {OUT}")


if __name__ == "__main__":
    main()
