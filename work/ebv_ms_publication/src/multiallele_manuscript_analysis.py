"""Reproducible multi-allele EBV--MS pMHC analysis helpers.

The module implements technical inventory, register eligibility, score-blind
control selection, fixed-seed robustness inputs, and within-allele P1--P9
geometry.  Outputs are computational pMHC descriptors only; they do not show
presentation, TCR binding, activation, cross-reactivity, or disease mechanism.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Optional, Sequence

import numpy as np

from analyze_af3_pmhc_downloads import (
    ca_coordinates,
    kabsch,
    parse_mmcif,
    peptide_hla_metrics,
    peptide_rmsd_after_hla_fit,
    request_details,
    sequence,
)
from premeeting_rigor import AMINO_ACIDS, binding_rank_bin, composition_distance
from register_aware_scoring import property_similarity


ALLELES = ("HLA-DRB1*13:03", "HLA-DRB1*03:01", "HLA-DRB1*08:01")
ALLELE_CODES = {
    "HLA-DRB1*13:03": "drb1303",
    "HLA-DRB1*03:01": "drb0301",
    "HLA-DRB1*08:01": "drb0801",
}
ANCHOR_EBV = "EBV_TCELL_950"
ANCHOR_HUMAN = "HUMAN_MYELIN_112214"
ANCHOR_POSITIONS = (1, 4, 6, 9)
EXPOSED_POSITIONS = (2, 3, 5, 7, 8)
ROBUSTNESS_SEEDS = (104729, 104759)
CLAIM_BOUNDARY = (
    "Computational pMHC geometry only; not evidence of presentation, TCR binding, "
    "activation, cross-reactivity, molecular mimicry, or MS mechanism."
)
IEDB_ENDPOINT = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
IEDB_METHOD = "recommended_binding"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: Optional[list[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fields:
        raise ValueError(f"refusing to write empty table without fields: {path}")
    fieldnames = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_request_name(name: str) -> str:
    """Collapse the observed `_2` download duplicate to its canonical request."""
    lowered = name.strip().lower()
    if lowered == "ebvms_drb1303_human_myelin_115622_2":
        return lowered[:-2]
    return lowered


def _truthy_clash(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "1.0", "true", "yes"}


def select_representative_sample(rows: Sequence[dict[str, Any]]) -> Optional[dict[str, Any]]:
    eligible = [
        row for row in rows
        if str(row.get("sequence_layout_status", "")).startswith("pass")
        and not _truthy_clash(row.get("has_clash", False))
    ]
    return max(eligible, key=lambda row: (float(row["ranking_score"]), -int(row["sample_index"]))) if eligible else None


def _sample_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.(?:json|cif)$", path.name)
    if not match:
        raise ValueError(f"cannot recover sample index from {path.name}")
    return int(match.group(1))


def _allele_for_request(name: str) -> str:
    canonical = canonical_request_name(name)
    for allele, code in ALLELE_CODES.items():
        if f"_{code}_" in f"_{canonical}_":
            return allele
    raise ValueError(f"request name lacks a supported allele code: {name}")


def discover_download_jobs(result_roots: Iterable[Path]) -> list[dict[str, Any]]:
    """Discover request-bearing folders across arbitrary AlphaFold result roots."""
    rows: list[dict[str, Any]] = []
    for result_root in sorted({Path(root).resolve() for root in result_roots}, key=str):
        if not result_root.exists():
            continue
        for directory in sorted((p for p in result_root.iterdir() if p.is_dir()), key=lambda p: p.name):
            request_files = list(directory.glob("*_job_request.json"))
            if len(request_files) != 1:
                continue
            details = request_details(json.loads(request_files[0].read_text(encoding="utf-8")))
            request_name = str(details["request_name"])
            rows.append({
                "result_root": str(result_root),
                "job_directory": directory.name,
                "job_directory_path": str(directory),
                "request_path": str(request_files[0]),
                "request_name": request_name,
                "canonical_request_name": canonical_request_name(request_name),
                "allele": _allele_for_request(request_name),
                "server_seed": details["server_seed"],
                "requested_dra": details["requested_dra"],
                "requested_drb": details["requested_drb"],
                "requested_peptide": details["requested_peptide"],
                "n_cif": len(list(directory.glob("*_model_*.cif"))),
                "n_confidence": len(list(directory.glob("*_summary_confidences_*.json"))),
                "n_full_data": len(list(directory.glob("*_full_data_*.json"))),
                "n_request": len(request_files),
            })
    return rows


def inventory_against_manifest(
    expected_rows: list[dict[str, str]], discovered: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Resolve canonical requests, retaining the known `_2` run as sensitivity."""
    by_canonical: dict[str, list[dict[str, Any]]] = {}
    for row in discovered:
        by_canonical.setdefault(str(row["canonical_request_name"]), []).append(row)
    inventory: list[dict[str, Any]] = []
    duplicates: list[dict[str, Any]] = []
    for expected in expected_rows:
        name = expected["job_name"].lower()
        matches = sorted(by_canonical.get(name, []), key=lambda row: (str(row["request_name"]).endswith("_2"), str(row["job_directory_path"])))
        primary = matches[0] if matches else None
        complete = bool(primary) and all(int(primary[field]) == required for field, required in (("n_cif", 5), ("n_confidence", 5), ("n_full_data", 5), ("n_request", 1)))
        inventory.append({
            **expected,
            "canonical_request_name": name,
            "download_status": "complete" if complete else "incomplete" if primary else "missing",
            "primary_job_directory_path": primary["job_directory_path"] if primary else "",
            "server_seed": primary["server_seed"] if primary else "",
            "n_cif": primary["n_cif"] if primary else 0,
            "n_confidence": primary["n_confidence"] if primary else 0,
            "n_full_data": primary["n_full_data"] if primary else 0,
            "n_request": primary["n_request"] if primary else 0,
            "duplicate_run_count": max(0, len(matches) - 1),
        })
        for duplicate in matches[1:]:
            duplicates.append({
                "canonical_request_name": name,
                "request_name": duplicate["request_name"],
                "job_directory_path": duplicate["job_directory_path"],
                "handling": "duplicate_run_sensitivity_only",
            })
    return inventory, duplicates


def analyze_job(
    inventory_row: dict[str, Any], expected_dra: str, expected_drb: str
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Extract five sample rows and one clash-aware job summary."""
    directory = Path(str(inventory_row["primary_job_directory_path"]))
    request_path = next(directory.glob("*_job_request.json"))
    details = request_details(json.loads(request_path.read_text(encoding="utf-8")))
    sample_rows: list[dict[str, Any]] = []
    models: dict[int, dict[str, list[dict[str, object]]]] = {}
    for summary_path in sorted(directory.glob("*_summary_confidences_*.json"), key=_sample_index):
        index = _sample_index(summary_path)
        cif_paths = list(directory.glob(f"*_model_{index}.cif"))
        full_paths = list(directory.glob(f"*_full_data_{index}.json"))
        if len(cif_paths) != 1 or len(full_paths) != 1:
            continue
        model = parse_mmcif(cif_paths[0])
        observed_dra = sequence(model.get("A", []))
        observed_drb = sequence(model.get("B", []))
        observed_peptide = sequence(model.get("C", []))
        layout_pass = set(model) == {"A", "B", "C"}
        exact = (
            layout_pass
            and observed_dra == details["requested_dra"] == expected_dra
            and observed_drb == details["requested_drb"] == expected_drb
            and observed_peptide == details["requested_peptide"] == inventory_row["peptide_sequence"]
        )
        status = "pass_exact_three_chain_sequence_match" if exact else "fail_input_model_sequence_mismatch"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        metrics = peptide_hla_metrics(model["C"], [model["A"], model["B"]]) if exact else {}
        row = {
            "allele": inventory_row["allele"],
            "candidate_id": inventory_row["candidate_id"],
            "request_name": details["request_name"],
            "canonical_request_name": inventory_row["canonical_request_name"],
            "server_seed": details["server_seed"],
            "sample_index": index,
            "canonical_model_key": f"{inventory_row['allele']}|{inventory_row['candidate_id']}|{details['server_seed']}|{index}",
            "sequence_layout_status": status,
            "requested_dra": details["requested_dra"],
            "observed_dra": observed_dra,
            "requested_drb": details["requested_drb"],
            "observed_drb": observed_drb,
            "requested_peptide": details["requested_peptide"],
            "observed_peptide": observed_peptide,
            "ranking_score": summary.get("ranking_score", ""),
            "iptm": summary.get("iptm", ""),
            "ptm": summary.get("ptm", ""),
            "has_clash": bool(_truthy_clash(summary.get("has_clash", False))),
            "cif_path": str(cif_paths[0]),
            "confidence_path": str(summary_path),
            "full_data_path": str(full_paths[0]),
            **metrics,
            "claim_boundary": CLAIM_BOUNDARY,
        }
        sample_rows.append(row)
        if exact:
            models[index] = model
    representative = select_representative_sample(sample_rows)
    valid_models = {
        int(row["sample_index"]): models[int(row["sample_index"])]
        for row in sample_rows
        if str(row["sequence_layout_status"]).startswith("pass")
        and not _truthy_clash(row["has_clash"])
        and int(row["sample_index"]) in models
    }
    rmsds: list[float] = []
    if representative is not None:
        selected_index = int(representative["sample_index"])
        rmsds = [
            peptide_rmsd_after_hla_fit(valid_models[selected_index], model)
            for index, model in valid_models.items() if index != selected_index
        ]
    job = {
        "allele": inventory_row["allele"],
        "candidate_id": inventory_row["candidate_id"],
        "request_name": details["request_name"],
        "canonical_request_name": inventory_row["canonical_request_name"],
        "server_seed": details["server_seed"],
        "technical_status": "geometry_evaluable" if representative is not None else "excluded_no_clash_free_exact_sample",
        "n_samples": len(sample_rows),
        "n_exact_samples": sum(str(row["sequence_layout_status"]).startswith("pass") for row in sample_rows),
        "n_clash_free_exact_samples": len(valid_models),
        "selected_sample_index": representative["sample_index"] if representative else "",
        "selected_ranking_score": representative["ranking_score"] if representative else "",
        "selected_iptm": representative["iptm"] if representative else "",
        "selected_ptm": representative["ptm"] if representative else "",
        "selected_has_clash": representative["has_clash"] if representative else "",
        "selected_peptide_mean_plddt": representative.get("peptide_mean_plddt", "") if representative else "",
        "selected_peptide_min_plddt": representative.get("peptide_min_plddt", "") if representative else "",
        "selected_peptide_residues_with_any_hla_contact": representative.get("peptide_residues_with_any_hla_contact", "") if representative else "",
        "selected_peptide_mean_hla_contacting_heavy_atoms_per_residue": representative.get("peptide_mean_hla_contacting_heavy_atoms_per_residue", "") if representative else "",
        "within_job_pose_rmsd_median_A": round(median(rmsds), 4) if rmsds else "",
        "within_job_pose_rmsd_max_A": round(max(rmsds), 4) if rmsds else "",
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return sample_rows, job


def _all_occurrences(sequence_text: str, core: str) -> list[int]:
    return [index + 1 for index in range(max(0, len(sequence_text) - len(core) + 1)) if sequence_text[index:index + len(core)] == core]


def register_record(
    *,
    candidate_id: str,
    modeled_peptide: str,
    prediction_input: str,
    predicted_core: str,
    percentile_rank: float,
    seq_num: int,
    original_start_in_prediction_1_based: int = 1,
) -> dict[str, Any]:
    """Classify a predicted 9-mer without allowing flank-dependent substitution."""
    input_starts = _all_occurrences(prediction_input, predicted_core)
    original_end = original_start_in_prediction_1_based + len(modeled_peptide) - 1
    contained_starts = [
        start - original_start_in_prediction_1_based + 1
        for start in input_starts
        if start >= original_start_in_prediction_1_based
        and start + len(predicted_core) - 1 <= original_end
    ]
    modeled_starts = _all_occurrences(modeled_peptide, predicted_core)
    if len(modeled_starts) > 1:
        status = "unresolved_tied_core_position"
        core_start: Any = ""
    elif len(modeled_starts) == 1 and len(contained_starts) == 1:
        status = "resolved_unique_fully_contained"
        core_start = modeled_starts[0]
    elif input_starts:
        status = "unresolved_flank_dependent"
        core_start = ""
    else:
        status = "unresolved_core_not_found"
        core_start = ""
    return {
        "candidate_id": candidate_id,
        "modeled_peptide": modeled_peptide,
        "prediction_input_peptide": prediction_input,
        "iedb_seq_num": seq_num,
        "predicted_core_peptide": predicted_core,
        "predicted_percentile_rank": percentile_rank,
        "binding_rank_bin": binding_rank_bin(float(percentile_rank)),
        "original_start_in_prediction_1_based": original_start_in_prediction_1_based,
        "predicted_core_input_start_positions_1_based": ";".join(map(str, input_starts)),
        "predicted_core_modeled_start_positions_1_based": ";".join(map(str, modeled_starts)),
        "predicted_core_fully_contained": status == "resolved_unique_fully_contained",
        "core_start_1_based": core_start,
        "register_status": status,
        "claim_boundary": "IEDB binding/register hypothesis only; not experimental presentation evidence.",
    }


def build_prediction_submissions(
    candidates: list[dict[str, Any]], flanks: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    """Build exactly one IEDB input per modeled peptide, using verified flanks."""
    submissions: list[dict[str, Any]] = []
    for seq_num, candidate in enumerate(candidates, start=1):
        candidate_id = str(candidate["candidate_id"])
        modeled = str(candidate.get("peptide_sequence", candidate.get("peptide", ""))).upper()
        if not modeled or not set(modeled).issubset(set(AMINO_ACIDS)):
            raise ValueError(f"invalid modeled peptide for {candidate_id}")
        if len(modeled) >= 11:
            prediction_input = modeled
            original_start = 1
            strategy = "direct_full_peptide"
        else:
            flank = flanks.get(candidate_id)
            if flank is None:
                raise ValueError(f"short peptide lacks a verified natural flank: {candidate_id}")
            prediction_input = str(flank["extended_sequence"]).upper()
            original_start = int(flank["original_start_in_extended_1_based"])
            original_end = int(flank["original_end_in_extended_1_based"])
            if prediction_input[original_start - 1:original_end] != modeled:
                raise ValueError(f"verified flank coordinates do not reproduce {candidate_id}")
            if not 11 <= len(prediction_input) <= 30:
                raise ValueError(f"flank input is outside the IEDB 11--30 aa range: {candidate_id}")
            strategy = "verified_natural_flank_extension"
        submissions.append({
            "seq_num": seq_num,
            "candidate_id": candidate_id,
            "modeled_peptide": modeled,
            "prediction_input_peptide": prediction_input,
            "original_start_in_prediction_1_based": original_start,
            "submission_strategy": strategy,
        })
    return submissions


def fetch_iedb_predictions(allele: str, submissions: list[dict[str, Any]]) -> str:
    if allele not in ALLELES:
        raise ValueError(f"unsupported frozen allele: {allele}")
    fasta = "\n".join(
        f">{row['candidate_id']}\n{row['prediction_input_peptide']}" for row in submissions
    )
    body = urllib.parse.urlencode({
        "method": IEDB_METHOD,
        "sequence_text": fasta,
        "allele": allele,
        "length": "asis",
    }).encode("utf-8")
    request = urllib.request.Request(
        IEDB_ENDPOINT,
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return response.read().decode("utf-8")


def prediction_records_from_tsv(
    allele: str,
    submissions: list[dict[str, Any]],
    raw_text: str,
    raw_response_path: str,
    retrieved_utc: Optional[str] = None,
) -> list[dict[str, Any]]:
    response_rows = list(csv.DictReader(io.StringIO(raw_text), delimiter="\t"))
    required = {"allele", "seq_num", "core_peptide", "peptide", "ic50", "rank"}
    if not response_rows or not required.issubset(response_rows[0]):
        raise ValueError(f"IEDB returned an incomplete table for {allele}")
    by_seq: dict[int, dict[str, str]] = {}
    for response in response_rows:
        seq_num = int(response["seq_num"])
        if seq_num in by_seq:
            raise ValueError(f"IEDB returned duplicate seq_num {seq_num} for {allele}")
        by_seq[seq_num] = response
    if len(by_seq) != len(submissions):
        raise ValueError(f"IEDB returned {len(by_seq)} rows for {len(submissions)} submissions to {allele}")
    retrieved = retrieved_utc or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    records: list[dict[str, Any]] = []
    for submission in submissions:
        seq_num = int(submission["seq_num"])
        response = by_seq.get(seq_num)
        if response is None:
            raise ValueError(f"IEDB omitted seq_num {seq_num} for {allele}")
        if response["allele"] != allele:
            raise ValueError(f"IEDB allele mismatch at seq_num {seq_num}: {response['allele']}")
        if response["peptide"] != submission["prediction_input_peptide"]:
            raise ValueError(f"IEDB peptide mismatch at seq_num {seq_num} for {allele}")
        record = register_record(
            candidate_id=str(submission["candidate_id"]),
            modeled_peptide=str(submission["modeled_peptide"]),
            prediction_input=response["peptide"],
            predicted_core=response["core_peptide"],
            percentile_rank=float(response["rank"]),
            seq_num=seq_num,
            original_start_in_prediction_1_based=int(submission["original_start_in_prediction_1_based"]),
        )
        records.append({
            "allele": allele,
            **record,
            "submission_strategy": submission["submission_strategy"],
            "predicted_ic50_nM": response["ic50"],
            "prediction_method_requested": IEDB_METHOD,
            "prediction_endpoint": IEDB_ENDPOINT,
            "prediction_retrieval_utc": retrieved,
            "raw_response_path": raw_response_path,
        })
    return records


def select_score_blind_controls(
    target_peptide: str,
    target_binding_bin: str,
    prediction_rows: list[dict[str, Any]],
    *,
    limit: int = 3,
) -> list[dict[str, Any]]:
    """Apply the frozen composition, length, numeric-ID order without AF scores."""
    if limit < 1:
        raise ValueError("control limit must be positive")
    candidates = [
        row for row in prediction_rows
        if abs(int(row["peptide_length"]) - len(target_peptide)) <= 1
        and str(row["binding_rank_bin"]) == target_binding_bin
        and str(row.get("register_status", "resolved_unique_fully_contained")) == "resolved_unique_fully_contained"
    ]
    ordered = sorted(candidates, key=lambda row: (
        composition_distance(target_peptide, str(row["peptide"])),
        int(row["peptide_length"]),
        int(row["iedb_epitope_id"]),
        str(row["candidate_id"]),
    ))
    fields = (
        "candidate_id", "iedb_epitope_id", "peptide", "peptide_length",
        "source_accession", "source_antigen_name", "source_start_1_based",
        "source_end_1_based", "binding_rank_bin", "predicted_percentile_rank",
        "predicted_core_peptide", "core_start_1_based", "register_status",
    )
    output: list[dict[str, Any]] = []
    for order, row in enumerate(ordered[:limit], start=1):
        selected = {field: row.get(field, "") for field in fields}
        selected.update({
            "selection_order": order,
            "composition_distance_to_mbp": round(composition_distance(target_peptide, str(row["peptide"])), 8),
            "length_distance_to_mbp": abs(int(row["peptide_length"]) - len(target_peptide)),
            "selection_rule": "composition_distance_then_peptide_length_then_numeric_iedb_identifier",
        })
        output.append(selected)
    return output


def build_pair_universe(
    allele: str,
    panel_rows: list[dict[str, Any]],
    registers: dict[tuple[str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    ebv = sorted((row for row in panel_rows if row["arm_group"] == "EBV"), key=lambda row: str(row["candidate_id"]))
    human = sorted((row for row in panel_rows if row["arm_group"] == "CNS/self"), key=lambda row: str(row["candidate_id"]))
    if len(ebv) != 25 or len(human) != 25:
        raise ValueError(f"panel must contain 25 EBV and 25 CNS/self peptides, got {len(ebv)} and {len(human)}")
    rows: list[dict[str, Any]] = []
    for left in ebv:
        for right in human:
            left_register = registers.get((allele, str(left["candidate_id"])), {})
            right_register = registers.get((allele, str(right["candidate_id"])), {})
            eligible = (
                left_register.get("register_status") == "resolved_unique_fully_contained"
                and right_register.get("register_status") == "resolved_unique_fully_contained"
            )
            rows.append({
                "allele": allele,
                "pair_id": f"{allele}::{left['candidate_id']}::{right['candidate_id']}",
                "ebv_candidate_id": left["candidate_id"],
                "human_candidate_id": right["candidate_id"],
                "analysis_role": "primary_anchor" if left["candidate_id"] == ANCHOR_EBV and right["candidate_id"] == ANCHOR_HUMAN else "exploratory",
                "ebv_register_status": left_register.get("register_status", "missing_prediction"),
                "human_register_status": right_register.get("register_status", "missing_prediction"),
                "register_eligible": eligible,
                "geometry_status": "pending_geometry" if eligible else "excluded_unresolved_register",
                "claim_boundary": CLAIM_BOUNDARY,
            })
    return rows


def _rmsd(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def same_register_geometry_from_coordinates(
    left_groove: np.ndarray,
    left_core: np.ndarray,
    right_groove: np.ndarray,
    right_core: np.ndarray,
) -> dict[str, float]:
    """Fit equivalent groove coordinates, then compare exact P1--P9 atoms."""
    if left_core.shape != (9, 3) or right_core.shape != (9, 3):
        raise ValueError("both peptide cores must contain exactly nine CA coordinates")
    if left_groove.shape != right_groove.shape or left_groove.ndim != 2 or left_groove.shape[1] != 3:
        raise ValueError("groove coordinate arrays must have equal N x 3 shapes")
    rotation, translation = kabsch(right_groove, left_groove)
    fitted = right_core @ rotation + translation
    anchor_indices = np.asarray([position - 1 for position in ANCHOR_POSITIONS])
    exposed_indices = np.asarray([position - 1 for position in EXPOSED_POSITIONS])
    return {
        "full_core_ca_rmsd_A": _rmsd(left_core, fitted),
        "anchor_ca_rmsd_A": _rmsd(left_core[anchor_indices], fitted[anchor_indices]),
        "exposed_ca_rmsd_A": _rmsd(left_core[exposed_indices], fitted[exposed_indices]),
    }


def same_register_geometry_from_models(
    left: dict[str, list[dict[str, object]]],
    right: dict[str, list[dict[str, object]]],
    left_core_start: int,
    right_core_start: int,
) -> dict[str, float]:
    left_groove = np.vstack([ca_coordinates(left[chain][:85]) for chain in ("A", "B")])
    right_groove = np.vstack([ca_coordinates(right[chain][:85]) for chain in ("A", "B")])
    left_core = ca_coordinates(left["C"])[left_core_start - 1:left_core_start + 8]
    right_core = ca_coordinates(right["C"])[right_core_start - 1:right_core_start + 8]
    return same_register_geometry_from_coordinates(left_groove, left_core, right_groove, right_core)


def direct_register_sequence_metrics(left_core: str, right_core: str) -> dict[str, float]:
    if len(left_core) != 9 or len(right_core) != 9:
        raise ValueError("sequence metrics require two 9-mer cores")
    def averaged_similarity(positions: Sequence[int]) -> float:
        descriptors = [property_similarity(left_core[p - 1], right_core[p - 1]) for p in positions]
        return float(np.mean([value for descriptor in descriptors for value in descriptor.values()]))
    return {
        "full_core_property_similarity_mean": round(averaged_similarity(tuple(range(1, 10))), 6),
        "anchor_property_similarity_mean": round(averaged_similarity(ANCHOR_POSITIONS), 6),
        "exposed_property_similarity_mean": round(averaged_similarity(EXPOSED_POSITIONS), 6),
        "exposed_exact_identity_count": sum(left_core[p - 1] == right_core[p - 1] for p in EXPOSED_POSITIONS),
    }


def build_robustness_jobs(
    entities: list[dict[str, Any]], seeds: Sequence[int] = ROBUSTNESS_SEEDS
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if tuple(seeds) != ROBUSTNESS_SEEDS:
        raise ValueError(f"robustness seeds are frozen as {ROBUSTNESS_SEEDS}")
    keys = [(str(row["allele"]), str(row["entity_id"])) for row in entities]
    if len(keys) != len(set(keys)):
        raise ValueError("robustness entities must be unique within allele")
    jobs: list[dict[str, Any]] = []
    manifest: list[dict[str, Any]] = []
    for entity in sorted(entities, key=lambda row: (ALLELES.index(str(row["allele"])), str(row["entity_id"]))):
        allele = str(entity["allele"])
        for seed in seeds:
            safe_id = re.sub(r"[^a-z0-9]+", "_", str(entity["entity_id"]).lower()).strip("_")
            name = f"ebvms_robust_{ALLELE_CODES[allele]}_{safe_id}_s{seed}"
            job = {
                "name": name,
                "modelSeeds": [int(seed)],
                "sequences": [
                    {"proteinChain": {"sequence": str(sequence_value), "count": 1}}
                    for sequence_value in (entity["dra_sequence"], entity["drb_sequence"], entity["peptide"])
                ],
                "dialect": "alphafoldserver",
                "version": 1,
            }
            jobs.append(job)
            manifest.append({
                "job_name": name,
                "allele": allele,
                "entity_id": entity["entity_id"],
                "entity_role": entity.get("entity_role", ""),
                "peptide": entity["peptide"],
                "server_seed": seed,
                "chain_order": f"HLA-DRA;{allele};peptide",
                "claim_boundary": CLAIM_BOUNDARY,
            })
    if len(jobs) > 30:
        raise ValueError("robustness batch exceeds the predeclared 30-job maximum")
    return jobs, manifest


def empirical_tail_fraction(anchor_value: float, control_values: Sequence[float]) -> float:
    """Exploratory lower-tail fraction with a +1 finite-control correction."""
    return (1 + sum(value <= anchor_value for value in control_values)) / (len(control_values) + 1)


def summarize_anchor_controls(anchor_values: Sequence[float], controls: dict[str, Sequence[float]]) -> dict[str, Any]:
    if not anchor_values or len(controls) != 3 or any(not values for values in controls.values()):
        raise ValueError("anchor summary requires anchor values and exactly three non-empty controls")
    anchor_median = median(anchor_values)
    control_medians = {key: median(values) for key, values in sorted(controls.items())}
    leave_one_out = {
        omitted: median([value for key, value in control_medians.items() if key != omitted])
        for omitted in control_medians
    }
    return {
        "anchor_exposed_rmsd_median_A": round(anchor_median, 6),
        "equal_weight_control_median_A": round(median(control_medians.values()), 6),
        "control_medians_A": json.dumps(control_medians, sort_keys=True),
        "leave_one_control_out_medians_A": json.dumps(leave_one_out, sort_keys=True),
        "empirical_lower_tail_fraction_exploratory": empirical_tail_fraction(anchor_median, list(control_medians.values())),
        "empirical_tail_interpretation": "Exploratory finite-control tail fraction; minimum 0.25 with three controls, not a p-value.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true", help="Run a lightweight deterministic contract check.")
    args = parser.parse_args()
    if args.self_check:
        panel = [
            *({"candidate_id": f"E{i}", "arm_group": "EBV"} for i in range(25)),
            *({"candidate_id": f"H{i}", "arm_group": "CNS/self"} for i in range(25)),
        ]
        assert len(build_pair_universe(ALLELES[0], panel, {})) == 625
        print("multi-allele analysis contracts: OK")
        return
    parser.error("use the dated orchestration script to generate project artifacts")


if __name__ == "__main__":
    main()
