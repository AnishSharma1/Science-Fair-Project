"""Analyze downloaded AlphaFold models for the held-out HLA-II controls.

The workflow is additive: it reads the frozen control package and model files,
writes a separate result package, and never changes discovery rankings.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from build_hla2_positive_control_benchmark import _hla_sequences, curated_registry
from hla2_positive_control_benchmark import (
    CLAIM_BOUNDARY,
    EXPOSED_INDICES,
    PmhcGeometry,
    build_trust_gate,
    geometry_from_mmcif,
    leave_one_system_out,
    pair_features,
    parse_mmcif_atoms,
    residue_sequence,
)


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "processed/hla2_positive_control_benchmark_2026-08-25"
DEFAULT_OUT = ROOT / "processed/hla2_positive_control_benchmark_results_2026-08-26"
V2_PACKAGE = ROOT / "processed/tcell_library_v2_2026-08-22"
OLD_MODEL_ANALYSIS = ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25"
DEFAULT_DOWNLOADS = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Downloads"
FEATURES = (
    "exposed_ca_rmsd_A",
    "exposed_sidechain_vector_rmsd_A",
    "tcr_face_physicochemical_mismatch",
    "anchor_ca_rmsd_A",
)


@dataclass(frozen=True)
class GeometrySample:
    sample_index: int
    geometry: PmhcGeometry


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] = ()
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or sorted({key for row in rows for key in row}))
    if not fieldnames:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _request_job(payload: Any) -> Mapping[str, Any]:
    job = payload[0] if isinstance(payload, list) else payload
    if not isinstance(job, Mapping):
        raise ValueError("AlphaFold request must contain one job object")
    return job


def _request_sequences(job: Mapping[str, Any]) -> list[str]:
    sequences = []
    for entry in job.get("sequences", []):
        if not isinstance(entry, Mapping) or "proteinChain" not in entry:
            raise ValueError("AlphaFold request must contain protein chains only")
        sequences.append(str(entry["proteinChain"]["sequence"]))
    if len(sequences) != 3:
        raise ValueError("AlphaFold request must contain exactly three protein chains")
    return sequences


def compare_request_to_expected(payload: Any, expected: Mapping[str, Any]) -> dict[str, Any]:
    observed = _request_job(payload)
    observed_seeds = list(observed.get("modelSeeds") or [])
    expected_seeds = list(expected.get("modelSeeds") or [])
    observed_seed = observed_seeds[0] if len(observed_seeds) == 1 else ""
    expected_seed = expected_seeds[0] if len(expected_seeds) == 1 else ""
    normalized_seed = int(observed_seed) if str(observed_seed).isdigit() else ""
    normalized_expected = int(expected_seed) if str(expected_seed).isdigit() else ""
    name_pass = str(observed.get("name", "")) == str(expected.get("name", ""))
    chain_sequences_pass = _request_sequences(observed) == _request_sequences(expected)
    seed_pass = normalized_seed != "" and normalized_seed == normalized_expected
    return {
        "request_identity_pass": name_pass and chain_sequences_pass and seed_pass,
        "request_name_pass": name_pass,
        "chain_sequences_pass": chain_sequences_pass,
        "seed_pass": seed_pass,
        "normalized_seed": normalized_seed,
        "seed_serialization": "string" if isinstance(observed_seed, str) else type(observed_seed).__name__,
    }


def _bundle_fingerprint(directory: Path) -> str:
    digest = hashlib.sha256()
    files = sorted([
        *directory.glob("*_job_request.json"),
        *directory.glob("*_model_*.cif"),
        *directory.glob("*_summary_confidences_*.json"),
        *directory.glob("*_full_data_*.json"),
    ], key=lambda path: path.name)
    for path in files:
        digest.update(path.name.rsplit("_", 3)[-1].encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    return digest.hexdigest()


def inventory_downloaded_jobs(
    manifest: Sequence[Mapping[str, Any]],
    batch_jobs: Sequence[Mapping[str, Any]],
    result_roots: Sequence[Path],
) -> list[dict[str, Any]]:
    expected_jobs = {str(job["name"]).lower(): job for job in batch_jobs}
    occurrences: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for root in sorted({Path(path).expanduser().resolve() for path in result_roots}, key=str):
        if not root.exists():
            continue
        for request_path in sorted(root.rglob("*_job_request.json"), key=str):
            directory = request_path.parent
            payload = json.loads(request_path.read_text(encoding="utf-8"))
            observed = _request_job(payload)
            name = str(observed.get("name", "")).lower()
            expected = expected_jobs.get(name)
            identity = compare_request_to_expected(payload, expected) if expected else {
                "request_identity_pass": False,
                "request_name_pass": False,
                "chain_sequences_pass": False,
                "seed_pass": False,
                "normalized_seed": "",
                "seed_serialization": "",
            }
            counts = {
                "model_cif_count": len(list(directory.glob("*_model_*.cif"))),
                "summary_confidence_count": len(list(directory.glob("*_summary_confidences_*.json"))),
                "full_data_count": len(list(directory.glob("*_full_data_*.json"))),
                "request_count": len(list(directory.glob("*_job_request.json"))),
            }
            complete = tuple(counts.values()) == (5, 5, 5, 1)
            occurrences[name].append({
                "directory": directory,
                "bundle_fingerprint": _bundle_fingerprint(directory),
                "complete": complete,
                **counts,
                **identity,
            })
    rows = []
    for expected_row in manifest:
        name = str(expected_row["job_name"]).lower()
        if name not in expected_jobs:
            raise ValueError(f"manifest job is absent from prepared batches: {name}")
        matches = occurrences.get(name, [])
        complete_exact = [
            row for row in matches if row["complete"] and row["request_identity_pass"]
        ]
        chosen = min(
            complete_exact,
            key=lambda row: (row["bundle_fingerprint"], str(row["directory"])),
            default=None,
        )
        base = chosen or (matches[0] if matches else {})
        rows.append({
            **dict(expected_row),
            "download_status": (
                "complete_exact" if chosen else "present_invalid_or_incomplete" if matches else "missing"
            ),
            "canonical_directory": str(chosen["directory"]) if chosen else "",
            "bundle_fingerprint": chosen["bundle_fingerprint"] if chosen else "",
            "observed_occurrence_count": len(matches),
            "complete_occurrence_count": len(complete_exact),
            "distinct_complete_bundle_count": len({row["bundle_fingerprint"] for row in complete_exact}),
            "duplicate_handling": "score_blind_minimum_bundle_fingerprint",
            "model_cif_count": base.get("model_cif_count", 0),
            "summary_confidence_count": base.get("summary_confidence_count", 0),
            "full_data_count": base.get("full_data_count", 0),
            "request_count": base.get("request_count", 0),
            "request_identity_pass": base.get("request_identity_pass", False),
            "request_name_pass": base.get("request_name_pass", False),
            "chain_sequences_pass": base.get("chain_sequences_pass", False),
            "seed_pass": base.get("seed_pass", False),
            "normalized_seed": base.get("normalized_seed", ""),
            "seed_serialization": base.get("seed_serialization", ""),
        })
    return rows


def choose_endpoint(
    composite_positive_ranks: Sequence[int], frozen_positive_ranks: Sequence[int]
) -> dict[str, Any]:
    if not composite_positive_ranks or not frozen_positive_ranks:
        return {
            "selected_endpoint": "frozen_exposed_ca",
            "selection_status": "not_evaluable_incomplete_outer_results",
        }
    composite_worst = max(int(value) for value in composite_positive_ranks)
    frozen_worst = max(int(value) for value in frozen_positive_ranks)
    selected = "candidate_composite" if composite_worst < frozen_worst else "frozen_exposed_ca"
    return {
        "selected_endpoint": selected,
        "selection_status": "complete_comparison",
        "candidate_composite_worst_rank": composite_worst,
        "frozen_exposed_ca_worst_rank": frozen_worst,
        "strict_worst_rank_improvement": composite_worst < frozen_worst,
    }


def sequence_identity(
    left_core: str, right_core: str, *, positions: Sequence[int] | None = None
) -> float:
    if len(left_core) != len(right_core) or not left_core:
        raise ValueError("sequence identity requires equal non-empty sequences")
    indices = tuple(range(len(left_core))) if positions is None else tuple(positions)
    if not indices or any(index < 0 or index >= len(left_core) for index in indices):
        raise ValueError("sequence identity positions are out of range")
    return sum(left_core[index] == right_core[index] for index in indices) / len(indices)


def summarize_feature_values(
    rows: Sequence[Mapping[str, Any]], features: Sequence[str] = FEATURES
) -> dict[str, float | str]:
    output: dict[str, float | str] = {}
    for feature in features:
        values = np.asarray([float(row[feature]) for row in rows], dtype=float)
        if values.size == 0:
            for suffix in ("min", "q25", "median", "q75", "max", "iqr"):
                output[f"{feature}_{suffix}"] = ""
            continue
        q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75])
        output.update({
            f"{feature}_min": round(float(values.min()), 6),
            f"{feature}_q25": round(float(q25), 6),
            f"{feature}_median": round(float(median), 6),
            f"{feature}_q75": round(float(q75), 6),
            f"{feature}_max": round(float(values.max()), 6),
            f"{feature}_iqr": round(float(q75 - q25), 6),
        })
    return output


def classify_panel(
    rows: Sequence[Mapping[str, Any]], *, expected_comparison_count: int = 26
) -> dict[str, Any]:
    complete = [row for row in rows if str(row.get("geometry_status")) == "complete"]
    ordered = sorted(complete, key=lambda row: (
        float(row["exposed_ca_rmsd_A_median"]),
        float(row["exposed_ca_rmsd_A_iqr"]),
        str(row["pair_id"]),
    ))
    positive_ranks = [
        index for index, row in enumerate(ordered, start=1)
        if str(row.get("pair_role")) == "positive"
    ]
    available_rank: int | str = max(positive_ranks) if positive_ranks else ""
    formal = (
        len(ordered) == expected_comparison_count
        and len(positive_ranks) == 1
        and sum(str(row.get("pair_role")) != "positive" for row in ordered)
        == expected_comparison_count - 1
    )
    return {
        "comparison_count": len(ordered),
        "expected_comparison_count": expected_comparison_count,
        "available_positive_rank": available_rank,
        "positive_rank": available_rank if formal else "",
        "capture_at_3": bool(int(available_rank) <= 3) if formal else "",
        "evaluation_status": "complete" if formal else "missing_required_comparisons",
        "ranking_endpoint": "frozen_exposed_ca_rmsd_A",
        "claim_boundary": CLAIM_BOUNDARY,
    }


def _sample_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.(?:json|cif)$", path.name)
    if match is None:
        raise ValueError(f"cannot recover model index from {path.name}")
    return int(match.group(1))


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "1.0", "true", "yes"}


def _batch_jobs(package: Path) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    for path in sorted((package / "alphafold_jobs").glob("hla2_controls_batch_*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, list):
            raise ValueError(f"batch is not an array: {path}")
        jobs.extend(payload)
    names = [str(job["name"]).lower() for job in jobs]
    if len(names) != len(set(names)):
        raise ValueError("prepared AlphaFold job names are not unique")
    return jobs


def _expected_request(spec: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "name": spec["job_name"],
        "modelSeeds": [int(spec["panel_seed"])],
        "sequences": [
            {"proteinChain": {"sequence": spec["mhc_alpha_sequence"]}},
            {"proteinChain": {"sequence": spec["mhc_beta_sequence"]}},
            {"proteinChain": {"sequence": spec["peptide_sequence"]}},
        ],
    }


def _new_job_specs(
    inventory: Sequence[Mapping[str, Any]], package: Path
) -> list[dict[str, Any]]:
    hla = _hla_sequences()
    specs = []
    for row in inventory:
        key = (str(row["mhc_alpha_allele"]), str(row["mhc_beta_allele"]))
        sequences = hla[key]
        specs.append({
            **dict(row),
            "source_cohort": "new_held_out_controls",
            "entity_id": row["ligand_id"],
            **sequences,
        })
    return specs


def _old_job_specs(package: Path) -> list[dict[str, Any]]:
    frozen_manifest = read_csv(package / "frozen_hy2e11_context/native_hla_calibration_manifest_24.csv")
    frozen_inventory = {
        row["job_name"].lower(): row
        for row in read_csv(OLD_MODEL_ANALYSIS / "inventory/calibration_model_inventory_24.csv")
    }
    predictions = [
        *read_csv(V2_PACKAGE / "calibration_control_binding_predictions.csv"),
        *read_csv(V2_PACKAGE / "allele_register_predictions_320.csv"),
    ]
    prediction_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for prediction in predictions:
        key = (prediction["allele"], prediction["candidate_id"])
        current = prediction_by_key.get(key)
        if current is None or current["sequence"] == prediction["sequence"]:
            prediction_by_key[key] = prediction
    hla = _hla_sequences()
    specs = []
    for row in frozen_manifest:
        old = frozen_inventory[row["job_name"].lower()]
        prediction = prediction_by_key[(row["allele"], row["entity_id"])]
        alpha, beta = "HLA-DRA*01:01", row["allele"]
        sequences = hla[(alpha, beta)]
        specs.append({
            **dict(row),
            "source_cohort": "frozen_hy2e11_reuse",
            "system_id": "SYS_BALF5_MBP_HY2E11",
            "panel_seed": int(row["server_seed"]),
            "mhc_alpha_allele": alpha,
            "mhc_beta_allele": beta,
            "core_sequence": prediction["predicted_core"],
            "core_start_1_based": prediction["core_start"],
            "register_resolution": prediction["register_resolution"],
            "canonical_directory": old["canonical_directory"],
            "download_status": old["download_status"],
            "bundle_fingerprint": old["bundle_fingerprint"],
            **sequences,
        })
    return specs


def analyze_job_bundle(
    spec: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], list[GeometrySample], dict[str, Any]]:
    directory_value = str(spec.get("canonical_directory", ""))
    if str(spec.get("download_status")) not in {"complete", "complete_exact"} or not directory_value:
        return [], [], {
            **dict(spec),
            "technical_status": "missing_preserved" if not directory_value else "incomplete",
            "observed_sample_count": 0,
            "exact_sample_count": 0,
            "clash_sample_count": 0,
            "valid_geometry_sample_count": 0,
            "claim_boundary": CLAIM_BOUNDARY,
        }
    directory = Path(directory_value)
    if not directory.is_dir():
        raise ValueError(f"frozen result directory is no longer available: {directory}")
    request_paths = list(directory.glob("*_job_request.json"))
    if len(request_paths) != 1:
        raise ValueError(f"expected one request file in {directory}")
    request_payload = json.loads(request_paths[0].read_text(encoding="utf-8"))
    request_check = compare_request_to_expected(request_payload, _expected_request(spec))
    sample_rows: list[dict[str, Any]] = []
    geometries: list[GeometrySample] = []
    summary_paths = sorted(directory.glob("*_summary_confidences_*.json"), key=_sample_index)
    for summary_path in summary_paths:
        index = _sample_index(summary_path)
        cif_paths = list(directory.glob(f"*_model_{index}.cif"))
        full_paths = list(directory.glob(f"*_full_data_{index}.json"))
        if len(cif_paths) != 1 or len(full_paths) != 1:
            sample_rows.append({
                "source_cohort": spec["source_cohort"], "job_name": spec["job_name"],
                "entity_id": spec["entity_id"], "panel_seed": spec["panel_seed"],
                "sample_index": index, "sample_status": "missing_model_or_full_data",
                "request_identity_pass": request_check["request_identity_pass"],
                "claim_boundary": CLAIM_BOUNDARY,
            })
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        model = parse_mmcif_atoms(cif_paths[0])
        observed_sequences = [residue_sequence(model.get(chain, [])) for chain in ("A", "B", "C")]
        expected_sequences = [
            str(spec["mhc_alpha_sequence"]),
            str(spec["mhc_beta_sequence"]),
            str(spec["peptide_sequence"]),
        ]
        exact = (
            request_check["request_identity_pass"]
            and set(model) == {"A", "B", "C"}
            and observed_sequences == expected_sequences
        )
        has_clash = _truth(summary.get("has_clash", False))
        sample_status = "fail_request_or_model_sequence_identity"
        if exact:
            sample_status = "excluded_has_clash" if has_clash else "pass_exact_clash_free"
            if not has_clash:
                try:
                    geometry = geometry_from_mmcif(
                        cif_paths[0],
                        ligand_id=str(spec["entity_id"]),
                        peptide_sequence=str(spec["peptide_sequence"]),
                        core_sequence=str(spec["core_sequence"]),
                        mhc_alpha_chain="A", mhc_beta_chain="B", peptide_chain="C",
                        mhc_alpha_reference_sequence=str(spec["mhc_alpha_sequence"]),
                        mhc_beta_reference_sequence=str(spec["mhc_beta_sequence"]),
                        mhc_alpha_reference_start=0, mhc_beta_reference_start=0,
                    )
                    geometries.append(GeometrySample(index, geometry))
                except (KeyError, StopIteration, ValueError) as error:
                    sample_status = f"fail_geometry_coordinates:{type(error).__name__}"
        sample_rows.append({
            "source_cohort": spec["source_cohort"],
            "system_id": spec["system_id"],
            "job_name": spec["job_name"],
            "entity_id": spec["entity_id"],
            "panel_seed": spec["panel_seed"],
            "sample_index": index,
            "sample_status": sample_status,
            "request_identity_pass": request_check["request_identity_pass"],
            "model_sequence_identity_pass": exact,
            "observed_mhc_alpha_length": len(observed_sequences[0]),
            "observed_mhc_beta_length": len(observed_sequences[1]),
            "observed_peptide": observed_sequences[2],
            "core_sequence": spec["core_sequence"],
            "core_start_1_based": spec["core_start_1_based"],
            "ranking_score": summary.get("ranking_score", ""),
            "iptm": summary.get("iptm", ""),
            "ptm": summary.get("ptm", ""),
            "has_clash": has_clash,
            "cif_path": str(cif_paths[0]),
            "summary_path": str(summary_path),
            "full_data_path": str(full_paths[0]),
            "claim_boundary": CLAIM_BOUNDARY,
        })
    job_summary = {
        **dict(spec),
        "technical_status": "geometry_evaluable" if geometries else "excluded_no_valid_geometry_samples",
        "request_identity_pass": request_check["request_identity_pass"],
        "observed_sample_count": len(sample_rows),
        "exact_sample_count": sum(
            str(row["sample_status"]).startswith(("pass", "excluded_has_clash"))
            for row in sample_rows
        ),
        "clash_sample_count": sum(bool(row.get("has_clash")) for row in sample_rows),
        "valid_geometry_sample_count": len(geometries),
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return sample_rows, geometries, job_summary


def _ligand_metadata(package: Path) -> dict[tuple[str, str], dict[str, Any]]:
    output: dict[tuple[str, str], dict[str, Any]] = {}
    positives = {
        row["ligand_id"]: {
            "ligand_id": row["ligand_id"],
            "sequence": row["sequence"],
            "core_sequence": row["core"],
            "mhc_alpha_allele": row["mhc_alpha_allele"],
            "mhc_beta_allele": row["mhc_beta_allele"],
        }
        for row in read_csv(package / "registry/control_ligand_registry.csv")
    }
    for pair in read_csv(package / "registry/positive_pair_registry.csv"):
        for field in ("left_ligand_id", "right_ligand_id"):
            ligand_id = pair[field]
            output[(pair["pair_id"], ligand_id)] = positives[ligand_id]
    for row in read_csv(package / "controls/control_decoy_registry.csv"):
        key = (row["positive_pair_id"], row["candidate_id"])
        current = output.get(key)
        value = {
            "ligand_id": row["candidate_id"],
            "sequence": row["sequence"],
            "core_sequence": row["predicted_core"],
            "mhc_alpha_allele": row["mhc_alpha_allele"],
            "mhc_beta_allele": row["mhc_beta_allele"],
        }
        if current is not None and current != value:
            raise ValueError(f"inconsistent comparator metadata for {key}")
        output[key] = value
    return output


def _summarize_comparison(
    base: Mapping[str, Any],
    left: Sequence[GeometrySample],
    right: Sequence[GeometrySample],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    combinations: list[dict[str, Any]] = []
    metric_rows: list[dict[str, float]] = []
    for left_sample in left:
        for right_sample in right:
            values = pair_features(left_sample.geometry, right_sample.geometry)
            metric_rows.append(values)
            combinations.append({
                "system_id": base["system_id"],
                "positive_pair_id": base["positive_pair_id"],
                "panel_seed": base["panel_seed"],
                "pair_id": base["pair_id"],
                "pair_role": base["pair_role"],
                "left_id": base["left_id"],
                "right_id": base["right_id"],
                "left_sample_index": left_sample.sample_index,
                "right_sample_index": right_sample.sample_index,
                **{key: round(value, 6) for key, value in values.items()},
                "interpretation": "Technical AlphaFold sample-combination sensitivity only; not biological replication.",
            })
    summary = {
        **dict(base),
        "left_valid_sample_count": len(left),
        "right_valid_sample_count": len(right),
        "model_combination_count": len(metric_rows),
        "geometry_status": "complete" if metric_rows else "missing_or_qc_excluded_model",
        **summarize_feature_values(metric_rows, FEATURES),
        "tcr_facing_sequence_identity": round(
            sequence_identity(str(base["left_core"]), str(base["right_core"]), positions=tuple(EXPOSED_INDICES)),
            9,
        ),
        "full_core_sequence_identity": round(
            sequence_identity(str(base["left_core"]), str(base["right_core"])), 9
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return summary, combinations


def _build_af3_feature_matrix(
    package: Path,
    specs: Sequence[Mapping[str, Any]],
    geometry_by_job: Mapping[str, Sequence[GeometrySample]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    jobs_by_sequence: dict[tuple[int, str, str, str], Mapping[str, Any]] = {}
    old_jobs_by_entity: dict[tuple[int, str], Mapping[str, Any]] = {}
    for spec in specs:
        sequence_key = (
            int(spec["panel_seed"]), str(spec["mhc_alpha_allele"]),
            str(spec["mhc_beta_allele"]), str(spec["peptide_sequence"]),
        )
        current = jobs_by_sequence.get(sequence_key)
        if current is not None and str(current["job_name"]) != str(spec["job_name"]):
            raise ValueError(f"multiple AlphaFold jobs for the same frozen ligand key: {sequence_key}")
        jobs_by_sequence[sequence_key] = spec
        if spec["source_cohort"] == "frozen_hy2e11_reuse":
            old_jobs_by_entity[(int(spec["panel_seed"]), str(spec["entity_id"]))] = spec
    ligand = _ligand_metadata(package)
    summaries: list[dict[str, Any]] = []
    ensembles: list[dict[str, Any]] = []
    for comparison in read_csv(package / "controls/new_control_comparison_universe.csv"):
        seed = int(comparison["panel_seed"])
        left_meta = ligand[(comparison["positive_pair_id"], comparison["left_id"])]
        right_meta = ligand[(comparison["positive_pair_id"], comparison["right_id"])]
        left_job = jobs_by_sequence[(
            seed, left_meta["mhc_alpha_allele"], left_meta["mhc_beta_allele"], left_meta["sequence"],
        )]
        right_job = jobs_by_sequence[(
            seed, right_meta["mhc_alpha_allele"], right_meta["mhc_beta_allele"], right_meta["sequence"],
        )]
        base = {
            "system_id": comparison["system_id"],
            "positive_pair_id": comparison["positive_pair_id"],
            "panel_seed": seed,
            "pair_id": comparison["pair_id"],
            "pair_role": "positive" if comparison["comparison_role"] == "positive" else "N3",
            "negative_tier": comparison["negative_tier"],
            "analysis_set": "primary_rank_of_26",
            "left_id": comparison["left_id"],
            "right_id": comparison["right_id"],
            "left_core": left_meta["core_sequence"],
            "right_core": right_meta["core_sequence"],
            "left_job_name": left_job["job_name"],
            "right_job_name": right_job["job_name"],
        }
        summary, combinations = _summarize_comparison(
            base,
            geometry_by_job.get(str(left_job["job_name"]).lower(), []),
            geometry_by_job.get(str(right_job["job_name"]).lower(), []),
        )
        summaries.append(summary)
        ensembles.extend(combinations)
    for comparison in read_csv(package / "frozen_hy2e11_context/calibration_comparison_universe_72.csv"):
        seed = int(comparison["seed"])
        left_job = old_jobs_by_entity[(seed, comparison["viral_entity_id"])]
        right_job = old_jobs_by_entity[(seed, comparison["self_entity_id"])]
        base = {
            "system_id": "SYS_BALF5_MBP_HY2E11",
            "positive_pair_id": "PAIR_HY2E11_BALF5_MBP",
            "panel_seed": seed,
            "pair_id": comparison["pair_id"],
            "pair_role": "positive" if comparison["pair_role"] == "E1_positive" else "N3",
            "negative_tier": "positive" if comparison["pair_role"] == "E1_positive" else "N3",
            "analysis_set": comparison["analysis_set"],
            "left_id": comparison["viral_entity_id"],
            "right_id": comparison["self_entity_id"],
            "left_core": left_job["core_sequence"],
            "right_core": right_job["core_sequence"],
            "left_job_name": left_job["job_name"],
            "right_job_name": right_job["job_name"],
        }
        summary, combinations = _summarize_comparison(
            base,
            geometry_by_job.get(str(left_job["job_name"]).lower(), []),
            geometry_by_job.get(str(right_job["job_name"]).lower(), []),
        )
        summaries.append(summary)
        ensembles.extend(combinations)
    return summaries, ensembles


def _af3_evaluation_rows(feature_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, int], list[Mapping[str, Any]]] = defaultdict(list)
    for row in feature_rows:
        if row["analysis_set"] == "primary_rank_of_26":
            grouped[(row["system_id"], row["positive_pair_id"], int(row["panel_seed"]))].append(row)
    output = []
    for (system_id, positive_pair_id, seed), rows in sorted(grouped.items()):
        result = classify_panel(rows)
        positive = next((row for row in rows if row["pair_role"] == "positive"), None)
        output.append({
            "system_id": system_id,
            "pair_id": positive_pair_id,
            "positive_pair_id": positive_pair_id,
            "layer": "af3",
            "panel_seed": seed,
            "required_for_system_pass": True,
            **result,
            **{
                feature: positive.get(f"{feature}_median", "") if positive else ""
                for feature in FEATURES
            },
        })
    return output


def _complete_scoring_rows(
    package: Path,
    af3_rows: Sequence[Mapping[str, Any]],
    af3_evaluations: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    complete_af3 = {
        (row["system_id"], row["pair_id"], str(row["panel_seed"]))
        for row in af3_evaluations if row["evaluation_status"] == "complete"
    }
    output = []
    for row in af3_rows:
        key = (row["system_id"], row["positive_pair_id"], str(row["panel_seed"]))
        if row["analysis_set"] != "primary_rank_of_26" or key not in complete_af3:
            continue
        output.append({
            "system_id": row["system_id"],
            "positive_pair_id": row["positive_pair_id"],
            "pair_id": row["pair_id"],
            "pair_role": row["pair_role"],
            "layer": "af3",
            "panel_seed": row["panel_seed"],
            **{feature: float(row[f"{feature}_median"]) for feature in FEATURES},
            "tcr_facing_sequence_identity": float(row["tcr_facing_sequence_identity"]),
            "full_core_sequence_identity": float(row["full_core_sequence_identity"]),
        })
    pdb_evaluations = {
        row["pair_id"]: row for row in read_csv(package / "benchmark/evaluation_status.csv")
        if row["layer"] == "pdb_oracle"
    }
    structural = {
        row["ligand_id"]: row
        for row in read_csv(package / "benchmark/pdb_structural_ligand_registry.csv")
    }
    for row in read_csv(package / "benchmark/pdb_oracle_feature_matrix.csv"):
        if pdb_evaluations[row["positive_pair_id"]]["evaluation_status"] != "complete":
            continue
        left_core = structural[row["left_ligand_id"]]["core_sequence"]
        right_core = structural[row["right_ligand_id"]]["core_sequence"]
        output.append({
            "system_id": row["system_id"],
            "positive_pair_id": row["positive_pair_id"],
            "pair_id": row["pair_id"],
            "pair_role": "positive" if row["pair_role"] == "positive" else "N3",
            "layer": "pdb_oracle",
            "panel_seed": "pdb",
            **{feature: float(row[feature]) for feature in FEATURES},
            "tcr_facing_sequence_identity": sequence_identity(
                left_core, right_core, positions=tuple(EXPOSED_INDICES)
            ),
            "full_core_sequence_identity": sequence_identity(left_core, right_core),
        })
    return output


def _rank_metric(group: Sequence[Mapping[str, Any]], metric: str, *, higher_better: bool) -> int:
    ordered = sorted(group, key=lambda row: (
        -float(row[metric]) if higher_better else float(row[metric]), str(row["pair_id"])
    ))
    ranks = [index for index, row in enumerate(ordered, start=1) if row["pair_role"] == "positive"]
    if len(ranks) != 1:
        raise ValueError("every complete panel must contain exactly one positive")
    return ranks[0]


def _baseline_and_permutation_rows(
    scoring_rows: Sequence[Mapping[str, Any]], selected_rank_by_panel: Mapping[tuple[str, str, str], int]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for row in scoring_rows:
        grouped[(row["system_id"], row["positive_pair_id"], row["layer"], str(row["panel_seed"]))].append(row)
    baselines: list[dict[str, Any]] = []
    permutations: list[dict[str, Any]] = []
    for (system_id, pair_id, layer, seed), group in sorted(grouped.items()):
        panel_key = (pair_id, layer, seed)
        ranks = {
            "frozen_exposed_ca_rmsd_A": _rank_metric(group, "exposed_ca_rmsd_A", higher_better=False),
            "tcr_facing_sequence_identity": _rank_metric(group, "tcr_facing_sequence_identity", higher_better=True),
            "full_core_sequence_identity": _rank_metric(group, "full_core_sequence_identity", higher_better=True),
        }
        for method, rank in ranks.items():
            baselines.append({
                "system_id": system_id, "positive_pair_id": pair_id, "layer": layer,
                "panel_seed": seed, "method": method, "positive_rank": rank,
                "comparison_count": len(group), "capture_at_3": rank <= 3,
            })
        baselines.append({
            "system_id": system_id, "positive_pair_id": pair_id, "layer": layer,
            "panel_seed": seed, "method": "random_rank_expectation",
            "positive_rank": round((len(group) + 1) / 2.0, 6),
            "comparison_count": len(group),
            "capture_at_3": round(min(3, len(group)) / len(group), 6),
        })
        observed_rank = int(selected_rank_by_panel[panel_key])
        seed_value = int(hashlib.sha256("|".join(panel_key).encode("utf-8")).hexdigest()[:16], 16)
        draws = np.random.default_rng(seed_value).integers(1, len(group) + 1, size=10000)
        permutations.append({
            "system_id": system_id, "positive_pair_id": pair_id, "layer": layer,
            "panel_seed": seed, "permutation_count": 10000,
            "observed_selected_endpoint_rank": observed_rank,
            "empirical_probability_random_rank_at_or_better": round(float(np.mean(draws <= observed_rank)), 6),
            "empirical_random_capture_at_3": round(float(np.mean(draws <= 3)), 6),
            "analytic_random_capture_at_3": round(min(3, len(group)) / len(group), 6),
            "random_seed_rule": "sha256(positive_pair_id|layer|panel_seed)",
        })
    return baselines, permutations


def _file_checksums(out: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted((value for value in out.rglob("*") if value.is_file() and value.name != "SHA256SUMS.csv"), key=str):
        rows.append({
            "relative_path": str(path.relative_to(out)),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        })
    return rows


def run_analysis(result_roots: Sequence[Path], out: Path = DEFAULT_OUT, package: Path = PACKAGE) -> dict[str, Any]:
    out.mkdir(parents=True, exist_ok=True)
    manifest = read_csv(package / "alphafold_jobs/job_manifest.csv")
    batch_jobs = _batch_jobs(package)
    new_inventory = inventory_downloaded_jobs(manifest, batch_jobs, result_roots)
    if any(row["download_status"] != "complete_exact" for row in new_inventory):
        missing = [row["job_name"] for row in new_inventory if row["download_status"] != "complete_exact"]
        raise ValueError(f"new held-out result set is incomplete or inexact: {missing}")
    observed_names = []
    for root in result_roots:
        for request_path in Path(root).rglob("*_job_request.json"):
            observed_names.append(str(_request_job(json.loads(request_path.read_text(encoding="utf-8")))["name"]).lower())
    expected_names = {str(row["job_name"]).lower() for row in manifest}
    unexpected = sorted(set(observed_names) - expected_names)
    if unexpected:
        raise ValueError(f"unexpected jobs in selected result roots: {unexpected}")
    write_csv(out / "inventory/new_job_inventory_48.csv", new_inventory)

    for source, destination in (
        ("registry/literature_and_structure_sources.csv", "provenance/literature_and_structure_sources.csv"),
        ("registry/control_system_registry.csv", "provenance/control_system_registry.csv"),
        ("registry/control_ligand_registry.csv", "provenance/control_ligand_registry.csv"),
        ("registry/positive_pair_registry.csv", "provenance/positive_pair_registry.csv"),
        ("controls/control_decoy_registry.csv", "provenance/control_decoy_registry.csv"),
        ("controls/new_control_comparison_universe.csv", "provenance/new_control_comparison_universe.csv"),
        ("benchmark/pdb_oracle_feature_matrix.csv", "benchmark/pdb_oracle_feature_matrix.csv"),
        ("benchmark/pdb_positive_feature_matrix.csv", "benchmark/pdb_positive_feature_matrix.csv"),
        ("benchmark/pdb_structural_ligand_registry.csv", "provenance/pdb_structural_ligand_registry.csv"),
    ):
        write_csv(out / destination, read_csv(package / source))

    old_specs = _old_job_specs(package)
    old_inventory_rows = [{
        key: value for key, value in spec.items()
        if key not in {"mhc_alpha_sequence", "mhc_beta_sequence"}
    } for spec in old_specs]
    write_csv(out / "inventory/reused_frozen_hy2e11_job_inventory_24.csv", old_inventory_rows)
    new_specs = _new_job_specs(new_inventory, package)
    all_specs = [*new_specs, *old_specs]
    sample_rows: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    geometry_by_job: dict[str, list[GeometrySample]] = {}
    for spec in all_specs:
        samples, geometries, job = analyze_job_bundle(spec)
        sample_rows.extend(samples)
        job_rows.append(job)
        geometry_by_job[str(spec["job_name"]).lower()] = geometries
    write_csv(out / "qc/model_sample_qc.csv", sample_rows)
    write_csv(out / "qc/job_qc_summary.csv", job_rows)

    feature_rows, ensemble_rows = _build_af3_feature_matrix(package, all_specs, geometry_by_job)
    write_csv(out / "benchmark/af3_model_ensemble.csv", ensemble_rows)
    write_csv(out / "benchmark/af3_pair_feature_matrix.csv", feature_rows)
    af3_evaluations = _af3_evaluation_rows(feature_rows)
    scoring_rows = _complete_scoring_rows(package, feature_rows, af3_evaluations)
    outer_results = leave_one_system_out(scoring_rows, FEATURES)
    frozen_rank_lookup = {}
    for row in af3_evaluations:
        if row["evaluation_status"] == "complete":
            frozen_rank_lookup[(row["pair_id"], "af3", str(row["panel_seed"]))] = int(row["positive_rank"])
    for row in read_csv(package / "benchmark/evaluation_status.csv"):
        if row["layer"] == "pdb_oracle" and row["evaluation_status"] == "complete":
            frozen_rank_lookup[(row["pair_id"], "pdb_oracle", "pdb")] = int(row["positive_rank"])
    for row in outer_results:
        key = (row["positive_pair_id"], row["layer"], str(row["panel_seed"]))
        row["frozen_exposed_ca_positive_rank"] = frozen_rank_lookup[key]
        row["candidate_composite_positive_rank"] = row["positive_rank"]
        row["weights_frozen"] = False
    endpoint = choose_endpoint(
        [int(row["candidate_composite_positive_rank"]) for row in outer_results],
        [int(row["frozen_exposed_ca_positive_rank"]) for row in outer_results],
    )
    selected_endpoint = endpoint["selected_endpoint"]
    selected_rank_by_panel = {
        (row["positive_pair_id"], row["layer"], str(row["panel_seed"])): int(
            row["candidate_composite_positive_rank"]
            if selected_endpoint == "candidate_composite"
            else row["frozen_exposed_ca_positive_rank"]
        )
        for row in outer_results
    }
    for row in outer_results:
        key = (row["positive_pair_id"], row["layer"], str(row["panel_seed"]))
        row["selected_endpoint"] = selected_endpoint
        row["selected_positive_rank"] = selected_rank_by_panel[key]
        row["selected_capture_at_3"] = selected_rank_by_panel[key] <= 3
    write_csv(out / "benchmark/outer_fold_results.csv", outer_results)
    weight_rows = []
    for held_out in sorted({row["held_out_system_id"] for row in outer_results}):
        row = next(value for value in outer_results if value["held_out_system_id"] == held_out)
        weight_rows.append({
            "held_out_system_id": held_out,
            "training_system_ids": row["training_system_ids"],
            **{f"weight_{feature}": row[f"weight_{feature}"] for feature in FEATURES},
            "weights_frozen": False,
            "freeze_block_reason": "formal_trust_gate_not_passed",
        })
    write_csv(out / "benchmark/selected_weights.csv", weight_rows)
    write_json(out / "benchmark/endpoint_selection.json", {
        **endpoint,
        "weights_frozen": False,
        "selection_rule": "Use candidate composite only for strict improvement in worst held-out rank; otherwise retain frozen exposed-CA.",
    })

    pdb_rows = [
        dict(row) for row in read_csv(package / "benchmark/evaluation_status.csv")
        if row["layer"] == "pdb_oracle"
    ]
    final_evaluations = [*pdb_rows, *af3_evaluations]
    for row in final_evaluations:
        key = (row["pair_id"], row["layer"], str(row["panel_seed"]))
        row["positive_pair_id"] = row["pair_id"]
        row["frozen_exposed_ca_positive_rank"] = row.get("available_positive_rank", "")
        row["decoy_count"] = max(0, int(row.get("comparison_count", 0)) - 1)
        if row["evaluation_status"] == "complete":
            row["positive_rank"] = selected_rank_by_panel[key]
            row["selected_positive_rank"] = selected_rank_by_panel[key]
            row["capture_at_3"] = selected_rank_by_panel[key] <= 3
            row["ranking_endpoint"] = selected_endpoint
        else:
            row["selected_positive_rank"] = ""
        row["weights_frozen"] = False
    final_evaluations.sort(key=lambda row: (
        row["system_id"], row["pair_id"], 0 if row["layer"] == "pdb_oracle" else 1,
        str(row["panel_seed"]),
    ))
    write_csv(out / "benchmark/evaluation_status.csv", final_evaluations)
    systems, _, _, _ = curated_registry()
    strict_system_ids = [row["system_id"] for row in systems if row["eligibility"] == "strict"]
    gate = build_trust_gate(final_evaluations, required_system_ids=strict_system_ids)
    gate.update({
        "selected_endpoint": selected_endpoint,
        "weights_frozen": False,
        "discovery_rankings_changed": False,
        "cross_allele_consensus_created": False,
        "blocking_evaluations": [
            {
                "system_id": row["system_id"], "pair_id": row["pair_id"],
                "layer": row["layer"], "panel_seed": row["panel_seed"],
                "evaluation_status": row["evaluation_status"],
            }
            for row in final_evaluations if row["evaluation_status"] != "complete"
        ],
    })
    write_json(out / "benchmark/trust_gate.json", gate)
    baselines, permutations = _baseline_and_permutation_rows(scoring_rows, selected_rank_by_panel)
    write_csv(out / "benchmark/baseline_comparisons.csv", baselines)
    write_csv(out / "benchmark/permutation_results.csv", permutations)

    missing_old = sorted(spec["job_name"] for spec in old_specs if not spec["canonical_directory"])
    verification = {
        "new_expected_job_count": len(manifest),
        "new_complete_exact_job_count": sum(row["download_status"] == "complete_exact" for row in new_inventory),
        "new_unexpected_job_count": len(unexpected),
        "new_duplicate_job_count": sum(int(row["observed_occurrence_count"]) > 1 for row in new_inventory),
        "new_all_five_model_bundles": all(
            (int(row["model_cif_count"]), int(row["summary_confidence_count"]), int(row["full_data_count"]), int(row["request_count"]))
            == (5, 5, 5, 1) for row in new_inventory
        ),
        "new_all_request_names_seeds_sequences_exact": all(row["request_identity_pass"] for row in new_inventory),
        "frozen_missing_jobs_preserved": missing_old,
        "frozen_missing_job_count": len(missing_old),
        "strict_biological_system_count": len(strict_system_ids),
        "required_positive_pair_count": 4,
        "evaluation_row_count": len(final_evaluations),
        "discovery_files_read": False,
        "discovery_files_written": False,
        "cross_allele_consensus_created": False,
        "overall_trust_status": gate["overall_trust_status"],
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "validation/package_verification_summary.json", verification)
    manifest_payload = {
        "analysis_version": "EBV_MS_HLA2_HELD_OUT_CONTROLS_RESULTS_V1",
        "frozen_control_package": str(package),
        "selected_result_roots": [str(Path(path).resolve()) for path in result_roots],
        "output_directory": str(out),
        "submission_state_transition": "prepared_not_submitted_to_downloaded_results_analyzed",
        "selected_endpoint": selected_endpoint,
        "weights_frozen": False,
        "overall_trust_status": gate["overall_trust_status"],
        "discovery_reranking_allowed": gate["discovery_reranking_allowed"],
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "analysis_manifest.json", manifest_payload)
    readme = f"""# Held-out HLA-II positive-control result analysis

All {len(manifest)} newly prepared jobs were recovered as exact five-model bundles. The two original Hy.2E11 missing jobs remain missing and were neither retried nor replaced.

Formal trust status: `{gate['overall_trust_status']}`. Selected descriptive endpoint: `{selected_endpoint}`. Candidate composite weights remain unfrozen. Discovery rankings are unchanged, and no cross-allele consensus output was created.

This package reports computational pMHC geometry prioritization only. It is not evidence of antigen presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or an MS mechanism.
"""
    (out / "README.md").write_text(readme, encoding="utf-8")
    checksums = _file_checksums(out)
    write_csv(out / "SHA256SUMS.csv", checksums)
    return {
        "out": str(out),
        "new_complete_exact_job_count": verification["new_complete_exact_job_count"],
        "frozen_missing_jobs": missing_old,
        "valid_geometry_sample_count": sum(int(row["valid_geometry_sample_count"]) for row in job_rows),
        "evaluation_rows": final_evaluations,
        "selected_endpoint": selected_endpoint,
        "trust_gate": gate,
        "package_checksum_rows": len(checksums),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--result-root", action="append", type=Path, default=[])
    args = parser.parse_args()
    roots = args.result_root or sorted(
        path for path in DEFAULT_DOWNLOADS.glob("folds_2026_08_26_*") if path.is_dir()
    )
    result = run_analysis(roots, out=args.out)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
