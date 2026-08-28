"""Auditable, lead-focused robustness calculations for structural analyses."""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from html import escape
from pathlib import Path
from statistics import median

import numpy as np

from analyze_af3_pmhc_downloads import (
    candidate_id_from_request_name,
    parse_mmcif,
    request_details,
    residue_plddt,
    sequence,
)
from same_register_af3_analysis import same_register_geometry


FIXED_LEADS = {
    1: {
        "pair_id": "EBV_TCELL_950::HUMAN_MYELIN_112214",
        "analysis_layer": "strict_primary_controls",
        "ebv_candidate_id": "EBV_TCELL_950",
        "target_candidate_id": "HUMAN_MYELIN_112214",
        "control_ids": (
            "HUMAN_BACKGROUND_115891",
            "HUMAN_BACKGROUND_118550",
            "HUMAN_BACKGROUND_119732",
        ),
    },
    2: {
        "pair_id": "EBV_TCELL_2268741::HUMAN_MYELIN_117032",
        "analysis_layer": "length_sensitivity_exact_bin_pm7",
        "ebv_candidate_id": "EBV_TCELL_2268741",
        "target_candidate_id": "HUMAN_MYELIN_117032",
        "control_ids": (
            "HUMAN_BACKGROUND_141561",
            "HUMAN_BACKGROUND_423369",
            "HUMAN_BACKGROUND_2258889",
        ),
    },
}

OUTPUT_DIRECTORY_NAME = "lead_focused_robustness_2026-08-15"
INPUT_TABLES = {
    "rank1_control_geometry": "processed/complete_model_pipeline_audit_2026-08-15/matched_background_structure_geometry.csv",
    "target_geometry": "processed/complete_model_pipeline_audit_2026-08-15/combined_same_register_geometry.csv",
    "legacy_job_summary": "processed/complete_model_pipeline_audit_2026-08-15/canonical_af3_job_summary.csv",
    "legacy_sample_metrics": "processed/complete_model_pipeline_audit_2026-08-15/canonical_af3_sample_metrics.csv",
    "score_sheet": "processed/structural_control_expansion_2026-08-15/master_pair_score_sheet_with_expanded_controls.csv",
    "rank2_control_geometry": "processed/structural_control_expansion_2026-08-15/complete_layered_control_geometry.csv",
    "new_control_sample_metrics": "processed/structural_control_expansion_2026-08-15/alphafold_control_sample_metrics.csv",
    "rank1_control_predictions": "processed/expanded_background/background_register_prediction_summary.csv",
    "rank2_control_predictions": "processed/structural_control_expansion_2026-08-15/control_binding_prediction_summary.csv",
}


def validate_lead_definition(lead):
    """Require the two frozen leads, controls, and non-poolable layers exactly."""
    rank = int(lead.get("rank", 0))
    if rank not in FIXED_LEADS:
        raise ValueError("Only frozen lead ranks 1 and 2 are permitted")
    expected = FIXED_LEADS[rank]
    for field in ("pair_id", "analysis_layer", "ebv_candidate_id", "target_candidate_id"):
        if lead.get(field) != expected[field]:
            raise ValueError(
                f"Rank {rank} requires {field}={expected[field]!r}; got {lead.get(field)!r}"
            )
    if set(lead.get("control_ids", ())) != set(expected["control_ids"]):
        raise ValueError(f"Rank {rank} requires the frozen three-control set")
    if rank == 2 and (
        lead.get("target_peptide_length") != 32 or lead.get("stratum_length") != 32
    ):
        raise ValueError(
            "Rank 2 sensitivity output requires an exactly 32-aa target and exactly 32-aa frozen stratum"
        )
    return expected


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _unique_core_start(peptide, core, candidate_id):
    starts = [
        index + 1
        for index in range(0, len(peptide) - len(core) + 1)
        if peptide[index:index + len(core)] == core
    ]
    if len(starts) != 1:
        raise ValueError(
            f"{candidate_id} core placement is ambiguous or absent: {core!r} in {peptide!r}"
        )
    return starts[0]


def _resolve_project_path(project_root, value):
    path = Path(value)
    return path if path.is_absolute() else Path(project_root) / path


def load_fixed_lead_metadata(project_root):
    """Load and strictly validate the frozen rows that define both analyses."""
    project_root = Path(project_root).resolve()
    paths = {name: project_root / relative for name, relative in INPUT_TABLES.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing required frozen input tables: {missing}")
    rows = {name: _read_csv(path) for name, path in paths.items()}
    score_by_pair = {}
    for rank, fixed in FIXED_LEADS.items():
        matches = [row for row in rows["score_sheet"] if row["pair_id"] == fixed["pair_id"]]
        if len(matches) != 1:
            raise ValueError(f"Expected one score-sheet row for {fixed['pair_id']}")
        score = matches[0]
        if int(score["discovery_priority_rank"]) != rank:
            raise ValueError(f"Frozen discovery rank changed for {fixed['pair_id']}")
        score_by_pair[fixed["pair_id"]] = score

    rank1 = FIXED_LEADS[1]
    rank1_geometry = [
        row for row in rows["rank1_control_geometry"]
        if row["pair_id"] == rank1["pair_id"]
    ]
    if not rank1_geometry or {row["background_candidate_id"] for row in rank1_geometry} != set(rank1["control_ids"]):
        raise ValueError("Rank 1 frozen strict-primary control set is incomplete or changed")

    rank2 = FIXED_LEADS[2]
    rank2_geometry = [
        row for row in rows["rank2_control_geometry"]
        if row["pair_id"] == rank2["pair_id"]
    ]
    if not rank2_geometry:
        raise ValueError("Rank 2 frozen length-sensitivity geometry is absent")
    require_single_analysis_layer(rank2_geometry)
    if {row["analysis_layer"] for row in rank2_geometry} != {rank2["analysis_layer"]}:
        raise ValueError(
            f"Rank 2 requires analysis layer {rank2['analysis_layer']}"
        )
    if {row["background_candidate_id"] for row in rank2_geometry} != set(rank2["control_ids"]):
        raise ValueError("Rank 2 frozen three-control set is incomplete or changed")
    try:
        rank2_strata = {int(row["stratum_length"]) for row in rank2_geometry}
    except (KeyError, ValueError) as error:
        raise ValueError("Rank 2 frozen stratum length is unavailable") from error
    if rank2_strata != {32}:
        raise ValueError("Rank 2 sensitivity output requires an exactly 32-aa frozen stratum")
    rank2_score = score_by_pair[rank2["pair_id"]]
    if len(rank2_score["human_peptide"]) != 32:
        raise ValueError("Rank 2 sensitivity output requires an exactly 32-aa target peptide")

    for fixed in FIXED_LEADS.values():
        target_rows = [
            row for row in rows["target_geometry"] if row["pair_id"] == fixed["pair_id"]
        ]
        if not target_rows:
            raise ValueError(f"Frozen target geometry is absent for {fixed['pair_id']}")

    candidate_metadata = {}
    for rank, fixed in FIXED_LEADS.items():
        score = score_by_pair[fixed["pair_id"]]
        for candidate_id, peptide_field, core_field, role in (
            (fixed["ebv_candidate_id"], "ebv_peptide", "ebv_p1_p9_core", "ebv"),
            (fixed["target_candidate_id"], "human_peptide", "human_p1_p9_core", "target"),
        ):
            peptide, core = score[peptide_field], score[core_field]
            if len(core) != 9:
                raise ValueError(f"{candidate_id} does not have an exact frozen P1-P9 core")
            candidate_metadata[candidate_id] = {
                "candidate_id": candidate_id,
                "peptide": peptide,
                "core": core,
                "core_start_1_based": _unique_core_start(peptide, core, candidate_id),
                "lead_rank": rank,
                "entity_role": role,
            }

    prediction_sources = {
        1: rows["rank1_control_predictions"],
        2: rows["rank2_control_predictions"],
    }
    geometry_sources = {1: rank1_geometry, 2: rank2_geometry}
    for rank, fixed in FIXED_LEADS.items():
        for control_id in fixed["control_ids"]:
            predictions = [
                row for row in prediction_sources[rank] if row["candidate_id"] == control_id
            ]
            if len(predictions) != 1:
                raise ValueError(f"Expected one frozen control prediction for {control_id}")
            prediction = predictions[0]
            if prediction.get("prediction_status") != "predicted":
                raise ValueError(f"Control prediction is not complete for {control_id}")
            peptide = prediction["peptide"]
            core = prediction["predicted_core_peptide"]
            start = _unique_core_start(peptide, core, control_id)
            stored_starts = prediction["predicted_core_start_positions_1_based"].split(";")
            if len(stored_starts) != 1 or not stored_starts[0].isdigit() or int(stored_starts[0]) != start:
                raise ValueError(f"Ambiguous or inconsistent stored control core placement for {control_id}")
            geometry_cores = {
                row["background_predicted_core"]
                for row in geometry_sources[rank]
                if row["background_candidate_id"] == control_id
            }
            if geometry_cores != {core}:
                raise ValueError(f"Geometry/prediction core mismatch for {control_id}")
            candidate_metadata[control_id] = {
                "candidate_id": control_id,
                "peptide": peptide,
                "core": core,
                "core_start_1_based": start,
                "lead_rank": rank,
                "entity_role": "control",
            }

    legacy_candidates = {
        identifier
        for rank, fixed in FIXED_LEADS.items()
        for identifier in (
            fixed["ebv_candidate_id"], fixed["target_candidate_id"],
            *(fixed["control_ids"] if rank == 1 else ()),
        )
    }
    legacy_jobs = [
        row for row in rows["legacy_job_summary"] if row["candidate_id"] in legacy_candidates
    ]
    if {row["candidate_id"] for row in legacy_jobs} != legacy_candidates:
        raise ValueError("Legacy canonical job summary is missing an involved entity")
    for row in legacy_jobs:
        if not Path(row["source_path"]).is_dir():
            raise FileNotFoundError(f"Saved legacy job path is unavailable: {row['source_path']}")

    legacy_samples = [
        row for row in rows["legacy_sample_metrics"] if row["candidate_id"] in legacy_candidates
    ]
    sample_indices_by_job = defaultdict(set)
    for row in legacy_samples:
        model_path = _resolve_project_path(project_root, row["model_path"])
        if not model_path.is_file():
            raise FileNotFoundError(f"Saved legacy model path is unavailable: {model_path}")
        sample_indices_by_job[row["canonical_job_key"]].add(int(row["sample_index"]))
    if any(indices != set(range(5)) for indices in sample_indices_by_job.values()):
        raise ValueError("Legacy canonical sample metrics do not contain exact model indices 0-4")

    new_controls = set(FIXED_LEADS[2]["control_ids"])
    new_samples = [
        row for row in rows["new_control_sample_metrics"]
        if row["candidate_id"] in new_controls
    ]
    new_indices_by_job = defaultdict(set)
    for row in new_samples:
        model_path = _resolve_project_path(project_root, row["model_path"])
        if not model_path.is_file():
            raise FileNotFoundError(f"Saved new-control model path is unavailable: {model_path}")
        new_indices_by_job[row["canonical_job_key"]].add(int(row["sample_index"]))
    if set(row["candidate_id"] for row in new_samples) != new_controls:
        raise ValueError("New-control sample metrics are missing an involved entity")
    if any(indices != set(range(5)) for indices in new_indices_by_job.values()):
        raise ValueError("New-control sample metrics do not contain exact model indices 0-4")

    lead_specs = []
    for rank, fixed in sorted(FIXED_LEADS.items()):
        lead_spec = {
            "rank": rank,
            **fixed,
            "control_ids": list(fixed["control_ids"]),
        }
        if rank == 2:
            lead_spec.update({
                "target_peptide_length": len(score_by_pair[fixed["pair_id"]]["human_peptide"]),
                "stratum_length": 32,
            })
        lead_specs.append(lead_spec)
    return {
        "lead_specs": lead_specs,
        "candidate_metadata": candidate_metadata,
        "input_paths": paths,
    }


def _indexed_artifacts(job_dir, glob_pattern, index_pattern, label):
    indexed = defaultdict(list)
    for path in job_dir.glob(glob_pattern):
        match = re.search(index_pattern, path.name)
        if match:
            indexed[int(match.group(1))].append(path)
    if set(indexed) != set(range(5)) or any(len(paths) != 1 for paths in indexed.values()):
        raise ValueError(
            f"{job_dir} requires exactly one {label} file for each model index 0-4"
        )
    return [indexed[index][0] for index in range(5)]


def discover_saved_jobs(project_root, candidate_metadata):
    """Rediscover and validate every involved saved AlphaFold job path."""
    project_root = Path(project_root).resolve()
    candidate_ids = set(candidate_metadata)
    saved_paths = []
    consumed_hashes = {}
    for request_path in sorted(project_root.rglob("*_job_request.json"), key=lambda path: str(path)):
        if OUTPUT_DIRECTORY_NAME in request_path.parts:
            continue
        request = json.loads(request_path.read_text(encoding="utf-8"))
        job = request[0] if isinstance(request, list) else request
        candidate_id = candidate_id_from_request_name(str(job.get("name", "")))
        if candidate_id not in candidate_ids:
            continue
        if not isinstance(job.get("sequences"), list) or len(job["sequences"]) != 3:
            raise ValueError(f"Expected exact three-chain request layout in {request_path}")
        chains = []
        for entry in job["sequences"]:
            protein_chain = entry.get("proteinChain") if isinstance(entry, dict) else None
            if not isinstance(protein_chain, dict) or protein_chain.get("count", 1) != 1:
                raise ValueError(f"Expected three single protein chains in {request_path}")
            chains.append(str(protein_chain.get("sequence", "")))
        seeds = job.get("modelSeeds") or []
        if len(seeds) != 1:
            raise ValueError(f"Expected exactly one model seed in {request_path}")
        details = request_details(request)
        expected_peptide = candidate_metadata[candidate_id]["peptide"]
        if chains[2] != expected_peptide or details["requested_peptide"] != expected_peptide:
            raise ValueError(f"Requested peptide mismatch for {candidate_id} in {request_path}")
        job_dir = request_path.parent
        if len(list(job_dir.glob("*_job_request.json"))) != 1:
            raise ValueError(f"Expected exactly one job request in {job_dir}")
        model_paths = _indexed_artifacts(
            job_dir, "*_model_*.cif", r"_model_(\d+)\.cif$", "CIF"
        )
        full_data_paths = _indexed_artifacts(
            job_dir, "*_full_data_*.json", r"_full_data_(\d+)\.json$", "full-data"
        )
        confidence_paths = _indexed_artifacts(
            job_dir,
            "*_summary_confidences_*.json",
            r"_summary_confidences_(\d+)\.json$",
            "confidence",
        )
        model_hashes = [_sha256_file(path) for path in model_paths]
        request_components = {
            "candidate_id": candidate_id,
            "ordered_chain_sequences": chains,
            "model_seed": str(seeds[0]),
        }
        request_identity = "sha256:" + hashlib.sha256(
            json.dumps(request_components, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        input_paths = [request_path, *model_paths, *full_data_paths, *confidence_paths]
        for path in input_paths:
            consumed_hashes[str(path.resolve())] = _sha256_file(path)
        saved_paths.append({
            "path": str(job_dir.resolve()),
            "candidate_id": candidate_id,
            "request_name": str(details["request_name"]),
            "request_identity": request_identity,
            "seed": str(seeds[0]),
            "chain_sequences": chains,
            "requested_peptide": expected_peptide,
            "model_hashes": model_hashes,
            "model_paths": model_paths,
            "full_data_paths": full_data_paths,
            "confidence_paths": confidence_paths,
            "request_path": request_path,
        })

    found_candidates = {job["candidate_id"] for job in saved_paths}
    if found_candidates != candidate_ids:
        raise FileNotFoundError(
            f"Saved jobs are missing involved entities: {sorted(candidate_ids - found_candidates)}"
        )
    hla_layouts = {(job["chain_sequences"][0], job["chain_sequences"][1]) for job in saved_paths}
    if len(hla_layouts) != 1:
        raise ValueError("Involved jobs do not share one exact ordered HLA-DRA/DRB sequence layout")

    saved_paths.sort(key=lambda job: job["path"])
    deduplicated = deduplicate_jobs(saved_paths)
    retained_jobs = []
    for job in deduplicated["retained_jobs"]:
        identity = (job["request_identity"], tuple(job["model_hashes"]))
        duplicate_paths = sorted(deduplicated["duplicate_paths_by_canonical"][identity])
        signature = hashlib.sha256(
            json.dumps(job["model_hashes"], separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        canonical_job_id = f"{job['candidate_id']}|seed-{job['seed']}|sig-{signature[:16]}"
        parsed_models = []
        for index, (model_path, confidence_path) in enumerate(
            zip(job["model_paths"], job["confidence_paths"])
        ):
            model = parse_mmcif(model_path)
            observed_sequences = [sequence(model.get(chain, [])) for chain in ("A", "B", "C")]
            if set(model) != {"A", "B", "C"} or observed_sequences != job["chain_sequences"]:
                raise ValueError(
                    f"Exact A/B/C sequence integrity failed for {model_path}: {observed_sequences}"
                )
            confidence = json.loads(confidence_path.read_text(encoding="utf-8"))
            if "iptm" not in confidence:
                raise ValueError(f"ipTM annotation is absent from {confidence_path}")
            peptide_plddt = float(np.mean([residue_plddt(residue) for residue in model["C"]]))
            parsed_models.append({
                "model_index": index,
                "model_id": f"{canonical_job_id}|model-{index}",
                "model_path": str(model_path.resolve()),
                "parsed_model": model,
                "model_plddt": peptide_plddt,
                "model_iptm": float(confidence["iptm"]),
            })
        retained_jobs.append({
            **job,
            "canonical_job_id": canonical_job_id,
            "complete_model_signature": "sha256:" + signature,
            "duplicate_paths": duplicate_paths,
            "all_saved_paths": [job["path"], *duplicate_paths],
            "unique_job_status": (
                "deduplicated_identical_request_and_complete_signature"
                if duplicate_paths else "unique_saved_path"
            ),
            "sequence_integrity_status": "pass_exact_three_chain_request_and_cif_sequences",
            "models": parsed_models,
            "model_ids": [model["model_id"] for model in parsed_models],
        })
    retained_jobs.sort(key=lambda job: (job["candidate_id"], job["canonical_job_id"]))
    return retained_jobs, consumed_hashes


def build_live_lead_inputs(lead_specs, candidate_metadata, jobs):
    """Parse the retained jobs into pure inputs and recompute all pose distances."""
    jobs_by_candidate = defaultdict(list)
    for job in jobs:
        jobs_by_candidate[job["candidate_id"]].append(job)
    lead_inputs = []
    for spec in sorted(lead_specs, key=lambda item: item["rank"]):
        validate_lead_definition(spec)
        required_candidates = [
            spec["ebv_candidate_id"],
            spec["target_candidate_id"],
            *spec["control_ids"],
        ]
        for candidate_id in required_candidates:
            if not jobs_by_candidate[candidate_id]:
                raise ValueError(f"No retained unique job for {candidate_id}")

        def public_job(job):
            return {
                "candidate_id": job["candidate_id"],
                "job_id": job["canonical_job_id"],
                "model_ids": list(job["model_ids"]),
            }

        entities = {
            "ebv": [public_job(job) for job in jobs_by_candidate[spec["ebv_candidate_id"]]],
            "target": [public_job(job) for job in jobs_by_candidate[spec["target_candidate_id"]]],
            "controls": {
                control_id: [public_job(job) for job in jobs_by_candidate[control_id]]
                for control_id in spec["control_ids"]
            },
        }
        involved_jobs = [
            job for candidate_id in required_candidates for job in jobs_by_candidate[candidate_id]
        ]
        models = [
            {**model, "candidate_id": job["candidate_id"]}
            for job in involved_jobs
            for model in job["models"]
        ]
        models.sort(key=lambda model: model["model_id"])
        geometry_lookup = {}
        for left_index, left in enumerate(models):
            left_start = candidate_metadata[left["candidate_id"]]["core_start_1_based"]
            for right in models[left_index + 1:]:
                right_start = candidate_metadata[right["candidate_id"]]["core_start_1_based"]
                metrics = same_register_geometry(
                    left["parsed_model"],
                    right["parsed_model"],
                    left_start,
                    right_start,
                )
                geometry_lookup[(left["model_id"], right["model_id"])] = float(
                    metrics["candidate_exposed_ca_rmsd_A"]
                )
        confidence_by_model = {
            model["model_id"]: {
                "model_plddt": model["model_plddt"],
                "model_iptm": model["model_iptm"],
            }
            for model in models
        }
        lead_inputs.append({
            **spec,
            "entities": entities,
            "geometry_lookup": geometry_lookup,
            "confidence_by_model": confidence_by_model,
        })
    return lead_inputs


def build_identity_manifest(lead_specs, jobs):
    """Build one deterministic manifest row per retained unique job."""
    placement = {}
    for spec in lead_specs:
        placement[spec["ebv_candidate_id"]] = (spec["rank"], "ebv")
        placement[spec["target_candidate_id"]] = (spec["rank"], "target")
        for control_id in spec["control_ids"]:
            placement[control_id] = (spec["rank"], "control")
    role_order = {"ebv": 0, "target": 1, "control": 2}
    rows = []
    for job in jobs:
        rank, role = placement[job["candidate_id"]]
        rows.append({
            "lead_rank": rank,
            "entity_role": role,
            "candidate_id": job["candidate_id"],
            "canonical_job_id": job["canonical_job_id"],
            "canonical_job_path": job["path"],
            "request_identity": job["request_identity"],
            "model_seed": job["seed"],
            "requested_dra_sha256": hashlib.sha256(job["chain_sequences"][0].encode("utf-8")).hexdigest(),
            "requested_drb_sha256": hashlib.sha256(job["chain_sequences"][1].encode("utf-8")).hexdigest(),
            "requested_peptide": job["requested_peptide"],
            "complete_model_signature": job["complete_model_signature"],
            "model_0_sha256": job["model_hashes"][0],
            "model_1_sha256": job["model_hashes"][1],
            "model_2_sha256": job["model_hashes"][2],
            "model_3_sha256": job["model_hashes"][3],
            "model_4_sha256": job["model_hashes"][4],
            "ordered_model_sha256": ";".join(job["model_hashes"]),
            "duplicate_path_count": len(job["duplicate_paths"]),
            "duplicate_paths": " || ".join(job["duplicate_paths"]),
            "all_saved_path_count": len(job["all_saved_paths"]),
            "all_saved_paths": " || ".join(job["all_saved_paths"]),
            "unique_job_status": job["unique_job_status"],
            "exact_sequence_integrity_status": job["sequence_integrity_status"],
        })
    rows.sort(
        key=lambda row: (
            row["lead_rank"], role_order[row["entity_role"]], row["candidate_id"], row["canonical_job_id"]
        )
    )
    return rows


def _format_metric(value):
    return f"{float(value):.3f}"


def build_findings(manifest_rows, tables, project_root):
    rank_rows = {row["lead_rank"]: row for row in tables["control_rank_and_leave_one_out"]}
    bootstrap_rows = {row["lead_rank"]: row for row in tables["technical_bootstrap_summary"]}
    rank1, rank2 = rank_rows[1], rank_rows[2]
    boot1, boot2 = bootstrap_rows[1], bootstrap_rows[2]
    duplicate_paths = sum(int(row["duplicate_path_count"]) for row in manifest_rows)
    return f"""# Lead-focused structural robustness findings

## Scope and methods

This audit keeps the frozen discovery ranking and analyzes the two leads separately. Rank 1 uses only its three strict primary controls. Rank 2 uses only `length_sensitivity_exact_bin_pm7` and is supplemental/sensitivity-only. The two layers were not pooled, averaged together, or given equal evidentiary weight.

Every saved path for the ten involved entities was rediscovered from its AlphaFold request and required exactly five CIFs, five full-data files, and five confidence files. Request identity used candidate ID, the exact ordered three-chain sequences, and model seed. Jobs were deduplicated only when that request identity and all five ordered CIF SHA-256 hashes matched. The manifest retains {len(manifest_rows)} distinct jobs and records {duplicate_paths} identical saved-path copies.

For every retained model, the exact frozen P1-P9 core had one unambiguous placement. CIFs were parsed with the existing project parser. Exposed-position RMSD used P2/P3/P5/P7/P8 after fitting the HLA groove. Confidence was retained only as an annotation and did not select models or influence the 2.0-A complete-linkage pose clusters.

The hierarchical bootstrap used 10,000 iterations per lead and seed 20260815. It resampled unique jobs and then the five technical models within jobs, retaining all three controls at equal top-level weight.

## Results

### Rank 1: strict primary-control lead

- Pair: `{rank1['pair_id']}`
- Classification: `{rank1['classification']}`
- Target median: {_format_metric(rank1['overall_target_median_A'])} A
- Equal-weight background median: {_format_metric(rank1['overall_background_median_A'])} A
- Background-minus-target delta: {_format_metric(rank1['background_minus_target_delta_A'])} A
- Leave-one-control-out delta range: {_format_metric(rank1['leave_one_out_delta_min_A'])} to {_format_metric(rank1['leave_one_out_delta_max_A'])} A
- Exploratory target rank: {rank1['target_rank']} of {int(rank1['control_count']) + 1}; empirical tail fraction {float(rank1['exploratory_empirical_tail_fraction']):.2f}
- Technical-stability delta interval: {_format_metric(boot1['delta_percentile_2_5_A'])} to {_format_metric(boot1['delta_percentile_97_5_A'])} A

### Rank 2: supplemental length sensitivity

- Pair: `{rank2['pair_id']}`
- Classification: `{rank2['classification']}`
- Target median: {_format_metric(rank2['overall_target_median_A'])} A
- Equal-weight background median: {_format_metric(rank2['overall_background_median_A'])} A
- Background-minus-target delta: {_format_metric(rank2['background_minus_target_delta_A'])} A
- Leave-one-control-out delta range: {_format_metric(rank2['leave_one_out_delta_min_A'])} to {_format_metric(rank2['leave_one_out_delta_max_A'])} A
- Exploratory target rank: {rank2['target_rank']} of {int(rank2['control_count']) + 1}; empirical tail fraction {float(rank2['exploratory_empirical_tail_fraction']):.2f}
- Technical-stability delta interval: {_format_metric(boot2['delta_percentile_2_5_A'])} to {_format_metric(boot2['delta_percentile_97_5_A'])} A

The empirical tail fractions are exploratory ranks, not p-values. With three controls, 0.25 is the smallest possible fraction. The bootstrap intervals quantify technical stability across saved AlphaFold jobs/models only; they are not biological replication or p-values.

## Limitations and claim boundary

This is descriptive computational pMHC geometry. AlphaFold jobs and models are technical samples, not biological replicates. The small frozen control sets limit empirical resolution. Rank 2 has a deliberate peptide-length mismatch and remains sensitivity-only. These outputs do not establish peptide presentation, TCR binding, activation, cross-reactivity, molecular mimicry, or an MS mechanism.

## Exact reproduce command

From the project root, with the output folder absent (the generator refuses overwrite):

```bash
PYTHONPATH=src python3 -m lead_focused_robustness --project-root "$PWD"
```

Project audited: `{Path(project_root).resolve()}`
"""


def collect_input_checksums(project_root, table_paths, consumed_hashes):
    """Return the deduplicated checksums for every consumed live-project input."""
    project_root = Path(project_root).resolve()
    rows = []
    combined = {
        str(path.resolve()): _sha256_file(path)
        for path in table_paths.values()
    }
    combined.update(consumed_hashes)
    table_set = {str(path.resolve()) for path in table_paths.values()}
    for path_string, digest in sorted(combined.items()):
        path = Path(path_string)
        try:
            display_path = str(path.relative_to(project_root))
        except ValueError:
            display_path = str(path)
        rows.append({
            "input_path": display_path,
            "input_role": "frozen_project_table" if path_string in table_set else "saved_job_artifact",
            "sha256": digest,
        })
    return rows


def generate_live_audit(project_root):
    """Run the fixed live-project audit with 10,000 iterations per lead."""
    project_root = Path(project_root).resolve()
    output_dir = project_root / "processed" / OUTPUT_DIRECTORY_NAME
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit folder: {output_dir}")
    metadata = load_fixed_lead_metadata(project_root)
    jobs, consumed_hashes = discover_saved_jobs(
        project_root, metadata["candidate_metadata"]
    )
    lead_inputs = build_live_lead_inputs(
        metadata["lead_specs"], metadata["candidate_metadata"], jobs
    )
    tables = build_audit_tables(lead_inputs, iterations=10_000, seed=20260815)
    expected_replicates = len(metadata["lead_specs"]) * 10_000
    if len(tables["technical_bootstrap_replicates"]) != expected_replicates:
        raise AssertionError(
            f"Audit requires exactly {expected_replicates} bootstrap replicate rows"
        )
    manifest_rows = build_identity_manifest(metadata["lead_specs"], jobs)
    rank_summaries = {
        row["lead_rank"]: row for row in tables["control_rank_and_leave_one_out"]
    }
    bootstrap_summaries = {
        row["lead_rank"]: row for row in tables["technical_bootstrap_summary"]
    }
    control_summaries = defaultdict(list)
    for row in tables["per_control_geometry_summary"]:
        control_summaries[row["lead_rank"]].append(row)
    target_job_values = defaultdict(list)
    for row in tables["job_pair_stability"]:
        if row["lead_rank"] == 2 and row["entity_role"] == "target":
            target_job_values[row["comparison_job_id"]].append(row["exposed_rmsd_median_A"])
    rank2_target_rows = [
        {
            "comparison_job_id": job_id,
            "exposed_rmsd_median_A": float(median(values)),
        }
        for job_id, values in sorted(target_job_values.items())
    ]
    rank1_svg = render_rank1_svg(
        rank_summaries[1], control_summaries[1], bootstrap_summaries[1]
    )
    rank2_svg = render_rank2_svg(
        rank_summaries[2],
        control_summaries[2],
        bootstrap_summaries[2],
        rank2_target_rows,
    )
    findings = build_findings(manifest_rows, tables, project_root)
    checksum_rows = collect_input_checksums(
        project_root, metadata["input_paths"], consumed_hashes
    )
    write_audit_outputs(
        output_dir,
        manifest_rows,
        tables,
        rank1_svg,
        rank2_svg,
        findings,
        checksum_rows,
    )
    return output_dir


def deduplicate_jobs(jobs):
    """Retain the first job for each exact request-plus-model-hash identity."""
    retained_jobs = []
    duplicate_paths_by_canonical = {}
    canonical_by_identity = {}
    for job in jobs:
        if len(job["model_hashes"]) != 5:
            raise ValueError("Each job requires ordered hashes for all five model files")
        identity = (job["request_identity"], tuple(job["model_hashes"]))
        if identity not in canonical_by_identity:
            canonical_by_identity[identity] = job
            retained_jobs.append(job)
            duplicate_paths_by_canonical[identity] = []
        else:
            duplicate_paths_by_canonical[identity].append(job["path"])
    return {
        "retained_jobs": retained_jobs,
        "duplicate_paths_by_canonical": duplicate_paths_by_canonical,
    }


def summarize_lead(target_job_medians, control_job_medians):
    """Summarize a target against controls without allowing job count to add weight."""
    if not target_job_medians or not control_job_medians:
        raise ValueError("A target and at least one control are required")
    if any(not job_medians for job_medians in control_job_medians.values()):
        raise ValueError("Every control requires at least one unique job median")
    target_median = float(median(target_job_medians))
    control_medians = {
        identifier: float(median(job_medians))
        for identifier, job_medians in sorted(control_job_medians.items())
    }
    background_median = float(median(control_medians.values()))
    return {
        "target_median": target_median,
        "target_job_medians": list(target_job_medians),
        "control_medians": control_medians,
        "background_median": background_median,
        "background_minus_target_delta": background_median - target_median,
        "leave_one_control_out_deltas": {
            identifier: float(median([value for other, value in control_medians.items() if other != identifier])) - target_median
            if len(control_medians) > 1
            else None
            for identifier in control_medians
        },
    }


def empirical_target_rank(target_identifier, target_median, control_medians):
    """Return a stable descriptive rank and conservative lower-tail fraction."""
    ordered = sorted(
        [(float(target_median), target_identifier)]
        + [(float(value), identifier) for identifier, value in control_medians.items()],
        key=lambda pair: (pair[0], pair[1]),
    )
    return {
        "sorted_identifiers": [identifier for _, identifier in ordered],
        "target_rank": [identifier for _, identifier in ordered].index(target_identifier) + 1,
        "one_sided_tail_fraction": (
            1 + sum(value <= target_median for value in control_medians.values())
        ) / (len(control_medians) + 1),
    }


def cluster_complete_linkage(model_ids, rmsd_matrix, threshold=2.0, confidence_by_model=None):
    """Cluster a symmetric RMSD matrix with deterministic complete linkage.

    ``confidence_by_model`` is deliberately accepted only as an annotation;
    it never contributes to distance or tie-breaking.
    """
    del confidence_by_model
    if len(set(model_ids)) != len(model_ids):
        raise ValueError("Model IDs must be unique")
    if len(rmsd_matrix) != len(model_ids) or any(len(row) != len(model_ids) for row in rmsd_matrix):
        raise ValueError("RMSD matrix dimensions must match model IDs")
    original_index = {model_id: index for index, model_id in enumerate(model_ids)}
    for left in model_ids:
        for right in model_ids:
            if abs(float(rmsd_matrix[original_index[left]][original_index[right]]) - float(rmsd_matrix[original_index[right]][original_index[left]])) > 1e-12:
                raise ValueError("RMSD matrix must be symmetric")
    clusters = [(model_id,) for model_id in sorted(model_ids)]

    def complete_distance(left_cluster, right_cluster):
        return max(
            float(rmsd_matrix[original_index[left]][original_index[right]])
            for left in left_cluster
            for right in right_cluster
        )

    while len(clusters) > 1:
        candidates = []
        for left_index in range(len(clusters)):
            for right_index in range(left_index + 1, len(clusters)):
                left, right = clusters[left_index], clusters[right_index]
                candidates.append((complete_distance(left, right), left, right, left_index, right_index))
        distance, left, right, left_index, right_index = min(candidates, key=lambda item: (item[0], item[1], item[2]))
        if distance > threshold:
            break
        merged = tuple(sorted(left + right))
        clusters = [cluster for index, cluster in enumerate(clusters) if index not in (left_index, right_index)]
        clusters.append(merged)
        clusters.sort()

    return {
        model_id: label
        for label, cluster in enumerate(sorted(clusters), start=1)
        for model_id in cluster
    }


def classify_separation(summary, bootstrap_summary, rank=1):
    """Classify descriptive robustness without treating it as biological proof."""
    delta = summary["background_minus_target_delta"]
    if delta <= 0:
        label = "no_positive_separation"
    else:
        leave_one_out = summary["leave_one_control_out_deltas"].values()
        all_leave_one_out_positive = bool(leave_one_out) and all(
            value is not None and value > 0 for value in leave_one_out
        )
        all_target_jobs_below_background = all(
            value < summary["background_median"] for value in summary["target_job_medians"]
        )
        if (
            all_leave_one_out_positive
            and all_target_jobs_below_background
            and bootstrap_summary["percentile_2_5"] > 0
        ):
            label = "consistent_positive"
        else:
            label = "mixed_positive"
    return "length_sensitivity_only__" + label if rank == 2 else label


def require_single_analysis_layer(records):
    """Reject a collection that would pool primary and sensitivity analyses."""
    layers = {record["analysis_layer"] for record in records}
    if len(layers) != 1:
        raise ValueError("Primary and sensitivity analysis layers cannot be pooled")
    return next(iter(layers))


def hierarchical_technical_bootstrap(
    entities, geometry_lookup, iterations=10_000, seed=20260815
):
    """Run the fixed 10,000-replicate job-then-model bootstrap.

    ``entities`` has ``ebv`` and ``target`` lists of unique jobs plus a
    ``controls`` mapping.  Every job has a stable ``job_id`` and exactly five
    ``model_ids``. ``geometry_lookup`` maps ordered model-ID pairs to RMSD;
    reverse-pair lookup is accepted for symmetric source tables.
    """
    ebv_jobs = entities["ebv"]
    target_jobs = entities["target"]
    controls = entities["controls"]
    if not ebv_jobs or not target_jobs or not controls:
        raise ValueError("EBV, target, and at least one control are required")

    def validate_jobs(jobs):
        if not jobs:
            raise ValueError("Every entity requires at least one unique job")
        if len({job["job_id"] for job in jobs}) != len(jobs):
            raise ValueError("Jobs must already be unique before bootstrap")
        for job in jobs:
            if len(job["model_ids"]) != 5:
                raise ValueError("Every job must expose exactly five model IDs")

    validate_jobs(ebv_jobs)
    validate_jobs(target_jobs)
    for jobs in controls.values():
        validate_jobs(jobs)

    def lookup(left, right):
        if (left, right) in geometry_lookup:
            return float(geometry_lookup[(left, right)])
        if (right, left) in geometry_lookup:
            return float(geometry_lookup[(right, left)])
        raise KeyError(f"No geometry for model pair {left!r}, {right!r}")

    def sample_jobs(jobs, rng):
        sampled = []
        for job_index in rng.integers(0, len(jobs), size=len(jobs)):
            job = jobs[int(job_index)]
            model_indices = rng.integers(0, 5, size=5)
            sampled.append([job["model_ids"][int(index)] for index in model_indices])
        return sampled

    def nested_job_median(left_jobs, right_jobs):
        comparison_job_medians = []
        for right_job in right_jobs:
            ebv_job_pair_medians = [
                float(np.median([
                    lookup(left_model, right_model)
                    for left_model in left_job
                    for right_model in right_job
                ]))
                for left_job in left_jobs
            ]
            comparison_job_medians.append(float(np.median(ebv_job_pair_medians)))
        return float(np.median(comparison_job_medians))

    if not isinstance(iterations, int) or iterations <= 0:
        raise ValueError("Bootstrap iterations must be a positive integer")
    rng = np.random.default_rng(seed)
    replicates = []
    replicate_rows = []
    for iteration in range(1, iterations + 1):
        sampled_ebv = sample_jobs(ebv_jobs, rng)
        target_median = nested_job_median(
            sampled_ebv, sample_jobs(target_jobs, rng)
        )
        control_medians = [
            nested_job_median(sampled_ebv, sample_jobs(jobs, rng))
            for _, jobs in sorted(controls.items())
        ]
        background_median = float(np.median(control_medians))
        delta = background_median - target_median
        replicates.append(delta)
        replicate_rows.append({
            "iteration": iteration,
            "target_median_A": target_median,
            "equal_weight_background_median_A": background_median,
            "delta_A": delta,
        })

    return {
        "replicates": replicates,
        "replicate_rows": replicate_rows,
        "median": float(np.median(replicates)),
        "percentile_2_5": float(np.percentile(replicates, 2.5)),
        "percentile_97_5": float(np.percentile(replicates, 97.5)),
        "fraction_positive": float(np.mean(np.asarray(replicates) > 0)),
    }


CLAIM_BOUNDARY = (
    "Descriptive computational pMHC geometry and technical stability only; "
    "not evidence of presentation, TCR binding, activation, cross-reactivity, "
    "molecular mimicry, or MS mechanism."
)
TECHNICAL_INTERVAL_LABEL = (
    "Technical-stability interval across saved AlphaFold jobs/models only; "
    "not biological replication and not a p-value."
)


def _geometry_value(geometry_lookup, left, right):
    if (left, right) in geometry_lookup:
        return float(geometry_lookup[(left, right)])
    if (right, left) in geometry_lookup:
        return float(geometry_lookup[(right, left)])
    raise KeyError(f"No geometry for model pair {left!r}, {right!r}")


def _job_pair_summary(left_job, right_job, geometry_lookup):
    values = [
        _geometry_value(geometry_lookup, left_model, right_model)
        for left_model in left_job["model_ids"]
        for right_model in right_job["model_ids"]
    ]
    return {
        "technical_comparison_count": len(values),
        "exposed_rmsd_median_A": float(np.median(values)),
        "below_2_A_fraction": sum(value < 2.0 for value in values) / len(values),
    }


def _relative_sign(value, reference):
    if value < reference:
        return "below_equal_weight_background_median"
    if value > reference:
        return "above_equal_weight_background_median"
    return "equal_to_equal_weight_background_median"


def _unique_comparison_job_medians(rows):
    """Collapse EBV-job pair medians to one value per comparison-side job."""
    values_by_job = defaultdict(list)
    for row in rows:
        values_by_job[row["comparison_job_id"]].append(row["exposed_rmsd_median_A"])
    return [
        float(median(values_by_job[job_id]))
        for job_id in sorted(values_by_job)
    ]


def build_audit_tables(leads, iterations=10_000, seed=20260815):
    """Build deterministic audit tables from injected, already parsed lead inputs."""
    tables = {
        "per_control_geometry_summary": [],
        "job_pair_stability": [],
        "control_rank_and_leave_one_out": [],
        "technical_bootstrap_replicates": [],
        "technical_bootstrap_summary": [],
        "pose_cluster_membership": [],
    }
    seen_ranks = set()
    role_order = {"ebv": 0, "target": 1, "control": 2}
    for lead in sorted(leads, key=lambda item: int(item["rank"])):
        validate_lead_definition(lead)
        rank = int(lead["rank"])
        if rank in seen_ranks:
            raise ValueError(f"Duplicate lead rank {rank}")
        seen_ranks.add(rank)
        entities = lead["entities"]
        geometry_lookup = lead["geometry_lookup"]
        control_ids = sorted(lead["control_ids"])
        if set(entities["controls"]) != set(control_ids):
            raise ValueError("Entity controls do not match the frozen lead definition")

        raw_job_pairs = []
        for ebv_job in sorted(entities["ebv"], key=lambda item: item["job_id"]):
            for target_job in sorted(entities["target"], key=lambda item: item["job_id"]):
                metrics = _job_pair_summary(ebv_job, target_job, geometry_lookup)
                raw_job_pairs.append({
                    "lead_rank": rank,
                    "pair_id": lead["pair_id"],
                    "analysis_layer": lead["analysis_layer"],
                    "ebv_job_id": ebv_job["job_id"],
                    "comparison_candidate_id": lead["target_candidate_id"],
                    "comparison_job_id": target_job["job_id"],
                    "entity_role": "target",
                    **metrics,
                })

        for control_id in control_ids:
            for ebv_job in sorted(entities["ebv"], key=lambda item: item["job_id"]):
                for control_job in sorted(
                    entities["controls"][control_id], key=lambda item: item["job_id"]
                ):
                    metrics = _job_pair_summary(ebv_job, control_job, geometry_lookup)
                    raw_job_pairs.append({
                        "lead_rank": rank,
                        "pair_id": lead["pair_id"],
                        "analysis_layer": lead["analysis_layer"],
                        "ebv_job_id": ebv_job["job_id"],
                        "comparison_candidate_id": control_id,
                        "comparison_job_id": control_job["job_id"],
                        "entity_role": "control",
                        **metrics,
                    })

        target_job_medians = _unique_comparison_job_medians([
            row for row in raw_job_pairs if row["entity_role"] == "target"
        ])
        control_job_medians = {
            control_id: _unique_comparison_job_medians([
                row for row in raw_job_pairs
                if row["entity_role"] == "control"
                and row["comparison_candidate_id"] == control_id
            ])
            for control_id in control_ids
        }
        summary = summarize_lead(target_job_medians, control_job_medians)
        for row in raw_job_pairs:
            tables["job_pair_stability"].append({
                **row,
                "relative_to_equal_weight_background_median": _relative_sign(
                    row["exposed_rmsd_median_A"], summary["background_median"]
                ),
                "claim_boundary": CLAIM_BOUNDARY,
            })

        for control_id in control_ids:
            control_rows = [
                row for row in raw_job_pairs
                if row["entity_role"] == "control"
                and row["comparison_candidate_id"] == control_id
            ]
            control_median = summary["control_medians"][control_id]
            tables["per_control_geometry_summary"].append({
                "lead_rank": rank,
                "pair_id": lead["pair_id"],
                "analysis_layer": lead["analysis_layer"],
                "control_candidate_id": control_id,
                "technical_geometry_count": sum(
                    row["technical_comparison_count"] for row in control_rows
                ),
                "job_count": len(entities["controls"][control_id]),
                "control_median_A": control_median,
                "target_median_A": summary["target_median"],
                "control_minus_target_delta_A": control_median - summary["target_median"],
                "claim_boundary": CLAIM_BOUNDARY,
            })

        rank_result = empirical_target_rank(
            lead["target_candidate_id"],
            summary["target_median"],
            summary["control_medians"],
        )
        bootstrap = hierarchical_technical_bootstrap(
            entities, geometry_lookup, iterations=iterations, seed=seed
        )
        classification = classify_separation(summary, bootstrap, rank=rank)
        leave_one_out = summary["leave_one_control_out_deltas"]
        leave_values = [value for value in leave_one_out.values() if value is not None]
        tables["control_rank_and_leave_one_out"].append({
            "lead_rank": rank,
            "pair_id": lead["pair_id"],
            "analysis_layer": lead["analysis_layer"],
            "target_rank": rank_result["target_rank"],
            "control_count": len(control_ids),
            "exploratory_empirical_tail_fraction": rank_result["one_sided_tail_fraction"],
            "overall_background_median_A": summary["background_median"],
            "overall_target_median_A": summary["target_median"],
            "background_minus_target_delta_A": summary["background_minus_target_delta"],
            "leave_one_out_deltas_A": ";".join(
                f"{identifier}={leave_one_out[identifier]:.6f}"
                for identifier in sorted(leave_one_out)
            ),
            "leave_one_out_delta_min_A": min(leave_values),
            "leave_one_out_delta_max_A": max(leave_values),
            "classification": classification,
            "inference_label": (
                "Exploratory empirical tail fraction, not a p-value; with three controls "
                "the minimum possible fraction is 0.25."
            ),
            "claim_boundary": CLAIM_BOUNDARY,
        })
        for row in bootstrap["replicate_rows"]:
            tables["technical_bootstrap_replicates"].append({
                "lead_rank": rank,
                "pair_id": lead["pair_id"],
                "analysis_layer": lead["analysis_layer"],
                **row,
                "claim_boundary": CLAIM_BOUNDARY,
            })
        tables["technical_bootstrap_summary"].append({
            "lead_rank": rank,
            "pair_id": lead["pair_id"],
            "analysis_layer": lead["analysis_layer"],
            "iterations": iterations,
            "seed": seed,
            "delta_median_A": bootstrap["median"],
            "delta_percentile_2_5_A": bootstrap["percentile_2_5"],
            "delta_percentile_97_5_A": bootstrap["percentile_97_5"],
            "fraction_positive": bootstrap["fraction_positive"],
            "inference_label": TECHNICAL_INTERVAL_LABEL,
            "claim_boundary": CLAIM_BOUNDARY,
        })

        model_metadata = {}
        entity_groups = [
            ("ebv", lead["ebv_candidate_id"], entities["ebv"]),
            ("target", lead["target_candidate_id"], entities["target"]),
        ] + [
            ("control", control_id, entities["controls"][control_id])
            for control_id in control_ids
        ]
        for role, candidate_id, jobs in entity_groups:
            for job in sorted(jobs, key=lambda item: item["job_id"]):
                if len(job["model_ids"]) != 5:
                    raise ValueError("Every pose-cluster job requires exactly five models")
                for model_id in job["model_ids"]:
                    if model_id in model_metadata:
                        raise ValueError("Model IDs must be unique across a lead")
                    model_metadata[model_id] = {
                        "entity_role": role,
                        "candidate_id": candidate_id,
                        "job_id": job["job_id"],
                    }
        model_ids = sorted(model_metadata)
        distance_matrix = [
            [
                0.0 if left == right else _geometry_value(geometry_lookup, left, right)
                for right in model_ids
            ]
            for left in model_ids
        ]
        labels = cluster_complete_linkage(
            model_ids,
            distance_matrix,
            threshold=2.0,
            confidence_by_model=lead.get("confidence_by_model", {}),
        )
        models_by_cluster = {
            cluster_id: [model for model in model_ids if labels[model] == cluster_id]
            for cluster_id in sorted(set(labels.values()))
        }
        for model_id in model_ids:
            metadata = model_metadata[model_id]
            cluster_models = models_by_cluster[labels[model_id]]
            confidence = lead.get("confidence_by_model", {}).get(model_id, {})
            tables["pose_cluster_membership"].append({
                "lead_rank": rank,
                "pair_id": lead["pair_id"],
                "analysis_layer": lead["analysis_layer"],
                "entity_role": metadata["entity_role"],
                "candidate_id": metadata["candidate_id"],
                "job_id": metadata["job_id"],
                "model_id": model_id,
                "cluster_id": labels[model_id],
                "cluster_size": len(cluster_models),
                "distinct_jobs_in_cluster": len({
                    model_metadata[model]["job_id"] for model in cluster_models
                }),
                "model_plddt": confidence.get("model_plddt", ""),
                "model_iptm": confidence.get("model_iptm", ""),
                "clustering_method_note": (
                    "Complete linkage at 2.0 A; confidence annotations did not affect clustering."
                ),
                "claim_boundary": CLAIM_BOUNDARY,
            })

    tables["per_control_geometry_summary"].sort(
        key=lambda row: (row["lead_rank"], row["control_candidate_id"])
    )
    tables["job_pair_stability"].sort(
        key=lambda row: (
            row["lead_rank"],
            role_order[row["entity_role"]],
            row["comparison_candidate_id"],
            row["ebv_job_id"],
            row["comparison_job_id"],
        )
    )
    tables["control_rank_and_leave_one_out"].sort(key=lambda row: row["lead_rank"])
    tables["technical_bootstrap_replicates"].sort(
        key=lambda row: (row["lead_rank"], row["iteration"])
    )
    tables["technical_bootstrap_summary"].sort(key=lambda row: row["lead_rank"])
    tables["pose_cluster_membership"].sort(
        key=lambda row: (
            row["lead_rank"],
            role_order[row["entity_role"]],
            row["candidate_id"],
            row["job_id"],
            row["model_id"],
        )
    )
    return tables


def _svg_number(value):
    return f"{float(value):.2f}"


def _svg_x(value, maximum, left=250.0, width=600.0):
    return left + (float(value) / maximum) * width


def render_rank1_svg(summary, controls, bootstrap, title="Rank 1 strict-primary robustness"):
    """Render an accessible SVG for the primary strict-control lead."""
    values = [
        float(summary["overall_target_median_A"]),
        float(summary["overall_background_median_A"]),
        *[float(row["control_median_A"]) for row in controls],
    ]
    maximum = max(values + [1.0]) * 1.15
    ticks = [maximum * index / 4 for index in range(5)]
    rows = [
        ("Target human", float(summary["overall_target_median_A"]), "#0072B2"),
        *[
            (str(row["control_candidate_id"]), float(row["control_median_A"]), "#D55E00")
            for row in sorted(controls, key=lambda item: str(item["control_candidate_id"]))
        ],
        ("Equal-weight background", float(summary["overall_background_median_A"]), "#009E73"),
    ]
    plot_rows = []
    for index, (label, value, color) in enumerate(rows):
        y = 180 + index * 54
        x = _svg_x(value, maximum)
        plot_rows.append(
            f'<text x="235" y="{y + 5}" text-anchor="end" class="label">{escape(label)}</text>'
            f'<line x1="250" y1="{y}" x2="850" y2="{y}" class="guide"/>'
            f'<circle cx="{x:.1f}" cy="{y}" r="8" fill="{color}"/>'
            f'<text x="{x + 14:.1f}" y="{y + 5}" class="value">{_svg_number(value)} A</text>'
        )
    tick_markup = "".join(
        f'<line x1="{_svg_x(value, maximum):.1f}" y1="145" x2="{_svg_x(value, maximum):.1f}" y2="420" class="grid"/>'
        f'<text x="{_svg_x(value, maximum):.1f}" y="445" text-anchor="middle" class="axis">{_svg_number(value)}</text>'
        for value in ticks
    )
    description = (
        "Rank 1 target and three strict-control medians, leave-one-control-out "
        "range, exploratory empirical rank, and technical bootstrap interval."
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="960" height="620" viewBox="0 0 960 620" role="img" aria-labelledby="title desc">
<title id="title">{escape(title)}</title>
<desc id="desc">{escape(description)}</desc>
<style>
text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill: #172B4D; }}
.title {{ font-size: 25px; font-weight: 700; }} .subtitle {{ font-size: 15px; fill: #42526E; }}
.label {{ font-size: 13px; }} .value {{ font-size: 13px; font-weight: 600; }} .axis {{ font-size: 12px; fill: #5E6C84; }}
.guide {{ stroke: #C1C7D0; stroke-width: 1; }} .grid {{ stroke: #DFE1E6; stroke-width: 1; }}
.callout {{ fill: #F4F5F7; stroke: #B3BAC5; }} .callout-title {{ font-size: 14px; font-weight: 700; }} .callout-text {{ font-size: 13px; }}
</style>
<rect width="960" height="620" fill="#FFFFFF"/>
<text x="48" y="50" class="title">{escape(title)}</text>
<text x="48" y="78" class="subtitle">Exposed-position RMSD after HLA-groove fitting; lower values are closer.</text>
{tick_markup}
{''.join(plot_rows)}
<text x="550" y="472" text-anchor="middle" class="axis">Exposed-position RMSD (A)</text>
<rect x="48" y="500" width="270" height="88" rx="8" class="callout"/>
<text x="64" y="524" class="callout-title">Exploratory empirical rank</text>
<text x="64" y="550" class="callout-text">Target rank: {int(summary['target_rank'])} of {int(summary['control_count']) + 1}</text>
<text x="64" y="574" class="callout-text">Minimum tail fraction with 3 controls: 0.25</text>
<rect x="334" y="500" width="290" height="88" rx="8" class="callout"/>
<text x="350" y="524" class="callout-title">Leave-one-control-out delta range</text>
<text x="350" y="553" class="callout-text">{_svg_number(summary['leave_one_out_delta_min_A'])} to {_svg_number(summary['leave_one_out_delta_max_A'])} A</text>
<rect x="640" y="500" width="272" height="88" rx="8" class="callout"/>
<text x="656" y="524" class="callout-title">Technical-stability interval (not a p-value)</text>
<text x="656" y="553" class="callout-text">{_svg_number(bootstrap['delta_percentile_2_5_A'])} to {_svg_number(bootstrap['delta_percentile_97_5_A'])} A</text>
</svg>
'''


def render_rank2_svg(summary, controls, bootstrap, target_job_rows):
    """Render a clearly supplemental SVG for rank 2 job dependence."""
    target_rows = sorted(target_job_rows, key=lambda row: str(row["comparison_job_id"]))
    plotted = [
        (
            f"Target human job {index}",
            float(row["exposed_rmsd_median_A"]),
            "#0072B2",
            str(row["comparison_job_id"]),
        )
        for index, row in enumerate(target_rows, start=1)
    ] + [
        (
            "Control " + str(row["control_candidate_id"]).removeprefix("HUMAN_BACKGROUND_"),
            float(row["control_median_A"]),
            "#D55E00",
            str(row["control_candidate_id"]),
        )
        for row in sorted(controls, key=lambda item: str(item["control_candidate_id"]))
    ]
    maximum = max([value for _, value, _, _ in plotted] + [1.0]) * 1.15
    ticks = [maximum * index / 4 for index in range(5)]
    plot_rows = []
    for index, (label, value, color, full_identifier) in enumerate(plotted):
        y = 185 + index * 52
        x = _svg_x(value, maximum)
        plot_rows.append(
            f'<g><title>{escape(full_identifier)}</title>'
            f'<text x="235" y="{y + 5}" text-anchor="end" class="label">{escape(label)}</text>'
            f'<line x1="250" y1="{y}" x2="850" y2="{y}" class="guide"/>'
            f'<circle cx="{x:.1f}" cy="{y}" r="8" fill="{color}"/>'
            f'<text x="{x + 14:.1f}" y="{y + 5}" class="value">{_svg_number(value)} A</text></g>'
        )
    tick_markup = "".join(
        f'<line x1="{_svg_x(value, maximum):.1f}" y1="150" x2="{_svg_x(value, maximum):.1f}" y2="430" class="grid"/>'
        f'<text x="{_svg_x(value, maximum):.1f}" y="455" text-anchor="middle" class="axis">{_svg_number(value)}</text>'
        for value in ticks
    )
    title = "Rank 2 target-job dependence"
    description = (
        "Supplemental length-sensitivity-only comparison of target-human job "
        "medians, three control medians, and a technical bootstrap interval."
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="960" height="620" viewBox="0 0 960 620" role="img" aria-labelledby="title desc">
<title id="title">{escape(title)}</title>
<desc id="desc">{escape(description)}</desc>
<style>
text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill: #172B4D; }}
.banner {{ fill: #FFF3CD; stroke: #FFAB00; }} .banner-text {{ font-size: 14px; font-weight: 700; fill: #7A4F01; }}
.title {{ font-size: 25px; font-weight: 700; }} .subtitle {{ font-size: 15px; fill: #42526E; }}
.label {{ font-size: 12px; }} .value {{ font-size: 13px; font-weight: 600; }} .axis {{ font-size: 12px; fill: #5E6C84; }}
.guide {{ stroke: #C1C7D0; stroke-width: 1; }} .grid {{ stroke: #DFE1E6; stroke-width: 1; }}
.callout {{ fill: #F4F5F7; stroke: #B3BAC5; }} .callout-title {{ font-size: 14px; font-weight: 700; }} .callout-text {{ font-size: 13px; }}
</style>
<rect width="960" height="620" fill="#FFFFFF"/>
<rect x="48" y="24" width="864" height="34" rx="6" class="banner"/>
<text x="480" y="46" text-anchor="middle" class="banner-text">Supplemental / length-sensitivity-only</text>
<text x="48" y="94" class="title">{escape(title)}</text>
<text x="48" y="120" class="subtitle">The 32-aa target is not given the evidentiary weight of rank 1.</text>
{tick_markup}
{''.join(plot_rows)}
<text x="550" y="482" text-anchor="middle" class="axis">Exposed-position RMSD (A)</text>
<rect x="48" y="510" width="410" height="78" rx="8" class="callout"/>
<text x="64" y="536" class="callout-title">Technical-stability interval (not a p-value)</text>
<text x="64" y="564" class="callout-text">{_svg_number(bootstrap['delta_percentile_2_5_A'])} to {_svg_number(bootstrap['delta_percentile_97_5_A'])} A</text>
<rect x="474" y="510" width="438" height="78" rx="8" class="callout"/>
<text x="490" y="536" class="callout-title">Interpretation boundary</text>
<text x="490" y="564" class="callout-text">Technical job dependence; no biological replication or p-value.</text>
</svg>
'''


def _write_csv(path, rows):
    if not rows:
        raise ValueError(f"Cannot write empty required table {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_audit_outputs(
    output_dir,
    manifest_rows,
    tables,
    rank1_svg,
    rank2_svg,
    findings,
    checksum_rows,
):
    """Write exactly the fixed 11-file audit, refusing any existing output."""
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing audit folder: {output_dir}")
    expected_tables = {
        "per_control_geometry_summary",
        "job_pair_stability",
        "control_rank_and_leave_one_out",
        "technical_bootstrap_replicates",
        "technical_bootstrap_summary",
        "pose_cluster_membership",
    }
    if set(tables) != expected_tables:
        raise ValueError("Audit tables do not match the fixed six generated analysis tables")
    output_dir.mkdir(parents=True)
    _write_csv(output_dir / "model_job_identity_manifest.csv", manifest_rows)
    for table_name in sorted(expected_tables):
        _write_csv(output_dir / f"{table_name}.csv", tables[table_name])
    (output_dir / "rank1_primary_control_robustness.svg").write_text(
        rank1_svg, encoding="utf-8"
    )
    (output_dir / "rank2_length_sensitivity_job_dependence.svg").write_text(
        rank2_svg, encoding="utf-8"
    )
    (output_dir / "LEAD_FOCUSED_FINDINGS.md").write_text(findings, encoding="utf-8")
    _write_csv(output_dir / "frozen_input_checksums.csv", checksum_rows)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the fixed lead-focused EBV-MS structural robustness audit."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="EBV-MS publication project root",
    )
    args = parser.parse_args(argv)
    output = generate_live_audit(args.project_root)
    print(output)


if __name__ == "__main__":
    main()
