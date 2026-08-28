"""Freeze and expand score-blind HLA-DRB1*15:01 structural controls.

The selection protocol uses source provenance, peptide length, stable IEDB IDs,
and predeclared binding-rank bins.  It never uses discovery ranks, structural
similarity, RMSD, or AlphaFold confidence to choose controls.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.request
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Optional

from analyze_af3_pmhc_downloads import analyze_complete_job, parse_mmcif
from build_expanded_background_inputs import hla_chains_from_existing_batch
from build_premeeting_rigor_artifacts import predicted_core_positions
from premeeting_rigor import binding_rank_bin
from same_register_af3_analysis import same_register_geometry


ALLELE = "HLA-DRB1*15:01"
HUMAN = "Homo sapiens (human)"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
MYELIN_ACCESSIONS = {"P02686", "P60201", "Q16653"}
MYELIN_NAME = re.compile(
    r"\b(?:myelin basic protein|MBP protein|myelin proteolipid protein|"
    r"proteolipid protein|myelin oligodendrocyte glycoprotein)\b",
    re.IGNORECASE,
)
TARGET_STRATA = (21, 23, 25, 32)
ALPHAFOLD_SEED = 104773
ALPHAFOLD_BATCH_FILE = "alphafold_server_expanded_controls_jobs.json"
LAYER_PREFIXES = {
    "primary_exact_bin_length_pm1": "expanded_primary",
    "binding_bin_sensitivity_length_pm1": "expanded_binding_sensitivity",
    "length_sensitivity_exact_bin_pm7": "expanded_length_sensitivity",
    "length_plus_binding_sensitivity_pm7": "expanded_length_plus_binding_sensitivity",
}
ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "structural_control_expansion_2026-08-15"
RAW_EXPORT = ROOT / "raw" / "tcell_human_drb1501.json"
STUDY_MANIFEST = PROC / "pmhc_candidate_manifest.csv"
BASELINE_SCORE_SHEET = PROC / "complete_model_pipeline_audit_2026-08-15" / "master_pair_score_sheet.csv"
BASELINE_SAMPLES = PROC / "complete_model_pipeline_audit_2026-08-15" / "canonical_af3_sample_metrics.csv"
PAIR_UNIVERSE = PROC / "register_aware_benchmark" / "benchmark_pair_universe.csv"
EXISTING_PMHCS = PROC / "pmhc_colabfold_batch.csv"
DEFAULT_DOWNLOADS = ROOT / "folds_structural_controls_2026_08_15"
IEDB_NEXTGEN_PIPELINE = "https://api-nextgen-tools.iedb.org/api/v1/pipeline"
IEDB_NEXTGEN_METHOD = "netmhciipan_ba"


def _integer(value: Any, default: int = 10**18) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _source(record: dict[str, Any]) -> dict[str, Any]:
    value = record.get("curated_source_antigen")
    return value if isinstance(value, dict) else {}


def _representative(records: Iterable[dict[str, Any]]) -> dict[str, Any]:
    return min(
        records,
        key=lambda row: (_integer(row.get("structure_id")), _integer(row.get("tcell_id"))),
    )


def freeze_control_registry(
    raw_rows: list[dict[str, Any]], study_peptides: set[str]
) -> list[dict[str, object]]:
    """Return one auditable row per exact peptide in the saved IEDB export."""
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in raw_rows:
        grouped[str(record.get("linear_sequence") or "").upper()].append(record)

    registry: list[dict[str, object]] = []
    for peptide, records in grouped.items():
        representative = _representative(records)
        source = _source(representative)
        epitope_id = _integer(representative.get("structure_id"))
        accession = str(source.get("accession") or "")
        accession_root = accession.split(".", 1)[0]
        antigen = str(source.get("name") or representative.get("parent_source_antigen_name") or "")
        start = source.get("starting_position")
        end = source.get("ending_position")
        valid_coordinates = (
            accession
            and isinstance(start, int)
            and isinstance(end, int)
            and start >= 1
            and end >= start
            and end - start + 1 == len(peptide)
        )
        if representative.get("mhc_restriction") != ALLELE:
            status = "excluded_wrong_allele"
            reason = f"Record is not restricted to {ALLELE}."
        elif source.get("source_organism_name") != HUMAN:
            status = "excluded_nonhuman_source"
            reason = "Curated source antigen is not Homo sapiens."
        elif not peptide or not set(peptide).issubset(AMINO_ACIDS):
            status = "excluded_invalid_peptide"
            reason = "Exact canonical amino-acid peptide sequence is unavailable."
        elif not valid_coordinates:
            status = "excluded_invalid_coordinates"
            reason = "Accession and one-based source coordinates must span the exact peptide."
        elif peptide in study_peptides:
            status = "excluded_study_candidate"
            reason = "Exact peptide already belongs to the EBV-MS study candidate set."
        elif accession_root in MYELIN_ACCESSIONS or MYELIN_NAME.search(antigen):
            status = "excluded_mbp_plp_mog"
            reason = "Original neutral-pool rule excludes MBP, PLP, and MOG peptides."
        elif not 20 <= len(peptide) <= 30:
            status = "excluded_outside_predeclared_submission_lengths"
            reason = "Peptide cannot serve the predeclared 21/23/25 primary or 32-mer sensitivity strata."
        else:
            status = "eligible_pre_prediction"
            reason = "Eligible by frozen provenance, length, and original MBP/PLP/MOG exclusion rule."
        registry.append({
            "candidate_id": f"HUMAN_BACKGROUND_{epitope_id}",
            "iedb_epitope_id": epitope_id,
            "peptide": peptide,
            "peptide_length": len(peptide),
            "source_accession": accession,
            "source_antigen_name": antigen,
            "source_start_1_based": start if isinstance(start, int) else "",
            "source_end_1_based": end if isinstance(end, int) else "",
            "source_record_count": len(records),
            "source_assay_ids": ";".join(
                str(row.get("tcell_id"))
                for row in sorted(records, key=lambda row: _integer(row.get("tcell_id")))
                if row.get("tcell_id") is not None
            ),
            "mhc_allele": str(representative.get("mhc_restriction") or ""),
            "source_organism": str(source.get("source_organism_name") or ""),
            "selection_status": status,
            "selection_reason": reason,
        })
    return sorted(registry, key=lambda row: (_integer(row["iedb_epitope_id"]), str(row["peptide"])))


def _control_order(row: dict[str, Any], target_length: int) -> tuple[int, int, int, str]:
    peptide_length = int(row["peptide_length"])
    return (
        abs(peptide_length - target_length),
        peptide_length,
        _integer(row.get("iedb_epitope_id")),
        str(row["candidate_id"]),
    )


def _selection_row(
    prediction: dict[str, Any],
    *,
    stratum_length: int,
    analysis_layer: str,
    selection_order: int,
    served_pair_ids: list[str],
) -> dict[str, object]:
    """Project only frozen selection fields; outcome/structure fields cannot leak."""
    return {
        "candidate_id": str(prediction["candidate_id"]),
        "iedb_epitope_id": str(prediction["iedb_epitope_id"]),
        "peptide": str(prediction["peptide"]),
        "peptide_length": int(prediction["peptide_length"]),
        "stratum_length": stratum_length,
        "analysis_layer": analysis_layer,
        "binding_bin": str(prediction["binding_rank_bin"]),
        "predicted_percentile_rank": str(prediction["predicted_percentile_rank"]),
        "selection_order": selection_order,
        "served_pair_ids": ";".join(sorted(served_pair_ids)),
    }


def select_layered_controls(
    prediction_rows: list[dict[str, Any]],
    served_pairs_by_stratum: dict[int, list[str]],
    *,
    limit: int = 5,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Select predeclared primary and sensitivity controls without score leakage."""
    if limit < 1:
        raise ValueError("control limit must be positive")
    allowed = [
        row for row in prediction_rows
        if str(row.get("binding_rank_bin")) in {"weak", "intermediate"}
        and 20 <= int(row["peptide_length"]) <= 30
    ]
    selected: list[dict[str, object]] = []
    feasibility: list[dict[str, object]] = []
    for target_length in TARGET_STRATA:
        served = served_pairs_by_stratum.get(target_length, [])
        if target_length == 32:
            eligible = [row for row in allowed if 25 <= int(row["peptide_length"]) <= 30]
            weak = sorted(
                (row for row in eligible if row["binding_rank_bin"] == "weak"),
                key=lambda row: _control_order(row, target_length),
            )[:limit]
            intermediate = sorted(
                (row for row in eligible if row["binding_rank_bin"] == "intermediate"),
                key=lambda row: _control_order(row, target_length),
            )[: max(0, limit - len(weak))]
            for order, row in enumerate(weak, start=1):
                selected.append(_selection_row(
                    row,
                    stratum_length=target_length,
                    analysis_layer="length_sensitivity_exact_bin_pm7",
                    selection_order=order,
                    served_pair_ids=served,
                ))
            for offset, row in enumerate(intermediate, start=len(weak) + 1):
                selected.append(_selection_row(
                    row,
                    stratum_length=target_length,
                    analysis_layer="length_plus_binding_sensitivity_pm7",
                    selection_order=offset,
                    served_pair_ids=served,
                ))
            chosen_count = len(weak) + len(intermediate)
            feasibility.append({
                "stratum_length": target_length,
                "served_pair_ids": ";".join(sorted(served)),
                "primary_assessment": "not_assessable_no_31_to_33aa_direct_controls",
                "primary_exact_bin_count": 0,
                "binding_sensitivity_count": 0,
                "length_sensitivity_exact_bin_count": len(weak),
                "length_plus_binding_sensitivity_count": len(intermediate),
                "selected_control_mapping_count": chosen_count,
                "control_shortfall": max(0, limit - chosen_count),
                "feasibility_status": "target_met" if chosen_count == limit else "partial_no_further_relaxation",
            })
            continue

        eligible = [
            row for row in allowed
            if abs(int(row["peptide_length"]) - target_length) <= 1
        ]
        weak = sorted(
            (row for row in eligible if row["binding_rank_bin"] == "weak"),
            key=lambda row: _control_order(row, target_length),
        )[:limit]
        intermediate = sorted(
            (row for row in eligible if row["binding_rank_bin"] == "intermediate"),
            key=lambda row: _control_order(row, target_length),
        )[: max(0, limit - len(weak))]
        for order, row in enumerate(weak, start=1):
            selected.append(_selection_row(
                row,
                stratum_length=target_length,
                analysis_layer="primary_exact_bin_length_pm1",
                selection_order=order,
                served_pair_ids=served,
            ))
        for offset, row in enumerate(intermediate, start=len(weak) + 1):
            selected.append(_selection_row(
                row,
                stratum_length=target_length,
                analysis_layer="binding_bin_sensitivity_length_pm1",
                selection_order=offset,
                served_pair_ids=served,
            ))
        chosen_count = len(weak) + len(intermediate)
        feasibility.append({
            "stratum_length": target_length,
            "served_pair_ids": ";".join(sorted(served)),
            "primary_assessment": "assessable_exact_length_and_binding_rule" if weak else "not_assessable_no_exact_bin_control",
            "primary_exact_bin_count": len(weak),
            "binding_sensitivity_count": len(intermediate),
            "length_sensitivity_exact_bin_count": 0,
            "length_plus_binding_sensitivity_count": 0,
            "selected_control_mapping_count": chosen_count,
            "control_shortfall": max(0, limit - chosen_count),
            "feasibility_status": "target_met" if chosen_count == limit else "partial_no_further_relaxation",
        })
    return selected, feasibility


def build_alphafold_jobs(
    selected_rows: list[dict[str, Any]], dra_sequence: str, drb_sequence: str
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Build one AF Server job per unique control and retain all layer mappings."""
    if not dra_sequence or not drb_sequence:
        raise ValueError("Both frozen HLA chain sequences are required")
    candidates: dict[str, tuple[str, int]] = {}
    for row in selected_rows:
        candidate_id = str(row["candidate_id"])
        value = (str(row["peptide"]), int(row["peptide_length"]))
        if candidate_id in candidates and candidates[candidate_id] != value:
            raise ValueError(f"Conflicting peptide mapping for {candidate_id}")
        candidates[candidate_id] = value
    if len(candidates) > 30:
        raise ValueError("Expanded control batch exceeds the predeclared 30-job limit")

    jobs: list[dict[str, object]] = []
    job_name_by_candidate: dict[str, str] = {}
    for candidate_id, (peptide, _) in sorted(candidates.items()):
        job_name = f"ebvms_bg_{candidate_id}_s04"
        job_name_by_candidate[candidate_id] = job_name
        jobs.append({
            "name": job_name,
            "modelSeeds": [ALPHAFOLD_SEED],
            "sequences": [
                {"proteinChain": {"sequence": chain, "count": 1}}
                for chain in (dra_sequence, drb_sequence, peptide)
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        })
    manifest = [
        {
            "batch_file": ALPHAFOLD_BATCH_FILE,
            "candidate_id": row["candidate_id"],
            "job_name": job_name_by_candidate[str(row["candidate_id"])],
            "chain_order": "HLA-DRA;HLA-DRB1*15:01;human-background peptide",
            "peptide": row["peptide"],
            "peptide_length": row["peptide_length"],
            "stratum_length": row["stratum_length"],
            "analysis_layer": row["analysis_layer"],
            "binding_bin": row["binding_bin"],
            "selection_order": row["selection_order"],
            "served_pair_ids": row["served_pair_ids"],
        }
        for row in selected_rows
    ]
    return jobs, manifest


def summarize_layer_geometry(
    geometry_rows: list[dict[str, Any]], target_median: Optional[float]
) -> dict[str, object]:
    """Equal-weight unique controls after collapsing their technical comparisons."""
    values_by_control: dict[str, list[float]] = defaultdict(list)
    for row in geometry_rows:
        values_by_control[str(row["background_candidate_id"])].append(
            float(row["candidate_exposed_ca_rmsd_A"])
        )
    if not values_by_control:
        return {
            "unique_control_count": 0,
            "technical_geometry_count": 0,
            "background_control_median_A": "",
            "background_minus_target_median_A": "",
            "p_value": "",
        }
    control_medians = [
        median(values) for _, values in sorted(values_by_control.items())
    ]
    background_median = median(control_medians)
    return {
        "unique_control_count": len(values_by_control),
        "technical_geometry_count": len(geometry_rows),
        "background_control_median_A": round(background_median, 6),
        "background_minus_target_median_A": (
            round(background_median - target_median, 6)
            if target_median is not None else ""
        ),
        "p_value": "",
    }


def extend_score_sheet_with_layered_controls(
    baseline_rows: list[dict[str, Any]],
    selected_rows: list[dict[str, Any]],
    geometry_rows: list[dict[str, Any]],
) -> list[dict[str, object]]:
    """Append separate control-layer fields without changing discovery ordering."""
    selected_by_pair_layer: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in selected_rows:
        for pair_id in str(row.get("served_pair_ids", "")).split(";"):
            if pair_id:
                selected_by_pair_layer[(pair_id, str(row["analysis_layer"]))].add(
                    str(row["candidate_id"])
                )
    geometry_by_pair_layer: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in geometry_rows:
        geometry_by_pair_layer[(str(row["pair_id"]), str(row["analysis_layer"]))].append(row)

    output: list[dict[str, object]] = []
    for baseline in baseline_rows:
        pair_id = str(baseline["pair_id"])
        row: dict[str, object] = dict(baseline)
        any_selected = False
        total_completed = 0
        target_text = str(
            baseline.get("target_candidate_exposed_rmsd_median_A", "")
            or baseline.get("candidate_exposed_ca_rmsd_A_median", "")
        )
        target_median = float(target_text) if target_text else None
        for layer, prefix in LAYER_PREFIXES.items():
            selected_ids = sorted(selected_by_pair_layer[(pair_id, layer)])
            any_selected = any_selected or bool(selected_ids)
            summary = summarize_layer_geometry(
                geometry_by_pair_layer[(pair_id, layer)], target_median
            )
            total_completed += int(summary["unique_control_count"])
            row.update({
                f"{prefix}_selected_control_count": len(selected_ids),
                f"{prefix}_selected_control_ids": ";".join(selected_ids),
                f"{prefix}_completed_control_count": summary["unique_control_count"],
                f"{prefix}_technical_geometry_count": summary["technical_geometry_count"],
                f"{prefix}_background_control_median_A": summary["background_control_median_A"],
                f"{prefix}_background_minus_target_median_A": summary["background_minus_target_median_A"],
                f"{prefix}_p_value": "",
            })
        if not any_selected:
            status = "not_selected_for_expansion"
        elif total_completed == 0:
            status = "awaiting_alphafold_downloads"
        else:
            status = "descriptive_layered_geometry_available"
        row["expanded_control_status"] = status
        row["expanded_control_claim_boundary"] = (
            "Descriptive pMHC structural context only; not evidence of TCR binding, "
            "cross-reactivity, activation, molecular mimicry, or MS mechanism."
        )
        output.append(row)
    return output


def uncovered_pairs_by_stratum(
    score_rows: list[dict[str, Any]],
) -> dict[int, list[str]]:
    """Freeze only uncovered, primary-allele-eligible pairs into target strata."""
    result = {length: [] for length in TARGET_STRATA}
    for row in score_rows:
        if str(row.get("register_eligible_primary_allele")) != "True":
            continue
        if _integer(row.get("structural_background_comparator_count"), 0) != 0:
            continue
        length = len(str(row.get("human_peptide", "")))
        if length in result:
            result[length].append(str(row["pair_id"]))
    return {length: sorted(pair_ids) for length, pair_ids in result.items()}


def _request_job(path: Path) -> dict[str, Any]:
    request = json.loads(path.read_text(encoding="utf-8"))
    return request[0] if isinstance(request, list) else request


def inventory_expected_control_jobs(
    downloads_root: Path, manifest_rows: list[dict[str, Any]]
) -> list[dict[str, object]]:
    """Inventory each unique expected job and validate its submitted peptide."""
    mappings: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in manifest_rows:
        mappings[str(row["job_name"])].append(row)
    directories = {
        path.name.lower(): path for path in downloads_root.iterdir() if path.is_dir()
    } if downloads_root.exists() else {}
    output: list[dict[str, object]] = []
    for job_name, rows in sorted(mappings.items()):
        candidate_ids = {str(row["candidate_id"]) for row in rows}
        peptides = {str(row["peptide"]) for row in rows}
        if len(candidate_ids) != 1 or len(peptides) != 1:
            raise ValueError(f"Conflicting manifest mappings for {job_name}")
        candidate_id = next(iter(candidate_ids))
        peptide = next(iter(peptides))
        directory = directories.get(job_name.lower())
        if directory is None:
            cifs = summaries = full_data = requests = 0
            status = "not_downloaded"
            request_name = observed_peptide = ""
        else:
            cifs = len(list(directory.glob("*_model_*.cif")))
            summaries = len(list(directory.glob("*_summary_confidences_*.json")))
            full_data = len(list(directory.glob("*_full_data_*.json")))
            request_paths = list(directory.glob("*_job_request.json"))
            requests = len(request_paths)
            request_name = observed_peptide = ""
            exact_request = False
            if requests == 1:
                job = _request_job(request_paths[0])
                chains = [
                    entry["proteinChain"]["sequence"]
                    for entry in job.get("sequences", [])
                    if "proteinChain" in entry
                ]
                request_name = str(job.get("name", ""))
                observed_peptide = chains[2] if len(chains) == 3 else ""
                exact_request = request_name == job_name and len(chains) == 3 and observed_peptide == peptide
            if (cifs, summaries, full_data, requests) == (5, 5, 5, 1) and exact_request:
                status = "complete_five_sample_exact_sequence"
            elif (cifs, summaries, full_data, requests) == (5, 5, 5, 1):
                status = "excluded_request_or_peptide_mismatch"
            else:
                status = "partial_download_excluded"
        output.append({
            "batch_folder": downloads_root.name,
            "batch_file": str(rows[0].get("batch_file", ALPHAFOLD_BATCH_FILE)),
            "candidate_id": candidate_id,
            "expected_job_name": job_name,
            "downloaded_job_directory": directory.name if directory else "",
            "expected_peptide": peptide,
            "request_name": request_name,
            "observed_peptide": observed_peptide,
            "served_stratum_lengths": ";".join(
                str(value) for value in sorted({int(row["stratum_length"]) for row in rows})
            ),
            "model_cif_count": cifs,
            "summary_confidences_count": summaries,
            "full_data_count": full_data,
            "job_request_count": requests,
            "completeness_status": status,
            "claim_boundary": "Technical availability and exact-sequence QA only; not biological evidence.",
        })
    return output


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: Optional[list[str]] = None,
) -> None:
    if fieldnames is None:
        if not rows:
            raise ValueError(f"Field names are required for empty table {path}")
        fieldnames = list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def build_prediction_inputs(
    registry: list[dict[str, object]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Create one direct IEDB submission per frozen eligible peptide."""
    candidates: list[dict[str, str]] = []
    submissions: list[dict[str, str]] = []
    for row in registry:
        if row["selection_status"] != "eligible_pre_prediction":
            continue
        candidate = {
            "candidate_id": str(row["candidate_id"]),
            "arm": "Expanded human background comparator",
            "evidence_tier": "Tier 4: source-coordinate-validated HLA-DRB1*15:01 record",
            "peptide": str(row["peptide"]),
            "peptide_length": str(row["peptide_length"]),
            "hla": ALLELE,
            "iedb_epitope_id": str(row["iedb_epitope_id"]),
            "source_accession": str(row["source_accession"]),
            "source_antigen_name": str(row["source_antigen_name"]),
            "source_start_1_based": str(row["source_start_1_based"]),
            "source_end_1_based": str(row["source_end_1_based"]),
        }
        candidates.append(candidate)
        submissions.append({
            "submission_id": f"{candidate['candidate_id']}__segment_001",
            "candidate_id": candidate["candidate_id"],
            "peptide": candidate["peptide"],
            "source_start_1_based": "1",
            "submission_strategy": "direct_full_peptide",
            "claim_boundary": "Computational binding/register hypothesis only; not experimental presentation evidence.",
        })
    return candidates, submissions


def attach_prediction_provenance(
    candidates: list[dict[str, str]], prediction_rows: list[dict[str, object]]
) -> list[dict[str, object]]:
    source = {row["candidate_id"]: row for row in candidates}
    output: list[dict[str, object]] = []
    for prediction in prediction_rows:
        candidate = source[str(prediction["candidate_id"])]
        output.append({
            **prediction,
            "iedb_epitope_id": candidate["iedb_epitope_id"],
            "source_accession": candidate["source_accession"],
            "source_antigen_name": candidate["source_antigen_name"],
            "source_start_1_based": candidate["source_start_1_based"],
            "source_end_1_based": candidate["source_end_1_based"],
        })
    return output


def prediction_rows_from_nextgen_result(
    candidates: list[dict[str, str]],
    payload: dict[str, Any],
    raw_response_path: str,
    retrieved_utc: str,
) -> list[dict[str, object]]:
    """Normalize the official NG IEDB peptide table in frozen input order."""
    result = payload.get("result", payload)
    if result.get("status") != "done":
        raise ValueError(f"IEDB Next-Generation result is not done: {result.get('status')}")
    data = result.get("data") or {}
    if data.get("errors"):
        raise ValueError(f"IEDB Next-Generation errors: {data['errors']}")
    tables = [table for table in data.get("results", []) if table.get("type") == "peptide_table"]
    if len(tables) != 1:
        raise ValueError("Expected exactly one IEDB Next-Generation peptide table")
    table = tables[0]
    column_names = [column["name"] for column in table.get("table_columns", [])]
    required = {
        "sequence_number", "peptide", "length", "allele",
        "netmhciipan_ba_core", "netmhciipan_ba_ic50",
        "netmhciipan_ba_percentile",
    }
    if not required.issubset(column_names):
        raise ValueError("IEDB Next-Generation peptide table lacks required binding fields")
    by_sequence_number: dict[int, dict[str, Any]] = {}
    for values in table.get("table_data", []):
        row = dict(zip(column_names, values))
        sequence_number = int(row["sequence_number"])
        if sequence_number in by_sequence_number:
            raise ValueError(f"Duplicate IEDB sequence number {sequence_number}")
        by_sequence_number[sequence_number] = row
    if set(by_sequence_number) != set(range(1, len(candidates) + 1)):
        raise ValueError("IEDB Next-Generation response does not cover every frozen input exactly once")

    output: list[dict[str, object]] = []
    for sequence_number, candidate in enumerate(candidates, start=1):
        prediction = by_sequence_number[sequence_number]
        peptide = candidate["peptide"]
        if prediction["peptide"] != peptide or int(prediction["length"]) != len(peptide):
            raise ValueError(f"IEDB peptide mismatch for {candidate['candidate_id']}")
        if prediction["allele"] != ALLELE:
            raise ValueError(f"IEDB allele mismatch for {candidate['candidate_id']}")
        core = str(prediction["netmhciipan_ba_core"])
        starts = predicted_core_positions(peptide, core)
        if not starts:
            raise ValueError(f"IEDB core is not contained in {candidate['candidate_id']}")
        rank = float(prediction["netmhciipan_ba_percentile"])
        output.append({
            "candidate_id": candidate["candidate_id"],
            "arm": candidate["arm"],
            "evidence_tier": candidate["evidence_tier"],
            "peptide": peptide,
            "peptide_length": candidate["peptide_length"],
            "hla": candidate["hla"],
            "prediction_method_requested": "NetMHCIIpan 4.1 BA (IEDB recommended binding predictor-2023.09)",
            "prediction_endpoint": IEDB_NEXTGEN_PIPELINE,
            "prediction_retrieval_utc": retrieved_utc,
            "raw_response_path": raw_response_path,
            "interpretation": "Computational binding/register hypothesis only; not experimental presentation evidence.",
            "prediction_status": "predicted",
            "submission_strategy": "direct_full_peptide_as_is",
            "submission_segment_count": 1,
            "iedb_seq_num": sequence_number,
            "predicted_core_peptide": core,
            "predicted_core_start_positions_1_based": starts,
            "prediction_input_peptide": peptide,
            "predicted_core_fully_contained_in_manifest_peptide": True,
            "predicted_ic50_nM": prediction["netmhciipan_ba_ic50"],
            "predicted_percentile_rank": prediction["netmhciipan_ba_percentile"],
            "binding_rank_bin": binding_rank_bin(rank),
            "iedb_epitope_id": candidate["iedb_epitope_id"],
            "source_accession": candidate["source_accession"],
            "source_antigen_name": candidate["source_antigen_name"],
            "source_start_1_based": candidate["source_start_1_based"],
            "source_end_1_based": candidate["source_end_1_based"],
        })
    return output


def fetch_nextgen_iedb_predictions(
    candidates: list[dict[str, str]],
    *,
    timeout_seconds: int = 600,
    poll_interval_seconds: float = 2.0,
) -> dict[str, Any]:
    """Submit and poll the official asynchronous NG IEDB MHC-II pipeline."""
    fasta = "\n".join(
        f">{candidate['candidate_id']}\n{candidate['peptide']}"
        for candidate in candidates
    )
    payload = {
        "pipeline_id": "",
        "pipeline_title": "EBV-MS frozen expanded structural controls",
        "run_stage_range": [1, 1],
        "stages": [{
            "stage_number": 1,
            "stage_type": "prediction",
            "tool_group": "mhcii",
            "input_sequence_text": fasta,
            "input_parameters": {
                "alleles": ALLELE,
                "peptide_length_range": None,
                "predictors": [{"type": "binding", "method": IEDB_NEXTGEN_METHOD}],
            },
        }],
    }
    request = urllib.request.Request(
        IEDB_NEXTGEN_PIPELINE,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        submission = json.load(response)
    if submission.get("errors"):
        raise ValueError(f"IEDB Next-Generation submission errors: {submission['errors']}")
    results_uri = submission.get("results_uri")
    if not results_uri:
        raise ValueError("IEDB Next-Generation submission omitted results_uri")
    deadline = time.monotonic() + timeout_seconds
    while True:
        with urllib.request.urlopen(str(results_uri), timeout=30) as response:
            result = json.load(response)
        if result.get("status") == "done":
            return {"submission": submission, "result": result}
        if result.get("status") in {"error", "failed"}:
            raise ValueError(f"IEDB Next-Generation job failed: {result}")
        if time.monotonic() >= deadline:
            raise TimeoutError("IEDB Next-Generation prediction job did not complete before timeout")
        time.sleep(poll_interval_seconds)


def extract_downloaded_control_geometry(
    downloads_root: Path,
    manifest_rows: list[dict[str, Any]],
    inventory_rows: list[dict[str, Any]],
    prediction_rows: list[dict[str, Any]],
    pair_rows: list[dict[str, str]],
    baseline_sample_rows: list[dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Validate completed jobs and compare them with existing EBV pMHC samples."""
    complete = {
        str(row["candidate_id"]): str(row["downloaded_job_directory"])
        for row in inventory_rows
        if row["completeness_status"] == "complete_five_sample_exact_sequence"
    }
    if not complete:
        return [], [], [], []
    metadata = {
        candidate_id: {"af3_cohort": "expanded_structural_control"}
        for candidate_id in complete
    }
    control_samples: dict[str, list[dict[str, object]]] = defaultdict(list)
    job_summaries: list[dict[str, object]] = []
    sample_summaries: list[dict[str, object]] = []
    for candidate_id, directory_name in sorted(complete.items()):
        directory = downloads_root / directory_name
        result = analyze_complete_job(directory, downloads_root.name, metadata)
        job_summaries.append(dict(result["job"]))
        for sample in result["samples"]:
            sample = dict(sample)
            sample_index = int(sample["sample_index"])
            matches = list(directory.glob(f"*_model_{sample_index}.cif"))
            if len(matches) != 1:
                raise FileNotFoundError(f"Expected one model {sample_index} in {directory}")
            sample["model_path"] = str(matches[0])
            sample["canonical_job_key"] = f"{candidate_id}|s04|{sample['requested_peptide']}"
            sample_summaries.append(sample)
            if sample["sequence_layout_status"] == "pass_exact_three_chain_peptide_match":
                control_samples[candidate_id].append(sample)

    predictions = {str(row["candidate_id"]): row for row in prediction_rows}
    pairs = {row["pair_id"]: row for row in pair_rows}
    ebv_samples: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in baseline_sample_rows:
        if (
            row.get("af3_cohort") == "legacy_candidate_pmhc"
            and row.get("sequence_layout_status") == "pass_exact_three_chain_peptide_match"
        ):
            ebv_samples[row["candidate_id"]].append(row)
    model_cache: dict[str, dict[str, list[dict[str, object]]]] = {}
    geometry: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    seen_mapping: set[tuple[str, str, str]] = set()
    for mapping in manifest_rows:
        candidate_id = str(mapping["candidate_id"])
        prediction = predictions[candidate_id]
        starts = str(prediction["predicted_core_start_positions_1_based"]).split(";")
        pair_ids = [value for value in str(mapping["served_pair_ids"]).split(";") if value]
        for pair_id in pair_ids:
            mapping_key = (pair_id, candidate_id, str(mapping["analysis_layer"]))
            if mapping_key in seen_mapping:
                continue
            seen_mapping.add(mapping_key)
            if candidate_id not in control_samples:
                continue
            if len(starts) != 1 or not starts[0].isdigit():
                exclusions.append({
                    "pair_id": pair_id,
                    "background_candidate_id": candidate_id,
                    "analysis_layer": mapping["analysis_layer"],
                    "exclusion_reason": "Predicted P1-P9 core is not uniquely positioned in the full peptide.",
                })
                continue
            pair = pairs[pair_id]
            for ebv_sample in ebv_samples[str(pair["ebv_candidate_id"])]:
                ebv_path = str(ebv_sample["model_path"])
                if ebv_path not in model_cache:
                    model_cache[ebv_path] = parse_mmcif(Path(ebv_path))
                for control_sample in control_samples[candidate_id]:
                    control_path = str(control_sample["model_path"])
                    if control_path not in model_cache:
                        model_cache[control_path] = parse_mmcif(Path(control_path))
                    metrics = same_register_geometry(
                        model_cache[ebv_path],
                        model_cache[control_path],
                        int(pair["ebv_top_core_start_1_based"]),
                        int(starts[0]),
                    )
                    geometry.append({
                        "pair_id": pair_id,
                        "ebv_candidate_id": pair["ebv_candidate_id"],
                        "target_human_candidate_id": pair["human_candidate_id"],
                        "background_candidate_id": candidate_id,
                        "stratum_length": mapping["stratum_length"],
                        "analysis_layer": mapping["analysis_layer"],
                        "binding_bin": mapping["binding_bin"],
                        "selection_order": mapping["selection_order"],
                        "background_predicted_core": prediction["predicted_core_peptide"],
                        "ebv_job_key": ebv_sample["canonical_job_key"],
                        "ebv_sample_index": ebv_sample["sample_index"],
                        "background_job_key": control_sample["canonical_job_key"],
                        "background_sample_index": control_sample["sample_index"],
                        **{key: round(value, 6) for key, value in metrics.items()},
                        "interpretation": "Frozen score-blind pMHC comparison; technical samples are not biological replicates.",
                    })
    return geometry, exclusions, job_summaries, sample_summaries


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_protocol(source_hash: str) -> None:
    protocol = {
        "protocol_frozen_before_structure_review": True,
        "source_export": str(RAW_EXPORT),
        "source_export_sha256": source_hash,
        "source_row_count_expected": 1139,
        "allele": ALLELE,
        "binding_predictor": "NetMHCIIpan 4.1 BA (IEDB recommended binding predictor-2023.09)",
        "binding_prediction_endpoint": IEDB_NEXTGEN_PIPELINE,
        "target_strata_aa": list(TARGET_STRATA),
        "controls_per_stratum_target": 5,
        "original_exclusion_rule": "exclude exact study peptides plus MBP, PLP, and MOG; other CNS proteins remain eligible",
        "primary_rule": "target length +/-1 aa and weak predicted binding bin",
        "binding_sensitivity_rule": "intermediate bin within target length +/-1; separate from primary",
        "long_peptide_rule": "32-aa primary remains not assessable; 25-30 aa controls form separate +/-7 length sensitivity",
        "tie_break_order": ["absolute_length_distance", "peptide_length", "numeric_iedb_epitope_id"],
        "post_hoc_relaxation": "forbidden",
        "p_values": "not calculated",
        "claim_boundary": "Computational pMHC structural context only; no TCR or mechanistic inference.",
    }
    with (OUT / "frozen_protocol.json").open("w", encoding="utf-8") as handle:
        json.dump(protocol, handle, indent=2)
        handle.write("\n")


def _write_findings(
    registry: list[dict[str, object]],
    predictions: list[dict[str, object]],
    selected: list[dict[str, object]],
    feasibility: list[dict[str, object]],
    jobs: list[dict[str, object]],
    inventory: list[dict[str, object]],
    geometry: list[dict[str, object]],
    expanded_scores: list[dict[str, object]],
) -> None:
    status_counts: dict[str, int] = defaultdict(int)
    for row in registry:
        status_counts[str(row["selection_status"])] += 1
    lines = [
        "# Expanded structural-control audit",
        "",
        "## Current state",
        "",
        f"- Saved IEDB assay rows scanned: **{sum(int(row['source_record_count']) for row in registry)}**",
        f"- Deduplicated exact peptides in the frozen registry: **{len(registry)}**",
        f"- Peptides predicted for HLA-DRB1*15:01: **{len(predictions)}**",
        f"- Selected stratum/layer mappings: **{len(selected)}**",
        f"- Unique AlphaFold Server jobs: **{len(jobs)}**",
        f"- Complete exact-sequence AlphaFold jobs currently available: **{sum(row['completeness_status'] == 'complete_five_sample_exact_sequence' for row in inventory)}**",
        f"- Layered structural geometry rows currently available: **{len(geometry)}**",
        "",
        "The discovery ranking is copied unchanged from the 2026-08-15 audit. Control selection did not inspect ranking, similarity, RMSD, AlphaFold confidence, or geometry fields.",
        "",
        "## Frozen strata",
        "",
    ]
    for row in feasibility:
        lines.append(
            f"- {row['stratum_length']} aa: {row['selected_control_mapping_count']} selected; "
            f"primary={row['primary_exact_bin_count']}, binding sensitivity={row['binding_sensitivity_count']}, "
            f"length sensitivity={row['length_sensitivity_exact_bin_count']}, "
            f"length+binding sensitivity={row['length_plus_binding_sensitivity_count']}, "
            f"shortfall={row['control_shortfall']} ({row['feasibility_status']})."
        )
    lines.extend([
        "",
        "The 32-aa primary stratum remains **not assessable** because the direct IEDB MHC-II range ends at 30 aa and no eligible 31-33 aa primary control exists. Its 25-30 aa controls remain a separate length-sensitivity layer.",
        "",
        "## Registry audit counts",
        "",
    ])
    lines.extend(f"- {status}: **{count}**" for status, count in sorted(status_counts.items()))
    completed_rows = [
        row for row in expanded_scores
        if row.get("expanded_control_status") == "descriptive_layered_geometry_available"
    ]
    lines.extend([
        "",
        "## Descriptive structural comparisons",
        "",
        "A positive background-minus-target delta means the candidate pair has a lower exposed-position RMSD than the median of its equal-weighted controls. A negative value means the controls have the lower median RMSD. These are descriptive comparisons only.",
        "",
        "| Discovery rank | Pair | Layer | Unique controls | Target median RMSD (A) | Background median RMSD (A) | Background - target (A) |",
        "|---:|---|---|---:|---:|---:|---:|",
    ])
    for row in completed_rows:
        if int(row["expanded_primary_completed_control_count"]):
            prefix = "expanded_primary"
            layer = "primary exact-bin, length +/-1"
        else:
            prefix = "expanded_length_sensitivity"
            layer = "32-aa length sensitivity only"
        lines.append(
            f"| {row['discovery_priority_rank']} | {row['pair_id']} | {layer} | "
            f"{row[f'{prefix}_completed_control_count']} | "
            f"{row['candidate_exposed_ca_rmsd_A_median']} | "
            f"{row[f'{prefix}_background_control_median_A']} | "
            f"{row[f'{prefix}_background_minus_target_median_A']} |"
        )
    lines.extend([
        "",
        "## Interpretation limit",
        "",
        "These controls provide descriptive pMHC structural context only. They do not establish peptide presentation, TCR binding, shared-TCR recognition, cross-reactivity, activation, molecular mimicry, or an MS mechanism. No p-value is reported.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "PYTHONPATH=src python3 src/build_structural_control_expansion.py --response-file processed/structural_control_expansion_2026-08-15/iedb_nextgen_mhcii_raw.json",
        "```",
        "",
    ])
    (OUT / "AUDIT_FINDINGS.md").write_text("\n".join(lines), encoding="utf-8")


GEOMETRY_FIELDS = [
    "pair_id", "ebv_candidate_id", "target_human_candidate_id",
    "background_candidate_id", "stratum_length", "analysis_layer", "binding_bin",
    "selection_order", "background_predicted_core", "ebv_job_key",
    "ebv_sample_index", "background_job_key", "background_sample_index",
    "hla_groove_ca_rmsd_A", "core_p1_p9_ca_rmsd_A", "anchor_ca_rmsd_A",
    "candidate_exposed_ca_rmsd_A", "interpretation",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--response-file", type=Path, help="Use a saved IEDB response instead of a live request.")
    parser.add_argument("--prepare-only", action="store_true", help="Freeze registry and IEDB submissions without selecting controls.")
    parser.add_argument("--downloads-root", type=Path, default=DEFAULT_DOWNLOADS)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    raw_rows = json.loads(RAW_EXPORT.read_text(encoding="utf-8"))
    if len(raw_rows) != 1139:
        raise ValueError(f"Frozen IEDB export row count changed: expected 1139, found {len(raw_rows)}")
    study_peptides = {row["peptide"].upper() for row in read_csv(STUDY_MANIFEST)}
    registry = freeze_control_registry(raw_rows, study_peptides)
    candidates, submissions = build_prediction_inputs(registry)
    write_csv(OUT / "frozen_control_universe.csv", registry)
    write_csv(OUT / "iedb_submission_manifest.csv", submissions)
    _write_protocol(_sha256(RAW_EXPORT))
    if args.prepare_only:
        print(f"Prepared {len(candidates)} frozen prediction inputs in {OUT}")
        return

    nextgen_payload = (
        json.loads(args.response_file.read_text(encoding="utf-8"))
        if args.response_file else fetch_nextgen_iedb_predictions(candidates)
    )
    raw_response = OUT / "iedb_nextgen_mhcii_raw.json"
    with raw_response.open("w", encoding="utf-8") as handle:
        json.dump(nextgen_payload, handle, indent=2)
        handle.write("\n")
    predictions = prediction_rows_from_nextgen_result(
        candidates,
        nextgen_payload,
        str(raw_response),
        time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    )
    write_csv(OUT / "control_binding_prediction_summary.csv", predictions)

    baseline_scores = read_csv(BASELINE_SCORE_SHEET)
    strata = uncovered_pairs_by_stratum(baseline_scores)
    selected, feasibility = select_layered_controls(predictions, strata, limit=5)
    dra, drb = hla_chains_from_existing_batch(read_csv(EXISTING_PMHCS))
    jobs, manifest = build_alphafold_jobs(selected, dra, drb)
    with (OUT / ALPHAFOLD_BATCH_FILE).open("w", encoding="utf-8") as handle:
        json.dump(jobs, handle, indent=2)
        handle.write("\n")
    write_csv(OUT / "selected_control_manifest.csv", manifest)
    write_csv(OUT / "control_selection_feasibility.csv", feasibility)

    inventory = inventory_expected_control_jobs(args.downloads_root, manifest)
    write_csv(OUT / "alphafold_download_inventory.csv", inventory)
    geometry, exclusions, job_summaries, sample_summaries = extract_downloaded_control_geometry(
        args.downloads_root,
        manifest,
        inventory,
        predictions,
        read_csv(PAIR_UNIVERSE),
        read_csv(BASELINE_SAMPLES),
    )
    write_csv(OUT / "complete_layered_control_geometry.csv", geometry, GEOMETRY_FIELDS)
    if exclusions:
        write_csv(OUT / "control_geometry_exclusions.csv", exclusions)
    if job_summaries:
        write_csv(OUT / "alphafold_control_job_summary.csv", job_summaries)
    if sample_summaries:
        write_csv(OUT / "alphafold_control_sample_metrics.csv", sample_summaries)
    expanded_scores = extend_score_sheet_with_layered_controls(
        baseline_scores, selected, geometry
    )
    write_csv(OUT / "master_pair_score_sheet_with_expanded_controls.csv", expanded_scores)
    _write_findings(
        registry, predictions, selected, feasibility, jobs, inventory, geometry,
        expanded_scores,
    )
    print(
        f"Prepared {len(jobs)} unique AlphaFold jobs across {len(selected)} layer mappings; "
        f"wrote {len(geometry)} available geometry rows to {OUT}"
    )


if __name__ == "__main__":
    main()
