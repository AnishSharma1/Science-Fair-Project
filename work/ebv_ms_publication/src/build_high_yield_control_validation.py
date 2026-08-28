"""Build the additive high-yield N3 ranking-context package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any, Mapping, Sequence

from high_yield_control_validation import (
    CLAIM_BOUNDARY,
    FROZEN_TARGETS,
    SEED,
    SURFACE_FEATURES,
    build_n3_panel,
    build_ranking_context_gate,
    build_specificity_gate,
    rank_panel_rows,
    select_comparator_arms,
    validate_frozen_targets,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V3 = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
DEFAULT_PANEL = ROOT / "processed/tcell_library_v2_2026-08-22/frozen_v2_80_peptide_panel.csv"
DEFAULT_OUT = ROOT / "processed/high_yield_control_validation_2026-08-28"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] = ()) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or sorted({key for row in rows for key in row}))
    if not fieldnames:
        raise ValueError(f"field names are required for empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="raise",
            lineterminator="\n",
        )
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


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def _normalized_identity(row: Mapping[str, Any]) -> tuple[str, str]:
    accession = re.sub(r"[^A-Z0-9]+", "", str(row.get("accession") or row.get("protein")).upper())
    core = re.sub(r"[^A-Z0-9]+", "", str(row.get("core")).upper())
    return accession, core


def _eligible_unique_count(catalog: Sequence[Mapping[str, Any]], excluded_ids: set[str]) -> int:
    return len({
        _normalized_identity(row)
        for row in catalog
        if str(row["candidate_id"]) not in excluded_ids
        and float(row["binding_percentile"]) <= 20
        and int(row["model_count"]) == 5
        and row["surface_status"] == "complete"
    })


def _panel_metadata(panel_path: Path) -> dict[str, dict[str, str]]:
    rows = read_csv(panel_path)
    mapping = {row["candidate_id"]: row for row in rows}
    if len(mapping) != 80:
        raise ValueError("the frozen V2 panel must contain 80 unique candidate IDs")
    return mapping


def _arm_catalog(
    rows: Sequence[Mapping[str, Any]],
    metadata: Mapping[str, Mapping[str, Any]],
    *,
    allele: str,
    side: str,
) -> list[dict[str, Any]]:
    prefix = "ebv" if side == "ebv" else "self"
    records: dict[str, dict[str, Any]] = {}
    for row in rows:
        if str(row["allele"]) != allele:
            continue
        candidate_id = str(row[f"{prefix}_candidate_id"])
        source = metadata[candidate_id]
        record = {
            "candidate_id": candidate_id,
            "allele": allele,
            "arm_class": side,
            "sequence": row[f"{prefix}_sequence"],
            "core": row[f"{prefix}_core_p1_p9"],
            "predicted_core": row[f"{prefix}_predicted_core"],
            "binding_percentile": row[f"{prefix}_binding_percentile_rank"],
            "protein": row[f"{prefix}_protein"],
            "accession": source.get("source_accession") or source.get("accession"),
            "source_certainty": row[f"{prefix}_source_certainty"],
        }
        model_count_value = row["left_model_count" if side == "ebv" else "right_model_count"]
        model_count = int(model_count_value) if str(model_count_value).strip() else 0
        if candidate_id in records and {
            key: value for key, value in records[candidate_id].items()
            if key not in {"surface_status", "model_count"}
        } != record:
            raise ValueError(f"inconsistent arm metadata for {allele} {candidate_id}")
        if candidate_id not in records:
            records[candidate_id] = {
                **record,
                "surface_status": "complete" if model_count == 5 else "missing",
                "model_count": model_count,
            }
        elif model_count > int(records[candidate_id]["model_count"]):
            records[candidate_id]["model_count"] = model_count
            records[candidate_id]["surface_status"] = "complete" if model_count == 5 else "missing"
    return [records[key] for key in sorted(records)]


def _control_census() -> list[dict[str, Any]]:
    common = {
        "census_date": "2026-08-28",
        "independent_v3_validation_vote": "false",
        "new_alphafold_modeling_allowed": "false",
    }
    return [
        {
            **common,
            "system_id": "SYS_BALF5_MBP_HY2E11",
            "tcr_id": "Hy.2E11",
            "hla_family": "DR",
            "hla_allotypes": "HLA-DRB5*01:01;HLA-DRB1*15:01",
            "ligand_sources": "EBV BALF5;human MBP",
            "admission_class": "strict_existing_development_control",
            "census_status": "not_new_and_already_used_in_v3_design",
            "distinct_source_rule": "pass",
            "exact_register_rule": "pass",
            "paired_structure_rule": "pass",
            "functional_both_arms_rule": "pass",
            "primary_source": "https://doi.org/10.1038/ni835",
            "pdb_ids": "1H15;1BX2",
            "exclusion_or_role_reason": "Valid strict system, but development use prevents an untouched V3 validation vote.",
        },
        {
            **common,
            "system_id": "SYS_ENGA_MBP_OB1A12",
            "tcr_id": "Ob.1A12",
            "hla_family": "DR",
            "hla_allotypes": "HLA-DRB1*15:01",
            "ligand_sources": "bacterial EngA;human MBP",
            "admission_class": "strict_existing_development_control",
            "census_status": "not_new_and_already_used_in_v3_design",
            "distinct_source_rule": "pass",
            "exact_register_rule": "pass",
            "paired_structure_rule": "pass",
            "functional_both_arms_rule": "pass",
            "primary_source": "https://doi.org/10.1016/j.immuni.2009.01.009",
            "pdb_ids": "2WBJ;1YMM",
            "exclusion_or_role_reason": "Valid strict system, but development use prevents an untouched V3 validation vote.",
        },
        {
            **common,
            "system_id": "SYS_MICROBIAL_MBP_HY1B11",
            "tcr_id": "Hy.1B11",
            "hla_family": "DQ",
            "hla_allotypes": "HLA-DQA1*01:02/HLA-DQB1*05:02",
            "ligand_sources": "HSV UL15;Pseudomonas PMM;human MBP",
            "admission_class": "strict_existing_development_control",
            "census_status": "not_new_and_already_used_in_v3_design",
            "distinct_source_rule": "pass",
            "exact_register_rule": "pass",
            "paired_structure_rule": "pass",
            "functional_both_arms_rule": "pass",
            "primary_source": "https://doi.org/10.1038/ncomms3623",
            "pdb_ids": "4MAY;4GRL;3PL6",
            "exclusion_or_role_reason": "Valid strict system with two required microbial pairs, but it contributes one development-system vote.",
        },
        {
            **common,
            "system_id": "LEAD_EBNA1_ANO2_2026",
            "tcr_id": "study clone set",
            "hla_family": "unresolved_for_exact_pair",
            "hla_allotypes": "unresolved",
            "ligand_sources": "EBV EBNA1;human ANO2",
            "admission_class": "prospective_lead",
            "census_status": "not_admitted",
            "distinct_source_rule": "pass",
            "exact_register_rule": "fail_unresolved",
            "paired_structure_rule": "fail_unresolved",
            "functional_both_arms_rule": "reported_but_exact_pair_unresolved",
            "primary_source": "https://pubmed.ncbi.nlm.nih.gov/41534529/",
            "pdb_ids": "",
            "exclusion_or_role_reason": "Exact paired peptides, registers, HLA complex, and structures are not yet resolved.",
        },
        {
            **common,
            "system_id": "LEAD_DQ8_A2_13_NATIVE_HIP",
            "tcr_id": "A2.13",
            "hla_family": "DQ",
            "hla_allotypes": "HLA-DQ8",
            "ligand_sources": "proinsulin C-peptide;hybrid insulin/NPY or IAPP peptide",
            "admission_class": "structure_resolved_lead",
            "census_status": "not_admitted",
            "distinct_source_rule": "fail_not_pathogen_vs_self",
            "exact_register_rule": "pass_structure_resolved",
            "paired_structure_rule": "pass",
            "functional_both_arms_rule": "pass",
            "primary_source": "https://doi.org/10.1016/j.jbc.2024.107612",
            "pdb_ids": "8VCX;8VCY;6XCP",
            "exclusion_or_role_reason": "Strong same-TCR cross-reactivity lead, but native and hybrid self ligands do not meet the frozen distinct-source rule.",
        },
        {
            **common,
            "system_id": "LEAD_DQ8_ET650_4_HIPS",
            "tcr_id": "ET650-4",
            "hla_family": "DQ",
            "hla_allotypes": "HLA-DQ8",
            "ligand_sources": "hybrid insulin/IAPP peptides",
            "admission_class": "structure_resolved_lead",
            "census_status": "not_admitted",
            "distinct_source_rule": "fail_not_pathogen_vs_self",
            "exact_register_rule": "pass_structure_resolved",
            "paired_structure_rule": "pass",
            "functional_both_arms_rule": "pass",
            "primary_source": "https://doi.org/10.1016/j.jbc.2024.107612",
            "pdb_ids": "8VD0;8VD2",
            "exclusion_or_role_reason": "Multiple recognized HIPs are informative, but they do not meet the frozen distinct-source rule.",
        },
        {
            **common,
            "system_id": "LEAD_DR15_TREE_NUT_ALLERGENS",
            "tcr_id": "multiple tetramer-isolated clones",
            "hla_family": "DR",
            "hla_allotypes": "includes HLA-DRB1*15:01",
            "ligand_sources": "cashew Ana o 1/Ana o 2;homologous tree-nut allergens",
            "admission_class": "functional_lead",
            "census_status": "not_admitted",
            "distinct_source_rule": "potential_pass",
            "exact_register_rule": "fail_not_locked",
            "paired_structure_rule": "fail_no_paired_structures",
            "functional_both_arms_rule": "pass_for_reported_clone_profiles",
            "primary_source": "https://pubmed.ncbi.nlm.nih.gov/27129138/",
            "pdb_ids": "",
            "exclusion_or_role_reason": "Potentially useful DR15 lead, but exact paired registers and matched ternary structures are unavailable.",
        },
    ]


def _control_search_audit() -> list[dict[str, Any]]:
    return [
        {
            "search_scope": "HLA-DRB1*03:01",
            "primary_record_checked": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9022012/",
            "result": "no_new_strict_system",
            "reason": "Proinsulin clones provide exact-allele functional data but not a same-TCR distinct-source, paired-structure system.",
        },
        {
            "search_scope": "HLA-DRB1*08:01",
            "primary_record_checked": "https://tcr3d.ibbr.umd.edu/class2",
            "result": "no_new_strict_system",
            "reason": "No exact-allele human multiligand system met all peptide, register, functional, and paired-structure fields.",
        },
        {
            "search_scope": "HLA-DRB1*13:03",
            "primary_record_checked": "https://doi.org/10.1038/s41467-025-66992-2",
            "result": "no_new_strict_system",
            "reason": "The inspected KSHV TCRs were not DRB1*13:03 cross-source systems; the paper explicitly reports lack of recognition for DRB1*13:03 in the relevant allotype test.",
        },
        {
            "search_scope": "HLA-DRB1*15:01",
            "primary_record_checked": "https://pubmed.ncbi.nlm.nih.gov/27129138/",
            "result": "functional_lead_only",
            "reason": "Tree-nut allergen cross-reactivity lacks locked exact registers and paired ternary structures.",
        },
        {
            "search_scope": "HLA-DQ family diversity",
            "primary_record_checked": "https://doi.org/10.1016/j.jbc.2024.107612",
            "result": "structure_resolved_leads_not_admitted",
            "reason": "A2.13 and ET650-4 are strong multiligand DQ8 systems but fail the frozen distinct-source requirement.",
        },
    ]


def _checksum_rows(out: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(out.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.csv":
            rows.append({
                "relative_path": str(path.relative_to(out)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            })
    return rows


def build_package(
    *,
    v3_dir: Path = DEFAULT_V3,
    panel_path: Path = DEFAULT_PANEL,
    out: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    target_summary = validate_frozen_targets(FROZEN_TARGETS)
    ranking_path = v3_dir / "v3_all_hla_ranked_pairs.csv"
    v3_manifest_path = v3_dir / "analysis_manifest.json"
    rows = read_csv(ranking_path)
    row_by_pair = {row["pair_id"]: row for row in rows}
    if len(row_by_pair) != len(rows):
        raise ValueError("V3 ranking input contains duplicate pair IDs")
    metadata = _panel_metadata(panel_path)
    target_rows = []
    for spec in FROZEN_TARGETS:
        if spec["pair_id"] not in row_by_pair:
            raise ValueError(f"missing frozen target {spec['pair_id']}")
        row = row_by_pair[spec["pair_id"]]
        if (
            row["allele"] != spec["allele"]
            or row["ebv_core_p1_p9"] != spec["ebv_core"]
            or row["self_core_p1_p9"] != spec["self_core"]
        ):
            raise ValueError(f"frozen target identity mismatch for {spec['target_id']}")
        if max(float(row["ebv_binding_percentile_rank"]), float(row["self_binding_percentile_rank"])) > 20:
            raise ValueError(f"frozen target fails binding eligibility for {spec['target_id']}")
        if row["surface_status"] != "complete" or int(row["left_model_count"]) != 5 or int(row["right_model_count"]) != 5:
            raise ValueError(f"frozen target lacks a complete five-model ensemble for {spec['target_id']}")
        target_rows.append({**row, **spec})

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    global_target_arm_ids = {
        str(row[field])
        for row in target_rows
        for field in ("ebv_candidate_id", "self_candidate_id")
    }
    confirmed_control_ids = {
        candidate_id
        for candidate_id, row in metadata.items()
        if _truth(row.get("required_for_confirmed_system"))
    }
    excluded_arm_ids = global_target_arm_ids | confirmed_control_ids

    panel_rows: list[dict[str, Any]] = []
    panel_summaries: list[dict[str, Any]] = []
    provenance_rows: list[dict[str, Any]] = []
    target_registry: list[dict[str, Any]] = []
    catalogs = {}
    for allele in sorted({row["allele"] for row in target_rows}):
        catalogs[(allele, "ebv")] = _arm_catalog(rows, metadata, allele=allele, side="ebv")
        catalogs[(allele, "self")] = _arm_catalog(rows, metadata, allele=allele, side="self")

    global_exclusion_feasibility = []
    for target in sorted(target_rows, key=lambda row: str(row["target_id"])):
        for side in ("ebv", "self"):
            count = _eligible_unique_count(catalogs[(target["allele"], side)], excluded_arm_ids)
            global_exclusion_feasibility.append({
                "target_id": target["target_id"],
                "allele": target["allele"],
                "arm_class": side,
                "eligible_unique_arm_count": count,
                "required_arm_count": 5,
                "global_exclusion_feasible": count >= 5,
                "assessment_uses_geometry": False,
                "resolution": "use_current_panel_target_exclusion_and_label_other_target_arm_overlap",
            })
    global_exclusion_feasible = all(
        int(row["eligible_unique_arm_count"]) >= 5 for row in global_exclusion_feasibility
    )

    for target in sorted(target_rows, key=lambda row: str(row["target_id"])):
        allele = str(target["allele"])
        ebv_target_arm = next(
            row for row in catalogs[(allele, "ebv")] if row["candidate_id"] == target["ebv_candidate_id"]
        )
        self_target_arm = next(
            row for row in catalogs[(allele, "self")] if row["candidate_id"] == target["self_candidate_id"]
        )
        panel_excluded_ids = confirmed_control_ids | {
            str(target["ebv_candidate_id"]),
            str(target["self_candidate_id"]),
        }
        ebv_selected, ebv_provenance = select_comparator_arms(
            ebv_target_arm,
            catalogs[(allele, "ebv")],
            allele=allele,
            arm_class="ebv",
            excluded_candidate_ids=panel_excluded_ids,
            count=5,
            seed=SEED,
        )
        self_selected, self_provenance = select_comparator_arms(
            self_target_arm,
            catalogs[(allele, "self")],
            allele=allele,
            arm_class="self",
            excluded_candidate_ids=panel_excluded_ids,
            count=5,
            seed=SEED,
        )
        for row in ebv_provenance + self_provenance:
            provenance_rows.append({
                **row,
                "target_id": target["target_id"],
                "panel_id": f"PANEL_{target['target_id']}",
                "target_pair_id": target["pair_id"],
                "overlaps_other_frozen_target_arm": (
                    str(row["candidate_id"]) in global_target_arm_ids
                    and str(row["candidate_id"]) not in {
                        str(target["ebv_candidate_id"]), str(target["self_candidate_id"])
                    }
                ),
            })
        pair_lookup = {
            (row["ebv_candidate_id"], row["self_candidate_id"]): row
            for row in rows
            if row["allele"] == allele
        }
        panel = build_n3_panel(target, ebv_selected, self_selected, pair_lookup)
        ranked = rank_panel_rows(panel, seed=SEED)
        for row in ranked:
            row.update({
                "panel_id": f"PANEL_{target['target_id']}",
                "target_id": target["target_id"],
                "target_lane": target["lane"],
                "target_pair_id": target["pair_id"],
                "register_interpretation": (
                    "register_robust_declared_geometry"
                    if target["lane"] == "register" and _truth(target["register_robust"])
                    else "sensitivity_only_register_uncertain"
                ),
                "contains_other_frozen_target_arm": (
                    row["row_role"] == "n3"
                    and (
                        str(row["ebv_candidate_id"]) in global_target_arm_ids
                        or str(row["self_candidate_id"]) in global_target_arm_ids
                    )
                ),
            })
            panel_rows.append(row)
        target_result = next(row for row in ranked if row["row_role"] == "target")
        best_n3 = min((row for row in ranked if row["row_role"] == "n3"), key=lambda row: int(row["panel_primary_rank"]))
        summary = {
            "panel_id": f"PANEL_{target['target_id']}",
            "target_id": target["target_id"],
            "allele": allele,
            "lane": target["lane"],
            "target_pair_id": target["pair_id"],
            "target_ebv_core": target["ebv_core"],
            "target_self_core": target["self_core"],
            "panel_status": "complete",
            "panel_pair_count": 26,
            "n3_pair_count": 25,
            "target_primary_rank": target_result["panel_primary_rank"],
            "target_local_surface_rank": target_result["panel_local_surface_rank"],
            "target_exposed_backbone_rank": target_result["panel_exposed_backbone_rank"],
            "target_full_core_rmsd_rank": target_result["panel_full_core_rmsd_rank"],
            "target_anchor_rmsd_rank": target_result["panel_anchor_rmsd_rank"],
            "target_physicochemical_rank": target_result["panel_physicochemical_rank"],
            "target_identity_rank": target_result["panel_identity_rank"],
            "target_random_rank": target_result["panel_random_rank"],
            "target_global_v3_rank": target["primary_rank"],
            "target_blosum": target["tcr_facing_blosum62_similarity"],
            "best_n3_pair_id": best_n3["pair_id"],
            "best_n3_blosum": best_n3["tcr_facing_blosum62_similarity"],
            "target_blosum_margin_over_best_n3": round(
                float(target["tcr_facing_blosum62_similarity"])
                - float(best_n3["tcr_facing_blosum62_similarity"]),
                12,
            ),
            "target_register_robust": _truth(target["register_robust"]),
            "register_interpretation": (
                "register_robust_declared_geometry"
                if target["lane"] == "register" and _truth(target["register_robust"])
                else "sensitivity_only_register_uncertain"
            ),
            "target_surface_ensemble_uncertainty": target["surface_ensemble_uncertainty"],
            "target_surface_ensemble_stable": target["surface_ensemble_stable"],
            "claim_boundary": CLAIM_BOUNDARY,
        }
        panel_summaries.append(summary)
        target_registry.append({
            **{key: target[key] for key in (
                "target_id", "allele", "lane", "pair_id", "ebv_core", "self_core",
                "ebv_candidate_id", "self_candidate_id", "ebv_protein", "self_protein",
                "ebv_sequence", "self_sequence", "ebv_binding_percentile_rank",
                "self_binding_percentile_rank", "register_robust", "surface_status",
                "left_model_count", "right_model_count",
            )},
            "binding_eligibility_rule": "both_arms_percentile_at_or_below_20_not_used_in_rank",
            "comparator_target_exclusion_scope": "current_panel_target_arms_plus_confirmed_control_ligands",
            "other_frozen_target_arms_may_be_labeled_n3": True,
            "claim_boundary": CLAIM_BOUNDARY,
        })

    ranking_gate = build_ranking_context_gate(panel_summaries)
    specificity_gate = build_specificity_gate([])
    census = _control_census()
    strict_registry = [row for row in census if row["admission_class"] == "strict_existing_development_control"]
    definitive_gate = {
        "status": "not_evaluable_insufficient_untouched_strict_systems",
        "current_independent_v3_validation_system_count": 0,
        "new_strict_system_count": 0,
        "required_independent_system_count": 6,
        "target_independent_system_count": 8,
        "required_hla_family_count": 2,
        "existing_three_system_role": "development_controls_only",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    af_status = {
        "status": "not_created_no_new_strict_systems",
        "new_strict_system_count": 0,
        "model_seeds_if_admitted": [271828, 314159],
        "maximum_jobs_per_batch": 30,
        "submission_status": "not_submitted",
        "prepared_not_submitted": False,
        "reason": "Only newly verified strict controls are eligible for new AlphaFold modeling.",
        "discovery_candidate_models_reused": True,
    }
    protocol = {
        "protocol_version": "HIGH_YIELD_CONTROL_VALIDATION_2026-08-28",
        "frozen_before_control_geometry_read": True,
        "seed": SEED,
        "target_count": 12,
        "targets_per_hla": 3,
        "comparators_per_arm": 5,
        "n3_pairs_per_panel": 25,
        "panel_pair_count": 26,
        "binding_percentile_eligibility_maximum": 20,
        "binding_percentile_used_in_rank": False,
        "comparator_order": [
            "peptide_length_difference_asc",
            "binding_percentile_difference_asc",
            "seeded_hash_asc",
            "candidate_id_asc",
        ],
        "ranking_order": [
            "tcr_facing_blosum62_desc",
            "tcr_face_physicochemical_mismatch_asc",
            "tcr_facing_sequence_identity_desc",
            "panel_local_surface_percentile_asc_when_register_robust",
            "pair_id_asc",
        ],
        "surface_features": list(SURFACE_FEATURES),
        "n3_specificity_role": "unknown_recognition_not_a_specificity_negative",
        "global_frozen_target_arm_exclusion_feasible": global_exclusion_feasible,
        "comparator_target_exclusion_scope": "current_panel_target_arms_plus_confirmed_control_ligands",
        "protocol_resolution_reason": (
            "Global frozen-target-arm exclusion leaves fewer than five eligible unique self arms; "
            "the current target arms and confirmed-control ligands remain excluded, while any other "
            "high-yield arm is explicitly labeled as unknown-recognition N3 overlap."
        ),
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "v3_input": {
            "relative_path": str(ranking_path.relative_to(ROOT)),
            "sha256": sha256_file(ranking_path),
            "bytes": ranking_path.stat().st_size,
        },
        "v3_manifest": {
            "relative_path": str(v3_manifest_path.relative_to(ROOT)),
            "sha256": sha256_file(v3_manifest_path),
            "bytes": v3_manifest_path.stat().st_size,
        },
        "frozen_panel_input": {
            "relative_path": str(panel_path.relative_to(ROOT)),
            "sha256": sha256_file(panel_path),
            "bytes": panel_path.stat().st_size,
        },
        "frozen_targets": list(FROZEN_TARGETS),
        "claim_boundary": CLAIM_BOUNDARY,
    }

    write_json(out / "protocol_lock.json", protocol)
    write_csv(out / "frozen_target_registry.csv", target_registry)
    write_csv(out / "comparator_provenance.csv", provenance_rows)
    write_csv(out / "global_exclusion_feasibility.csv", global_exclusion_feasibility)
    write_csv(out / "panel_feature_matrix.csv", panel_rows)
    write_csv(out / "panel_rank_summary.csv", panel_summaries)
    write_csv(out / "control_system_census.csv", census)
    write_csv(out / "control_search_audit.csv", _control_search_audit())
    write_csv(out / "strict_control_registry.csv", strict_registry)
    write_json(out / "ranking_context_gate.json", ranking_gate)
    write_json(out / "specificity_gate.json", specificity_gate)
    write_json(out / "definitive_validation_gate.json", definitive_gate)
    write_json(out / "alphafold_manifest_status.json", af_status)

    supportive = sum(int(row["target_primary_rank"]) <= 3 for row in panel_summaries)
    readme = f"""# High-Yield Control Validation

This additive package evaluates 12 frozen discovery candidates against 25
score-blind, exact-HLA N3 comparison pairs each. It reuses the existing V3
model-derived features and recomputes all panel-level percentiles and ranks.

## Result

- Complete panels: 12/12
- Total rows: 312 (12 targets and 300 N3 comparisons)
- Targets ranked in the top three: {supportive}/12
- New untouched strict positive-control systems admitted: 0
- Specificity gate: not evaluable; no explicit N1/N2 registry
- Discovery unlock: false

N3 means recognition is unknown. These rows are fair computational ranking
comparators, not biological negative controls. Sequence-lane target structure
is retained only as register-sensitivity diagnostics and abstains from the V3
primary structural tie-break.

The originally proposed global exclusion of every frozen target arm was not
feasible: fewer than five unique eligible self arms remained. The frozen
resolution excludes the current panel's target arms and all confirmed-control
ligands, while explicitly flagging N3 rows that contain an arm from another
high-yield target. The feasibility audit uses no geometry.

## Primary-source census

The existing Hy.2E11, Ob.1A12, and Hy.1B11 systems remain development
controls. DQ8 A2.13/ET650-4 systems are retained as structure-resolved leads
but fail the frozen distinct-source rule. Exact-allele searches for the four
studied DR alleles produced no newly admissible system with all strict fields.

{CLAIM_BOUNDARY}
"""
    (out / "README.md").write_text(readme, encoding="utf-8")
    results_lines = [
        "# Panel Results",
        "",
        "| HLA | Lane | EBV core | Self core | V3 panel rank | Surface rank | Register interpretation |",
        "|---|---|---|---|---:|---:|---|",
    ]
    for row in sorted(panel_summaries, key=lambda value: (value["allele"], value["lane"], value["target_id"])):
        results_lines.append(
            f"| {row['allele']} | {row['lane']} | {row['target_ebv_core']} | {row['target_self_core']} | "
            f"{row['target_primary_rank']} | {row['target_local_surface_rank']} | {row['register_interpretation']} |"
        )
    results_lines.extend(["", CLAIM_BOUNDARY, ""])
    (out / "RESULTS_SUMMARY.md").write_text("\n".join(results_lines), encoding="utf-8")

    manifest = {
        "package": "high_yield_control_validation_2026-08-28",
        "status": "complete_additive_n3_ranking_context",
        "target_count": target_summary["target_count"],
        "targets_per_hla": target_summary["targets_per_hla"],
        "panel_count": len(panel_summaries),
        "panel_row_count": len(panel_rows),
        "target_row_count": sum(row["row_role"] == "target" for row in panel_rows),
        "n3_row_count": sum(row["row_role"] == "n3" for row in panel_rows),
        "supportive_rank_context_count": supportive,
        "new_strict_control_count": 0,
        "new_alphafold_job_count": 0,
        "existing_v1_v2_v3_outputs_modified": False,
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "specificity_claim_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "analysis_manifest.json", manifest)
    write_csv(out / "SHA256SUMS.csv", _checksum_rows(out))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v3-dir", type=Path, default=DEFAULT_V3)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    manifest = build_package(v3_dir=args.v3_dir, panel_path=args.panel, out=args.out)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
