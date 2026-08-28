"""Build the frozen, prepared-not-submitted HLA-II benchmark v2 pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
from typing import Any, Mapping, Sequence

import numpy as np

from build_hla2_positive_control_benchmark import _hla_sequences, curated_registry
from hla2_positive_control_benchmark import (
    build_af3_job_batches,
    build_pdb_oracle_pairings,
    validate_af3_job_package,
    validate_comparator_registry,
    validate_control_registry,
)
from hla2_positive_control_benchmark_v2 import (
    CLAIM_BOUNDARY_V2,
    PILOT_SEEDS,
    build_definitive_ranking_gate,
    build_oracle_availability,
    build_pilot_attribution_gate,
    build_protocol_lock,
    validate_specificity_registry,
)


ROOT = Path(__file__).resolve().parents[1]
V1_PACKAGE = ROOT / "processed/hla2_positive_control_benchmark_2026-08-25"
V1_RESULTS = ROOT / "processed/hla2_positive_control_benchmark_results_2026-08-26"
TCELL_PACKAGE = ROOT / "processed/tcell_library_v2_2026-08-22"
DEFAULT_OUT = ROOT / "processed/hla2_positive_control_benchmark_v2_pilot_2026-08-26"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] = (),
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


def stable_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _family(beta_allele: str) -> str:
    if "DRB" in beta_allele:
        return "DR"
    if "DQB" in beta_allele:
        return "DQ"
    if "DPB" in beta_allele:
        return "DP"
    raise ValueError(f"unsupported HLA-II beta allotype: {beta_allele}")


def _hy2_comparators(out: Path) -> list[dict[str, Any]]:
    controls = read_csv(TCELL_PACKAGE / "frozen_native_hla_controls.csv")
    predictions = {
        (row["allele"], row["candidate_id"]): row
        for row in read_csv(TCELL_PACKAGE / "calibration_control_binding_predictions.csv")
    }
    raw_files = {
        "viral": "iedb_calibration_ebv_drb5_0101.tsv",
        "self": "iedb_calibration_human_non_cns_drb1_1501.tsv",
    }
    for filename in raw_files.values():
        _copy_file(
            TCELL_PACKAGE / "raw_responses" / filename,
            out / "sources/raw_responses" / filename,
        )
    rows = []
    for control in controls:
        arm = control["arm"]
        alpha = "HLA-DRA*01:01"
        beta = "HLA-DRB5*01:01" if arm == "viral" else "HLA-DRB1*15:01"
        prediction = predictions[(beta, control["control_candidate_id"])]
        rows.append({
            "system_id": "SYS_BALF5_MBP_HY2E11",
            "positive_pair_id": "PAIR_HY2E11_BALF5_MBP",
            "candidate_id": control["control_candidate_id"],
            "comparator_arm": "microbial" if arm == "viral" else "self",
            "sequence": control["control_sequence"],
            "predicted_core": prediction["predicted_core"],
            "core_start_1_based": int(prediction["core_start"]),
            "register_resolution": prediction["register_resolution"],
            "seq_num": int(prediction["seq_num"]),
            "raw_response_file": f"sources/raw_responses/{raw_files[arm]}",
            "prediction_method": prediction["prediction_method"],
            "prediction_status": prediction["prediction_status"],
            "binding_percentile": prediction["percentile_rank"],
            "mhc_alpha_allele": alpha,
            "mhc_beta_allele": beta,
            "negative_tier": "N3",
            "recognition_status": "unknown_not_specificity_negative",
            "source_protein": control["control_source"],
            "source_accession": control["control_accession"],
            "selection_length_difference": control["length_difference"],
            "selection_binding_percentile_difference": "frozen_rank_bin_match",
            "selection_seeded_hash": hashlib.sha256(
                f"HY2_FROZEN_SCORE_BLIND|{control['control_candidate_id']}".encode("utf-8")
            ).hexdigest(),
            "selection_is_score_blind": True,
            "selection_provenance": "v1_frozen_before_geometry_reused_as_identity_only",
        })
    return rows


def _other_comparators(out: Path) -> list[dict[str, Any]]:
    rows = read_csv(V1_PACKAGE / "controls/control_decoy_registry.csv")
    referenced = sorted({row["raw_response_file"] for row in rows})
    for relative in referenced:
        _copy_file(V1_PACKAGE / relative, out / relative)
    for row in rows:
        row["selection_provenance"] = "v1_score_blind_identity_reused_before_v2_geometry"
    return rows


def _comparison_universe(
    pairs: Sequence[Mapping[str, Any]],
    comparators: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows = []
    for pair in sorted(pairs, key=lambda row: str(row["pair_id"])):
        pair_id = str(pair["pair_id"])
        left = sorted(
            (row for row in comparators if row["positive_pair_id"] == pair_id and row["comparator_arm"] == "microbial"),
            key=lambda row: str(row["candidate_id"]),
        )
        right = sorted(
            (row for row in comparators if row["positive_pair_id"] == pair_id and row["comparator_arm"] == "self"),
            key=lambda row: str(row["candidate_id"]),
        )
        if len(left) != 5 or len(right) != 5:
            raise ValueError(f"pair {pair_id} does not have five comparators per arm")
        for seed in PILOT_SEEDS:
            rows.append({
                "system_id": pair["system_id"],
                "positive_pair_id": pair_id,
                "panel_seed": seed,
                "pair_id": f"s{seed}|{pair['left_ligand_id']}|{pair['right_ligand_id']}",
                "comparison_role": "positive",
                "negative_tier": "positive",
                "left_id": pair["left_ligand_id"],
                "right_id": pair["right_ligand_id"],
                "selection_is_score_blind": True,
            })
            for left_row in left:
                for right_row in right:
                    rows.append({
                        "system_id": pair["system_id"],
                        "positive_pair_id": pair_id,
                        "panel_seed": seed,
                        "pair_id": f"s{seed}|{left_row['candidate_id']}|{right_row['candidate_id']}",
                        "comparison_role": "N3_pair_decoy",
                        "negative_tier": "N3",
                        "left_id": left_row["candidate_id"],
                        "right_id": right_row["candidate_id"],
                        "selection_is_score_blind": True,
                    })
    return rows


def _unique_af_ligands(
    ligands: Sequence[Mapping[str, Any]],
    comparators: Sequence[Mapping[str, Any]],
    comparisons: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    ligand_by_id = {str(row["ligand_id"]): row for row in ligands}
    comparator_by_pair_id = {
        (str(row["positive_pair_id"]), str(row["candidate_id"])): row for row in comparators
    }
    candidates = []
    for comparison in comparisons:
        pair_id = str(comparison["positive_pair_id"])
        for side in ("left", "right"):
            entity_id = str(comparison[f"{side}_id"])
            if entity_id in ligand_by_id:
                source = ligand_by_id[entity_id]
                value = {
                    "source_entity_id": entity_id,
                    "system_id": source["system_id"],
                    "ligand_role": "positive_control",
                    "sequence": source["sequence"],
                    "mhc_alpha_allele": source["mhc_alpha_allele"],
                    "mhc_beta_allele": source["mhc_beta_allele"],
                    "core_sequence": source["core"],
                    "core_start_1_based": source["core_start_1_based"],
                    "register_resolution": "experimentally_resolved_structure",
                    "register_source": f"PDB_{source['pdb_id']}",
                    "seq_num": "",
                    "raw_response_file": "",
                }
            else:
                source = comparator_by_pair_id[(pair_id, entity_id)]
                value = {
                    "source_entity_id": entity_id,
                    "system_id": source["system_id"],
                    "ligand_role": "N3_comparator",
                    "sequence": source["sequence"],
                    "mhc_alpha_allele": source["mhc_alpha_allele"],
                    "mhc_beta_allele": source["mhc_beta_allele"],
                    "core_sequence": source["predicted_core"],
                    "core_start_1_based": source["core_start_1_based"],
                    "register_resolution": source["register_resolution"],
                    "register_source": source["prediction_method"],
                    "seq_num": source["seq_num"],
                    "raw_response_file": source["raw_response_file"],
                }
            value["pmhc_key"] = "|".join((
                str(value["mhc_alpha_allele"]), str(value["mhc_beta_allele"]), str(value["sequence"]),
            ))
            candidates.append(value)
    by_key: dict[str, dict[str, Any]] = {}
    memberships: dict[str, set[str]] = {}
    for row in sorted(candidates, key=lambda value: (
        value["pmhc_key"], value["ligand_role"] != "positive_control", value["source_entity_id"]
    )):
        key = str(row["pmhc_key"])
        memberships.setdefault(key, set()).add(str(row["system_id"]))
        current = by_key.get(key)
        if current is None or (
            current["ligand_role"] != "positive_control" and row["ligand_role"] == "positive_control"
        ):
            by_key[key] = dict(row)
    output = []
    for key, row in sorted(by_key.items()):
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:16]
        row["ligand_id"] = f"V2PMHC_{digest}"
        row["system_id"] = ";".join(sorted(memberships[key]))
        row["panel_membership_count"] = sum(
            entity["pmhc_key"] == key for entity in candidates
        )
        output.append(row)
    pmhc_by_key = {row["pmhc_key"]: row["ligand_id"] for row in output}
    for comparison in comparisons:
        pair_id = str(comparison["positive_pair_id"])
        for side in ("left", "right"):
            entity_id = str(comparison[f"{side}_id"])
            if entity_id in ligand_by_id:
                source = ligand_by_id[entity_id]
                key = "|".join((
                    str(source["mhc_alpha_allele"]), str(source["mhc_beta_allele"]), str(source["sequence"]),
                ))
            else:
                source = comparator_by_pair_id[(pair_id, entity_id)]
                key = "|".join((
                    str(source["mhc_alpha_allele"]), str(source["mhc_beta_allele"]), str(source["sequence"]),
                ))
            comparison[f"{side}_pmhc_id"] = pmhc_by_key[key]
    return output


def _frozen_oracle_pairings(
    positive_pairs: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Rebuild score-blind pairing identities without opening any geometry table."""
    structural = read_csv(V1_PACKAGE / "benchmark/pdb_structural_ligand_registry.csv")
    structural_by_id = {row["ligand_id"]: row for row in structural}
    fields = (
        "system_id", "positive_pair_id", "pair_id", "pair_role", "left_ligand_id",
        "right_ligand_id", "left_pdb_id", "right_pdb_id", "selection_is_score_blind",
    )
    output = []
    for pair in positive_pairs:
        positive_ids = {str(pair["left_ligand_id"]), str(pair["right_ligand_id"])}
        unique_pmhc: dict[tuple[str, str, str], Mapping[str, Any]] = {}
        for row in sorted(structural, key=lambda value: (
            str(value["ligand_id"]) not in positive_ids,
            float(value["resolution_A"]),
            str(value["ligand_id"]),
        )):
            key = (
                str(row["mhc_alpha_allele"]), str(row["mhc_beta_allele"]),
                str(row["core_sequence"]),
            )
            unique_pmhc.setdefault(key, row)
        left = structural_by_id[str(pair["left_ligand_id"])]
        right = structural_by_id[str(pair["right_ligand_id"])]
        pairings, _summary = build_pdb_oracle_pairings({
            **dict(pair),
            "left_mhc_alpha_allele": left["mhc_alpha_allele"],
            "left_mhc_beta_allele": left["mhc_beta_allele"],
            "right_mhc_alpha_allele": right["mhc_alpha_allele"],
            "right_mhc_beta_allele": right["mhc_beta_allele"],
        }, list(unique_pmhc.values()))
        for row in pairings:
            value = {
                **dict(row),
                "left_pdb_id": structural_by_id[str(row["left_ligand_id"])]["pdb_id"],
                "right_pdb_id": structural_by_id[str(row["right_ligand_id"])]["pdb_id"],
            }
            output.append({field: value[field] for field in fields})
    return output


def _oracle_availability(pairings: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    by_pair: dict[str, list[Mapping[str, Any]]] = {}
    for row in pairings:
        by_pair.setdefault(str(row["positive_pair_id"]), []).append(row)
    rows = []
    for pair_id, pair_rows in sorted(by_pair.items()):
        source = pair_rows[0]
        rows.append({
            "system_id": source["system_id"],
            "positive_pair_id": pair_id,
            "eligible_decoy_count": sum(row["pair_role"] != "positive" for row in pair_rows),
            "availability_source": "pre_v2_score_blind_exact_hla_pairing_identities",
            "prior_rank_reused_for_v2_tuning": False,
            "prior_geometry_reused_for_v2_tuning": False,
        })
    return build_oracle_availability(rows)


def _checksums(out: Path) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(
        (value for value in out.rglob("*") if value.is_file() and value.name != "SHA256SUMS.csv"),
        key=str,
    ):
        rows.append({
            "relative_path": str(path.relative_to(out)),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        })
    return rows


def build_pilot_package(out: Path = DEFAULT_OUT) -> dict[str, Any]:
    systems, ligands, pairs, sources = curated_registry()
    ligand_by_system: dict[str, list[Mapping[str, Any]]] = {}
    for ligand in ligands:
        ligand_by_system.setdefault(str(ligand["system_id"]), []).append(ligand)
    for system in systems:
        if system["eligibility"] == "strict":
            families = {_family(str(row["mhc_beta_allele"])) for row in ligand_by_system[system["system_id"]]}
            system["hla_family"] = ";".join(sorted(families))
            system["distinct_biological_sources_verified"] = True
        else:
            system["hla_family"] = "unresolved"
            system["distinct_biological_sources_verified"] = False
    registry_validation = validate_control_registry(systems, ligands, pairs)

    comparators = [*_hy2_comparators(out), *_other_comparators(out)]
    comparators.sort(key=lambda row: (
        str(row["positive_pair_id"]), str(row["comparator_arm"]), str(row["candidate_id"])
    ))
    comparator_validation = validate_comparator_registry(
        comparators, expected_pair_ids=[row["pair_id"] for row in pairs]
    )
    comparisons = _comparison_universe(pairs, comparators)
    oracle_pairings = _frozen_oracle_pairings(pairs)
    af_ligands = _unique_af_ligands(ligands, comparators, comparisons)
    jobs, job_manifest, batches = build_af3_job_batches(
        af_ligands, _hla_sequences(), panel_seeds=PILOT_SEEDS, batch_size=30
    )
    job_validation = validate_af3_job_package(batches, job_manifest, _hla_sequences())

    registry_payload = {"systems": systems, "ligands": ligands, "pairs": pairs, "sources": sources}
    protocol = build_protocol_lock(
        strict_system_ids=registry_validation["strict_system_ids"],
        positive_pair_ids=[row["pair_id"] for row in pairs],
        registry_sha256=stable_sha256(registry_payload),
        comparator_sha256=stable_sha256(comparators),
        oracle_pairings_sha256=stable_sha256(oracle_pairings),
        software_versions={
            "python": platform.python_version(),
            "numpy": np.__version__,
            "alphafold_server_json_schema": "dialect_alphafoldserver_version_1",
            "geometry_module_sha256": sha256_file(ROOT / "src/hla2_positive_control_benchmark.py"),
            "v2_contract_module_sha256": sha256_file(ROOT / "src/hla2_positive_control_benchmark_v2.py"),
        },
    )

    write_json(out / "protocol/protocol_lock.json", protocol)
    write_csv(out / "registry/control_system_registry.csv", systems)
    write_csv(out / "registry/control_ligand_registry.csv", ligands)
    write_csv(out / "registry/positive_pair_registry.csv", pairs)
    write_csv(out / "registry/literature_and_structure_sources.csv", sources)
    write_csv(out / "registry/specificity_negative_registry.csv", [], fields=(
        "negative_id", "tcr_system", "peptide", "exact_hla", "assay", "tested_condition",
        "outcome", "negative_tier", "source_location", "citation",
    ))
    write_csv(out / "controls/control_decoy_registry.csv", comparators)
    write_csv(out / "controls/comparison_universe.csv", comparisons)
    write_csv(out / "benchmark/pdb_oracle_frozen_pairings.csv", oracle_pairings)
    write_csv(out / "alphafold_jobs/job_manifest.csv", job_manifest)
    write_csv(out / "alphafold_jobs/unique_pmhc_inventory.csv", af_ligands)
    for index, batch in enumerate(batches, start=1):
        write_json(
            out / f"alphafold_jobs/hla2_v2_pilot_batch_{index:02d}_{len(batch)}_jobs.json",
            batch,
        )

    oracle = _oracle_availability(oracle_pairings)
    write_csv(out / "benchmark/pdb_oracle_availability.csv", oracle)
    pilot_gate = build_pilot_attribution_gate(
        [], required_system_ids=registry_validation["strict_system_ids"]
    )
    definitive_gate = build_definitive_ranking_gate([], systems)
    specificity_gate = validate_specificity_registry([])
    write_json(out / "benchmark/pilot_attribution_gate.json", pilot_gate)
    write_json(out / "benchmark/definitive_ranking_gate.json", definitive_gate)
    write_json(out / "benchmark/specificity_gate.json", specificity_gate)

    v1_snapshot = {
        "v1_files_modified": False,
        "v1_missing_jobs_preserved": True,
        "v1_manifest_sha256": sha256_file(V1_PACKAGE / "analysis_manifest.json"),
        "v1_results_manifest_sha256": sha256_file(V1_RESULTS / "analysis_manifest.json"),
        "v1_trust_gate_sha256": sha256_file(V1_RESULTS / "benchmark/trust_gate.json"),
    }
    write_json(out / "validation/v1_immutability_snapshot.json", v1_snapshot)
    verification = {
        "benchmark_version": protocol["benchmark_version"],
        "strict_system_count": len(registry_validation["strict_system_ids"]),
        "positive_pair_count": len(pairs),
        "comparator_row_count": len(comparators),
        "comparators_per_pair_arm": 5,
        "comparison_row_count": len(comparisons),
        "pair_decoys_per_panel": 25,
        "frozen_pdb_oracle_pairing_count": len(oracle_pairings),
        "unique_pmhc_count": len(af_ligands),
        "job_count": len(jobs),
        "batch_sizes": [len(batch) for batch in batches],
        "all_jobs_prepared_not_submitted": all(
            row["status"] == "prepared_not_submitted" for row in job_manifest
        ),
        "all_job_names_unique": len({row["job_name"] for row in job_manifest}) == len(job_manifest),
        "job_package_validation": job_validation,
        "comparator_registry_validation": comparator_validation,
        "discovery_files_read": False,
        "discovery_files_written": False,
        "cross_allele_consensus_created": False,
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY_V2,
    }
    write_json(out / "validation/package_verification_summary.json", verification)
    write_json(out / "analysis_manifest.json", {
        "benchmark_version": protocol["benchmark_version"],
        "protocol_sha256": protocol["protocol_sha256"],
        "status": "prepared_not_submitted",
        "strict_system_count": verification["strict_system_count"],
        "positive_pair_count": verification["positive_pair_count"],
        "job_count": verification["job_count"],
        "batch_sizes": verification["batch_sizes"],
        "pilot_attribution_status": pilot_gate["pilot_attribution_status"],
        "definitive_status": definitive_gate["definitive_status"],
        "specificity_status": specificity_gate["specificity_status"],
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "claim_boundary": CLAIM_BOUNDARY_V2,
    })
    (out / "README.md").write_text(
        "# HLA-II benchmark v2 attribution pilot\n\n"
        "This is a fresh, frozen, prepared-not-submitted three-system attribution pilot. "
        "It does not reuse v1 AlphaFold jobs, cannot freeze weights, and cannot unlock discovery.\n\n"
        "The definitive benchmark remains blocked until at least six independent strict systems "
        "spanning at least two HLA-II families are admitted.\n\n"
        f"{CLAIM_BOUNDARY_V2}\n",
        encoding="utf-8",
    )
    (out / "METHODS.md").write_text(
        "# Frozen methods\n\n"
        "Each required positive is ranked against 25 score-blind N3 pair decoys assembled from "
        "five exact-HLA comparators per arm. Fresh pMHC-only AlphaFold jobs use fixed seeds "
        "271828 and 314159. TCR-facing features use P2/P3/P5/P7/P8; anchors use P1/P4/P6/P9. "
        "Binding percentile and length/register agreement are matching diagnostics only.\n\n"
        "The pilot requires top-3 recovery in every panel, improvement over the strongest "
        "training-selected nonstructural baseline in a majority of systems, no worsened system, "
        "structural weight of at least 0.25 on every credited improvement, and ablation removal "
        "of that improvement. A supportive pilot still cannot unlock discovery.\n",
        encoding="utf-8",
    )
    write_csv(out / "SHA256SUMS.csv", _checksums(out))
    return {
        "output_directory": str(out),
        "strict_system_count": verification["strict_system_count"],
        "positive_pair_count": verification["positive_pair_count"],
        "comparator_row_count": verification["comparator_row_count"],
        "comparison_row_count": verification["comparison_row_count"],
        "unique_pmhc_count": verification["unique_pmhc_count"],
        "job_count": verification["job_count"],
        "batch_sizes": verification["batch_sizes"],
        "pilot_attribution_status": pilot_gate["pilot_attribution_status"],
        "definitive_status": definitive_gate["definitive_status"],
        "specificity_status": specificity_gate["specificity_status"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(build_pilot_package(args.out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
