"""Audit every saved pMHC modeling phase and rebuild the same-register score sheet.

The audit keeps modeling phases separate, removes exact duplicated AF3 jobs by
candidate/peptide/seed, and evaluates only frozen shortlist pairs with resolved
primary-allele P1--P9 hypotheses.  Outputs are prioritization evidence, not
evidence of TCR binding, cross-reactivity, or disease mechanism.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any

from analyze_af3_pmhc_downloads import (
    analyze_complete_job,
    build_candidate_metadata,
    candidate_id_from_request_name,
    parse_mmcif,
)
from premeeting_rigor import binding_rank_bin
from same_register_af3_analysis import (
    VALID_REGISTER_STATUSES,
    _metric_distribution,
    direct_register_sequence_metrics,
    matched_background_feasibility,
    same_register_geometry,
)


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "complete_model_pipeline_audit_2026-08-15"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def request_details(path: Path) -> dict[str, Any]:
    request = json.loads(path.read_text(encoding="utf-8"))
    job = request[0] if isinstance(request, list) else request
    chains = [
        entry["proteinChain"]["sequence"]
        for entry in job.get("sequences", [])
        if "proteinChain" in entry
    ]
    seed = str((job.get("modelSeeds") or [""])[0])
    name = str(job.get("name", ""))
    return {
        "request_name": name,
        "candidate_id": candidate_id_from_request_name(name),
        "server_seed": seed,
        "protein_chain_count": len(chains),
        "requested_peptide": chains[2] if len(chains) == 3 else "",
    }


def monday_root(project_root: Path = ROOT) -> Path:
    candidates = [
        path
        for path in project_root.glob("Alphafold3_pMHC_Folds*")
        if path.is_dir() and len(list(path.glob("folds_2026_08_10*"))) >= 4
    ]
    if not candidates:
        return project_root / "Alphafold3_pMHC_Folds _v1"
    if len(candidates) != 1:
        raise FileNotFoundError(f"Expected one Monday AF3 root, found {candidates}")
    return candidates[0]


def af3_discovery_roots(project_root: Path) -> list[tuple[str, Path, str, int]]:
    """Return the preserved AF3 phases, including the Seed 03 completion batch."""
    return [
        ("monday_af3_full", monday_root(project_root), "folds_2026_08_10*/*", 0),
        ("focused_af3_rerun", project_root / "Alphafold3_pMHCs", "json[12]folds/*", 1),
        ("new_background_af3", project_root / "Alphafold3_pMHCs", "json3folds/*", 2),
        (
            "new_background_af3_seed03_completion",
            project_root / "folds_2026_08_15_01_53",
            "*",
            3,
        ),
    ]


def discover_jobs() -> list[dict[str, Any]]:
    roots = af3_discovery_roots(ROOT)
    metadata = build_candidate_metadata()
    discovered: list[dict[str, Any]] = []
    for phase, root, pattern, phase_priority in roots:
        for directory in sorted(root.glob(pattern)):
            if not directory.is_dir():
                continue
            requests = list(directory.glob("*_job_request.json"))
            if len(requests) != 1:
                continue
            details = request_details(requests[0])
            candidate_id = str(details["candidate_id"])
            if int(details["protein_chain_count"]) != 3:
                cohort = "excluded_non_pmhc_chain_layout"
            elif candidate_id in metadata:
                cohort = str(metadata[candidate_id]["af3_cohort"])
            elif candidate_id.startswith("GLIALCAM_"):
                cohort = "special_glialcam_control"
            else:
                cohort = "excluded_unmapped_decoy_or_other"
            discovered.append({
                "phase": phase,
                "phase_priority": phase_priority,
                "source_group": directory.parent.name,
                "job_directory": directory.name,
                "request_name": details["request_name"],
                "candidate_id": candidate_id,
                "server_seed": details["server_seed"],
                "protein_chain_count": details["protein_chain_count"],
                "requested_peptide": details["requested_peptide"],
                "model_cif_count": len(list(directory.glob("*_model_*.cif"))),
                "confidence_json_count": len(list(directory.glob("*_summary_confidences_*.json"))),
                "full_data_json_count": len(list(directory.glob("*_full_data_*.json"))),
                "cohort": cohort,
                "source_path": str(directory),
            })

    # The same saved AF3 job exists in several copied folders. Retain exactly one
    # canonical copy, preferring the original Monday download for identical jobs.
    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in discovered:
        groups[(str(row["candidate_id"]), str(row["requested_peptide"]), str(row["server_seed"]))].append(row)
    for group in groups.values():
        group.sort(key=lambda row: (int(row["phase_priority"]), str(row["source_path"])))
        canonical = group[0]
        canonical_key = f"{canonical['candidate_id']}|{canonical['server_seed']}|{canonical['requested_peptide']}"
        for index, row in enumerate(group):
            row["canonical_job_key"] = canonical_key
            row["canonical_for_analysis"] = index == 0
            row["duplicate_of_source_path"] = "" if index == 0 else canonical["source_path"]
            row["integrity_status"] = (
                "complete_5_of_5"
                if row["model_cif_count"] == row["confidence_json_count"] == row["full_data_json_count"] == 5
                else "incomplete"
            )
    return sorted(discovered, key=lambda row: (str(row["phase"]), str(row["source_group"]), str(row["job_directory"])))


def canonical_sample_rows(inventory: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    metadata = build_candidate_metadata()
    samples: list[dict[str, Any]] = []
    jobs: list[dict[str, Any]] = []
    for row in inventory:
        if not row["canonical_for_analysis"] or row["integrity_status"] != "complete_5_of_5":
            continue
        if row["candidate_id"] not in metadata:
            continue
        directory = Path(str(row["source_path"]))
        result = analyze_complete_job(directory, str(row["phase"]), metadata)
        job = dict(result["job"])
        job.update({
            "canonical_job_key": row["canonical_job_key"],
            "source_path": row["source_path"],
        })
        jobs.append(job)
        for sample in result["samples"]:
            sample = dict(sample)
            index = int(sample["sample_index"])
            matches = list(directory.glob(f"*_model_{index}.cif"))
            if len(matches) != 1:
                raise FileNotFoundError(f"Expected one model {index} in {directory}")
            sample.update({
                "canonical_job_key": row["canonical_job_key"],
                "model_path": str(matches[0]),
                "source_path": row["source_path"],
            })
            samples.append(sample)
    return samples, jobs


def phase_candidate_ids(inventory: list[dict[str, Any]], phase: str) -> set[str]:
    return {
        str(row["candidate_id"])
        for row in inventory
        if row["phase"] == phase
        and row["cohort"] == "legacy_candidate_pmhc"
        and row["integrity_status"] == "complete_5_of_5"
    }


def structural_consistency_class(lt2_fraction: float, ge5_fraction: float, median_rmsd: float) -> str:
    """Label technical pose consistency without converting it to biology."""
    if lt2_fraction >= 0.8 and ge5_fraction <= 0.2 and median_rmsd < 2.0:
        return "tier_A_robust"
    if lt2_fraction >= 0.5 and median_rmsd < 2.0:
        return "tier_B_mixed"
    if lt2_fraction > 0:
        return "tier_C_unstable_partial_pose"
    return "tier_D_no_consistent_pose"


def summarize_controlled_comparison(
    target_geometry: list[dict[str, Any]],
    background_geometry: list[dict[str, Any]],
) -> dict[str, Any]:
    """Summarize the frozen target-versus-background geometry descriptively."""
    if not background_geometry:
        return {
            "controlled_comparison_status": "not_evaluable_no_completed_matched_background",
            "structural_background_comparator_count": 0,
            "structural_background_geometry_count": 0,
            "target_candidate_exposed_rmsd_median_A": "",
            "background_candidate_exposed_rmsd_median_A": "",
            "background_minus_target_exposed_rmsd_median_A": "",
            "controlled_comparison_p_value": "",
        }
    target_values = [float(row["candidate_exposed_ca_rmsd_A"]) for row in target_geometry]
    values_by_candidate: dict[str, list[float]] = defaultdict(list)
    for row in background_geometry:
        values_by_candidate[str(row["background_candidate_id"])].append(
            float(row["candidate_exposed_ca_rmsd_A"])
        )
    background_candidate_medians = [
        median(values) for _, values in sorted(values_by_candidate.items())
    ]
    target_median = median(target_values)
    background_median = median(background_candidate_medians)
    return {
        "controlled_comparison_status": "evaluable_descriptive_matched_background",
        "structural_background_comparator_count": len(values_by_candidate),
        "structural_background_geometry_count": len(background_geometry),
        "target_candidate_exposed_rmsd_median_A": round(target_median, 6),
        "background_candidate_exposed_rmsd_median_A": round(background_median, 6),
        "background_minus_target_exposed_rmsd_median_A": round(
            background_median - target_median, 6
        ),
        "controlled_comparison_p_value": "",
    }


def build_pair_score_sheet(
    inventory: list[dict[str, Any]], samples: list[dict[str, Any]], jobs: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    shortlist = read_csv(PROC / "fullscreen_tier1_ebv_myelin_shortlist.csv")
    universe = {
        row["pair_id"]: row
        for row in read_csv(PROC / "register_aware_benchmark" / "benchmark_pair_universe.csv")
    }
    backgrounds = read_csv(PROC / "expanded_background" / "background_register_prediction_summary.csv")
    completed_background_ids = {
        str(row["candidate_id"])
        for row in jobs
        if row["af3_cohort"] == "new_human_background_pmhc"
        and row["sequence_layout_status"] == "pass_exact_three_chain_peptide_match"
    }
    monday_ids = phase_candidate_ids(inventory, "monday_af3_full")
    focused_ids = phase_candidate_ids(inventory, "focused_af3_rerun")
    samples_by_candidate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    background_samples_by_candidate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in samples:
        if row["af3_cohort"] == "legacy_candidate_pmhc" and row["sequence_layout_status"] == "pass_exact_three_chain_peptide_match":
            samples_by_candidate[str(row["candidate_id"])].append(row)
        elif row["af3_cohort"] == "new_human_background_pmhc" and row["sequence_layout_status"] == "pass_exact_three_chain_peptide_match":
            background_samples_by_candidate[str(row["candidate_id"])].append(row)

    model_cache: dict[str, dict[str, list[dict[str, object]]]] = {}
    score_rows: list[dict[str, Any]] = []
    geometry_rows: list[dict[str, Any]] = []
    controlled_geometry_rows: list[dict[str, Any]] = []
    for shortlist_rank, prior in enumerate(shortlist, start=1):
        pair_id = f"{prior['ebv_candidate_id']}::{prior['human_candidate_id']}"
        pair = universe[pair_id]
        ebv_id, human_id = pair["ebv_candidate_id"], pair["human_candidate_id"]
        register_eligible = pair["ebv_register_status"] in VALID_REGISTER_STATUSES and pair["human_register_status"] in VALID_REGISTER_STATUSES
        monday_both = ebv_id in monday_ids and human_id in monday_ids
        focused_both = ebv_id in focused_ids and human_id in focused_ids
        available = ebv_id in samples_by_candidate and human_id in samples_by_candidate
        sequence_metrics: dict[str, Any] = {}
        if register_eligible:
            sequence_metrics = direct_register_sequence_metrics(pair["ebv_top_core_peptide"], pair["human_top_core_peptide"])
        pair_geometry: list[dict[str, float]] = []
        if register_eligible and available:
            for ebv_sample in samples_by_candidate[ebv_id]:
                ebv_path = str(ebv_sample["model_path"])
                if ebv_path not in model_cache:
                    model_cache[ebv_path] = parse_mmcif(Path(ebv_path))
                for human_sample in samples_by_candidate[human_id]:
                    human_path = str(human_sample["model_path"])
                    if human_path not in model_cache:
                        model_cache[human_path] = parse_mmcif(Path(human_path))
                    metrics = same_register_geometry(
                        model_cache[ebv_path],
                        model_cache[human_path],
                        int(pair["ebv_top_core_start_1_based"]),
                        int(pair["human_top_core_start_1_based"]),
                    )
                    pair_geometry.append(metrics)
                    geometry_rows.append({
                        "pair_id": pair_id,
                        "shortlist_rank": shortlist_rank,
                        "ebv_candidate_id": ebv_id,
                        "human_candidate_id": human_id,
                        "ebv_job_key": ebv_sample["canonical_job_key"],
                        "ebv_sample_index": ebv_sample["sample_index"],
                        "human_job_key": human_sample["canonical_job_key"],
                        "human_sample_index": human_sample["sample_index"],
                        **{key: round(value, 6) for key, value in metrics.items()},
                    })

        human_bin = binding_rank_bin(float(pair["human_binding_rank"]))
        controls = matched_background_feasibility(len(pair["human_peptide"]), human_bin, backgrounds, completed_background_ids)
        matched_background_geometry: list[dict[str, Any]] = []
        completed_control_ids = [
            identifier
            for identifier in str(controls["completed_matched_background_ids"]).split(";")
            if identifier
        ]
        if pair_geometry:
            for background_id in completed_control_ids:
                background = next(
                    row for row in backgrounds if row["candidate_id"] == background_id
                )
                core_starts = str(background["predicted_core_start_positions_1_based"]).split(";")
                if (
                    background.get("predicted_core_fully_contained_in_manifest_peptide") != "True"
                    or len(core_starts) != 1
                    or not core_starts[0].isdigit()
                ):
                    continue
                background_start = int(core_starts[0])
                for ebv_sample in samples_by_candidate[ebv_id]:
                    ebv_path = str(ebv_sample["model_path"])
                    if ebv_path not in model_cache:
                        model_cache[ebv_path] = parse_mmcif(Path(ebv_path))
                    for background_sample in background_samples_by_candidate[background_id]:
                        background_path = str(background_sample["model_path"])
                        if background_path not in model_cache:
                            model_cache[background_path] = parse_mmcif(Path(background_path))
                        metrics = same_register_geometry(
                            model_cache[ebv_path],
                            model_cache[background_path],
                            int(pair["ebv_top_core_start_1_based"]),
                            background_start,
                        )
                        comparison = {
                            "pair_id": pair_id,
                            "shortlist_rank": shortlist_rank,
                            "ebv_candidate_id": ebv_id,
                            "target_human_candidate_id": human_id,
                            "background_candidate_id": background_id,
                            "background_predicted_core": background["predicted_core_peptide"],
                            "ebv_job_key": ebv_sample["canonical_job_key"],
                            "ebv_sample_index": ebv_sample["sample_index"],
                            "background_job_key": background_sample["canonical_job_key"],
                            "background_sample_index": background_sample["sample_index"],
                            **{key: round(value, 6) for key, value in metrics.items()},
                            "interpretation": "Frozen score-blind human-background pMHC comparison; technical samples are not biological replicates.",
                        }
                        matched_background_geometry.append(comparison)
                        controlled_geometry_rows.append(comparison)
        controlled_summary = summarize_controlled_comparison(
            pair_geometry, matched_background_geometry
        )
        row: dict[str, Any] = {
            "discovery_priority_rank": "",
            "pair_id": pair_id,
            "original_colabfold_shortlist_rank": shortlist_rank,
            "ebv_candidate_id": ebv_id,
            "human_candidate_id": human_id,
            "ebv_peptide": pair["ebv_peptide"],
            "human_peptide": pair["human_peptide"],
            "ebv_p1_p9_core": pair["ebv_top_core_peptide"],
            "human_p1_p9_core": pair["human_top_core_peptide"],
            "ebv_register_status": pair["ebv_register_status"],
            "human_register_status": pair["human_register_status"],
            "register_eligible_primary_allele": register_eligible,
            "monday_af3_both_models_available": monday_both,
            "focused_rerun_both_models_available": focused_both,
            "combined_af3_both_models_available": available,
            "ebv_unique_af3_job_count": len({row["canonical_job_key"] for row in samples_by_candidate.get(ebv_id, [])}),
            "human_unique_af3_job_count": len({row["canonical_job_key"] for row in samples_by_candidate.get(human_id, [])}),
            "ebv_af3_sample_count": len(samples_by_candidate.get(ebv_id, [])),
            "human_af3_sample_count": len(samples_by_candidate.get(human_id, [])),
            "af3_cross_sample_geometry_count": len(pair_geometry),
            "colabfold_local_property_similarity": prior["property_similarity"],
            "colabfold_local_peptide_ca_rmsd_A": prior["local_peptide_ca_rmsd_after_hla_fit"],
            "colabfold_review_priority_heuristic": prior["review_priority_heuristic"],
            "same_register_property_similarity": sequence_metrics.get("same_register_property_similarity", ""),
            "candidate_exposed_identity_fraction": sequence_metrics.get("candidate_exposed_identity_fraction", ""),
            "candidate_exposed_property_similarity": sequence_metrics.get("candidate_exposed_property_similarity", ""),
            "core_p1_p9_ca_rmsd_A_median": "",
            "core_p1_p9_ca_rmsd_A_min": "",
            "core_p1_p9_ca_rmsd_A_max": "",
            "candidate_exposed_ca_rmsd_A_median": "",
            "candidate_exposed_ca_rmsd_A_min": "",
            "candidate_exposed_ca_rmsd_A_max": "",
            "candidate_exposed_rmsd_lt_1A_fraction": "",
            "candidate_exposed_rmsd_lt_2A_fraction": "",
            "candidate_exposed_rmsd_ge_5A_fraction": "",
            "structural_consistency_tier": "",
            "target_human_binding_rank_bin": human_bin,
            **controls,
            **controlled_summary,
            "audit_status": "",
            "claim_boundary": "Computational pMHC prioritization only; not evidence of shared-TCR recognition or molecular mimicry.",
        }
        if pair_geometry:
            row.update(_metric_distribution(pair_geometry, "core_p1_p9_ca_rmsd_A"))
            row.update(_metric_distribution(pair_geometry, "candidate_exposed_ca_rmsd_A"))
            exposed = [float(item["candidate_exposed_ca_rmsd_A"]) for item in pair_geometry]
            row.update({
                "candidate_exposed_rmsd_lt_1A_fraction": round(sum(value < 1.0 for value in exposed) / len(exposed), 6),
                "candidate_exposed_rmsd_lt_2A_fraction": round(sum(value < 2.0 for value in exposed) / len(exposed), 6),
                "candidate_exposed_rmsd_ge_5A_fraction": round(sum(value >= 5.0 for value in exposed) / len(exposed), 6),
                "audit_status": "eligible_and_structurally_evaluated",
            })
            row["structural_consistency_tier"] = structural_consistency_class(
                float(row["candidate_exposed_rmsd_lt_2A_fraction"]),
                float(row["candidate_exposed_rmsd_ge_5A_fraction"]),
                float(row["candidate_exposed_ca_rmsd_A_median"]),
            )
        elif not register_eligible:
            row["audit_status"] = "not_discovery_rankable_register_or_allele_unresolved"
        else:
            row["audit_status"] = "eligible_but_missing_complete_af3_partner"
        score_rows.append(row)

    # This is a transparent prioritization order, not a synthetic biological score:
    # robustness first, then lower median exposed RMSD, then exposed chemistry.
    rankable = [row for row in score_rows if row["audit_status"] == "eligible_and_structurally_evaluated"]
    rankable.sort(key=lambda row: (
        -float(row["candidate_exposed_rmsd_lt_2A_fraction"]),
        float(row["candidate_exposed_ca_rmsd_A_median"]),
        -float(row["candidate_exposed_property_similarity"]),
        int(row["original_colabfold_shortlist_rank"]),
    ))
    for rank, row in enumerate(rankable, start=1):
        row["discovery_priority_rank"] = rank
    score_rows.sort(key=lambda row: (row["discovery_priority_rank"] == "", int(row["discovery_priority_rank"] or 999), int(row["original_colabfold_shortlist_rank"])))
    return score_rows, geometry_rows, controlled_geometry_rows


def render_findings(inventory: list[dict[str, Any]], jobs: list[dict[str, Any]], scores: list[dict[str, Any]]) -> str:
    phase_counts = Counter(str(row["phase"]) for row in inventory)
    canonical_counts = Counter(str(row["phase"]) for row in inventory if row["canonical_for_analysis"])
    eligible = [row for row in scores if row["register_eligible_primary_allele"]]
    monday_unlocked = [
        row for row in scores
        if row["register_eligible_primary_allele"]
        and row["monday_af3_both_models_available"]
        and not row["focused_rerun_both_models_available"]
    ]
    top = [row for row in scores if row["discovery_priority_rank"] != ""][:10]
    controlled = [
        row for row in eligible
        if row["controlled_comparison_status"] == "evaluable_descriptive_matched_background"
    ]
    lines = [
        "# Complete pMHC modeling audit",
        "",
        "## Direct answer",
        "",
        "The Monday AF3 pMHC models were preserved, but the focused same-register analysis did not read their folder. It read the later `json1folds/json2folds` collection instead. This omitted seven otherwise eligible frozen-shortlist pairs from the nine-pair result table. The combined audit now includes the Monday models and removes copied duplicate jobs before scoring.",
        "",
        "## Reconciled model inventory",
        "",
        f"- Saved job directories discovered: **{len(inventory)}** ({', '.join(f'{key}: {value}' for key, value in sorted(phase_counts.items()))}).",
        f"- Unique candidate/peptide/seed jobs after deduplication: **{sum(canonical_counts.values())}** ({', '.join(f'{key}: {value}' for key, value in sorted(canonical_counts.items()))}).",
        f"- Canonical study jobs with complete parsed pMHC outputs: **{len(jobs)}**.",
        f"- Frozen shortlist pairs: **{len(scores)}**; primary-allele register-eligible: **{len(eligible)}**; structurally evaluated after the audit: **{sum(row['audit_status'] == 'eligible_and_structurally_evaluated' for row in scores)}**.",
        "",
        "## What Monday changes",
        "",
        f"Monday unlocks **{len(monday_unlocked)}** eligible pairs that were absent from the focused rerun table:",
        "",
    ]
    lines.extend(f"- Rank {row['original_colabfold_shortlist_rank']}: `{row['pair_id']}`" for row in monday_unlocked)
    lines.extend([
        "",
        "## Best current pairings",
        "",
        "The audit rank is lexicographic: fraction of AF3 comparisons below 2 A (higher first), median exposed-position RMSD (lower first), then exposed-position property similarity (higher first). It is a transparent prioritization order, not a biological probability.",
        "",
        "| Audit rank | Pair | Robustness tier | <2 A fraction | >=5 A fraction | Median exposed RMSD (A) | Exposed property similarity |",
        "|---:|---|---|---:|---:|---:|---:|",
    ])
    for row in top:
        lines.append(
            f"| {row['discovery_priority_rank']} | `{row['pair_id']}` | {row['structural_consistency_tier']} | {float(row['candidate_exposed_rmsd_lt_2A_fraction']):.2f} | {float(row['candidate_exposed_rmsd_ge_5A_fraction']):.2f} | {float(row['candidate_exposed_ca_rmsd_A_median']):.3f} | {float(row['candidate_exposed_property_similarity']):.3f} |"
        )
    lines.extend([
        "",
        "## Claim boundary and remaining limitation",
        "",
        f"Length/bin-matched structural background comparisons are now available for **{len(controlled)}** of the **{len(eligible)}** primary-allele register-eligible shortlist pairs. Background summaries weight each unique comparator candidate once; repeated server jobs and model samples remain technical sensitivity analyses rather than biological replicates. These are descriptive controls with no p-value. The score sheet therefore prioritizes candidates for a computational paper and future testing; it does not establish presentation, shared-TCR recognition, cross-reactivity, molecular mimicry, or an MS mechanism.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inventory = discover_jobs()
    samples, jobs = canonical_sample_rows(inventory)
    scores, geometry, controlled_geometry = build_pair_score_sheet(inventory, samples, jobs)
    write_csv(OUT / "complete_af3_job_inventory.csv", inventory)
    write_csv(OUT / "canonical_af3_job_summary.csv", jobs)
    write_csv(OUT / "canonical_af3_sample_metrics.csv", samples)
    write_csv(OUT / "master_pair_score_sheet.csv", scores)
    write_csv(OUT / "combined_same_register_geometry.csv", geometry)
    write_csv(OUT / "matched_background_structure_geometry.csv", controlled_geometry)
    (OUT / "AUDIT_FINDINGS.md").write_text(render_findings(inventory, jobs, scores), encoding="utf-8")
    print(
        f"Audited {len(inventory)} saved jobs, retained {sum(bool(row['canonical_for_analysis']) for row in inventory)} unique jobs, "
        f"and evaluated {sum(row['audit_status'] == 'eligible_and_structurally_evaluated' for row in scores)} of {len(scores)} shortlist pairs."
    )


if __name__ == "__main__":
    main()
