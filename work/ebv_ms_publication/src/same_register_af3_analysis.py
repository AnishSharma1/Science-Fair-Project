"""Direct P1--P9 analysis of predeclared AF3 EBV--myelin pMHC pairs.

The analysis compares equivalent HLA-II register positions without requiring a
local sequence alignment. It is restricted to the frozen 32-pair shortlist and
pre-score human-background comparator universe. Outputs are computational pMHC
descriptors, not evidence of presentation, TCR binding, cross-reactivity, or MS
mechanism.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

import numpy as np

from analyze_af3_pmhc_downloads import AF3_ROOT, ca_coordinates, kabsch, parse_mmcif
from premeeting_rigor import binding_rank_bin
from register_aware_scoring import ANCHOR_POSITIONS, CANDIDATE_EXPOSED_POSITIONS, property_similarity


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "same_register_af3_analysis"
VALID_REGISTER_STATUSES = {"iedb_top_core_hypothesis", "experimental_primary_allele_reference"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _position_summary(left: str, right: str, positions: list[int]) -> dict[str, object]:
    properties = [property_similarity(left[position - 1], right[position - 1]) for position in positions]
    components = {
        key: round(mean(item[key] for item in properties), 6)
        for key in properties[0]
    }
    return {
        "positions": ";".join(f"P{position}" for position in positions),
        "count": len(positions),
        "identity_count": sum(left[position - 1] == right[position - 1] for position in positions),
        "identity_fraction": round(sum(left[position - 1] == right[position - 1] for position in positions) / len(positions), 6),
        **components,
        "property_similarity": round(mean(components.values()), 6),
    }


def direct_register_sequence_metrics(ebv_core: str, human_core: str) -> dict[str, object]:
    """Compare every equivalent position in two exact nine-residue cores."""
    if len(ebv_core) != 9 or len(human_core) != 9:
        raise ValueError("Direct register comparison requires two exact 9-mer cores")
    all_positions = list(range(1, 10))
    anchors = sorted(ANCHOR_POSITIONS)
    exposed = sorted(CANDIDATE_EXPOSED_POSITIONS)
    all_summary = _position_summary(ebv_core, human_core, all_positions)
    anchor_summary = _position_summary(ebv_core, human_core, anchors)
    exposed_summary = _position_summary(ebv_core, human_core, exposed)
    return {
        "same_register_positions": all_summary["positions"],
        "same_register_position_count": all_summary["count"],
        "same_register_identity_count": all_summary["identity_count"],
        "same_register_identity_fraction": all_summary["identity_fraction"],
        "same_register_property_similarity": all_summary["property_similarity"],
        "anchor_positions": anchor_summary["positions"],
        "anchor_identity_count": anchor_summary["identity_count"],
        "anchor_identity_fraction": anchor_summary["identity_fraction"],
        "anchor_property_similarity": anchor_summary["property_similarity"],
        "candidate_exposed_positions": exposed_summary["positions"],
        "candidate_exposed_identity_count": exposed_summary["identity_count"],
        "candidate_exposed_identity_fraction": exposed_summary["identity_fraction"],
        "candidate_exposed_property_similarity": exposed_summary["property_similarity"],
    }


def _rmsd(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape:
        raise ValueError(f"Coordinate shape mismatch: {left.shape} versus {right.shape}")
    return float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def same_register_geometry(
    ebv_model: dict[str, list[dict[str, object]]],
    human_model: dict[str, list[dict[str, object]]],
    ebv_core_start_1_based: int,
    human_core_start_1_based: int,
) -> dict[str, float]:
    """Fit the HLA grooves, then compare C-alpha positions at matching P1--P9."""
    for model in (ebv_model, human_model):
        if set(model) != {"A", "B", "C"}:
            raise ValueError("Expected exact A/B/C pMHC chain layout")
    ebv_hla = np.vstack([ca_coordinates(ebv_model[chain][:85]) for chain in ("A", "B")])
    human_hla = np.vstack([ca_coordinates(human_model[chain][:85]) for chain in ("A", "B")])
    rotation, translation = kabsch(human_hla, ebv_hla)
    fitted_hla = human_hla @ rotation + translation
    ebv_start, human_start = ebv_core_start_1_based - 1, human_core_start_1_based - 1
    ebv_core = ca_coordinates(ebv_model["C"][ebv_start:ebv_start + 9])
    human_core = ca_coordinates(human_model["C"][human_start:human_start + 9]) @ rotation + translation
    if ebv_core.shape != (9, 3) or human_core.shape != (9, 3):
        raise ValueError("Resolved register is not fully contained in model peptide")
    anchor_indices = np.asarray([position - 1 for position in sorted(ANCHOR_POSITIONS)])
    exposed_indices = np.asarray([position - 1 for position in sorted(CANDIDATE_EXPOSED_POSITIONS)])
    return {
        "hla_groove_ca_rmsd_A": _rmsd(ebv_hla, fitted_hla),
        "core_p1_p9_ca_rmsd_A": _rmsd(ebv_core, human_core),
        "anchor_ca_rmsd_A": _rmsd(ebv_core[anchor_indices], human_core[anchor_indices]),
        "candidate_exposed_ca_rmsd_A": _rmsd(ebv_core[exposed_indices], human_core[exposed_indices]),
    }


def matched_background_feasibility(
    target_human_length: int,
    target_human_binding_bin: str,
    backgrounds: list[dict[str, str]],
    completed_background_ids: set[str],
) -> dict[str, object]:
    """Count frozen human backgrounds meeting the predeclared length/bin rule."""
    planned = sorted(
        row["candidate_id"] for row in backgrounds
        if abs(len(row["peptide"]) - target_human_length) <= 1
        and row["binding_rank_bin"] == target_human_binding_bin
    )
    completed = [identifier for identifier in planned if identifier in completed_background_ids]
    return {
        "planned_matched_background_count": len(planned),
        "planned_matched_background_ids": ";".join(planned),
        "completed_matched_background_count": len(completed),
        "completed_matched_background_ids": ";".join(completed),
    }


def evaluate_sequence_decoys(
    target_pair_id: str,
    target_metrics: dict[str, object],
    decoy_metrics: list[dict[str, object]],
) -> dict[str, object]:
    """Compare exposed-position chemistry to frozen decoys, descriptively only."""
    target_score = float(target_metrics["candidate_exposed_property_similarity"])
    if not decoy_metrics:
        return {
            "sequence_decoy_evaluation_status": "not_evaluable_no_matched_sequence_decoys",
            "sequence_decoy_count": 0,
            "sequence_decoy_ids": "",
            "sequence_decoy_exposed_property_scores": "",
            "target_exposed_property_similarity": target_score,
            "decoy_exposed_property_similarity_median": "",
            "target_minus_decoy_median": "",
            "target_rank_among_target_plus_decoys": "",
            "p_value": "",
        }
    ordered_decoys = sorted(decoy_metrics, key=lambda row: str(row["candidate_id"]))
    decoy_scores = [float(row["candidate_exposed_property_similarity"]) for row in ordered_decoys]
    ranked = sorted(
        [(target_score, target_pair_id), *[(float(row["candidate_exposed_property_similarity"]), str(row["candidate_id"])) for row in ordered_decoys]],
        key=lambda item: (-item[0], item[1]),
    )
    target_rank = next(index for index, (_, identifier) in enumerate(ranked, start=1) if identifier == target_pair_id)
    return {
        "sequence_decoy_evaluation_status": "evaluable_descriptive_sequence_only",
        "sequence_decoy_count": len(ordered_decoys),
        "sequence_decoy_ids": ";".join(str(row["candidate_id"]) for row in ordered_decoys),
        "sequence_decoy_exposed_property_scores": ";".join(f"{float(row['candidate_exposed_property_similarity']):.6f}" for row in ordered_decoys),
        "target_exposed_property_similarity": target_score,
        "decoy_exposed_property_similarity_median": round(median(decoy_scores), 6),
        "target_minus_decoy_median": round(target_score - median(decoy_scores), 6),
        "target_rank_among_target_plus_decoys": target_rank,
        "p_value": "",
    }


def _model_path(sample: dict[str, str]) -> Path:
    job = AF3_ROOT / sample["source_folder"] / sample["job_directory"]
    matches = list(job.glob(f"*_model_{int(sample['sample_index'])}.cif"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one model for {sample['candidate_id']} sample {sample['sample_index']}")
    return matches[0]


def _metric_distribution(rows: list[dict[str, float]], field: str) -> dict[str, float]:
    values = [row[field] for row in rows]
    return {
        f"{field}_median": round(median(values), 3),
        f"{field}_min": round(min(values), 3),
        f"{field}_max": round(max(values), 3),
    }


def summarize_job_pair_geometry(geometry_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Expose seed/job sensitivity without treating model samples as replicates."""
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in geometry_rows:
        groups[(str(row["pair_id"]), str(row["ebv_job_directory"]), str(row["human_job_directory"]))].append(row)
    output = []
    for (pair_id, ebv_job, human_job), rows in sorted(groups.items()):
        core = [float(row["core_p1_p9_ca_rmsd_A"]) for row in rows]
        exposed = [float(row["candidate_exposed_ca_rmsd_A"]) for row in rows]
        output.append({
            "pair_id": pair_id,
            "shortlist_rank": rows[0]["shortlist_rank"],
            "ebv_job_directory": ebv_job,
            "human_job_directory": human_job,
            "cross_sample_comparison_count": len(rows),
            "core_p1_p9_ca_rmsd_A_median": round(median(core), 3),
            "core_p1_p9_ca_rmsd_A_min": round(min(core), 3),
            "core_p1_p9_ca_rmsd_A_max": round(max(core), 3),
            "candidate_exposed_ca_rmsd_A_median": round(median(exposed), 3),
            "candidate_exposed_ca_rmsd_A_min": round(min(exposed), 3),
            "candidate_exposed_ca_rmsd_A_max": round(max(exposed), 3),
            "interpretation": "Seed/job sensitivity only; model samples are not independent biological replicates.",
        })
    return output


def build_analysis() -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Analyze only shortlist pairs with AF3 samples and resolved primary cores."""
    shortlist = read_csv(PROC / "fullscreen_tier1_ebv_myelin_shortlist.csv")
    pair_universe = {
        row["pair_id"]: row
        for row in read_csv(PROC / "register_aware_benchmark" / "benchmark_pair_universe.csv")
    }
    samples = [
        row for row in read_csv(PROC / "alphafold_server_pmhc_descriptive_analysis" / "af3_pmhc_sample_metrics.csv")
        if row["af3_cohort"] == "legacy_candidate_pmhc"
        and row["sequence_layout_status"] == "pass_exact_three_chain_peptide_match"
    ]
    samples_by_candidate: dict[str, list[dict[str, str]]] = {}
    for sample in samples:
        samples_by_candidate.setdefault(sample["candidate_id"], []).append(sample)
    completed_background_ids = {
        row["candidate_id"]
        for row in read_csv(PROC / "alphafold_server_pmhc_descriptive_analysis" / "af3_pmhc_job_summary.csv")
        if row["af3_cohort"] == "new_human_background_pmhc"
    }
    backgrounds = read_csv(PROC / "expanded_background" / "background_register_prediction_summary.csv")
    background_by_id = {row["candidate_id"]: row for row in backgrounds}
    model_cache: dict[Path, dict[str, list[dict[str, object]]]] = {}
    pair_rows, geometry_rows, sequence_decoy_rows = [], [], []
    for shortlist_rank, shortlist_row in enumerate(shortlist, start=1):
        pair_id = f"{shortlist_row['ebv_candidate_id']}::{shortlist_row['human_candidate_id']}"
        pair = pair_universe[pair_id]
        if pair["ebv_candidate_id"] not in samples_by_candidate or pair["human_candidate_id"] not in samples_by_candidate:
            continue
        if pair["ebv_register_status"] not in VALID_REGISTER_STATUSES or pair["human_register_status"] not in VALID_REGISTER_STATUSES:
            continue
        ebv_start, human_start = int(pair["ebv_top_core_start_1_based"]), int(pair["human_top_core_start_1_based"])
        ebv_core, human_core = pair["ebv_top_core_peptide"], pair["human_top_core_peptide"]
        sequence_metrics = direct_register_sequence_metrics(ebv_core, human_core)
        pair_geometry = []
        for ebv_sample in samples_by_candidate[pair["ebv_candidate_id"]]:
            ebv_path = _model_path(ebv_sample)
            if ebv_path not in model_cache:
                model_cache[ebv_path] = parse_mmcif(ebv_path)
            for human_sample in samples_by_candidate[pair["human_candidate_id"]]:
                human_path = _model_path(human_sample)
                if human_path not in model_cache:
                    model_cache[human_path] = parse_mmcif(human_path)
                metrics = same_register_geometry(model_cache[ebv_path], model_cache[human_path], ebv_start, human_start)
                combination = {
                    "pair_id": pair_id,
                    "shortlist_rank": shortlist_rank,
                    "ebv_job_directory": ebv_sample["job_directory"],
                    "ebv_sample_index": ebv_sample["sample_index"],
                    "human_job_directory": human_sample["job_directory"],
                    "human_sample_index": human_sample["sample_index"],
                    **{key: round(value, 6) for key, value in metrics.items()},
                    "interpretation": "Within-pair AF3 sample sensitivity only; samples are not independent biological replicates.",
                }
                geometry_rows.append(combination)
                pair_geometry.append(metrics)
        human_bin = binding_rank_bin(float(pair["human_binding_rank"]))
        feasibility = matched_background_feasibility(len(pair["human_peptide"]), human_bin, backgrounds, completed_background_ids)
        planned_ids = [identifier for identifier in str(feasibility["planned_matched_background_ids"]).split(";") if identifier]
        decoy_metrics = []
        for identifier in planned_ids:
            background = background_by_id[identifier]
            core = background["predicted_core_peptide"]
            if len(core) != 9 or background.get("predicted_core_fully_contained_in_manifest_peptide") != "True":
                continue
            decoy_metrics.append({
                "candidate_id": identifier,
                **direct_register_sequence_metrics(ebv_core, core),
            })
        sequence_evaluation = evaluate_sequence_decoys(pair_id, sequence_metrics, decoy_metrics)
        sequence_decoy_rows.append({
            "pair_id": pair_id,
            "shortlist_rank": shortlist_rank,
            "ebv_p1_p9_core": ebv_core,
            "target_human_p1_p9_core": human_core,
            **sequence_evaluation,
            "statistical_status": "descriptive_only_no_p_value_small_nonindependent_control_set",
            "interpretation": "Frozen length/bin-matched sequence chemistry control only; no structural comparator or biological inference.",
        })
        pair_rows.append({
            "pair_id": pair_id,
            "shortlist_rank": shortlist_rank,
            "ebv_candidate_id": pair["ebv_candidate_id"],
            "human_candidate_id": pair["human_candidate_id"],
            "ebv_peptide": pair["ebv_peptide"],
            "human_peptide": pair["human_peptide"],
            "ebv_core_start_1_based": ebv_start,
            "human_core_start_1_based": human_start,
            "ebv_p1_p9_core": ebv_core,
            "human_p1_p9_core": human_core,
            "register_evidence_level": "computational_unique_IEDB_top_core_hypotheses",
            "ebv_af3_sample_count": len(samples_by_candidate[pair["ebv_candidate_id"]]),
            "human_af3_sample_count": len(samples_by_candidate[pair["human_candidate_id"]]),
            "af3_cross_sample_geometry_count": len(pair_geometry),
            **sequence_metrics,
            **_metric_distribution(pair_geometry, "hla_groove_ca_rmsd_A"),
            **_metric_distribution(pair_geometry, "core_p1_p9_ca_rmsd_A"),
            **_metric_distribution(pair_geometry, "anchor_ca_rmsd_A"),
            **_metric_distribution(pair_geometry, "candidate_exposed_ca_rmsd_A"),
            "target_human_binding_rank_bin": human_bin,
            **feasibility,
            **sequence_evaluation,
            "controlled_comparison_status": "not_evaluable_no_completed_matched_background",
            "interpretation": "Predeclared-pair direct P1-P9 pMHC comparison only; no TCR, cross-reactivity, or MS-mechanism inference.",
        })
    return pair_rows, geometry_rows, sequence_decoy_rows


def render_readme(pair_rows: list[dict[str, object]], geometry_rows: list[dict[str, object]]) -> str:
    return f"""# Direct same-register AF3 analysis

This analysis compares the frozen 32-pair shortlist only. A pair enters the
table when both pMHCs have complete AF3 samples and each peptide has a unique,
fully contained primary-allele computational P1--P9 core hypothesis.

Unlike the earlier local-alignment diagnostic, this endpoint compares all nine
equivalent register positions directly. HLA-DRA/DRB groove C-alpha atoms are
fitted first; P1--P9, anchor-position (P1/P4/P6/P9), and candidate-exposed-
position (P2/P3/P5/P7/P8) C-alpha RMSDs are then reported across every
available AF3 sample combination.

- Eligible predeclared pairs: **{len(pair_rows)}**
- AF3 cross-sample geometry comparisons: **{len(geometry_rows)}**
- Pairs with at least one completed matched background comparator: **{sum(int(row['completed_matched_background_count']) > 0 for row in pair_rows)}**

The background rule was frozen before structural inspection: human peptide
length within one residue and the same IEDB predicted-binding-rank bin. No
completed matched comparator exists for any eligible pair, so no controlled
structural effect size, enrichment test, or p-value is computed. Where frozen
background cores exist despite failed structure generation, the sequence-only
table reports a descriptive same-register chemistry rank; these small control
sets do not support inferential statistics.

AlphaFold samples and seed repeats are technical sensitivity analyses, not
independent biological replicates. These files do not establish peptide
presentation, a shared TCR surface, T-cell activation, cross-reactivity,
molecular mimicry, or an MS mechanism.
"""


def main() -> None:
    pair_rows, geometry_rows, sequence_decoy_rows = build_analysis()
    if not pair_rows or not geometry_rows:
        raise ValueError("No eligible same-register AF3 pairs")
    write_csv(OUT / "same_register_af3_pair_results.csv", pair_rows)
    write_csv(OUT / "same_register_af3_geometry_ensemble.csv", geometry_rows)
    write_csv(OUT / "same_register_af3_job_pair_sensitivity.csv", summarize_job_pair_geometry(geometry_rows))
    write_csv(OUT / "same_register_sequence_decoy_comparison.csv", sequence_decoy_rows)
    (OUT / "README.md").write_text(render_readme(pair_rows, geometry_rows), encoding="utf-8")
    print(f"Wrote {len(pair_rows)} pair rows and {len(geometry_rows)} geometry rows to {OUT}")


if __name__ == "__main__":
    main()
