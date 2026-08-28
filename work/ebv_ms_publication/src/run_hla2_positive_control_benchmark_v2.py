"""Analyze fresh AlphaFold results for the frozen HLA-II benchmark v2 pilot.

The command is deliberately safe to run before submission: missing jobs produce a
complete, machine-readable not-evaluable package. Pilot results can never freeze
weights or unlock discovery rankings.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
from typing import Any, Mapping, Sequence

from build_hla2_positive_control_benchmark import _hla_sequences
from build_hla2_positive_control_benchmark_v2 import (
    DEFAULT_OUT as DEFAULT_PACKAGE,
    V1_PACKAGE,
    _checksums,
    read_csv,
    write_csv,
    write_json,
)
from hla2_positive_control_benchmark import (
    EXPOSED_INDICES,
    generate_weight_grid,
    pair_features,
    rank_feature_percentiles,
)
from hla2_positive_control_benchmark_v2 import (
    CLAIM_BOUNDARY_V2,
    FULL_COMPOSITE_FEATURES,
    PILOT_SEEDS,
    STRUCTURAL_FEATURES,
    TCR_FACING_INDICES,
    aggregate_system_results,
    blosum62_similarity,
    build_definitive_ranking_gate,
    build_pilot_attribution_gate,
    validate_specificity_registry,
)
from run_hla2_positive_control_results_analysis import (
    GeometrySample,
    analyze_job_bundle,
    inventory_downloaded_jobs,
    sequence_identity,
    summarize_feature_values,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_OUT = ROOT / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26"
DEFAULT_DOWNLOADS = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Downloads"

METHOD_NAMES = (
    "full_structural_composite",
    "physicochemical_only",
    "tcr_facing_identity",
    "full_core_identity",
    "tcr_facing_blosum62",
    "full_core_blosum62",
    "structural_only",
    "frozen_exposed_ca",
    "random_ranking",
)
NONSTRUCTURAL_METHODS = (
    "physicochemical_only",
    "tcr_facing_identity",
    "full_core_identity",
    "tcr_facing_blosum62",
    "full_core_blosum62",
)
HIGHER_IS_BETTER = {
    "tcr_facing_identity": "tcr_facing_sequence_identity",
    "full_core_identity": "full_core_sequence_identity",
    "tcr_facing_blosum62": "tcr_facing_blosum62_similarity",
    "full_core_blosum62": "full_core_blosum62_similarity",
}


def _panel_key(row: Mapping[str, Any]) -> tuple[str, str, int]:
    return (
        str(row["system_id"]),
        str(row["positive_pair_id"]),
        int(row["panel_seed"]),
    )


def _panels(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, str, int], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[_panel_key(row)].append(dict(row))
    for key, group in grouped.items():
        if sum(str(row["pair_role"]) == "positive" for row in group) != 1:
            raise ValueError(f"panel {key} must contain exactly one positive")
    return grouped


def _rank_metric(group: Sequence[Mapping[str, Any]], metric: str, *, higher: bool = False) -> int:
    ordered = sorted(
        group,
        key=lambda row: (
            -float(row[metric]) if higher else float(row[metric]),
            str(row["pair_id"]),
        ),
    )
    return next(
        index for index, row in enumerate(ordered, start=1)
        if str(row["pair_role"]) == "positive"
    )


def _rank_weighted(
    group: Sequence[Mapping[str, Any]],
    features: Sequence[str],
    weights: Mapping[str, float],
) -> int:
    ranked = rank_feature_percentiles(group, features)
    ordered = sorted(
        ranked,
        key=lambda row: (
            sum(float(weights[feature]) * float(row[f"{feature}_percentile"])
                for feature in features),
            str(row["pair_id"]),
        ),
    )
    return next(
        index for index, row in enumerate(ordered, start=1)
        if str(row["pair_role"]) == "positive"
    )


def _objective(ranks: Sequence[int], lexical: str = "") -> tuple[Any, ...]:
    if not ranks:
        return (math.inf, math.inf, math.inf, lexical)
    return (
        -sum(rank <= 3 for rank in ranks),
        max(ranks),
        -sum(1.0 / rank for rank in ranks) / len(ranks),
        lexical,
    )


def _system_worst_ranks(
    panels: Mapping[tuple[str, str, int], Sequence[Mapping[str, Any]]],
    ranker: Any,
) -> list[int]:
    by_system: dict[str, list[int]] = defaultdict(list)
    for (system_id, _pair_id, _seed), group in panels.items():
        by_system[system_id].append(int(ranker(group)))
    return [max(ranks) for _system, ranks in sorted(by_system.items())]


def _select_weights(
    rows: Sequence[Mapping[str, Any]], features: Sequence[str]
) -> dict[str, float]:
    panels = _panels(rows)
    choices = []
    for weights in generate_weight_grid(features):
        ranks = _system_worst_ranks(
            panels, lambda group: _rank_weighted(group, features, weights)
        )
        nonzero = sum(value > 0 for value in weights.values())
        exposed = float(weights.get("exposed_ca_rmsd_A", 0.0))
        lexical = ";".join(f"{feature}={weights[feature]:.2f}" for feature in features)
        choices.append((
            (
                -sum(rank <= 3 for rank in ranks),
                max(ranks),
                -sum(1.0 / rank for rank in ranks) / len(ranks),
                nonzero,
                -exposed,
                lexical,
            ),
            dict(weights),
        ))
    return min(choices, key=lambda value: value[0])[1]


def _fixed_method_rank(group: Sequence[Mapping[str, Any]], method: str) -> int:
    if method == "physicochemical_only":
        return _rank_metric(group, "tcr_face_physicochemical_mismatch")
    if method in HIGHER_IS_BETTER:
        return _rank_metric(group, HIGHER_IS_BETTER[method], higher=True)
    if method == "frozen_exposed_ca":
        return _rank_metric(group, "exposed_ca_rmsd_A")
    if method == "random_ranking":
        identity = "|".join(str(value) for value in _panel_key(group[0]))
        seed = int(hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16], 16)
        ordered = list(group)
        random.Random(seed).shuffle(ordered)
        return next(
            index for index, row in enumerate(ordered, start=1)
            if str(row["pair_role"]) == "positive"
        )
    raise ValueError(f"unknown fixed method: {method}")


def select_strongest_nonstructural_baseline(rows: Sequence[Mapping[str, Any]]) -> str:
    panels = _panels(rows)
    choices = []
    for method in NONSTRUCTURAL_METHODS:
        ranks = _system_worst_ranks(
            panels, lambda group: _fixed_method_rank(group, method)
        )
        choices.append((_objective(ranks, method), method))
    return min(choices, key=lambda value: value[0])[1]


def evaluate_outer_folds(rows: Sequence[Mapping[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Tune all tunable models on training systems and score only the held-out system."""
    systems = sorted({str(row["system_id"]) for row in rows})
    if len(systems) < 2:
        raise ValueError("outer-fold analysis requires at least two biological systems")
    panel_method_ranks = []
    method_rank_long = []
    fold_weights = []
    ablations = []
    for held_out in systems:
        training = [row for row in rows if str(row["system_id"]) != held_out]
        testing = [row for row in rows if str(row["system_id"]) == held_out]
        training_ids = sorted({str(row["system_id"]) for row in training})
        full_weights = _select_weights(training, FULL_COMPOSITE_FEATURES)
        structural_weights = _select_weights(training, STRUCTURAL_FEATURES)
        baseline = select_strongest_nonstructural_baseline(training)
        structural_weight = sum(full_weights[feature] for feature in STRUCTURAL_FEATURES)
        fold_weights.append({
            "held_out_system_id": held_out,
            "training_system_ids": ";".join(training_ids),
            "selected_nonstructural_baseline": baseline,
            **{f"full_weight_{feature}": full_weights[feature]
               for feature in FULL_COMPOSITE_FEATURES},
            **{f"structural_only_weight_{feature}": structural_weights[feature]
               for feature in STRUCTURAL_FEATURES},
            "full_structural_weight": structural_weight,
            "weights_frozen": False,
        })
        for (system_id, pair_id, seed), group in sorted(_panels(testing).items()):
            ranks = {
                "full_structural_composite": _rank_weighted(
                    group, FULL_COMPOSITE_FEATURES, full_weights
                ),
                "structural_only": _rank_weighted(group, STRUCTURAL_FEATURES, structural_weights),
            }
            for method in METHOD_NAMES:
                if method not in ranks:
                    ranks[method] = _fixed_method_rank(group, method)
                method_rank_long.append({
                    "held_out_system_id": held_out,
                    "training_system_ids": ";".join(training_ids),
                    "system_id": system_id,
                    "positive_pair_id": pair_id,
                    "panel_seed": seed,
                    "method": method,
                    "positive_rank": ranks[method],
                    "capture_at_3": ranks[method] <= 3,
                    "comparison_count": len(group),
                })
            composite_rank = ranks["full_structural_composite"]
            baseline_rank = ranks[baseline]
            ablated_rank = ranks["physicochemical_only"]
            panel_method_ranks.append({
                "held_out_system_id": held_out,
                "training_system_ids": ";".join(training_ids),
                "system_id": system_id,
                "positive_pair_id": pair_id,
                "panel_seed": seed,
                "evaluation_status": "complete",
                "comparison_count": len(group),
                "selected_nonstructural_baseline": baseline,
                "composite_rank": composite_rank,
                "baseline_rank": baseline_rank,
                "ablated_rank": ablated_rank,
                "structural_weight": structural_weight,
                "capture_at_3": composite_rank <= 3,
            })
            ablations.append({
                "held_out_system_id": held_out,
                "system_id": system_id,
                "positive_pair_id": pair_id,
                "panel_seed": seed,
                "full_composite_rank": composite_rank,
                "structural_ablation_rank": ablated_rank,
                "selected_nonstructural_baseline": baseline,
                "selected_nonstructural_rank": baseline_rank,
                "full_structural_weight": structural_weight,
                "credited_improvement": composite_rank < baseline_rank,
                "ablation_removes_improvement": (
                    composite_rank >= baseline_rank or ablated_rank >= baseline_rank
                ),
            })
    return {
        "panel_method_ranks": panel_method_ranks,
        "method_rank_long": method_rank_long,
        "fold_weights": fold_weights,
        "ablations": ablations,
    }


def _batch_jobs(package: Path) -> list[dict[str, Any]]:
    jobs = []
    for path in sorted((package / "alphafold_jobs").glob("hla2_v2_pilot_batch_*_jobs.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, list):
            raise ValueError(f"AlphaFold batch is not an array: {path}")
        jobs.extend(payload)
    if len({str(job["name"]).lower() for job in jobs}) != len(jobs):
        raise ValueError("prepared v2 job names are not unique")
    return jobs


def _feature_matrices(
    package: Path,
    specs: Sequence[Mapping[str, Any]],
    geometry_by_job: Mapping[str, Sequence[GeometrySample]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    pmhc = {row["ligand_id"]: row for row in read_csv(package / "alphafold_jobs/unique_pmhc_inventory.csv")}
    spec_by_key = {
        (str(row["ligand_id"]), int(row["panel_seed"])): row for row in specs
    }
    summaries = []
    ensembles = []
    for comparison in read_csv(package / "controls/comparison_universe.csv"):
        seed = int(comparison["panel_seed"])
        left_meta = pmhc[comparison["left_pmhc_id"]]
        right_meta = pmhc[comparison["right_pmhc_id"]]
        left_spec = spec_by_key[(comparison["left_pmhc_id"], seed)]
        right_spec = spec_by_key[(comparison["right_pmhc_id"], seed)]
        left_samples = geometry_by_job.get(str(left_spec["job_name"]).lower(), ())
        right_samples = geometry_by_job.get(str(right_spec["job_name"]).lower(), ())
        metrics = []
        for left in left_samples:
            for right in right_samples:
                values = pair_features(left.geometry, right.geometry)
                metrics.append(values)
                ensembles.append({
                    "system_id": comparison["system_id"],
                    "positive_pair_id": comparison["positive_pair_id"],
                    "panel_seed": seed,
                    "pair_id": comparison["pair_id"],
                    "left_sample_index": left.sample_index,
                    "right_sample_index": right.sample_index,
                    **{key: round(value, 9) for key, value in values.items()},
                    "interpretation": "Technical AlphaFold ensemble sensitivity, not biological replication.",
                })
        left_core = str(left_meta["core_sequence"])
        right_core = str(right_meta["core_sequence"])
        summary = {
            "system_id": comparison["system_id"],
            "positive_pair_id": comparison["positive_pair_id"],
            "panel_seed": seed,
            "pair_id": comparison["pair_id"],
            "pair_role": "positive" if comparison["comparison_role"] == "positive" else "N3",
            "negative_tier": comparison["negative_tier"],
            "left_id": comparison["left_id"],
            "right_id": comparison["right_id"],
            "left_pmhc_id": comparison["left_pmhc_id"],
            "right_pmhc_id": comparison["right_pmhc_id"],
            "geometry_status": "complete" if metrics else "missing_or_qc_excluded_model",
            "model_combination_count": len(metrics),
            **summarize_feature_values(metrics, FULL_COMPOSITE_FEATURES),
            "tcr_facing_sequence_identity": sequence_identity(
                left_core, right_core, positions=TCR_FACING_INDICES
            ),
            "full_core_sequence_identity": sequence_identity(left_core, right_core),
            "tcr_facing_blosum62_similarity": blosum62_similarity(
                left_core, right_core, positions=TCR_FACING_INDICES
            ),
            "full_core_blosum62_similarity": blosum62_similarity(left_core, right_core),
            "claim_boundary": CLAIM_BOUNDARY_V2,
        }
        summaries.append(summary)
    return summaries, ensembles


def _scoring_rows(feature_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for row in feature_rows:
        if str(row["geometry_status"]) != "complete":
            continue
        output.append({
            "system_id": row["system_id"],
            "positive_pair_id": row["positive_pair_id"],
            "panel_seed": int(row["panel_seed"]),
            "pair_id": row["pair_id"],
            "pair_role": row["pair_role"],
            **{feature: float(row[f"{feature}_median"])
               for feature in FULL_COMPOSITE_FEATURES},
            **{feature: float(row[feature]) for feature in (
                "tcr_facing_sequence_identity", "full_core_sequence_identity",
                "tcr_facing_blosum62_similarity", "full_core_blosum62_similarity",
            )},
        })
    return output


def _permutation_and_uncertainty(
    system_rows: Sequence[Mapping[str, Any]], draws: int = 10000
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    complete = [row for row in system_rows if row.get("evaluation_status") == "complete"]
    if not complete:
        return ([{"evaluation_status": "not_evaluable", "reason": "no_complete_system_results"}],
                [{"evaluation_status": "not_evaluable", "reason": "no_complete_system_results"}])
    differences = [
        1.0 / int(row["system_score"]) - 1.0 / int(row["baseline_system_score"])
        for row in complete
    ]
    observed = sum(differences) / len(differences)
    permutations = [
        sum(sign * value for sign, value in zip(signs, differences)) / len(differences)
        for signs in itertools.product((-1, 1), repeat=len(differences))
    ]
    permutation_rows = [{
        "evaluation_status": "complete",
        "independent_system_count": len(complete),
        "permutation_count": len(permutations),
        "observed_mean_reciprocal_rank_difference": round(observed, 9),
        "two_sided_sign_flip_p_value": round(
            sum(abs(value) >= abs(observed) for value in permutations) / len(permutations), 9
        ),
        "small_system_caveat": "descriptive_not_universal_statistical_proof",
    }]
    generator = random.Random(271828314159)
    bootstraps = []
    for _ in range(draws):
        sample = [generator.choice(differences) for _ in differences]
        bootstraps.append(sum(sample) / len(sample))
    bootstraps.sort()
    lower = bootstraps[int(0.025 * (draws - 1))]
    upper = bootstraps[int(0.975 * (draws - 1))]
    uncertainty_rows = [{
        "evaluation_status": "complete",
        "independent_system_count": len(complete),
        "bootstrap_draw_count": draws,
        "mean_reciprocal_rank_difference": round(observed, 9),
        "bootstrap_95_percent_lower": round(lower, 9),
        "bootstrap_95_percent_upper": round(upper, 9),
        "resampling_unit": "biological_system",
        "small_system_caveat": "descriptive_not_universal_statistical_proof",
    }]
    return permutation_rows, uncertainty_rows


def _write_empty_analysis_tables(out: Path) -> None:
    tables = {
        "qc/model_sample_qc.csv": ("job_name", "sample_index", "sample_status"),
        "qc/job_qc_summary.csv": ("job_name", "technical_status"),
        "benchmark/af3_model_ensemble.csv": ("system_id", "positive_pair_id", "panel_seed", "pair_id"),
        "benchmark/af3_pair_feature_matrix.csv": (
            "system_id", "positive_pair_id", "panel_seed", "pair_id", "pair_role", "geometry_status"
        ),
        "benchmark/panel_method_ranks.csv": (
            "held_out_system_id", "system_id", "positive_pair_id", "panel_seed", "evaluation_status"
        ),
        "benchmark/method_rank_long.csv": (
            "held_out_system_id", "positive_pair_id", "panel_seed", "method", "positive_rank"
        ),
        "benchmark/fold_weights.csv": (
            "held_out_system_id", "training_system_ids", "selected_nonstructural_baseline", "weights_frozen"
        ),
        "benchmark/structural_ablation.csv": (
            "held_out_system_id", "positive_pair_id", "panel_seed", "full_composite_rank", "structural_ablation_rank"
        ),
        "benchmark/system_level_comparisons.csv": (
            "system_id", "evaluation_status", "system_score", "baseline_system_score"
        ),
        "benchmark/pdb_oracle_results.csv": (
            "system_id", "positive_pair_id", "eligible_decoy_count", "oracle_status", "positive_rank"
        ),
    }
    for relative, fields in tables.items():
        write_csv(out / relative, [], fields=fields)


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]], fallback: Sequence[str]) -> None:
    write_csv(path, rows, fields=() if rows else fallback)


def _truth(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "true", "yes"}


def _score_pdb_oracles(
    package: Path,
    availability: Sequence[Mapping[str, Any]],
    fold_weights: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    frozen = {
        (row["positive_pair_id"], row["pair_id"], row["left_ligand_id"], row["right_ligand_id"])
        for row in read_csv(package / "benchmark/pdb_oracle_frozen_pairings.csv")
    }
    prior_geometry = read_csv(V1_PACKAGE / "benchmark/pdb_oracle_feature_matrix.csv")
    observed = {
        (row["positive_pair_id"], row["pair_id"], row["left_ligand_id"], row["right_ligand_id"])
        for row in prior_geometry
    }
    if frozen != observed:
        raise ValueError("PDB oracle geometry rows do not match the pre-geometry frozen pairings")
    weights_by_system = {
        row["held_out_system_id"]: {
            feature: float(row[f"full_weight_{feature}"])
            for feature in FULL_COMPOSITE_FEATURES
        }
        for row in fold_weights
    }
    geometry_by_pair: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in prior_geometry:
        geometry_by_pair[row["positive_pair_id"]].append({
            **dict(row),
            "pair_role": "positive" if row["pair_role"] == "positive" else "N3",
            **{feature: float(row[feature]) for feature in FULL_COMPOSITE_FEATURES},
        })
    output = []
    for source in availability:
        row = dict(source)
        row["mandatory_if_scored"] = _truth(row.get("mandatory_if_scored"))
        if not row["mandatory_if_scored"]:
            row["oracle_status"] = "not_evaluable_availability"
            row["positive_rank"] = ""
        else:
            group = geometry_by_pair.get(str(row["positive_pair_id"]), [])
            weights = weights_by_system.get(str(row["system_id"]))
            if not group or weights is None:
                row["oracle_status"] = "required_pending_results"
                row["positive_rank"] = ""
            else:
                rank = _rank_weighted(group, FULL_COMPOSITE_FEATURES, weights)
                row["positive_rank"] = rank
                row["oracle_status"] = "pass" if rank <= 3 else "fail"
        output.append(row)
    return output


def run_v2_analysis(
    result_roots: Sequence[Path],
    *,
    package: Path = DEFAULT_PACKAGE,
    out: Path = DEFAULT_RESULTS_OUT,
) -> dict[str, Any]:
    """Build a result package from zero, partial, or complete downloaded jobs."""
    out.mkdir(parents=True, exist_ok=True)
    manifest = read_csv(package / "alphafold_jobs/job_manifest.csv")
    batch_jobs = _batch_jobs(package)
    inventory = inventory_downloaded_jobs(manifest, batch_jobs, result_roots)
    write_csv(out / "inventory/job_inventory.csv", inventory)
    complete_count = sum(row["download_status"] == "complete_exact" for row in inventory)
    all_complete = complete_count == len(manifest)
    systems = read_csv(package / "registry/control_system_registry.csv")
    strict_ids = sorted(row["system_id"] for row in systems if row["eligibility"] == "strict")
    required_pairs = defaultdict(list)
    for row in read_csv(package / "registry/positive_pair_registry.csv"):
        required_pairs[row["system_id"]].append(row["pair_id"])

    feature_rows: list[dict[str, Any]] = []
    ensemble_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    job_rows: list[dict[str, Any]] = []
    folds = {"panel_method_ranks": [], "method_rank_long": [], "fold_weights": [], "ablations": []}
    system_rows: list[dict[str, Any]] = []
    if all_complete:
        hla = _hla_sequences()
        specs = []
        for row in inventory:
            alpha = str(row["mhc_alpha_allele"])
            beta = str(row["mhc_beta_allele"])
            specs.append({
                **dict(row),
                "source_cohort": "v2_fresh_attribution_pilot",
                "entity_id": row["ligand_id"],
                **hla[(alpha, beta)],
            })
        geometry_by_job = {}
        for spec in specs:
            samples, geometries, job = analyze_job_bundle(spec)
            sample_rows.extend(samples)
            job_rows.append(job)
            geometry_by_job[str(spec["job_name"]).lower()] = geometries
        feature_rows, ensemble_rows = _feature_matrices(package, specs, geometry_by_job)
        expected_panel_sizes = defaultdict(int)
        for row in read_csv(package / "controls/comparison_universe.csv"):
            expected_panel_sizes[_panel_key(row)] += 1
        scoring = _scoring_rows(feature_rows)
        observed_panel_sizes: dict[tuple[str, str, int], int] = defaultdict(int)
        for row in scoring:
            observed_panel_sizes[_panel_key(row)] += 1
        geometry_complete = expected_panel_sizes == observed_panel_sizes and all(
            size == 26 for size in observed_panel_sizes.values()
        )
        if geometry_complete:
            folds = evaluate_outer_folds(scoring)
            system_rows = aggregate_system_results(
                folds["panel_method_ranks"],
                required_seeds=PILOT_SEEDS,
                required_pairs_by_system=required_pairs,
            )

    if all_complete:
        _write_rows(out / "qc/model_sample_qc.csv", sample_rows, (
            "job_name", "sample_index", "sample_status"
        ))
        _write_rows(out / "qc/job_qc_summary.csv", job_rows, ("job_name", "technical_status"))
        _write_rows(out / "benchmark/af3_model_ensemble.csv", ensemble_rows, (
            "system_id", "positive_pair_id", "panel_seed", "pair_id"
        ))
        _write_rows(out / "benchmark/af3_pair_feature_matrix.csv", feature_rows, (
            "system_id", "positive_pair_id", "panel_seed", "pair_id", "pair_role", "geometry_status"
        ))
        _write_rows(out / "benchmark/panel_method_ranks.csv", folds["panel_method_ranks"], (
            "held_out_system_id", "system_id", "positive_pair_id", "panel_seed", "evaluation_status"
        ))
        _write_rows(out / "benchmark/method_rank_long.csv", folds["method_rank_long"], (
            "held_out_system_id", "positive_pair_id", "panel_seed", "method", "positive_rank"
        ))
        _write_rows(out / "benchmark/fold_weights.csv", folds["fold_weights"], (
            "held_out_system_id", "training_system_ids", "selected_nonstructural_baseline", "weights_frozen"
        ))
        _write_rows(out / "benchmark/structural_ablation.csv", folds["ablations"], (
            "held_out_system_id", "positive_pair_id", "panel_seed", "full_composite_rank", "structural_ablation_rank"
        ))
        _write_rows(out / "benchmark/system_level_comparisons.csv", system_rows, (
            "system_id", "evaluation_status", "system_score", "baseline_system_score"
        ))
    else:
        _write_empty_analysis_tables(out)

    availability = read_csv(package / "benchmark/pdb_oracle_availability.csv")
    oracle_rows = _score_pdb_oracles(package, availability, folds["fold_weights"]) if system_rows else [
        {**row, "mandatory_if_scored": _truth(row.get("mandatory_if_scored"))}
        for row in availability
    ]
    if all_complete:
        _write_rows(out / "benchmark/pdb_oracle_results.csv", oracle_rows, (
            "system_id", "positive_pair_id", "eligible_decoy_count", "oracle_status", "positive_rank"
        ))
    pilot_gate = build_pilot_attribution_gate(
        system_rows, required_system_ids=strict_ids, oracle_rows=oracle_rows
    )
    definitive_gate = build_definitive_ranking_gate(
        system_rows, systems, oracle_rows=oracle_rows
    )
    specificity_rows = read_csv(package / "registry/specificity_negative_registry.csv")
    specificity_gate = validate_specificity_registry(specificity_rows)
    write_json(out / "benchmark/pilot_attribution_gate.json", pilot_gate)
    write_json(out / "benchmark/definitive_ranking_gate.json", definitive_gate)
    write_json(out / "benchmark/specificity_gate.json", specificity_gate)
    write_csv(out / "benchmark/pdb_oracle_availability.csv", oracle_rows)
    permutations, uncertainty = _permutation_and_uncertainty(system_rows)
    write_csv(out / "benchmark/permutation_results.csv", permutations)
    write_csv(out / "benchmark/uncertainty_intervals.csv", uncertainty)

    analysis_manifest = {
        "benchmark_version": "EBV_MS_HLA2_BENCHMARK_V2_PILOT_RESULTS",
        "prepared_job_count": len(manifest),
        "downloaded_complete_exact_job_count": complete_count,
        "all_required_jobs_complete": all_complete,
        "pilot_attribution_status": pilot_gate["pilot_attribution_status"],
        "definitive_status": definitive_gate["definitive_status"],
        "specificity_status": specificity_gate["specificity_status"],
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "discovery_files_read": False,
        "discovery_files_written": False,
        "cross_allele_consensus_created": False,
        "claim_boundary": CLAIM_BOUNDARY_V2,
    }
    write_json(out / "analysis_manifest.json", analysis_manifest)
    (out / "README.md").write_text(
        "# HLA-II benchmark v2 result analysis\n\n"
        f"Pilot attribution status: `{pilot_gate['pilot_attribution_status']}`. "
        f"Exact downloaded jobs: {complete_count}/{len(manifest)}. "
        "Pilot weights are never frozen and discovery remains locked.\n\n"
        f"{CLAIM_BOUNDARY_V2}\n",
        encoding="utf-8",
    )
    write_csv(out / "SHA256SUMS.csv", _checksums(out))
    return {
        "output_directory": str(out),
        "prepared_job_count": len(manifest),
        "downloaded_complete_exact_job_count": complete_count,
        "pilot_attribution_status": pilot_gate["pilot_attribution_status"],
        "definitive_status": definitive_gate["definitive_status"],
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--out", type=Path, default=DEFAULT_RESULTS_OUT)
    parser.add_argument("--result-root", action="append", type=Path, default=[])
    args = parser.parse_args()
    roots = args.result_root or sorted(
        path for path in DEFAULT_DOWNLOADS.glob("*") if path.is_dir()
    )
    print(json.dumps(
        run_v2_analysis(roots, package=args.package, out=args.out),
        indent=2,
        sort_keys=True,
    ))


if __name__ == "__main__":
    main()
