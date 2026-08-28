"""Build additive RMSD sensitivity rankings for V2 and the combined DRB1*15:01 universe."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from analyze_af3_pmhc_downloads import parse_mmcif
from build_combined_drb1501_same_register_rankings import (
    DEFAULT_LEGACY,
    audit_legacy_eligibility,
)
from build_same_register_hla_rankings_v2 import (
    ALLELES,
    CLAIM_BOUNDARY,
    DEFAULT_PANEL,
    read_csv,
    sha256_file,
    write_csv,
    write_json,
)
from same_register_af3_analysis import same_register_geometry


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V2_GEOMETRY = (
    ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25/discovery/pair_summary_6400.csv"
)
DEFAULT_V2_SEQUENCE_RANKS = (
    ROOT / "processed/tcell_library_v2_same_register_hla_rankings_2026-08-27/all_hla_ranked_pairs.csv"
)
DEFAULT_COMBINED = (
    ROOT / "processed/drb1501_combined_same_register_rankings_2026-08-27/combined_ranked_pairs.csv"
)
DEFAULT_LEGACY_SAMPLES = (
    ROOT / "processed/complete_model_pipeline_audit_2026-08-15/canonical_af3_sample_metrics.csv"
)
DEFAULT_BENCHMARK_METHOD_RANKS = (
    ROOT
    / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26/benchmark/"
    "method_rank_long.csv"
)
DEFAULT_OUT = ROOT / "processed/rmsd_rerankings_2026-08-27"
ALLELE_SLUGS = {
    "HLA-DRB1*15:01": "hla_drb1_15_01",
    "HLA-DRB1*13:03": "hla_drb1_13_03",
    "HLA-DRB1*03:01": "hla_drb1_03_01",
    "HLA-DRB1*08:01": "hla_drb1_08_01",
}
RMSD_ENDPOINT = "median P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove fit"
RMSD_CLAIM = (
    "Exploratory structural sensitivity ranking; the frozen exposed-C-alpha RMSD endpoint "
    "did not pass the three-system control benchmark and is not the control-supported primary ranking."
)


def _distribution(values: Sequence[float], prefix: str) -> dict[str, Any]:
    if not values:
        return {f"{prefix}_{field}": "" for field in ("min", "q25", "median", "q75", "max", "iqr")}
    array = np.asarray(values, dtype=float)
    q25, q50, q75 = np.quantile(array, [0.25, 0.5, 0.75])
    return {
        f"{prefix}_min": round(float(np.min(array)), 6),
        f"{prefix}_q25": round(float(q25), 6),
        f"{prefix}_median": round(float(q50), 6),
        f"{prefix}_q75": round(float(q75), 6),
        f"{prefix}_max": round(float(np.max(array)), 6),
        f"{prefix}_iqr": round(float(q75 - q25), 6),
    }


def rank_rmsd_rows(
    rows: Sequence[Mapping[str, Any]], *, group_field: str = "allele"
) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    missing = []
    for source in rows:
        row = dict(source)
        if str(row.get("rmsd_status", "")) == "complete" and str(row.get("exposed_ca_rmsd_A_median", "")):
            groups[str(row[group_field])].append(row)
        else:
            row["rmsd_rank_status"] = "not_ranked_missing_comparable_rmsd"
            missing.append(row)

    ranked: dict[str, list[dict[str, Any]]] = {}
    for group, values in sorted(groups.items()):
        ordered = sorted(
            values,
            key=lambda row: (
                float(row["exposed_ca_rmsd_A_median"]),
                float(row.get("exposed_ca_rmsd_A_iqr", 0.0) or 0.0),
                str(row["pair_id"]),
            ),
        )
        position = 0
        while position < len(ordered):
            end = position + 1
            score = float(ordered[position]["exposed_ca_rmsd_A_median"])
            while end < len(ordered) and float(ordered[end]["exposed_ca_rmsd_A_median"]) == score:
                end += 1
            for row in ordered[position:end]:
                row["rmsd_score_rank"] = position + 1
                row["rmsd_score_tie_size"] = end - position
            position = end
        count = len(ordered)
        for display_rank, row in enumerate(ordered, start=1):
            row["rmsd_rank"] = display_rank
            row["rmsd_evaluable_pair_count"] = count
            row["rmsd_percentile"] = round((int(row["rmsd_score_rank"]) - 1) / max(1, count - 1), 8)
            row["rmsd_rank_status"] = "ranked_comparable_exposed_ca_rmsd"
            row["rmsd_endpoint"] = RMSD_ENDPOINT
            row["rmsd_control_status"] = "failed_three_system_control_benchmark"
        ranked[group] = ordered
    return ranked, sorted(missing, key=lambda row: (str(row.get(group_field, "")), str(row.get("pair_id", ""))))


def _coordinate_label(
    row: Mapping[str, Any], registry: Mapping[str, Mapping[str, Any]]
) -> str:
    left = registry[str(row["ebv_candidate_id"])]
    right = registry[str(row["self_candidate_id"])]
    ebv_protein = str(row["ebv_protein"]).replace("_", "/")
    self_protein = str(row["self_protein"]).replace("_", "/")
    return (
        f"{ebv_protein} {left['source_start_1_based']}-{left['source_end_1_based']} / "
        f"{self_protein} {right['source_start_1_based']}-{right['source_end_1_based']}*"
    )


def build_v2_rmsd_rows(
    geometry_rows: Sequence[Mapping[str, Any]],
    sequence_rows: Sequence[Mapping[str, Any]],
    registry_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    sequence = {str(row["pair_id"]): row for row in sequence_rows}
    registry = {str(row["candidate_id"]): row for row in registry_rows}
    output = []
    for source in geometry_rows:
        row = dict(source)
        sequence_row = sequence[str(row["pair_id"])]
        output.append({
            "allele": row["allele"],
            "pair_id": row["pair_id"],
            "pair_coordinate_label": _coordinate_label(row, registry),
            "ebv_candidate_id": row["ebv_candidate_id"],
            "ebv_protein": row["ebv_protein"],
            "ebv_sequence": row["ebv_sequence"],
            "ebv_core_p1_p9": row["ebv_predicted_core"],
            "self_candidate_id": row["self_candidate_id"],
            "self_protein": row["self_protein"],
            "self_sequence": row["self_sequence"],
            "self_core_p1_p9": row["self_predicted_core"],
            "sequence_hla_rank": sequence_row["hla_rank"],
            "sequence_primary_score": sequence_row["primary_score"],
            "rmsd_status": "complete" if row["geometry_status"] == "complete" else "missing_v2_geometry",
            "rmsd_source": "v2_af3_full_5_by_5_model_ensemble",
            "model_combination_count": row["model_combination_count"],
            "exposed_ca_rmsd_A_min": row["exposed_ca_rmsd_A_min"],
            "exposed_ca_rmsd_A_q25": row["exposed_ca_rmsd_A_q25"],
            "exposed_ca_rmsd_A_median": row["exposed_ca_rmsd_A_median"],
            "exposed_ca_rmsd_A_q75": row["exposed_ca_rmsd_A_q75"],
            "exposed_ca_rmsd_A_max": row["exposed_ca_rmsd_A_max"],
            "exposed_ca_rmsd_A_iqr": row["exposed_ca_rmsd_A_iqr"],
            "full_core_ca_rmsd_A_median": row["full_core_ca_rmsd_A_median"],
            "anchor_ca_rmsd_A_median": row["anchor_ca_rmsd_A_median"],
            "computational_pair_marker": "*",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return output


def compute_legacy_rmsd(
    eligible_pairs: Sequence[Mapping[str, Any]],
    sample_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    samples: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for source in sample_rows:
        path = Path(str(source.get("model_path", "")))
        if (
            str(source.get("af3_cohort", "")) == "legacy_candidate_pmhc"
            and str(source.get("sequence_layout_status", "")) == "pass_exact_three_chain_peptide_match"
            and path.is_file()
        ):
            samples[str(source["candidate_id"])].append(dict(source))
    for candidate_id in samples:
        samples[candidate_id].sort(
            key=lambda row: (str(row["canonical_job_key"]), int(row["sample_index"]), str(row["model_path"]))
        )

    model_cache: dict[str, Any] = {}
    summaries = []
    ensemble = []
    for pair in eligible_pairs:
        ebv_id = str(pair["ebv_candidate_id"])
        self_id = str(pair["human_candidate_id"])
        ebv_samples = samples.get(ebv_id, [])
        self_samples = samples.get(self_id, [])
        values: dict[str, list[float]] = defaultdict(list)
        if ebv_samples and self_samples:
            for left in ebv_samples:
                left_path = str(left["model_path"])
                if left_path not in model_cache:
                    model_cache[left_path] = parse_mmcif(Path(left_path))
                for right in self_samples:
                    right_path = str(right["model_path"])
                    if right_path not in model_cache:
                        model_cache[right_path] = parse_mmcif(Path(right_path))
                    metrics = same_register_geometry(
                        model_cache[left_path],
                        model_cache[right_path],
                        int(pair["ebv_top_core_start_1_based"]),
                        int(pair["human_top_core_start_1_based"]),
                    )
                    for metric, value in metrics.items():
                        values[metric].append(float(value))
                    ensemble.append({
                        "legacy_pair_id": pair["pair_id"],
                        "ebv_candidate_id": ebv_id,
                        "self_candidate_id": self_id,
                        "ebv_job_key": left["canonical_job_key"],
                        "ebv_sample_index": left["sample_index"],
                        "self_job_key": right["canonical_job_key"],
                        "self_sample_index": right["sample_index"],
                        **{key: round(float(value), 6) for key, value in metrics.items()},
                    })
        complete = bool(values)
        summary = {
            "legacy_pair_id": pair["pair_id"],
            "ebv_candidate_id": ebv_id,
            "self_candidate_id": self_id,
            "ebv_core_p1_p9": pair["ebv_top_core_peptide"],
            "self_core_p1_p9": pair["human_top_core_peptide"],
            "rmsd_status": "complete" if complete else "missing_complete_legacy_af3_partner",
            "rmsd_source": "legacy_af3_all_canonical_job_sample_combinations",
            "ebv_unique_job_count": len({row["canonical_job_key"] for row in ebv_samples}),
            "self_unique_job_count": len({row["canonical_job_key"] for row in self_samples}),
            "ebv_sample_count": len(ebv_samples),
            "self_sample_count": len(self_samples),
            "model_combination_count": len(values.get("candidate_exposed_ca_rmsd_A", [])),
        }
        summary.update(_distribution(values.get("candidate_exposed_ca_rmsd_A", []), "exposed_ca_rmsd_A"))
        summary.update(_distribution(values.get("core_p1_p9_ca_rmsd_A", []), "full_core_ca_rmsd_A"))
        summary.update(_distribution(values.get("anchor_ca_rmsd_A", []), "anchor_ca_rmsd_A"))
        summaries.append(summary)
    return summaries, ensemble


def build_combined_rmsd_rows(
    combined_rows: Sequence[Mapping[str, Any]],
    v2_rows: Sequence[Mapping[str, Any]],
    legacy_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    v2 = {str(row["pair_id"]): row for row in v2_rows}
    legacy = {str(row["legacy_pair_id"]): row for row in legacy_rows}
    output = []
    for source in combined_rows:
        row = dict(source)
        geometry: Mapping[str, Any] | None = None
        if str(row.get("v2_pair_id", "")):
            candidate = v2.get(str(row["v2_pair_id"]))
            if candidate and candidate["rmsd_status"] == "complete":
                geometry = candidate
        if geometry is None and str(row.get("legacy_pair_id", "")):
            candidate = legacy.get(str(row["legacy_pair_id"]))
            if candidate and candidate["rmsd_status"] == "complete":
                geometry = candidate
        output.append({
            "allele": row["allele"],
            "pair_id": row["pair_id"],
            "pair_coordinate_label": row["pair_coordinate_label"],
            "source_membership": row["source_membership"],
            "v2_pair_id": row["v2_pair_id"],
            "legacy_pair_id": row["legacy_pair_id"],
            "sequence_combined_rank": row["combined_rank"],
            "sequence_primary_score": row["primary_score"],
            "ebv_sequence": row["ebv_sequence"],
            "ebv_core_p1_p9": row["ebv_core_p1_p9"],
            "self_sequence": row["self_sequence"],
            "self_core_p1_p9": row["self_core_p1_p9"],
            "rmsd_status": geometry["rmsd_status"] if geometry else "missing_comparable_geometry",
            "rmsd_source": geometry["rmsd_source"] if geometry else "",
            "model_combination_count": geometry["model_combination_count"] if geometry else "",
            "exposed_ca_rmsd_A_min": geometry["exposed_ca_rmsd_A_min"] if geometry else "",
            "exposed_ca_rmsd_A_q25": geometry["exposed_ca_rmsd_A_q25"] if geometry else "",
            "exposed_ca_rmsd_A_median": geometry["exposed_ca_rmsd_A_median"] if geometry else "",
            "exposed_ca_rmsd_A_q75": geometry["exposed_ca_rmsd_A_q75"] if geometry else "",
            "exposed_ca_rmsd_A_max": geometry["exposed_ca_rmsd_A_max"] if geometry else "",
            "exposed_ca_rmsd_A_iqr": geometry["exposed_ca_rmsd_A_iqr"] if geometry else "",
            "full_core_ca_rmsd_A_median": geometry["full_core_ca_rmsd_A_median"] if geometry else "",
            "anchor_ca_rmsd_A_median": geometry["anchor_ca_rmsd_A_median"] if geometry else "",
            "computational_pair_marker": "*",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return output


def _control_audit(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    method_rows = [row for row in rows if str(row["method"]) == "frozen_exposed_ca"]
    ranks = [int(row["positive_rank"]) for row in method_rows]
    by_system: dict[str, list[int]] = defaultdict(list)
    for row in method_rows:
        by_system[str(row["system_id"])].append(int(row["positive_rank"]))
    return {
        "method": "frozen_exposed_ca",
        "panel_count": len(ranks),
        "panel_positive_ranks": ranks,
        "panel_capture_at_3_count": sum(rank <= 3 for rank in ranks),
        "system_worst_ranks": {system: max(values) for system, values in sorted(by_system.items())},
        "control_benchmark_status": "fail",
        "use_as_primary_ranking_allowed": False,
        "use_as_exploratory_sensitivity_ranking_allowed": True,
        "claim_boundary": RMSD_CLAIM,
    }


def _results_summary(
    v2_ranked: Mapping[str, Sequence[Mapping[str, Any]]],
    v2_missing: Sequence[Mapping[str, Any]],
    combined: Sequence[Mapping[str, Any]],
    combined_missing: Sequence[Mapping[str, Any]],
    control: Mapping[str, Any],
) -> str:
    lines = [
        "# RMSD sensitivity rankings",
        "",
        f"Endpoint: **{RMSD_ENDPOINT}**. Lower is better.",
        "",
        f"Control benchmark: **{control['panel_capture_at_3_count']}/{control['panel_count']} panels captured at rank 3**; status **fail**. These rankings are structural sensitivity analyses, not replacements for the control-supported sequence ranking.",
        "",
    ]
    for allele in ALLELES:
        rows = list(v2_ranked.get(allele, ()))
        lines.extend([
            f"## {allele}", "",
            f"RMSD-ranked pairs: **{len(rows)}**.", "",
            "| Rank | Epitope pair | Exposed RMSD (A) | Sequence rank |",
            "|---:|---|---:|---:|",
        ])
        for row in rows[:10]:
            lines.append(
                f"| {row['rmsd_rank']} | {row['pair_coordinate_label']} | "
                f"{float(row['exposed_ca_rmsd_A_median']):.3f} | {row['sequence_hla_rank']} |"
            )
        lines.append("")
    lines.extend([
        "## Combined DRB1*15:01 universe", "",
        f"RMSD-ranked pairs: **{len(combined)}**; missing comparable RMSD: **{len(combined_missing)}**.", "",
        "| Rank | Epitope pair | Exposed RMSD (A) | Sequence rank | Source |",
        "|---:|---|---:|---:|---|",
    ])
    for row in combined[:25]:
        lines.append(
            f"| {row['rmsd_rank']} | {row['pair_coordinate_label']} | "
            f"{float(row['exposed_ca_rmsd_A_median']):.3f} | {row['sequence_combined_rank']} | {row['rmsd_source']} |"
        )
    lines.extend([
        "", f"V2 missing comparable RMSD rows across all HLAs: **{len(v2_missing)}**.", "",
        f"> {RMSD_CLAIM}", "", f"> {CLAIM_BOUNDARY}", "",
    ])
    return "\n".join(lines)


def run(
    *,
    v2_geometry_path: Path = DEFAULT_V2_GEOMETRY,
    v2_sequence_path: Path = DEFAULT_V2_SEQUENCE_RANKS,
    panel_path: Path = DEFAULT_PANEL,
    combined_path: Path = DEFAULT_COMBINED,
    legacy_path: Path = DEFAULT_LEGACY,
    legacy_samples_path: Path = DEFAULT_LEGACY_SAMPLES,
    benchmark_path: Path = DEFAULT_BENCHMARK_METHOD_RANKS,
    out: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    v2_rows = build_v2_rmsd_rows(
        read_csv(v2_geometry_path), read_csv(v2_sequence_path), read_csv(panel_path)
    )
    v2_ranked, v2_missing = rank_rmsd_rows(v2_rows)

    legacy_eligible, _legacy_audit = audit_legacy_eligibility(read_csv(legacy_path))
    legacy_summaries, legacy_ensemble = compute_legacy_rmsd(
        legacy_eligible, read_csv(legacy_samples_path)
    )
    combined_rows = build_combined_rmsd_rows(
        read_csv(combined_path), v2_rows, legacy_summaries
    )
    combined_ranked_groups, combined_missing = rank_rmsd_rows(combined_rows)
    combined_ranked = combined_ranked_groups.get("HLA-DRB1*15:01", [])
    control = _control_audit(read_csv(benchmark_path))

    out.mkdir(parents=True, exist_ok=True)
    all_v2_ranked = [row for allele in ALLELES for row in v2_ranked.get(allele, [])]
    write_csv(out / "v2_all_hla_rmsd_ranked_pairs.csv", all_v2_ranked)
    for allele in ALLELES:
        rows = v2_ranked.get(allele, [])
        write_csv(out / f"v2_rankings/{ALLELE_SLUGS[allele]}_rmsd_ranked_pairs.csv", rows)
    write_csv(out / "v2_top_10_rmsd_by_hla.csv", [row for allele in ALLELES for row in v2_ranked.get(allele, [])[:10]])
    write_csv(out / "v2_missing_rmsd.csv", v2_missing)
    write_csv(out / "legacy_rmsd_feature_matrix.csv", legacy_summaries)
    write_csv(out / "legacy_rmsd_model_ensemble.csv", legacy_ensemble)
    write_csv(out / "combined_drb1501_rmsd_ranked_pairs.csv", combined_ranked)
    write_csv(out / "combined_drb1501_top_25_rmsd.csv", combined_ranked[:25], list(combined_ranked[0]))
    write_csv(out / "combined_drb1501_missing_rmsd.csv", combined_missing)
    write_json(out / "rmsd_control_audit.json", control)

    manifest = {
        "analysis_version": "EBV_MS_RMSD_SENSITIVITY_RANKINGS_2026-08-27",
        "endpoint": RMSD_ENDPOINT,
        "endpoint_direction": "lower_is_better",
        "ranking_role": "exploratory_structural_sensitivity_not_primary",
        "control_benchmark_status": control["control_benchmark_status"],
        "v2_ranked_pair_count": len(all_v2_ranked),
        "v2_missing_pair_count": len(v2_missing),
        "v2_ranked_by_hla": {allele: len(v2_ranked.get(allele, [])) for allele in ALLELES},
        "legacy_eligible_pair_count": len(legacy_eligible),
        "legacy_rmsd_complete_pair_count": sum(row["rmsd_status"] == "complete" for row in legacy_summaries),
        "legacy_rmsd_missing_pair_count": sum(row["rmsd_status"] != "complete" for row in legacy_summaries),
        "combined_drb1501_pair_count": len(combined_rows),
        "combined_drb1501_rmsd_ranked_pair_count": len(combined_ranked),
        "combined_drb1501_missing_rmsd_pair_count": len(combined_missing),
        "combined_rmsd_source_counts": dict(sorted(Counter(row["rmsd_source"] for row in combined_ranked).items())),
        "rank_1_pair_by_hla": {
            allele: v2_ranked[allele][0]["pair_id"] for allele in ALLELES
        },
        "combined_drb1501_rank_1_pair_id": combined_ranked[0]["pair_id"],
        "input_sha256": {
            "v2_geometry": sha256_file(v2_geometry_path),
            "v2_sequence_ranks": sha256_file(v2_sequence_path),
            "v2_candidate_registry": sha256_file(panel_path),
            "combined_sequence_rank": sha256_file(combined_path),
            "legacy_pair_universe": sha256_file(legacy_path),
            "legacy_af3_sample_inventory": sha256_file(legacy_samples_path),
            "held_out_control_method_ranks": sha256_file(benchmark_path),
        },
        "claim_boundary": RMSD_CLAIM,
    }
    write_json(out / "analysis_manifest.json", manifest)
    (out / "RESULTS_SUMMARY.md").write_text(
        _results_summary(v2_ranked, v2_missing, combined_ranked, combined_missing, control),
        encoding="utf-8",
    )
    (out / "README.md").write_text(
        """# RMSD sensitivity rankings

This additive package ranks every pair with a comparable structural measurement by median exposed-position P2/P3/P5/P7/P8 C-alpha RMSD after HLA-groove alignment. Lower RMSD ranks first. V2 alleles remain separate; the expanded DRB1*15:01 universe also receives its own structural rank.

Missing RMSDs are never imputed. Legacy whole-local-alignment RMSDs are not used. Eligible legacy RMSDs are recomputed from saved AF3 structures using the same register positions and atom-level endpoint.

The frozen RMSD endpoint failed the completed three-system control benchmark, capturing only 2 of 8 positive panels within the top three. Therefore these outputs are sensitivity analyses and do not replace the control-supported same-register BLOSUM62 rankings.
""" + "\n" + RMSD_CLAIM + "\n" + CLAIM_BOUNDARY + "\n",
        encoding="utf-8",
    )
    checksums = [
        {
            "relative_path": str(path.relative_to(out)),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in sorted(out.rglob("*"), key=str)
        if path.is_file() and path.name != "SHA256SUMS.csv"
    ]
    write_csv(out / "SHA256SUMS.csv", checksums)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(run(out=args.out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
