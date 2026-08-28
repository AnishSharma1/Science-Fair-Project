"""Build HLA-specific sequence-structure concordance rankings.

The ranking combines three same-register sequence features and two structural
features without allowing one evidence family to compensate for poor evidence
in the other. Existing sequence and RMSD packages remain unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from build_same_register_hla_rankings_v2 import (
    ALLELES,
    ALLELE_SLUGS,
    read_csv,
    sha256_file,
    write_csv,
    write_json,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V2_SEQUENCE = (
    ROOT / "processed/tcell_library_v2_same_register_hla_rankings_2026-08-27/all_hla_ranked_pairs.csv"
)
DEFAULT_V2_RMSD = (
    ROOT / "processed/rmsd_rerankings_2026-08-27/v2_all_hla_rmsd_ranked_pairs.csv"
)
DEFAULT_COMBINED_SEQUENCE = (
    ROOT / "processed/drb1501_combined_same_register_rankings_2026-08-27/combined_ranked_pairs.csv"
)
DEFAULT_COMBINED_RMSD = (
    ROOT / "processed/rmsd_rerankings_2026-08-27/combined_drb1501_rmsd_ranked_pairs.csv"
)
DEFAULT_CONTROL_FEATURES = (
    ROOT
    / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26/benchmark/"
    "af3_pair_feature_matrix.csv"
)
DEFAULT_OUT = ROOT / "processed/multifeature_concordance_rankings_2026-08-27"

FEATURE_SPECS = (
    ("tcr_facing_blosum62_similarity", "tcr_facing_blosum62_percentile", True),
    ("tcr_face_physicochemical_mismatch", "physicochemical_mismatch_percentile", False),
    ("tcr_facing_sequence_identity", "same_register_identity_percentile", True),
    ("exposed_ca_rmsd_A_median", "exposed_register_rmsd_percentile", False),
    ("full_core_ca_rmsd_A_median", "full_core_rmsd_percentile", False),
)
SEQUENCE_PERCENTILES = (
    "tcr_facing_blosum62_percentile",
    "physicochemical_mismatch_percentile",
    "same_register_identity_percentile",
)
STRUCTURE_PERCENTILES = (
    "exposed_register_rmsd_percentile",
    "full_core_rmsd_percentile",
)
CLAIM_BOUNDARY = (
    "Descriptive same-register pMHC sequence-structure concordance prioritization only; "
    "not evidence of presentation, TCR binding, activation, cross-reactivity, molecular "
    "mimicry, MS mechanism, probability, or false-discovery rate."
)
CONTROL_BOUNDARY = (
    "Retrospective exploratory audit on the three existing positive-control systems; it "
    "does not freeze weights, satisfy the six-system definitive gate, or unlock discovery."
)


def _float(row: Mapping[str, Any], field: str) -> float:
    value = row.get(field, "")
    if value is None or str(value).strip() == "":
        raise ValueError(f"missing required numeric feature {field} for {row.get('pair_id', '')}")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite feature {field} for {row.get('pair_id', '')}")
    return number


def _average_tie_percentiles(
    rows: Sequence[Mapping[str, Any]], field: str, *, higher_is_better: bool
) -> dict[str, float]:
    """Return average-rank percentiles from 0 (best) to 1 (worst)."""
    ordered = sorted(
        rows,
        key=lambda row: (
            -_float(row, field) if higher_is_better else _float(row, field),
            str(row["pair_id"]),
        ),
    )
    output: dict[str, float] = {}
    position = 0
    denominator = max(1, len(ordered) - 1)
    while position < len(ordered):
        end = position + 1
        value = _float(ordered[position], field)
        while end < len(ordered) and _float(ordered[end], field) == value:
            end += 1
        average_zero_based_rank = (position + end - 1) / 2.0
        percentile = round(average_zero_based_rank / denominator, 12)
        for row in ordered[position:end]:
            output[str(row["pair_id"])] = percentile
        position = end
    return output


def _evidence_band(sequence_percentile: float, structure_percentile: float) -> str:
    if sequence_percentile <= 0.10 and structure_percentile <= 0.10:
        return "both_families_top_10_percent"
    if sequence_percentile <= 0.25 and structure_percentile <= 0.25:
        return "both_families_top_25_percent"
    if sequence_percentile <= 0.25:
        return "sequence_stronger_than_structure"
    if structure_percentile <= 0.25:
        return "structure_stronger_than_sequence"
    return "neither_family_top_25_percent"


def build_concordance_ranking(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Rank one HLA/evaluable universe using balanced evidence families."""
    if len({str(row.get("pair_id", "")) for row in rows}) != len(rows):
        raise ValueError("pair IDs must be unique within one concordance universe")

    complete: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    required = [field for field, _percentile, _higher in FEATURE_SPECS]
    for source in rows:
        row = dict(source)
        status = str(row.get("multifeature_input_status", ""))
        has_all_features = all(str(row.get(field, "")).strip() for field in required)
        if status == "complete" and has_all_features:
            for field in required:
                _float(row, field)
            complete.append(row)
        else:
            row["concordance_rank_status"] = "not_ranked_missing_comparable_structure"
            missing.append(row)

    percentile_maps = {
        output_field: _average_tie_percentiles(
            complete, source_field, higher_is_better=higher_is_better
        )
        for source_field, output_field, higher_is_better in FEATURE_SPECS
    }
    scored = []
    for row in complete:
        pair_id = str(row["pair_id"])
        for _source_field, output_field, _higher in FEATURE_SPECS:
            row[output_field] = percentile_maps[output_field][pair_id]
        sequence = sum(float(row[field]) for field in SEQUENCE_PERCENTILES) / 3.0
        structure = sum(float(row[field]) for field in STRUCTURE_PERCENTILES) / 2.0
        row["sequence_family_percentile"] = round(sequence, 12)
        row["structure_family_percentile"] = round(structure, 12)
        row["worst_family_percentile"] = round(max(sequence, structure), 12)
        row["balanced_multifeature_percentile"] = round((sequence + structure) / 2.0, 12)
        row["concordance_evidence_band"] = _evidence_band(sequence, structure)
        scored.append(row)

    ordered = sorted(
        scored,
        key=lambda row: (
            float(row["worst_family_percentile"]),
            float(row["balanced_multifeature_percentile"]),
            str(row["pair_id"]),
        ),
    )
    position = 0
    while position < len(ordered):
        end = position + 1
        score = (
            float(ordered[position]["worst_family_percentile"]),
            float(ordered[position]["balanced_multifeature_percentile"]),
        )
        while end < len(ordered) and (
            float(ordered[end]["worst_family_percentile"]),
            float(ordered[end]["balanced_multifeature_percentile"]),
        ) == score:
            end += 1
        for row in ordered[position:end]:
            row["concordance_score_rank"] = position + 1
            row["concordance_score_tie_size"] = end - position
        position = end

    count = len(ordered)
    for display_rank, row in enumerate(ordered, start=1):
        row["concordance_rank"] = display_rank
        row["concordance_evaluable_pair_count"] = count
        row["concordance_percentile"] = round(
            (int(row["concordance_score_rank"]) - 1) / max(1, count - 1), 12
        )
        row["concordance_rank_status"] = "ranked_complete_sequence_and_structure"
        row["primary_objective"] = "minimize_worst_sequence_or_structure_family_percentile"
        row["secondary_objective"] = "minimize_balanced_sequence_structure_mean_percentile"
        row["deterministic_tie_break"] = "lexical_pair_id"
    return ordered, sorted(missing, key=lambda row: str(row.get("pair_id", "")))


def _join_sequence_and_rmsd(
    sequence_rows: Sequence[Mapping[str, Any]],
    rmsd_rows: Sequence[Mapping[str, Any]],
    *, sequence_rank_field: str,
) -> list[dict[str, Any]]:
    rmsd = {str(row["pair_id"]): row for row in rmsd_rows}
    output = []
    for source in sequence_rows:
        row = dict(source)
        geometry = rmsd.get(str(row["pair_id"]))
        complete = bool(
            geometry
            and str(geometry.get("rmsd_status", "")) == "complete"
            and str(geometry.get("exposed_ca_rmsd_A_median", "")).strip()
            and str(geometry.get("full_core_ca_rmsd_A_median", "")).strip()
        )
        output.append({
            "allele": row["allele"],
            "pair_id": row["pair_id"],
            "pair_coordinate_label": (
                geometry.get("pair_coordinate_label", "") if geometry else row.get("pair_coordinate_label", "")
            ),
            "source_membership": row.get("source_membership", "v2"),
            "v2_pair_id": row.get("v2_pair_id", row["pair_id"]),
            "legacy_pair_id": row.get("legacy_pair_id", ""),
            "ebv_candidate_id": row.get("ebv_candidate_id", ""),
            "ebv_protein": row.get("ebv_protein", ""),
            "ebv_sequence": row.get("ebv_sequence", ""),
            "ebv_core_p1_p9": row.get("ebv_core_p1_p9", row.get("ebv_predicted_core", "")),
            "self_candidate_id": row.get("self_candidate_id", ""),
            "self_protein": row.get("self_protein", ""),
            "self_sequence": row.get("self_sequence", ""),
            "self_core_p1_p9": row.get("self_core_p1_p9", row.get("self_predicted_core", "")),
            "sequence_rank": row.get(sequence_rank_field, ""),
            "tcr_face_physicochemical_mismatch": row["tcr_face_physicochemical_mismatch"],
            "tcr_facing_blosum62_similarity": row["tcr_facing_blosum62_similarity"],
            "tcr_facing_sequence_identity": row["tcr_facing_sequence_identity"],
            "full_core_sequence_identity": row["full_core_sequence_identity"],
            "full_core_blosum62_similarity": row.get("full_core_blosum62_similarity", ""),
            "rmsd_rank": geometry.get("rmsd_rank", "") if geometry else "",
            "rmsd_source": geometry.get("rmsd_source", "") if geometry else "",
            "model_combination_count": geometry.get("model_combination_count", "") if geometry else "",
            "exposed_ca_rmsd_A_median": geometry.get("exposed_ca_rmsd_A_median", "") if geometry else "",
            "exposed_ca_rmsd_A_iqr": geometry.get("exposed_ca_rmsd_A_iqr", "") if geometry else "",
            "full_core_ca_rmsd_A_median": geometry.get("full_core_ca_rmsd_A_median", "") if geometry else "",
            "anchor_ca_rmsd_A_median": geometry.get("anchor_ca_rmsd_A_median", "") if geometry else "",
            "multifeature_input_status": "complete" if complete else "missing_comparable_structure",
            "computational_pair_marker": "*",
            "pair_evidence_status": "computational_pair_no_exact_paired_recognition_evidence",
            "claim_boundary": CLAIM_BOUNDARY,
        })
    return output


def _rank_groups(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row["allele"])].append(dict(row))
    ranked: dict[str, list[dict[str, Any]]] = {}
    missing = []
    for allele, values in sorted(groups.items()):
        ranked[allele], group_missing = build_concordance_ranking(values)
        missing.extend(group_missing)
    return ranked, sorted(missing, key=lambda row: (str(row["allele"]), str(row["pair_id"])))


def _control_panel_audit(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    panels: dict[tuple[str, str, int], list[dict[str, Any]]] = defaultdict(list)
    for source in rows:
        if str(source.get("geometry_status", "")) != "complete":
            continue
        exposed = _float(source, "exposed_ca_rmsd_A_median")
        anchor = _float(source, "anchor_ca_rmsd_A_median")
        full_core = math.sqrt((5.0 * exposed * exposed + 4.0 * anchor * anchor) / 9.0)
        row = {
            "allele": "control_exact_hla_panel",
            "pair_id": source["pair_id"],
            "pair_role": source["pair_role"],
            "system_id": source["system_id"],
            "positive_pair_id": source["positive_pair_id"],
            "panel_seed": int(source["panel_seed"]),
            "tcr_facing_blosum62_similarity": source["tcr_facing_blosum62_similarity"],
            "tcr_face_physicochemical_mismatch": source["tcr_face_physicochemical_mismatch_median"],
            "tcr_facing_sequence_identity": source["tcr_facing_sequence_identity"],
            "full_core_sequence_identity": source["full_core_sequence_identity"],
            "exposed_ca_rmsd_A_median": exposed,
            "full_core_ca_rmsd_A_median": round(full_core, 12),
            "anchor_ca_rmsd_A_median": anchor,
            "multifeature_input_status": "complete",
        }
        key = (str(source["system_id"]), str(source["positive_pair_id"]), int(source["panel_seed"]))
        panels[key].append(row)

    panel_results = []
    for (system_id, positive_pair_id, seed), panel_rows in sorted(panels.items()):
        ranked, missing = build_concordance_ranking(panel_rows)
        positives = [row for row in ranked if str(row["pair_role"]) == "positive"]
        if len(panel_rows) != 26 or missing or len(positives) != 1:
            raise ValueError(
                f"control panel {system_id}/{positive_pair_id}/{seed} is not one complete 26-pair panel"
            )
        positive = positives[0]
        positive_row_id = str(positive["pair_id"])
        panel_results.append({
            "system_id": system_id,
            "positive_pair_id": positive_pair_id,
            "panel_seed": seed,
            "comparison_count": len(panel_rows),
            "positive_concordance_rank": positive["concordance_score_rank"],
            "positive_display_rank": positive["concordance_rank"],
            "positive_worst_family_percentile": positive["worst_family_percentile"],
            "positive_balanced_multifeature_percentile": positive["balanced_multifeature_percentile"],
            "positive_sequence_family_percentile": positive["sequence_family_percentile"],
            "positive_structure_family_percentile": positive["structure_family_percentile"],
            "positive_blosum62_feature_rank": _feature_rank(panel_rows, positive_row_id, "tcr_facing_blosum62_similarity", True),
            "positive_exposed_rmsd_feature_rank": _feature_rank(panel_rows, positive_row_id, "exposed_ca_rmsd_A_median", False),
            "capture_at_3": int(positive["concordance_score_rank"]) <= 3,
            "audit_role": "retrospective_exploratory_only",
        })

    by_system: dict[str, list[int]] = defaultdict(list)
    for row in panel_results:
        by_system[str(row["system_id"])].append(int(row["positive_concordance_rank"]))
    system_worst = {system: max(values) for system, values in sorted(by_system.items())}
    audit = {
        "analysis": "five_feature_sequence_structure_concordance",
        "panel_count": len(panel_results),
        "independent_system_count": len(system_worst),
        "panel_capture_at_3_count": sum(bool(row["capture_at_3"]) for row in panel_results),
        "all_panels_capture_at_3": all(bool(row["capture_at_3"]) for row in panel_results),
        "system_worst_ranks": system_worst,
        "all_systems_capture_at_3": all(rank <= 3 for rank in system_worst.values()),
        "validation_status": "retrospective_exploratory_not_definitive",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "control_boundary": CONTROL_BOUNDARY,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return panel_results, audit


def _feature_rank(
    rows: Sequence[Mapping[str, Any]], pair_id: str, field: str, higher_is_better: bool
) -> int:
    ordered = sorted(
        rows,
        key=lambda row: (
            -_float(row, field) if higher_is_better else _float(row, field),
            str(row["pair_id"]),
        ),
    )
    target = next(_float(row, field) for row in ordered if str(row["pair_id"]) == pair_id)
    return 1 + sum(
        (_float(row, field) > target if higher_is_better else _float(row, field) < target)
        for row in ordered
    )


def _results_summary(
    v2_ranked: Mapping[str, Sequence[Mapping[str, Any]]],
    v2_missing: Sequence[Mapping[str, Any]],
    combined: Sequence[Mapping[str, Any]],
    combined_missing: Sequence[Mapping[str, Any]],
    control: Mapping[str, Any],
) -> str:
    lines = [
        "# Multi-feature sequence-structure concordance rankings",
        "",
        "Five displayed-pMHC proxies are shown and combined: TCR-facing physicochemical mismatch, "
        "TCR-facing BLOSUM62, TCR-facing same-register identity, exposed-register RMSD, and full "
        "P1-P9 core RMSD. Lower percentiles and lower concordance ranks are better.",
        "",
        "The primary objective minimizes the worse of the sequence-family and structure-family "
        "percentiles. This rewards pairs that are jointly good in both families instead of allowing "
        "one excellent family to hide one poor family.",
        "",
        f"Retrospective control audit: **{control['panel_capture_at_3_count']}/{control['panel_count']} "
        "panels captured at rank 3**. This audit is exploratory and cannot freeze weights or unlock discovery.",
        "",
    ]
    for allele in ALLELES:
        rows = list(v2_ranked.get(allele, ()))
        lines.extend([
            f"## {allele}",
            "",
            f"Ranked: **{len(rows)}**.",
            "",
            "| Rank | Epitope pair | Sequence family | Structure family | Worst family | Exposed RMSD (A) | Full-core RMSD (A) |",
            "|---:|---|---:|---:|---:|---:|---:|",
        ])
        for row in rows[:10]:
            lines.append(
                f"| {row['concordance_rank']} | {row['pair_coordinate_label']} | "
                f"{float(row['sequence_family_percentile']):.3f} | "
                f"{float(row['structure_family_percentile']):.3f} | "
                f"{float(row['worst_family_percentile']):.3f} | "
                f"{float(row['exposed_ca_rmsd_A_median']):.3f} | "
                f"{float(row['full_core_ca_rmsd_A_median']):.3f} |"
            )
        lines.append("")
    lines.extend([
        "## Combined DRB1*15:01 universe",
        "",
        f"Ranked: **{len(combined)}**; missing comparable structure: **{len(combined_missing)}**.",
        "",
        "| Rank | Epitope pair | Sequence family | Structure family | Worst family | Sequence rank | RMSD rank |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ])
    for row in combined[:25]:
        lines.append(
            f"| {row['concordance_rank']} | {row['pair_coordinate_label']} | "
            f"{float(row['sequence_family_percentile']):.3f} | "
            f"{float(row['structure_family_percentile']):.3f} | "
            f"{float(row['worst_family_percentile']):.3f} | {row['sequence_rank']} | {row['rmsd_rank']} |"
        )
    lines.extend([
        "",
        f"V2 missing comparable structure across all HLAs: **{len(v2_missing)}**.",
        "",
        "> Full-core RMSD means P1-P9 peptide-core C-alpha RMSD after HLA-groove alignment. It is not whole-source-protein RMSD and not whole-pMHC RMSD.",
        "",
        f"> {CLAIM_BOUNDARY}",
        "",
    ])
    return "\n".join(lines)


def run(
    *,
    v2_sequence_path: Path = DEFAULT_V2_SEQUENCE,
    v2_rmsd_path: Path = DEFAULT_V2_RMSD,
    combined_sequence_path: Path = DEFAULT_COMBINED_SEQUENCE,
    combined_rmsd_path: Path = DEFAULT_COMBINED_RMSD,
    control_features_path: Path = DEFAULT_CONTROL_FEATURES,
    out: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    v2_inputs = _join_sequence_and_rmsd(
        read_csv(v2_sequence_path), read_csv(v2_rmsd_path), sequence_rank_field="hla_rank"
    )
    v2_ranked, v2_missing = _rank_groups(v2_inputs)
    combined_inputs = _join_sequence_and_rmsd(
        read_csv(combined_sequence_path),
        read_csv(combined_rmsd_path),
        sequence_rank_field="combined_rank",
    )
    combined_groups, combined_missing = _rank_groups(combined_inputs)
    combined = combined_groups.get("HLA-DRB1*15:01", [])
    control_rows, control = _control_panel_audit(read_csv(control_features_path))

    out.mkdir(parents=True, exist_ok=True)
    all_v2 = [row for allele in ALLELES for row in v2_ranked.get(allele, [])]
    write_csv(out / "v2_all_hla_multifeature_ranked_pairs.csv", all_v2)
    for allele in ALLELES:
        write_csv(
            out / f"v2_rankings/{ALLELE_SLUGS[allele]}_multifeature_ranked_pairs.csv",
            v2_ranked.get(allele, []),
        )
    write_csv(
        out / "v2_top_10_multifeature_by_hla.csv",
        [row for allele in ALLELES for row in v2_ranked.get(allele, [])[:10]],
    )
    write_csv(out / "v2_missing_multifeature.csv", v2_missing)
    write_csv(out / "combined_drb1501_multifeature_ranked_pairs.csv", combined)
    write_csv(out / "combined_drb1501_top_25_multifeature.csv", combined[:25])
    write_csv(out / "combined_drb1501_missing_multifeature.csv", combined_missing)
    write_csv(out / "control_panel_concordance_ranks.csv", control_rows)
    write_json(out / "control_concordance_audit.json", control)

    protocol = {
        "analysis_version": "EBV_MS_MULTIFEATURE_CONCORDANCE_2026-08-27",
        "rank_scope": "within_each_HLA_or_declared_DRB1501_universe_only",
        "tcr_facing_positions": ["P2", "P3", "P5", "P7", "P8"],
        "anchor_positions": ["P1", "P4", "P6", "P9"],
        "feature_directions": {
            "tcr_facing_blosum62_similarity": "higher_is_better",
            "tcr_face_physicochemical_mismatch": "lower_is_better",
            "tcr_facing_sequence_identity": "higher_is_better",
            "exposed_ca_rmsd_A_median": "lower_is_better",
            "full_core_ca_rmsd_A_median": "lower_is_better",
        },
        "feature_normalization": "average_tie_rank_percentile_within_evaluable_universe",
        "sequence_family": "equal_mean_of_blosum62_physicochemical_and_identity_percentiles",
        "structure_family": "equal_mean_of_exposed_register_and_full_core_rmsd_percentiles",
        "primary_objective": "minimize_max_sequence_family_and_structure_family_percentile",
        "secondary_objective": "minimize_equal_mean_of_sequence_and_structure_family_percentiles",
        "full_core_rmsd_definition": "P1-P9_peptide_CA_RMSD_after_HLA_groove_alignment",
        "whole_source_protein_rmsd_available": False,
        "binding_percentile_used": False,
        "missing_structure_imputed": False,
        "cross_allele_consensus_used": False,
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "control_audit_role": "retrospective_exploratory_only",
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "scoring_protocol.json", protocol)

    manifest = {
        **protocol,
        "v2_input_pair_count": len(v2_inputs),
        "v2_ranked_pair_count": len(all_v2),
        "v2_missing_structure_pair_count": len(v2_missing),
        "v2_ranked_by_hla": {allele: len(v2_ranked.get(allele, [])) for allele in ALLELES},
        "combined_drb1501_input_pair_count": len(combined_inputs),
        "combined_drb1501_ranked_pair_count": len(combined),
        "combined_drb1501_missing_structure_pair_count": len(combined_missing),
        "combined_source_membership_counts": dict(sorted(Counter(str(row["source_membership"]) for row in combined).items())),
        "rank_1_pair_by_hla": {
            allele: v2_ranked[allele][0]["pair_id"] for allele in ALLELES
        },
        "combined_drb1501_rank_1_pair_id": combined[0]["pair_id"],
        "control_audit": control,
        "input_sha256": {
            "v2_sequence_rankings": sha256_file(v2_sequence_path),
            "v2_rmsd_rankings": sha256_file(v2_rmsd_path),
            "combined_drb1501_sequence_rankings": sha256_file(combined_sequence_path),
            "combined_drb1501_rmsd_rankings": sha256_file(combined_rmsd_path),
            "positive_control_feature_matrix": sha256_file(control_features_path),
        },
    }
    write_json(out / "analysis_manifest.json", manifest)
    (out / "RESULTS_SUMMARY.md").write_text(
        _results_summary(v2_ranked, v2_missing, combined, combined_missing, control),
        encoding="utf-8",
    )
    (out / "README.md").write_text(
        """# Multi-feature sequence-structure concordance rankings

This additive package places five interpretable measurements beside every pair and ranks only pairs with complete same-register sequence and comparable structural data. Each HLA remains separate. The expanded DRB1*15:01 universe is also ranked separately.

The three sequence metrics are converted to within-universe percentiles and averaged. The two structural RMSDs are converted the same way and averaged. The primary rank minimizes the worse family percentile, so a pair must be reasonably strong in both sequence and structure. The balanced family mean is the secondary score and lexical pair ID is the deterministic final display tie-break.

"Whole-protein RMSD" is not available because the AlphaFold models contain pMHC complexes, not the complete EBV and human source proteins. The scientifically relevant broad structural endpoint here is full P1-P9 peptide-core C-alpha RMSD after aligning equivalent HLA-groove atoms. Whole-pMHC RMSD is not used because the shared HLA scaffold would dominate it.

The five features describe complementary proxies of the displayed pMHC surface. They are not a mechanistic model of what a TCR "considers," and no TCR sequence, TCR contact map, binding affinity, or activation assay enters this ranking. The control audit is retrospective and exploratory; it cannot freeze weights or unlock discovery.

Reproduce with:

```bash
PYTHONPATH=src python3 src/build_multifeature_concordance_rankings.py
```
"""
        + "\n"
        + CONTROL_BOUNDARY
        + "\n"
        + CLAIM_BOUNDARY
        + "\n",
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
