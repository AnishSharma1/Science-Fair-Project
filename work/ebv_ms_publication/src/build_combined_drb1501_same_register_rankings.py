"""Combine V2 and eligible legacy DRB1*15:01 pairs under one locked score."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any, Mapping, Sequence

from build_same_register_hla_rankings_v2 import (
    CLAIM_BOUNDARY,
    DEFAULT_BENCHMARK,
    read_csv,
    select_control_supported_method,
    sequence_metrics,
    sha256_file,
    write_csv,
    write_json,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_V2 = (
    ROOT
    / "processed/tcell_library_v2_same_register_hla_rankings_2026-08-27/"
    "rankings/hla_drb1_15_01_ranked_pairs.csv"
)
DEFAULT_LEGACY = ROOT / "processed/register_aware_scoring/register_aware_pair_scores.csv"
DEFAULT_V2_REGISTRY = ROOT / "processed/tcell_library_v2_2026-08-22/frozen_v2_80_peptide_panel.csv"
DEFAULT_LEGACY_ANNOTATIONS = (
    ROOT / "processed/protein_region_annotations/candidate_protein_region_annotations.csv"
)
DEFAULT_OUT = ROOT / "processed/drb1501_combined_same_register_rankings_2026-08-27"
ALLELE = "HLA-DRB1*15:01"
VALID_LEGACY_SELF_REGISTERS = {
    "iedb_top_core_hypothesis",
    "experimental_primary_allele_reference",
}


def audit_legacy_eligibility(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    eligible = []
    audit = []
    for source in rows:
        row = dict(source)
        ebv_status = str(row.get("ebv_register_status", ""))
        self_status = str(row.get("human_register_status", ""))
        ebv_core = str(row.get("ebv_top_core_peptide", ""))
        self_core = str(row.get("human_top_core_peptide", ""))
        if ebv_status == "calibration_only_nonprimary_allele":
            status = "excluded_nonprimary_hla"
        elif ebv_status != "iedb_top_core_hypothesis" or len(ebv_core) != 9:
            status = "excluded_unresolved_ebv_register"
        elif self_status not in VALID_LEGACY_SELF_REGISTERS or len(self_core) != 9:
            status = "excluded_unresolved_self_register"
        else:
            status = "eligible_primary_drb1501_resolved_registers"
            eligible.append(row)
        audit.append({
            "pair_id": row.get("pair_id", ""),
            "ebv_candidate_id": row.get("ebv_candidate_id", ""),
            "human_candidate_id": row.get("human_candidate_id", ""),
            "ebv_register_status": ebv_status,
            "human_register_status": self_status,
            "ebv_core_p1_p9": ebv_core,
            "self_core_p1_p9": self_core,
            "legacy_score_coverage_status": row.get("score_coverage_status", ""),
            "combined_eligibility_status": status,
        })
    return eligible, audit


def _coordinate_label(
    protein: str,
    source: Mapping[str, Any] | None,
    *,
    start_field: str,
    end_field: str,
) -> tuple[str, str, str]:
    display_protein = protein.replace("_", "/")
    if source is None:
        return display_protein, "", ""
    start = str(source.get(start_field, ""))
    end = str(source.get(end_field, ""))
    return display_protein, start, end


def _v2_record(
    row: Mapping[str, Any], registry: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    ebv_source = registry.get(str(row["ebv_candidate_id"]))
    self_source = registry.get(str(row["self_candidate_id"]))
    ebv_protein, ebv_start, ebv_end = _coordinate_label(
        str(row["ebv_protein"]), ebv_source,
        start_field="source_start_1_based", end_field="source_end_1_based",
    )
    self_protein, self_start, self_end = _coordinate_label(
        str(row["self_protein"]), self_source,
        start_field="source_start_1_based", end_field="source_end_1_based",
    )
    metrics = sequence_metrics(str(row["ebv_predicted_core"]), str(row["self_predicted_core"]))
    return {
        "allele": ALLELE,
        "pair_id": row["pair_id"],
        "source_membership": "v2_only",
        "v2_pair_id": row["pair_id"],
        "legacy_pair_id": "",
        "v2_original_hla_rank": row.get("hla_rank", ""),
        "legacy_original_score_coverage_status": "",
        "legacy_pair_validation": "",
        "ebv_candidate_id": row["ebv_candidate_id"],
        "ebv_legacy_candidate_id": "",
        "ebv_protein": ebv_protein,
        "ebv_start_1_based": ebv_start,
        "ebv_end_1_based": ebv_end,
        "ebv_sequence": row["ebv_sequence"],
        "ebv_core_p1_p9": row["ebv_predicted_core"],
        "ebv_binding_percentile_rank": row.get("ebv_binding_percentile_rank", ""),
        "self_candidate_id": row["self_candidate_id"],
        "self_legacy_candidate_id": "",
        "self_protein": self_protein,
        "self_start_1_based": self_start,
        "self_end_1_based": self_end,
        "self_sequence": row["self_sequence"],
        "self_core_p1_p9": row["self_predicted_core"],
        "self_binding_percentile_rank": row.get("self_binding_percentile_rank", ""),
        **metrics,
    }


def _legacy_record(
    row: Mapping[str, Any], annotations: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    ebv_source = annotations.get(str(row["ebv_candidate_id"]))
    self_source = annotations.get(str(row["human_candidate_id"]))
    if ebv_source is None or self_source is None:
        raise ValueError(f"legacy candidate lacks protein-coordinate annotation: {row['pair_id']}")
    if str(ebv_source.get("peptide", "")) != str(row["ebv_peptide"]):
        raise ValueError(f"legacy EBV annotation sequence mismatch: {row['pair_id']}")
    if str(self_source.get("peptide", "")) != str(row["human_peptide"]):
        raise ValueError(f"legacy self annotation sequence mismatch: {row['pair_id']}")
    ebv_protein, ebv_start, ebv_end = _coordinate_label(
        str(ebv_source["short_protein_name"]), ebv_source,
        start_field="parent_residue_start_1_based", end_field="parent_residue_end_1_based",
    )
    self_protein, self_start, self_end = _coordinate_label(
        str(self_source["short_protein_name"]), self_source,
        start_field="parent_residue_start_1_based", end_field="parent_residue_end_1_based",
    )
    metrics = sequence_metrics(str(row["ebv_top_core_peptide"]), str(row["human_top_core_peptide"]))
    return {
        "allele": ALLELE,
        "pair_id": f"{ALLELE}|LEGACY|{row['pair_id']}",
        "source_membership": "legacy_only",
        "v2_pair_id": "",
        "legacy_pair_id": row["pair_id"],
        "v2_original_hla_rank": "",
        "legacy_original_score_coverage_status": row.get("score_coverage_status", ""),
        "legacy_pair_validation": row.get("pair_validation", ""),
        "ebv_candidate_id": row["ebv_candidate_id"],
        "ebv_legacy_candidate_id": row["ebv_candidate_id"],
        "ebv_protein": ebv_protein,
        "ebv_start_1_based": ebv_start,
        "ebv_end_1_based": ebv_end,
        "ebv_sequence": row["ebv_peptide"],
        "ebv_core_p1_p9": row["ebv_top_core_peptide"],
        "ebv_binding_percentile_rank": row.get("ebv_binding_rank", ""),
        "self_candidate_id": row["human_candidate_id"],
        "self_legacy_candidate_id": row["human_candidate_id"],
        "self_protein": self_protein,
        "self_start_1_based": self_start,
        "self_end_1_based": self_end,
        "self_sequence": row["human_peptide"],
        "self_core_p1_p9": row["human_top_core_peptide"],
        "self_binding_percentile_rank": row.get("human_binding_rank", ""),
        **metrics,
    }


def combine_and_rank(
    v2_rows: Sequence[Mapping[str, Any]],
    legacy_rows: Sequence[Mapping[str, Any]],
    *,
    v2_registry: Sequence[Mapping[str, Any]] = (),
    legacy_annotations: Sequence[Mapping[str, Any]] = (),
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    v2_reg = {str(row.get("candidate_id", "")): row for row in v2_registry}
    legacy_reg = {str(row.get("candidate_id", "")): row for row in legacy_annotations}
    canonical: dict[tuple[str, str], dict[str, Any]] = {}
    for source in v2_rows:
        record = _v2_record(source, v2_reg)
        key = (str(record["ebv_sequence"]), str(record["self_sequence"]))
        if key in canonical:
            raise ValueError("duplicate exact peptide pair within V2")
        canonical[key] = record

    overlaps = []
    for source in legacy_rows:
        legacy = _legacy_record(source, legacy_reg)
        key = (str(legacy["ebv_sequence"]), str(legacy["self_sequence"]))
        if key not in canonical:
            canonical[key] = legacy
            continue
        current = canonical[key]
        if (
            current["ebv_core_p1_p9"] != legacy["ebv_core_p1_p9"]
            or current["self_core_p1_p9"] != legacy["self_core_p1_p9"]
        ):
            raise ValueError(f"exact duplicate register disagreement: {source['pair_id']}")
        current["source_membership"] = "v2_and_legacy_exact_duplicate"
        current["legacy_pair_id"] = legacy["legacy_pair_id"]
        current["legacy_original_score_coverage_status"] = legacy[
            "legacy_original_score_coverage_status"
        ]
        current["legacy_pair_validation"] = legacy["legacy_pair_validation"]
        current["ebv_legacy_candidate_id"] = legacy["ebv_candidate_id"]
        current["self_legacy_candidate_id"] = legacy["self_candidate_id"]
        overlaps.append({
            "v2_pair_id": current["v2_pair_id"],
            "legacy_pair_id": legacy["legacy_pair_id"],
            "ebv_sequence": key[0],
            "self_sequence": key[1],
            "ebv_core_p1_p9": current["ebv_core_p1_p9"],
            "self_core_p1_p9": current["self_core_p1_p9"],
            "registers_agree": True,
            "deduplication_action": "retain_v2_canonical_record_link_legacy_provenance",
        })

    ordered = sorted(
        canonical.values(),
        key=lambda row: (-float(row["tcr_facing_blosum62_similarity"]), str(row["pair_id"])),
    )
    count = len(ordered)
    position = 0
    while position < count:
        end = position + 1
        value = float(ordered[position]["tcr_facing_blosum62_similarity"])
        while end < count and float(ordered[end]["tcr_facing_blosum62_similarity"]) == value:
            end += 1
        for row in ordered[position:end]:
            row["combined_score_rank"] = position + 1
            row["combined_score_tie_size"] = end - position
        position = end
    for rank, row in enumerate(ordered, start=1):
        row["combined_rank"] = rank
        row["combined_pair_count"] = count
        row["combined_percentile"] = round((int(row["combined_score_rank"]) - 1) / max(1, count - 1), 8)
        old_rank = str(row.get("v2_original_hla_rank", ""))
        row["rank_shift_vs_v2"] = rank - int(old_rank) if old_rank else ""
        ebv_label = (
            f"{row['ebv_protein']} {row['ebv_start_1_based']}-{row['ebv_end_1_based']}"
            if row["ebv_start_1_based"] and row["ebv_end_1_based"] else row["ebv_protein"]
        )
        self_label = (
            f"{row['self_protein']} {row['self_start_1_based']}-{row['self_end_1_based']}"
            if row["self_start_1_based"] and row["self_end_1_based"] else row["self_protein"]
        )
        row["pair_coordinate_label"] = f"{ebv_label} / {self_label}*"
        row["primary_method"] = "tcr_facing_blosum62"
        row["primary_score"] = round(float(row["tcr_facing_blosum62_similarity"]), 12)
        row["rank_scope"] = "combined_unique_drb1501_pair_universe"
        row["computational_pair_marker"] = "*"
        row["pair_evidence_status"] = "computational_pair_no_exact_paired_recognition_evidence"
        row["claim_boundary"] = CLAIM_BOUNDARY
    return ordered, sorted(overlaps, key=lambda row: (row["v2_pair_id"], row["legacy_pair_id"]))


def _summary(rows: Sequence[Mapping[str, Any]], eligible_count: int, excluded_count: int) -> str:
    lines = [
        "# Combined HLA-DRB1*15:01 same-register ranking",
        "",
        f"Unique ranked pairs: **{len(rows)}**.",
        f"Legacy pairs admitted: **{eligible_count}**; excluded with reasons: **{excluded_count}**.",
        "",
        "The V2 and legacy universes were rescored together using the same control-selected P2/P3/P5/P7/P8 BLOSUM62 endpoint. Original rank numbers were not averaged.",
        "",
        "| Rank | Epitope pair | EBV core | Self core | Score | Source |",
        "|---:|---|---|---|---:|---|",
    ]
    for row in rows[:25]:
        lines.append(
            f"| {row['combined_rank']} | {row['pair_coordinate_label']} | `{row['ebv_core_p1_p9']}` | "
            f"`{row['self_core_p1_p9']}` | {float(row['primary_score']):.3f} | {row['source_membership']} |"
        )
    lines.extend([
        "",
        "\\* Computationally prioritized exact peptide pair; exact paired recognition has not been experimentally confirmed.",
        "",
        f"> {CLAIM_BOUNDARY}",
        "",
    ])
    return "\n".join(lines)


def run(
    *,
    v2_path: Path = DEFAULT_V2,
    legacy_path: Path = DEFAULT_LEGACY,
    benchmark_path: Path = DEFAULT_BENCHMARK,
    v2_registry_path: Path = DEFAULT_V2_REGISTRY,
    legacy_annotations_path: Path = DEFAULT_LEGACY_ANNOTATIONS,
    out: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    selected_method, control_evidence = select_control_supported_method(read_csv(benchmark_path))
    if selected_method != "tcr_facing_blosum62":
        raise ValueError(f"combined rank expected the locked control winner, got {selected_method}")
    legacy_all = read_csv(legacy_path)
    legacy_eligible, legacy_audit = audit_legacy_eligibility(legacy_all)
    ranked, overlaps = combine_and_rank(
        read_csv(v2_path),
        legacy_eligible,
        v2_registry=read_csv(v2_registry_path),
        legacy_annotations=read_csv(legacy_annotations_path),
    )
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "combined_ranked_pairs.csv", ranked)
    write_csv(out / "top_25_combined.csv", ranked[:25], list(ranked[0]))
    write_csv(out / "legacy_eligibility_audit_636.csv", legacy_audit)
    write_csv(out / "exact_duplicate_crosswalk.csv", overlaps)
    write_csv(out / "control_method_comparison.csv", control_evidence)
    write_json(out / "ranking_basis.json", {
        "allele": ALLELE,
        "primary_method": selected_method,
        "primary_metric": "normalized BLOSUM62 similarity at P2/P3/P5/P7/P8",
        "metric_direction": "higher_is_better",
        "method_selection": "held_out_positive_control_systems_only",
        "control_system_count": 3,
        "control_panel_capture": "8_of_8_rank_1",
        "geometry_used": False,
        "binding_percentile_used": False,
        "original_rank_numbers_combined": False,
        "combination_rule": "rescore_all_eligible_pairs_once_then_rank_one_deduplicated_universe",
        "tie_rule": "shared_score_rank_and_lexical_pair_id_display_order",
        "universe_limitation": (
            "legacy_and_v2_candidate_selection_rules_differ;combined_percentiles_describe_"
            "this_expanded_frozen_library_not_a_uniform_biological_background"
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    })
    excluded_count = len(legacy_all) - len(legacy_eligible)
    manifest = {
        "analysis_version": "EBV_MS_DRB1501_COMBINED_SAME_REGISTER_RANKING_2026-08-27",
        "allele": ALLELE,
        "v2_input_pair_count": len(read_csv(v2_path)),
        "legacy_input_pair_count": len(legacy_all),
        "legacy_eligible_pair_count": len(legacy_eligible),
        "legacy_excluded_pair_count": excluded_count,
        "exact_duplicate_pair_count": len(overlaps),
        "combined_unique_ranked_pair_count": len(ranked),
        "source_membership_counts": dict(sorted(Counter(row["source_membership"] for row in ranked).items())),
        "selected_primary_method": selected_method,
        "definitive_validation_complete": False,
        "ranking_status": "provisional_three_system_control_supported",
        "combined_universe_is_heterogeneous": True,
        "input_sha256": {
            "v2_drb1501_ranked_pairs": sha256_file(v2_path),
            "legacy_register_aware_pairs": sha256_file(legacy_path),
            "held_out_control_method_ranks": sha256_file(benchmark_path),
            "v2_candidate_registry": sha256_file(v2_registry_path),
            "legacy_protein_coordinate_annotations": sha256_file(legacy_annotations_path),
        },
        "rank_1_pair_id": ranked[0]["pair_id"],
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "analysis_manifest.json", manifest)
    (out / "RESULTS_SUMMARY.md").write_text(
        _summary(ranked, len(legacy_eligible), excluded_count), encoding="utf-8"
    )
    (out / "README.md").write_text(
        """# Combined DRB1*15:01 same-register ranking

This additive package unifies the V2 DRB1*15:01 screen with every eligible pair from the earlier 636-pair register-aware universe. It does not average old and new rank numbers. Every admitted pair is rescored using the control-selected TCR-facing BLOSUM62 method, exact peptide-pair duplicates are collapsed, and one new rank is assigned across the combined universe.

- `combined_ranked_pairs.csv`: complete combined ranking
- `top_25_combined.csv`: compact review table
- `legacy_eligibility_audit_636.csv`: inclusion or exclusion reason for all old pairs
- `exact_duplicate_crosswalk.csv`: deduplicated old/new pair identities
- `ranking_basis.json`: machine-readable ranking rule

An asterisk marks a computationally prioritized pair without exact experimental paired-recognition evidence. The result remains provisional because the selected endpoint has three strict control systems, below the six-system definitive gate.

The legacy and V2 candidate-selection rules were not identical. Therefore, the combined rank is an expanded-library prioritization, and its percentile describes only this frozen 2,043-pair universe rather than a uniformly sampled biological background.
""" + "\n" + CLAIM_BOUNDARY + "\n",
        encoding="utf-8",
    )
    (out / "METHODS.md").write_text(
        """# Methods

## Ranking endpoint

Both input universes were rescored with normalized BLOSUM62 similarity at predicted HLA-II core positions P2/P3/P5/P7/P8. Higher values rank first. This endpoint was selected before the discovery merge because it placed every declared positive first in all eight panels of the three-system held-out control benchmark. Binding percentile, AlphaFold geometry, original local alignments, and original rank numbers do not enter the combined score.

## Legacy eligibility

An older pair enters only when both P1-P9 cores contain exactly nine residues, the EBV core is a primary HLA-DRB1*15:01 prediction, and the self core is either a primary HLA-DRB1*15:01 prediction or the experimental primary-allele reference. DRB5 calibration rows and unresolved or sensitivity-only registers remain in the eligibility audit but cannot enter this DRB1*15:01 rank.

## Merge and ties

Exact full-peptide duplicates are retained once after proving that both predicted cores agree; V2 is the canonical record and the legacy identifier remains linked. Every remaining unique pair is ranked once. Equal primary scores share `combined_score_rank`; lexical pair ID provides only a deterministic display order in `combined_rank`.

## Limitation

The legacy and V2 libraries were assembled under different candidate-selection rules. The resulting percentile is descriptive of this frozen expanded universe, not a population probability, false-discovery rate, or uniformly sampled biological background.

""" + CLAIM_BOUNDARY + "\n",
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
    parser.add_argument("--v2", type=Path, default=DEFAULT_V2)
    parser.add_argument("--legacy", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK)
    parser.add_argument("--v2-registry", type=Path, default=DEFAULT_V2_REGISTRY)
    parser.add_argument("--legacy-annotations", type=Path, default=DEFAULT_LEGACY_ANNOTATIONS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(run(
        v2_path=args.v2,
        legacy_path=args.legacy,
        benchmark_path=args.benchmark,
        v2_registry_path=args.v2_registry,
        legacy_annotations_path=args.legacy_annotations,
        out=args.out,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
