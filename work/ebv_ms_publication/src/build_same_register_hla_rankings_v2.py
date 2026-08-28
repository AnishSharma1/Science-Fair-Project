"""Build additive, HLA-specific same-register rankings for T-cell Library V2.

The primary nonstructural method is selected only from the completed held-out
positive-control benchmark. Discovery rows are read after method selection.
Existing structural rankings and benchmark packages are never modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from hla2_positive_control_benchmark import physicochemical_mismatch
from hla2_positive_control_benchmark_v2 import TCR_FACING_INDICES, blosum62_similarity


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DISCOVERY = (
    ROOT / "processed/tcell_library_v2_model_analysis_2026-08-25/discovery/pair_summary_6400.csv"
)
DEFAULT_BENCHMARK = (
    ROOT
    / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26/benchmark/method_rank_long.csv"
)
DEFAULT_PANEL = ROOT / "processed/tcell_library_v2_2026-08-22/frozen_v2_80_peptide_panel.csv"
DEFAULT_OUT = ROOT / "processed/tcell_library_v2_same_register_hla_rankings_2026-08-27"

ALLELES = (
    "HLA-DRB1*15:01",
    "HLA-DRB1*13:03",
    "HLA-DRB1*03:01",
    "HLA-DRB1*08:01",
)
ALLELE_SLUGS = {
    "HLA-DRB1*15:01": "hla_drb1_15_01",
    "HLA-DRB1*13:03": "hla_drb1_13_03",
    "HLA-DRB1*03:01": "hla_drb1_03_01",
    "HLA-DRB1*08:01": "hla_drb1_08_01",
}
NONSTRUCTURAL_METHODS = (
    "physicochemical_only",
    "tcr_facing_identity",
    "full_core_identity",
    "tcr_facing_blosum62",
    "full_core_blosum62",
)
METHOD_SPECS = {
    "physicochemical_only": ("tcr_face_physicochemical_mismatch", False),
    "tcr_facing_identity": ("tcr_facing_sequence_identity", True),
    "full_core_identity": ("full_core_sequence_identity", True),
    "tcr_facing_blosum62": ("tcr_facing_blosum62_similarity", True),
    "full_core_blosum62": ("full_core_blosum62_similarity", True),
}
CLAIM_BOUNDARY = (
    "Descriptive same-register pMHC sequence prioritization only; not evidence of "
    "presentation, TCR binding, activation, cross-reactivity, molecular mimicry, "
    "MS mechanism, probability, or false-discovery rate."
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str] | None = None,
) -> None:
    fieldnames = list(fields or (list(rows[0]) if rows else []))
    if not fieldnames:
        raise ValueError(f"cannot write empty CSV without fields: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _system_objective(system_ranks: Sequence[int], method: str) -> tuple[Any, ...]:
    if not system_ranks:
        return (math.inf, math.inf, math.inf, method)
    return (
        -sum(rank <= 3 for rank in system_ranks),
        max(system_ranks),
        -sum(1.0 / rank for rank in system_ranks) / len(system_ranks),
        method,
    )


def select_control_supported_method(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[str, list[dict[str, Any]]]:
    """Select a fixed sequence method using one conservative vote per system."""
    observed_methods = {str(row.get("method", "")) for row in rows}
    missing = set(NONSTRUCTURAL_METHODS) - observed_methods
    if missing:
        raise ValueError(f"missing nonstructural methods: {sorted(missing)}")

    expected_panels = None
    evidence = []
    choices = []
    for method in NONSTRUCTURAL_METHODS:
        method_rows = [row for row in rows if str(row["method"]) == method]
        panel_keys = {
            (
                str(row["system_id"]),
                str(row["positive_pair_id"]),
                int(row["panel_seed"]),
            )
            for row in method_rows
        }
        if len(panel_keys) != len(method_rows):
            raise ValueError(f"duplicate control panel rows for method {method}")
        if expected_panels is None:
            expected_panels = panel_keys
        elif panel_keys != expected_panels:
            raise ValueError("nonstructural methods do not cover identical control panels")
        if any(int(row["comparison_count"]) != 26 for row in method_rows):
            raise ValueError(f"method {method} includes an incomplete 26-comparison panel")

        by_system: dict[str, list[int]] = defaultdict(list)
        for row in method_rows:
            by_system[str(row["system_id"])].append(int(row["positive_rank"]))
        system_worst = {system: max(ranks) for system, ranks in sorted(by_system.items())}
        system_ranks = list(system_worst.values())
        row = {
            "method": method,
            "control_panel_count": len(method_rows),
            "independent_system_count": len(system_ranks),
            "system_capture_at_3_count": sum(rank <= 3 for rank in system_ranks),
            "all_panels_capture_at_3": all(int(item["positive_rank"]) <= 3 for item in method_rows),
            "worst_system_rank": max(system_ranks),
            "system_weighted_mrr": round(sum(1.0 / rank for rank in system_ranks) / len(system_ranks), 8),
            "system_worst_ranks": ";".join(f"{system}:{rank}" for system, rank in system_worst.items()),
            "selection_scope": "held_out_control_systems_only",
        }
        evidence.append(row)
        choices.append((_system_objective(system_ranks, method), method))

    selected = min(choices, key=lambda item: item[0])[1]
    for row in evidence:
        row["selected_primary_method"] = row["method"] == selected
    return selected, evidence


def _identity(left: str, right: str, indices: Sequence[int]) -> float:
    return sum(left[index] == right[index] for index in indices) / len(indices)


def sequence_metrics(left_core: str, right_core: str) -> dict[str, float]:
    if len(left_core) != 9 or len(right_core) != 9:
        raise ValueError("same-register ranking requires two exact nine-residue cores")
    return {
        "tcr_face_physicochemical_mismatch": physicochemical_mismatch(left_core, right_core),
        "tcr_facing_sequence_identity": _identity(left_core, right_core, TCR_FACING_INDICES),
        "full_core_sequence_identity": _identity(left_core, right_core, tuple(range(9))),
        "tcr_facing_blosum62_similarity": blosum62_similarity(
            left_core, right_core, positions=TCR_FACING_INDICES
        ),
        "full_core_blosum62_similarity": blosum62_similarity(left_core, right_core),
    }


def _ordered(rows: Sequence[Mapping[str, Any]], method: str) -> list[dict[str, Any]]:
    metric, higher_is_better = METHOD_SPECS[method]
    return sorted(
        (dict(row) for row in rows),
        key=lambda row: (
            -float(row[metric]) if higher_is_better else float(row[metric]),
            str(row["pair_id"]),
        ),
    )


def _score_rank_maps(
    rows: Sequence[Mapping[str, Any]], method: str
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    metric, _higher = METHOD_SPECS[method]
    ordered = _ordered(rows, method)
    display = {str(row["pair_id"]): index for index, row in enumerate(ordered, start=1)}
    score_rank: dict[str, int] = {}
    tie_size: dict[str, int] = {}
    position = 0
    while position < len(ordered):
        end = position + 1
        value = float(ordered[position][metric])
        while end < len(ordered) and float(ordered[end][metric]) == value:
            end += 1
        for row in ordered[position:end]:
            pair_id = str(row["pair_id"])
            score_rank[pair_id] = position + 1
            tie_size[pair_id] = end - position
        position = end
    return display, score_rank, tie_size


def rank_within_hla(
    pair_rows: Sequence[Mapping[str, Any]], selected_method: str
) -> dict[str, list[dict[str, Any]]]:
    if selected_method not in NONSTRUCTURAL_METHODS:
        raise ValueError(f"unsupported selected method: {selected_method}")
    if len({str(row.get("pair_id", "")) for row in pair_rows}) != len(pair_rows):
        raise ValueError("discovery pair IDs must be unique")

    by_allele: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for source in pair_rows:
        row = dict(source)
        allele = str(row.get("allele", ""))
        if allele not in ALLELES:
            raise ValueError(f"unexpected discovery allele: {allele}")
        if allele != str(row.get("ebv_allele", allele)) or allele != str(row.get("self_allele", allele)):
            raise ValueError(f"pair {row.get('pair_id')} is not an exact-HLA comparison")
        left = str(row.get("ebv_predicted_core", ""))
        right = str(row.get("self_predicted_core", ""))
        row.update(sequence_metrics(left, right))
        by_allele[allele].append(row)

    ranked: dict[str, list[dict[str, Any]]] = {}
    primary_metric = METHOD_SPECS[selected_method][0]
    for allele in ALLELES:
        rows = by_allele.get(allele, [])
        rank_maps = {
            method: _score_rank_maps(rows, method)
            for method in NONSTRUCTURAL_METHODS
        }
        ordered = _ordered(rows, selected_method)
        count = len(ordered)
        output = []
        for display_rank, row in enumerate(ordered, start=1):
            pair_id = str(row["pair_id"])
            primary_score_rank = rank_maps[selected_method][1][pair_id]
            primary_tie_size = rank_maps[selected_method][2][pair_id]
            output.append({
                "allele": allele,
                "hla_rank": display_rank,
                "hla_score_rank": primary_score_rank,
                "hla_score_tie_size": primary_tie_size,
                "hla_pair_count": count,
                "hla_percentile": round((primary_score_rank - 1) / max(1, count - 1), 8),
                "rank_scope": "within_hla_only",
                "pair_id": pair_id,
                "ebv_candidate_id": row["ebv_candidate_id"],
                "ebv_protein": row["ebv_protein"],
                "ebv_sequence": row["ebv_sequence"],
                "ebv_predicted_core": row["ebv_predicted_core"],
                "ebv_binding_percentile_rank": row["ebv_binding_percentile_rank"],
                "ebv_source_certainty": row["ebv_source_certainty"],
                "self_candidate_id": row["self_candidate_id"],
                "self_protein": row["self_protein"],
                "self_sequence": row["self_sequence"],
                "self_predicted_core": row["self_predicted_core"],
                "self_binding_percentile_rank": row["self_binding_percentile_rank"],
                "self_source_certainty": row["self_source_certainty"],
                "primary_method": selected_method,
                "primary_metric": primary_metric,
                "primary_score": round(float(row[primary_metric]), 12),
                "tcr_face_physicochemical_mismatch": round(float(row["tcr_face_physicochemical_mismatch"]), 12),
                "tcr_facing_sequence_identity": round(float(row["tcr_facing_sequence_identity"]), 12),
                "full_core_sequence_identity": round(float(row["full_core_sequence_identity"]), 12),
                "tcr_facing_blosum62_similarity": round(float(row["tcr_facing_blosum62_similarity"]), 12),
                "full_core_blosum62_similarity": round(float(row["full_core_blosum62_similarity"]), 12),
                **{
                    f"{method}_hla_rank": rank_maps[method][0][pair_id]
                    for method in NONSTRUCTURAL_METHODS
                },
                "register_comparison": "direct_P1_to_P9",
                "tcr_facing_positions": "P2;P3;P5;P7;P8",
                "ranking_input_status": "complete_sequence_register",
                "legacy_geometry_status": row.get("geometry_status", "not_available"),
                "geometry_used_in_ranking": False,
                "binding_percentile_used_in_ranking": False,
                "deterministic_tie_break": "lexical_pair_id",
                "ranking_validation_status": "provisional_three_system_control_supported",
                "pair_evidence_status": "computational_pair_no_exact_paired_recognition_evidence",
                "computational_pair_marker": "*",
                "claim_boundary": CLAIM_BOUNDARY,
            })
        ranked[allele] = output
    return ranked


def select_exact_target_epitopes(
    ranked: Mapping[str, Sequence[Mapping[str, Any]]],
    candidate_registry: Sequence[Mapping[str, Any]] = (),
) -> list[dict[str, Any]]:
    """Name exact EBNA1-ANO2 peptides only where that pair is HLA rank 1."""
    registry = {str(row.get("candidate_id", "")): row for row in candidate_registry}
    selected = []
    for allele in ALLELES:
        matches = [
            row for row in ranked.get(allele, ())
            if int(row["hla_rank"]) == 1
            and str(row["ebv_protein"]).upper() == "EBNA1"
            and str(row["self_protein"]).upper() == "ANO2"
        ]
        for row in matches:
            ebv_source = registry.get(str(row["ebv_candidate_id"]), {})
            self_source = registry.get(str(row["self_candidate_id"]), {})
            if ebv_source and str(ebv_source.get("sequence")) != str(row["ebv_sequence"]):
                raise ValueError("selected EBNA1 sequence does not match the frozen candidate registry")
            if self_source and str(self_source.get("sequence")) != str(row["self_sequence"]):
                raise ValueError("selected ANO2 sequence does not match the frozen candidate registry")
            selected.append({
                "allele": allele,
                "hla_rank": int(row["hla_rank"]),
                "pair_id": row["pair_id"],
                "ebv_candidate_id": row["ebv_candidate_id"],
                "ebv_epitope_15mer": row["ebv_sequence"],
                "ebv_core_p1_p9": row["ebv_predicted_core"],
                "ebv_source_accession": ebv_source.get("source_accession", ""),
                "ebv_source_start_1_based": ebv_source.get("source_start_1_based", ""),
                "ebv_source_end_1_based": ebv_source.get("source_end_1_based", ""),
                "ebv_evidence_status": ebv_source.get("source_certainty", row["ebv_source_certainty"]),
                "self_candidate_id": row["self_candidate_id"],
                "self_epitope_candidate_15mer": row["self_sequence"],
                "self_core_p1_p9": row["self_predicted_core"],
                "self_source_accession": self_source.get("source_accession", ""),
                "self_source_start_1_based": self_source.get("source_start_1_based", ""),
                "self_source_end_1_based": self_source.get("source_end_1_based", ""),
                "self_evidence_status": self_source.get("source_certainty", row["self_source_certainty"]),
                "tcr_facing_identity": row["tcr_facing_sequence_identity"],
                "tcr_facing_blosum62_similarity": row["tcr_facing_blosum62_similarity"],
                "selection_status": "selected_rank_1_within_hla",
                "prospective_not_experimentally_confirmed_pair": True,
                "cross_allele_consensus_used": False,
                "claim_boundary": CLAIM_BOUNDARY,
            })
    return selected


def build_exact_top_10_rows(
    ranked: Mapping[str, Sequence[Mapping[str, Any]]],
    candidate_registry: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Return presentation-ready exact peptide identities for each HLA top 10."""
    registry = {str(row.get("candidate_id", "")): row for row in candidate_registry}
    output = []
    for allele in ALLELES:
        for row in ranked.get(allele, ())[:10]:
            ebv_source = registry.get(str(row["ebv_candidate_id"]))
            self_source = registry.get(str(row["self_candidate_id"]))
            if not ebv_source or not self_source:
                raise ValueError(f"top-10 pair is missing frozen coordinate provenance: {row['pair_id']}")
            for label, source, sequence in (
                ("EBV", ebv_source, row["ebv_sequence"]),
                ("self", self_source, row["self_sequence"]),
            ):
                if str(source.get("sequence", "")) != str(sequence):
                    raise ValueError(f"{label} top-10 sequence does not match frozen candidate registry")
                if not str(source.get("source_start_1_based", "")) or not str(source.get("source_end_1_based", "")):
                    raise ValueError(f"{label} top-10 candidate lacks source coordinates")
            ebv_protein = str(row["ebv_protein"]).replace("_", "/")
            self_protein = str(row["self_protein"]).replace("_", "/")
            ebv_label = f"{ebv_protein} {ebv_source['source_start_1_based']}-{ebv_source['source_end_1_based']}"
            self_label = f"{self_protein} {self_source['source_start_1_based']}-{self_source['source_end_1_based']}"
            output.append({
                "allele": allele,
                "hla_rank": row["hla_rank"],
                "rank_display": str(row["hla_rank"]),
                "pair_id": row["pair_id"],
                "pair_coordinate_label": f"{ebv_label} / {self_label}*",
                "ebv_candidate_id": row["ebv_candidate_id"],
                "ebv_protein": ebv_protein,
                "ebv_source_start_1_based": ebv_source["source_start_1_based"],
                "ebv_source_end_1_based": ebv_source["source_end_1_based"],
                "ebv_epitope_label": ebv_label,
                "ebv_epitope_sequence": row["ebv_sequence"],
                "ebv_core_p1_p9": row["ebv_predicted_core"],
                "self_candidate_id": row["self_candidate_id"],
                "self_protein": self_protein,
                "self_source_start_1_based": self_source["source_start_1_based"],
                "self_source_end_1_based": self_source["source_end_1_based"],
                "self_epitope_label": self_label,
                "self_epitope_sequence": row["self_sequence"],
                "self_core_p1_p9": row["self_predicted_core"],
                "tcr_facing_blosum62_similarity": row["tcr_facing_blosum62_similarity"],
                "tcr_facing_sequence_identity": row["tcr_facing_sequence_identity"],
                "computational_pair_marker": "*",
                "pair_evidence_status": "computational_pair_no_exact_paired_recognition_evidence",
                "asterisk_definition": (
                    "computationally_prioritized_exact_peptide_pair_not_experimentally_confirmed"
                ),
                "claim_boundary": CLAIM_BOUNDARY,
            })
    return output


def _results_markdown(
    ranked: Mapping[str, Sequence[Mapping[str, Any]]],
    selected: str,
    evidence: Sequence[Mapping[str, Any]],
    selected_targets: Sequence[Mapping[str, Any]],
    selected_top_10: Sequence[Mapping[str, Any]],
) -> str:
    selected_row = next(row for row in evidence if row["method"] == selected)
    lines = [
        "# V2 same-register HLA-specific results",
        "",
        f"Primary method: **{selected}**, selected using held-out positive-control systems only.",
        "",
        f"Control result: **{selected_row['system_capture_at_3_count']}/{selected_row['independent_system_count']} systems captured at rank 3 or better**; worst system rank **{selected_row['worst_system_rank']}**.",
        "",
        "These are provisional sequence-based rankings. They do not satisfy the predeclared six-system definitive-validation requirement.",
        "",
    ]
    if selected_targets:
        first = selected_targets[0]
        alleles = ", ".join(str(row["allele"]) for row in selected_targets)
        lines.extend([
            "## Selected exact prospective epitope pair",
            "",
            f"- EBNA1 {first['ebv_source_start_1_based']}-{first['ebv_source_end_1_based']}: `{first['ebv_epitope_15mer']}`; predicted P1-P9 core `{first['ebv_core_p1_p9']}`.",
            f"- ANO2 {first['self_source_start_1_based']}-{first['self_source_end_1_based']}: `{first['self_epitope_candidate_15mer']}`; predicted P1-P9 core `{first['self_core_p1_p9']}`.",
            f"- Rank 1 independently in: **{alleles}**.",
            "- The EBNA1 peptide has exact positive T-cell evidence. The ANO2 sequence remains a region-derived candidate tile, not an experimentally confirmed epitope.",
            "",
        ])
    for allele in ALLELES:
        rows = list(ranked[allele])
        lines.extend([
            f"## {allele}",
            "",
            f"Ranked pairs: **{len(rows)}**.",
            "",
            "| Rank | Epitope pair | EBV core | Self core | P2/P3/P5/P7/P8 BLOSUM62 | Score rank | Tie size |",
            "|---:|---|---|---|---:|---:|---:|",
        ])
        exact_rows = [row for row in selected_top_10 if row["allele"] == allele]
        for row, exact in zip(rows[:10], exact_rows):
            lines.append(
                f"| {row['hla_rank']} | {exact['pair_coordinate_label']} | `{row['ebv_predicted_core']}` | "
                f"`{row['self_predicted_core']}` | "
                f"{float(row['tcr_facing_blosum62_similarity']):.3f} | "
                f"{row['hla_score_rank']} | {row['hla_score_tie_size']} |"
            )
        lines.extend([
            "",
            "\\* Computationally prioritized exact peptide pair; exact paired recognition has not been experimentally confirmed. Individual peptide evidence may differ.",
            "",
        ])
    lines.extend([f"> {CLAIM_BOUNDARY}", ""])
    return "\n".join(lines)


def run(
    *,
    discovery_path: Path = DEFAULT_DISCOVERY,
    benchmark_path: Path = DEFAULT_BENCHMARK,
    panel_path: Path = DEFAULT_PANEL,
    out: Path = DEFAULT_OUT,
) -> dict[str, Any]:
    benchmark_rows = read_csv(benchmark_path)
    selected, evidence = select_control_supported_method(benchmark_rows)
    discovery_rows = read_csv(discovery_path)
    ranked = rank_within_hla(discovery_rows, selected)
    panel_rows = read_csv(panel_path)
    selected_targets = select_exact_target_epitopes(ranked, panel_rows)
    exact_top_10 = build_exact_top_10_rows(ranked, panel_rows)
    all_ranked = [row for allele in ALLELES for row in ranked[allele]]
    fields = list(all_ranked[0]) if all_ranked else []
    if not fields:
        raise ValueError("same-register ranking requires at least one discovery pair")

    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "control_method_comparison.csv", evidence)
    write_csv(out / "all_hla_ranked_pairs.csv", all_ranked, fields)
    for allele in ALLELES:
        write_csv(out / f"rankings/{ALLELE_SLUGS[allele]}_ranked_pairs.csv", ranked[allele], fields)
    top_rows = [row for allele in ALLELES for row in ranked[allele][:25]]
    write_csv(out / "top_25_by_hla.csv", top_rows, fields)
    write_csv(out / "top_10_exact_epitopes_by_hla.csv", exact_top_10)
    selected_target_fields = list(selected_targets[0]) if selected_targets else (
        "allele", "hla_rank", "pair_id", "ebv_candidate_id", "ebv_epitope_15mer",
        "ebv_core_p1_p9", "ebv_source_accession", "ebv_source_start_1_based",
        "ebv_source_end_1_based", "ebv_evidence_status", "self_candidate_id",
        "self_epitope_candidate_15mer", "self_core_p1_p9", "self_source_accession",
        "self_source_start_1_based", "self_source_end_1_based", "self_evidence_status",
        "tcr_facing_identity", "tcr_facing_blosum62_similarity", "selection_status",
        "prospective_not_experimentally_confirmed_pair", "cross_allele_consensus_used",
        "claim_boundary",
    )
    write_csv(out / "selected_epitope_pair.csv", selected_targets, selected_target_fields)
    write_json(out / "selected_epitope_pair.json", {
        "selection_status": "selected_exact_epitope_pair" if selected_targets else "no_rank_1_ebna1_ano2_pair",
        "selected_hla_count": len(selected_targets),
        "selected_rows": selected_targets,
        "cross_allele_consensus_used": False,
        "claim_boundary": CLAIM_BOUNDARY,
    })

    selected_evidence = next(row for row in evidence if row["method"] == selected)
    systems = int(selected_evidence["independent_system_count"])
    gate = {
        "analysis_version": "EBV_MS_TCELL_V2_SAME_REGISTER_HLA_RANKINGS_2026-08-27",
        "control_support_status": "supportive_provisional",
        "selected_primary_method": selected,
        "selection_scope": "held_out_control_systems_only_before_discovery_read",
        "strict_independent_system_count": systems,
        "selected_method_system_capture_at_3_count": int(selected_evidence["system_capture_at_3_count"]),
        "selected_method_worst_system_rank": int(selected_evidence["worst_system_rank"]),
        "selected_method_system_weighted_mrr": float(selected_evidence["system_weighted_mrr"]),
        "selected_method_all_control_panels_capture_at_3": bool(selected_evidence["all_panels_capture_at_3"]),
        "minimum_systems_for_definitive_validation": 6,
        "definitive_validation_complete": systems >= 6,
        "definitive_status": "blocked_registry_size" if systems < 6 else "requires_full_gate_evaluation",
        "weights_frozen": False,
        "discovery_unlock_allowed": False,
        "exploratory_hla_specific_ranking_emitted": True,
        "cross_allele_consensus_used": False,
        "geometry_used_in_ranking": False,
        "binding_percentile_used_in_ranking": False,
        "claim_boundary": CLAIM_BOUNDARY,
    }
    write_json(out / "ranking_gate.json", gate)

    counts = {allele: len(ranked[allele]) for allele in ALLELES}
    legacy_missing = Counter(str(row.get("geometry_status", "not_available")) for row in discovery_rows)
    manifest = {
        **gate,
        "ranked_pair_count": len(all_ranked),
        "ranked_pairs_by_hla": counts,
        "legacy_geometry_status_counts": dict(sorted(legacy_missing.items())),
        "register_comparison": "direct_P1_to_P9",
        "tcr_facing_positions": ["P2", "P3", "P5", "P7", "P8"],
        "primary_metric_direction": "higher_is_better",
        "deterministic_tie_break": "lexical_pair_id",
        "input_sha256": {
            "discovery_pair_summary": sha256_file(discovery_path),
            "held_out_control_method_ranks": sha256_file(benchmark_path),
            "frozen_v2_peptide_panel": sha256_file(panel_path),
        },
        "selected_exact_epitope_pair_count_by_hla": len(selected_targets),
        "exact_top_10_row_count": len(exact_top_10),
        "all_top_10_pairs_marked_computational": all(
            row["computational_pair_marker"] == "*" for row in exact_top_10
        ),
        "rank_1_pair_by_hla": {
            allele: ranked[allele][0]["pair_id"] if ranked[allele] else None
            for allele in ALLELES
        },
    }
    write_json(out / "analysis_manifest.json", manifest)
    (out / "RESULTS_SUMMARY.md").write_text(
        _results_markdown(ranked, selected, evidence, selected_targets, exact_top_10),
        encoding="utf-8",
    )
    (out / "METHODS.md").write_text(
        """# Methods

## Control-locked method selection

Only the five nonstructural methods from the completed held-out HLA-II benchmark v2 were eligible. Each biological system contributed one conservative score: its worst required positive rank across ligands and seeds. Methods were selected by system capture at 3, then worst system rank, system-weighted reciprocal rank, and a fixed lexical tie-break. Discovery files were read only after this selection. The selected method was TCR-facing BLOSUM62 similarity.

## Same-register ranking

Each pair's two predicted nine-residue HLA-II binding cores were compared directly P1-to-P1 through P9-to-P9. The primary score uses only P2/P3/P5/P7/P8. Every HLA allele is a separate 40-by-40 EBV-self ranking. Higher primary similarity is better; exact score ties retain a shared score rank and tie size, while lexical pair ID supplies a reproducible display order. Geometry and binding percentile do not enter the ranking.

## Validation status

The selected method recovered every positive panel at rank 1 across three independent control systems, but this is below the predeclared minimum of six systems for definitive validation. Results are therefore provisional and do not overwrite earlier structural rankings.

## Interpretation

""" + CLAIM_BOUNDARY + "\n",
        encoding="utf-8",
    )
    (out / "README.md").write_text(
        """# T-cell Library V2 same-register rankings

This additive package ranks all V2 EBV-self pairs separately within each HLA using the control-selected P2/P3/P5/P7/P8 BLOSUM62 method. Existing structural packages remain unchanged.

- Full table: `all_hla_ranked_pairs.csv`
- Separate HLA tables: `rankings/`
- Compact view: `top_25_by_hla.csv`
- Exact top 10 for every HLA: `top_10_exact_epitopes_by_hla.csv`
- Exact selected candidate: `selected_epitope_pair.csv` and `selected_epitope_pair.json`
- Control comparison: `control_method_comparison.csv`
- Validation status: `ranking_gate.json`

The ranking is provisional because the method has three strict control systems rather than the six required for definitive validation. Score ties are explicit; `hla_rank` is a deterministic display order and `hla_score_rank` is the tied scientific rank. An asterisk marks a computationally prioritized exact peptide pair whose paired recognition has not been experimentally confirmed; it does not erase evidence attached to either individual peptide.

Reproduce with:

```bash
PYTHONPATH=src python3 src/build_same_register_hla_rankings_v2.py
```
""" + "\n" + CLAIM_BOUNDARY + "\n",
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
    parser.add_argument("--discovery", type=Path, default=DEFAULT_DISCOVERY)
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(run(discovery_path=args.discovery, benchmark_path=args.benchmark, panel_path=args.panel, out=args.out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
