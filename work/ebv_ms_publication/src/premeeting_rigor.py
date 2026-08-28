"""Deterministic helpers for pre-meeting pMHC rigor artifacts.

These utilities organize computational hypotheses and matching covariates. They
do not infer biological presentation, TCR recognition, or molecular mimicry.
"""

from __future__ import annotations

import csv
import io
from collections import Counter
from typing import Any


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def iedb_mhcii_eligible(peptide: str) -> bool:
    """Return whether a peptide fits the IEDB MHC-II API's 11--30 aa range."""
    return 11 <= len(peptide) <= 30


def iedb_submission_segments(candidate: dict[str, str]) -> list[dict[str, object]]:
    """Split long peptides into overlapping 30-mers that cover every 9-mer core."""
    peptide = candidate["peptide"]
    if len(peptide) < 11:
        return []
    if len(peptide) <= 30:
        return [{
            "submission_id": f"{candidate['candidate_id']}__segment_001",
            "candidate_id": candidate["candidate_id"],
            "peptide": peptide,
            "source_start_1_based": 1,
            "submission_strategy": "direct_full_peptide",
        }]
    starts = set(range(0, len(peptide) - 30 + 1, 22))
    starts.add(len(peptide) - 30)
    return [
        {
            "submission_id": f"{candidate['candidate_id']}__segment_{ordinal:03d}",
            "candidate_id": candidate["candidate_id"],
            "peptide": peptide[start:start + 30],
            "source_start_1_based": start + 1,
            "submission_strategy": "overlapping_30mer_scan",
        }
        for ordinal, start in enumerate(sorted(starts), start=1)
    ]


def natural_flank_submission_segment(
    candidate: dict[str, str], flank: dict[str, str]
) -> dict[str, object]:
    """Build an IEDB-safe segment from a source-coordinate-verified flank window."""
    sequence = flank["extended_sequence"]
    start = int(flank["original_start_in_extended_1_based"])
    end = int(flank["original_end_in_extended_1_based"])
    if not iedb_mhcii_eligible(sequence):
        raise ValueError("Natural-flank extension must be 11--30 residues for IEDB MHC-II")
    if sequence[start - 1:end] != candidate["peptide"]:
        raise ValueError("Natural-flank extension does not contain the manifest peptide at recorded positions")
    return {
        "submission_id": f"{candidate['candidate_id']}__natural_flank",
        "candidate_id": candidate["candidate_id"],
        "peptide": sequence,
        "source_start_1_based": 1,
        "original_start_in_submission_1_based": start,
        "submission_strategy": "verified_natural_flank_extension",
    }


def enumerate_core_windows(peptide: str) -> list[dict[str, object]]:
    """Return every inclusive 9-mer window with one-based peptide positions."""
    return [
        {"start": index + 1, "end": index + 9, "core_peptide": peptide[index:index + 9]}
        for index in range(max(0, len(peptide) - 8))
    ]


def parse_iedb_mhcii_tsv(text: str) -> list[dict[str, str]]:
    """Parse a tabular IEDB MHC-II API response and reject non-tabular errors."""
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    required = {"seq_num", "core_peptide", "peptide", "rank"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("IEDB MHC-II response is not a complete prediction table")
    return rows


def binding_rank_bin(rank: float) -> str:
    """Use fixed, predeclared bins for matching; lower percentile rank is better."""
    if rank <= 2.0:
        return "strong"
    if rank <= 10.0:
        return "intermediate"
    return "weak"


def map_prediction_rows(
    candidates: list[dict[str, str]], response_rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    """Attach API prediction rows to manifest candidates with sequence checks."""
    if len(candidates) != len(response_rows):
        raise ValueError("IEDB response row count does not match manifest candidate count")
    response_by_sequence_number: dict[int, dict[str, str]] = {}
    for response in response_rows:
        sequence_number = int(response["seq_num"])
        if sequence_number in response_by_sequence_number:
            raise ValueError(f"IEDB response contains duplicate seq_num {sequence_number}")
        response_by_sequence_number[sequence_number] = response
    mapped = []
    for index, candidate in enumerate(candidates, start=1):
        response = response_by_sequence_number.get(index)
        if response is None:
            raise ValueError(f"IEDB response is missing seq_num {index}")
        if candidate["peptide"] != response["peptide"]:
            raise ValueError(
                f"IEDB response peptide {response['peptide']} does not match candidate peptide {candidate['peptide']}"
            )
        mapped.append({**candidate, **response})
    return mapped


def composition_distance(left: str, right: str) -> float:
    """Return normalized L1 residue-composition distance on a 0--2 scale."""
    if not left or not right:
        raise ValueError("Peptides must be non-empty")
    left_counts, right_counts = Counter(left), Counter(right)
    return sum(
        abs(left_counts[residue] / len(left) - right_counts[residue] / len(right))
        for residue in AMINO_ACIDS
    )


def _numeric(record: dict[str, Any], field: str) -> float:
    return float(record[field])


def ordered_decoys(
    target: dict[str, Any], candidates: list[dict[str, Any]], limit: int
) -> list[dict[str, object]]:
    """Order candidate decoys only by predeclared nuisance-variable matching.

    The returned records deliberately omit screen priority fields. Selection uses
    paired peptide lengths, residue composition, peptide-model confidence, and
    predicted-binding rank bins. It never inspects a pMHC similarity/rank score.
    """
    if limit < 1:
        raise ValueError("limit must be positive")

    target_ebv_bin = binding_rank_bin(_numeric(target, "ebv_binding_rank"))
    target_human_bin = binding_rank_bin(_numeric(target, "human_binding_rank"))
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        if candidate["pair_id"] == target["pair_id"]:
            continue
        ebv_length_distance = abs(len(str(target["ebv_peptide"])) - len(str(candidate["ebv_peptide"])))
        human_length_distance = abs(len(str(target["human_peptide"])) - len(str(candidate["human_peptide"])))
        ebv_composition = composition_distance(str(target["ebv_peptide"]), str(candidate["ebv_peptide"]))
        human_composition = composition_distance(str(target["human_peptide"]), str(candidate["human_peptide"]))
        confidence_distance = abs(_numeric(target, "ebv_plddt") - _numeric(candidate, "ebv_plddt"))
        confidence_distance += abs(_numeric(target, "human_plddt") - _numeric(candidate, "human_plddt"))
        candidate_ebv_bin = binding_rank_bin(_numeric(candidate, "ebv_binding_rank"))
        candidate_human_bin = binding_rank_bin(_numeric(candidate, "human_binding_rank"))
        binding_bin_mismatches = int(target_ebv_bin != candidate_ebv_bin) + int(target_human_bin != candidate_human_bin)
        rows.append({
            "pair_id": str(candidate["pair_id"]),
            "ebv_length_distance": ebv_length_distance,
            "human_length_distance": human_length_distance,
            "total_length_distance": ebv_length_distance + human_length_distance,
            "ebv_composition_distance": round(ebv_composition, 6),
            "human_composition_distance": round(human_composition, 6),
            "total_composition_distance": round(ebv_composition + human_composition, 6),
            "model_confidence_distance": round(confidence_distance, 6),
            "target_ebv_binding_rank_bin": target_ebv_bin,
            "target_human_binding_rank_bin": target_human_bin,
            "candidate_ebv_binding_rank_bin": candidate_ebv_bin,
            "candidate_human_binding_rank_bin": candidate_human_bin,
            "binding_rank_bin_mismatches": binding_bin_mismatches,
            "meets_length_tolerance": ebv_length_distance <= 1 and human_length_distance <= 1,
        })
    rows.sort(key=lambda row: (
        int(row["binding_rank_bin_mismatches"]),
        int(row["total_length_distance"]),
        float(row["total_composition_distance"]),
        float(row["model_confidence_distance"]),
        str(row["pair_id"]),
    ))
    return rows[:limit]


def eligible_decoys(
    target: dict[str, Any], candidates: list[dict[str, Any]], limit: int
) -> tuple[list[dict[str, object]], int]:
    """Return only decoys meeting the predeclared length and binding-bin rules."""
    ranked = ordered_decoys(target, candidates, limit=len(candidates))
    eligible = [
        row
        for row in ranked
        if bool(row["meets_length_tolerance"]) and int(row["binding_rank_bin_mismatches"]) == 0
    ]
    return eligible[:limit], len(eligible)
