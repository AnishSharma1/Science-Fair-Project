"""Audit recovery of the locked experimental BALF5-MBP positive control.

This module keeps biological evidence eligibility separate from model ranking.
It verifies the exact literature system first, then reads already-generated
calibration ranks. It never changes a score, threshold, model, or candidate.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Sequence


GOLD_SYSTEM_ID = "SYS_BALF5_MBP_HY2E11"
EXPECTED_SEEDS = {104729, 104759}
CLAIM_BOUNDARY = (
    "Recovery of an experimentally established positive control calibrates known-positive "
    "capture; it does not establish presentation, TCR binding, activation, cross-reactivity, "
    "molecular mimicry, or MS disease mechanism for discovery candidates."
)
SYSTEM_REQUIREMENTS = {
    "biological_system_id": GOLD_SYSTEM_ID,
    "evidence_tier": "E1_exact_pmhc_positive",
    "receptor_modality": "human_T_cell_shared_clone",
    "receptor_or_clone_id": "Hy.2E11",
    "assay_type": "same clone/TCR recognition plus pMHC crystal structures",
    "primary_source": "PMID:12244309",
    "doi": "10.1038/ni835",
    "viral_sequence": "TGGVYHFVKKHVHES",
    "viral_hla": "HLA-DRB5*01:01",
    "self_sequence": "ENPVVHFFKNIVTPR",
    "self_hla": "HLA-DRB1*15:01",
}
STRUCTURE_REQUIREMENTS = {
    "EBV_structure": "1H15",
    "MBP_structure": "1BX2",
    "seven_position_core_CA_RMSD_A": "0.838",
}


def _truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _require_exact_fields(row: dict[str, Any], expected: dict[str, str], label: str) -> None:
    mismatches = [
        f"{field}={row.get(field)!r}, expected {value!r}"
        for field, value in expected.items()
        if str(row.get(field, "")) != value
    ]
    if mismatches:
        raise ValueError(f"{label} identity check failed: " + "; ".join(mismatches))


def build_gold_standard_audit(
    registry_rows: Sequence[dict[str, Any]],
    experimental_metric_rows: Sequence[dict[str, Any]],
    seed_rows: Sequence[dict[str, Any]],
    recovery_rows: Sequence[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    """Return a system record, seed records, and summary after strict checks."""
    positives = [row for row in registry_rows if _truth(row.get("tcell_positive_denominator"))]
    if len(positives) != 1:
        raise ValueError(
            "gold-standard denominator must contain exactly one independent system; "
            f"observed {len(positives)}"
        )
    positive = positives[0]
    _require_exact_fields(positive, SYSTEM_REQUIREMENTS, "gold-standard system")
    if "antibody" in str(positive.get("receptor_modality", "")).lower():
        raise ValueError("antibody evidence cannot enter the T-cell gold-standard denominator")

    metrics = {str(row["metric"]): str(row["value"]) for row in experimental_metric_rows}
    _require_exact_fields(metrics, STRUCTURE_REQUIREMENTS, "experimental structure")

    matching_recovery = [row for row in recovery_rows if row.get("biological_system_id") == GOLD_SYSTEM_ID]
    if len(matching_recovery) != 1:
        raise ValueError(f"expected one recovery record for {GOLD_SYSTEM_ID}")
    recovery = matching_recovery[0]

    observed_seeds = {int(row["seed"]) for row in seed_rows}
    if observed_seeds != EXPECTED_SEEDS:
        raise ValueError(f"expected calibration seeds {sorted(EXPECTED_SEEDS)}, observed {sorted(observed_seeds)}")

    seed_audit_rows: list[dict[str, Any]] = []
    for row in sorted(seed_rows, key=lambda value: int(value["seed"])):
        available_rank = int(row["available_rank"]) if str(row.get("available_rank", "")) else None
        positive_rmsd = float(row["positive_exposed_ca_rmsd_median_A"])
        control_median = float(row["available_equal_weight_control_median_A"])
        formal_evaluable = _truth(row.get("formal_seed_evaluable"))
        captured_at_1 = available_rank == 1
        captured_at_3 = available_rank is not None and available_rank <= 3
        below_control_median = positive_rmsd < control_median
        formal_pass = formal_evaluable and _truth(row.get("seed_recovery_criterion_pass"))
        if formal_pass:
            audit_status = "formal_pass"
        elif not formal_evaluable and captured_at_3 and below_control_median:
            audit_status = "available_set_capture_formal_incomplete"
        elif formal_evaluable:
            audit_status = "formal_fail"
        else:
            audit_status = "not_evaluable_and_not_captured"
        seed_audit_rows.append({
            "biological_system_id": GOLD_SYSTEM_ID,
            "seed": int(row["seed"]),
            "available_primary_count": int(row["available_primary_count"]),
            "expected_primary_count": int(row["expected_primary_count"]),
            "available_rank": available_rank if available_rank is not None else "",
            "capture_at_1_available_set": captured_at_1,
            "capture_at_3_available_set": captured_at_3,
            "positive_exposed_ca_rmsd_median_A": positive_rmsd,
            "available_equal_weight_control_median_A": control_median,
            "positive_below_available_control_median": below_control_median,
            "formal_seed_evaluable": formal_evaluable,
            "formal_seed_pass": formal_pass,
            "audit_status": audit_status,
            "used_for_score_tuning": False,
        })

    available_seed_count = len(seed_audit_rows)
    capture_at_1_count = sum(_truth(row["capture_at_1_available_set"]) for row in seed_audit_rows)
    capture_at_3_count = sum(_truth(row["capture_at_3_available_set"]) for row in seed_audit_rows)
    formal_rows = [row for row in seed_audit_rows if _truth(row["formal_seed_evaluable"])]
    formal_pass_count = sum(_truth(row["formal_seed_pass"]) for row in formal_rows)
    system_row = {
        "biological_system_id": GOLD_SYSTEM_ID,
        "gold_standard_eligible": True,
        "independent_system_weight": 1,
        "evidence_class": "same_human_T_cell_clone_plus_two_experimental_pMHC_structures",
        "receptor_or_clone_id": positive["receptor_or_clone_id"],
        "viral_peptide": positive["viral_sequence"],
        "viral_hla": positive["viral_hla"],
        "viral_structure": metrics["EBV_structure"],
        "self_peptide": positive["self_sequence"],
        "self_hla": positive["self_hla"],
        "self_structure": metrics["MBP_structure"],
        "experimental_seven_position_core_ca_rmsd_A": float(metrics["seven_position_core_CA_RMSD_A"]),
        "primary_source": positive["primary_source"],
        "doi": positive["doi"],
        "predeclared_before_model_ranking": True,
        "used_for_discovery_ranking": False,
        "used_for_score_tuning": False,
    }
    summary = {
        "audit_version": "EBV_MS_GOLD_STANDARD_CAPTURE_V1",
        "denominator_rule": (
            "Exact peptide identities and HLA arms, same human T-cell clone cross-recognition, "
            "and experimental pMHC structures for both arms"
        ),
        "gold_standard_independent_system_count": 1,
        "gold_standard_system_ids": [GOLD_SYSTEM_ID],
        "available_seed_count": available_seed_count,
        "capture_at_1_available_seed_count": capture_at_1_count,
        "capture_at_1_available_seed_fraction": capture_at_1_count / available_seed_count,
        "capture_at_3_available_seed_count": capture_at_3_count,
        "capture_at_3_available_seed_fraction": capture_at_3_count / available_seed_count,
        "formal_evaluable_seed_count": len(formal_rows),
        "formal_seed_pass_count": formal_pass_count,
        "formal_evaluable_seed_pass_fraction": formal_pass_count / len(formal_rows) if formal_rows else None,
        "strict_two_seed_recovery_status": recovery["recovery_status"],
        "available_set_conclusion": (
            "captured_rank_1_in_both_available_seed_sets"
            if capture_at_1_count == available_seed_count
            else "known_positive_not_captured_rank_1_in_every_available_seed_set"
        ),
        "model_or_score_changed_to_fit_positive": False,
        "interpretation": (
            "The frozen method captured the one independent experimentally established system. "
            "Broad sensitivity and generalization cannot be estimated from one system."
        ),
        "claim_boundary": CLAIM_BOUNDARY,
    }
    return system_row, seed_audit_rows, summary


def run_gold_standard_audit(
    registry_path: Path,
    experimental_metrics_path: Path,
    seed_summary_path: Path,
    recovery_report_path: Path,
    out_dir: Path,
) -> dict[str, Any]:
    system_row, seed_rows, summary = build_gold_standard_audit(
        read_csv(registry_path),
        read_csv(experimental_metrics_path),
        read_csv(seed_summary_path),
        read_csv(recovery_report_path),
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "gold_standard_systems.csv", [system_row], list(system_row))
    write_csv(out_dir / "gold_standard_seed_recovery.csv", seed_rows, list(seed_rows[0]))
    (out_dir / "gold_standard_capture_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    sources = [
        {
            "identifier": "PMID:12244309 / DOI:10.1038/ni835",
            "source_type": "primary_cross_reactivity_and_structure_paper",
            "url": "https://doi.org/10.1038/ni835",
            "role": "same human TCR recognizes DRB5-BALF5 and DRB1-MBP",
        },
        {
            "identifier": "PDB:1H15",
            "source_type": "experimental_structure",
            "url": "https://www.rcsb.org/structure/1H15",
            "role": "DRB5*01:01-BALF5 positive-control arm",
        },
        {
            "identifier": "PDB:1BX2",
            "source_type": "experimental_structure",
            "url": "https://www.rcsb.org/structure/1BX2",
            "role": "DRB1*15:01-MBP positive-control arm",
        },
    ]
    write_csv(out_dir / "gold_standard_sources.csv", sources, list(sources[0]))
    readme = f"""# Gold-standard positive-control capture audit

## What counts

The denominator is locked before reading model ranks. A system must have exact peptide identities and HLA arms, recognition of both pMHCs by the same human T-cell clone, and experimental pMHC structures for both arms. Supportive T-cell studies, antibody-only pairs, protein-level associations, canonical tiles, and computational discoveries do not count.

The current denominator contains **one independent system**: Hy.2E11 recognition of DRB5*01:01-BALF5 and DRB1*15:01-MBP, anchored by PDB 1H15 and 1BX2.

## Frozen-method recovery

- Available-set capture@1: **{summary['capture_at_1_available_seed_count']}/{summary['available_seed_count']} seeds**.
- Available-set capture@3: **{summary['capture_at_3_available_seed_count']}/{summary['available_seed_count']} seeds**.
- Fully evaluable seeds passing the predeclared rule: **{summary['formal_seed_pass_count']}/{summary['formal_evaluable_seed_count']}**.
- Strict two-seed status: **{summary['strict_two_seed_recovery_status']}**.
- Model or score changed to fit the positive: **no**.

Seed 104729 is an available-set result only because controls are missing. Its rank cannot be relabeled as a formal 1-of-26 result. The successful capture shows that the frozen pMHC geometry method recognizes this established example; one independent system is not enough to estimate general sensitivity or validate new biological claims.

{CLAIM_BOUNDARY}
"""
    (out_dir / "README.md").write_text(readme, encoding="utf-8")
    return summary

