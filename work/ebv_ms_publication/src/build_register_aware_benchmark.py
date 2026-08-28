"""Build a complete-universe, register-aware pMHC decoy benchmark.

This builder performs computational bookkeeping only. Its outputs organize
pMHC candidate-prioritization hypotheses; they are not evidence of peptide
presentation, shared-TCR binding, activation, cross-reactivity, or an MS
mechanism.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from register_aware_benchmark import (
    is_assessable_same_register_pair,
    strict_eligible_decoys,
)
from register_aware_diagnostic import (
    parse_local_alignment_positions,
    same_register_alignment_count,
)


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "register_aware_benchmark"
PRIMARY_ALLELE = "HLA-DRB1*15:01"
REGISTER_SOURCE = "IEDB recommended_binding HLA-DRB1*15:01 top-core hypothesis"
STRICT_MATCHING_RULE = (
    "Both EBV and human peptide lengths within one residue and zero IEDB "
    "predicted-binding-rank-bin mismatches; composition distance and peptide "
    "pLDDT order only already-eligible decoys."
)
SELECTION_BOUNDARY = (
    "Decoys were selected from assessable background rows using only peptide "
    "length, amino-acid composition, peptide pLDDT, and IEDB binding-rank bins; "
    "never a pMHC priority score, rank, or structural similarity value."
)


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a UTF-8 CSV into dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    """Write rows with an explicit field order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def unique_contained_core_start(prediction: dict[str, str]) -> int | None:
    """Return a unique manifest-contained top-core start, otherwise ``None``."""
    if prediction.get("predicted_core_fully_contained_in_manifest_peptide") != "True":
        return None
    starts = [
        value
        for value in prediction.get("predicted_core_start_positions_1_based", "").split(";")
        if value
    ]
    if len(starts) != 1:
        return None
    return int(starts[0])


def resolve_candidate_register(
    candidate_id: str,
    peptide: str,
    prediction: dict[str, str],
    override: dict[str, str] | None = None,
) -> dict[str, object]:
    """Resolve a candidate through reviewed evidence before an IEDB fallback.

    An override never mutates a stored predictor result. It records the source
    hierarchy used for this benchmark and makes nonprimary or unresolved
    records ineligible for same-allele comparison.
    """
    if override is not None:
        if override.get("candidate_id") != candidate_id:
            raise ValueError(f"Override candidate mismatch for {candidate_id}")
        role = override["analysis_role"]
        source = override["register_source"]
        core = override.get("core_peptide", "")
        start_text = override.get("core_start_1_based", "")
        if role == "sensitivity_only_unresolved":
            if core or start_text:
                raise ValueError(f"Sensitivity-only override must not set a core for {candidate_id}")
            return {
                "register_status": "sensitivity_only_unresolved",
                "register_source": source,
                "core_start_1_based": None,
                "core_peptide": "",
            }
        if role not in {
            "primary_experimental_reference",
            "calibration_only_nonprimary_allele",
        }:
            raise ValueError(f"Unknown override role {role} for {candidate_id}")
        if len(core) != 9 or not start_text:
            raise ValueError(f"Override must supply a 9-mer core and start for {candidate_id}")
        start = int(start_text)
        if peptide[start - 1:start + 8] != core:
            raise ValueError(f"Override core does not match peptide for {candidate_id}")
        if role == "primary_experimental_reference":
            if override["presenting_allele"] != PRIMARY_ALLELE:
                raise ValueError(f"Primary reference has wrong allele for {candidate_id}")
            status = "experimental_primary_allele_reference"
        else:
            if override["presenting_allele"] == PRIMARY_ALLELE:
                raise ValueError(f"Calibration override must be nonprimary for {candidate_id}")
            status = "calibration_only_nonprimary_allele"
        return {
            "register_status": status,
            "register_source": source,
            "core_start_1_based": start,
            "core_peptide": core,
        }

    start = unique_contained_core_start(prediction)
    return {
        "register_status": "iedb_top_core_hypothesis" if start is not None else "unresolved_or_flank_dependent_core",
        "register_source": REGISTER_SOURCE,
        "core_start_1_based": start,
        "core_peptide": prediction.get("predicted_core_peptide", ""),
    }


def load_register_overrides() -> dict[str, dict[str, str]]:
    """Load reviewed register decisions and reject duplicate candidate rows."""
    overrides: dict[str, dict[str, str]] = {}
    for row in read_csv(PROC / "register_sensitivity" / "experimental_register_overrides.csv"):
        candidate_id = row["candidate_id"]
        if candidate_id in overrides:
            raise ValueError(f"Duplicate register override for {candidate_id}")
        overrides[candidate_id] = row
    return overrides


def classify_pair_validation(
    ebv_candidate_id: str,
    human_candidate_id: str,
    validation_groups_by_candidate: dict[str, set[str]],
) -> str:
    """Apply the project’s existing source/context overlay labels exactly."""
    ebv_groups = validation_groups_by_candidate.get(ebv_candidate_id, set())
    human_groups = validation_groups_by_candidate.get(human_candidate_id, set())
    groups = ebv_groups | human_groups
    classic = "classic_BALF5_MBP_structural_positive"
    drosu = "drosu_2024_DRB1501_EBV_glycoprotein"
    wang = "wang_2026_MBP90_region"
    if ebv_candidate_id == "EBV_TCELL_63843" and classic in human_groups:
        return "classic_BALF5_MBP_pair"
    if drosu in groups and wang in groups:
        return "combined_new_literature_overlay"
    if drosu in groups:
        return "drosu_2024_EBV_glycoprotein"
    if wang in groups:
        return "wang_2026_MBP90_region"
    if classic in groups:
        return "classic_component_only"
    return "background"


def build_pair_universe(
    geometry_rows: list[dict[str, str]],
    prediction_by_candidate: dict[str, dict[str, str]],
    manifest_by_candidate: dict[str, dict[str, str]],
    validation_groups_by_candidate: dict[str, set[str]] | None = None,
    overrides_by_candidate: dict[str, dict[str, str]] | None = None,
) -> list[dict[str, object]]:
    """Return every PASS geometry pair with an auditable register decision."""
    validation_groups_by_candidate = validation_groups_by_candidate or {}
    overrides_by_candidate = overrides_by_candidate or {}
    rows: list[dict[str, object]] = []
    for geometry_row in geometry_rows:
        if geometry_row.get("status") != "PASS":
            continue
        ebv_id = geometry_row["ebv_candidate_id"]
        human_id = geometry_row["human_candidate_id"]
        ebv_prediction = prediction_by_candidate.get(ebv_id)
        human_prediction = prediction_by_candidate.get(human_id)
        if ebv_prediction is None or human_prediction is None:
            raise ValueError(f"Missing register prediction for {ebv_id}::{human_id}")
        ebv_manifest = manifest_by_candidate.get(ebv_id)
        human_manifest = manifest_by_candidate.get(human_id)
        if ebv_manifest is None or human_manifest is None:
            raise ValueError(f"Missing manifest record for {ebv_id}::{human_id}")

        ebv_register = resolve_candidate_register(
            ebv_id,
            ebv_manifest["peptide"],
            ebv_prediction,
            overrides_by_candidate.get(ebv_id),
        )
        human_register = resolve_candidate_register(
            human_id,
            human_manifest["peptide"],
            human_prediction,
            overrides_by_candidate.get(human_id),
        )
        ebv_start = ebv_register["core_start_1_based"]
        human_start = human_register["core_start_1_based"]
        alignment = parse_local_alignment_positions(
            geometry_row.get("aligned_positions_ebv_to_human", "")
        )
        same_count: int | str = ""
        statuses = {ebv_register["register_status"], human_register["register_status"]}
        if "sensitivity_only_unresolved" in statuses:
            assessment = "sensitivity_only_unresolved_register"
        elif "calibration_only_nonprimary_allele" in statuses:
            assessment = "calibration_only_nonprimary_allele"
        elif ebv_start is None or human_start is None:
            assessment = "unresolved_or_flank_dependent_core"
        else:
            if not isinstance(ebv_start, int) or not isinstance(human_start, int):
                raise ValueError(f"Resolved core start must be integer for {ebv_id}::{human_id}")
            same_count = same_register_alignment_count(alignment, ebv_start, human_start)
            assessment = (
                "assessable_register_hypothesis"
                if is_assessable_same_register_pair(alignment, ebv_start, human_start)
                else "no_same_register_local_alignment"
            )
        pair_validation = classify_pair_validation(
            ebv_id, human_id, validation_groups_by_candidate
        )
        rows.append({
            "pair_id": f"{ebv_id}::{human_id}",
            "ebv_candidate_id": ebv_id,
            "human_candidate_id": human_id,
            "pair_validation": pair_validation,
            "register_source": (
                f"EBV: {ebv_register['register_source']} | "
                f"Human: {human_register['register_source']}"
            ),
            "register_assessment": assessment,
            "ebv_register_status": ebv_register["register_status"],
            "human_register_status": human_register["register_status"],
            "ebv_register_source": ebv_register["register_source"],
            "human_register_source": human_register["register_source"],
            "ebv_top_core_start_1_based": ebv_start if ebv_start is not None else "",
            "human_top_core_start_1_based": human_start if human_start is not None else "",
            "ebv_top_core_peptide": ebv_register["core_peptide"],
            "human_top_core_peptide": human_register["core_peptide"],
            "same_register_alignment_count": same_count,
            "original_local_alignment_coordinates": geometry_row.get(
                "aligned_positions_ebv_to_human", ""
            ),
            "ebv_peptide": ebv_manifest["peptide"],
            "human_peptide": human_manifest["peptide"],
            "ebv_binding_rank": float(ebv_prediction["predicted_percentile_rank"]),
            "human_binding_rank": float(human_prediction["predicted_percentile_rank"]),
            "ebv_plddt": float(geometry_row["ebv_peptide_mean_plddt"]),
            "human_plddt": float(geometry_row["human_peptide_mean_plddt"]),
            "decoy_background_eligible": (
                pair_validation == "background"
                and assessment == "assessable_register_hypothesis"
            ),
            "interpretation": (
                "Computational register-aware pMHC bookkeeping only; not evidence of "
                "presentation, shared-TCR binding, activation, cross-reactivity, or MS mechanism."
            ),
        })
    return rows


def load_validation_groups() -> dict[str, set[str]]:
    """Map each candidate to its pre-existing literature/context group(s)."""
    groups: dict[str, set[str]] = {}
    for row in read_csv(PROC / "external_validation_panel.csv"):
        groups.setdefault(row["candidate_id"], set()).add(row["validation_group"])
    return groups


def build_decoy_benchmark(
    pair_rows: list[dict[str, Any]], target_decoy_count: int = 5
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Select strict decoys and retain all target feasibility decisions."""
    if target_decoy_count < 1:
        raise ValueError("target_decoy_count must be positive")
    background = [
        row
        for row in pair_rows
        if row["pair_validation"] == "background" and row["decoy_background_eligible"]
    ]
    targets = [row for row in pair_rows if row["pair_validation"] != "background"]
    decoy_rows: list[dict[str, object]] = []
    feasibility_rows: list[dict[str, object]] = []
    for target in targets:
        assessment = str(target["register_assessment"])
        base = {
            "target_pair_id": target["pair_id"],
            "target_validation_label": target["pair_validation"],
            "target_register_assessment": assessment,
            "target_decoy_count": target_decoy_count,
            "matching_rule": STRICT_MATCHING_RULE,
            "selection_boundary": SELECTION_BOUNDARY,
            "interpretation": (
                "Computational pMHC candidate-prioritization feasibility only; "
                "not evidence of presentation, shared-TCR binding, activation, "
                "cross-reactivity, or MS mechanism."
            ),
        }
        if assessment != "assessable_register_hypothesis":
            feasibility_rows.append({
                **base,
                "eligible_decoy_count": 0,
                "selected_decoy_count": 0,
                "readiness_status": "not assessable",
            })
            continue
        selected, available = strict_eligible_decoys(
            target, background, limit=target_decoy_count
        )
        readiness = (
            "ready"
            if available >= target_decoy_count
            else "partial"
            if available > 0
            else "no eligible decoy"
        )
        feasibility_rows.append({
            **base,
            "eligible_decoy_count": available,
            "selected_decoy_count": len(selected),
            "readiness_status": readiness,
        })
        for ordinal, selected_row in enumerate(selected, start=1):
            decoy_rows.append({
                "target_pair_id": target["pair_id"],
                "target_validation_label": target["pair_validation"],
                "decoy_ordinal": ordinal,
                "decoy_pair_id": selected_row["pair_id"],
                "ebv_length_distance": selected_row["ebv_length_distance"],
                "human_length_distance": selected_row["human_length_distance"],
                "total_length_distance": selected_row["total_length_distance"],
                "ebv_composition_distance": selected_row["ebv_composition_distance"],
                "human_composition_distance": selected_row["human_composition_distance"],
                "total_composition_distance": selected_row["total_composition_distance"],
                "model_confidence_distance": selected_row["model_confidence_distance"],
                "target_ebv_binding_rank_bin": selected_row["target_ebv_binding_rank_bin"],
                "target_human_binding_rank_bin": selected_row["target_human_binding_rank_bin"],
                "candidate_ebv_binding_rank_bin": selected_row["candidate_ebv_binding_rank_bin"],
                "candidate_human_binding_rank_bin": selected_row["candidate_human_binding_rank_bin"],
                "binding_rank_bin_mismatches": selected_row["binding_rank_bin_mismatches"],
                "meets_length_tolerance": selected_row["meets_length_tolerance"],
                "matching_rule": STRICT_MATCHING_RULE,
                "selection_boundary": SELECTION_BOUNDARY,
                "interpretation": (
                    "Strict matched decoy for computational pMHC candidate "
                    "prioritization; not a validated biological negative."
                ),
            })
    return decoy_rows, feasibility_rows


def summarize_benchmark(
    pair_rows: list[dict[str, object]], feasibility_rows: list[dict[str, object]]
) -> dict[str, int]:
    """Count assessability and strict-decoy readiness without scoring pairs."""
    return {
        "pass_pair_universe": len(pair_rows),
        "assessable_pairs": sum(
            row["register_assessment"] == "assessable_register_hypothesis"
            for row in pair_rows
        ),
        "eligible_background_pairs": sum(
            bool(row["decoy_background_eligible"]) for row in pair_rows
        ),
        "annotated_target_records": len(feasibility_rows),
        "ready_targets": sum(row["readiness_status"] == "ready" for row in feasibility_rows),
        "partial_targets": sum(row["readiness_status"] == "partial" for row in feasibility_rows),
        "no_eligible_decoy_targets": sum(
            row["readiness_status"] == "no eligible decoy" for row in feasibility_rows
        ),
        "not_assessable_targets": sum(
            row["readiness_status"] == "not assessable" for row in feasibility_rows
        ),
    }


def render_benchmark_readme(summary: dict[str, int]) -> str:
    """Render a concise, claim-safe report for the benchmark outputs."""
    return f"""# Register-aware matched-decoy benchmark

## What this artifact does

This is a complete-universe feasibility and pMHC candidate-prioritization
artifact. It does not establish shared-TCR binding, peptide presentation,
T-cell activation, cross-reactivity, or an MS mechanism.

## Frozen computational boundary

- Register hierarchy: `processed/register_sensitivity/experimental_register_overrides.csv`
  is applied before the stored IEDB top-core hypotheses.
- PDB 1BX2 makes `VHFFKNIVT` at positions 5--13 the exact MBP(85--99)
  DRB1*15:01 reference register. PDB 1H15 makes BALF5 `YHFVKKHVH` a
  DRB5*01:01 calibration-only register, not a primary-screen override.
- Exact gH `EKQLFYYIGTMLPN` and candidate-MBP `QRPGFGYGGRASDYKSAHK`
  records are sensitivity-only because no exact DRB1*15:01 experimental
  register was established.
- All candidates without a registry decision retain `{REGISTER_SOURCE}` as a
  computational hypothesis.
- A pair is assessable only when both resolved cores are unique and contained
  in the manifest peptide, neither arm is calibration-only or sensitivity-only,
  and at least one pre-existing local alignment maps to the same P1--P9 position.
- Decoy matching: {STRICT_MATCHING_RULE}
- Selection boundary: {SELECTION_BOUNDARY}

## Results counts

- pass_pair_universe: {summary.get('pass_pair_universe', 0)}
- assessable_pairs: {summary.get('assessable_pairs', 0)}
- eligible_background_pairs: {summary.get('eligible_background_pairs', 0)}
- annotated_target_records: {summary.get('annotated_target_records', 0)}
- ready_targets: {summary.get('ready_targets', 0)}
- partial_targets: {summary.get('partial_targets', 0)}
- no_eligible_decoy_targets: {summary.get('no_eligible_decoy_targets', 0)}
- not_assessable_targets: {summary.get('not_assessable_targets', 0)}

Annotated target records are source/context overlays, not independent
biological positives. In particular, BALF5--MBP overlapping records remain one
DR15-haplotype calibration system, not replicated validation systems.

## Interpretation

The output measures whether the current computational inputs can support the
predeclared strict comparison. It does not calculate a new pMHC effect size or
rerank candidates. A target labelled `ready` has five covariate-matched,
score-blind decoys available; this is a design-readiness property, not support
for molecular mimicry.

## Reproduce

```bash
PYTHONPATH=src python3 src/build_register_aware_benchmark.py
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py' -v
```

## Remaining primary-evidence gates

1. Establish allele-resolved experimental registers for the exact gH and
   candidate-MBP peptides before they enter a same-register benchmark.
2. Obtain direct DRB1*15:01 evidence for BALF5 before promoting it beyond
   DR15-haplotype calibration.
3. Use exact pMHC structures, not the register table alone, to assign any
   peptide-specific solvent exposure or candidate TCR-facing positions.
"""


def main() -> None:
    """Write the complete pair universe and strict matched-decoy outputs."""
    geometry_rows = read_csv(PROC / "colabfold_tier1_ebv_myelin_geometry_matrix.csv")
    prediction_by_candidate = {
        row["candidate_id"]: row
        for row in read_csv(PROC / "register_sensitivity" / "register_prediction_summary.csv")
    }
    manifest_by_candidate = {
        row["candidate_id"]: row
        for row in read_csv(PROC / "pmhc_candidate_manifest.csv")
    }
    rows = build_pair_universe(
        geometry_rows,
        prediction_by_candidate,
        manifest_by_candidate,
        load_validation_groups(),
        load_register_overrides(),
    )
    if not rows:
        raise ValueError("No PASS geometry rows available for the benchmark universe")
    write_csv(OUT / "benchmark_pair_universe.csv", rows, list(rows[0]))
    decoy_rows, feasibility_rows = build_decoy_benchmark(rows)
    feasibility_fields = [
        "target_pair_id", "target_validation_label", "target_register_assessment",
        "eligible_decoy_count", "selected_decoy_count", "target_decoy_count",
        "readiness_status", "matching_rule", "selection_boundary", "interpretation",
    ]
    decoy_fields = [
        "target_pair_id", "target_validation_label", "decoy_ordinal", "decoy_pair_id",
        "ebv_length_distance", "human_length_distance", "total_length_distance",
        "ebv_composition_distance", "human_composition_distance", "total_composition_distance",
        "model_confidence_distance", "target_ebv_binding_rank_bin",
        "target_human_binding_rank_bin", "candidate_ebv_binding_rank_bin",
        "candidate_human_binding_rank_bin", "binding_rank_bin_mismatches",
        "meets_length_tolerance", "matching_rule", "selection_boundary", "interpretation",
    ]
    write_csv(OUT / "target_feasibility.csv", feasibility_rows, feasibility_fields)
    write_csv(OUT / "matched_decoy_sets.csv", decoy_rows, decoy_fields)
    summary = summarize_benchmark(rows, feasibility_rows)
    (OUT / "README.md").write_text(render_benchmark_readme(summary), encoding="utf-8")
    print(f"Wrote {OUT / 'benchmark_pair_universe.csv'} ({len(rows)} PASS pairs)")
    print(f"Wrote {OUT / 'target_feasibility.csv'} ({len(feasibility_rows)} targets)")
    print(f"Wrote {OUT / 'matched_decoy_sets.csv'} ({len(decoy_rows)} strict decoys)")
    print(f"Wrote {OUT / 'README.md'}")


if __name__ == "__main__":
    main()
