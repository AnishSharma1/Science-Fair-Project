"""Create factual provenance tables and claim-bounded paper-result wording."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
OUT = PROC / "register_aware_scoring"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def build_evidence_hierarchy_rows(
    overrides: list[dict[str, Any]]
) -> list[dict[str, object]]:
    """Make reviewed experimental and exclusion decisions explicit."""
    use_by_role = {
        "primary_experimental_reference": (
            "same_allele_register_reference", True,
            "Same-allele pMHC register reference; not a TCR or cross-reactivity result.",
        ),
        "calibration_only_nonprimary_allele": (
            "cross_allotype_calibration_only", False,
            "DR15-haplotype calibration only; not same-allele benchmark evidence.",
        ),
        "sensitivity_only_unresolved": (
            "sensitivity_appendix_only", False,
            "Unresolved same-allele register; excluded from primary scoring.",
        ),
    }
    rows: list[dict[str, object]] = []
    for override in overrides:
        role = str(override["analysis_role"])
        if role not in use_by_role:
            raise ValueError(f"unknown override analysis role: {role}")
        eligible_use, primary_eligible, default_boundary = use_by_role[role]
        rows.append({
            "candidate_id": override["candidate_id"],
            "analysis_role": role,
            "presenting_allele": override.get("presenting_allele", ""),
            "core_peptide": override.get("core_peptide", ""),
            "core_start_1_based": override.get("core_start_1_based", ""),
            "register_source": override.get("register_source", ""),
            "eligible_use": eligible_use,
            "primary_analysis_eligible": primary_eligible,
            "claim_boundary": override.get("claim_boundary", default_boundary),
        })
    return rows


def build_sensitivity_rows(score_rows: list[dict[str, Any]]) -> list[dict[str, object]]:
    """Retain non-primary register outcomes instead of silently discarding them."""
    primary = "assessable_register_hypothesis"
    category_by_assessment = {
        "no_same_register_local_alignment": (
            "same_allele_hypothesis_no_matched_register",
            "No aligned residue occupied the same P1-P9 index under the retained same-allele register hypotheses.",
        ),
        "sensitivity_only_unresolved_register": (
            "unresolved_register_sensitivity_only",
            "An exact same-allele experimental P1-P9 register was not established; retained for sensitivity only.",
        ),
        "calibration_only_nonprimary_allele": (
            "cross_allotype_calibration_only",
            "Register support is from a nonprimary HLA allele and cannot enter the DRB1*15:01 primary analysis.",
        ),
        "unresolved_or_flank_dependent_core": (
            "unresolved_or_flank_dependent_register",
            "The retained predictor record did not provide a unique manifest-contained core for primary same-register comparison.",
        ),
    }
    appendix_rows: list[dict[str, object]] = []
    for row in score_rows:
        assessment = str(row.get("register_assessment", ""))
        if assessment == primary:
            continue
        category, reason = category_by_assessment.get(
            assessment,
            ("other_primary_exclusion", "Excluded from primary analysis; see the preserved register-assessment field."),
        )
        appendix_rows.append({
            "pair_id": row["pair_id"],
            "register_assessment": assessment,
            "register_source": row.get("register_source", ""),
            "ebv_register_status": row.get("ebv_register_status", ""),
            "human_register_status": row.get("human_register_status", ""),
            "score_coverage_status": row.get("score_coverage_status", ""),
            "geometry_context_status": row.get("geometry_context_status", ""),
            "appendix_category": category,
            "exclusion_reason": reason,
        })
    return appendix_rows


def render_paper_result_branch(decoy_summary: dict[str, Any]) -> str:
    """Render only the currently supported result branch from machine outputs."""
    evaluable = int(decoy_summary.get("evaluable_target_count", 0))
    if evaluable > 0 and decoy_summary.get("global_inference_status") == "predeclared_positive_gate_met":
        branch = "positive_computational_prioritization"
        body = (
            "The predefined score and strict-decoy conditions were met for the stated "
            "independent systems. This prioritizes pMHC candidate pairs for future testing "
            "under the defined computational model."
        )
    else:
        branch = "negative_or_mixed_method_result"
        body = (
            "The frozen inputs do not provide a robust evaluable strict-decoy comparison "
            "for the annotated targets. Therefore, the current evidence does not support "
            "robust register-aware computational prioritization beyond the reported "
            "hypothesis-generating candidates."
        )
    return "\n".join([
        f"# {branch}",
        "",
        body,
        "",
        "This result is limited to computational pMHC candidate prioritization and "
        "does not establish peptide presentation, T-cell activation, cross-reactivity, "
        "or an EBV-driven multiple-sclerosis mechanism.",
        "",
        "Prospective wet-lab work, outside this paper, could test final candidates with "
        "peptide-HLA binding/register assays and clone-defined functional assays.",
        "",
    ])


def main() -> None:
    overrides = read_csv(PROC / "register_sensitivity" / "experimental_register_overrides.csv")
    score_rows = read_csv(OUT / "register_aware_pair_scores.csv")
    decoy_summary = read_csv(OUT / "decoy_evaluation_summary.csv")
    if len(decoy_summary) != 1:
        raise ValueError("expected exactly one decoy evaluation summary row")
    write_csv(OUT / "evidence_hierarchy_table.csv", build_evidence_hierarchy_rows(overrides))
    write_csv(OUT / "register_sensitivity_appendix.csv", build_sensitivity_rows(score_rows))
    (OUT / "paper_result_branch.md").write_text(
        render_paper_result_branch(decoy_summary[0]), encoding="utf-8"
    )
    print("wrote evidence hierarchy, sensitivity appendix, and result branch")


if __name__ == "__main__":
    main()
