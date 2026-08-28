"""Tests for paper-ready register provenance reporting."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_register_aware_reporting_tables import (  # noqa: E402
    build_evidence_hierarchy_rows,
    build_sensitivity_rows,
    render_paper_result_branch,
)


class RegisterAwareReportingTests(unittest.TestCase):
    def test_keeps_same_allele_reference_and_nonprimary_calibration_distinct(self):
        overrides = [
            {
                "candidate_id": "HUMAN_MYELIN_13572",
                "analysis_role": "primary_experimental_reference",
                "presenting_allele": "HLA-DRB1*15:01",
                "core_peptide": "VHFFKNIVT",
                "core_start_1_based": "5",
                "register_source": "PDB 1BX2",
                "claim_boundary": "Same-allele pMHC register reference only.",
            },
            {
                "candidate_id": "EBV_TCELL_63843",
                "analysis_role": "calibration_only_nonprimary_allele",
                "presenting_allele": "HLA-DRB5*01:01",
                "core_peptide": "YHFVKKHVH",
                "core_start_1_based": "5",
                "register_source": "PDB 1H15",
                "claim_boundary": "Calibration only.",
            },
        ]

        rows = build_evidence_hierarchy_rows(overrides)

        self.assertEqual(rows[0]["primary_analysis_eligible"], True)
        self.assertEqual(rows[1]["primary_analysis_eligible"], False)
        self.assertEqual(rows[1]["eligible_use"], "cross_allotype_calibration_only")

    def test_reports_sensitivity_records_and_negative_method_branch(self):
        sensitivity = build_sensitivity_rows([{
            "pair_id": "gH::MBP",
            "register_assessment": "sensitivity_only_unresolved_register",
            "register_source": "IEDB hypothesis",
            "score_coverage_status": "excluded_nonprimary_register_status",
            "geometry_context_status": "whole_original_local_alignment_only",
        }])
        branch = render_paper_result_branch({
            "evaluable_target_count": 0,
            "global_inference_status": "insufficient_independent_systems",
        })

        self.assertEqual(sensitivity[0]["pair_id"], "gH::MBP")
        self.assertIn("negative_or_mixed_method_result", branch)
        self.assertIn("prospective", branch.lower())
        self.assertNotIn("shared-TCR", branch)

    def test_distinguishes_missing_same_register_alignment_from_sensitivity_only(self):
        rows = build_sensitivity_rows([{
            "pair_id": "EBV::myelin",
            "register_assessment": "no_same_register_local_alignment",
            "register_source": "IEDB hypothesis",
            "score_coverage_status": "excluded_nonprimary_register_status",
            "geometry_context_status": "whole_original_local_alignment_only",
        }])

        self.assertEqual(rows[0]["appendix_category"], "same_allele_hypothesis_no_matched_register")
        self.assertIn("No aligned residue", rows[0]["exclusion_reason"])


if __name__ == "__main__":
    unittest.main()
