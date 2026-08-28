"""Tests for register-filtered sequence and chemistry scoring."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_register_aware_score_table import build_score_rows, rendered_readme  # noqa: E402
from register_aware_scoring import score_same_register_alignment  # noqa: E402


class RegisterAwareScoringTests(unittest.TestCase):
    def test_keeps_only_positions_at_the_same_p1_p9_register_index(self):
        row = {
            "register_assessment": "assessable_register_hypothesis",
            "ebv_peptide": "MAVILKDFGQS",
            "human_peptide": "MAVILKDFGQS",
            "ebv_top_core_start_1_based": "3",
            "human_top_core_start_1_based": "3",
            "original_local_alignment_coordinates": (
                "2A:2A;3V:3V;4I:4I;5L:5L;6K:6K;7D:7D;"
                "8F:8F;9G:9G;10Q:10Q;11S:11S"
            ),
        }

        score = score_same_register_alignment(row)

        self.assertEqual(score["same_register_alignment_count"], 9)
        self.assertEqual(
            score["all_same_register_positions"], "P1;P2;P3;P4;P5;P6;P7;P8;P9"
        )
        self.assertEqual(score["anchor_same_register_positions"], "P1;P4;P6;P9")
        self.assertEqual(
            score["candidate_exposed_same_register_positions"], "P2;P3;P5;P7;P8"
        )
        self.assertEqual(score["score_coverage_status"], "robust_primary_ranking_eligible")

    def test_reports_component_similarity_for_a_known_charged_mismatch(self):
        row = {
            "register_assessment": "assessable_register_hypothesis",
            "ebv_peptide": "MKE",
            "human_peptide": "MEE",
            "ebv_top_core_start_1_based": "2",
            "human_top_core_start_1_based": "2",
            "original_local_alignment_coordinates": "2K:2E;3E:3E",
        }

        score = score_same_register_alignment(row)

        self.assertEqual(score["same_register_positions"], "P1;P2")
        self.assertEqual(score["same_register_identity_count"], 1)
        self.assertEqual(score["same_register_charge_similarity"], 0.5)
        self.assertLess(score["same_register_property_similarity"], 1.0)
        self.assertEqual(score["score_coverage_status"], "limited_coverage_report_only")

    def test_excludes_nonprimary_registers_before_feature_scoring(self):
        score = score_same_register_alignment(
            {
                "register_assessment": "calibration_only_nonprimary_allele",
                "ebv_peptide": "MAV",
                "human_peptide": "MAV",
                "ebv_top_core_start_1_based": "1",
                "human_top_core_start_1_based": "1",
                "original_local_alignment_coordinates": "1M:1M;2A:2A;3V:3V",
            }
        )

        self.assertEqual(score["score_coverage_status"], "excluded_nonprimary_register_status")
        self.assertEqual(score["same_register_alignment_count"], 0)

    def test_joins_geometry_as_context_without_making_it_a_primary_score(self):
        universe = [{
            "pair_id": "E::H",
            "ebv_candidate_id": "E",
            "human_candidate_id": "H",
            "register_assessment": "assessable_register_hypothesis",
            "ebv_peptide": "MAV",
            "human_peptide": "MAV",
            "ebv_top_core_start_1_based": "1",
            "human_top_core_start_1_based": "1",
            "original_local_alignment_coordinates": "1M:1M;2A:2A;3V:3V",
        }]
        geometry = [{
            "ebv_candidate_id": "E",
            "human_candidate_id": "H",
            "local_peptide_ca_rmsd_after_hla_fit": "1.234",
        }]

        rows = build_score_rows(universe, geometry)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["pair_id"], "E::H")
        self.assertEqual(rows[0]["whole_local_alignment_geometry_rmsd_context"], "1.234")
        self.assertEqual(rows[0]["geometry_context_status"], "whole_original_local_alignment_only")
        self.assertFalse(rows[0]["geometry_primary_score_eligible"])

    def test_readme_states_the_strict_decoy_decision_rule(self):
        readme = rendered_readme()

        self.assertIn("complete frozen strict decoys", readme)
        self.assertIn("negative/mixed", readme)


if __name__ == "__main__":
    unittest.main()
