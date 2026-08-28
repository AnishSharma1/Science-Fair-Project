"""Tests for direct P1--P9 AF3 pMHC comparison of predeclared pairs."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from same_register_af3_analysis import (  # noqa: E402
    direct_register_sequence_metrics,
    evaluate_sequence_decoys,
    matched_background_feasibility,
    same_register_geometry,
    summarize_job_pair_geometry,
)


def residue(x, y, z):
    return {"aa": "A", "atoms": [{"name": "CA", "element": "C", "xyz": (x, y, z)}], "bfactors": [90.0]}


class SameRegisterAf3AnalysisTests(unittest.TestCase):
    def test_summarizes_seed_job_pair_sensitivity_without_treating_samples_as_replicates(self):
        rows = [
            {"pair_id": "P", "shortlist_rank": 1, "ebv_job_directory": "E", "human_job_directory": "H", "core_p1_p9_ca_rmsd_A": 1.0, "candidate_exposed_ca_rmsd_A": 2.0},
            {"pair_id": "P", "shortlist_rank": 1, "ebv_job_directory": "E", "human_job_directory": "H", "core_p1_p9_ca_rmsd_A": 3.0, "candidate_exposed_ca_rmsd_A": 4.0},
        ]

        summary = summarize_job_pair_geometry(rows)

        self.assertEqual(summary[0]["cross_sample_comparison_count"], 2)
        self.assertEqual(summary[0]["core_p1_p9_ca_rmsd_A_median"], 2.0)
        self.assertEqual(summary[0]["candidate_exposed_ca_rmsd_A_median"], 3.0)

    def test_sequence_decoy_evaluation_ranks_target_without_computing_a_p_value(self):
        result = evaluate_sequence_decoys(
            "TARGET",
            {"candidate_exposed_property_similarity": 0.9},
            [
                {"candidate_id": "D1", "candidate_exposed_property_similarity": 0.7},
                {"candidate_id": "D2", "candidate_exposed_property_similarity": 0.8},
            ],
        )

        self.assertEqual(result["target_rank_among_target_plus_decoys"], 1)
        self.assertAlmostEqual(result["target_minus_decoy_median"], 0.15)
        self.assertEqual(result["p_value"], "")

    def test_scores_all_equivalent_register_positions_without_local_alignment_filter(self):
        metrics = direct_register_sequence_metrics("ACDEFGHIK", "ACDEYGHIK")

        self.assertEqual(metrics["same_register_position_count"], 9)
        self.assertEqual(metrics["anchor_identity_count"], 4)
        self.assertEqual(metrics["candidate_exposed_identity_count"], 4)
        self.assertEqual(metrics["candidate_exposed_positions"], "P2;P3;P5;P7;P8")

    def test_geometry_fits_hla_then_compares_matching_p1_to_p9_coordinates(self):
        hla = [residue(0, 0, 0), residue(1, 0, 0), residue(0, 1, 0)]
        peptide = [residue(float(index), 2, 0) for index in range(9)]
        reference = {"A": hla, "B": [residue(0, 0, 1), residue(1, 0, 1), residue(0, 1, 1)], "C": peptide}
        translated_hla = [residue(item["atoms"][0]["xyz"][0] + 10, item["atoms"][0]["xyz"][1], item["atoms"][0]["xyz"][2]) for item in hla]
        translated_beta = [residue(10, 0, 1), residue(11, 0, 1), residue(10, 1, 1)]
        shifted_peptide = [residue(float(index) + 10, 2, 2) for index in range(9)]
        other = {"A": translated_hla, "B": translated_beta, "C": shifted_peptide}

        metrics = same_register_geometry(reference, other, 1, 1)

        self.assertAlmostEqual(metrics["core_p1_p9_ca_rmsd_A"], 2.0, places=6)
        self.assertAlmostEqual(metrics["candidate_exposed_ca_rmsd_A"], 2.0, places=6)

    def test_matched_backgrounds_require_length_bin_and_completed_model(self):
        backgrounds = [
            {"candidate_id": "A", "peptide": "ACDEFGHIKLM", "binding_rank_bin": "weak"},
            {"candidate_id": "B", "peptide": "ACDEFGHIKLMN", "binding_rank_bin": "strong"},
            {"candidate_id": "C", "peptide": "ACDEFGHIKLMNPQ", "binding_rank_bin": "weak"},
        ]

        result = matched_background_feasibility(12, "weak", backgrounds, {"A"})

        self.assertEqual(result["planned_matched_background_count"], 1)
        self.assertEqual(result["completed_matched_background_count"], 1)
        self.assertEqual(result["completed_matched_background_ids"], "A")


if __name__ == "__main__":
    unittest.main()
