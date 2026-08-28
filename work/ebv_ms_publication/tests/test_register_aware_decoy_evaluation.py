"""Tests for score-blind strict-decoy evaluation."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from run_register_aware_decoy_evaluation import evaluate_targets  # noqa: E402


class RegisterAwareDecoyEvaluationTests(unittest.TestCase):
    def test_evaluates_only_a_complete_five_decoy_robust_set(self):
        scores = [
            {
                "pair_id": "target",
                "score_coverage_status": "robust_primary_ranking_eligible",
                "same_register_property_similarity": "0.90",
                "same_register_identity_fraction": "0.80",
                "anchor_same_register_property_similarity": "0.95",
                "candidate_exposed_same_register_property_similarity": "0.85",
            },
            *[
                {
                    "pair_id": f"decoy_{index}",
                    "score_coverage_status": "robust_primary_ranking_eligible",
                    "same_register_property_similarity": str(value),
                    "same_register_identity_fraction": "0.10",
                    "anchor_same_register_property_similarity": "0.20",
                    "candidate_exposed_same_register_property_similarity": "0.30",
                }
                for index, value in enumerate([0.1, 0.2, 0.3, 0.4, 0.5], start=1)
            ],
        ]
        feasibility = [{
            "target_pair_id": "target",
            "target_validation_label": "system_a",
            "target_register_assessment": "assessable_register_hypothesis",
            "selected_decoy_count": "5",
        }]
        decoys = [
            {"target_pair_id": "target", "decoy_pair_id": f"decoy_{index}", "decoy_ordinal": str(index)}
            for index in range(1, 6)
        ]

        rows, summary = evaluate_targets(scores, decoys, feasibility)

        self.assertEqual(rows[0]["evaluation_status"], "evaluable_descriptive_only")
        self.assertEqual(rows[0]["target_minus_decoy_median"], 0.6)
        self.assertEqual(rows[0]["target_rank_among_six"], 1)
        self.assertEqual(rows[0]["descriptive_within_set_rank_fraction"], 1 / 6)
        self.assertEqual(summary["global_inference_status"], "insufficient_independent_systems")

    def test_does_not_backfill_or_score_a_limited_coverage_decoy(self):
        scores = [
            {
                "pair_id": "target",
                "score_coverage_status": "robust_primary_ranking_eligible",
                "same_register_property_similarity": "0.90",
            },
            {
                "pair_id": "decoy_1",
                "score_coverage_status": "limited_coverage_report_only",
                "same_register_property_similarity": "0.80",
            },
        ]
        feasibility = [{
            "target_pair_id": "target",
            "target_validation_label": "system_a",
            "target_register_assessment": "assessable_register_hypothesis",
            "selected_decoy_count": "1",
        }]
        decoys = [{"target_pair_id": "target", "decoy_pair_id": "decoy_1", "decoy_ordinal": "1"}]

        rows, _ = evaluate_targets(scores, decoys, feasibility)

        self.assertEqual(rows[0]["evaluation_status"], "not_evaluable_incomplete_strict_decoy_set")
        self.assertEqual(rows[0]["target_score"], "")


if __name__ == "__main__":
    unittest.main()
