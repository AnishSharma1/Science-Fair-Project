"""Tests for auditable IEDB prediction records for expanded comparators."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from run_expanded_background_iedb import build_prediction_rows  # noqa: E402


class ExpandedBackgroundIedbTests(unittest.TestCase):
    def test_records_the_best_mapped_prediction_without_claiming_presentation(self):
        candidates = [{
            "candidate_id": "HUMAN_BACKGROUND_101",
            "arm": "Human background comparator",
            "evidence_tier": "Tier 4",
            "peptide": "TGGVYHFVKKHVHES",
            "peptide_length": "15",
            "hla": "HLA-DRB1*15:01",
        }]
        mapped = [{
            "candidate_id": "HUMAN_BACKGROUND_101",
            "seq_num": "1",
            "peptide": "TGGVYHFVKKHVHES",
            "core_peptide": "YHFVKKHVH",
            "rank": "1.8",
            "ic50": "55.0",
            "submission_strategy": "direct_full_peptide",
        }]

        rows = build_prediction_rows(candidates, mapped, "raw.tsv", "2026-08-11T00:00:00Z")

        self.assertEqual(rows[0]["prediction_status"], "predicted")
        self.assertEqual(rows[0]["predicted_core_start_positions_1_based"], "5")
        self.assertEqual(rows[0]["binding_rank_bin"], "strong")
        self.assertIn("not experimental presentation", rows[0]["interpretation"])


if __name__ == "__main__":
    unittest.main()
