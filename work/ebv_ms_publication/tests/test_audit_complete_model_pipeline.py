import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import audit_complete_model_pipeline as audit_pipeline


request_details = audit_pipeline.request_details
structural_consistency_class = audit_pipeline.structural_consistency_class


class RequestDetailsTests(unittest.TestCase):
    def test_discovery_roots_include_the_new_seed03_completion_folder(self):
        self.assertTrue(hasattr(audit_pipeline, "af3_discovery_roots"))
        project_root = Path("/project-without-local-model-downloads")

        try:
            roots = audit_pipeline.af3_discovery_roots(project_root)
        except FileNotFoundError as error:
            self.fail(f"Discovery roots must remain declarative when model downloads are absent: {error}")

        self.assertIn(
            (
                "new_background_af3_seed03_completion",
                project_root / "folds_2026_08_15_01_53",
                "*",
                3,
            ),
            roots,
        )

    def test_non_pmhc_request_is_inventoried_without_being_parsed_as_three_chain_pmhc(self):
        request = [{
            "name": "DECOY_02_HY_ENGA_DRB1_S101",
            "modelSeeds": ["101"],
            "sequences": [
                {"proteinChain": {"sequence": sequence}}
                for sequence in ("AAA", "BBB", "CCC", "DDD", "EEE")
            ],
        }]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fold_decoy_job_request.json"
            path.write_text(json.dumps(request), encoding="utf-8")

            try:
                observed = request_details(path)
            except ValueError as error:
                self.fail(f"Non-pMHC requests must remain inventory records: {error}")

        self.assertEqual(observed["candidate_id"], "DECOY_02_HY_ENGA_DRB1_S101")
        self.assertEqual(observed["protein_chain_count"], 5)
        self.assertEqual(observed["requested_peptide"], "")

    def test_structural_consistency_class_separates_robust_mixed_and_unstable_results(self):
        self.assertEqual(structural_consistency_class(0.90, 0.10, 0.758), "tier_A_robust")
        self.assertEqual(structural_consistency_class(0.60, 0.40, 1.208), "tier_B_mixed")
        self.assertEqual(structural_consistency_class(0.26, 0.60, 8.809), "tier_C_unstable_partial_pose")
        self.assertEqual(structural_consistency_class(0.00, 1.00, 5.630), "tier_D_no_consistent_pose")

    def test_controlled_summary_compares_target_with_unique_background_candidates(self):
        self.assertTrue(hasattr(audit_pipeline, "summarize_controlled_comparison"))
        target = [
            {"candidate_exposed_ca_rmsd_A": 1.0},
            {"candidate_exposed_ca_rmsd_A": 2.0},
        ]
        background = [
            {"background_candidate_id": "D1", "candidate_exposed_ca_rmsd_A": 3.0},
            {"background_candidate_id": "D1", "candidate_exposed_ca_rmsd_A": 4.0},
            {"background_candidate_id": "D2", "candidate_exposed_ca_rmsd_A": 5.0},
        ]

        summary = audit_pipeline.summarize_controlled_comparison(target, background)

        self.assertEqual(summary["controlled_comparison_status"], "evaluable_descriptive_matched_background")
        self.assertEqual(summary["structural_background_comparator_count"], 2)
        self.assertEqual(summary["structural_background_geometry_count"], 3)
        self.assertEqual(summary["target_candidate_exposed_rmsd_median_A"], 1.5)
        self.assertEqual(summary["background_candidate_exposed_rmsd_median_A"], 4.25)
        self.assertEqual(summary["background_minus_target_exposed_rmsd_median_A"], 2.75)
        self.assertEqual(summary["controlled_comparison_p_value"], "")


if __name__ == "__main__":
    unittest.main()
