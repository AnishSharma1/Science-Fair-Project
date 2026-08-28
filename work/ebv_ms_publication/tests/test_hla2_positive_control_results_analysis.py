import json
import tempfile
import unittest
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from run_hla2_positive_control_results_analysis import (  # noqa: E402
    PACKAGE,
    _ligand_metadata,
    choose_endpoint,
    classify_panel,
    compare_request_to_expected,
    inventory_downloaded_jobs,
    sequence_identity,
    summarize_feature_values,
)


FEATURES = (
    "exposed_ca_rmsd_A",
    "exposed_sidechain_vector_rmsd_A",
    "tcr_face_physicochemical_mismatch",
    "anchor_ca_rmsd_A",
)


class ResultRequestTests(unittest.TestCase):
    def test_string_seed_is_accepted_only_after_exact_numeric_normalization(self):
        expected = {
            "name": "control_job",
            "modelSeeds": [104729],
            "sequences": [
                {"proteinChain": {"sequence": "ALPHA"}},
                {"proteinChain": {"sequence": "BETA"}},
                {"proteinChain": {"sequence": "PEPTIDE"}},
            ],
        }
        downloaded = json.loads(json.dumps([{**expected, "modelSeeds": ["104729"]}]))
        result = compare_request_to_expected(downloaded, expected)
        self.assertTrue(result["request_identity_pass"])
        self.assertEqual(result["normalized_seed"], 104729)
        self.assertEqual(result["seed_serialization"], "string")

    def test_request_identity_rejects_any_chain_sequence_change(self):
        expected = {
            "name": "control_job",
            "modelSeeds": [104729],
            "sequences": [
                {"proteinChain": {"sequence": "ALPHA"}},
                {"proteinChain": {"sequence": "BETA"}},
                {"proteinChain": {"sequence": "PEPTIDE"}},
            ],
        }
        downloaded = json.loads(json.dumps([expected]))
        downloaded[0]["sequences"][2]["proteinChain"]["sequence"] = "CHANGED"
        result = compare_request_to_expected(downloaded, expected)
        self.assertFalse(result["request_identity_pass"])
        self.assertFalse(result["chain_sequences_pass"])

    def test_inventory_requires_one_exact_five_model_bundle(self):
        expected = {
            "name": "control_job",
            "modelSeeds": [104729],
            "sequences": [
                {"proteinChain": {"sequence": "ALPHA"}},
                {"proteinChain": {"sequence": "BETA"}},
                {"proteinChain": {"sequence": "PEPTIDE"}},
            ],
        }
        manifest = [{"job_name": "control_job", "panel_seed": "104729"}]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = root / "folds_test" / "control_job"
            bundle.mkdir(parents=True)
            (bundle / "fold_control_job_job_request.json").write_text(
                json.dumps([{**expected, "modelSeeds": ["104729"]}]), encoding="utf-8"
            )
            for index in range(5):
                (bundle / f"fold_control_job_model_{index}.cif").write_text("model", encoding="utf-8")
                (bundle / f"fold_control_job_summary_confidences_{index}.json").write_text("{}", encoding="utf-8")
                (bundle / f"fold_control_job_full_data_{index}.json").write_text("{}", encoding="utf-8")
            rows = inventory_downloaded_jobs(manifest, [expected], [root])
        self.assertEqual(rows[0]["download_status"], "complete_exact")
        self.assertEqual(rows[0]["model_cif_count"], 5)
        self.assertEqual(rows[0]["complete_occurrence_count"], 1)


class ResultSummaryTests(unittest.TestCase):
    def test_reused_candidate_ids_are_scoped_to_their_hla_specific_pair(self):
        metadata = _ligand_metadata(PACKAGE)
        dr = metadata[("PAIR_OB1A12_ENGA_MBP", "HUMAN_BACKGROUND_186917")]
        dq = metadata[("PAIR_HY1B11_UL15_MBP", "HUMAN_BACKGROUND_186917")]
        self.assertEqual(dr["mhc_beta_allele"], "HLA-DRB1*15:01")
        self.assertEqual(dq["mhc_beta_allele"], "HLA-DQB1*05:02")

    def test_feature_summary_uses_median_and_interquartile_range(self):
        rows = [
            {feature: float(index) for feature in FEATURES}
            for index in (1, 2, 8, 9)
        ]
        summary = summarize_feature_values(rows, FEATURES)
        self.assertEqual(summary["exposed_ca_rmsd_A_median"], 5.0)
        self.assertEqual(summary["exposed_ca_rmsd_A_iqr"], 6.5)

    def test_incomplete_panel_keeps_available_rank_but_has_no_formal_rank(self):
        rows = []
        for index in range(17):
            rows.append({
                "pair_id": f"pair_{index:02d}",
                "pair_role": "positive" if index == 0 else "N3",
                "geometry_status": "complete",
                "exposed_ca_rmsd_A_median": float(index),
                "exposed_ca_rmsd_A_iqr": 0.0,
            })
        result = classify_panel(rows, expected_comparison_count=26)
        self.assertEqual(result["available_positive_rank"], 1)
        self.assertEqual(result["positive_rank"], "")
        self.assertEqual(result["evaluation_status"], "missing_required_comparisons")

    def test_complete_rank_above_three_is_a_completed_failure(self):
        rows = []
        for index in range(26):
            rows.append({
                "pair_id": f"pair_{index:02d}",
                "pair_role": "positive" if index == 4 else "N3",
                "geometry_status": "complete",
                "exposed_ca_rmsd_A_median": float(index),
                "exposed_ca_rmsd_A_iqr": 0.0,
            })
        result = classify_panel(rows, expected_comparison_count=26)
        self.assertEqual(result["positive_rank"], 5)
        self.assertFalse(result["capture_at_3"])
        self.assertEqual(result["evaluation_status"], "complete")

    def test_endpoint_choice_retains_exposed_ca_without_strict_worst_rank_improvement(self):
        self.assertEqual(choose_endpoint([1, 3], [1, 3])["selected_endpoint"], "frozen_exposed_ca")
        self.assertEqual(choose_endpoint([1, 2], [1, 3])["selected_endpoint"], "candidate_composite")
        self.assertEqual(choose_endpoint([1, 4], [1, 3])["selected_endpoint"], "frozen_exposed_ca")

    def test_sequence_identity_supports_exposed_and_full_core_positions(self):
        self.assertEqual(sequence_identity("ABCDEFGHI", "ABXDEYQHI"), 6 / 9)
        self.assertEqual(sequence_identity("ABCDEFGHI", "ABXDEYQHI", positions=(1, 2, 4, 6, 7)), 3 / 5)


if __name__ == "__main__":
    unittest.main()
