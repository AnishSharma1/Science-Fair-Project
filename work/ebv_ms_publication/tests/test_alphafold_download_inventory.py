"""Tests for non-destructive AlphaFold Server download inventorying."""

import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_alphafold_download_inventory import inventory_expected_seed_jobs  # noqa: E402


class AlphaFoldDownloadInventoryTests(unittest.TestCase):
    def test_marks_absent_expected_jobs_as_no_result_without_retry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            folder = root / "json3folds"
            folder.mkdir()
            manifest = [{
                "batch_file": "alphafold_server_seed_03_jobs.json",
                "seed_label": "03",
                "model_seed": "104761",
                "candidate_id": "HUMAN_BACKGROUND_101",
                "job_name": "ebvms_bg_HUMAN_BACKGROUND_101_s03",
            }]

            rows = inventory_expected_seed_jobs(folder, manifest)

        self.assertEqual(rows[0]["completeness_status"], "not_downloaded_or_failed_no_retry")
        self.assertEqual(rows[0]["availability_interpretation"], "No model result; exclude from structural summaries without retrying.")

    def test_marks_a_job_complete_only_with_all_five_output_triplets(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary) / "json3folds"
            job = folder / "ebvms_bg_human_background_101_s03"
            job.mkdir(parents=True)
            prefix = "fold_ebvms_bg_human_background_101_s03"
            (job / f"{prefix}_job_request.json").touch()
            for index in range(5):
                (job / f"{prefix}_model_{index}.cif").touch()
                (job / f"{prefix}_summary_confidences_{index}.json").touch()
                (job / f"{prefix}_full_data_{index}.json").touch()
            manifest = [{
                "batch_file": "alphafold_server_seed_03_jobs.json",
                "seed_label": "03",
                "model_seed": "104761",
                "candidate_id": "HUMAN_BACKGROUND_101",
                "job_name": "ebvms_bg_HUMAN_BACKGROUND_101_s03",
            }]

            rows = inventory_expected_seed_jobs(folder, manifest)

        self.assertEqual(rows[0]["completeness_status"], "complete_five_sample_result")
        self.assertEqual(rows[0]["model_cif_count"], 5)


if __name__ == "__main__":
    unittest.main()
