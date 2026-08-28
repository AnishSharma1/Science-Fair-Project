"""Tests for AlphaFold Server three-seed pMHC import batches."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import build_alphafold_server_seed_jobs as seed_jobs  # noqa: E402


build_server_jobs = seed_jobs.build_server_jobs


class AlphaFoldServerSeedJobTests(unittest.TestCase):
    def test_builds_server_dialect_jobs_with_one_predeclared_seed(self):
        rows = [{"id": "HUMAN_BACKGROUND_101", "sequence": "DRASEQ:DRBSEQ:ACDEFGHIKLM"}]

        jobs = build_server_jobs(rows, seed=104729, seed_label="01")

        self.assertEqual(jobs, [{
            "name": "ebvms_bg_HUMAN_BACKGROUND_101_s01",
            "modelSeeds": [104729],
            "sequences": [
                {"proteinChain": {"sequence": "DRASEQ", "count": 1}},
                {"proteinChain": {"sequence": "DRBSEQ", "count": 1}},
                {"proteinChain": {"sequence": "ACDEFGHIKLM", "count": 1}},
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }])

    def test_rejects_a_batch_row_that_is_not_exactly_three_chains(self):
        rows = [{"id": "HUMAN_BACKGROUND_101", "sequence": "DRASEQ:DRBSEQ"}]

        with self.assertRaisesRegex(ValueError, "exactly three chains"):
            build_server_jobs(rows, seed=104729, seed_label="01")

    def test_selects_only_jobs_without_complete_five_sample_results(self):
        self.assertTrue(hasattr(seed_jobs, "select_incomplete_jobs"))
        jobs = [
            {"name": "ebvms_bg_HUMAN_BACKGROUND_101_s03"},
            {"name": "ebvms_bg_HUMAN_BACKGROUND_102_s03"},
        ]
        inventory = [
            {
                "expected_job_name": "ebvms_bg_HUMAN_BACKGROUND_101_s03",
                "completeness_status": "complete_five_sample_result",
            },
            {
                "expected_job_name": "ebvms_bg_HUMAN_BACKGROUND_102_s03",
                "completeness_status": "not_downloaded_or_failed_no_retry",
            },
        ]

        selected = seed_jobs.select_incomplete_jobs(jobs, inventory)

        self.assertEqual(
            [job["name"] for job in selected],
            ["ebvms_bg_HUMAN_BACKGROUND_102_s03"],
        )


if __name__ == "__main__":
    unittest.main()
