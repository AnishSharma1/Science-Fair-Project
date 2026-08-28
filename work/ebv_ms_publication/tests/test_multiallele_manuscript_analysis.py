"""Contract tests for the fixed multi-allele EBV--MS analysis."""

import random
import sys
import unittest
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from multiallele_manuscript_analysis import (  # noqa: E402
    build_prediction_submissions,
    build_pair_universe,
    build_robustness_jobs,
    canonical_request_name,
    direct_register_sequence_metrics,
    prediction_records_from_tsv,
    register_record,
    same_register_geometry_from_coordinates,
    select_score_blind_controls,
    select_representative_sample,
)


class MultiAlleleAnalysisTests(unittest.TestCase):
    def test_sequence_property_summary_collapses_descriptor_records_to_numbers(self):
        result = direct_register_sequence_metrics("ACDEFGHIK", "ACDEYGHIK")
        self.assertIsInstance(result["full_core_property_similarity_mean"], float)
        self.assertGreaterEqual(result["exposed_property_similarity_mean"], 0.0)
        self.assertLessEqual(result["exposed_property_similarity_mean"], 1.0)

    def test_prediction_submissions_preserve_50_records_and_use_verified_short_flank(self):
        panel = [
            {"candidate_id": f"P{i}", "peptide_sequence": "ACDEFGHIKLM"}
            for i in range(49)
        ] + [{"candidate_id": "SHORT", "peptide_sequence": "CDEFGHIKLM"}]
        flanks = {
            "SHORT": {
                "original_peptide": "CDEFGHIKLM",
                "extended_sequence": "AACDEFGHIKLMZZ",
                "original_start_in_extended_1_based": "3",
                "original_end_in_extended_1_based": "12",
            }
        }
        submissions = build_prediction_submissions(panel, flanks)
        self.assertEqual(len(submissions), 50)
        short = submissions[-1]
        self.assertEqual(short["prediction_input_peptide"], "AACDEFGHIKLMZZ")
        self.assertEqual(short["submission_strategy"], "verified_natural_flank_extension")
        self.assertEqual(short["original_start_in_prediction_1_based"], 3)

    def test_prediction_tsv_mapping_preserves_seq_num_and_register_status(self):
        submissions = [{
            "seq_num": 1,
            "candidate_id": "P1",
            "modeled_peptide": "ABCDEFGHIJKLM",
            "prediction_input_peptide": "ABCDEFGHIJKLM",
            "original_start_in_prediction_1_based": 1,
            "submission_strategy": "direct_full_peptide",
        }]
        raw = "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tic50\trank\nHLA-DRB1*13:03\t1\t1\t13\t13\tCDEFGHIJK\tABCDEFGHIJKLM\t25.0\t1.2\n"
        records = prediction_records_from_tsv("HLA-DRB1*13:03", submissions, raw, "raw.tsv")
        self.assertEqual(records[0]["iedb_seq_num"], 1)
        self.assertEqual(records[0]["register_status"], "resolved_unique_fully_contained")
        self.assertEqual(records[0]["predicted_ic50_nM"], "25.0")

    def test_canonical_request_name_only_collapses_download_duplicate_suffix(self):
        self.assertEqual(
            canonical_request_name("ebvms_drb1303_human_myelin_115622_2"),
            "ebvms_drb1303_human_myelin_115622",
        )
        self.assertEqual(
            canonical_request_name("ebvms_drb0301_ebv_tcell_950"),
            "ebvms_drb0301_ebv_tcell_950",
        )

    def test_representative_is_highest_ranked_clash_free_sample(self):
        rows = [
            {"sample_index": 0, "ranking_score": 0.95, "has_clash": True, "sequence_layout_status": "pass"},
            {"sample_index": 1, "ranking_score": 0.80, "has_clash": False, "sequence_layout_status": "pass"},
            {"sample_index": 2, "ranking_score": 0.85, "has_clash": False, "sequence_layout_status": "pass"},
        ]
        self.assertEqual(select_representative_sample(rows)["sample_index"], 2)
        self.assertIsNone(select_representative_sample([{**rows[0]}]))

    def test_register_requires_unique_fully_contained_core(self):
        direct = register_record(
            candidate_id="P1",
            modeled_peptide="ABCDEFGHIJKLM",
            prediction_input="ABCDEFGHIJKLM",
            predicted_core="CDEFGHIJK",
            percentile_rank=1.2,
            seq_num=7,
        )
        self.assertEqual(direct["register_status"], "resolved_unique_fully_contained")
        self.assertEqual(direct["core_start_1_based"], 3)

        flank = register_record(
            candidate_id="P2",
            modeled_peptide="CDEFGHIJKL",
            prediction_input="AACDEFGHIJKLZZ",
            predicted_core="AACDEFGHI",
            percentile_rank=4.0,
            seq_num=8,
            original_start_in_prediction_1_based=3,
        )
        self.assertEqual(flank["register_status"], "unresolved_flank_dependent")
        self.assertEqual(flank["core_start_1_based"], "")

        tied = register_record(
            candidate_id="P3",
            modeled_peptide="AAAAAAAAAA",
            prediction_input="AAAAAAAAAA",
            predicted_core="AAAAAAAAA",
            percentile_rank=20.0,
            seq_num=9,
        )
        self.assertEqual(tied["register_status"], "unresolved_tied_core_position")

    def test_control_selection_is_score_blind_and_deterministic(self):
        target = "LSRFSWGAEGQRPGFGYGG"
        rows = [
            {"candidate_id": "HUMAN_BACKGROUND_20", "iedb_epitope_id": "20", "peptide": "A" * 20, "peptide_length": 20, "binding_rank_bin": "weak", "geometry": 0.1},
            {"candidate_id": "HUMAN_BACKGROUND_10", "iedb_epitope_id": "10", "peptide": "LSRFSWGAEGQRPGFGYGGA", "peptide_length": 20, "binding_rank_bin": "weak", "geometry": 99.0},
            {"candidate_id": "HUMAN_BACKGROUND_30", "iedb_epitope_id": "30", "peptide": "C" * 20, "peptide_length": 20, "binding_rank_bin": "weak", "geometry": 1.0},
            {"candidate_id": "HUMAN_BACKGROUND_40", "iedb_epitope_id": "40", "peptide": "D" * 20, "peptide_length": 20, "binding_rank_bin": "strong", "geometry": 0.0},
        ]
        selected = select_score_blind_controls(target, "weak", rows, limit=3)
        shuffled = list(rows)
        random.Random(17).shuffle(shuffled)
        for row in shuffled:
            row["geometry"] = random.Random(row["iedb_epitope_id"]).random()
        repeated = select_score_blind_controls(target, "weak", shuffled, limit=3)
        self.assertEqual([r["candidate_id"] for r in selected], [r["candidate_id"] for r in repeated])
        self.assertEqual(selected[0]["candidate_id"], "HUMAN_BACKGROUND_10")
        self.assertTrue(all("geometry" not in row for row in selected))

    def test_pair_universe_is_625_per_allele_and_never_crosses_alleles(self):
        panel = [
            *({"candidate_id": f"E{i}", "arm_group": "EBV"} for i in range(25)),
            *({"candidate_id": f"H{i}", "arm_group": "CNS/self"} for i in range(25)),
        ]
        rows = build_pair_universe("HLA-DRB1*13:03", panel, {})
        self.assertEqual(len(rows), 625)
        self.assertEqual({row["allele"] for row in rows}, {"HLA-DRB1*13:03"})
        self.assertTrue(all(row["pair_id"].startswith("HLA-DRB1*13:03::") for row in rows))

    def test_geometry_is_rigid_body_invariant_and_uses_declared_positions(self):
        rng = np.random.default_rng(9)
        groove = rng.normal(size=(24, 3))
        core = rng.normal(size=(9, 3))
        theta = 0.71
        rotation = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        translation = np.array([4.0, -7.0, 2.5])
        moved_groove = groove @ rotation + translation
        moved_core = core @ rotation + translation
        result = same_register_geometry_from_coordinates(groove, core, moved_groove, moved_core)
        self.assertAlmostEqual(result["full_core_ca_rmsd_A"], 0.0, places=10)
        self.assertAlmostEqual(result["anchor_ca_rmsd_A"], 0.0, places=10)
        self.assertAlmostEqual(result["exposed_ca_rmsd_A"], 0.0, places=10)

    def test_robustness_jobs_use_two_fixed_seeds_and_maximum_30_jobs(self):
        entities = []
        for allele in ("HLA-DRB1*13:03", "HLA-DRB1*03:01", "HLA-DRB1*08:01"):
            for index in range(5):
                entities.append({"allele": allele, "entity_id": f"P{index}", "peptide": "ACDEFGHIKLM", "dra_sequence": "DRA", "drb_sequence": "DRB"})
        jobs, manifest = build_robustness_jobs(entities, seeds=(104729, 104759))
        self.assertEqual(len(jobs), 30)
        self.assertEqual(len(manifest), 30)
        self.assertEqual({job["modelSeeds"][0] for job in jobs}, {104729, 104759})
        with self.assertRaises(ValueError):
            build_robustness_jobs(entities + [entities[0]], seeds=(104729, 104759))


if __name__ == "__main__":
    unittest.main()
