import copy
import csv
import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from run_tcell_library_v2_model_analysis import (  # noqa: E402
    Occurrence,
    SampleGeometry,
    _bundle_fingerprint,
    choose_occurrence,
    classify_recovery,
    pair_geometry,
    rank_pair_rows,
    summarize_cross_allele,
    write_csv,
)


class V2ModelAnalysisTests(unittest.TestCase):
    def test_bundle_fingerprint_ignores_finder_renamed_prefixes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left, right = root / "left", root / "right"
            left.mkdir()
            right.mkdir()
            for role, content in (
                ("job_request.json", "request"),
                ("model_0.cif", "model"),
                ("summary_confidences_0.json", "summary"),
            ):
                (left / f"fold_original_{role}").write_text(content, encoding="utf-8")
                (right / f"fold_original_2_{role}").write_text(content, encoding="utf-8")
            self.assertEqual(_bundle_fingerprint(left), _bundle_fingerprint(right))

    def test_csv_writer_uses_union_for_mixed_inventory_schemas(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mixed.csv"
            write_csv(path, [{"shared": 1, "discovery": 2}, {"shared": 3, "calibration": 4}])
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(set(rows[0]), {"shared", "discovery", "calibration"})
        self.assertEqual(rows[1]["calibration"], "4")

    def test_occurrence_selection_is_score_blind_and_fingerprint_deterministic(self):
        rows = [
            Occurrence("job", Path("/z"), "bbb", 5, 5, 5, 1),
            Occurrence("job", Path("/a"), "aaa", 5, 5, 5, 1),
            Occurrence("job", Path("/incomplete"), "000", 4, 5, 5, 1),
        ]
        self.assertEqual(choose_occurrence(rows).directory, Path("/a"))

    def test_pair_geometry_is_rigid_body_invariant(self):
        rng = np.random.default_rng(21)
        groove = rng.normal(size=(170, 3))
        core = rng.normal(size=(9, 3))
        theta = 0.53
        rotation = np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ])
        translation = np.array([7.0, -2.0, 4.0])
        left = SampleGeometry("d", "l", "a", "e", "", 0, groove, core, core, core, core)
        right = SampleGeometry(
            "d", "r", "a", "s", "", 0,
            groove @ rotation + translation,
            core @ rotation + translation,
            core @ rotation + translation,
            core @ rotation + translation,
            core @ rotation + translation,
        )
        result = pair_geometry(left, right)
        for value in result.values():
            self.assertAlmostEqual(value, 0.0, places=9)

    def test_pair_ranking_is_within_allele_and_deterministic(self):
        rows = [
            {"allele": "A", "pair_id": "p2", "geometry_status": "complete", "exposed_ca_rmsd_A_median": 1.0, "exposed_ca_rmsd_A_iqr": 0.2},
            {"allele": "A", "pair_id": "p1", "geometry_status": "complete", "exposed_ca_rmsd_A_median": 1.0, "exposed_ca_rmsd_A_iqr": 0.1},
            {"allele": "B", "pair_id": "p3", "geometry_status": "complete", "exposed_ca_rmsd_A_median": 9.0, "exposed_ca_rmsd_A_iqr": 1.0},
            {"allele": "A", "pair_id": "missing", "geometry_status": "missing", "exposed_ca_rmsd_A_median": "", "exposed_ca_rmsd_A_iqr": ""},
        ]
        rank_pair_rows(rows)
        self.assertEqual(next(row for row in rows if row["pair_id"] == "p1")["within_allele_rank"], 1)
        self.assertEqual(next(row for row in rows if row["pair_id"] == "p2")["within_allele_rank"], 2)
        self.assertEqual(next(row for row in rows if row["pair_id"] == "p3")["within_allele_rank"], 1)
        self.assertEqual(next(row for row in rows if row["pair_id"] == "missing")["within_allele_rank"], "")

    def test_cross_allele_consensus_excludes_incomplete_pairs(self):
        rows = []
        specifications = {
            "stable": [0.04, 0.05, 0.06, 0.07],
            "mixed": [0.01, 0.01, 0.01, 0.60],
            "partial": [0.01, 0.01, 0.01],
        }
        for pair, percentiles in specifications.items():
            for index, percentile in enumerate(percentiles):
                rows.append({
                    "allele": f"A{index}", "ebv_candidate_id": pair, "self_candidate_id": "S",
                    "geometry_status": "complete", "within_allele_percentile": percentile,
                    "exposed_ca_rmsd_A_median": 1.0,
                })
        result = summarize_cross_allele(rows)
        stable = next(row for row in result if row["ebv_candidate_id"] == "stable")
        mixed = next(row for row in result if row["ebv_candidate_id"] == "mixed")
        partial = next(row for row in result if row["ebv_candidate_id"] == "partial")
        self.assertEqual(stable["cross_allele_consensus_rank"], 1)
        self.assertEqual(mixed["cross_allele_consensus_rank"], 2)
        self.assertEqual(partial["cross_allele_consensus_rank"], "")

    def test_recovery_is_not_evaluable_when_one_seed_is_incomplete(self):
        rows = [
            {"seed": 104729, "formal_seed_evaluable": False, "seed_recovery_criterion_pass": ""},
            {"seed": 104759, "formal_seed_evaluable": True, "seed_recovery_criterion_pass": True},
        ]
        self.assertEqual(classify_recovery(rows), "not_evaluable_incomplete_calibration")
        complete = copy.deepcopy(rows)
        complete[0].update(formal_seed_evaluable=True, seed_recovery_criterion_pass=True)
        self.assertEqual(classify_recovery(complete), "recovered")
        complete[1]["seed_recovery_criterion_pass"] = False
        self.assertEqual(classify_recovery(complete), "failed_calibration")


if __name__ == "__main__":
    unittest.main()
