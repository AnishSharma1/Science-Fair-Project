import unittest
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_control_calibrated_hla_rankings import (  # noqa: E402
    ALLELES,
    control_reference_metrics,
    derive_formal_control_reference,
    rank_within_hla,
)


def calibration_rows(seed=104759):
    rows = [{
        "seed": str(seed), "analysis_set": "primary_rank_of_26", "geometry_status": "complete",
        "pair_role": "E1_positive", "exposed_ca_rmsd_A_median": "0.5",
    }]
    rows.extend({
        "seed": str(seed), "analysis_set": "primary_rank_of_26", "geometry_status": "complete",
        "pair_role": "full_decoy", "exposed_ca_rmsd_A_median": str(value),
    } for value in range(1, 26))
    return rows


def seed_rows(formal=True):
    return [{
        "seed": "104759", "formal_seed_evaluable": str(formal),
    }]


def pair(allele, pair_id, value, iqr=0.2, status="complete"):
    return {
        "allele": allele,
        "geometry_status": status,
        "pair_id": pair_id,
        "ebv_candidate_id": f"E_{pair_id}",
        "ebv_protein": "EBV",
        "ebv_sequence": "ABCDEFGHIJKLMNO",
        "ebv_predicted_core": "BCDEFGHIJ",
        "ebv_binding_percentile_rank": "1.0",
        "ebv_source_certainty": "exact",
        "self_candidate_id": f"S_{pair_id}",
        "self_protein": "SELF",
        "self_sequence": "PQRSTUVWXYZABCD",
        "self_predicted_core": "QRSTUVWXY",
        "self_binding_percentile_rank": "1.0",
        "self_source_certainty": "exact",
        "exposed_ca_rmsd_A_median": str(value),
        "exposed_ca_rmsd_A_iqr": str(iqr),
        "exposed_ca_rmsd_A_q25": str(value - 0.1),
        "exposed_ca_rmsd_A_q75": str(value + 0.1),
        "model_combination_count": "25",
        "primary_endpoint": "exposed RMSD",
    }


class ControlCalibratedHlaRankingTests(unittest.TestCase):
    def test_formal_reference_requires_complete_positive_plus_25_decoys(self):
        reference = derive_formal_control_reference(calibration_rows(), seed_rows())
        self.assertEqual(reference["positive_available_rank"], 1)
        self.assertEqual(reference["available_primary_count"], 26)
        self.assertEqual(reference["available_decoy_count"], 25)
        self.assertEqual(reference["decoy_exposed_ca_rmsd_median_A"], 13)

    def test_incomplete_seed_cannot_set_formal_control_index(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            derive_formal_control_reference(calibration_rows(), seed_rows(formal=False))

    def test_control_index_has_locked_positive_and_decoy_landmarks(self):
        reference = derive_formal_control_reference(calibration_rows(), seed_rows())
        positive = control_reference_metrics(0.5, reference)
        decoy = control_reference_metrics(13, reference)
        self.assertEqual(positive["control_separation_index"], 1.0)
        self.assertEqual(decoy["control_separation_index"], 0.0)
        self.assertEqual(positive["control_geometry_band"], "at_or_below_gold_positive_median")
        self.assertEqual(positive["control_metric_interpretation"], "descriptive_method_reference_not_probability")

    def test_each_hla_is_ranked_independently_and_missing_pairs_are_retained(self):
        reference = derive_formal_control_reference(calibration_rows(), seed_rows())
        rows = [
            pair(ALLELES[0], "a_slow", 5.0),
            pair(ALLELES[0], "a_fast", 1.0),
            pair(ALLELES[1], "b_only", 20.0),
            pair(ALLELES[0], "a_missing", 0.0, status="missing"),
        ]
        ranked, missing = rank_within_hla(rows, reference)
        self.assertEqual([row["pair_id"] for row in ranked[ALLELES[0]]], ["a_fast", "a_slow"])
        self.assertEqual(ranked[ALLELES[1]][0]["hla_rank"], 1)
        self.assertEqual(ranked[ALLELES[1]][0]["hla_percentile"], 0.0)
        self.assertEqual(ranked[ALLELES[0]][0]["rank_scope"], "within_hla_only")
        self.assertNotIn("cross_allele_consensus_rank", ranked[ALLELES[0]][0])
        self.assertEqual(missing[0]["ranking_status"], "not_ranked_missing_geometry")


if __name__ == "__main__":
    unittest.main()

