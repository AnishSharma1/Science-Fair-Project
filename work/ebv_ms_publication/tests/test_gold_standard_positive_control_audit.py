import copy
import unittest
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from gold_standard_positive_control_audit import build_gold_standard_audit  # noqa: E402


def valid_registry():
    return [
        {
            "biological_system_id": "SYS_BALF5_MBP_HY2E11",
            "evidence_tier": "E1_exact_pmhc_positive",
            "receptor_modality": "human_T_cell_shared_clone",
            "receptor_or_clone_id": "Hy.2E11",
            "assay_type": "same clone/TCR recognition plus pMHC crystal structures",
            "primary_source": "PMID:12244309",
            "doi": "10.1038/ni835",
            "viral_sequence": "TGGVYHFVKKHVHES",
            "viral_hla": "HLA-DRB5*01:01",
            "self_sequence": "ENPVVHFFKNIVTPR",
            "self_hla": "HLA-DRB1*15:01",
            "tcell_positive_denominator": "True",
        },
        {
            "biological_system_id": "SUPPORTIVE_ONLY",
            "evidence_tier": "E3_supportive_tcell",
            "receptor_modality": "human_T_cell_clones_pooled_targets",
            "tcell_positive_denominator": "False",
        },
    ]


def valid_metrics():
    return [
        {"metric": "EBV_structure", "value": "1H15"},
        {"metric": "MBP_structure", "value": "1BX2"},
        {"metric": "seven_position_core_CA_RMSD_A", "value": "0.838"},
    ]


def valid_seeds():
    return [
        {
            "seed": "104729", "available_primary_count": "17", "expected_primary_count": "26",
            "available_rank": "1", "positive_exposed_ca_rmsd_median_A": "0.490404",
            "available_equal_weight_control_median_A": "5.042155", "formal_seed_evaluable": "False",
            "seed_recovery_criterion_pass": "", "formal_seed_status": "not_evaluable_incomplete_seed",
        },
        {
            "seed": "104759", "available_primary_count": "26", "expected_primary_count": "26",
            "available_rank": "1", "positive_exposed_ca_rmsd_median_A": "0.495705",
            "available_equal_weight_control_median_A": "6.581933", "formal_seed_evaluable": "True",
            "seed_recovery_criterion_pass": "True", "formal_seed_status": "pass",
        },
    ]


def valid_recovery():
    return [{
        "biological_system_id": "SYS_BALF5_MBP_HY2E11",
        "recovery_status": "not_evaluable_incomplete_calibration",
    }]


class GoldStandardPositiveControlAuditTests(unittest.TestCase):
    def test_captures_locked_positive_without_calling_incomplete_seed_formal(self):
        system, seeds, summary = build_gold_standard_audit(
            valid_registry(), valid_metrics(), valid_seeds(), valid_recovery()
        )
        self.assertTrue(system["gold_standard_eligible"])
        self.assertFalse(system["used_for_score_tuning"])
        self.assertEqual(summary["gold_standard_independent_system_count"], 1)
        self.assertEqual(summary["capture_at_1_available_seed_fraction"], 1.0)
        self.assertEqual(summary["formal_evaluable_seed_count"], 1)
        self.assertEqual(summary["formal_seed_pass_count"], 1)
        self.assertEqual(seeds[0]["audit_status"], "available_set_capture_formal_incomplete")
        self.assertEqual(seeds[1]["audit_status"], "formal_pass")

    def test_rejects_wrong_hla_arm_even_when_pair_is_marked_positive(self):
        registry = valid_registry()
        registry[0]["viral_hla"] = "HLA-DRB1*15:01"
        with self.assertRaisesRegex(ValueError, "viral_hla"):
            build_gold_standard_audit(registry, valid_metrics(), valid_seeds(), valid_recovery())

    def test_supportive_records_cannot_inflate_gold_standard_denominator(self):
        registry = valid_registry()
        registry[1]["tcell_positive_denominator"] = "True"
        with self.assertRaisesRegex(ValueError, "exactly one independent system"):
            build_gold_standard_audit(registry, valid_metrics(), valid_seeds(), valid_recovery())

    def test_requires_both_fixed_seeds(self):
        with self.assertRaisesRegex(ValueError, "expected calibration seeds"):
            build_gold_standard_audit(
                valid_registry(), valid_metrics(), valid_seeds()[1:], valid_recovery()
            )


if __name__ == "__main__":
    unittest.main()

