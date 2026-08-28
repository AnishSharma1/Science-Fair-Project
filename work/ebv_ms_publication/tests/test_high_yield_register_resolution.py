import csv
import json
from pathlib import Path
import tempfile
import unittest


from high_yield_register_resolution import (
    build_experimental_peptide_panel,
    build_register_resolution_gate,
    enumerate_windows,
    evaluate_target_windows,
    prioritize_register_confirmation,
    summarize_target_windows,
)
from build_high_yield_register_resolution import build_package


def panel_row(pair_id, blosum, physicochemical, identity, *, role="n3"):
    return {
        "pair_id": pair_id,
        "row_role": role,
        "tcr_facing_blosum62_similarity": blosum,
        "tcr_face_physicochemical_mismatch": physicochemical,
        "tcr_facing_sequence_identity": identity,
        "local_surface_score": 0.5,
        "exposure_weighted_backbone_rmsd_A_q75": 1.0,
        "full_core_ca_rmsd_A_q75": 1.0,
        "anchor_ca_rmsd_A_q75": 1.0,
        "register_robust": False,
    }


class WindowEnumerationTests(unittest.TestCase):
    def test_enumerates_every_fully_contained_nine_residue_window_with_one_based_starts(self):
        windows = enumerate_windows("ABCDEFGHIJK", width=9)
        self.assertEqual(
            windows,
            [
                {"start_1_based": 1, "core": "ABCDEFGHI"},
                {"start_1_based": 2, "core": "BCDEFGHIJ"},
                {"start_1_based": 3, "core": "CDEFGHIJK"},
            ],
        )

    def test_short_sequence_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "shorter than window width"):
            enumerate_windows("ABCDEFGH", width=9)


class SensitivityTests(unittest.TestCase):
    def test_recomputes_sequence_features_and_ranks_each_window_against_frozen_n3_rows(self):
        target = {
            "target_id": "T1",
            "allele": "HLA-DRB1*15:01",
            "pair_id": "TARGET",
            "ebv_sequence": "AAACDEFGHIK",
            "self_sequence": "AAACDEFGHIK",
            "ebv_core_p1_p9": "ACDEFGHIK",
            "self_core_p1_p9": "ACDEFGHIK",
            "ebv_declared_core_start_1_based": 3,
            "self_declared_core_start_1_based": 3,
        }
        panel = [
            panel_row("TARGET", 0.0, 1.0, 0.0, role="target"),
            panel_row("N3_A", 0.8, 0.2, 0.8),
            panel_row("N3_B", 0.7, 0.3, 0.6),
            panel_row("N3_C", 0.6, 0.4, 0.4),
        ]
        rows = evaluate_target_windows(target, panel, local_shift=1)
        declared = next(row for row in rows if row["is_declared_window_pair"])
        shifted = next(
            row
            for row in rows
            if row["ebv_window_start_1_based"] == 1
            and row["self_window_start_1_based"] == 3
        )
        self.assertEqual(declared["panel_primary_rank"], 1)
        self.assertTrue(declared["capture_at_3"])
        self.assertTrue(declared["is_local_shift_window_pair"])
        self.assertEqual(declared["tcr_facing_sequence_identity"], 1.0)
        self.assertTrue(declared["structure_abstained_for_alternate_register"])
        self.assertFalse(shifted["is_local_shift_window_pair"])
        self.assertGreater(shifted["panel_primary_rank"], 1)

    def test_panel_must_retain_one_target_and_at_least_one_n3_row(self):
        target = {
            "target_id": "T1",
            "allele": "HLA-DRB1*15:01",
            "pair_id": "TARGET",
            "ebv_sequence": "ABCDEFGHI",
            "self_sequence": "ABCDEFGHI",
            "ebv_core_p1_p9": "ABCDEFGHI",
            "self_core_p1_p9": "ABCDEFGHI",
            "ebv_declared_core_start_1_based": 1,
            "self_declared_core_start_1_based": 1,
        }
        with self.assertRaisesRegex(ValueError, "one target and at least one N3"):
            evaluate_target_windows(target, [panel_row("TARGET", 0, 1, 0, role="target")])


class SummaryAndGateTests(unittest.TestCase):
    def test_summary_distinguishes_all_window_local_and_declared_only_support(self):
        base = {
            "target_id": "T1",
            "is_declared_window_pair": False,
            "is_local_shift_window_pair": False,
            "capture_at_3": False,
            "panel_primary_rank": 8,
            "ebv_window_core": "AAAAAAAAA",
            "self_window_core": "CCCCCCCCC",
        }
        rows = [
            {**base, "is_declared_window_pair": True, "is_local_shift_window_pair": True, "capture_at_3": True, "panel_primary_rank": 1},
            {**base, "is_local_shift_window_pair": True, "capture_at_3": True, "panel_primary_rank": 2},
            {**base, "capture_at_3": False, "panel_primary_rank": 8},
        ]
        summary = summarize_target_windows(rows)
        self.assertEqual(summary["register_resolution_status"], "local_shift_robust_only")
        self.assertEqual(summary["worst_local_shift_rank"], 2)
        self.assertEqual(summary["worst_all_window_rank"], 8)

        rows[1]["capture_at_3"] = False
        rows[1]["panel_primary_rank"] = 5
        summary = summarize_target_windows(rows)
        self.assertEqual(summary["register_resolution_status"], "declared_window_only")

    def test_gate_never_unlocks_discovery_or_freezes_weights(self):
        gate = build_register_resolution_gate([
            {"target_id": "T1", "register_resolution_status": "all_window_robust"},
            {"target_id": "T2", "register_resolution_status": "declared_window_only"},
        ])
        self.assertEqual(gate["status"], "declared_register_dependent")
        self.assertFalse(gate["discovery_unlock_allowed"])
        self.assertFalse(gate["weights_frozen"])
        self.assertFalse(gate["specificity_claim_allowed"])

    def test_experimental_priority_favors_more_local_captures_then_declared_rank(self):
        rows = [
            {"target_id": "B", "local_shift_capture_fraction": 0.2, "declared_window_rank": 1, "worst_local_shift_rank": 10},
            {"target_id": "A", "local_shift_capture_fraction": 0.4, "declared_window_rank": 2, "worst_local_shift_rank": 20},
            {"target_id": "C", "local_shift_capture_fraction": 0.4, "declared_window_rank": 1, "worst_local_shift_rank": 22},
        ]
        prioritized = prioritize_register_confirmation(rows)
        self.assertEqual([row["target_id"] for row in prioritized], ["C", "A", "B"])
        self.assertEqual([row["experimental_priority_rank"] for row in prioritized], [1, 2, 3])


class PeptidePanelTests(unittest.TestCase):
    def test_panel_includes_parent_peptides_and_declared_plus_local_shift_cores(self):
        target = {
            "target_id": "T1",
            "allele": "HLA-DRB1*15:01",
            "ebv_sequence": "ABCDEFGHIJKLMNO",
            "self_sequence": "PQRSTUVWXYZABCD",
            "ebv_declared_core_start_1_based": 4,
            "self_declared_core_start_1_based": 4,
        }
        rows = build_experimental_peptide_panel([target], local_shift=1)
        self.assertEqual(sum(row["peptide_role"] == "parent_assay_peptide" for row in rows), 2)
        self.assertEqual(sum(row["peptide_role"] == "register_discrimination_core" for row in rows), 6)
        self.assertEqual(len({row["panel_peptide_id"] for row in rows}), len(rows))
        self.assertTrue(all(row["proposed_not_ordered"] for row in rows))


class PackageIntegrationTests(unittest.TestCase):
    def test_real_package_has_eight_targets_and_deterministic_checksums(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            first_manifest = build_package(output_dir=Path(first))
            second_manifest = build_package(output_dir=Path(second))
            self.assertEqual(first_manifest["target_count"], 8)
            self.assertEqual(first_manifest["n3_pairs_per_target"], 25)
            self.assertEqual(first_manifest["file_checksums"], second_manifest["file_checksums"])
            with (Path(first) / "register_resolution_gate.json").open(encoding="utf-8") as handle:
                gate = json.load(handle)
            self.assertFalse(gate["discovery_unlock_allowed"])
            with (Path(first) / "target_register_summary.csv").open(newline="", encoding="utf-8") as handle:
                summaries = list(csv.DictReader(handle))
            self.assertEqual(len(summaries), 8)


if __name__ == "__main__":
    unittest.main()
