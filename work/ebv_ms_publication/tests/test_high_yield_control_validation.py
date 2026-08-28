import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest


from high_yield_control_validation import (
    CLAIM_BOUNDARY,
    FROZEN_TARGETS,
    build_n3_panel,
    build_ranking_context_gate,
    build_specificity_gate,
    rank_panel_rows,
    select_comparator_arms,
    validate_frozen_targets,
)
from build_high_yield_control_validation import build_package


def arm(candidate_id, core, percentile, *, allele="HLA-DRB1*15:01", protein="P"):
    return {
        "candidate_id": candidate_id,
        "allele": allele,
        "sequence": f"AA{core}AAAA",
        "core": core,
        "binding_percentile": percentile,
        "protein": protein,
        "accession": candidate_id.split("_")[0],
        "surface_status": "complete",
        "model_count": 5,
    }


def pair(pair_id, blosum, physicochemical, identity, surface):
    return {
        "allele": "HLA-DRB1*15:01",
        "pair_id": pair_id,
        "tcr_facing_blosum62_similarity": blosum,
        "tcr_face_physicochemical_mismatch": physicochemical,
        "tcr_facing_sequence_identity": identity,
        "local_surface_score": surface,
        "exposure_weighted_backbone_rmsd_A_q75": surface + 0.1,
        "full_core_ca_rmsd_A_q75": surface + 0.2,
        "anchor_ca_rmsd_A_q75": surface + 0.3,
    }


class FrozenTargetTests(unittest.TestCase):
    def test_registry_has_three_unique_targets_per_hla(self):
        summary = validate_frozen_targets(FROZEN_TARGETS)
        self.assertEqual(summary["target_count"], 12)
        self.assertEqual(summary["targets_per_hla"], {
            "HLA-DRB1*03:01": 3,
            "HLA-DRB1*08:01": 3,
            "HLA-DRB1*13:03": 3,
            "HLA-DRB1*15:01": 3,
        })

    def test_duplicate_target_is_rejected(self):
        duplicate = list(FROZEN_TARGETS) + [copy.deepcopy(FROZEN_TARGETS[0])]
        with self.assertRaisesRegex(ValueError, "duplicate frozen target"):
            validate_frozen_targets(duplicate)


class ComparatorSelectionTests(unittest.TestCase):
    def test_selection_is_exact_hla_eligible_complete_and_score_blind(self):
        target = arm("TARGET", "ABCDEFGHI", 10.0)
        eligible = [
            arm(f"A{i}_ID", f"ABCDE{i:04d}"[-9:], 10.0 + i / 10, protein=f"P{i}")
            for i in range(8)
        ]
        candidates = eligible + [
            arm("WRONG_HLA", "LMNOPQRST", 10.0, allele="HLA-DRB1*03:01"),
            arm("WEAK_BIND", "QRSTUVWXZ", 20.1),
            {**arm("MISSING_MODELS", "STUVWXYZA", 10.0), "model_count": 4},
            {**arm("INCOMPLETE", "TUVWXYZAB", 10.0), "surface_status": "missing"},
            arm("BANNED", "UVWXYZABC", 10.0),
        ]
        selected, provenance = select_comparator_arms(
            target,
            candidates,
            allele="HLA-DRB1*15:01",
            arm_class="ebv",
            excluded_candidate_ids={"BANNED"},
            count=5,
            seed=271828,
        )
        self.assertEqual(len(selected), 5)
        self.assertTrue(all(row["allele"] == target["allele"] for row in selected))
        self.assertTrue(all(float(row["binding_percentile"]) <= 20 for row in selected))
        self.assertTrue(all(row["model_count"] == 5 for row in selected))
        self.assertNotIn("BANNED", {row["candidate_id"] for row in selected})
        reasons = {row["candidate_id"]: row["eligibility_reason"] for row in provenance}
        self.assertEqual(reasons["WRONG_HLA"], "wrong_hla")
        self.assertEqual(reasons["WEAK_BIND"], "binding_percentile_above_20")
        self.assertEqual(reasons["MISSING_MODELS"], "incomplete_model_ensemble")
        self.assertEqual(reasons["INCOMPLETE"], "incomplete_surface_features")
        self.assertEqual(reasons["BANNED"], "excluded_frozen_or_control_arm")
        selected_again, _ = select_comparator_arms(
            target,
            list(reversed(candidates)),
            allele="HLA-DRB1*15:01",
            arm_class="ebv",
            excluded_candidate_ids={"BANNED"},
            count=5,
            seed=271828,
        )
        self.assertEqual(
            [row["candidate_id"] for row in selected],
            [row["candidate_id"] for row in selected_again],
        )

    def test_insufficient_comparators_do_not_relax_rules(self):
        target = arm("TARGET", "ABCDEFGHI", 10.0)
        with self.assertRaisesRegex(ValueError, "insufficient eligible ebv comparators"):
            select_comparator_arms(
                target,
                [arm("A", "BCDEFGHIJ", 10.0)],
                allele=target["allele"],
                arm_class="ebv",
                excluded_candidate_ids=set(),
                count=5,
                seed=271828,
            )


class PanelTests(unittest.TestCase):
    def test_panel_cross_product_has_target_plus_25_unique_n3_pairs(self):
        target = pair("TARGET", 9, 0.1, 0.8, 0.2)
        left = [arm(f"E{i}", f"AAAAAAA{i}A", 10.0) for i in range(5)]
        right = [arm(f"S{i}", f"CCCCCCC{i}C", 10.0) for i in range(5)]
        lookup = {}
        for left_arm in left:
            for right_arm in right:
                pid = f"{left_arm['candidate_id']}|{right_arm['candidate_id']}"
                lookup[(left_arm["candidate_id"], right_arm["candidate_id"])] = pair(
                    pid, 1, 1, 0, 1
                )
        panel = build_n3_panel(target, left, right, lookup)
        self.assertEqual(len(panel), 26)
        self.assertEqual(sum(row["row_role"] == "target" for row in panel), 1)
        self.assertEqual(sum(row["row_role"] == "n3" for row in panel), 25)
        self.assertEqual(len({row["pair_id"] for row in panel}), 26)

    def test_rank_uses_frozen_v3_order_and_reports_diagnostic_ranks(self):
        rows = [
            {**pair("target", 10, 0.5, 0.2, 0.9), "register_robust": True},
            {**pair("better_structure_lower_blosum", 9, 0.0, 1.0, 0.0), "register_robust": True},
            {**pair("same_blosum_better_physchem", 10, 0.4, 0.0, 1.0), "register_robust": True},
            {**pair("same_upstream_better_surface", 10, 0.5, 0.2, 0.1), "register_robust": True},
        ]
        ranked = rank_panel_rows(rows, seed=271828)
        ordered = sorted(ranked, key=lambda row: row["panel_primary_rank"])
        self.assertEqual(
            [row["pair_id"] for row in ordered],
            ["same_blosum_better_physchem", "same_upstream_better_surface", "target", "better_structure_lower_blosum"],
        )
        target = next(row for row in ranked if row["pair_id"] == "target")
        self.assertEqual(target["panel_primary_rank"], 3)
        self.assertEqual(target["panel_local_surface_rank"], 3)
        self.assertIn("panel_random_rank", target)

    def test_register_uncertain_rows_abstain_from_structural_primary_tie_break(self):
        rows = [
            {**pair("z_uncertain_better_surface", 10, 0.5, 0.2, 0.0), "register_robust": False},
            {**pair("a_uncertain_worse_surface", 10, 0.5, 0.2, 1.0), "register_robust": False},
        ]
        ranked = rank_panel_rows(rows, seed=271828)
        ordered = sorted(ranked, key=lambda row: row["panel_primary_rank"])
        self.assertEqual(
            [row["pair_id"] for row in ordered],
            ["a_uncertain_worse_surface", "z_uncertain_better_surface"],
        )
        self.assertTrue(all(row["panel_primary_structure_abstained"] for row in ranked))


class GateTests(unittest.TestCase):
    def test_top_three_is_supportive_rank_context_not_validation(self):
        summaries = [
            {"target_id": "A", "panel_status": "complete", "target_primary_rank": 3},
            {"target_id": "B", "panel_status": "complete", "target_primary_rank": 4},
            {"target_id": "C", "panel_status": "not_evaluable", "target_primary_rank": ""},
        ]
        gate = build_ranking_context_gate(summaries)
        self.assertEqual(gate["targets"][0]["ranking_context_status"], "rank_context_supportive")
        self.assertEqual(gate["targets"][1]["ranking_context_status"], "rank_context_not_supportive")
        self.assertEqual(gate["targets"][2]["ranking_context_status"], "not_evaluable")
        self.assertFalse(gate["discovery_unlock_allowed"])
        self.assertFalse(gate["weights_frozen"])
        self.assertEqual(gate["claim_boundary"], CLAIM_BOUNDARY)

    def test_n3_never_counts_as_specificity(self):
        gate = build_specificity_gate([])
        self.assertEqual(gate["status"], "not_evaluable_no_explicit_n1_n2")
        self.assertTrue(gate["n3_excluded_from_specificity"])
        self.assertFalse(gate["specificity_claim_allowed"])


class RealPackageTests(unittest.TestCase):
    def test_real_package_has_frozen_counts_and_locked_claim_boundaries(self):
        root = Path(__file__).resolve().parents[1]
        v3 = root / "processed/literature_grounded_hla2_rankings_v3_2026-08-27"
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "package"
            manifest = build_package(v3_dir=v3, out=out)
            self.assertEqual(manifest["target_count"], 12)
            self.assertEqual(manifest["panel_row_count"], 312)
            with (out / "panel_feature_matrix.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 312)
            self.assertEqual(sum(row["row_role"] == "target" for row in rows), 12)
            self.assertEqual(sum(row["row_role"] == "n3" for row in rows), 300)
            with (out / "panel_rank_summary.csv").open(newline="") as handle:
                summaries = list(csv.DictReader(handle))
            self.assertEqual(len(summaries), 12)
            gate = json.loads((out / "ranking_context_gate.json").read_text())
            self.assertFalse(gate["discovery_unlock_allowed"])
            self.assertFalse(gate["specificity_claim_allowed"])
            specificity = json.loads((out / "specificity_gate.json").read_text())
            self.assertTrue(specificity["n3_excluded_from_specificity"])
            protocol = json.loads((out / "protocol_lock.json").read_text())
            self.assertFalse(protocol["global_frozen_target_arm_exclusion_feasible"])
            self.assertEqual(
                protocol["comparator_target_exclusion_scope"],
                "current_panel_target_arms_plus_confirmed_control_ligands",
            )
            with (out / "global_exclusion_feasibility.csv").open(newline="") as handle:
                feasibility = list(csv.DictReader(handle))
            self.assertEqual(len(feasibility), 24)
            self.assertTrue(any(int(row["eligible_unique_arm_count"]) < 5 for row in feasibility))
            with (out / "control_system_census.csv").open(newline="") as handle:
                census = list(csv.DictReader(handle))
            self.assertGreaterEqual(len(census), 4)
            self.assertEqual(
                sum(row["independent_v3_validation_vote"] == "true" for row in census),
                0,
            )
            with (out / "SHA256SUMS.csv").open(newline="") as handle:
                checksums = list(csv.DictReader(handle))
            self.assertGreaterEqual(len(checksums), 10)
            for path in out.glob("*.csv"):
                self.assertNotIn(b"\r\n", path.read_bytes(), path.name)


if __name__ == "__main__":
    unittest.main()
