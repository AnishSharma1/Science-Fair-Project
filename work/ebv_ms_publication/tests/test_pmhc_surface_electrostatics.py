import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pmhc_surface_electrostatics import (
    APBSParameters,
    GridSpec,
    assign_electrostatic_context,
    build_apbs_input,
    build_shared_grid,
    build_common_accessible_field_patch,
    candidate_exposed_histidines,
    carbo_similarity,
    dielectric_robustness,
    hodgkin_similarity,
    parse_open_dx,
    potential_rmse,
    rank_panel,
    select_supported_comparators,
    sign_agreement_fraction,
    summarize_electrostatic_ensemble,
    align_model_to_reference,
    surface_patch_points,
    write_model_pdb,
    trilinear_sample,
)


class MetricTests(unittest.TestCase):
    def test_histidine_positions_use_one_based_register_coordinates(self):
        self.assertEqual(candidate_exposed_histidines("YHFVKKHVH"), (2, 7))
        self.assertEqual(candidate_exposed_histidines("VYHFVKKHV"), (3, 8))

    def test_hodgkin_identity_inversion_and_zero_handling(self):
        values = np.asarray([1.0, -2.0, 3.0])
        self.assertAlmostEqual(hodgkin_similarity(values, values), 1.0)
        self.assertAlmostEqual(hodgkin_similarity(values, -values), -1.0)
        self.assertAlmostEqual(hodgkin_similarity(np.zeros(3), np.zeros(3)), 1.0)
        self.assertAlmostEqual(hodgkin_similarity(values, np.zeros(3)), 0.0)

    def test_secondary_metrics_are_deterministic(self):
        left = np.asarray([1.0, -2.0, 3.0, 0.0])
        right = np.asarray([1.0, -1.0, -3.0, 0.0])
        expected_carbo = float(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)))
        self.assertAlmostEqual(carbo_similarity(left, right), expected_carbo)
        self.assertAlmostEqual(sign_agreement_fraction(left, right), 0.75)
        self.assertAlmostEqual(potential_rmse(left, right), math.sqrt(37.0 / 4.0))

    def test_ensemble_uses_conservative_similarity_and_error_quantiles(self):
        rows = []
        for index in range(25):
            rows.append(
                {
                    "hodgkin_similarity": index / 24.0,
                    "carbo_similarity": 1.0 - index / 48.0,
                    "sign_agreement_fraction": 0.5 + index / 100.0,
                    "potential_rmse": float(index),
                }
            )
        summary = summarize_electrostatic_ensemble(rows)
        self.assertEqual(summary["model_combination_count"], 25)
        self.assertAlmostEqual(summary["hodgkin_similarity_q25"], 0.25)
        self.assertAlmostEqual(summary["potential_rmse_q75"], 18.0)

    def test_ensemble_rejects_incomplete_five_by_five_comparison(self):
        with self.assertRaisesRegex(ValueError, "exactly 25"):
            summarize_electrostatic_ensemble(
                [{"hodgkin_similarity": 1, "carbo_similarity": 1,
                  "sign_agreement_fraction": 1, "potential_rmse": 0}] * 24
            )


class GridAndDxTests(unittest.TestCase):
    @staticmethod
    def _model(offset=(0.0, 0.0, 0.0)):
        ox, oy, oz = offset
        def residue(aa, x):
            return {
                "aa": aa,
                "atoms": [
                    {"name": "N", "element": "N", "xyz": (x + ox, oy, oz)},
                    {"name": "CA", "element": "C", "xyz": (x + ox, 1.0 + oy, oz)},
                    {"name": "C", "element": "C", "xyz": (x + ox, 2.0 + oy, oz)},
                    {"name": "CB", "element": "C", "xyz": (x + ox, 1.0 + oy, 1.5 + oz)},
                ],
                "bfactors": [90.0] * 4,
            }
        return {
            "A": [residue("A", float(i)) for i in range(3)],
            "B": [residue("V", float(i + 4)) for i in range(3)],
            "C": [residue("H", float(i + 8)) for i in range(9)],
        }

    def test_alignment_is_rotation_translation_invariant(self):
        reference = self._model()
        moving = self._model(offset=(11.0, -4.0, 2.0))
        aligned, fit_rmsd = align_model_to_reference(moving, reference, groove_residue_count=3)
        self.assertLess(fit_rmsd, 1e-8)
        np.testing.assert_allclose(
            [atom["xyz"] for atom in aligned["C"][0]["atoms"]],
            [atom["xyz"] for atom in reference["C"][0]["atoms"]],
            atol=1e-8,
        )

    def test_surface_patch_uses_requested_core_positions_and_is_deterministic(self):
        model = self._model()
        first = surface_patch_points(model, core_start_1_based=1, samples_per_atom=24)
        second = surface_patch_points(model, core_start_1_based=1, samples_per_atom=24)
        self.assertGreater(len(first), 0)
        np.testing.assert_allclose(first, second)

    def test_common_field_patch_is_position_matched_and_solvent_accessible(self):
        core = np.column_stack((np.arange(9, dtype=float), np.zeros(9), np.ones(9) * 3.0))
        groove = np.column_stack((np.arange(18, dtype=float) / 2.0, np.zeros(18), np.zeros(18)))
        atoms = [
            (
                np.asarray([[4.0, 0.0, 3.0], [4.0, 1.5, 3.0]]),
                np.asarray([1.7, 1.7]),
            )
        ]
        patch, metadata = build_common_accessible_field_patch(
            core,
            groove,
            atoms,
            probe_radius_A=1.4,
            minimum_clearance_A=0.25,
            height_step_A=0.25,
            maximum_height_A=12.0,
        )
        self.assertEqual(patch.shape, (25, 3))
        self.assertEqual({row["core_position"] for row in metadata}, {2, 3, 5, 7, 8})
        self.assertTrue(all(row["minimum_clearance_A"] >= 0.25 - 1e-8 for row in metadata))
        self.assertTrue(all(row["height_A"] >= 2.0 for row in metadata))

    def test_pdb_export_preserves_chain_order_and_can_drop_peptide(self):
        with tempfile.TemporaryDirectory() as directory:
            full = Path(directory) / "full.pdb"
            hla = Path(directory) / "hla.pdb"
            write_model_pdb(self._model(), full, include_peptide=True)
            write_model_pdb(self._model(), hla, include_peptide=False)
            full_text = full.read_text()
            hla_text = hla.read_text()
        self.assertIn(" A   1", full_text)
        self.assertIn(" B   1", full_text)
        self.assertIn(" C   1", full_text)
        self.assertNotIn(" C   1", hla_text)

    def test_shared_grid_uses_union_bounds_and_multigrid_dimensions(self):
        bounds = [
            (np.asarray([-10.0, -5.0, 0.0]), np.asarray([10.0, 5.0, 20.0])),
            (np.asarray([-12.0, -4.0, -1.0]), np.asarray([9.0, 8.0, 18.0])),
        ]
        patch = np.asarray([[0.0, 0.0, 9.0], [2.0, 1.0, 11.0]])
        grid = build_shared_grid(bounds, patch, padding_A=12.0, maximum_spacing_A=0.5)
        self.assertTrue(all((dimension - 1) % 32 == 0 for dimension in grid.dime))
        self.assertTrue(all(spacing <= 0.5 for spacing in grid.fine_spacing_A))
        self.assertEqual(grid.coarse_min_A, (-24.0, -17.0, -13.0))
        self.assertEqual(grid.coarse_max_A, (22.0, 20.0, 32.0))

    def test_shared_grid_coarse_bounds_enclose_surface_fine_grid(self):
        bounds = [(np.asarray([0.0, 0.0, 0.0]), np.asarray([10.0, 10.0, 10.0]))]
        patch = np.asarray([[12.5, 5.0, 5.0], [13.0, 6.0, 5.0]])
        grid = build_shared_grid(bounds, patch, padding_A=2.0, maximum_spacing_A=0.5)
        fine_min = np.asarray(grid.fine_center_A) - np.asarray(grid.fine_length_A) / 2.0
        fine_max = np.asarray(grid.fine_center_A) + np.asarray(grid.fine_length_A) / 2.0
        self.assertTrue(np.all(np.asarray(grid.coarse_min_A) <= fine_min))
        self.assertTrue(np.all(np.asarray(grid.coarse_max_A) >= fine_max))

    def test_open_dx_parse_and_trilinear_sampling(self):
        dx = """object 1 class gridpositions counts 2 2 2
origin 0 0 0
delta 1 0 0
delta 0 1 0
delta 0 0 1
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
0 1 1 2 1 2 2 3
attribute \"dep\" string \"positions\"
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.dx"
            path.write_text(dx, encoding="ascii")
            grid = parse_open_dx(path)
        self.assertEqual(grid.values.shape, (2, 2, 2))
        sampled = trilinear_sample(grid, np.asarray([[0.5, 0.5, 0.5]]))
        self.assertAlmostEqual(float(sampled[0]), 1.5)

    def test_apbs_input_freezes_primary_conditions_and_grid(self):
        grid = GridSpec(
            dime=(65, 65, 65),
            coarse_length_A=(80.0, 80.0, 80.0),
            coarse_center_A=(0.0, 0.0, 0.0),
            fine_length_A=(30.0, 30.0, 30.0),
            fine_center_A=(1.0, 2.0, 3.0),
            coarse_min_A=(-40.0, -40.0, -40.0),
            coarse_max_A=(40.0, 40.0, 40.0),
            fine_spacing_A=(30.0 / 64, 30.0 / 64, 30.0 / 64),
        )
        text = build_apbs_input(Path("model.pqr"), Path("potential"), grid, APBSParameters())
        for expected in ("pdie 2", "sdie 78.5", "temp 298.15", "srad 1.4"):
            self.assertIn(expected, text)
        self.assertEqual(text.count("ion charge"), 2)
        self.assertIn("dime 65 65 65", text)


class SelectionAndGateTests(unittest.TestCase):
    @staticmethod
    def _candidate(identifier, protein, core, sequence, rank, *, allele="HLA-DRB1*15:01"):
        return {
            "candidate_id": identifier,
            "accession": protein,
            "protein": protein,
            "core": core,
            "predicted_core": core,
            "sequence": sequence,
            "allele": allele,
            "netmhciipan_el_percentile": rank,
            "mixmhc2pred_percentile": rank,
            "binding_consensus": True,
            "register_consensus": True,
            "declared_core_match": True,
            "model_count": 5,
            "surface_status": "complete",
        }

    def test_comparator_selection_is_exact_hla_supported_unique_and_score_blind(self):
        target = self._candidate("TARGET", "BALF5", "VYHFVKKHV", "TGGVYHFVKKHVHES", 2.0)
        candidates = []
        for index in range(1, 8):
            core = f"AAAAAAAA{chr(65 + index)}"
            candidates.append(
                self._candidate(f"C{index}", f"P{index}", core, f"GGG{core}GGG", float(index))
            )
        candidates.append(self._candidate("WRONG", "PX", "CCCCCCCCC", "M" * 15, 1.0, allele="HLA-DRB1*13:03"))
        candidates.append({**self._candidate("BADREG", "PY", "DDDDDDDDD", "M" * 15, 1.0), "register_consensus": False})
        selected, provenance = select_supported_comparators(
            target, candidates, excluded_candidate_ids={"C1"}, count=5, seed=271828
        )
        self.assertEqual(len(selected), 5)
        self.assertTrue(all(row["allele"] == target["allele"] for row in selected))
        self.assertNotIn("C1", {row["candidate_id"] for row in selected})
        reasons = {row["candidate_id"]: row["eligibility_reason"] for row in provenance}
        self.assertEqual(reasons["WRONG"], "wrong_hla")
        self.assertEqual(reasons["BADREG"], "register_not_supported")

    def test_panel_rank_uses_q25_hodgkin_then_pair_id(self):
        rows = [
            {"pair_id": "B", "hodgkin_similarity_q25": 0.8},
            {"pair_id": "A", "hodgkin_similarity_q25": 0.8},
            {"pair_id": "C", "hodgkin_similarity_q25": 0.2},
        ]
        ranked = rank_panel(rows)
        self.assertEqual([(row["pair_id"], row["electrostatic_rank"]) for row in ranked], [("A", 1), ("B", 2), ("C", 3)])

    def test_context_requires_top_three_and_complete_register_qc(self):
        self.assertEqual(assign_electrostatic_context(3, True, True), "electrostatic_context_supportive")
        self.assertEqual(assign_electrostatic_context(4, True, True), "electrostatic_context_not_supportive")
        self.assertEqual(assign_electrostatic_context(1, False, True), "not_evaluable")
        self.assertEqual(assign_electrostatic_context(1, True, False), "not_evaluable")

    def test_dielectric_robustness_requires_identical_classification(self):
        self.assertTrue(dielectric_robustness({2.0: "supportive", 4.0: "supportive", 8.0: "supportive"}))
        self.assertFalse(dielectric_robustness({2.0: "supportive", 4.0: "supportive", 8.0: "not_supportive"}))
        self.assertFalse(dielectric_robustness({2.0: "supportive", 4.0: "not_evaluable", 8.0: "supportive"}))


if __name__ == "__main__":
    unittest.main()
