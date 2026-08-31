import tempfile
import unittest
from pathlib import Path

import numpy as np

from pmhc_surface_electrostatics import OpenDXGrid
from pmhc_surface_electrostatics_v2 import (
    align_local_groove,
    build_apbs_surface_input,
    build_control_gate,
    canonical_groove_frame,
    dense_lateral_grid,
    hierarchical_pair_summary,
    label_surface_regions,
    sample_outer_surface,
    standardize_pmhc_chains,
    validate_control_only_paths,
)
from run_pmhc_surface_electrostatics_v2_controls import (
    _apbs_command,
    _finder_duplicate_paths,
    _gate_requirements,
    _gate_invariants,
    _grid_definition,
    _normalize_lateral_arrays,
    _readme_rank_table,
    _write_figures,
    _write_svg_rank_figure,
    _write_matrix_archive,
)


class SurfaceGeometryTests(unittest.TestCase):
    @staticmethod
    def _model(offset=(0.0, 0.0, 0.0), extra_chain=False):
        shift = np.asarray(offset, dtype=float)

        def residue(aa, xyz):
            base = np.asarray(xyz, dtype=float) + shift
            return {
                "aa": aa,
                "atoms": [
                    {"name": "N", "element": "N", "xyz": tuple(base + [-0.2, 0.0, 0.0])},
                    {"name": "CA", "element": "C", "xyz": tuple(base)},
                    {"name": "C", "element": "C", "xyz": tuple(base + [0.2, 0.0, 0.0])},
                    {"name": "O", "element": "O", "xyz": tuple(base + [0.3, 0.1, 0.0])},
                ],
                "bfactors": [90.0] * 4,
            }

        model = {
            "X": [residue("A", (i, -3.0, 0.0)) for i in range(12)],
            "Y": [residue("V", (i, 3.0, 0.0)) for i in range(12)],
            "P": [residue("G", (i + 1.0, 0.0, 3.0)) for i in range(9)],
        }
        if extra_chain:
            model["T"] = [residue("L", (i, 0.0, 10.0)) for i in range(4)]
        return model

    def test_chain_standardization_removes_tcr_and_records_it(self):
        standardized, provenance = standardize_pmhc_chains(
            self._model(extra_chain=True), alpha_chain="X", beta_chain="Y", peptide_chain="P"
        )
        self.assertEqual(set(standardized), {"A", "B", "C"})
        self.assertEqual(provenance["excluded_chain_ids"], ["T"])
        self.assertTrue(provenance["tcr_or_other_protein_chains_removed"])

    def test_local_groove_alignment_uses_backbone_atoms_near_reference_core(self):
        reference, _ = standardize_pmhc_chains(
            self._model(), alpha_chain="X", beta_chain="Y", peptide_chain="P"
        )
        moving, _ = standardize_pmhc_chains(
            self._model(offset=(7.0, -2.0, 4.0)), alpha_chain="X", beta_chain="Y", peptide_chain="P"
        )
        aligned, qc = align_local_groove(moving, reference, reference_core_start_1_based=1)
        self.assertGreater(qc["matched_backbone_atom_count"], 20)
        self.assertLess(qc["fit_rmsd_A"], 1e-8)
        np.testing.assert_allclose(
            aligned["C"][0]["atoms"][1]["xyz"],
            reference["C"][0]["atoms"][1]["xyz"],
            atol=1e-8,
        )

    def test_canonical_frame_has_locked_orientation(self):
        core = np.column_stack((np.arange(9, dtype=float), np.zeros(9), np.ones(9) * 3.0))
        alpha = np.asarray([[0.0, -3.0, 0.0], [8.0, -3.0, 0.0]])
        beta = np.asarray([[0.0, 3.0, 0.0], [8.0, 3.0, 0.0]])
        frame = canonical_groove_frame(core, alpha, beta)
        np.testing.assert_allclose(frame.longitudinal, [1.0, 0.0, 0.0], atol=1e-8)
        np.testing.assert_allclose(frame.outward, [0.0, 0.0, 1.0], atol=1e-8)
        self.assertGreater(float(np.dot(frame.transverse, beta.mean(0) - alpha.mean(0))), 0.0)

    def test_dense_grid_meets_minimum_composite_density(self):
        core = np.column_stack((np.arange(9, dtype=float) * 3.5, np.zeros(9), np.ones(9) * 3.0))
        alpha = np.asarray([[0.0, -3.0, 0.0], [28.0, -3.0, 0.0]])
        beta = np.asarray([[0.0, 3.0, 0.0], [28.0, 3.0, 0.0]])
        frame = canonical_groove_frame(core, alpha, beta)
        points, metadata = dense_lateral_grid(core, frame, spacing_A=1.0)
        self.assertGreaterEqual(len(points), 500)
        self.assertEqual(len(points), len(metadata))
        self.assertLessEqual(min(row["longitudinal_A"] for row in metadata), -4.0)
        self.assertGreaterEqual(max(row["longitudinal_A"] for row in metadata), 32.0)
        self.assertEqual(min(row["transverse_A"] for row in metadata), -14.0)
        self.assertEqual(max(row["transverse_A"] for row in metadata), 14.0)

    def test_surface_crossing_and_offset_follow_accessibility_normal(self):
        axis = np.arange(7, dtype=float)
        values = np.zeros((7, 3, 3), dtype=float)
        values[3:, :, :] = 1.0
        accessibility = OpenDXGrid(
            origin=np.asarray([0.0, -1.0, -1.0]),
            deltas=np.ones(3),
            values=values,
        )
        potential = OpenDXGrid(
            origin=np.asarray([0.0, -1.0, -1.0]),
            deltas=np.ones(3),
            values=np.broadcast_to(axis[:, None, None], (7, 3, 3)).copy(),
        )
        result = sample_outer_surface(
            accessibility,
            potential,
            np.asarray([[0.0, 0.0, 0.0]]),
            np.asarray([1.0, 0.0, 0.0]),
            search_min_A=0.0,
            search_max_A=6.0,
            search_step_A=0.25,
            offset_A=0.5,
        )
        self.assertTrue(result[0]["covered"])
        self.assertAlmostEqual(result[0]["surface_height_A"], 2.5, places=6)
        np.testing.assert_allclose(result[0]["normal"], [1.0, 0.0, 0.0], atol=1e-6)
        self.assertAlmostEqual(result[0]["potential_kT_per_e"], 3.0, places=6)

    def test_apbs_input_requests_nonlinear_potential_and_accessibility(self):
        text = build_apbs_surface_input(
            Path("model.pqr"),
            Path("potential"),
            Path("accessibility"),
            dime=(65, 65, 65),
            lengths_A=(32.0, 32.0, 32.0),
            center_A=(0.0, 0.0, 0.0),
            solute_dielectric=4.0,
            linear=False,
            write_accessibility=True,
        )
        self.assertIn("npbe", text)
        self.assertIn("pdie 4", text)
        self.assertIn("write pot dx potential", text)
        self.assertIn("write smol dx accessibility", text)

    def test_apbs_coarse_box_encloses_actual_fine_mesh(self):
        atoms = np.asarray([
            [-18.0, -12.0, -10.0],
            [18.0, 12.0, 8.0],
        ])
        base_points = np.asarray([
            [-15.0, -14.0, 3.0],
            [15.0, 14.0, 3.0],
        ])
        grid = _grid_definition(
            atoms,
            base_points,
            np.asarray([0.0, 0.0, 1.0]),
            maximum_spacing_A=0.5,
        )
        fine_center = np.asarray(grid["fine_center_A"])
        fine_half = np.asarray(grid["fine_lengths_A"]) / 2.0
        coarse_center = np.asarray(grid["coarse_center_A"])
        coarse_half = np.asarray(grid["coarse_lengths_A"]) / 2.0
        self.assertTrue(np.all(coarse_center - coarse_half <= fine_center - fine_half))
        self.assertTrue(np.all(coarse_center + coarse_half >= fine_center + fine_half))

    def test_native_apbs_command_uses_host_input_path(self):
        output = Path("/tmp/control_package")
        input_path = output / "raw_calculations/apbs/model/apbs.in"
        command = _apbs_command(output, input_path, None)
        self.assertEqual(Path(command[0]).name, "apbs")
        self.assertEqual(command[1], str(input_path))

    def test_lateral_normalization_maps_variable_p1_p9_spans_to_one_shape(self):
        def arrays(longitudinal_count):
            count = longitudinal_count * 57
            return {
                "coverage": np.ones(count, dtype=bool),
                "height": np.arange(count, dtype=float),
                "normal": np.tile(np.asarray([0.0, 0.0, 1.0]), (count, 1)),
                "potential": np.linspace(-1.0, 1.0, count),
                "label": np.full(count, "peptide"),
            }

        short = _normalize_lateral_arrays(arrays(55), spacing_A=0.5)
        long = _normalize_lateral_arrays(arrays(67), spacing_A=0.5)
        self.assertEqual(len(short["coverage"]), 73 * 57)
        self.assertEqual(len(short["coverage"]), len(long["coverage"]))
        self.assertTrue(np.allclose(np.linalg.norm(short["normal"], axis=1), 1.0))

    def test_region_labels_use_nearest_contributing_atom(self):
        points = np.asarray([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        atoms = {
            "peptide": np.asarray([[0.1, 0.0, 0.0]]),
            "hla_alpha": np.asarray([[5.1, 0.0, 0.0]]),
            "hla_beta": np.asarray([[10.1, 0.0, 0.0]]),
        }
        self.assertEqual(label_surface_regions(points, atoms), ["peptide", "hla_alpha", "hla_beta"])


class EnsembleAndGateTests(unittest.TestCase):
    def test_gate_requirement_abstains_when_any_panel_decoy_is_incomplete(self):
        common = {
            "layer": "alphafold",
            "positive_pair_id": "PAIR",
            "panel_seed": "271828",
            "variant": "npbe_eps4_grid0.50",
            "density": "fine",
            "offset_A": 0.5,
        }
        rows = [
            {
                **common,
                "pair_id": "positive",
                "comparison_role": "positive",
                "status": "complete",
                "peptide_rank": 1,
                "composite_rank": 1,
                "shape_rank": 1,
            },
            {
                **common,
                "pair_id": "decoy",
                "comparison_role": "N3_pair_decoy",
                "status": "not_evaluable_surface_qc",
                "peptide_rank": "",
                "composite_rank": "",
                "shape_rank": "",
            },
        ]
        requirements = _gate_requirements(rows, [])
        self.assertEqual({row["status"] for row in requirements}, {"missing"})

    def test_matrix_archive_is_byte_deterministic(self):
        matrices = {
            ("alphafold", "PAIR", "271828", "DECOY", "primary", "fine", 0.5): {
                "peptide_hodgkin": np.arange(25, dtype=float).reshape(5, 5),
                "composite_hodgkin": np.ones((5, 5), dtype=float),
                "helix_hodgkin": np.zeros((5, 5), dtype=float),
                "surface_height_rmse_A": np.full((5, 5), 2.0),
                "surface_normal_mean_dot": np.full((5, 5), 0.75),
            }
        }
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            _write_matrix_archive(Path(first), matrices)
            _write_matrix_archive(Path(second), matrices)
            self.assertEqual(
                (Path(first) / "region_score_matrices.npy").read_bytes(),
                (Path(second) / "region_score_matrices.npy").read_bytes(),
            )
            self.assertEqual(
                (Path(first) / "region_score_matrix_manifest.csv").read_bytes(),
                (Path(second) / "region_score_matrix_manifest.csv").read_bytes(),
            )

    def test_hierarchical_summary_uses_lower_marginal_median(self):
        matrix = np.arange(25, dtype=float).reshape(5, 5) / 24.0
        summary = hierarchical_pair_summary(matrix)
        expected_left = np.median(matrix, axis=1)
        expected_right = np.median(matrix, axis=0)
        self.assertEqual(summary["left_marginal_medians"], expected_left.tolist())
        self.assertEqual(summary["right_marginal_medians"], expected_right.tolist())
        self.assertEqual(
            summary["conservative_score"],
            min(float(np.median(expected_left)), float(np.median(expected_right))),
        )

    def test_gate_fails_on_any_completed_required_rank_above_three(self):
        requirements = []
        for layer in ("pdb", "af_271828", "af_314159"):
            for pair_id in ("HY2", "OB1"):
                for endpoint in ("peptide", "composite", "shape"):
                    requirements.append({
                        "layer": layer,
                        "pair_id": pair_id,
                        "endpoint": endpoint,
                        "status": "complete",
                        "rank": 1,
                        "sensitivity_top3": True,
                        "resampling_top3_fraction": 0.9,
                    })
        requirements[-1]["rank"] = 4
        gate = build_control_gate(requirements)
        self.assertEqual(gate["status"], "fail")
        self.assertFalse(gate["candidate_evaluation_allowed"])
        self.assertTrue(gate["electrostatics_retired_from_candidate_ranking"])

    def test_gate_reports_missing_required_calculation_as_not_evaluable(self):
        gate = build_control_gate([{
            "layer": "af_271828",
            "pair_id": "HY2",
            "endpoint": "peptide",
            "status": "missing",
            "rank": None,
            "sensitivity_top3": False,
            "resampling_top3_fraction": None,
        }])
        self.assertEqual(gate["status"], "not_evaluable")
        self.assertFalse(gate["candidate_evaluation_allowed"])

    def test_gate_invariants_reject_inconsistent_unlock_flags(self):
        requirements = [{
            "layer": "af_271828", "pair_id": "HY2", "endpoint": "peptide",
            "status": "complete", "rank": 4, "sensitivity_top3": True,
            "resampling_top3_fraction": 0.9,
        }]
        gate = build_control_gate(requirements)
        self.assertTrue(_gate_invariants(gate, requirements))
        gate["candidate_evaluation_allowed"] = True
        self.assertFalse(_gate_invariants(gate, requirements))

    def test_readme_table_distinguishes_failure_from_missing_panel(self):
        requirements = [
            {"layer": "pdb", "pair_id": "PAIR_HY2E11_BALF5_MBP", "endpoint": "peptide", "status": "complete", "rank": 5, "sensitivity_top3": False, "resampling_top3_fraction": 1.0},
            {"layer": "pdb", "pair_id": "PAIR_HY2E11_BALF5_MBP", "endpoint": "composite", "status": "missing", "rank": None, "sensitivity_top3": False, "resampling_top3_fraction": 1.0},
        ]
        table = _readme_rank_table(requirements)
        self.assertIn("rank 5 (fail)", table)
        self.assertIn("not evaluable", table)

    def test_finder_suffix_duplicate_detection_is_narrow(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "normal_2026.csv").write_text("ok")
            (root / "stale 2.csv").write_text("duplicate")
            self.assertEqual(_finder_duplicate_paths(root), ["stale 2.csv"])

    def test_dependency_free_svg_figure_is_deterministic(self):
        rows = [
            {"layer": "pdb", "pair_id": "PAIR_HY2", "endpoint": "peptide", "rank": 4},
        ]
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            self.assertEqual(_write_svg_rank_figure(Path(first), rows), "complete_svg_fallback")
            self.assertEqual(_write_svg_rank_figure(Path(second), rows), "complete_svg_fallback")
            first_svg = Path(first) / "figures/primary_control_ranks.svg"
            second_svg = Path(second) / "figures/primary_control_ranks.svg"
            self.assertEqual(first_svg.read_bytes(), second_svg.read_bytes())
            self.assertIn(b"raw rank above 3", first_svg.read_bytes())

    def test_figure_writer_ignores_csv_empty_missing_rank(self):
        rows = [
            {"layer": "pdb", "pair_id": "PAIR_HY2", "endpoint": "peptide", "rank": "4"},
            {"layer": "af_271828", "pair_id": "PAIR_OB1", "endpoint": "shape", "rank": ""},
        ]
        with tempfile.TemporaryDirectory() as directory:
            self.assertIn(_write_figures(Path(directory), rows), {"complete", "complete_svg_fallback"})

    def test_candidate_and_discovery_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            allowed = root / "processed/hla2_positive_control_benchmark_v2_results_2026-08-26"
            forbidden = root / "processed/high_yield_candidate_evidence_2026-08-28"
            allowed.mkdir(parents=True)
            forbidden.mkdir(parents=True)
            validate_control_only_paths([allowed], root)
            with self.assertRaisesRegex(ValueError, "candidate/discovery"):
                validate_control_only_paths([forbidden], root)


if __name__ == "__main__":
    unittest.main()
