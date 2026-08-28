import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hla2_positive_control_benchmark import (  # noqa: E402
    PmhcGeometry,
    _evaluate_weights,
    _sidechain_vector,
    aligned_chain_ca,
    build_af3_job_batches,
    build_pdb_oracle_pairings,
    build_trust_gate,
    generate_weight_grid,
    leave_one_system_out,
    pair_features,
    physicochemical_mismatch,
    rank_feature_percentiles,
    select_score_blind_comparators,
    validate_comparator_registry,
    validate_control_registry,
    validate_af3_job_package,
)
from build_hla2_positive_control_benchmark import (  # noqa: E402
    build_evaluation_skeleton,
    curated_registry,
)


FEATURES = (
    "exposed_ca_rmsd_A",
    "exposed_sidechain_vector_rmsd_A",
    "tcr_face_physicochemical_mismatch",
    "anchor_ca_rmsd_A",
)


def registry_rows():
    systems = [
        {"system_id": "HY2", "tcr_id": "Hy.2E11", "eligibility": "strict", "independent_system_weight": 1},
        {"system_id": "OB1", "tcr_id": "Ob.1A12", "eligibility": "strict", "independent_system_weight": 1},
        {"system_id": "HY1", "tcr_id": "Hy.1B11", "eligibility": "strict", "independent_system_weight": 1},
        {"system_id": "ANO2", "tcr_id": "study clone set", "eligibility": "prospective", "independent_system_weight": 0},
    ]
    ligands = [
        {"ligand_id": "HY2_V", "system_id": "HY2", "sequence": "TGGVYHFVKKHVHES", "core": "YHFVKKHVH", "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB5*01:01", "pdb_id": "1H15"},
        {"ligand_id": "HY2_S", "system_id": "HY2", "sequence": "ENPVVHFFKNIVTPR", "core": "VHFFKNIVT", "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01", "pdb_id": "1BX2"},
        {"ligand_id": "OB1_V", "system_id": "OB1", "sequence": "DFARVHFISALHGSG", "core": "VHFISALHG", "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01", "pdb_id": "2WBJ"},
        {"ligand_id": "OB1_S", "system_id": "OB1", "sequence": "ENPVVHFFKNIVTPR", "core": "VHFFKNIVT", "mhc_alpha_allele": "HLA-DRA*01:01", "mhc_beta_allele": "HLA-DRB1*15:01", "pdb_id": "1YMM"},
        {"ligand_id": "HY1_U", "system_id": "HY1", "sequence": "QLVHFVRDFAQL", "core": "VHFVRDFAQ", "mhc_alpha_allele": "HLA-DQA1*01:02", "mhc_beta_allele": "HLA-DQB1*05:02", "pdb_id": "4MAY"},
        {"ligand_id": "HY1_P", "system_id": "HY1", "sequence": "RLLMLFAKDVVSRN", "core": "MLFAKDVVS", "mhc_alpha_allele": "HLA-DQA1*01:02", "mhc_beta_allele": "HLA-DQB1*05:02", "pdb_id": "4GRL"},
        {"ligand_id": "HY1_S", "system_id": "HY1", "sequence": "ENPVVHFFKNIVTPR", "core": "VHFFKNIVT", "mhc_alpha_allele": "HLA-DQA1*01:02", "mhc_beta_allele": "HLA-DQB1*05:02", "pdb_id": "3PL6"},
    ]
    pairs = [
        {"pair_id": "HY2_PAIR", "system_id": "HY2", "left_ligand_id": "HY2_V", "right_ligand_id": "HY2_S", "required_for_system_pass": True},
        {"pair_id": "OB1_PAIR", "system_id": "OB1", "left_ligand_id": "OB1_V", "right_ligand_id": "OB1_S", "required_for_system_pass": True},
        {"pair_id": "HY1_UL15", "system_id": "HY1", "left_ligand_id": "HY1_U", "right_ligand_id": "HY1_S", "required_for_system_pass": True},
        {"pair_id": "HY1_PMM", "system_id": "HY1", "left_ligand_id": "HY1_P", "right_ligand_id": "HY1_S", "required_for_system_pass": True},
    ]
    return systems, ligands, pairs


class RegistryTests(unittest.TestCase):
    def test_curated_registry_locks_three_systems_four_pairs_and_exact_structures(self):
        systems, ligands, pairs, sources = curated_registry()
        summary = validate_control_registry(systems, ligands, pairs)
        self.assertEqual(summary["strict_independent_system_count"], 3)
        self.assertEqual(summary["strict_positive_pair_count"], 4)
        self.assertEqual({row["pdb_id"] for row in ligands}, {"1H15", "1BX2", "2WBJ", "1YMM", "3PL6", "4MAY", "4GRL"})
        self.assertTrue(all(row["functional_evidence_for_both_arms"] for row in systems if row["eligibility"] == "strict"))
        self.assertGreaterEqual(len(sources), 10)

    def test_registry_counts_tcr_systems_not_positive_pairs(self):
        summary = validate_control_registry(*registry_rows())
        self.assertEqual(summary["strict_independent_system_count"], 3)
        self.assertEqual(summary["strict_positive_pair_count"], 4)
        self.assertEqual(summary["prospective_system_ids"], ["ANO2"])

    def test_registry_rejects_inexact_core_or_unknown_ligand(self):
        systems, ligands, pairs = registry_rows()
        ligands[0]["core"] = "NOTINPEPT"
        with self.assertRaisesRegex(ValueError, "core"):
            validate_control_registry(systems, ligands, pairs)


class FeatureTests(unittest.TestCase):
    def test_chain_alignment_uses_reference_positions_despite_missing_prefix(self):
        reference = "XX" + "ACDEFGHIKLMNPQRSTVWY" * 5
        observed = reference[2:]
        residues = []
        for index, aa in enumerate(observed):
            residues.append({
                "aa": aa,
                "atoms": [{"name": "CA", "element": "C", "xyz": (float(index), 0.0, 0.0)}],
            })
        coordinates = aligned_chain_ca(residues, reference, reference_start=2, count=85)
        self.assertEqual(coordinates.shape, (85, 3))
        self.assertEqual(coordinates[0, 0], 0.0)
        self.assertEqual(coordinates[-1, 0], 84.0)

    def test_chain_alignment_is_deterministic_and_rejects_missing_groove_residue(self):
        reference = "ACDEFGHIKLMNPQRSTVWY" * 5
        residues = [{
            "aa": aa,
            "atoms": [{"name": "CA", "element": "C", "xyz": (float(index), 0.0, 0.0)}],
        } for index, aa in enumerate(reference)]
        first = aligned_chain_ca(residues, reference, reference_start=0, count=85)
        second = aligned_chain_ca(residues, reference, reference_start=0, count=85)
        np.testing.assert_array_equal(first, second)
        with self.assertRaisesRegex(ValueError, "missing"):
            aligned_chain_ca(residues[:40] + residues[41:], reference, reference_start=0, count=85)

    def test_glycine_sidechain_fallback_is_zero_vector(self):
        residue = {
            "aa": "G",
            "atoms": [
                {"name": "N", "element": "N", "xyz": (0.0, 0.0, 0.0)},
                {"name": "CA", "element": "C", "xyz": (1.0, 2.0, 3.0)},
                {"name": "C", "element": "C", "xyz": (2.0, 2.0, 3.0)},
            ],
        }
        np.testing.assert_array_equal(_sidechain_vector(residue), np.zeros(3))

    def test_physicochemical_mismatch_is_zero_for_identical_core(self):
        self.assertEqual(physicochemical_mismatch("VVHFFKNIV", "VVHFFKNIV"), 0.0)
        self.assertGreater(physicochemical_mismatch("VVHFFKNIV", "VVDFFDNIV"), 0.0)

    def test_geometry_features_are_rigid_body_invariant(self):
        rng = np.random.default_rng(41)
        groove = rng.normal(size=(170, 3))
        core_ca = rng.normal(size=(9, 3))
        vectors = rng.normal(size=(9, 3))
        theta = 0.37
        rotation = np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ])
        translation = np.array([4.0, -3.0, 8.0])
        left = PmhcGeometry("L", "VVHFFKNIV", groove, core_ca, vectors)
        right = PmhcGeometry(
            "R", "VVHFFKNIV", groove @ rotation + translation,
            core_ca @ rotation + translation, vectors @ rotation,
        )
        result = pair_features(left, right)
        for feature in FEATURES:
            self.assertAlmostEqual(result[feature], 0.0, places=9)

    def test_feature_percentiles_use_average_ties(self):
        rows = [
            {"pair_id": "a", "x": 1.0},
            {"pair_id": "b", "x": 1.0},
            {"pair_id": "c", "x": 3.0},
        ]
        ranked = rank_feature_percentiles(rows, ("x",))
        by_id = {row["pair_id"]: row for row in ranked}
        self.assertEqual(by_id["a"]["x_percentile"], 0.25)
        self.assertEqual(by_id["b"]["x_percentile"], 0.25)
        self.assertEqual(by_id["c"]["x_percentile"], 1.0)


class SelectionAndValidationTests(unittest.TestCase):
    def test_pdb_oracle_uses_all_exact_hla_pairings_and_requires_five_decoys(self):
        structural_ligands = [
            {"ligand_id": "L1", "mhc_alpha_allele": "A", "mhc_beta_allele": "B1"},
            {"ligand_id": "L2", "mhc_alpha_allele": "A", "mhc_beta_allele": "B1"},
            {"ligand_id": "R1", "mhc_alpha_allele": "A", "mhc_beta_allele": "B2"},
            {"ligand_id": "R2", "mhc_alpha_allele": "A", "mhc_beta_allele": "B2"},
            {"ligand_id": "R3", "mhc_alpha_allele": "A", "mhc_beta_allele": "B2"},
        ]
        cross_hla = {
            "pair_id": "CROSS", "system_id": "SYS", "left_ligand_id": "L1",
            "right_ligand_id": "R1", "left_mhc_alpha_allele": "A",
            "left_mhc_beta_allele": "B1", "right_mhc_alpha_allele": "A",
            "right_mhc_beta_allele": "B2",
        }
        rows, summary = build_pdb_oracle_pairings(cross_hla, structural_ligands)
        self.assertEqual(len(rows), 6)
        self.assertEqual(summary["decoy_count"], 5)
        self.assertEqual(summary["evaluation_status"], "complete")

        same_hla = {
            **cross_hla, "pair_id": "SAME", "left_ligand_id": "R1",
            "right_ligand_id": "R2", "left_mhc_beta_allele": "B2",
        }
        rows, summary = build_pdb_oracle_pairings(same_hla, structural_ligands)
        self.assertEqual(len(rows), 3)
        self.assertEqual(summary["decoy_count"], 2)
        self.assertEqual(summary["evaluation_status"], "not_evaluable_insufficient_exact_hla_decoys")

    def test_comparator_registry_requires_five_exact_registers_per_arm(self):
        rows = []
        for arm in ("microbial", "self"):
            for index in range(5):
                rows.append({
                    "positive_pair_id": "PAIR",
                    "comparator_arm": arm,
                    "candidate_id": f"{arm}_{index}",
                    "sequence": "ACDEFGHIKLMNPQR",
                    "predicted_core": "DEFGHIKLM",
                    "core_start_1_based": 3,
                    "register_resolution": "resolved_unique_fully_contained",
                    "seq_num": index + 1,
                    "raw_response_file": "raw.tsv",
                    "negative_tier": "N3",
                    "recognition_status": "unknown_not_specificity_negative",
                    "selection_is_score_blind": True,
                })
        summary = validate_comparator_registry(rows, expected_pair_ids=("PAIR",))
        self.assertEqual(summary["comparison_pair_count"], 1)
        self.assertEqual(summary["unique_comparator_count"], 10)
        rows[0]["predicted_core"] = "NOTACOREX"
        with self.assertRaisesRegex(ValueError, "core"):
            validate_comparator_registry(rows, expected_pair_ids=("PAIR",))

    def test_evaluation_skeleton_requires_pdb_and_both_seeds_for_every_pair(self):
        _, _, pairs, _ = curated_registry()
        rows = build_evaluation_skeleton(pairs, panel_seeds=(104729, 104759))
        self.assertEqual(len(rows), 12)
        for pair_id in {row["pair_id"] for row in rows}:
            pair_rows = [row for row in rows if row["pair_id"] == pair_id]
            self.assertEqual(
                {(row["layer"], str(row["panel_seed"])) for row in pair_rows},
                {("pdb_oracle", "pdb"), ("af3", "104729"), ("af3", "104759")},
            )
            self.assertTrue(all(row["evaluation_status"] == "pending" for row in pair_rows))

    def test_weight_grid_is_nonnegative_and_sums_to_one(self):
        grid = generate_weight_grid(FEATURES)
        self.assertEqual(len(grid), 35)
        self.assertTrue(all(abs(sum(row.values()) - 1.0) < 1e-12 for row in grid))
        self.assertTrue(all(value >= 0 for row in grid for value in row.values()))

    def test_score_blind_selection_ignores_geometry_fields(self):
        target = {"sequence": "ABCDEFGHIJKLMNO", "binding_percentile": 2.0}
        candidates = [
            {"candidate_id": "z", "sequence": "A" * 15, "binding_percentile": 2.0, "geometry_score": 0.0},
            {"candidate_id": "a", "sequence": "C" * 15, "binding_percentile": 2.0, "geometry_score": 99.0},
        ]
        selected = select_score_blind_comparators(target, candidates, count=2, seed=104759)
        reversed_scores = [dict(row, geometry_score=100 - row["geometry_score"]) for row in candidates]
        selected_again = select_score_blind_comparators(target, reversed_scores, count=2, seed=104759)
        self.assertEqual([row["candidate_id"] for row in selected], [row["candidate_id"] for row in selected_again])

    def test_leave_one_system_out_never_trains_on_held_out_system(self):
        rows = []
        for system_index, system_id in enumerate(("HY2", "OB1", "HY1"), start=1):
            for pair_index in range(4):
                rows.append({
                    "system_id": system_id,
                    "pair_id": f"{system_id}_{pair_index}",
                    "pair_role": "positive" if pair_index == 0 else "N3",
                    "layer": "pdb_oracle",
                    "panel_seed": "pdb",
                    **{feature: float(pair_index + system_index) for feature in FEATURES},
                })
        results = leave_one_system_out(rows, FEATURES)
        self.assertEqual({row["held_out_system_id"] for row in results}, {"HY2", "OB1", "HY1"})
        for row in results:
            self.assertNotIn(row["held_out_system_id"], row["training_system_ids"].split(";"))

    def test_multi_ligand_panels_are_ranked_separately_with_one_system_vote(self):
        rows = []
        for positive_pair_id in ("UL15", "PMM"):
            for pair_index in range(4):
                rows.append({
                    "system_id": "HY1",
                    "positive_pair_id": positive_pair_id,
                    "pair_id": f"{positive_pair_id}_{pair_index}",
                    "pair_role": "positive" if pair_index == 0 else "N3",
                    "layer": "af3",
                    "panel_seed": "104729",
                    **{feature: float(pair_index) for feature in FEATURES},
                })
        weights = {feature: float(feature == "exposed_ca_rmsd_A") for feature in FEATURES}
        results = _evaluate_weights(rows, FEATURES, weights)
        self.assertEqual(len(results), 2)
        self.assertEqual({row["positive_pair_id"] for row in results}, {"UL15", "PMM"})
        self.assertTrue(all(row["comparison_count"] == 4 for row in results))

    def test_missing_required_seed_blocks_overall_pass(self):
        results = [
            {"system_id": "HY2", "pair_id": "HY2_PAIR", "layer": "pdb_oracle", "panel_seed": "pdb", "evaluation_status": "complete", "positive_rank": 1},
            {"system_id": "HY2", "pair_id": "HY2_PAIR", "layer": "af3", "panel_seed": "104759", "evaluation_status": "complete", "positive_rank": 1},
            {"system_id": "HY2", "pair_id": "HY2_PAIR", "layer": "af3", "panel_seed": "104729", "evaluation_status": "missing", "positive_rank": ""},
        ]
        gate = build_trust_gate(results, required_system_ids=("HY2",))
        self.assertEqual(gate["overall_trust_status"], "not_evaluable")
        self.assertEqual(gate["failed_system_count"], 0)

    def test_multi_ligand_system_gets_one_vote_and_every_positive_is_required(self):
        results = []
        for pair_id, rank in (("UL15", 1), ("PMM", 3)):
            for layer, seed in (("pdb_oracle", "pdb"), ("af3", "104729"), ("af3", "104759")):
                results.append({
                    "system_id": "HY1", "pair_id": pair_id, "layer": layer,
                    "panel_seed": seed, "evaluation_status": "complete", "positive_rank": rank,
                })
        gate = build_trust_gate(results, required_system_ids=("HY1",))
        self.assertEqual(gate["passed_system_count"], 1)
        self.assertEqual(len(gate["system_statuses"]), 1)
        results[-1]["positive_rank"] = 4
        gate = build_trust_gate(results, required_system_ids=("HY1",))
        self.assertEqual(gate["overall_trust_status"], "fail")


class AlphaFoldJobTests(unittest.TestCase):
    def test_jobs_use_generic_alpha_beta_peptide_roles_and_batches_of_30(self):
        ligands = [
            {
                "ligand_id": f"L{index:02d}",
                "sequence": "ACDEFGHIKLMNPQR",
                "mhc_alpha_allele": "HLA-DQA1*01:02",
                "mhc_beta_allele": "HLA-DQB1*05:02",
            }
            for index in range(31)
        ]
        hla = {
            ("HLA-DQA1*01:02", "HLA-DQB1*05:02"): {
                "mhc_alpha_sequence": "A" * 180,
                "mhc_beta_sequence": "B" * 190,
            }
        }
        jobs, manifest, batches = build_af3_job_batches(ligands, hla, panel_seeds=(104729,))
        self.assertEqual([len(batch) for batch in batches], [30, 1])
        self.assertEqual(manifest[0]["chain_roles"], "mhc_alpha;mhc_beta;peptide")
        self.assertEqual(len(jobs[0]["sequences"]), 3)
        self.assertEqual(jobs[0]["modelSeeds"], [104729])

    def test_job_package_validator_checks_chain_order_exact_sequence_and_status(self):
        ligands = [{
            "ligand_id": "L", "system_id": "S", "ligand_role": "positive",
            "sequence": "ACDEFGHIKLMNPQR", "core_sequence": "DEFGHIKLM",
            "core_start_1_based": 3, "register_resolution": "experimentally_resolved",
            "register_source": "PDB_TEST", "mhc_alpha_allele": "A", "mhc_beta_allele": "B",
        }]
        hla = {("A", "B"): {"mhc_alpha_sequence": "A" * 180, "mhc_beta_sequence": "B" * 190}}
        jobs, manifest, batches = build_af3_job_batches(ligands, hla, panel_seeds=(1,))
        summary = validate_af3_job_package(batches, manifest, hla)
        self.assertEqual(summary["job_count"], 1)
        self.assertEqual(summary["batch_sizes"], [1])
        batches[0][0]["sequences"].reverse()
        with self.assertRaisesRegex(ValueError, "chain"):
            validate_af3_job_package(batches, manifest, hla)


if __name__ == "__main__":
    unittest.main()
