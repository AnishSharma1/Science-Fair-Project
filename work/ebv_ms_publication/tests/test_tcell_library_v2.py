import copy
import math
import unittest

import numpy as np

from tcell_library_v2 import (
    ALLELES,
    CALIBRATION_SEEDS,
    build_af3_jobs,
    build_calibration_comparison_universe,
    build_discovery_pairs,
    build_native_calibration_jobs,
    build_prediction_manifest,
    classify_positive_recovery,
    freeze_panel,
    groove_fitted_rmsd,
    select_native_controls,
    validate_registry,
)


def candidate(candidate_id, kingdom, protein, sequence, priority=1, required=False):
    return {
        "candidate_id": candidate_id,
        "kingdom": kingdom,
        "protein_symbol": protein,
        "sequence": sequence,
        "sequence_length": len(sequence),
        "evidence_priority": priority,
        "native_hla_evidence": priority == 1,
        "source_certainty": "exact_primary_source" if priority == 1 else "canonical_exploratory",
        "accession": f"ACC_{protein}",
        "start": priority,
        "end": priority + len(sequence) - 1,
        "required_for_confirmed_system": required,
        "proposed_core": sequence[priority - 1:priority + 8],
    }


def synthetic_candidates():
    rows = []
    ebv_proteins = [f"EBV{i:02d}" for i in range(20)]
    self_proteins = [f"SELF{i:02d}" for i in range(11)]
    for p_i, protein in enumerate(ebv_proteins):
        for rank in range(1, 4):
            seq = ("ACDEFGHIKLMNPQRSTVWY" * 2)[p_i % 10 : p_i % 10 + 14 + rank]
            rows.append(candidate(f"E_{p_i:02d}_{rank}", "EBV", protein, seq, rank, p_i == 0 and rank == 1))
    for p_i, protein in enumerate(self_proteins):
        for rank in range(1, 6):
            seq = ("YWVTSRQPNMLKIHGFEDCA" * 2)[p_i % 10 : p_i % 10 + 14 + rank]
            rows.append(candidate(f"S_{p_i:02d}_{rank}", "human_self", protein, seq, rank, p_i == 0 and rank == 1))
    return rows


class RegistryContractTests(unittest.TestCase):
    def test_registry_requires_unique_systems_and_excludes_antibody_from_denominator(self):
        rows = [
            {
                "biological_system_id": "SYS_BALF5_MBP",
                "evidence_tier": "E1_exact_pmhc_positive",
                "receptor_modality": "T_cell",
                "tcell_positive_denominator": True,
                "primary_source": "PMID:12244309",
            },
            {
                "biological_system_id": "SYS_GLIALCAM",
                "evidence_tier": "context_only",
                "receptor_modality": "antibody",
                "tcell_positive_denominator": False,
                "primary_source": "PMID:35073561",
            },
        ]
        summary = validate_registry(rows)
        self.assertEqual(summary["strict_e1_count"], 1)
        self.assertEqual(summary["tcell_denominator_ids"], ["SYS_BALF5_MBP"])

        bad = copy.deepcopy(rows)
        bad[1]["tcell_positive_denominator"] = True
        with self.assertRaisesRegex(ValueError, "antibody"):
            validate_registry(bad)

        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_registry(rows + [copy.deepcopy(rows[0])])


class PanelContractTests(unittest.TestCase):
    def test_panel_is_balanced_diverse_deterministic_and_score_blind(self):
        candidates = synthetic_candidates()
        panel = freeze_panel(candidates)
        self.assertEqual(len(panel), 80)
        self.assertEqual(sum(row["kingdom"] == "EBV" for row in panel), 40)
        self.assertEqual(sum(row["kingdom"] == "human_self" for row in panel), 40)
        self.assertGreaterEqual(len({row["protein_symbol"] for row in panel if row["kingdom"] == "EBV"}), 18)
        self.assertEqual(len({row["protein_symbol"] for row in panel if row["kingdom"] == "human_self"}), 11)
        self.assertIn("E_00_1", {row["candidate_id"] for row in panel})
        self.assertIn("S_00_1", {row["candidate_id"] for row in panel})

        shuffled_scores = copy.deepcopy(candidates)
        for i, row in enumerate(reversed(shuffled_scores)):
            row["structural_score"] = i * 1000
        self.assertEqual(
            [row["candidate_id"] for row in panel],
            [row["candidate_id"] for row in freeze_panel(shuffled_scores)],
        )

        ebv_counts = {}
        self_counts = {}
        for row in panel:
            target = ebv_counts if row["kingdom"] == "EBV" else self_counts
            target[row["protein_symbol"]] = target.get(row["protein_symbol"], 0) + 1
        self.assertLessEqual(max(ebv_counts.values()), 3)
        self.assertLessEqual(max(self_counts.values()), 4)

    def test_same_accession_and_proposed_core_collapses_to_one_panel_primary(self):
        candidates = synthetic_candidates()
        duplicate = copy.deepcopy(candidates[0])
        duplicate["candidate_id"] = "E_DUPLICATE_LENGTH_OR_COORDINATE_RECORD"
        duplicate["sequence"] = duplicate["sequence"] + "A"
        duplicate["sequence_length"] = len(duplicate["sequence"])
        duplicate["start"] = 999
        duplicate["end"] = 999 + len(duplicate["sequence"]) - 1
        panel = freeze_panel(candidates + [duplicate])
        matching = [
            row for row in panel
            if row["accession"] == candidates[0]["accession"]
            and row["proposed_core"] == candidates[0]["proposed_core"]
        ]
        self.assertEqual(len(matching), 1)

    def test_prediction_and_pair_contracts(self):
        panel = freeze_panel(synthetic_candidates())
        predictions = build_prediction_manifest(panel)
        self.assertEqual(len(predictions), 320)
        self.assertEqual({row["allele"] for row in predictions}, set(ALLELES))
        for allele in ALLELES:
            allele_rows = [row for row in predictions if row["allele"] == allele]
            self.assertEqual(len(allele_rows), 80)
            self.assertEqual([row["seq_num"] for row in allele_rows], list(range(1, 81)))

        pairs = build_discovery_pairs(panel)
        self.assertEqual(len(pairs), 6400)
        for allele in ALLELES:
            self.assertEqual(sum(row["allele"] == allele for row in pairs), 1600)
        self.assertTrue(all(row["ebv_allele"] == row["self_allele"] == row["allele"] for row in pairs))


class ModelContractTests(unittest.TestCase):
    def setUp(self):
        self.panel = freeze_panel(synthetic_candidates())
        self.hla = {
            allele: {"dra_sequence": "A" * 178, "drb_sequence": chr(66 + i) * 189}
            for i, allele in enumerate(ALLELES)
        }

    def test_general_jobs_are_320_and_packaged_30_x_10_plus_20(self):
        jobs, inventory, batches = build_af3_jobs(self.panel, self.hla)
        self.assertEqual(len(jobs), 320)
        self.assertEqual(len(inventory), 320)
        self.assertEqual([len(batch) for batch in batches], [30] * 10 + [20])
        self.assertTrue(all(len(job["sequences"]) == 3 for job in jobs))
        self.assertEqual(len({job["name"] for job in jobs}), 320)

    def test_native_controls_are_deterministic_and_ignore_geometry(self):
        target = {"candidate_id": "TARGET", "sequence": "ACDEFGHIKLMNPQR", "binding_rank_bin": "weak"}
        pool = []
        for i in range(12):
            pool.append({
                "candidate_id": f"CTRL_{i:02d}",
                "sequence": ("ACDEYGHIKLMNPQR" if i == 0 else "YYYYYGHIKLMNPQR"),
                "binding_rank_bin": "weak",
                "protein_symbol": f"P{i}",
                "excluded_source": False,
                "geometry_rmsd": 100 - i,
            })
        chosen = select_native_controls(target, pool, count=5)
        changed = copy.deepcopy(pool)
        for i, row in enumerate(changed):
            row["geometry_rmsd"] = i * 100
        self.assertEqual(
            [row["candidate_id"] for row in chosen],
            [row["candidate_id"] for row in select_native_controls(target, changed, count=5)],
        )
        self.assertEqual(len(chosen), 5)
        self.assertTrue(all(abs(len(row["sequence"]) - len(target["sequence"])) <= 1 for row in chosen))

    def test_native_calibration_is_24_fixed_seed_native_hla_jobs(self):
        entities = []
        for arm, allele in (("viral", "HLA-DRB5*01:01"), ("self", "HLA-DRB1*15:01")):
            for i in range(6):
                entities.append({
                    "entity_id": f"{arm}_{i}",
                    "entity_role": f"E1_positive_{arm}" if i == 0 else "control",
                    "arm": arm,
                    "allele": allele,
                    "sequence": "ACDEFGHIKLMNPQR",
                    "dra_sequence": "A" * 178,
                    "drb_sequence": "B" * 189,
                })
        jobs, manifest = build_native_calibration_jobs(entities)
        self.assertEqual(len(jobs), 24)
        self.assertEqual(len(manifest), 24)
        self.assertEqual({job["modelSeeds"][0] for job in jobs}, set(CALIBRATION_SEEDS))
        self.assertEqual(sum(row["allele"] == "HLA-DRB5*01:01" for row in manifest), 12)
        self.assertEqual(sum(row["allele"] == "HLA-DRB1*15:01" for row in manifest), 12)

        comparisons = build_calibration_comparison_universe(manifest)
        self.assertEqual(len(comparisons), 72)
        for seed in CALIBRATION_SEEDS:
            seed_rows = [row for row in comparisons if row["seed"] == seed]
            self.assertEqual(sum(row["analysis_set"] == "primary_rank_of_26" for row in seed_rows), 26)
            self.assertEqual(sum(row["analysis_set"] == "single_arm_sensitivity" for row in seed_rows), 10)
            self.assertEqual(sum(row["pair_role"] == "E1_positive" for row in seed_rows), 1)


class RecoveryAndGeometryTests(unittest.TestCase):
    def test_recovery_requires_top_three_both_seeds_and_below_control_median(self):
        good = [
            {"seed": 104729, "rank_of_26": 2, "positive_rmsd": 0.7, "equal_weight_control_median": 1.0},
            {"seed": 104759, "rank_of_26": 3, "positive_rmsd": 0.8, "equal_weight_control_median": 1.1},
        ]
        self.assertEqual(classify_positive_recovery(good), "recovered")
        bad = copy.deepcopy(good)
        bad[1]["rank_of_26"] = 4
        self.assertEqual(classify_positive_recovery(bad), "failed_calibration")
        self.assertEqual(classify_positive_recovery([]), "not_evaluable")

    def test_cross_hla_groove_fit_is_rigid_body_invariant(self):
        groove = np.array([[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]], dtype=float)
        core = np.array([[i, math.sin(i), math.cos(i)] for i in range(9)], dtype=float)
        theta = math.pi / 3
        rotation = np.array([
            [math.cos(theta), -math.sin(theta), 0],
            [math.sin(theta), math.cos(theta), 0],
            [0, 0, 1],
        ])
        translation = np.array([11.0, -7.0, 4.5])
        moved_groove = groove @ rotation.T + translation
        moved_core = core @ rotation.T + translation
        result = groove_fitted_rmsd(groove, core, moved_groove, moved_core)
        self.assertAlmostEqual(result["full_core_rmsd"], 0.0, places=7)
        self.assertAlmostEqual(result["exposed_rmsd"], 0.0, places=7)


if __name__ == "__main__":
    unittest.main()
