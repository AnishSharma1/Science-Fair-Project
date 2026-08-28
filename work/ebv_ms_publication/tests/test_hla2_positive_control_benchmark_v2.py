import unittest
import csv
import hashlib
import json
import tempfile
from pathlib import Path


from hla2_positive_control_benchmark_v2 import (
    PILOT_SEEDS,
    aggregate_system_results,
    blosum62_similarity,
    build_definitive_ranking_gate,
    build_oracle_availability,
    build_pilot_attribution_gate,
    build_protocol_lock,
    validate_specificity_registry,
)
from build_hla2_positive_control_benchmark_v2 import build_pilot_package
from run_hla2_positive_control_benchmark_v2 import (
    METHOD_NAMES,
    evaluate_outer_folds,
    run_v2_analysis,
    select_strongest_nonstructural_baseline,
)


class ProtocolLockTests(unittest.TestCase):
    def test_protocol_is_fresh_frozen_and_never_submitted(self):
        protocol = build_protocol_lock(
            strict_system_ids=("HY2", "OB1", "HY1"),
            positive_pair_ids=("HY2_PAIR", "OB1_PAIR", "HY1_A", "HY1_B"),
            registry_sha256="a" * 64,
            comparator_sha256="b" * 64,
            software_versions={"python": "test"},
        )
        self.assertEqual(PILOT_SEEDS, (271828, 314159))
        self.assertEqual(protocol["status"], "prepared_not_submitted")
        self.assertEqual(protocol["benchmark_stage"], "three_system_attribution_pilot")
        self.assertFalse(protocol["reuses_v1_jobs"])
        self.assertIn("python", protocol["software_versions"])
        self.assertEqual(protocol["comparators_per_arm"], 5)
        self.assertEqual(protocol["pair_decoys_per_panel"], 25)
        self.assertEqual(protocol["tcr_facing_positions"], ["P2", "P3", "P5", "P7", "P8"])
        self.assertEqual(protocol["anchor_positions"], ["P1", "P4", "P6", "P9"])


class BaselineTests(unittest.TestCase):
    def test_blosum_can_be_scoped_to_tcr_face_or_full_core(self):
        self.assertGreater(
            blosum62_similarity("VHFFKNIVT", "VHFIKNIVT"),
            blosum62_similarity("VHFFKNIVT", "AAAAAAAAA"),
        )
        self.assertEqual(
            blosum62_similarity("ACDEFGHIK", "ACDEFGHIK", positions=(1, 2, 4, 6, 7)),
            blosum62_similarity("ACDEFGHIK", "ACDEFGHIK"),
        )


class AggregationTests(unittest.TestCase):
    def test_multi_ligand_system_uses_worst_required_rank_and_one_vote(self):
        rows = [
            {"system_id": "HY1", "positive_pair_id": "UL15", "panel_seed": 271828,
             "evaluation_status": "complete", "composite_rank": 1, "baseline_rank": 2,
             "structural_weight": 0.5, "ablated_rank": 2},
            {"system_id": "HY1", "positive_pair_id": "PMM", "panel_seed": 271828,
             "evaluation_status": "complete", "composite_rank": 3, "baseline_rank": 4,
             "structural_weight": 0.5, "ablated_rank": 4},
            {"system_id": "HY1", "positive_pair_id": "UL15", "panel_seed": 314159,
             "evaluation_status": "complete", "composite_rank": 2, "baseline_rank": 2,
             "structural_weight": 0.5, "ablated_rank": 2},
            {"system_id": "HY1", "positive_pair_id": "PMM", "panel_seed": 314159,
             "evaluation_status": "complete", "composite_rank": 3, "baseline_rank": 4,
             "structural_weight": 0.5, "ablated_rank": 4},
        ]
        systems = aggregate_system_results(rows, required_seeds=PILOT_SEEDS)
        self.assertEqual(len(systems), 1)
        self.assertEqual(systems[0]["system_score"], 3)
        self.assertEqual(systems[0]["baseline_system_score"], 4)
        self.assertEqual(systems[0]["independent_system_vote"], 1)

    def test_every_credited_panel_improvement_must_be_removed_by_ablation(self):
        rows = [
            {"system_id": "HY1", "positive_pair_id": pair, "panel_seed": seed,
             "evaluation_status": "complete", "composite_rank": 1, "baseline_rank": 3,
             "structural_weight": 0.5, "ablated_rank": 3}
            for pair in ("UL15", "PMM") for seed in PILOT_SEEDS
        ]
        rows[-1]["ablated_rank"] = 2
        systems = aggregate_system_results(
            rows,
            required_seeds=PILOT_SEEDS,
            required_pairs_by_system={"HY1": ("UL15", "PMM")},
        )
        self.assertFalse(systems[0]["structural_ablation_removes_improvement"])


class PilotGateTests(unittest.TestCase):
    @staticmethod
    def supportive_rows():
        return [
            {"system_id": "HY2", "evaluation_status": "complete", "system_score": 2,
             "baseline_system_score": 3, "minimum_structural_weight_on_improvements": 0.25,
             "structural_ablation_removes_improvement": True},
            {"system_id": "OB1", "evaluation_status": "complete", "system_score": 1,
             "baseline_system_score": 2, "minimum_structural_weight_on_improvements": 0.5,
             "structural_ablation_removes_improvement": True},
            {"system_id": "HY1", "evaluation_status": "complete", "system_score": 3,
             "baseline_system_score": 3, "minimum_structural_weight_on_improvements": 0.0,
             "structural_ablation_removes_improvement": True},
        ]

    def test_majority_better_none_worse_is_supportive_but_never_unlocks(self):
        gate = build_pilot_attribution_gate(
            self.supportive_rows(), required_system_ids=("HY2", "OB1", "HY1")
        )
        self.assertEqual(gate["pilot_attribution_status"], "supportive")
        self.assertFalse(gate["weights_frozen"])
        self.assertFalse(gate["discovery_unlock_allowed"])

    def test_any_worsened_system_fails(self):
        rows = self.supportive_rows()
        rows[-1]["system_score"] = 4
        gate = build_pilot_attribution_gate(rows, required_system_ids=("HY2", "OB1", "HY1"))
        self.assertEqual(gate["pilot_attribution_status"], "fail")

    def test_zero_structural_weight_cannot_support_attribution(self):
        rows = self.supportive_rows()
        rows[0]["minimum_structural_weight_on_improvements"] = 0.0
        gate = build_pilot_attribution_gate(rows, required_system_ids=("HY2", "OB1", "HY1"))
        self.assertEqual(gate["pilot_attribution_status"], "fail")

    def test_incomplete_system_is_not_evaluable(self):
        rows = self.supportive_rows()[:-1]
        gate = build_pilot_attribution_gate(rows, required_system_ids=("HY2", "OB1", "HY1"))
        self.assertEqual(gate["pilot_attribution_status"], "not_evaluable")

    def test_mandatory_pdb_oracle_pending_blocks_supportive_alpha_result(self):
        gate = build_pilot_attribution_gate(
            self.supportive_rows(),
            required_system_ids=("HY2", "OB1", "HY1"),
            oracle_rows=[{"mandatory_if_scored": True, "oracle_status": "required_pending_results"}],
        )
        self.assertEqual(gate["pilot_attribution_status"], "not_evaluable")

    def test_unavailable_pdb_oracle_does_not_block_supportive_alpha_result(self):
        gate = build_pilot_attribution_gate(
            self.supportive_rows(),
            required_system_ids=("HY2", "OB1", "HY1"),
            oracle_rows=[{"mandatory_if_scored": False, "oracle_status": "not_evaluable_availability"}],
        )
        self.assertEqual(gate["pilot_attribution_status"], "supportive")


class OracleAndSpecificityTests(unittest.TestCase):
    def test_unavailable_oracle_is_not_a_pass_and_does_not_block_af(self):
        availability = build_oracle_availability([
            {"positive_pair_id": "DQ_PAIR", "eligible_decoy_count": 2},
            {"positive_pair_id": "DR_PAIR", "eligible_decoy_count": 5},
        ])
        by_pair = {row["positive_pair_id"]: row for row in availability}
        self.assertEqual(by_pair["DQ_PAIR"]["oracle_status"], "not_evaluable_availability")
        self.assertFalse(by_pair["DQ_PAIR"]["mandatory_if_scored"])
        self.assertEqual(by_pair["DR_PAIR"]["oracle_status"], "required_pending_results")
        self.assertTrue(by_pair["DR_PAIR"]["mandatory_if_scored"])

    def test_only_explicit_n1_n2_negatives_enter_specificity_gate(self):
        gate = validate_specificity_registry([
            {"negative_id": "N1_OK", "negative_tier": "N1", "peptide": "ABCDEFGHI",
             "exact_hla": "DRA*01:01/DRB1*15:01", "assay": "IL-2", "tested_condition": "10 uM",
             "outcome": "no response", "source_location": "Table S1"},
            {"negative_id": "N3", "negative_tier": "N3", "peptide": "JKLMNOPQR"},
        ])
        self.assertEqual(gate["specificity_status"], "prepared")
        self.assertEqual(gate["admitted_negative_count"], 1)
        self.assertEqual(gate["excluded_n3_count"], 1)


class DefinitiveGateTests(unittest.TestCase):
    def test_fewer_than_six_systems_is_blocked_even_if_pilot_is_supportive(self):
        systems = PilotGateTests.supportive_rows()
        registry = [
            {"system_id": row["system_id"], "eligibility": "strict", "hla_family": "DR"}
            for row in systems
        ]
        gate = build_definitive_ranking_gate(systems, registry)
        self.assertEqual(gate["definitive_status"], "blocked_registry_size")
        self.assertFalse(gate["weights_frozen"])
        self.assertFalse(gate["discovery_unlock_allowed"])

    def test_six_systems_two_families_majority_better_none_worse_can_pass(self):
        systems = []
        registry = []
        for index in range(6):
            system_id = f"SYS{index}"
            systems.append({
                "system_id": system_id,
                "evaluation_status": "complete",
                "system_score": 1 if index < 4 else 2,
                "baseline_system_score": 2,
                "minimum_structural_weight_on_improvements": 0.25 if index < 4 else 0.0,
                "structural_ablation_removes_improvement": True,
            })
            registry.append({
                "system_id": system_id,
                "eligibility": "strict",
                "hla_family": "DR" if index < 3 else "DQ",
                "distinct_biological_sources_verified": True,
            })
        gate = build_definitive_ranking_gate(systems, registry)
        self.assertEqual(gate["definitive_status"], "pass")
        self.assertTrue(gate["weights_frozen"])
        self.assertTrue(gate["discovery_unlock_allowed"])


class PilotPackageTests(unittest.TestCase):
    @staticmethod
    def tree_digest(root: Path) -> str:
        digest = hashlib.sha256()
        for path in sorted(value for value in root.rglob("*") if value.is_file()):
            if path.name == "SHA256SUMS.csv":
                continue
            digest.update(str(path.relative_to(root)).encode("utf-8"))
            digest.update(b"\0")
            digest.update(path.read_bytes())
        return digest.hexdigest()

    def test_builder_prepares_fresh_complete_pilot_without_unlocking(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "pilot"
            result = build_pilot_package(out)
            protocol = json.loads((out / "protocol/protocol_lock.json").read_text())
            pilot_gate = json.loads((out / "benchmark/pilot_attribution_gate.json").read_text())
            definitive_gate = json.loads((out / "benchmark/definitive_ranking_gate.json").read_text())
            specificity = json.loads((out / "benchmark/specificity_gate.json").read_text())
            with (out / "controls/control_decoy_registry.csv").open() as handle:
                comparators = list(csv.DictReader(handle))
            with (out / "controls/comparison_universe.csv").open() as handle:
                comparisons = list(csv.DictReader(handle))
            with (out / "alphafold_jobs/job_manifest.csv").open() as handle:
                jobs = list(csv.DictReader(handle))
            oracle_pairings_exist = (out / "benchmark/pdb_oracle_frozen_pairings.csv").exists()
        self.assertEqual(result["strict_system_count"], 3)
        self.assertEqual(result["positive_pair_count"], 4)
        self.assertEqual(len(comparators), 40)
        self.assertEqual(len(comparisons), 208)
        self.assertTrue(all(row["status"] == "prepared_not_submitted" for row in jobs))
        self.assertTrue(all(size <= 30 for size in result["batch_sizes"]))
        self.assertFalse(protocol["reuses_v1_jobs"])
        self.assertIn("python", protocol["software_versions"])
        self.assertTrue(oracle_pairings_exist)
        self.assertEqual(pilot_gate["pilot_attribution_status"], "not_evaluable")
        self.assertEqual(definitive_gate["definitive_status"], "blocked_registry_size")
        self.assertEqual(specificity["specificity_status"], "not_evaluable_no_verified_negatives")
        self.assertFalse(pilot_gate["discovery_unlock_allowed"])

    def test_rebuild_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first"
            second = Path(temporary) / "second"
            build_pilot_package(first)
            build_pilot_package(second)
            self.assertEqual(self.tree_digest(first), self.tree_digest(second))


class V2AnalysisTests(unittest.TestCase):
    @staticmethod
    def synthetic_rows():
        rows = []
        for system_index, system_id in enumerate(("SYS_A", "SYS_B", "SYS_C")):
            for pair_index in range(2):
                pair_id = f"{system_id}_PAIR_{pair_index}"
                for comparison_index in range(6):
                    positive = comparison_index == 0
                    rows.append({
                        "system_id": system_id,
                        "positive_pair_id": pair_id,
                        "panel_seed": 271828 + pair_index,
                        "pair_id": f"{pair_id}_{comparison_index}",
                        "pair_role": "positive" if positive else "N3",
                        "exposed_ca_rmsd_A": 0.1 + comparison_index + system_index * 0.01,
                        "exposed_sidechain_vector_rmsd_A": 0.2 + comparison_index,
                        "tcr_face_physicochemical_mismatch": 0.3 + comparison_index,
                        "anchor_ca_rmsd_A": 0.4 + comparison_index,
                        "tcr_facing_sequence_identity": 1.0 if positive else 0.5,
                        "full_core_sequence_identity": 1.0 if positive else 0.5,
                        "tcr_facing_blosum62_similarity": 1.0 if positive else 0.2,
                        "full_core_blosum62_similarity": 1.0 if positive else 0.2,
                    })
        return rows

    def test_nonstructural_baseline_selection_reads_training_systems_only(self):
        rows = self.synthetic_rows()
        training = [row for row in rows if row["system_id"] != "SYS_C"]
        selected = select_strongest_nonstructural_baseline(training)
        self.assertIn(selected, {
            "physicochemical_only", "tcr_facing_identity", "full_core_identity",
            "tcr_facing_blosum62", "full_core_blosum62",
        })
        folds = evaluate_outer_folds(rows)
        held_out = [row for row in folds["panel_method_ranks"] if row["held_out_system_id"] == "SYS_C"]
        self.assertTrue(held_out)
        self.assertTrue(all("SYS_C" not in row["training_system_ids"].split(";") for row in held_out))

    def test_every_locked_method_is_reported_inside_each_outer_fold(self):
        folds = evaluate_outer_folds(self.synthetic_rows())
        first_panel = folds["panel_method_ranks"][0]
        observed = {
            row["method"] for row in folds["method_rank_long"]
            if row["held_out_system_id"] == first_panel["held_out_system_id"]
            and row["positive_pair_id"] == first_panel["positive_pair_id"]
            and row["panel_seed"] == first_panel["panel_seed"]
        }
        self.assertEqual(observed, set(METHOD_NAMES))

    def test_prepared_package_analysis_is_not_evaluable_without_downloads(self):
        with tempfile.TemporaryDirectory() as temporary:
            package = Path(temporary) / "pilot"
            out = Path(temporary) / "results"
            build_pilot_package(package)
            result = run_v2_analysis([], package=package, out=out)
            pilot_gate = json.loads((out / "benchmark/pilot_attribution_gate.json").read_text())
            manifest = json.loads((out / "analysis_manifest.json").read_text())
            self.assertEqual(result["pilot_attribution_status"], "not_evaluable")
            self.assertEqual(pilot_gate["pilot_attribution_status"], "not_evaluable")
            self.assertEqual(manifest["downloaded_complete_exact_job_count"], 0)
            self.assertFalse(manifest["weights_frozen"])
            self.assertFalse(manifest["discovery_unlock_allowed"])
            self.assertTrue((out / "benchmark/panel_method_ranks.csv").exists())


if __name__ == "__main__":
    unittest.main()
