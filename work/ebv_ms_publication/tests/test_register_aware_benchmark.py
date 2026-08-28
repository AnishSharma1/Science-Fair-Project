"""Tests for register-aware, score-blind benchmark helpers."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from register_aware_benchmark import (  # noqa: E402
    is_assessable_same_register_pair,
    strict_eligible_decoys,
)
from build_register_aware_benchmark import (  # noqa: E402
    build_decoy_benchmark,
    build_pair_universe,
    render_benchmark_readme,
    resolve_candidate_register,
)


class RegisterAwareBenchmarkTests(unittest.TestCase):
    @staticmethod
    def assessable_row(pair_id: str, validation: str, eligible: bool) -> dict[str, object]:
        return {
            "pair_id": pair_id,
            "pair_validation": validation,
            "register_assessment": "assessable_register_hypothesis",
            "decoy_background_eligible": eligible,
            "ebv_peptide": "ABCDEFGHIJK",
            "human_peptide": "LMNOPQRSTUV",
            "ebv_plddt": 80.0,
            "human_plddt": 80.0,
            "ebv_binding_rank": 1.0,
            "human_binding_rank": 1.0,
        }

    def test_assessable_pair_requires_same_register_aligned_residue(self):
        alignment = [(4, "Y", 3, "H"), (7, "F", 6, "F")]
        self.assertTrue(is_assessable_same_register_pair(alignment, 4, 3))
        self.assertFalse(is_assessable_same_register_pair(alignment, 4, 5))

    def test_strict_decoys_reject_binding_bin_mismatches(self):
        target = {
            "pair_id": "target",
            "ebv_peptide": "ABCDEFGHIJK",
            "human_peptide": "LMNOPQRSTUV",
            "ebv_plddt": 80,
            "human_plddt": 80,
            "ebv_binding_rank": 1,
            "human_binding_rank": 1,
        }
        mismatched = {**target, "pair_id": "bad", "ebv_binding_rank": 11}

        selected, available = strict_eligible_decoys(target, [mismatched])

        self.assertEqual(selected, [])
        self.assertEqual(available, 0)

    def test_pair_universe_excludes_zero_same_register_background(self):
        rows = build_pair_universe(
            geometry_rows=[
                {
                    "status": "PASS",
                    "ebv_candidate_id": "E",
                    "human_candidate_id": "H",
                    "aligned_positions_ebv_to_human": "4Y:3H",
                    "ebv_peptide_mean_plddt": "80",
                    "human_peptide_mean_plddt": "81",
                }
            ],
            prediction_by_candidate={
                "E": {
                    "predicted_core_start_positions_1_based": "4",
                    "predicted_core_fully_contained_in_manifest_peptide": "True",
                    "predicted_percentile_rank": "1",
                },
                "H": {
                    "predicted_core_start_positions_1_based": "5",
                    "predicted_core_fully_contained_in_manifest_peptide": "True",
                    "predicted_percentile_rank": "1",
                },
            },
            manifest_by_candidate={
                "E": {"peptide": "ABCDEFGHIJK"},
                "H": {"peptide": "LMNOPQRSTUV"},
            },
        )

        self.assertEqual(rows[0]["register_assessment"], "no_same_register_local_alignment")
        self.assertFalse(rows[0]["decoy_background_eligible"])

    def test_pair_universe_excludes_nonprimary_calibration_override(self):
        geometry = [{
            "status": "PASS",
            "ebv_candidate_id": "E",
            "human_candidate_id": "H",
            "aligned_positions_ebv_to_human": "5Y:5V",
            "ebv_peptide_mean_plddt": "80",
            "human_peptide_mean_plddt": "81",
        }]
        predictions = {
            "E": {
                "predicted_core_start_positions_1_based": "4",
                "predicted_core_peptide": "VYHFVKKHV",
                "predicted_core_fully_contained_in_manifest_peptide": "True",
                "predicted_percentile_rank": "1",
            },
            "H": {
                "predicted_core_start_positions_1_based": "5",
                "predicted_core_peptide": "VHFFKNIVT",
                "predicted_core_fully_contained_in_manifest_peptide": "True",
                "predicted_percentile_rank": "1",
            },
        }
        manifest = {
            "E": {"peptide": "TGGVYHFVKKHVHES"},
            "H": {"peptide": "ENPVVHFFKNIVTPR"},
        }
        overrides = {
            "E": {
                "candidate_id": "E",
                "analysis_role": "calibration_only_nonprimary_allele",
                "presenting_allele": "HLA-DRB5*01:01",
                "core_peptide": "YHFVKKHVH",
                "core_start_1_based": "5",
                "register_source": "PDB 1H15",
            }
        }

        rows = build_pair_universe(
            geometry, predictions, manifest, overrides_by_candidate=overrides
        )

        self.assertEqual(rows[0]["register_assessment"], "calibration_only_nonprimary_allele")
        self.assertFalse(rows[0]["decoy_background_eligible"])
        self.assertEqual(rows[0]["ebv_register_source"], "PDB 1H15")

    def test_decoy_selection_uses_only_assessable_background_rows(self):
        target = self.assessable_row("target", "classic_component_only", False)
        valid = self.assessable_row("valid", "background", True)
        invalid = self.assessable_row("invalid", "background", False)

        decoys, feasibility = build_decoy_benchmark(
            [target, valid, invalid], target_decoy_count=1
        )

        self.assertEqual([row["decoy_pair_id"] for row in decoys], ["valid"])
        self.assertEqual(feasibility[0]["eligible_decoy_count"], 1)

    def test_unassessable_target_is_reported_without_decoy_substitution(self):
        target = {
            **self.assessable_row("target", "classic_component_only", False),
            "register_assessment": "no_same_register_local_alignment",
        }
        decoys, feasibility = build_decoy_benchmark([target], target_decoy_count=5)

        self.assertEqual(decoys, [])
        self.assertEqual(feasibility[0]["readiness_status"], "not assessable")

    def test_report_does_not_overstate_tcr_or_disease_evidence(self):
        report = render_benchmark_readme(
            {"ready_targets": 0, "partial_targets": 1, "not_assessable_targets": 2}
        )

        self.assertIn("does not establish shared-TCR binding", report)
        self.assertIn("MS mechanism", report)
        self.assertIn("ready_targets: 0", report)
        self.assertIn("PDB 1BX2", report)
        self.assertIn("sensitivity-only", report)

    def test_primary_experimental_override_has_priority_over_predictor(self):
        prediction = {
            "predicted_core_start_positions_1_based": "4",
            "predicted_core_peptide": "NPVVHFFKN",
            "predicted_core_fully_contained_in_manifest_peptide": "True",
        }
        override = {
            "candidate_id": "HUMAN_MYELIN_13572",
            "analysis_role": "primary_experimental_reference",
            "presenting_allele": "HLA-DRB1*15:01",
            "core_peptide": "VHFFKNIVT",
            "core_start_1_based": "5",
            "register_source": "PDB 1BX2",
        }

        resolved = resolve_candidate_register(
            "HUMAN_MYELIN_13572", "ENPVVHFFKNIVTPR", prediction, override
        )

        self.assertEqual(resolved["register_status"], "experimental_primary_allele_reference")
        self.assertEqual(resolved["core_start_1_based"], 5)
        self.assertEqual(resolved["core_peptide"], "VHFFKNIVT")

    def test_nonprimary_and_sensitivity_overrides_cannot_be_benchmarked(self):
        prediction = {
            "predicted_core_start_positions_1_based": "4",
            "predicted_core_peptide": "VYHFVKKHV",
            "predicted_core_fully_contained_in_manifest_peptide": "True",
        }
        calibration = {
            "candidate_id": "EBV_TCELL_63843",
            "analysis_role": "calibration_only_nonprimary_allele",
            "presenting_allele": "HLA-DRB5*01:01",
            "core_peptide": "YHFVKKHVH",
            "core_start_1_based": "5",
            "register_source": "PDB 1H15",
        }
        sensitivity = {
            "candidate_id": "EBV_TCELL_2268683",
            "analysis_role": "sensitivity_only_unresolved",
            "presenting_allele": "HLA-DRB1*15:01",
            "core_peptide": "",
            "core_start_1_based": "",
            "register_source": "No exact experimental register located",
        }

        calibration_resolved = resolve_candidate_register(
            "EBV_TCELL_63843", "TGGVYHFVKKHVHES", prediction, calibration
        )
        sensitivity_resolved = resolve_candidate_register(
            "EBV_TCELL_2268683", "EKQLFYYIGTMLPN", prediction, sensitivity
        )

        self.assertEqual(
            calibration_resolved["register_status"], "calibration_only_nonprimary_allele"
        )
        self.assertEqual(
            sensitivity_resolved["register_status"], "sensitivity_only_unresolved"
        )


if __name__ == "__main__":
    unittest.main()
