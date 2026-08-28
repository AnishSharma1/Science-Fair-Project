"""Tests for the predeclared structural-background control expansion."""

import json
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_structural_control_expansion import (  # noqa: E402
    build_alphafold_jobs,
    extend_score_sheet_with_layered_controls,
    freeze_control_registry,
    inventory_expected_control_jobs,
    prediction_rows_from_nextgen_result,
    select_layered_controls,
    summarize_layer_geometry,
    uncovered_pairs_by_stratum,
)


def raw_record(
    epitope_id: int,
    peptide: str,
    *,
    accession: str = "P12345.1",
    antigen: str = "Comparator protein",
    start: int = 10,
) -> dict[str, object]:
    return {
        "tcell_id": epitope_id + 1000,
        "structure_id": epitope_id,
        "linear_sequence": peptide,
        "mhc_restriction": "HLA-DRB1*15:01",
        "curated_source_antigen": {
            "accession": accession,
            "name": antigen,
            "starting_position": start,
            "ending_position": start + len(peptide) - 1,
            "source_organism_name": "Homo sapiens (human)",
        },
    }


class StructuralControlExpansionTests(unittest.TestCase):
    def test_registry_deduplicates_peptides_and_audits_exclusions(self):
        eligible = raw_record(101, "ACDEFGHIKLMNPQRSTVWY")
        duplicate = raw_record(102, "ACDEFGHIKLMNPQRSTVWY")
        mbp = raw_record(
            103,
            "CDEFGHIKLMNPQRSTVWYA",
            accession="P02686.3",
            antigen="Myelin basic protein",
        )
        study = raw_record(104, "DEFGHIKLMNPQRSTVWYAC")
        bad_coordinates = raw_record(105, "EFGHIKLMNPQRSTVWYACD")
        bad_coordinates["curated_source_antigen"]["ending_position"] += 1

        rows = freeze_control_registry(
            [duplicate, mbp, study, bad_coordinates, eligible],
            study_peptides={"DEFGHIKLMNPQRSTVWYAC"},
        )

        self.assertEqual(len(rows), 4)
        by_id = {row["iedb_epitope_id"]: row for row in rows}
        self.assertEqual(by_id[101]["selection_status"], "eligible_pre_prediction")
        self.assertEqual(by_id[101]["source_record_count"], 2)
        self.assertEqual(by_id[103]["selection_status"], "excluded_mbp_plp_mog")
        self.assertEqual(by_id[104]["selection_status"], "excluded_study_candidate")
        self.assertEqual(by_id[105]["selection_status"], "excluded_invalid_coordinates")

    def test_selection_ignores_discovery_and_structure_fields_and_separates_layers(self):
        predictions = [
            {
                "candidate_id": "HUMAN_BACKGROUND_20",
                "iedb_epitope_id": "20",
                "peptide": "A" * 21,
                "peptide_length": "21",
                "binding_rank_bin": "weak",
                "predicted_percentile_rank": "40.0",
                "discovery_priority_rank": "1",
                "candidate_exposed_ca_rmsd_A": "0.1",
            },
            {
                "candidate_id": "HUMAN_BACKGROUND_10",
                "iedb_epitope_id": "10",
                "peptide": "C" * 21,
                "peptide_length": "21",
                "binding_rank_bin": "weak",
                "predicted_percentile_rank": "60.0",
                "discovery_priority_rank": "999",
                "candidate_exposed_ca_rmsd_A": "99.0",
            },
            {
                "candidate_id": "HUMAN_BACKGROUND_30",
                "iedb_epitope_id": "30",
                "peptide": "D" * 22,
                "peptide_length": "22",
                "binding_rank_bin": "intermediate",
                "predicted_percentile_rank": "8.0",
            },
            {
                "candidate_id": "HUMAN_BACKGROUND_40",
                "iedb_epitope_id": "40",
                "peptide": "E" * 25,
                "peptide_length": "25",
                "binding_rank_bin": "weak",
                "predicted_percentile_rank": "55.0",
            },
            {
                "candidate_id": "HUMAN_BACKGROUND_50",
                "iedb_epitope_id": "50",
                "peptide": "F" * 29,
                "peptide_length": "29",
                "binding_rank_bin": "intermediate",
                "predicted_percentile_rank": "7.0",
            },
        ]
        served = {21: ["PAIR_21"], 23: [], 25: ["PAIR_25"], 32: ["PAIR_32_A", "PAIR_32_B"]}

        selected, feasibility = select_layered_controls(predictions, served, limit=3)

        rows_21 = [row for row in selected if row["stratum_length"] == 21]
        self.assertEqual(
            [(row["candidate_id"], row["analysis_layer"]) for row in rows_21],
            [
                ("HUMAN_BACKGROUND_10", "primary_exact_bin_length_pm1"),
                ("HUMAN_BACKGROUND_20", "primary_exact_bin_length_pm1"),
                ("HUMAN_BACKGROUND_30", "binding_bin_sensitivity_length_pm1"),
            ],
        )
        rows_32 = [row for row in selected if row["stratum_length"] == 32]
        self.assertEqual(
            [(row["candidate_id"], row["analysis_layer"]) for row in rows_32],
            [
                ("HUMAN_BACKGROUND_40", "length_sensitivity_exact_bin_pm7"),
                ("HUMAN_BACKGROUND_50", "length_plus_binding_sensitivity_pm7"),
            ],
        )
        status_32 = next(row for row in feasibility if row["stratum_length"] == 32)
        self.assertEqual(status_32["primary_assessment"], "not_assessable_no_31_to_33aa_direct_controls")
        self.assertEqual(status_32["control_shortfall"], 1)
        self.assertNotIn("discovery_priority_rank", selected[0])
        self.assertNotIn("candidate_exposed_ca_rmsd_A", selected[0])

    def test_alphafold_jobs_deduplicate_reused_controls_and_keep_manifest_mappings(self):
        selected = [
            {
                "candidate_id": "HUMAN_BACKGROUND_10",
                "iedb_epitope_id": "10",
                "peptide": "A" * 21,
                "peptide_length": 21,
                "stratum_length": 21,
                "analysis_layer": "primary_exact_bin_length_pm1",
                "binding_bin": "weak",
                "selection_order": 1,
                "served_pair_ids": "PAIR_21",
            },
            {
                "candidate_id": "HUMAN_BACKGROUND_10",
                "iedb_epitope_id": "10",
                "peptide": "A" * 21,
                "peptide_length": 21,
                "stratum_length": 32,
                "analysis_layer": "length_sensitivity_exact_bin_pm7",
                "binding_bin": "weak",
                "selection_order": 1,
                "served_pair_ids": "PAIR_32",
            },
        ]

        jobs, manifest = build_alphafold_jobs(selected, "DRASEQ", "DRBSEQ")

        self.assertEqual(len(jobs), 1)
        self.assertEqual([chain["proteinChain"]["sequence"] for chain in jobs[0]["sequences"]], ["DRASEQ", "DRBSEQ", "A" * 21])
        self.assertEqual(len(manifest), 2)
        self.assertEqual({row["stratum_length"] for row in manifest}, {21, 32})
        self.assertTrue(all(row["job_name"] == jobs[0]["name"] for row in manifest))

    def test_layer_summary_equal_weights_unique_controls_despite_repeat_rows(self):
        geometry = [
            {"background_candidate_id": "D1", "candidate_exposed_ca_rmsd_A": 2.0},
            {"background_candidate_id": "D1", "candidate_exposed_ca_rmsd_A": 4.0},
            {"background_candidate_id": "D1", "candidate_exposed_ca_rmsd_A": 100.0},
            {"background_candidate_id": "D2", "candidate_exposed_ca_rmsd_A": 8.0},
        ]

        summary = summarize_layer_geometry(geometry, target_median=1.0)

        self.assertEqual(summary["unique_control_count"], 2)
        self.assertEqual(summary["technical_geometry_count"], 4)
        self.assertEqual(summary["background_control_median_A"], 6.0)
        self.assertEqual(summary["background_minus_target_median_A"], 5.0)
        self.assertEqual(summary["p_value"], "")

    def test_extended_score_sheet_keeps_ranking_and_layers_separate_when_pending(self):
        baseline = [{
            "discovery_priority_rank": "7",
            "pair_id": "PAIR_21",
            "human_peptide": "A" * 21,
            "target_candidate_exposed_rmsd_median_A": "1.25",
            "claim_boundary": "Existing claim boundary.",
        }]
        selected = [
            {
                "candidate_id": "D_PRIMARY",
                "analysis_layer": "primary_exact_bin_length_pm1",
                "served_pair_ids": "PAIR_21",
            },
            {
                "candidate_id": "D_SENSITIVITY",
                "analysis_layer": "binding_bin_sensitivity_length_pm1",
                "served_pair_ids": "PAIR_21",
            },
        ]

        rows = extend_score_sheet_with_layered_controls(baseline, selected, [])

        self.assertEqual(rows[0]["discovery_priority_rank"], "7")
        self.assertEqual(rows[0]["expanded_primary_selected_control_ids"], "D_PRIMARY")
        self.assertEqual(rows[0]["expanded_binding_sensitivity_selected_control_ids"], "D_SENSITIVITY")
        self.assertEqual(rows[0]["expanded_primary_completed_control_count"], 0)
        self.assertEqual(rows[0]["expanded_primary_background_control_median_A"], "")
        self.assertEqual(rows[0]["expanded_control_status"], "awaiting_alphafold_downloads")

    def test_extended_score_sheet_uses_existing_pair_geometry_for_new_control_delta(self):
        baseline = [{
            "discovery_priority_rank": "7",
            "pair_id": "PAIR_21",
            "human_peptide": "A" * 21,
            "candidate_exposed_ca_rmsd_A_median": "1.25",
            "target_candidate_exposed_rmsd_median_A": "",
        }]
        selected = [{
            "candidate_id": "D_PRIMARY",
            "analysis_layer": "primary_exact_bin_length_pm1",
            "served_pair_ids": "PAIR_21",
        }]
        geometry = [{
            "pair_id": "PAIR_21",
            "analysis_layer": "primary_exact_bin_length_pm1",
            "background_candidate_id": "D_PRIMARY",
            "candidate_exposed_ca_rmsd_A": 4.25,
        }]

        rows = extend_score_sheet_with_layered_controls(baseline, selected, geometry)

        self.assertEqual(rows[0]["expanded_primary_background_control_median_A"], 4.25)
        self.assertEqual(rows[0]["expanded_primary_background_minus_target_median_A"], 3.0)

    def test_only_uncovered_primary_allele_pairs_define_control_strata(self):
        rows = [
            {
                "pair_id": "PAIR_21",
                "human_peptide": "A" * 21,
                "register_eligible_primary_allele": "True",
                "structural_background_comparator_count": "0",
            },
            {
                "pair_id": "ALREADY_DONE",
                "human_peptide": "C" * 21,
                "register_eligible_primary_allele": "True",
                "structural_background_comparator_count": "2",
            },
            {
                "pair_id": "INELIGIBLE",
                "human_peptide": "D" * 23,
                "register_eligible_primary_allele": "False",
                "structural_background_comparator_count": "0",
            },
        ]

        strata = uncovered_pairs_by_stratum(rows)

        self.assertEqual(strata, {21: ["PAIR_21"], 23: [], 25: [], 32: []})

    def test_download_inventory_deduplicates_manifest_and_requires_exact_request(self):
        manifest = [
            {
                "batch_file": "jobs.json",
                "candidate_id": "HUMAN_BACKGROUND_10",
                "job_name": "ebvms_bg_HUMAN_BACKGROUND_10_s04",
                "peptide": "A" * 21,
                "stratum_length": 21,
            },
            {
                "batch_file": "jobs.json",
                "candidate_id": "HUMAN_BACKGROUND_10",
                "job_name": "ebvms_bg_HUMAN_BACKGROUND_10_s04",
                "peptide": "A" * 21,
                "stratum_length": 32,
            },
        ]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            job = root / "ebvms_bg_human_background_10_s04"
            job.mkdir()
            prefix = "fold_ebvms_bg_human_background_10_s04"
            request = [{
                "name": "ebvms_bg_HUMAN_BACKGROUND_10_s04",
                "modelSeeds": [104773],
                "sequences": [
                    {"proteinChain": {"sequence": sequence, "count": 1}}
                    for sequence in ("DRA", "DRB", "A" * 21)
                ],
            }]
            (job / f"{prefix}_job_request.json").write_text(json.dumps(request), encoding="utf-8")
            for index in range(5):
                for suffix in (
                    f"model_{index}.cif",
                    f"summary_confidences_{index}.json",
                    f"full_data_{index}.json",
                ):
                    (job / f"{prefix}_{suffix}").touch()

            inventory = inventory_expected_control_jobs(root, manifest)

        self.assertEqual(len(inventory), 1)
        self.assertEqual(inventory[0]["completeness_status"], "complete_five_sample_exact_sequence")
        self.assertEqual(inventory[0]["served_stratum_lengths"], "21;32")

    def test_nextgen_predictions_map_by_sequence_number_and_preserve_core(self):
        candidates = [
            {
                "candidate_id": "HUMAN_BACKGROUND_10",
                "arm": "Expanded human background comparator",
                "evidence_tier": "Tier 4",
                "peptide": "A" * 20,
                "peptide_length": "20",
                "hla": "HLA-DRB1*15:01",
                "iedb_epitope_id": "10",
                "source_accession": "P1",
                "source_antigen_name": "Protein 1",
                "source_start_1_based": "1",
                "source_end_1_based": "20",
            },
            {
                "candidate_id": "HUMAN_BACKGROUND_20",
                "arm": "Expanded human background comparator",
                "evidence_tier": "Tier 4",
                "peptide": "C" * 21,
                "peptide_length": "21",
                "hla": "HLA-DRB1*15:01",
                "iedb_epitope_id": "20",
                "source_accession": "P2",
                "source_antigen_name": "Protein 2",
                "source_start_1_based": "2",
                "source_end_1_based": "22",
            },
        ]
        columns = [
            "sequence_number", "peptide", "start", "end", "length", "allele",
            "peptide_index", "median_percentile", "netmhciipan_ba_core",
            "netmhciipan_ba_ic50", "netmhciipan_ba_percentile",
        ]
        result = {
            "status": "done",
            "data": {
                "errors": [],
                "results": [{
                    "type": "peptide_table",
                    "table_columns": [{"name": name} for name in columns],
                    "table_data": [
                        [2, "C" * 21, 1, 21, 21, "HLA-DRB1*15:01", 2, 20.0, "C" * 9, 900.0, 20.0],
                        [1, "A" * 20, 1, 20, 20, "HLA-DRB1*15:01", 1, 5.8, "A" * 9, 82.92, 5.8],
                    ],
                }],
            },
        }

        rows = prediction_rows_from_nextgen_result(
            candidates, result, "raw.json", "2026-08-15T00:00:00Z"
        )

        self.assertEqual([row["candidate_id"] for row in rows], ["HUMAN_BACKGROUND_10", "HUMAN_BACKGROUND_20"])
        self.assertEqual(rows[0]["predicted_core_peptide"], "A" * 9)
        self.assertEqual(rows[0]["binding_rank_bin"], "intermediate")
        self.assertEqual(rows[1]["binding_rank_bin"], "weak")
        self.assertEqual(rows[0]["iedb_epitope_id"], "10")


if __name__ == "__main__":
    unittest.main()
