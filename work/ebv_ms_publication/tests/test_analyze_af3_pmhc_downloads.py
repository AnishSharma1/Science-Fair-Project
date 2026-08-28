"""Tests for descriptive, source-traceable AF3 pMHC extraction."""

import json
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analyze_af3_pmhc_downloads import analyze_complete_job, request_details, study_candidate_id, summarize_cohorts, background_descriptive_rows  # noqa: E402


class AnalyzeAf3PmhcDownloadsTests(unittest.TestCase):
    def test_keeps_binding_prediction_and_structure_fields_descriptive_in_background_table(self):
        jobs = [{"candidate_id": "HUMAN_BACKGROUND_101", "af3_cohort": "new_human_background_pmhc", "selected_model_iptm": "0.9", "selected_model_peptide_mean_plddt": "50.0", "selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_median_A": "1.2"}]
        metadata = {"HUMAN_BACKGROUND_101": {"peptide": "ACDEFGHIKLM", "source_antigen": "Comparator", "binding_rank_bin": "weak", "predicted_core_peptide": "DEFGHIKLM"}}

        rows = background_descriptive_rows(jobs, metadata)

        self.assertEqual(rows, [{
            "candidate_id": "HUMAN_BACKGROUND_101", "peptide": "ACDEFGHIKLM", "source_antigen": "Comparator",
            "iedb_predicted_binding_rank_bin": "weak", "iedb_predicted_core": "DEFGHIKLM",
            "selected_model_iptm": "0.9", "selected_model_peptide_mean_plddt": "50.0",
            "within_job_peptide_pose_rmsd_median_A": "1.2",
            "interpretation": "Descriptive pMHC context only; not a decoy decision or biological ranking.",
        }])

    def test_summarizes_cohorts_without_using_scores_to_rank_candidates(self):
        rows = [
            {"af3_cohort": "background", "selected_model_iptm": "0.90", "selected_model_peptide_mean_plddt": "40.0", "selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_median_A": "1.0"},
            {"af3_cohort": "background", "selected_model_iptm": "0.94", "selected_model_peptide_mean_plddt": "60.0", "selected_to_other_sample_peptide_ca_rmsd_after_hla_groove_fit_median_A": "3.0"},
        ]

        summary = summarize_cohorts(rows)

        self.assertEqual(summary, [{
            "af3_cohort": "background", "completed_jobs": 2,
            "selected_model_iptm_median": 0.92,
            "selected_model_peptide_mean_plddt_median": 50.0,
            "within_job_peptide_pose_rmsd_median_A": 2.0,
            "interpretation": "Descriptive cohort summary only; no between-peptide biological ranking.",
        }])

    def test_identifies_an_unmapped_job_before_pmhc_chain_validation(self):
        request = [{"name": "decoy_02_hy_enga_drb1_s101", "sequences": [{"proteinChain": {"sequence": "A"}}]}]

        candidate_id = study_candidate_id(request, {"HUMAN_BACKGROUND_101": {}})

        self.assertIsNone(candidate_id)

    def test_uses_the_third_submitted_chain_as_the_expected_peptide(self):
        request = [{
            "name": "HUMAN_BACKGROUND_101",
            "modelSeeds": [7],
            "sequences": [
                {"proteinChain": {"sequence": "DRA"}},
                {"proteinChain": {"sequence": "DRB"}},
                {"proteinChain": {"sequence": "ACDEFGHIKLM"}},
            ],
        }]

        details = request_details(request)

        self.assertEqual(details["requested_peptide"], "ACDEFGHIKLM")
        self.assertEqual(details["server_seed"], 7)

    def test_selects_highest_server_ranking_score_and_rejects_wrong_peptide_sequence(self):
        with tempfile.TemporaryDirectory() as temporary:
            job = Path(temporary) / "HUMAN_BACKGROUND_101"
            job.mkdir()
            request = [{
                "name": "HUMAN_BACKGROUND_101",
                "modelSeeds": [7],
                "sequences": [
                    {"proteinChain": {"sequence": "DRA"}},
                    {"proteinChain": {"sequence": "DRB"}},
                    {"proteinChain": {"sequence": "ACDEFGHIKLM"}},
                ],
            }]
            (job / "fold_HUMAN_BACKGROUND_101_job_request.json").write_text(json.dumps(request))
            for index, score, peptide in [(0, 0.10, "ACDEFGHIKLM"), (1, 0.90, "ACDEFGHIKLM")]:
                (job / f"fold_HUMAN_BACKGROUND_101_summary_confidences_{index}.json").write_text(json.dumps({
                    "ranking_score": score, "iptm": 0.8, "ptm": 0.9,
                    "fraction_disordered": 0.0, "has_clash": 0.0,
                }))
                (job / f"fold_HUMAN_BACKGROUND_101_model_{index}.cif").write_text(_cif(peptide))

            result = analyze_complete_job(job, "json3folds", {"HUMAN_BACKGROUND_101": {}})

        self.assertEqual(result["job"]["selected_model_index"], 1)
        self.assertEqual(result["job"]["sequence_layout_status"], "pass_exact_three_chain_peptide_match")

        with tempfile.TemporaryDirectory() as temporary:
            job = Path(temporary) / "HUMAN_BACKGROUND_101"
            job.mkdir()
            (job / "fold_HUMAN_BACKGROUND_101_job_request.json").write_text(json.dumps(request))
            (job / "fold_HUMAN_BACKGROUND_101_summary_confidences_0.json").write_text(json.dumps({
                "ranking_score": 0.9, "iptm": 0.8, "ptm": 0.9,
                "fraction_disordered": 0.0, "has_clash": 0.0,
            }))
            (job / "fold_HUMAN_BACKGROUND_101_model_0.cif").write_text(_cif("ACDEFGHIKLN"))

            result = analyze_complete_job(job, "json3folds", {"HUMAN_BACKGROUND_101": {}})

        self.assertEqual(result["job"]["sequence_layout_status"], "fail_peptide_sequence_mismatch")


def _cif(peptide: str) -> str:
    aa3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN"}
    rows = []
    atom_id = 1
    for chain, sequence in [("A", "AAAA"), ("B", "AAAA"), ("C", peptide)]:
        for position, residue in enumerate(sequence, start=1):
            rows.append(f"ATOM {atom_id} C CA . {aa3[residue]} {chain} 1 {position} ? {atom_id}.0 0.0 0.0 1.00 80.0 {position} {chain} 1")
            atom_id += 1
    return "\n".join([
        "data_test", "#", "loop_", "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
        "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
        "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x",
        "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id", "_atom_site.auth_asym_id", "_atom_site.pdbx_PDB_model_num", *rows, "#",
    ])


if __name__ == "__main__":
    unittest.main()
