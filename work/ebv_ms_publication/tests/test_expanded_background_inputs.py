"""Tests for source-traceable human-background pMHC batch preparation."""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_expanded_background_inputs import (  # noqa: E402
    build_background_candidates,
    build_colabfold_batch_rows,
    potential_target_coverage,
    render_readme,
)


class ExpandedBackgroundInputTests(unittest.TestCase):
    def test_selects_only_coordinate_validated_direct_mhcii_backgrounds(self):
        source_rows = [
            {
                "iedb_epitope_id": "101",
                "peptide": "ACDEFGHIKLM",
                "candidate_class": "human_background",
                "provenance_status": "coordinate_validated",
                "source_antigen_name": "Comparator protein",
                "mhc_allele": "HLA-DRB1*15:01",
            },
            {
                "iedb_epitope_id": "102",
                "peptide": "ACDEFGHIK",
                "candidate_class": "human_background",
                "provenance_status": "coordinate_validated",
                "source_antigen_name": "Too short without a verified flank",
                "mhc_allele": "HLA-DRB1*15:01",
            },
            {
                "iedb_epitope_id": "103",
                "peptide": "ACDEFGHIKLM",
                "candidate_class": "myelin_candidate",
                "provenance_status": "coordinate_validated",
                "source_antigen_name": "Not a comparator",
                "mhc_allele": "HLA-DRB1*15:01",
            },
        ]

        registry, candidates = build_background_candidates(source_rows)

        self.assertEqual(len(registry), 3)
        self.assertEqual([candidate["candidate_id"] for candidate in candidates], ["HUMAN_BACKGROUND_101"])
        self.assertEqual(registry[1]["selection_status"], "retained_not_modeled_missing_verified_flank")
        self.assertEqual(registry[2]["selection_status"], "excluded_not_human_background")

    def test_uses_the_existing_hla_chains_and_preserves_full_peptides(self):
        candidates = [{
            "candidate_id": "HUMAN_BACKGROUND_101",
            "peptide": "ACDEFGHIKLM",
        }]

        batch = build_colabfold_batch_rows(candidates, "DRASEQ", "DRBSEQ")

        self.assertEqual(batch, [{
            "id": "HUMAN_BACKGROUND_101",
            "sequence": "DRASEQ:DRBSEQ:ACDEFGHIKLM",
        }])

    def test_reports_only_length_feasibility_before_binding_and_geometry_are_known(self):
        targets = [{
            "target_pair_id": "EBV::MYELIN",
            "target_register_assessment": "assessable_register_hypothesis",
            "readiness_status": "partial",
        }]
        universe = [{
            "pair_id": "EBV::MYELIN",
            "human_peptide": "ABCDEFGHIJKLMNOP",
            "ebv_peptide": "ABCDEFGHIJKLMN",
        }]
        candidates = [
            {"candidate_id": "HUMAN_BACKGROUND_A", "peptide": "ABCDEFGHIJKLMNO"},
            {"candidate_id": "HUMAN_BACKGROUND_B", "peptide": "ABCDEFGHIJKL"},
        ]

        rows = potential_target_coverage(targets, universe, candidates)

        self.assertEqual(rows[0]["direct_background_length_match_count"], 1)
        self.assertEqual(rows[0]["coverage_interpretation"], "length_feasible_pending_iedb_and_structure")

    def test_readme_counts_background_records_not_all_source_rows(self):
        registry = [
            {"candidate_class": "myelin_candidate", "selection_status": "excluded_not_human_background"},
            {"candidate_class": "human_background", "selection_status": "selected_for_direct_iedb_and_pmhc_batch"},
            {"candidate_class": "human_background", "selection_status": "retained_not_modeled_missing_verified_flank"},
        ]

        readme = render_readme(registry, [{"candidate_id": "HUMAN_BACKGROUND_1"}])

        self.assertIn("Human-background records reviewed: **2**", readme)
        self.assertIn("Source-table records scanned: **3**", readme)


if __name__ == "__main__":
    unittest.main()
