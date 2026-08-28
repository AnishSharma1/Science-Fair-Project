import unittest

from premeeting_rigor import (
    binding_rank_bin,
    eligible_decoys,
    enumerate_core_windows,
    iedb_mhcii_eligible,
    iedb_submission_segments,
    natural_flank_submission_segment,
    map_prediction_rows,
    ordered_decoys,
    parse_iedb_mhcii_tsv,
)
from run_external_validation_benchmark import format_validation_summary_line


class PremeetingRigorTests(unittest.TestCase):
    def test_enumerate_core_windows_includes_every_nine_mer_position(self):
        """A one-residue off-by-one bug would silently omit a register hypothesis."""
        self.assertEqual(
            enumerate_core_windows("ABCDEFGHIJK"),
            [
                {"start": 1, "end": 9, "core_peptide": "ABCDEFGHI"},
                {"start": 2, "end": 10, "core_peptide": "BCDEFGHIJ"},
                {"start": 3, "end": 11, "core_peptide": "CDEFGHIJK"},
            ],
        )

    def test_parse_iedb_mhcii_tsv_preserves_predicted_core_and_rank(self):
        """A parser that drops core/rank fields would make register matching unauditable."""
        text = (
            "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tic50\trank\n"
            "HLA-DRB1*15:01\t1\t1\t15\t15\tVYHFVKKHV\tTGGVYHFVKKHVHES\t314.19\t20.0\n"
        )
        self.assertEqual(
            parse_iedb_mhcii_tsv(text),
            [{
                "allele": "HLA-DRB1*15:01",
                "seq_num": "1",
                "start": "1",
                "end": "15",
                "length": "15",
                "core_peptide": "VYHFVKKHV",
                "peptide": "TGGVYHFVKKHVHES",
                "ic50": "314.19",
                "rank": "20.0",
            }],
        )

    def test_ordered_decoys_ignore_priority_heuristic(self):
        """A decoy selector that reads the screen priority score would leak the outcome."""
        target = {
            "pair_id": "target",
            "ebv_peptide": "ACDEFGHIKL",
            "human_peptide": "MNPQRSTVWY",
            "ebv_plddt": 80.0,
            "human_plddt": 85.0,
            "ebv_binding_rank": 4.0,
            "human_binding_rank": 5.0,
        }
        candidates = [
            {
                "pair_id": "matching_covariates",
                "ebv_peptide": "ACDEFGHIKL",
                "human_peptide": "MNPQRSTVWY",
                "ebv_plddt": 80.0,
                "human_plddt": 85.0,
                "ebv_binding_rank": 4.0,
                "human_binding_rank": 5.0,
                "review_priority_heuristic": 0.001,
            },
            {
                "pair_id": "high_priority_bad_match",
                "ebv_peptide": "ACDEFGHIKLMN",
                "human_peptide": "MNPQRSTVWYAC",
                "ebv_plddt": 30.0,
                "human_plddt": 35.0,
                "ebv_binding_rank": 30.0,
                "human_binding_rank": 30.0,
                "review_priority_heuristic": 0.999,
            },
        ]
        result = ordered_decoys(target, candidates, limit=1)
        self.assertEqual(result[0]["pair_id"], "matching_covariates")
        self.assertNotIn("review_priority_heuristic", result[0])

    def test_binding_rank_bin_has_fixed_boundary_behavior(self):
        """A rank-bin boundary error would make matched-decoy strata inconsistent."""
        self.assertEqual(binding_rank_bin(2.0), "strong")
        self.assertEqual(binding_rank_bin(2.1), "intermediate")
        self.assertEqual(binding_rank_bin(10.0), "intermediate")
        self.assertEqual(binding_rank_bin(10.1), "weak")

    def test_iedb_mhcii_eligibility_respects_supported_length_range(self):
        """Submitting a 10-mer or 31-mer silently loses a candidate in the API response."""
        self.assertFalse(iedb_mhcii_eligible("A" * 10))
        self.assertTrue(iedb_mhcii_eligible("A" * 11))
        self.assertTrue(iedb_mhcii_eligible("A" * 30))
        self.assertFalse(iedb_mhcii_eligible("A" * 31))

    def test_long_peptide_segments_cover_every_possible_nine_mer_core(self):
        """A gap between 30-mer tiles would erase valid register hypotheses from long peptides."""
        segments = iedb_submission_segments({"candidate_id": "LONG", "peptide": "A" * 60})
        self.assertEqual(
            [(segment["source_start_1_based"], len(segment["peptide"])) for segment in segments],
            [(1, 30), (23, 30), (31, 30)],
        )
        covered_core_starts = set()
        for segment in segments:
            covered_core_starts.update(
                range(segment["source_start_1_based"], segment["source_start_1_based"] + 22)
            )
        self.assertEqual(covered_core_starts, set(range(1, 53)))

    def test_natural_flank_submission_preserves_the_original_short_peptide(self):
        """A shifted flank window would attach a predicted core to the wrong source epitope."""
        segment = natural_flank_submission_segment(
            {"candidate_id": "SHORT", "peptide": "VLRYHVLLEE"},
            {
                "extended_sequence": "EDTVVLRYHVLLEEIIER",
                "original_start_in_extended_1_based": "5",
                "original_end_in_extended_1_based": "14",
            },
        )
        self.assertEqual(segment["peptide"], "EDTVVLRYHVLLEEIIER")
        self.assertEqual(segment["original_start_in_submission_1_based"], 5)
        self.assertEqual(segment["submission_strategy"], "verified_natural_flank_extension")

    def test_eligible_decoys_excludes_candidates_outside_prespecified_tolerances(self):
        """A nearly matched row must not be mislabeled as an analysis-ready decoy."""
        target = {
            "pair_id": "target",
            "ebv_peptide": "ACDEFGHIKL",
            "human_peptide": "MNPQRSTVWY",
            "ebv_plddt": 80.0,
            "human_plddt": 85.0,
            "ebv_binding_rank": 4.0,
            "human_binding_rank": 5.0,
        }
        candidates = [
            {
                "pair_id": "eligible",
                "ebv_peptide": "ACDEFGHIKL",
                "human_peptide": "MNPQRSTVWY",
                "ebv_plddt": 80.0,
                "human_plddt": 85.0,
                "ebv_binding_rank": 4.0,
                "human_binding_rank": 5.0,
            },
            {
                "pair_id": "wrong_binding_bin",
                "ebv_peptide": "ACDEFGHIKL",
                "human_peptide": "MNPQRSTVWY",
                "ebv_plddt": 80.0,
                "human_plddt": 85.0,
                "ebv_binding_rank": 50.0,
                "human_binding_rank": 50.0,
            },
        ]
        selected, available = eligible_decoys(target, candidates, limit=5)
        self.assertEqual(available, 1)
        self.assertEqual([row["pair_id"] for row in selected], ["eligible"])

    def test_map_prediction_rows_rejects_a_sequence_mismatch(self):
        """A reordered API response must not silently attach a core to the wrong peptide."""
        candidates = [
            {"candidate_id": "EBV_1", "peptide": "TGGVYHFVKKHVHES"},
            {"candidate_id": "HUMAN_1", "peptide": "ENPVVHFFKNIVTPR"},
        ]
        response_rows = [
            {"seq_num": "1", "peptide": "TGGVYHFVKKHVHES", "core_peptide": "VYHFVKKHV", "rank": "20.0"},
            {"seq_num": "2", "peptide": "WRONGSEQUENCE", "core_peptide": "WRONGSEQU", "rank": "40.0"},
        ]
        with self.assertRaisesRegex(ValueError, "does not match candidate peptide"):
            map_prediction_rows(candidates, response_rows)

    def test_map_prediction_rows_uses_seq_num_not_response_order(self):
        """IEDB may return prediction rows out of FASTA order, so zip mapping would mislabel cores."""
        candidates = [
            {"candidate_id": "EBV_1", "peptide": "TGGVYHFVKKHVHES"},
            {"candidate_id": "HUMAN_1", "peptide": "ENPVVHFFKNIVTPR"},
        ]
        response_rows = [
            {"seq_num": "2", "peptide": "ENPVVHFFKNIVTPR", "core_peptide": "VHFFKNIVT", "rank": "0.08"},
            {"seq_num": "1", "peptide": "TGGVYHFVKKHVHES", "core_peptide": "VYHFVKKHV", "rank": "20.0"},
        ]
        mapped = map_prediction_rows(candidates, response_rows)
        self.assertEqual(mapped[0]["candidate_id"], "EBV_1")
        self.assertEqual(mapped[0]["core_peptide"], "VYHFVKKHV")
        self.assertEqual(mapped[1]["candidate_id"], "HUMAN_1")
        self.assertEqual(mapped[1]["core_peptide"], "VHFFKNIVT")

    def test_validation_summary_describes_records_not_independent_positives(self):
        """Overlap in the BALF5--MBP family makes 'positive pairs' scientifically misleading."""
        line = format_validation_summary_line({
            "test": "classic_BALF5_MBP_pair_recovery",
            "positive_pair_count": 10,
            "scored_pair_universe_n": 32,
            "observed_mean_rank": 7.1,
            "observed_top10_positive_count": 8,
            "empirical_p_top10_count_ge_observed": 0.0001,
        })
        self.assertIn("annotated pair records", line)
        self.assertNotIn("positive pairs", line)


if __name__ == "__main__":
    unittest.main()
