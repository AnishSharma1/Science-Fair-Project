import unittest


from build_combined_drb1501_same_register_rankings import (
    audit_legacy_eligibility,
    combine_and_rank,
)


def legacy_row(pair_id, left="LRALLARSH", right="LEARLSRMH"):
    return {
        "pair_id": pair_id,
        "ebv_candidate_id": f"OLD_E_{pair_id}",
        "human_candidate_id": f"OLD_S_{pair_id}",
        "ebv_peptide": f"AAA{left}AAA",
        "human_peptide": f"AAA{right}AAA",
        "ebv_top_core_peptide": left,
        "human_top_core_peptide": right,
        "ebv_register_status": "iedb_top_core_hypothesis",
        "human_register_status": "iedb_top_core_hypothesis",
        "pair_validation": "background",
        "score_coverage_status": "excluded_nonprimary_register_status",
        "ebv_binding_rank": "1.0",
        "human_binding_rank": "2.0",
    }


def v2_row(pair_id, left="LRALLARSH", right="LEARLSRMH"):
    return {
        "pair_id": pair_id,
        "ebv_candidate_id": f"NEW_E_{pair_id}",
        "self_candidate_id": f"NEW_S_{pair_id}",
        "ebv_protein": "EBNA1",
        "self_protein": "ANO2",
        "ebv_sequence": f"AAA{left}AAA",
        "self_sequence": f"AAA{right}AAA",
        "ebv_predicted_core": left,
        "self_predicted_core": right,
        "ebv_binding_percentile_rank": "1.0",
        "self_binding_percentile_rank": "2.0",
        "ebv_source_certainty": "exact",
        "self_source_certainty": "candidate",
        "hla_rank": "1",
    }


def annotation(candidate_id, protein, start, end, sequence):
    return {
        "candidate_id": candidate_id,
        "short_protein_name": protein,
        "parent_residue_start_1_based": str(start),
        "parent_residue_end_1_based": str(end),
        "resolved_parent_accession": "ACC",
        "peptide": sequence,
    }


class LegacyEligibilityTests(unittest.TestCase):
    def test_only_resolved_primary_drb1501_rows_are_admitted(self):
        valid = legacy_row("valid")
        nonprimary = legacy_row("nonprimary")
        nonprimary["ebv_register_status"] = "calibration_only_nonprimary_allele"
        unresolved_ebv = legacy_row("unresolved_ebv")
        unresolved_ebv["ebv_register_status"] = "sensitivity_only_unresolved"
        unresolved_ebv["ebv_top_core_peptide"] = ""
        unresolved_self = legacy_row("unresolved_self")
        unresolved_self["human_register_status"] = "unresolved_or_flank_dependent_core"
        eligible, audit = audit_legacy_eligibility(
            [valid, nonprimary, unresolved_ebv, unresolved_self]
        )
        self.assertEqual([row["pair_id"] for row in eligible], ["valid"])
        status = {row["pair_id"]: row["combined_eligibility_status"] for row in audit}
        self.assertEqual(status["valid"], "eligible_primary_drb1501_resolved_registers")
        self.assertEqual(status["nonprimary"], "excluded_nonprimary_hla")
        self.assertEqual(status["unresolved_ebv"], "excluded_unresolved_ebv_register")
        self.assertEqual(status["unresolved_self"], "excluded_unresolved_self_register")


class CombinationTests(unittest.TestCase):
    def test_exact_sequence_duplicate_is_one_pair_and_registers_must_agree(self):
        new = v2_row("new")
        old = legacy_row("old")
        annotations = [
            annotation(old["ebv_candidate_id"], "EBNA1", 1, 15, old["ebv_peptide"]),
            annotation(old["human_candidate_id"], "ANO2", 20, 34, old["human_peptide"]),
        ]
        combined, overlaps = combine_and_rank([new], [old], legacy_annotations=annotations)
        self.assertEqual(len(combined), 1)
        self.assertEqual(len(overlaps), 1)
        self.assertEqual(combined[0]["source_membership"], "v2_and_legacy_exact_duplicate")
        self.assertEqual(combined[0]["legacy_pair_id"], "old")

        old["human_top_core_peptide"] = "AAAAAAAAA"
        with self.assertRaisesRegex(ValueError, "register disagreement"):
            combine_and_rank([new], [old], legacy_annotations=annotations)

    def test_all_unique_pairs_receive_one_common_blosum_ranking(self):
        new = v2_row("new_best", "AAAAAAAAA", "AAAAAAAAA")
        old = legacy_row("old_worse", "RRRRRRRRR", "AAAAAAAAA")
        annotations = [
            annotation(old["ebv_candidate_id"], "EBV", 1, 15, old["ebv_peptide"]),
            annotation(old["human_candidate_id"], "SELF", 20, 34, old["human_peptide"]),
        ]
        combined, overlaps = combine_and_rank([new], [old], legacy_annotations=annotations)
        self.assertEqual(overlaps, [])
        self.assertEqual([row["combined_rank"] for row in combined], [1, 2])
        self.assertEqual([row["source_membership"] for row in combined], ["v2_only", "legacy_only"])
        self.assertEqual(combined[1]["pair_coordinate_label"], "EBV 1-15 / SELF 20-34*")
        self.assertTrue(all(row["primary_method"] == "tcr_facing_blosum62" for row in combined))


if __name__ == "__main__":
    unittest.main()
