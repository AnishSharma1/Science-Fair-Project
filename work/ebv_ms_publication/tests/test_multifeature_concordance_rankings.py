import unittest


from build_multifeature_concordance_rankings import (
    _control_panel_audit,
    build_concordance_ranking,
)


def row(pair_id, blosum, chemistry, identity, exposed, full_core, status="complete"):
    return {
        "allele": "HLA-DRB1*15:01",
        "pair_id": pair_id,
        "pair_coordinate_label": f"{pair_id}*",
        "tcr_facing_blosum62_similarity": blosum,
        "tcr_face_physicochemical_mismatch": chemistry,
        "tcr_facing_sequence_identity": identity,
        "full_core_sequence_identity": identity,
        "exposed_ca_rmsd_A_median": exposed,
        "full_core_ca_rmsd_A_median": full_core,
        "anchor_ca_rmsd_A_median": full_core,
        "multifeature_input_status": status,
    }


class ConcordanceRankingTests(unittest.TestCase):
    def test_jointly_good_pair_beats_one_sided_pairs(self):
        rows = [
            row("joint", 0.7, 0.1, 0.8, 0.5, 0.5),
            row("sequence_only", 0.9, 0.0, 1.0, 5.0, 5.0),
            row("structure_only", 0.0, 0.5, 0.0, 0.1, 0.1),
        ]
        ranked, missing = build_concordance_ranking(rows)
        self.assertEqual(missing, [])
        self.assertEqual(ranked[0]["pair_id"], "joint")
        self.assertEqual(ranked[0]["concordance_rank"], 1)
        self.assertLessEqual(ranked[0]["sequence_family_percentile"], 0.5)
        self.assertLessEqual(ranked[0]["structure_family_percentile"], 0.5)

    def test_feature_directions_and_ties_are_explicit(self):
        rows = [
            row("a", 1.0, 0.0, 1.0, 1.0, 1.0),
            row("b", 1.0, 0.0, 1.0, 1.0, 1.0),
            row("c", 0.0, 1.0, 0.0, 2.0, 2.0),
        ]
        ranked, _missing = build_concordance_ranking(rows)
        by_id = {item["pair_id"]: item for item in ranked}
        self.assertEqual(by_id["a"]["tcr_facing_blosum62_percentile"], 0.25)
        self.assertEqual(by_id["b"]["tcr_facing_blosum62_percentile"], 0.25)
        self.assertEqual(by_id["c"]["tcr_facing_blosum62_percentile"], 1.0)
        self.assertEqual([item["pair_id"] for item in ranked[:2]], ["a", "b"])

    def test_missing_structure_cannot_enter_composite_rank(self):
        ranked, missing = build_concordance_ranking([
            row("complete", 0.5, 0.2, 0.4, 1.0, 1.0),
            row("missing", 1.0, 0.0, 1.0, "", "", status="missing_structure"),
        ])
        self.assertEqual([item["pair_id"] for item in ranked], ["complete"])
        self.assertEqual(missing[0]["concordance_rank_status"], "not_ranked_missing_comparable_structure")

    def test_control_audit_distinguishes_declared_and_seed_qualified_positive_ids(self):
        rows = []
        for index in range(26):
            is_positive = index == 0
            rows.append({
                "system_id": "system",
                "positive_pair_id": "PAIR_DECLARED_POSITIVE",
                "panel_seed": "271828",
                "pair_id": "s271828|positive" if is_positive else f"s271828|decoy_{index:02d}",
                "pair_role": "positive" if is_positive else "N3",
                "geometry_status": "complete",
                "tcr_facing_blosum62_similarity": 1.0 if is_positive else 0.0,
                "tcr_face_physicochemical_mismatch_median": 0.0 if is_positive else 1.0,
                "tcr_facing_sequence_identity": 1.0 if is_positive else 0.0,
                "full_core_sequence_identity": 1.0 if is_positive else 0.0,
                "exposed_ca_rmsd_A_median": 0.1 if is_positive else 2.0,
                "anchor_ca_rmsd_A_median": 0.1 if is_positive else 2.0,
            })
        panel_rows, audit = _control_panel_audit(rows)
        self.assertEqual(panel_rows[0]["positive_concordance_rank"], 1)
        self.assertEqual(audit["panel_capture_at_3_count"], 1)


if __name__ == "__main__":
    unittest.main()
