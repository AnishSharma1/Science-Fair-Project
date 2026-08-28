import unittest


from build_rmsd_rerankings import rank_rmsd_rows


def row(group, pair_id, median="", iqr="", status="complete"):
    return {
        "allele": group,
        "pair_id": pair_id,
        "rmsd_status": status,
        "exposed_ca_rmsd_A_median": median,
        "exposed_ca_rmsd_A_iqr": iqr,
    }


class RmsdRankingTests(unittest.TestCase):
    def test_lower_rmsd_ranks_first_within_each_hla(self):
        ranked, missing = rank_rmsd_rows([
            row("A", "a2", "2.0", "0.1"),
            row("A", "a1", "1.0", "0.2"),
            row("B", "b1", "5.0", "0.4"),
            row("A", "missing", status="missing_geometry"),
        ])
        self.assertEqual([item["pair_id"] for item in ranked["A"]], ["a1", "a2"])
        self.assertEqual(ranked["B"][0]["rmsd_rank"], 1)
        self.assertEqual(missing[0]["pair_id"], "missing")
        self.assertEqual(missing[0]["rmsd_rank_status"], "not_ranked_missing_comparable_rmsd")

    def test_equal_medians_share_score_rank_but_iqr_orders_display(self):
        ranked, _missing = rank_rmsd_rows([
            row("A", "wide", "1.0", "0.5"),
            row("A", "tight", "1.0", "0.1"),
        ])
        self.assertEqual([item["pair_id"] for item in ranked["A"]], ["tight", "wide"])
        self.assertEqual([item["rmsd_score_rank"] for item in ranked["A"]], [1, 1])
        self.assertEqual([item["rmsd_score_tie_size"] for item in ranked["A"]], [2, 2])


if __name__ == "__main__":
    unittest.main()
