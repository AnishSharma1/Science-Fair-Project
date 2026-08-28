import csv
import tempfile
import unittest
from pathlib import Path


from build_same_register_hla_rankings_v2 import (
    NONSTRUCTURAL_METHODS,
    rank_within_hla,
    run,
    select_exact_target_epitopes,
    select_control_supported_method,
)


def control_row(method, system, rank, pair="PAIR", seed=271828):
    return {
        "method": method,
        "system_id": system,
        "positive_pair_id": pair,
        "panel_seed": str(seed),
        "positive_rank": str(rank),
        "comparison_count": "26",
    }


def pair(allele, pair_id, left, right, geometry_status="complete"):
    return {
        "allele": allele,
        "ebv_allele": allele,
        "self_allele": allele,
        "pair_id": pair_id,
        "ebv_candidate_id": f"E_{pair_id}",
        "ebv_protein": "EBV",
        "ebv_sequence": f"AAA{left}AAA",
        "ebv_predicted_core": left,
        "ebv_binding_percentile_rank": "1.0",
        "ebv_source_certainty": "exact",
        "self_candidate_id": f"S_{pair_id}",
        "self_protein": "SELF",
        "self_sequence": f"AAA{right}AAA",
        "self_predicted_core": right,
        "self_binding_percentile_rank": "2.0",
        "self_source_certainty": "exact",
        "geometry_status": geometry_status,
    }


class MethodSelectionTests(unittest.TestCase):
    def test_control_selection_uses_one_worst_rank_per_system(self):
        rows = []
        ranks = {
            "physicochemical_only": {"A": (1, 1), "B": (2, 2), "C": (1, 1)},
            "tcr_facing_identity": {"A": (1, 1), "B": (2, 2), "C": (2, 2)},
            "full_core_identity": {"A": (2, 2), "B": (4, 4), "C": (3, 3)},
            "tcr_facing_blosum62": {"A": (1, 1), "B": (1, 1), "C": (1, 1)},
            "full_core_blosum62": {"A": (2, 2), "B": (1, 1), "C": (2, 2)},
        }
        for method, systems in ranks.items():
            for system, seed_ranks in systems.items():
                for seed, rank in zip((271828, 314159), seed_ranks):
                    rows.append(control_row(method, system, rank, seed=seed))
        selected, evidence = select_control_supported_method(rows)
        self.assertEqual(selected, "tcr_facing_blosum62")
        selected_row = next(row for row in evidence if row["method"] == selected)
        self.assertEqual(selected_row["system_worst_ranks"], "A:1;B:1;C:1")
        self.assertEqual(selected_row["system_capture_at_3_count"], 3)

    def test_selection_rejects_missing_methods(self):
        rows = [control_row("tcr_facing_blosum62", "A", 1)]
        with self.assertRaisesRegex(ValueError, "missing nonstructural methods"):
            select_control_supported_method(rows)


class RankingTests(unittest.TestCase):
    def test_hla_ranking_is_separate_and_includes_geometry_missing_pairs(self):
        rows = [
            pair("HLA-DRB1*15:01", "a_best", "AAAAAAAAA", "AAAAAAAAA"),
            pair("HLA-DRB1*15:01", "a_worse", "AAAAAAAAA", "RRRRRRRRR", "missing_or_qc_excluded_model"),
            pair("HLA-DRB1*03:01", "b_only", "RRRRRRRRR", "RRRRRRRRR"),
        ]
        ranked = rank_within_hla(rows, "tcr_facing_blosum62")
        self.assertEqual([row["pair_id"] for row in ranked["HLA-DRB1*15:01"]], ["a_best", "a_worse"])
        self.assertEqual(ranked["HLA-DRB1*03:01"][0]["hla_rank"], 1)
        missing_geometry = ranked["HLA-DRB1*15:01"][1]
        self.assertEqual(missing_geometry["legacy_geometry_status"], "missing_or_qc_excluded_model")
        self.assertEqual(missing_geometry["ranking_input_status"], "complete_sequence_register")
        self.assertFalse(missing_geometry["geometry_used_in_ranking"])
        self.assertEqual(missing_geometry["rank_scope"], "within_hla_only")
        self.assertEqual(missing_geometry["computational_pair_marker"], "*")
        self.assertEqual(
            missing_geometry["pair_evidence_status"],
            "computational_pair_no_exact_paired_recognition_evidence",
        )

    def test_equal_primary_scores_report_ties_with_deterministic_display_order(self):
        rows = [
            pair("HLA-DRB1*15:01", "z_pair", "AAAAAAAAA", "AAAAAAAAA"),
            pair("HLA-DRB1*15:01", "a_pair", "AAAAAAAAA", "AAAAAAAAA"),
        ]
        ranked = rank_within_hla(rows, "tcr_facing_blosum62")["HLA-DRB1*15:01"]
        self.assertEqual([row["pair_id"] for row in ranked], ["a_pair", "z_pair"])
        self.assertEqual([row["hla_score_rank"] for row in ranked], [1, 1])
        self.assertEqual([row["hla_score_tie_size"] for row in ranked], [2, 2])

    def test_invalid_register_or_cross_hla_pair_is_rejected(self):
        invalid = pair("HLA-DRB1*15:01", "bad", "AAAAAAAA", "AAAAAAAAA")
        with self.assertRaisesRegex(ValueError, "nine-residue"):
            rank_within_hla([invalid], "tcr_facing_blosum62")
        cross = pair("HLA-DRB1*15:01", "cross", "AAAAAAAAA", "AAAAAAAAA")
        cross["self_allele"] = "HLA-DRB1*03:01"
        with self.assertRaisesRegex(ValueError, "exact-HLA"):
            rank_within_hla([cross], "tcr_facing_blosum62")

    def test_selected_target_names_exact_peptides_and_registers(self):
        rows = [
            {
                **pair("HLA-DRB1*15:01", "target", "LRALLARSH", "LEARLSRMH"),
                "ebv_protein": "EBNA1",
                "self_protein": "ANO2",
            },
            {
                **pair("HLA-DRB1*15:01", "other", "RRRRRRRRR", "AAAAAAAAA"),
                "ebv_protein": "OTHER",
                "self_protein": "OTHER",
            },
        ]
        ranked = rank_within_hla(rows, "tcr_facing_blosum62")
        selected = select_exact_target_epitopes(ranked)
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["ebv_epitope_15mer"], "AAALRALLARSHAAA")
        self.assertEqual(selected[0]["ebv_core_p1_p9"], "LRALLARSH")
        self.assertEqual(selected[0]["self_core_p1_p9"], "LEARLSRMH")
        self.assertEqual(selected[0]["selection_status"], "selected_rank_1_within_hla")
        self.assertTrue(selected[0]["prospective_not_experimentally_confirmed_pair"])


class PackageTests(unittest.TestCase):
    def test_run_emits_provisional_reproducible_hla_specific_package(self):
        control_rows = []
        for method in NONSTRUCTURAL_METHODS:
            for system in ("A", "B", "C"):
                rank = 1 if method == "tcr_facing_blosum62" else 2
                for seed in (271828, 314159):
                    control_rows.append(control_row(method, system, rank, seed=seed))
        discovery_rows = [
            pair("HLA-DRB1*15:01", "p1", "AAAAAAAAA", "AAAAAAAAA"),
            pair("HLA-DRB1*15:01", "p2", "AAAAAAAAA", "RRRRRRRRR"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            discovery = root / "pair_summary.csv"
            controls = root / "method_rank_long.csv"
            panel = root / "frozen_panel.csv"
            out = root / "out"
            for path, rows in ((discovery, discovery_rows), (controls, control_rows)):
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
            panel_rows = []
            for index, row in enumerate(discovery_rows, start=1):
                panel_rows.extend([
                    {
                        "candidate_id": row["ebv_candidate_id"],
                        "sequence": row["ebv_sequence"],
                        "source_start_1_based": str(index * 10),
                        "source_end_1_based": str(index * 10 + len(row["ebv_sequence"]) - 1),
                    },
                    {
                        "candidate_id": row["self_candidate_id"],
                        "sequence": row["self_sequence"],
                        "source_start_1_based": str(index * 20),
                        "source_end_1_based": str(index * 20 + len(row["self_sequence"]) - 1),
                    },
                ])
            with panel.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(panel_rows[0]))
                writer.writeheader()
                writer.writerows(panel_rows)
            first = run(
                discovery_path=discovery, benchmark_path=controls, panel_path=panel, out=out
            )
            first_checksums = (out / "SHA256SUMS.csv").read_bytes()
            second = run(
                discovery_path=discovery, benchmark_path=controls, panel_path=panel, out=out
            )
            self.assertEqual(first, second)
            self.assertEqual(first_checksums, (out / "SHA256SUMS.csv").read_bytes())
            self.assertEqual(first["selected_primary_method"], "tcr_facing_blosum62")
            self.assertEqual(first["ranked_pair_count"], 2)
            self.assertFalse(first["definitive_validation_complete"])
            self.assertFalse(first["discovery_unlock_allowed"])
            self.assertTrue((out / "rankings/hla_drb1_15_01_ranked_pairs.csv").exists())
            with (out / "top_10_exact_epitopes_by_hla.csv").open(newline="", encoding="utf-8") as handle:
                exact_top = list(csv.DictReader(handle))
            self.assertEqual(len(exact_top), 2)
            self.assertTrue(all(row["computational_pair_marker"] == "*" for row in exact_top))
            self.assertEqual(exact_top[0]["ebv_epitope_sequence"], "AAAAAAAAAAAAAAA")
            self.assertEqual(exact_top[0]["ebv_core_p1_p9"], "AAAAAAAAA")
            self.assertEqual(exact_top[0]["ebv_epitope_label"], "EBV 10-24")
            self.assertEqual(exact_top[0]["self_epitope_label"], "SELF 20-34")
            self.assertEqual(exact_top[0]["pair_coordinate_label"], "EBV 10-24 / SELF 20-34*")


if __name__ == "__main__":
    unittest.main()
