import csv
import hashlib
import unittest
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED = PROJECT_ROOT / "processed"
AUDIT = PROCESSED / "lead_focused_robustness_2026-08-15"
LEDGER = PROCESSED / "publication_evidence_ledger.md"
CLAIM_MATRIX = PROCESSED / "publication_claim_matrix.csv"
READINESS = PROCESSED / "same_register_af3_analysis" / "PUBLICATION_READINESS.md"

EXPECTED_AUDIT_FILES = {
    "LEAD_FOCUSED_FINDINGS.md",
    "control_rank_and_leave_one_out.csv",
    "frozen_input_checksums.csv",
    "job_pair_stability.csv",
    "model_job_identity_manifest.csv",
    "per_control_geometry_summary.csv",
    "pose_cluster_membership.csv",
    "rank1_primary_control_robustness.svg",
    "rank2_length_sensitivity_job_dependence.svg",
    "technical_bootstrap_replicates.csv",
    "technical_bootstrap_summary.csv",
}
FROZEN_HASHES = {
    "complete_model_pipeline_audit_2026-08-12": "cf02fe6d2111dda133eb6e9cd5d942e74cd08f70ad8b47399cf1561e60bf0be4",
    "complete_model_pipeline_audit_2026-08-15": "198e195d4352884d7873daf8810f1353666b014112590c7b6abd9d1571d7d769",
    "structural_control_expansion_2026-08-15": "a7218643435537113d05ce27f5286dcecbc506add039fb8b891a37ff6a028324",
}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def directory_sha256(directory):
    digest = hashlib.sha256()
    for path in sorted(item for item in directory.rglob("*") if item.is_file()):
        digest.update(path.relative_to(directory).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


class LeadFocusedAuditOutputTests(unittest.TestCase):
    def test_audit_has_exact_output_set_and_twenty_thousand_bootstrap_rows(self):
        self.assertEqual({path.name for path in AUDIT.iterdir()}, EXPECTED_AUDIT_FILES)

        rows = read_csv(AUDIT / "technical_bootstrap_replicates.csv")
        self.assertEqual(len(rows), 20_000)
        self.assertEqual(Counter(row["lead_rank"] for row in rows), {"1": 10_000, "2": 10_000})

    def test_rank_one_is_strict_primary_and_rank_two_is_unpooled_length_sensitivity(self):
        rows = read_csv(AUDIT / "control_rank_and_leave_one_out.csv")
        self.assertEqual(rows[0]["analysis_layer"], "strict_primary_controls")
        self.assertEqual(rows[0]["classification"], "consistent_positive")
        self.assertEqual(rows[1]["analysis_layer"], "length_sensitivity_exact_bin_pm7")
        self.assertTrue(rows[1]["classification"].startswith("length_sensitivity_only__"))

        findings = (AUDIT / "LEAD_FOCUSED_FINDINGS.md").read_text(encoding="utf-8")
        self.assertIn("The two layers were not pooled", findings)

    def test_tail_fractions_are_exploratory_and_output_headers_never_name_a_p_value(self):
        rank_rows = read_csv(AUDIT / "control_rank_and_leave_one_out.csv")
        self.assertTrue(all(float(row["exploratory_empirical_tail_fraction"]) >= 0.25 for row in rank_rows))
        self.assertTrue(all("not a p-value" in row["inference_label"] for row in rank_rows))

        for csv_path in AUDIT.glob("*.csv"):
            with csv_path.open(newline="", encoding="utf-8") as handle:
                headers = next(csv.reader(handle))
            self.assertFalse(any("p_value" in header.lower() or "p-value" in header.lower() for header in headers))

    def test_all_materials_keep_the_pmhc_only_claim_boundary(self):
        materials = [
            AUDIT / "LEAD_FOCUSED_FINDINGS.md",
            LEDGER,
            READINESS,
            CLAIM_MATRIX,
        ]
        for path in materials:
            with self.subTest(path=path.name):
                text = path.read_text(encoding="utf-8").lower()
                self.assertIn("pmhc", text)
                self.assertTrue(
                    "does not establish" in text or "do not establish" in text or "not a standalone proof" in text
                )
                self.assertNotIn("demonstrates tcr binding", text)
                self.assertNotIn("demonstrates activation", text)
                self.assertNotIn("demonstrates cross-reactivity", text)
                self.assertNotIn("demonstrates molecular mimicry", text)
                self.assertNotIn("demonstrates an ms mechanism", text)

    def test_claim_matrix_inserts_c06b_preserves_discovery_context_and_routes_c08b_to_supplement(self):
        rows = read_csv(CLAIM_MATRIX)
        claim_ids = [row["claim_id"] for row in rows]
        c06_index = claim_ids.index("C06")
        self.assertEqual(claim_ids[c06_index + 1], "C06B")

        c06b = rows[c06_index + 1]
        self.assertEqual(c06b["manuscript_section"], "Results 2")
        self.assertEqual(c06b["evidence_class"], "computational_structural_control")
        self.assertEqual(c06b["status"], "publication_ready_with_caveat")
        self.assertIn("strict", c06b["allowed_wording"].lower())
        self.assertIn("length", c06b["allowed_wording"].lower())
        self.assertIn("lead_focused_robustness_2026-08-15", c06b["key_artifacts"])
        self.assertIn("TCR binding", c06b["prohibited_overclaim"])

        c06 = rows[c06_index]
        self.assertIn("BALF5-family", c06["proposed_claim"])
        c08b = next(row for row in rows if row["claim_id"] == "C08B")
        self.assertIn("supplement", c08b["next_action"].lower())
        self.assertNotIn("Figure 4", c08b["next_action"])

    def test_publication_materials_report_exact_lead_results_and_retire_absent_controls_language(self):
        ledger = " ".join(LEDGER.read_text(encoding="utf-8").split())
        readiness = " ".join(READINESS.read_text(encoding="utf-8").split())
        bootstrap_rows = {
            row["lead_rank"]: row
            for row in read_csv(AUDIT / "technical_bootstrap_summary.csv")
        }
        self.assertAlmostEqual(float(bootstrap_rows["1"]["delta_percentile_2_5_A"]), 1.620233109253708)
        self.assertAlmostEqual(float(bootstrap_rows["1"]["delta_percentile_97_5_A"]), 12.758541333935227)
        self.assertAlmostEqual(float(bootstrap_rows["2"]["delta_percentile_2_5_A"]), -0.5167436997405864)
        self.assertAlmostEqual(float(bootstrap_rows["2"]["delta_percentile_97_5_A"]), 16.879650260161725)
        interval_phrases = {
            rank: (
                f'{float(row["delta_percentile_2_5_A"]):.3f} to '
                f'{float(row["delta_percentile_97_5_A"]):.3f} A'
            )
            for rank, row in bootstrap_rows.items()
        }
        c06b_wording = next(
            row["allowed_wording"]
            for row in read_csv(CLAIM_MATRIX)
            if row["claim_id"] == "C06B"
        )
        for interval in interval_phrases.values():
            self.assertIn(interval, ledger)
            self.assertIn(interval, readiness)
            self.assertIn(interval, c06b_wording)

        expected_phrases = (
            "7.321 A",
            "0.25",
            "length_sensitivity_only__mixed_positive",
            "Pose/job consistency",
        )
        for phrase in expected_phrases:
            self.assertIn(phrase, ledger)
        self.assertIn("rank-1 primary-control panel in the main results", readiness.lower())
        self.assertIn("rank-2 sensitivity/job-dependence and pytorch classifier in the supplement", readiness.lower())
        self.assertNotIn("no completed structural controls", readiness.lower())
        self.assertNotIn("completed matched structural controls are absent", readiness.lower())

    def test_frozen_discovery_order_and_frozen_directory_hashes_are_unchanged(self):
        frozen_rows = read_csv(PROCESSED / "complete_model_pipeline_audit_2026-08-15" / "master_pair_score_sheet.csv")
        expanded_rows = read_csv(PROCESSED / "structural_control_expansion_2026-08-15" / "master_pair_score_sheet_with_expanded_controls.csv")
        frozen_bytes = b"".join(
            f'{row["discovery_priority_rank"]},{row["pair_id"]}\n'.encode("utf-8") for row in frozen_rows
        )
        expanded_bytes = b"".join(
            f'{row["discovery_priority_rank"]},{row["pair_id"]}\n'.encode("utf-8") for row in expanded_rows
        )
        self.assertEqual(expanded_bytes, frozen_bytes)

        for directory_name, expected_hash in FROZEN_HASHES.items():
            with self.subTest(directory=directory_name):
                self.assertEqual(directory_sha256(PROCESSED / directory_name), expected_hash)


if __name__ == "__main__":
    unittest.main()
