import csv
import json
from pathlib import Path
import tempfile
import unittest


from high_yield_candidate_evidence import (
    build_stage2_gate,
    classify_assay_evidence,
    classify_ligand_hit,
    classify_stage1,
    normalize_hla,
    presentation_conditioned_rarity,
    scan_similarity_rarity,
    scan_similarity_rarity_fast,
    sequence_relation,
    summarize_conservation,
    summarize_predictor_evidence,
)
from build_high_yield_candidate_evidence import (
    analyze_package,
    parse_netmhcii_batch,
    prepare_package,
)


class ProvenanceClassificationTests(unittest.TestCase):
    def test_hla_normalization_preserves_exact_class_ii_allotype(self):
        self.assertEqual(normalize_hla("DRB1_15_01"), "HLA-DRB1*15:01")
        self.assertEqual(normalize_hla("HLA-DRB1*08:01"), "HLA-DRB1*08:01")

    def test_sequence_relation_distinguishes_exact_overlap_and_none(self):
        self.assertEqual(sequence_relation("ACDEFGHIK", "ACDEFGHIK"), "exact")
        self.assertEqual(sequence_relation("ACDEFGHIK", "YYACDEFGHIKWW"), "overlap")
        self.assertEqual(sequence_relation("ACDEFGHIK", "YYYYYYYYY"), "none")

    def test_exact_sequence_and_exact_hla_is_the_only_exact_hla_class(self):
        base = {
            "epitope_sequence": "ACDEFGHIKLMNPQR",
            "mhc_allele": "HLA-DRB1*15:01",
            "host_organism": "Homo sapiens",
        }
        self.assertEqual(
            classify_assay_evidence(base, "ACDEFGHIKLMNPQR", "HLA-DRB1*15:01"),
            "exact_sequence_exact_hla",
        )
        wrong_hla = {**base, "mhc_allele": "HLA-DRB1*04:01"}
        self.assertEqual(
            classify_assay_evidence(wrong_hla, "ACDEFGHIKLMNPQR", "HLA-DRB1*15:01"),
            "other_human_hla",
        )
        class_i = {**base, "mhc_allele": "HLA-A*02:01"}
        self.assertEqual(
            classify_assay_evidence(class_i, "ACDEFGHIKLMNPQR", "HLA-DRB1*15:01"),
            "class_i",
        )
        mouse = {**base, "host_organism": "Mus musculus", "mhc_allele": "H-2-IAb"}
        self.assertEqual(
            classify_assay_evidence(mouse, "ACDEFGHIKLMNPQR", "HLA-DRB1*15:01"),
            "nonhuman",
        )


class PredictorEvidenceTests(unittest.TestCase):
    def test_batched_netmhcii_response_maps_seq_num_back_to_exact_arm(self):
        text = (
            "allele\tseq_num\tstart\tend\tlength\tcore_peptide\tpeptide\tscore\trank\n"
            "HLA-DRB1*15:01\t1\t1\t15\t15\tACDEFGHIK\tACDEFGHIKLMNPQR\t0.2\t2.0\n"
            "HLA-DRB1*15:01\t2\t1\t15\t15\tCDEFGHIKL\tCDEFGHIKLMNPQRA\t0.1\t3.0\n"
        )
        arms = [
            {"arm_id": "A", "target_id": "T1", "allele": "HLA-DRB1*15:01", "sequence": "ACDEFGHIKLMNPQR"},
            {"arm_id": "B", "target_id": "T2", "allele": "HLA-DRB1*15:01", "sequence": "CDEFGHIKLMNPQRA"},
        ]
        rows = parse_netmhcii_batch(text, "netmhciipan_4_3_el", arms, "raw.tsv")
        self.assertEqual([row["arm_id"] for row in rows], ["A", "B"])
        self.assertEqual([row["core"] for row in rows], ["ACDEFGHIK", "CDEFGHIKL"])

    def test_binding_and_register_consensus_require_both_independent_tools(self):
        rows = [
            {"predictor": "netmhciipan_4_3_el", "percentile_rank": 1.2, "core": "ACDEFGHIK"},
            {"predictor": "netmhciipan_4_3_ba", "percentile_rank": 2.3, "core": "ACDEFGHIK"},
            {"predictor": "mixmhc2pred_2_1_context", "percentile_rank": 3.4, "core": "ACDEFGHIK", "orientation": "-1"},
            {"predictor": "mixmhc2pred_2_1_no_context", "percentile_rank": 4.0, "core": "CDEFGHIKL"},
        ]
        summary = summarize_predictor_evidence(rows, "ACDEFGHIK")
        self.assertTrue(summary["binding_consensus"])
        self.assertTrue(summary["register_consensus"])
        self.assertTrue(summary["register_consensus_matches_declared"])
        self.assertEqual(summary["predictor_status"], "complete")
        self.assertTrue(summary["mixmhc2pred_reverse_orientation"])

    def test_disagreement_and_missing_tools_cannot_become_consensus(self):
        disagreement = summarize_predictor_evidence(
            [
                {"predictor": "netmhciipan_4_3_el", "percentile_rank": 1.0, "core": "ACDEFGHIK"},
                {"predictor": "mixmhc2pred_2_1_context", "percentile_rank": 2.0, "core": "CDEFGHIKL"},
            ],
            "ACDEFGHIK",
        )
        self.assertTrue(disagreement["binding_consensus"])
        self.assertFalse(disagreement["register_consensus"])
        missing = summarize_predictor_evidence(
            [{"predictor": "netmhciipan_4_3_el", "percentile_rank": 1.0, "core": "ACDEFGHIK"}],
            "ACDEFGHIK",
        )
        self.assertEqual(missing["predictor_status"], "not_evaluable")
        self.assertFalse(missing["binding_consensus"])


class LigandAndRarityTests(unittest.TestCase):
    def test_ligand_hits_keep_exact_nested_region_and_absent_distinct(self):
        self.assertEqual(
            classify_ligand_hit("ACDEFGHIK", "ACDEFGHIK", exact_hla=True, monoallelic=True),
            "exact_sequence_monoallelic_exact_hla",
        )
        self.assertEqual(
            classify_ligand_hit("ACDEFGHIK", "YYACDEFGHIKWW", exact_hla=True, monoallelic=False),
            "nested_overlap_multiallelic_compatible",
        )
        self.assertEqual(
            classify_ligand_hit("ACDEFGHIK", "YYYYYYYYY", exact_hla=False, monoallelic=False),
            "not_a_hit",
        )
        self.assertEqual(
            classify_ligand_hit(
                "XXACDEFGHIKYY", "ZZACDEFGHIKZZ", exact_hla=True,
                monoallelic=False, target_core="ACDEFGHIK",
            ),
            "core_overlap_multiallelic_compatible",
        )

    def test_rarity_scan_reports_unrelated_windows_at_least_as_good(self):
        database = [
            {"accession": "P1", "protein": "one", "sequence": "YYACDEFGHIKWW"},
            {"accession": "P2", "protein": "two", "sequence": "YYACDEYGHIKWW"},
            {"accession": "P3", "protein": "three", "sequence": "YYYYYYYYYYYYY"},
        ]
        result = scan_similarity_rarity(
            query_core="ACDEFGHIK",
            paired_core="ACDEFGHIK",
            database_records=database,
            exclude_accession="P1",
            exclude_core="ACDEFGHIK",
            top_n=2,
        )
        self.assertEqual(result["excluded_target_window_count"], 1)
        self.assertGreater(result["evaluated_window_count"], 0)
        self.assertEqual(len(result["nearest_neighbors"]), 2)
        self.assertEqual(result["at_least_as_good_count"], 0)
        self.assertEqual(result["empirical_percentile"], 0.0)

    def test_fast_rarity_scan_matches_reference_implementation(self):
        database = [
            {"accession": "P1", "protein": "one", "sequence": "YYACDEFGHIKWW"},
            {"accession": "P2", "protein": "two", "sequence": "YYACDEYGHIKWW"},
            {"accession": "P3", "protein": "three", "sequence": "YYYYYYYYYYYYY"},
        ]
        arguments = {
            "query_core": "ACDEFGHIK",
            "paired_core": "ACDEYGHIK",
            "database_records": database,
            "exclude_accession": "P2",
            "exclude_core": "ACDEYGHIK",
            "top_n": 3,
        }
        reference = scan_similarity_rarity(**arguments)
        fast = scan_similarity_rarity_fast(**arguments)
        for key in (
            "evaluated_window_count",
            "excluded_target_window_count",
            "at_least_as_good_count",
            "empirical_percentile",
        ):
            self.assertEqual(fast[key], reference[key])
        self.assertEqual(
            [(row["accession"], row["start_1_based"], row["core"]) for row in fast["nearest_neighbors"]],
            [(row["accession"], row["start_1_based"], row["core"]) for row in reference["nearest_neighbors"]],
        )

    def test_presentation_conditioned_rarity_uses_exact_hla_and_excludes_target(self):
        rows = [
            {"allele": "HLA-DRB1*03:01", "pair_id": "better", "tcr_facing_blosum62_similarity": 0.8,
             "tcr_face_physicochemical_mismatch": 0.2, "tcr_facing_sequence_identity": 0.2},
            {"allele": "HLA-DRB1*03:01", "pair_id": "target", "tcr_facing_blosum62_similarity": 0.5,
             "tcr_face_physicochemical_mismatch": 0.2, "tcr_facing_sequence_identity": 0.2},
            {"allele": "HLA-DRB1*03:01", "pair_id": "worse", "tcr_facing_blosum62_similarity": 0.1,
             "tcr_face_physicochemical_mismatch": 0.2, "tcr_facing_sequence_identity": 0.2},
            {"allele": "HLA-DRB1*08:01", "pair_id": "wrong_hla", "tcr_facing_blosum62_similarity": 1.0,
             "tcr_face_physicochemical_mismatch": 0.0, "tcr_facing_sequence_identity": 1.0},
        ]
        result = presentation_conditioned_rarity(
            target_pair_id="target",
            target_allele="DRB1*03:01",
            candidate_rows=rows,
        )
        self.assertEqual(result["status"], "evaluable")
        self.assertEqual(result["candidate_pair_count"], 3)
        self.assertEqual(result["at_least_as_good_count"], 1)
        self.assertEqual(result["empirical_percentile"], 50.0)

    def test_conservation_keeps_missing_coverage_explicit(self):
        complete = summarize_conservation(
            "XXACDEFGHIKYY", "ACDEFGHIK",
            ["ZZXXACDEFGHIKYYZZ", "ZZXXACDEYGHIKYYZZ"],
        )
        self.assertEqual(complete["conservation_status"], "evaluable")
        self.assertEqual(complete["exact_peptide_count"], 1)
        self.assertEqual(complete["exact_core_count"], 1)
        missing = summarize_conservation("XXACDEFGHIKYY", "ACDEFGHIK", [])
        self.assertEqual(missing["conservation_status"], "not_evaluable_missing_sequence_coverage")


class ExperimentalGateTests(unittest.TestCase):
    def arm(self, **overrides):
        base = {
            "predictor_status": "complete",
            "binding_consensus": True,
            "binding_supported": True,
            "both_predictors_above_20": False,
            "register_consensus_matches_declared": True,
            "identity_or_hla_conflict": False,
        }
        return {**base, **overrides}

    def test_stage1_high_medium_and_hold_rules(self):
        self.assertEqual(
            classify_stage1(self.arm(), self.arm(), rarity_percentile=0.4),
            "stage1_high_priority",
        )
        self.assertEqual(
            classify_stage1(
                self.arm(register_consensus_matches_declared=False),
                self.arm(),
                rarity_percentile=0.4,
            ),
            "stage1_medium_priority",
        )
        self.assertEqual(
            classify_stage1(
                self.arm(predictor_status="not_evaluable"),
                self.arm(),
                rarity_percentile=0.4,
            ),
            "stage1_hold",
        )
        self.assertEqual(
            classify_stage1(
                self.arm(both_predictors_above_20=True),
                self.arm(),
                rarity_percentile=0.4,
            ),
            "stage1_hold",
        )

    def test_stage2_gate_is_locked_without_experimental_binding_and_register(self):
        gate = build_stage2_gate([{"target_id": "T1"}, {"target_id": "T2"}])
        self.assertEqual(
            gate["status"],
            "not_evaluable_pending_experimental_binding_and_register",
        )
        self.assertFalse(gate["tcell_assay_recommendation_allowed"])
        self.assertFalse(gate["specificity_claim_allowed"])


class PackageWorkflowTests(unittest.TestCase):
    def test_prepare_and_offline_analyze_build_locked_eight_candidate_package(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            out = Path(temp_dir)
            prepared = prepare_package(output_dir=out)
            self.assertEqual(prepared["target_count"], 8)
            self.assertEqual(prepared["arm_count"], 16)
            self.assertTrue((out / "protocol_lock.json").exists())
            analyzed = analyze_package(output_dir=out)
            self.assertEqual(analyzed["target_count"], 8)
            with (out / "candidate_evidence_matrix.csv").open(newline="", encoding="utf-8") as handle:
                candidates = list(csv.DictReader(handle))
            self.assertEqual(len(candidates), 8)
            with (out / "peptide_arm_evidence.csv").open(newline="", encoding="utf-8") as handle:
                arms = list(csv.DictReader(handle))
            self.assertEqual(len(arms), 16)
            with (out / "stage2_tcell_gate.json").open(encoding="utf-8") as handle:
                gate = json.load(handle)
            self.assertFalse(gate["tcell_assay_recommendation_allowed"])
            self.assertTrue((out / "SHA256SUMS.csv").exists())
            self.assertTrue((out / "raw_response_manifest.csv").exists())


if __name__ == "__main__":
    unittest.main()
