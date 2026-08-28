import unittest

from lead_focused_robustness import (
    cluster_complete_linkage,
    classify_separation,
    deduplicate_jobs,
    empirical_target_rank,
    hierarchical_technical_bootstrap,
    require_single_analysis_layer,
    summarize_lead,
)


class DeduplicateJobsTests(unittest.TestCase):
    def test_duplicate_request_and_ordered_model_hashes_retains_one_canonical_job(self):
        jobs = [
            {
                "path": "/download/one",
                "request_identity": "target-A",
                "model_hashes": ["a", "b", "c", "d", "e"],
            },
            {
                "path": "/download/copied",
                "request_identity": "target-A",
                "model_hashes": ["a", "b", "c", "d", "e"],
            },
        ]

        result = deduplicate_jobs(jobs)

        self.assertEqual([job["path"] for job in result["retained_jobs"]], ["/download/one"])
        self.assertEqual(
            result["duplicate_paths_by_canonical"],
            {("target-A", ("a", "b", "c", "d", "e")): ["/download/copied"]},
        )

    def test_distinct_model_hashes_survive_even_for_the_same_control(self):
        jobs = [
            {"path": "/download/one", "request_identity": "control-A", "model_hashes": ["a", "b", "c", "d", "e"]},
            {"path": "/download/two", "request_identity": "control-A", "model_hashes": ["a", "b", "c", "d", "changed"]},
        ]

        result = deduplicate_jobs(jobs)

        self.assertEqual([job["path"] for job in result["retained_jobs"]], ["/download/one", "/download/two"])

    def test_same_path_distinct_jobs_keep_separate_duplicate_audit_buckets(self):
        first_hashes = ["a", "b", "c", "d", "e"]
        second_hashes = ["f", "g", "h", "i", "j"]
        jobs = [
            {"path": "/download/shared", "request_identity": "control-A", "model_hashes": first_hashes},
            {"path": "/download/shared", "request_identity": "control-A", "model_hashes": second_hashes},
            {"path": "/download/copy-one", "request_identity": "control-A", "model_hashes": first_hashes},
            {"path": "/download/copy-two", "request_identity": "control-A", "model_hashes": second_hashes},
        ]

        result = deduplicate_jobs(jobs)

        self.assertEqual(len(result["retained_jobs"]), 2)
        self.assertEqual(
            result["duplicate_paths_by_canonical"],
            {
                ("control-A", tuple(first_hashes)): ["/download/copy-one"],
                ("control-A", tuple(second_hashes)): ["/download/copy-two"],
            },
        )

    def test_duplicate_control_job_cannot_add_background_weight_after_deduplication(self):
        jobs = [
            {"path": "/target", "request_identity": "target", "model_hashes": ["t"] * 5, "entity": "target", "median": 1.0},
            {"path": "/control-original", "request_identity": "control-A", "model_hashes": ["a"] * 5, "entity": "control-A", "median": 9.0},
            {"path": "/control-copy", "request_identity": "control-A", "model_hashes": ["a"] * 5, "entity": "control-A", "median": 9.0},
            {"path": "/other-control", "request_identity": "control-B", "model_hashes": ["b"] * 5, "entity": "control-B", "median": 5.0},
        ]
        retained = deduplicate_jobs(jobs)["retained_jobs"]
        control_jobs = {"control-A": [], "control-B": []}
        target_jobs = []
        for job in retained:
            if job["entity"] == "target":
                target_jobs.append(job["median"])
            else:
                control_jobs[job["entity"]].append(job["median"])

        summary = summarize_lead(target_jobs, control_jobs)

        self.assertEqual(summary["background_median"], 7.0)

    def test_request_identity_requires_all_five_ordered_model_hashes(self):
        incomplete_job = [{"path": "/partial", "request_identity": "control-A", "model_hashes": ["a", "b", "c", "d"]}]

        with self.assertRaisesRegex(ValueError, "five"):
            deduplicate_jobs(incomplete_job)


class LeadSummaryTests(unittest.TestCase):
    def test_controls_receive_one_equal_weight_top_level_median(self):
        summary = summarize_lead(
            target_job_medians=[1.0, 3.0],
            control_job_medians={"control_many_jobs": [3.0, 3.0, 3.0, 3.0, 99.0], "control_one_job": [7.0]},
        )

        self.assertEqual(summary["target_median"], 2.0)
        self.assertEqual(summary["control_medians"], {"control_many_jobs": 3.0, "control_one_job": 7.0})
        self.assertEqual(summary["background_median"], 5.0)
        self.assertEqual(summary["background_minus_target_delta"], 3.0)
        self.assertEqual(summary["leave_one_control_out_deltas"], {"control_many_jobs": 5.0, "control_one_job": 1.0})


class EmpiricalRankTests(unittest.TestCase):
    def test_ties_count_against_lower_is_better_target(self):
        result = empirical_target_rank(
            target_identifier="lead",
            target_median=2.0,
            control_medians={"z_control": 2.0, "a_control": 1.0, "b_control": 4.0},
        )

        self.assertEqual(result["sorted_identifiers"], ["a_control", "lead", "z_control", "b_control"])
        self.assertEqual(result["target_rank"], 2)
        self.assertEqual(result["one_sided_tail_fraction"], 0.75)


class CompleteLinkageTests(unittest.TestCase):
    def test_distance_exactly_at_threshold_is_in_the_same_cluster(self):
        labels = cluster_complete_linkage(["model-b", "model-a"], [[0.0, 2.0], [2.0, 0.0]], threshold=2.0)

        self.assertEqual(labels, {"model-a": 1, "model-b": 1})

    def test_threshold_and_labels_are_stable_and_ignore_confidence_annotations(self):
        model_ids = ["model-c", "model-a", "model-b"]
        matrix = [
            [0.0, 2.0, 2.1],
            [2.0, 0.0, 1.0],
            [2.1, 1.0, 0.0],
        ]

        without_confidence = cluster_complete_linkage(model_ids, matrix, threshold=2.0)
        with_confidence = cluster_complete_linkage(
            model_ids,
            matrix,
            threshold=2.0,
            confidence_by_model={"model-a": 1.0, "model-b": 99.0, "model-c": 50.0},
        )

        self.assertEqual(without_confidence, {"model-a": 1, "model-b": 1, "model-c": 2})
        self.assertEqual(with_confidence, without_confidence)


class ClassificationTests(unittest.TestCase):
    def test_consistent_and_nonpositive_rules_require_their_declared_conditions(self):
        consistent_summary = {
            "background_minus_target_delta": 3.0,
            "leave_one_control_out_deltas": {"control-a": 1.0, "control-b": 2.0},
            "target_job_medians": [1.0, 2.0],
            "background_median": 4.0,
        }

        self.assertEqual(classify_separation(consistent_summary, {"percentile_2_5": 0.01}), "consistent_positive")
        self.assertEqual(classify_separation(consistent_summary, {"percentile_2_5": 0.0}), "mixed_positive")
        self.assertEqual(
            classify_separation({**consistent_summary, "background_minus_target_delta": 0.0}, {"percentile_2_5": 1.0}),
            "no_positive_separation",
        )

    def test_rank_two_mixed_result_has_the_required_sensitivity_prefix(self):
        summary = {
            "background_minus_target_delta": 3.0,
            "leave_one_control_out_deltas": {"control-a": 1.0, "control-b": -0.1},
            "target_job_medians": [1.0, 1.5],
            "background_median": 4.0,
        }

        label = classify_separation(summary, {"percentile_2_5": 0.5}, rank=2)

        self.assertEqual(label, "length_sensitivity_only__mixed_positive")

    def test_primary_and_sensitivity_records_cannot_be_pooled(self):
        records = [
            {"analysis_layer": "primary", "job_id": "one"},
            {"analysis_layer": "sensitivity_rank_2", "job_id": "two"},
        ]

        with self.assertRaisesRegex(ValueError, "cannot be pooled"):
            require_single_analysis_layer(records)


class HierarchicalBootstrapTests(unittest.TestCase):
    def test_each_comparison_job_occurrence_gets_a_nested_ebv_job_median(self):
        def job(job_id):
            return {
                "job_id": job_id,
                "model_ids": [f"{job_id}-model-{index}" for index in range(5)],
            }

        ebv_jobs = [job("ebv-1"), job("ebv-2")]
        target_jobs = [job("target-1"), job("target-2")]
        control_jobs = [job("control-1"), job("control-2")]
        entities = {
            "ebv": ebv_jobs,
            "target": target_jobs,
            "controls": {"control": control_jobs},
        }
        job_pair_values = {
            ("ebv-1", "target-1"): 0.0,
            ("ebv-2", "target-1"): 100.0,
            ("ebv-1", "target-2"): 1.0,
            ("ebv-2", "target-2"): 2.0,
        }
        geometry = {}
        for ebv_job in ebv_jobs:
            for comparison_job in [*target_jobs, *control_jobs]:
                value = job_pair_values.get(
                    (ebv_job["job_id"], comparison_job["job_id"]), 200.0
                )
                for ebv_model in ebv_job["model_ids"]:
                    for comparison_model in comparison_job["model_ids"]:
                        geometry[(ebv_model, comparison_model)] = value

        result = hierarchical_technical_bootstrap(
            entities, geometry, iterations=1, seed=1
        )

        self.assertEqual(
            result["replicate_rows"],
            [{
                "iteration": 1,
                "target_median_A": 25.75,
                "equal_weight_background_median_A": 200.0,
                "delta_A": 174.25,
            }],
        )
        self.assertNotEqual(result["replicate_rows"][0]["target_median_A"], 1.5)

    def test_bootstrap_resamples_five_models_with_replacement_inside_a_job(self):
        entities = {
            "ebv": [{"job_id": "ebv", "model_ids": ["e1", "e2", "e3", "e4", "e5"]}],
            "target": [{"job_id": "target", "model_ids": ["t1", "t2", "t3", "t4", "t5"]}],
            "controls": {"control": [{"job_id": "control", "model_ids": ["c1", "c2", "c3", "c4", "c5"]}]},
        }
        geometry = {}
        for index, ebv_model in enumerate(entities["ebv"][0]["model_ids"], start=1):
            for target_model in entities["target"][0]["model_ids"]:
                geometry[(ebv_model, target_model)] = float(index)
            for control_model in entities["controls"]["control"][0]["model_ids"]:
                geometry[(ebv_model, control_model)] = float(index * 2)

        result = hierarchical_technical_bootstrap(entities, geometry)

        self.assertIn(3.0, result["replicates"])
        self.assertTrue(any(delta != 3.0 for delta in result["replicates"]))

    def test_bootstrap_is_deterministic_and_resamples_jobs_before_models(self):
        entities = {
            "ebv": [
                {"job_id": "ebv-one", "model_ids": ["e1a", "e1b", "e1c", "e1d", "e1e"]},
                {"job_id": "ebv-two", "model_ids": ["e2a", "e2b", "e2c", "e2d", "e2e"]},
            ],
            "target": [{"job_id": "target", "model_ids": ["ta", "tb", "tc", "td", "te"]}],
            "controls": {"control": [{"job_id": "control", "model_ids": ["ca", "cb", "cc", "cd", "ce"]}]},
        }
        geometry = {}
        for ebv_model in entities["ebv"][0]["model_ids"]:
            for target_model in entities["target"][0]["model_ids"]:
                geometry[(ebv_model, target_model)] = 1.0
            for control_model in entities["controls"]["control"][0]["model_ids"]:
                geometry[(ebv_model, control_model)] = 10.0
        for ebv_model in entities["ebv"][1]["model_ids"]:
            for target_model in entities["target"][0]["model_ids"]:
                geometry[(ebv_model, target_model)] = 9.0
            for control_model in entities["controls"]["control"][0]["model_ids"]:
                geometry[(ebv_model, control_model)] = 20.0

        first = hierarchical_technical_bootstrap(entities, geometry)
        second = hierarchical_technical_bootstrap(entities, geometry)

        self.assertEqual(first, second)
        self.assertEqual(len(first["replicates"]), 10_000)
        self.assertEqual(set(first["replicates"]), {9.0, 10.0, 11.0})
        self.assertEqual(first["fraction_positive"], 1.0)
