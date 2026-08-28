import csv
import hashlib
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from lead_focused_robustness import (
    build_audit_tables,
    build_identity_manifest,
    build_live_lead_inputs,
    collect_input_checksums,
    discover_saved_jobs,
    generate_live_audit,
    hierarchical_technical_bootstrap,
    load_fixed_lead_metadata,
    render_rank1_svg,
    render_rank2_svg,
    validate_lead_definition,
    write_audit_outputs,
)


_TEST_CANDIDATES = {
    "EBV_TCELL_950": ("GGACDEFGHIKLMGG", "ACDEFGHIK"),
    "HUMAN_MYELIN_112214": ("TTLMNPQRSTVWYTT", "LMNPQRSTV"),
    "HUMAN_BACKGROUND_115891": ("VVDEFGHIKLMNPVV", "DEFGHIKLM"),
    "HUMAN_BACKGROUND_118550": ("WWPQRSTVWYACDWW", "PQRSTVWYA"),
    "HUMAN_BACKGROUND_119732": ("YYHIKLMNPQRSTYY", "HIKLMNPQR"),
    "EBV_TCELL_2268741": ("GGFGHIKLMNPQRSS", "FGHIKLMNP"),
    "HUMAN_MYELIN_117032": ("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN", "RSTVWYACD"),
    "HUMAN_BACKGROUND_141561": ("AANPQRSTVWYACAA", "NPQRSTVWY"),
    "HUMAN_BACKGROUND_423369": ("CCKLMNPQRSTVWCC", "KLMNPQRST"),
    "HUMAN_BACKGROUND_2258889": ("DDRSTVWYACDEFGDD", "RSTVWYACD"),
}
_AA3 = {
    one: three for one, three in zip(
        "ACDEFGHIKLMNPQRSTVWY",
        ("ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"),
    )
}
_HLA_SEQUENCE = "A" * 85


def _write_fixture_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _fixture_cif(peptide, candidate_offset, model_index, signature_marker=""):
    lines = [
        "data_fixture",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.B_iso_or_equiv",
    ]
    for chain, y_base, modulus in (("A", 0.0, 5), ("B", 10.0, 7)):
        for index in range(1, 86):
            lines.append(
                f"ATOM C CA ALA {chain} {index} {index:.3f} {y_base:.3f} {(index % modulus) * 0.1:.3f} 70.0"
            )
    for index, amino_acid in enumerate(peptide, start=1):
        offset = candidate_offset + model_index * 0.05
        lines.append(
            f"ATOM C CA {_AA3[amino_acid]} C {index} {index + 5.0:.3f} {5.0 + offset:.3f} {(index % 3) + offset * 0.1:.3f} 60.0"
        )
    lines.extend(["#", f"# {signature_marker}"])
    return "\n".join(lines) + "\n"


def _request_name(candidate_id):
    if candidate_id.startswith("HUMAN_BACKGROUND_"):
        return f"ebvms_bg_{candidate_id}_s04"
    return candidate_id


def _write_fixture_job(project_root, candidate_id, directory_name, signature_marker=""):
    job_dir = project_root / "jobs" / directory_name
    job_dir.mkdir(parents=True)
    peptide = _TEST_CANDIDATES[candidate_id][0]
    seed = str(sorted(_TEST_CANDIDATES).index(candidate_id) + 1)
    request = [{
        "name": _request_name(candidate_id),
        "modelSeeds": [seed],
        "sequences": [
            {"proteinChain": {"sequence": _HLA_SEQUENCE, "count": 1}},
            {"proteinChain": {"sequence": _HLA_SEQUENCE, "count": 1}},
            {"proteinChain": {"sequence": peptide, "count": 1}},
        ],
    }]
    prefix = f"fold_{directory_name}"
    (job_dir / f"{prefix}_job_request.json").write_text(json.dumps(request), encoding="utf-8")
    candidate_offset = sorted(_TEST_CANDIDATES).index(candidate_id) * 0.7
    for model_index in range(5):
        marker = signature_marker if model_index == 0 else ""
        (job_dir / f"{prefix}_model_{model_index}.cif").write_text(
            _fixture_cif(peptide, candidate_offset, model_index, marker), encoding="utf-8"
        )
        (job_dir / f"{prefix}_full_data_{model_index}.json").write_text("{}", encoding="utf-8")
        (job_dir / f"{prefix}_summary_confidences_{model_index}.json").write_text(
            json.dumps({"iptm": 0.8 + model_index * 0.01}), encoding="utf-8"
        )
    return job_dir


def _model_path(job_dir, index):
    return next(job_dir.glob(f"*_model_{index}.cif"))


def _build_live_project_fixture(project_root, rank2_stratum="32", rank2_target_peptide=None):
    primary_jobs = {}
    for candidate_id in _TEST_CANDIDATES:
        primary_jobs[candidate_id] = _write_fixture_job(
            project_root, candidate_id, f"00_{candidate_id.lower()}"
        )
    _write_fixture_job(project_root, "EBV_TCELL_950", "01_ebv_tcell_950_duplicate")
    _write_fixture_job(
        project_root, "EBV_TCELL_950", "02_ebv_tcell_950_distinct", signature_marker="distinct-signature"
    )

    rank1_pair = "EBV_TCELL_950::HUMAN_MYELIN_112214"
    rank2_pair = "EBV_TCELL_2268741::HUMAN_MYELIN_117032"
    rank1_controls = ("HUMAN_BACKGROUND_115891", "HUMAN_BACKGROUND_118550", "HUMAN_BACKGROUND_119732")
    rank2_controls = ("HUMAN_BACKGROUND_141561", "HUMAN_BACKGROUND_423369", "HUMAN_BACKGROUND_2258889")
    score_rows = [
        {
            "pair_id": rank1_pair,
            "discovery_priority_rank": "1",
            "ebv_peptide": _TEST_CANDIDATES["EBV_TCELL_950"][0],
            "ebv_p1_p9_core": _TEST_CANDIDATES["EBV_TCELL_950"][1],
            "human_peptide": _TEST_CANDIDATES["HUMAN_MYELIN_112214"][0],
            "human_p1_p9_core": _TEST_CANDIDATES["HUMAN_MYELIN_112214"][1],
        },
        {
            "pair_id": rank2_pair,
            "discovery_priority_rank": "2",
            "ebv_peptide": _TEST_CANDIDATES["EBV_TCELL_2268741"][0],
            "ebv_p1_p9_core": _TEST_CANDIDATES["EBV_TCELL_2268741"][1],
            "human_peptide": rank2_target_peptide or _TEST_CANDIDATES["HUMAN_MYELIN_117032"][0],
            "human_p1_p9_core": _TEST_CANDIDATES["HUMAN_MYELIN_117032"][1],
        },
    ]
    _write_fixture_csv(
        project_root / "processed/complete_model_pipeline_audit_2026-08-15/matched_background_structure_geometry.csv",
        [
            {"pair_id": rank1_pair, "background_candidate_id": candidate_id, "background_predicted_core": _TEST_CANDIDATES[candidate_id][1]}
            for candidate_id in rank1_controls
        ],
    )
    _write_fixture_csv(
        project_root / "processed/complete_model_pipeline_audit_2026-08-15/combined_same_register_geometry.csv",
        [{"pair_id": rank1_pair}, {"pair_id": rank2_pair}],
    )
    legacy_candidates = {
        "EBV_TCELL_950", "HUMAN_MYELIN_112214", *rank1_controls,
        "EBV_TCELL_2268741", "HUMAN_MYELIN_117032",
    }
    _write_fixture_csv(
        project_root / "processed/complete_model_pipeline_audit_2026-08-15/canonical_af3_job_summary.csv",
        [
            {"candidate_id": candidate_id, "source_path": str(primary_jobs[candidate_id])}
            for candidate_id in sorted(legacy_candidates)
        ],
    )
    _write_fixture_csv(
        project_root / "processed/complete_model_pipeline_audit_2026-08-15/canonical_af3_sample_metrics.csv",
        [
            {
                "candidate_id": candidate_id,
                "model_path": str(_model_path(primary_jobs[candidate_id], model_index)),
                "canonical_job_key": f"{candidate_id}|fixture",
                "sample_index": model_index,
            }
            for candidate_id in sorted(legacy_candidates)
            for model_index in range(5)
        ],
    )
    _write_fixture_csv(
        project_root / "processed/structural_control_expansion_2026-08-15/master_pair_score_sheet_with_expanded_controls.csv",
        score_rows,
    )
    _write_fixture_csv(
        project_root / "processed/structural_control_expansion_2026-08-15/complete_layered_control_geometry.csv",
        [
            {
                "pair_id": rank2_pair,
                "analysis_layer": "length_sensitivity_exact_bin_pm7",
                "background_candidate_id": candidate_id,
                "background_predicted_core": _TEST_CANDIDATES[candidate_id][1],
                "stratum_length": rank2_stratum,
            }
            for candidate_id in rank2_controls
        ],
    )
    _write_fixture_csv(
        project_root / "processed/structural_control_expansion_2026-08-15/alphafold_control_sample_metrics.csv",
        [
            {
                "candidate_id": candidate_id,
                "model_path": str(_model_path(primary_jobs[candidate_id], model_index).relative_to(project_root)),
                "canonical_job_key": f"{candidate_id}|fixture",
                "sample_index": model_index,
            }
            for candidate_id in rank2_controls
            for model_index in range(5)
        ],
    )
    for relative_path, controls in (
        ("processed/expanded_background/background_register_prediction_summary.csv", rank1_controls),
        ("processed/structural_control_expansion_2026-08-15/control_binding_prediction_summary.csv", rank2_controls),
    ):
        _write_fixture_csv(
            project_root / relative_path,
            [
                {
                    "candidate_id": candidate_id,
                    "prediction_status": "predicted",
                    "peptide": _TEST_CANDIDATES[candidate_id][0],
                    "predicted_core_peptide": _TEST_CANDIDATES[candidate_id][1],
                    "predicted_core_start_positions_1_based": _TEST_CANDIDATES[candidate_id][0].index(_TEST_CANDIDATES[candidate_id][1]) + 1,
                }
                for candidate_id in controls
            ],
        )
    return primary_jobs


def _lead_fixture(rank):
    if rank == 1:
        pair_id = "EBV_TCELL_950::HUMAN_MYELIN_112214"
        layer = "strict_primary_controls"
        ebv_id = "EBV_TCELL_950"
        target_id = "HUMAN_MYELIN_112214"
        control_ids = [
            "HUMAN_BACKGROUND_119732",
            "HUMAN_BACKGROUND_118550",
            "HUMAN_BACKGROUND_115891",
        ]
        target_distance, control_distances = 1.0, [5.0, 4.0, 3.0]
    else:
        pair_id = "EBV_TCELL_2268741::HUMAN_MYELIN_117032"
        layer = "length_sensitivity_exact_bin_pm7"
        ebv_id = "EBV_TCELL_2268741"
        target_id = "HUMAN_MYELIN_117032"
        control_ids = [
            "HUMAN_BACKGROUND_2258889",
            "HUMAN_BACKGROUND_423369",
            "HUMAN_BACKGROUND_141561",
        ]
        target_distance, control_distances = 2.0, [8.0, 7.0, 6.0]

    def job(candidate_id, suffix):
        return {
            "candidate_id": candidate_id,
            "job_id": f"{candidate_id}::{suffix}",
            "model_ids": [f"{candidate_id}::{suffix}::model-{index}" for index in range(5)],
        }

    ebv_job = job(ebv_id, "job")
    target_job = job(target_id, "job")
    control_jobs = {identifier: [job(identifier, "job")] for identifier in control_ids}
    entities = {"ebv": [ebv_job], "target": [target_job], "controls": control_jobs}
    groups = [("ebv", ebv_job["model_ids"]), ("target", target_job["model_ids"])] + [
        (identifier, control_jobs[identifier][0]["model_ids"]) for identifier in control_ids
    ]
    group_distance = {"target": target_distance, **dict(zip(control_ids, control_distances))}
    geometry = {}
    all_models = [(group, model) for group, models in groups for model in models]
    for left_index, (left_group, left_model) in enumerate(all_models):
        for right_group, right_model in all_models[left_index + 1:]:
            if left_group == right_group:
                value = 0.5
            elif left_group == "ebv":
                value = group_distance[right_group]
            elif right_group == "ebv":
                value = group_distance[left_group]
            else:
                value = 9.0
            geometry[(left_model, right_model)] = value
    confidence = {
        model: {"model_plddt": float(index + 10), "model_iptm": float(index) / 10}
        for index, (_, model) in enumerate(all_models)
    }
    return {
        "rank": rank,
        "pair_id": pair_id,
        "analysis_layer": layer,
        "ebv_candidate_id": ebv_id,
        "target_candidate_id": target_id,
        "control_ids": control_ids,
        "entities": entities,
        "geometry_lookup": geometry,
        "confidence_by_model": confidence,
        **(
            {"target_peptide_length": 32, "stratum_length": 32}
            if rank == 2 else {}
        ),
    }


class LeadDefinitionIntegrationTests(unittest.TestCase):
    def test_rank_two_rejects_any_layer_other_than_the_frozen_length_sensitivity_layer(self):
        lead = {
            "rank": 2,
            "pair_id": "EBV_TCELL_2268741::HUMAN_MYELIN_117032",
            "analysis_layer": "primary_exact_bin_length_pm1",
            "ebv_candidate_id": "EBV_TCELL_2268741",
            "target_candidate_id": "HUMAN_MYELIN_117032",
            "control_ids": [
                "HUMAN_BACKGROUND_141561",
                "HUMAN_BACKGROUND_423369",
                "HUMAN_BACKGROUND_2258889",
            ],
        }

        with self.assertRaisesRegex(ValueError, "length_sensitivity_exact_bin_pm7"):
            validate_lead_definition(lead)

    def test_rank_two_requires_an_exact_32_aa_target_and_frozen_32_aa_stratum(self):
        lead = {
            "rank": 2,
            "pair_id": "EBV_TCELL_2268741::HUMAN_MYELIN_117032",
            "analysis_layer": "length_sensitivity_exact_bin_pm7",
            "ebv_candidate_id": "EBV_TCELL_2268741",
            "target_candidate_id": "HUMAN_MYELIN_117032",
            "control_ids": [
                "HUMAN_BACKGROUND_141561",
                "HUMAN_BACKGROUND_423369",
                "HUMAN_BACKGROUND_2258889",
            ],
            "target_peptide_length": 32,
            "stratum_length": 32,
        }

        for field, wrong_value in (("target_peptide_length", 31), ("stratum_length", 31)):
            with self.subTest(field=field):
                invalid = {**lead, field: wrong_value}
                with self.assertRaisesRegex(ValueError, "exactly 32"):
                    validate_lead_definition(invalid)

        self.assertEqual(validate_lead_definition(lead)["analysis_layer"], "length_sensitivity_exact_bin_pm7")


class ConfigurableBootstrapIntegrationTests(unittest.TestCase):
    def test_injected_iteration_count_controls_full_replicate_rows(self):
        entities = {
            "ebv": [{"job_id": "ebv-job", "model_ids": [f"e{i}" for i in range(5)]}],
            "target": [{"job_id": "target-job", "model_ids": [f"t{i}" for i in range(5)]}],
            "controls": {
                "control": [{"job_id": "control-job", "model_ids": [f"c{i}" for i in range(5)]}]
            },
        }
        geometry = {}
        for ebv_model in entities["ebv"][0]["model_ids"]:
            for target_model in entities["target"][0]["model_ids"]:
                geometry[(ebv_model, target_model)] = 1.0
            for control_model in entities["controls"]["control"][0]["model_ids"]:
                geometry[(ebv_model, control_model)] = 4.0

        result = hierarchical_technical_bootstrap(entities, geometry, iterations=3, seed=7)

        self.assertEqual(len(result["replicates"]), 3)
        self.assertEqual(
            result["replicate_rows"],
            [
                {"iteration": 1, "target_median_A": 1.0, "equal_weight_background_median_A": 4.0, "delta_A": 3.0},
                {"iteration": 2, "target_median_A": 1.0, "equal_weight_background_median_A": 4.0, "delta_A": 3.0},
                {"iteration": 3, "target_median_A": 1.0, "equal_weight_background_median_A": 4.0, "delta_A": 3.0},
            ],
        )


class AuditTableIntegrationTests(unittest.TestCase):
    def test_summary_collapses_ebv_pair_medians_to_one_value_per_unique_target_job(self):
        lead = _lead_fixture(1)

        def job(candidate_id, suffix):
            return {
                "candidate_id": candidate_id,
                "job_id": f"{candidate_id}::{suffix}",
                "model_ids": [f"{candidate_id}::{suffix}::model-{index}" for index in range(5)],
            }

        ebv_jobs = [job(lead["ebv_candidate_id"], "ebv-1"), job(lead["ebv_candidate_id"], "ebv-2")]
        target_jobs = [job(lead["target_candidate_id"], "target-1"), job(lead["target_candidate_id"], "target-2")]
        lead["entities"]["ebv"] = ebv_jobs
        lead["entities"]["target"] = target_jobs
        grouped_jobs = [("ebv-1", ebv_jobs[0]), ("ebv-2", ebv_jobs[1]), ("target-1", target_jobs[0]), ("target-2", target_jobs[1])]
        grouped_jobs.extend(
            (control_id, jobs[0]) for control_id, jobs in lead["entities"]["controls"].items()
        )
        target_pair_values = {
            ("ebv-1", "target-1"): 1.0,
            ("ebv-2", "target-1"): 100.0,
            ("ebv-1", "target-2"): 2.0,
            ("ebv-2", "target-2"): 3.0,
        }
        geometry = {}
        flattened = [
            (group, model)
            for group, grouped_job in grouped_jobs
            for model in grouped_job["model_ids"]
        ]
        for left_index, (left_group, left_model) in enumerate(flattened):
            for right_group, right_model in flattened[left_index + 1:]:
                pair = (left_group, right_group)
                reverse = (right_group, left_group)
                if left_group == right_group:
                    value = 0.5
                elif pair in target_pair_values:
                    value = target_pair_values[pair]
                elif reverse in target_pair_values:
                    value = target_pair_values[reverse]
                elif "ebv-1" in pair or "ebv-2" in pair:
                    value = 60.0
                else:
                    value = 9.0
                geometry[(left_model, right_model)] = value
        lead["geometry_lookup"] = geometry
        lead["confidence_by_model"] = {
            model: {"model_plddt": 50.0, "model_iptm": 0.8}
            for _, model in flattened
        }

        tables = build_audit_tables([lead], iterations=1)

        summary = tables["control_rank_and_leave_one_out"][0]
        self.assertEqual(summary["overall_target_median_A"], 26.5)
        self.assertNotEqual(summary["overall_target_median_A"], 2.5)

    def test_two_small_injected_leads_emit_exactly_twenty_thousand_bootstrap_rows(self):
        tables = build_audit_tables([_lead_fixture(2), _lead_fixture(1)], iterations=10_000)

        rows = tables["technical_bootstrap_replicates"]
        self.assertEqual(len(rows), 20_000)
        self.assertEqual([row["lead_rank"] for row in rows[:10_000]], [1] * 10_000)
        self.assertEqual([row["iteration"] for row in rows[:10_000]], list(range(1, 10_001)))
        self.assertEqual([row["lead_rank"] for row in rows[10_000:]], [2] * 10_000)

    def test_every_output_table_has_deterministic_row_order(self):
        rank_one, rank_two = _lead_fixture(1), _lead_fixture(2)

        forward = build_audit_tables([rank_one, rank_two], iterations=3)
        reverse = build_audit_tables([rank_two, rank_one], iterations=3)

        self.assertEqual(forward, reverse)
        self.assertEqual(
            [row["control_candidate_id"] for row in forward["per_control_geometry_summary"][:3]],
            ["HUMAN_BACKGROUND_115891", "HUMAN_BACKGROUND_118550", "HUMAN_BACKGROUND_119732"],
        )

    def test_confidence_changes_annotations_but_cannot_change_integrated_pose_clusters(self):
        lead = _lead_fixture(1)
        baseline = build_audit_tables([lead], iterations=2)["pose_cluster_membership"]
        for values in lead["confidence_by_model"].values():
            values["model_plddt"] = 999.0 - values["model_plddt"]
            values["model_iptm"] = 1.0 - values["model_iptm"]
        changed = build_audit_tables([lead], iterations=2)["pose_cluster_membership"]

        cluster_fields = ("model_id", "cluster_id", "cluster_size", "distinct_jobs_in_cluster")
        self.assertEqual(
            [tuple(row[field] for field in cluster_fields) for row in baseline],
            [tuple(row[field] for field in cluster_fields) for row in changed],
        )
        self.assertNotEqual(
            [row["model_plddt"] for row in baseline],
            [row["model_plddt"] for row in changed],
        )


class FigureIntegrationTests(unittest.TestCase):
    def test_rank_two_uses_compact_visible_job_labels_and_accessible_full_ids(self):
        summary = {
            "overall_target_median_A": 1.0,
            "overall_background_median_A": 4.0,
            "background_minus_target_delta_A": 3.0,
            "leave_one_out_delta_min_A": 2.0,
            "leave_one_out_delta_max_A": 4.0,
            "target_rank": 1,
            "control_count": 3,
            "classification": "length_sensitivity_only__mixed_positive",
        }
        controls = [
            {"control_candidate_id": f"HUMAN_BACKGROUND_{identifier}", "control_median_A": value}
            for identifier, value in ((141561, 3.0), (423369, 4.0), (2258889, 5.0))
        ]
        bootstrap = {"delta_percentile_2_5_A": -1.5, "delta_percentile_97_5_A": 4.5}
        full_ids = [
            "HUMAN_MYELIN_117032|seed-179807928|sig-e2e8bb38a67a3c46",
            "HUMAN_MYELIN_117032|seed-27|sig-2f9a0e012ee11aa4",
        ]

        svg = render_rank2_svg(
            summary,
            controls,
            bootstrap,
            target_job_rows=[
                {"comparison_job_id": full_ids[0], "exposed_rmsd_median_A": 1.0},
                {"comparison_job_id": full_ids[1], "exposed_rmsd_median_A": 2.0},
            ],
        )

        root = ET.fromstring(svg)
        visible_labels = [
            "".join(element.itertext())
            for element in root.iter()
            if element.tag.endswith("text") and element.attrib.get("class") == "label"
        ]
        accessible_titles = [
            "".join(element.itertext())
            for element in root.iter()
            if element.tag.endswith("title")
        ]
        self.assertIn("Target human job 1", visible_labels)
        self.assertIn("Target human job 2", visible_labels)
        self.assertTrue(all(identifier not in visible_labels for identifier in full_ids))
        self.assertTrue(all(identifier in accessible_titles for identifier in full_ids))

    def test_rank_figures_escape_dynamic_labels_and_explain_their_intervals(self):
        summary = {
            "overall_target_median_A": 1.0,
            "overall_background_median_A": 4.0,
            "background_minus_target_delta_A": 3.0,
            "leave_one_out_delta_min_A": 2.0,
            "leave_one_out_delta_max_A": 4.0,
            "target_rank": 1,
            "control_count": 3,
            "classification": "consistent_positive",
        }
        controls = [
            {"control_candidate_id": "control&A", "control_median_A": 3.0},
            {"control_candidate_id": "control<B>", "control_median_A": 4.0},
            {"control_candidate_id": "control-C", "control_median_A": 5.0},
        ]
        bootstrap = {"delta_percentile_2_5_A": 1.5, "delta_percentile_97_5_A": 4.5}

        primary = render_rank1_svg(summary, controls, bootstrap, title="Rank 1 <primary> & strict")
        sensitivity = render_rank2_svg(
            summary,
            controls,
            bootstrap,
            target_job_rows=[
                {"comparison_job_id": "human&job-1", "exposed_rmsd_median_A": 1.0},
                {"comparison_job_id": "human<job-2>", "exposed_rmsd_median_A": 2.0},
            ],
        )

        self.assertIn("Rank 1 &lt;primary&gt; &amp; strict", primary)
        self.assertIn("control&amp;A", primary)
        self.assertIn("control&lt;B&gt;", primary)
        self.assertNotIn("control<B>", primary)
        self.assertIn("Target rank: 1 of 4", primary)
        self.assertIn("Leave-one-control-out delta range", primary)
        self.assertIn("Technical-stability interval (not a p-value)", primary)
        self.assertIn("Supplemental / length-sensitivity-only", sensitivity)
        self.assertIn("human&amp;job-1", sensitivity)
        self.assertIn("human&lt;job-2&gt;", sensitivity)


class LiveProjectPathIntegrationTests(unittest.TestCase):
    def test_temp_project_exercises_dedup_geometry_manifest_checksums_and_generator(self):
        with tempfile.TemporaryDirectory() as directory:
            project_root = Path(directory)
            _build_live_project_fixture(project_root)

            metadata = load_fixed_lead_metadata(project_root)
            self.assertEqual(metadata["lead_specs"][1]["target_peptide_length"], 32)
            self.assertEqual(metadata["lead_specs"][1]["stratum_length"], 32)

            jobs, consumed_hashes = discover_saved_jobs(
                project_root, metadata["candidate_metadata"]
            )
            ebv_950_jobs = [job for job in jobs if job["candidate_id"] == "EBV_TCELL_950"]
            self.assertEqual(len(jobs), 11)
            self.assertEqual(len(ebv_950_jobs), 2)
            self.assertEqual(sum(len(job["duplicate_paths"]) for job in ebv_950_jobs), 1)
            self.assertEqual(len({job["request_identity"] for job in ebv_950_jobs}), 1)
            self.assertEqual({job["request_name"] for job in ebv_950_jobs}, {"EBV_TCELL_950"})
            self.assertEqual({job["seed"] for job in ebv_950_jobs}, {"2"})
            self.assertEqual(
                {tuple(job["chain_sequences"]) for job in ebv_950_jobs},
                {(_HLA_SEQUENCE, _HLA_SEQUENCE, _TEST_CANDIDATES["EBV_TCELL_950"][0])},
            )
            self.assertEqual(len({job["complete_model_signature"] for job in ebv_950_jobs}), 2)
            self.assertTrue(all(job["sequence_integrity_status"].startswith("pass_exact") for job in jobs))

            lead_inputs = build_live_lead_inputs(
                metadata["lead_specs"], metadata["candidate_metadata"], jobs
            )
            self.assertEqual([len(lead["geometry_lookup"]) for lead in lead_inputs], [435, 300])

            manifest = build_identity_manifest(metadata["lead_specs"], jobs)
            self.assertEqual(len(manifest), 11)
            self.assertEqual(sum(int(row["duplicate_path_count"]) for row in manifest), 1)
            self.assertEqual(
                len({row["complete_model_signature"] for row in manifest if row["candidate_id"] == "EBV_TCELL_950"}),
                2,
            )

            checksum_rows = collect_input_checksums(
                project_root, metadata["input_paths"], consumed_hashes
            )
            self.assertEqual(len(checksum_rows), 201)
            for row in checksum_rows:
                path = project_root / row["input_path"]
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row["sha256"])

            output = generate_live_audit(project_root)
            self.assertEqual(
                sorted(path.name for path in output.iterdir()),
                sorted([
                    "model_job_identity_manifest.csv",
                    "per_control_geometry_summary.csv",
                    "job_pair_stability.csv",
                    "control_rank_and_leave_one_out.csv",
                    "technical_bootstrap_replicates.csv",
                    "technical_bootstrap_summary.csv",
                    "pose_cluster_membership.csv",
                    "rank1_primary_control_robustness.svg",
                    "rank2_length_sensitivity_job_dependence.svg",
                    "LEAD_FOCUSED_FINDINGS.md",
                    "frozen_input_checksums.csv",
                ]),
            )
            with (output / "technical_bootstrap_replicates.csv").open(newline="", encoding="utf-8") as handle:
                bootstrap_rows = list(csv.DictReader(handle))
            self.assertEqual(len(bootstrap_rows), 20_000)
            with (output / "technical_bootstrap_summary.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                bootstrap_summary_rows = list(csv.DictReader(handle))
            self.assertEqual(
                {row["iterations"] for row in bootstrap_summary_rows}, {"10000"}
            )
            with (output / "control_rank_and_leave_one_out.csv").open(newline="", encoding="utf-8") as handle:
                rank_rows = list(csv.DictReader(handle))
            self.assertTrue(rank_rows[1]["classification"].startswith("length_sensitivity_only__"))
            with self.assertRaises(FileExistsError):
                generate_live_audit(project_root)
            with self.assertRaises(TypeError):
                generate_live_audit(project_root, iterations=3)

    def test_metadata_rejects_non32_rank2_target_stratum_ambiguous_core_and_wrong_layer(self):
        mutations = (
            ("stratum", "Rank 2 sensitivity output requires an exactly 32-aa frozen stratum"),
            ("target", "Rank 2 sensitivity output requires an exactly 32-aa target peptide"),
            ("core", "ambiguous"),
            ("layer", "length_sensitivity_exact_bin_pm7"),
        )
        for mutation, message in mutations:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                project_root = Path(directory)
                _build_live_project_fixture(
                    project_root,
                    rank2_stratum="31" if mutation == "stratum" else "32",
                    rank2_target_peptide="RSTVWYACDRSTVWYACD" if mutation == "target" else None,
                )
                if mutation == "core":
                    prediction_path = project_root / "processed/structural_control_expansion_2026-08-15/control_binding_prediction_summary.csv"
                    with prediction_path.open(newline="", encoding="utf-8") as handle:
                        rows = list(csv.DictReader(handle))
                    rows[0]["peptide"] = rows[0]["predicted_core_peptide"] * 2
                    rows[0]["predicted_core_start_positions_1_based"] = "1;10"
                    _write_fixture_csv(prediction_path, rows)
                if mutation == "layer":
                    geometry_path = project_root / "processed/structural_control_expansion_2026-08-15/complete_layered_control_geometry.csv"
                    with geometry_path.open(newline="", encoding="utf-8") as handle:
                        rows = list(csv.DictReader(handle))
                    for row in rows:
                        row["analysis_layer"] = "primary_exact_bin_length_pm1"
                    _write_fixture_csv(geometry_path, rows)
                with self.assertRaisesRegex(ValueError, message):
                    load_fixed_lead_metadata(project_root)

    def test_discovery_rejects_wrong_request_layout_peptide_and_incomplete_file_set(self):
        for mutation, message in (
            ("two_chains", "three-chain"),
            ("wrong_peptide", "Requested peptide mismatch"),
            ("missing_confidence", "exactly one confidence file"),
        ):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                project_root = Path(directory)
                primary_jobs = _build_live_project_fixture(project_root)
                metadata = load_fixed_lead_metadata(project_root)
                job_dir = primary_jobs["EBV_TCELL_950"]
                request_path = next(job_dir.glob("*_job_request.json"))
                if mutation in {"two_chains", "wrong_peptide"}:
                    request = json.loads(request_path.read_text(encoding="utf-8"))
                    if mutation == "two_chains":
                        request[0]["sequences"] = request[0]["sequences"][:2]
                    else:
                        request[0]["sequences"][2]["proteinChain"]["sequence"] = "ACDEFGHIK"
                    request_path.write_text(json.dumps(request), encoding="utf-8")
                else:
                    next(job_dir.glob("*_summary_confidences_4.json")).unlink()
                with self.assertRaisesRegex(ValueError, message):
                    discover_saved_jobs(project_root, metadata["candidate_metadata"])


class AuditWriterIntegrationTests(unittest.TestCase):
    def test_writer_creates_only_the_eleven_required_files_and_refuses_overwrite(self):
        tables = {
            "per_control_geometry_summary": [{"lead_rank": 1}],
            "job_pair_stability": [{"lead_rank": 1}],
            "control_rank_and_leave_one_out": [{"lead_rank": 1}],
            "technical_bootstrap_replicates": [{"lead_rank": 1}],
            "technical_bootstrap_summary": [{"lead_rank": 1}],
            "pose_cluster_membership": [{"lead_rank": 1}],
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "processed" / "lead_focused_robustness_2026-08-15"

            write_audit_outputs(
                output,
                manifest_rows=[{"lead_rank": 1}],
                tables=tables,
                rank1_svg="<svg/>",
                rank2_svg="<svg/>",
                findings="# Findings\n",
                checksum_rows=[{"input_path": "input.csv", "sha256": "abc"}],
            )

            self.assertEqual(
                sorted(path.name for path in output.iterdir()),
                sorted([
                    "model_job_identity_manifest.csv",
                    "per_control_geometry_summary.csv",
                    "job_pair_stability.csv",
                    "control_rank_and_leave_one_out.csv",
                    "technical_bootstrap_replicates.csv",
                    "technical_bootstrap_summary.csv",
                    "pose_cluster_membership.csv",
                    "rank1_primary_control_robustness.svg",
                    "rank2_length_sensitivity_job_dependence.svg",
                    "LEAD_FOCUSED_FINDINGS.md",
                    "frozen_input_checksums.csv",
                ]),
            )
            with self.assertRaises(FileExistsError):
                write_audit_outputs(
                    output,
                    manifest_rows=[{"lead_rank": 1}],
                    tables=tables,
                    rank1_svg="<svg/>",
                    rank2_svg="<svg/>",
                    findings="# Findings\n",
                    checksum_rows=[{"input_path": "input.csv", "sha256": "abc"}],
                )


if __name__ == "__main__":
    unittest.main()
