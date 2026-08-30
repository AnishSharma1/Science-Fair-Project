import csv
import tempfile
import unittest
from pathlib import Path

from run_pmhc_surface_electrostatics_pilot import (
    REMAINING_SEQUENCE_TARGET_IDS,
    TARGET_IDS,
    _apbs_command,
    _panel_rows,
    _package_target_ids,
    _patch_protocol_amendment,
    resolve_target_set,
)


class SequenceExpansionTargetTests(unittest.TestCase):
    def test_remaining_sequence_targets_are_the_six_untested_candidates(self):
        self.assertEqual(
            REMAINING_SEQUENCE_TARGET_IDS,
            (
                "HY03_SEQ_01",
                "HY03_SEQ_02",
                "HY08_SEQ_01",
                "HY08_SEQ_02",
                "HY13_SEQ_01",
                "HY15_SEQ_01",
            ),
        )
        self.assertTrue(set(REMAINING_SEQUENCE_TARGET_IDS).isdisjoint(TARGET_IDS))

    def test_target_set_resolution_keeps_pilot_and_expansion_separate(self):
        pilot = resolve_target_set("pilot")
        expansion = resolve_target_set("remaining-sequence")

        self.assertEqual(pilot["target_ids"], TARGET_IDS)
        self.assertEqual(pilot["protocol_version"], "pmhc_surface_electrostatics_pilot_v1")
        self.assertEqual(expansion["target_ids"], REMAINING_SEQUENCE_TARGET_IDS)
        self.assertEqual(
            expansion["protocol_version"],
            "pmhc_surface_electrostatics_sequence_expansion_v1",
        )
        self.assertNotEqual(pilot["output_dir"], expansion["output_dir"])

    def test_remaining_sequence_panels_are_complete_and_exact_hla(self):
        panels = _panel_rows(REMAINING_SEQUENCE_TARGET_IDS)

        self.assertEqual(tuple(panels), REMAINING_SEQUENCE_TARGET_IDS)
        for target_id, rows in panels.items():
            self.assertEqual(len(rows), 26, target_id)
            self.assertEqual(sum(row["row_role"] == "target" for row in rows), 1)
            self.assertEqual(sum(row["row_role"] == "n3" for row in rows), 25)
            self.assertEqual(len({row["allele"] for row in rows}), 1)

    def test_package_stages_recover_target_ids_from_the_frozen_registry(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with (output / "frozen_target_registry.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=("target_id",))
                writer.writeheader()
                for target_id in REMAINING_SEQUENCE_TARGET_IDS:
                    writer.writerow({"target_id": target_id})

            self.assertEqual(_package_target_ids(output), REMAINING_SEQUENCE_TARGET_IDS)

    def test_persistent_apbs_command_uses_the_same_pinned_container_runtime(self):
        command = _apbs_command(
            Path("/tmp/output"),
            Path("/tmp/output/raw_calculations/apbs/test.in"),
            container_name="pmhc-apbs-test",
        )

        self.assertEqual(command[:4], ["docker", "exec", "-e", "LD_LIBRARY_PATH=/apbs/lib"])
        self.assertEqual(command[4], "pmhc-apbs-test")
        self.assertEqual(command[-2:], ["/apbs/bin/apbs", "/work/raw_calculations/apbs/test.in"])
        self.assertNotIn("--platform", command)

    def test_expansion_records_inherited_geometry_rule_without_discarded_analysis(self):
        amendment = _patch_protocol_amendment(
            {"protocol_version": "pmhc_surface_electrostatics_sequence_expansion_v1", "frozen_at": "2026-08-30"}
        )

        self.assertEqual(amendment["date"], "2026-08-30")
        self.assertEqual(amendment["initial_analysis_status"], "not_run")
        self.assertFalse(amendment["uses_electrostatic_scores"])


if __name__ == "__main__":
    unittest.main()
