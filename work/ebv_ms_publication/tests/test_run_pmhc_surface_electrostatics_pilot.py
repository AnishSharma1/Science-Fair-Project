import tempfile
import unittest
from pathlib import Path

from run_pmhc_surface_electrostatics_pilot import (
    TARGET_IDS,
    build_lead_gate,
    parse_propka_histidines,
    validate_panel_rows,
)


class PilotWorkflowTests(unittest.TestCase):
    def test_exactly_two_frozen_leads_are_in_scope(self):
        self.assertEqual(TARGET_IDS, ("HY13_SEQ_02", "HY15_SEQ_02"))

    def test_panel_requires_one_target_and_twenty_five_n3_rows(self):
        rows = [
            {"row_role": "target" if index == 0 else "n3", "pair_id": f"P{index}", "allele": "HLA-DRB1*15:01"}
            for index in range(26)
        ]
        validate_panel_rows("PANEL", rows)
        with self.assertRaisesRegex(ValueError, "26 rows"):
            validate_panel_rows("PANEL", rows[:-1])
        with self.assertRaisesRegex(ValueError, "one target"):
            validate_panel_rows(
                "PANEL",
                [{"row_role": "n3", "pair_id": f"P{index}", "allele": "HLA-DRB1*15:01"} for index in range(26)],
            )

    def test_propka_parser_preserves_peptide_histidine_positions(self):
        text = """
   HIS   6 C     6.64       6.50
   HIS  11 C     2.69       6.50
   HIS  13 C     1.11       6.50
   HIS  80 B     6.04       6.50
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.log"
            path.write_text(text)
            rows = parse_propka_histidines(path)
        self.assertEqual([row["sequence_position_1_based"] for row in rows], [6, 11, 13])
        self.assertEqual([row["predicted_pka"] for row in rows], [6.64, 2.69, 1.11])

    def test_gate_never_unlocks_discovery_or_specificity(self):
        gate = build_lead_gate(
            target_id="HY15_SEQ_02",
            rank=1,
            register_qc=False,
            model_qc=True,
            dielectric_robust=True,
        )
        self.assertEqual(gate["status"], "not_evaluable")
        self.assertEqual(gate["rank_only_context"], "electrostatic_context_supportive")
        self.assertFalse(gate["weights_frozen"])
        self.assertFalse(gate["discovery_unlock_allowed"])
        self.assertFalse(gate["specificity_claim_allowed"])


if __name__ == "__main__":
    unittest.main()
