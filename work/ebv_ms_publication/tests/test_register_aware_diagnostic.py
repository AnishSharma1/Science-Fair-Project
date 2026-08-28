"""Tests for register-position diagnostic helpers.

The diagnostic is intentionally descriptive: it must preserve every possible
window and never silently convert an unassessable short/flanked peptide into a
favorable register claim.
"""

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from register_aware_diagnostic import (  # noqa: E402
    parse_local_alignment_positions,
    same_register_alignment_count,
    window_pair_sensitivity,
)


class RegisterAwareDiagnosticTests(unittest.TestCase):
    def test_parse_local_alignment_positions_retains_residue_coordinates(self):
        self.assertEqual(
            parse_local_alignment_positions("4Y:3H;7F:6F"),
            [(4, "Y", 3, "H"), (7, "F", 6, "F")],
        )

    def test_same_register_count_requires_identical_p_positions(self):
        alignment = [(4, "Y", 3, "H"), (7, "F", 6, "F")]
        self.assertEqual(
            same_register_alignment_count(alignment, ebv_core_start=4, human_core_start=3),
            2,
        )
        self.assertEqual(
            same_register_alignment_count(alignment, ebv_core_start=4, human_core_start=5),
            0,
        )

    def test_window_sensitivity_retains_all_manifest_core_combinations(self):
        # Two 10-mers each contain two possible 9-mer windows, so the output
        # must retain four pre-review hypotheses rather than one best row.
        rows = window_pair_sensitivity(
            ebv_peptide="ABCDEFGHIJ",
            human_peptide="KLMNOPQRST",
            alignment=[(1, "A", 1, "K")],
        )
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["same_register_alignment_count"] for row in rows}, {0, 1})
        self.assertEqual(
            sum(row["same_register_alignment_count"] == 1 for row in rows), 1
        )


if __name__ == "__main__":
    unittest.main()
