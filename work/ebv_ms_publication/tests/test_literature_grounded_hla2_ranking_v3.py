import math
import unittest

import numpy as np


from build_literature_grounded_hla2_rankings_v3 import (
    SurfaceModel,
    allowed_register_starts,
    assign_evidence_tier,
    enumerate_register_windows,
    rank_v3_rows,
    sidechain_centroid,
    surface_pair_metrics,
)


def surface_model(offset=0.0):
    groove = np.asarray(
        [[float(i), math.sin(i), math.cos(i)] for i in range(12)], dtype=float
    )
    ca = np.asarray([[float(i), 2.0, 0.2 * i] for i in range(9)], dtype=float)
    centroids = ca + np.asarray([[0.0, 1.0, 0.5]] * 9)
    orientations = centroids - ca
    exposure = np.asarray([1, 2, 4, 7, 9, 7, 4, 2, 1], dtype=float)
    return SurfaceModel(
        sequence="ACDEFGHIK",
        groove_ca=groove + offset,
        peptide_ca=ca + offset,
        sidechain_centroids=centroids + offset,
        sidechain_orientations=orientations,
        sidechain_sasa=exposure,
    )


class SurfaceGeometryTests(unittest.TestCase):
    def test_pair_metrics_are_rotation_and_translation_invariant(self):
        reference = surface_model()
        theta = math.pi / 3
        rotation = np.asarray(
            [
                [math.cos(theta), -math.sin(theta), 0.0],
                [math.sin(theta), math.cos(theta), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        translation = np.asarray([13.0, -7.0, 4.5])
        transformed = SurfaceModel(
            sequence=reference.sequence,
            groove_ca=reference.groove_ca @ rotation + translation,
            peptide_ca=reference.peptide_ca @ rotation + translation,
            sidechain_centroids=reference.sidechain_centroids @ rotation + translation,
            sidechain_orientations=reference.sidechain_orientations @ rotation,
            sidechain_sasa=reference.sidechain_sasa.copy(),
        )
        metrics = surface_pair_metrics(reference, transformed)
        for key in (
            "exposure_weighted_backbone_rmsd_A",
            "exposure_weighted_centroid_rmsd_A",
            "exposure_weighted_orientation_rmsd_A",
            "exposed_distance_matrix_rmsd_A",
        ):
            self.assertAlmostEqual(metrics[key], 0.0, places=8)

    def test_glycine_sidechain_falls_back_to_ca(self):
        residue = {
            "aa": "G",
            "atoms": [
                {"name": "N", "element": "N", "xyz": (0.0, 0.0, 0.0)},
                {"name": "CA", "element": "C", "xyz": (1.0, 2.0, 3.0)},
                {"name": "C", "element": "C", "xyz": (2.0, 2.0, 3.0)},
                {"name": "O", "element": "O", "xyz": (3.0, 2.0, 3.0)},
            ],
        }
        np.testing.assert_allclose(sidechain_centroid(residue), [1.0, 2.0, 3.0])


class RegisterTests(unittest.TestCase):
    def test_all_windows_and_local_allowed_starts_are_deterministic(self):
        windows = enumerate_register_windows("ABCDEFGHIJKLMNO")
        self.assertEqual(len(windows), 7)
        self.assertEqual(windows[0], (1, "ABCDEFGHI"))
        self.assertEqual(windows[-1], (7, "GHIJKLMNO"))
        self.assertEqual(allowed_register_starts(15, 4), (3, 4, 5))
        self.assertEqual(allowed_register_starts(9, 1), (1,))


class RankingAndTierTests(unittest.TestCase):
    def test_blosum_is_primary_and_structure_only_breaks_equal_sequence_scores(self):
        rows = [
            {
                "allele": "HLA-DRB1*15:01",
                "pair_id": "lower_blosum_better_structure",
                "tcr_facing_blosum62_similarity": 0.6,
                "tcr_face_physicochemical_mismatch": 0.01,
                "tcr_facing_sequence_identity": 1.0,
                "local_surface_percentile": 0.01,
            },
            {
                "allele": "HLA-DRB1*15:01",
                "pair_id": "higher_blosum",
                "tcr_facing_blosum62_similarity": 0.7,
                "tcr_face_physicochemical_mismatch": 0.50,
                "tcr_facing_sequence_identity": 0.0,
                "local_surface_percentile": 0.90,
            },
            {
                "allele": "HLA-DRB1*15:01",
                "pair_id": "tie_worse_structure",
                "tcr_facing_blosum62_similarity": 0.6,
                "tcr_face_physicochemical_mismatch": 0.01,
                "tcr_facing_sequence_identity": 1.0,
                "local_surface_percentile": 0.10,
            },
            {
                "allele": "HLA-DQB1*06:02",
                "pair_id": "different_hla",
                "tcr_facing_blosum62_similarity": 0.1,
                "tcr_face_physicochemical_mismatch": 1.0,
                "tcr_facing_sequence_identity": 0.0,
                "local_surface_percentile": 1.0,
            },
        ]
        ranked = rank_v3_rows(rows)
        dr = [row for row in ranked if row["allele"] == "HLA-DRB1*15:01"]
        self.assertEqual(
            [row["pair_id"] for row in dr],
            ["higher_blosum", "lower_blosum_better_structure", "tie_worse_structure"],
        )
        self.assertEqual([row["primary_rank"] for row in dr], [1, 2, 3])
        dq = next(row for row in ranked if row["allele"] == "HLA-DQB1*06:02")
        self.assertEqual(dq["primary_rank"], 1)

    def test_missing_or_register_uncertain_abstains_with_m_tier(self):
        self.assertEqual(assign_evidence_tier(0.05, None, True), "M")
        self.assertEqual(assign_evidence_tier(0.05, 0.05, False), "M")
        self.assertEqual(assign_evidence_tier(0.05, 0.05, True), "A")
        self.assertEqual(assign_evidence_tier(0.20, 0.20, True), "B")
        self.assertEqual(assign_evidence_tier(0.20, 0.80, True), "C")
        self.assertEqual(assign_evidence_tier(0.80, 0.20, True), "D")
        self.assertEqual(assign_evidence_tier(0.80, 0.80, True), "E")


if __name__ == "__main__":
    unittest.main()
