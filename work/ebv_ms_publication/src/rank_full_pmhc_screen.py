"""Rank the full Tier-1 EBV x myelin pMHC screen for structural review.

This is a transparent shortlist heuristic, not a cross-reactivity predictor.
It combines provenance-preserving sequence descriptors with peptide-level
model confidence and a pMHC-groove-normalized local pose metric.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"


def read_csv(name: str) -> list[dict[str, str]]:
    with (PROC / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    manifest = read_csv("pmhc_candidate_manifest.csv")
    geometry = read_csv("colabfold_tier1_ebv_myelin_geometry_matrix.csv")
    similarity = read_csv("ebv_physicochemical_pairwise.csv")

    # An EBV peptide can appear in both the Tier-1 and MHC-ligand arms.  Keep
    # the Tier-1 identity explicitly rather than silently overwriting it.
    tier1_id_by_peptide = {
        row["peptide"]: row["candidate_id"]
        for row in manifest
        if row["candidate_id"].startswith("EBV_TCELL_")
    }
    myelin_ids_by_peptide: dict[str, list[str]] = {}
    for row in manifest:
        if row["candidate_id"].startswith("HUMAN_MYELIN_"):
            myelin_ids_by_peptide.setdefault(row["peptide"], []).append(row["candidate_id"])

    descriptors: dict[tuple[str, str], dict[str, str]] = {}
    for row in similarity:
        if row["target_group"] != "myelin":
            continue
        ebv_id = tier1_id_by_peptide.get(row["ebv_peptide"])
        for human_id in myelin_ids_by_peptide.get(row["human_peptide"], []):
            if ebv_id:
                descriptors[(ebv_id, human_id)] = row

    rows = []
    for row in geometry:
        descriptor = descriptors.get((row["ebv_candidate_id"], row["human_candidate_id"]))
        if not descriptor or row["status"] != "PASS":
            continue
        aligned = int(descriptor["aligned_length"])
        if aligned < 5:
            continue
        ebv_confidence = float(row["ebv_peptide_mean_plddt"])
        human_confidence = float(row["human_peptide_mean_plddt"])
        local_rmsd = float(row["local_peptide_ca_rmsd_after_hla_fit"])
        property_similarity = float(descriptor["property_similarity"])
        # Penalize low-confidence peptide models and divergent poses.  The
        # score orders candidates for review only; it is never effect size.
        review_priority = property_similarity * min(ebv_confidence, human_confidence) / 100 * math.exp(-local_rmsd / 5)
        rows.append({
            "ebv_candidate_id": row["ebv_candidate_id"],
            "human_candidate_id": row["human_candidate_id"],
            "locally_aligned_residues": aligned,
            "property_similarity": property_similarity,
            "ebv_peptide_mean_plddt": ebv_confidence,
            "human_peptide_mean_plddt": human_confidence,
            "hla_groove_ca_rmsd_after_fit": row["hla_groove_ca_rmsd_after_fit"],
            "local_peptide_ca_rmsd_after_hla_fit": local_rmsd,
            "review_priority_heuristic": round(review_priority, 6),
            "interpretation": "hypothesis-ranking metric only; not evidence of TCR cross-reactivity",
        })

    rows.sort(key=lambda row: row["review_priority_heuristic"], reverse=True)
    out = PROC / "fullscreen_tier1_ebv_myelin_shortlist.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out} ({len(rows)} pairs with >=5 locally aligned residues)")


if __name__ == "__main__":
    main()
