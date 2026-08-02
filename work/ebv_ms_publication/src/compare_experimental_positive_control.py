"""Compare the experimental DR2a--BALF5 and DR2b--MBP pMHC structures.

This is a provenance-anchored positive control, not a discovery analysis.
Inputs are RCSB PDB 1H15 (DRB5*01:01--BALF5) and 1BX2
(DRB1*15:01--MBP).  Each PDB has two crystallographic copies; chain triplets
A/B/C are used consistently.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from triage_colabfold_pmhc import ca_coordinates, kabsch, parse_pdb, sequence


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "experimental_positive_control"


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: compare_experimental_positive_control.py 1H15.pdb 1BX2.pdb")
    ebv = parse_pdb(Path(sys.argv[1]))
    mbp = parse_pdb(Path(sys.argv[2]))
    # The first 85 residues of each alpha/beta chain define the peptide-binding
    # platform.  This avoids allowing distal domains to dominate the fit.
    ebv_frame = np.vstack((ca_coordinates(ebv["A"][:85]), ca_coordinates(ebv["B"][:85])))
    mbp_frame = np.vstack((ca_coordinates(mbp["A"][:85]), ca_coordinates(mbp["B"][:85])))
    rotation, translation, groove_rmsd = kabsch(ebv_frame, mbp_frame)
    ebv_peptide = ca_coordinates(ebv["C"])
    mbp_peptide = ca_coordinates(mbp["C"])
    ebv_fitted = ebv_peptide @ rotation + translation

    # The published TCR-facing equivalent core is the H/F/V/K region.  Compare
    # those positional residues after fitting the experimentally determined HLA
    # grooves. This metric describes geometric similarity only.
    ebv_core = "YHFVKKH"
    mbp_core = "VHFFKNI"
    ebv_seq, mbp_seq = sequence(ebv["C"]), sequence(mbp["C"])
    ebv_start, mbp_start = ebv_seq.index(ebv_core), mbp_seq.index(mbp_core)
    pairs = [(ebv_start + i, mbp_start + i) for i in range(len(ebv_core))]
    # There are different residues in the homologous core, so retain all seven
    # aligned positions rather than selecting only exact matches.
    distances = [float(np.linalg.norm(ebv_fitted[i] - mbp_peptide[j])) for i, j in pairs]
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "experimental_core_position_distances.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "core_position", "balf5_residue", "mbp_residue", "residue_pair",
            "ca_distance_after_hla_fit_a"
        ])
        writer.writeheader()
        for position, ((i, j), value) in enumerate(zip(pairs, distances), start=1):
            writer.writerow({
                "core_position": position,
                "balf5_residue": ebv_seq[i],
                "mbp_residue": mbp_seq[j],
                "residue_pair": f"{ebv_seq[i]}/{mbp_seq[j]}",
                "ca_distance_after_hla_fit_a": f"{value:.4f}",
            })
    rows = [
        {"metric": "EBV_structure", "value": "1H15", "interpretation": "Experimental DRB5*01:01--BALF5 pMHC"},
        {"metric": "MBP_structure", "value": "1BX2", "interpretation": "Experimental DRB1*15:01--MBP pMHC"},
        {"metric": "EBV_modeled_peptide", "value": ebv_seq, "interpretation": "Resolved residues; N-terminal Thr is not modeled"},
        {"metric": "MBP_modeled_peptide", "value": mbp_seq, "interpretation": "Resolved residues"},
        {"metric": "HLA_groove_CA_RMSD_A", "value": f"{groove_rmsd:.3f}", "interpretation": "After fitting alpha1/beta1 peptide-binding platforms"},
        {"metric": "seven_position_core_CA_RMSD_A", "value": f"{np.sqrt(np.mean(np.square(distances))):.3f}", "interpretation": "After HLA-groove fit; geometric similarity, not TCR binding evidence"},
    ]
    with (OUT / "experimental_drb2_positive_control_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)
    (OUT / "README.md").write_text(
        "# Experimental positive control\n\n"
        "This comparison uses published experimental structures, not ColabFold. "
        "PDB 1H15 is HLA-DRA1*01:01/DRB5*01:01 with EBV BALF5(627-641); "
        "PDB 1BX2 is HLA-DRA1*01:01/DRB1*15:01 with MBP(85-99). "
        "They are the paired pMHC ligands for the Hy.2E11 cross-reactivity system. "
        "Metrics document pMHC surface geometry only and cannot establish TCR recognition.\n"
    )


if __name__ == "__main__":
    main()
