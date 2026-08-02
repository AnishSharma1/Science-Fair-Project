"""Scaffold a completed pMHC peptide onto experimental DR15--Ob.1A12 geometry.

This is a geometry-compatibility check.  The target pMHC's DRB1 chain is fitted
to the DRB1 chain in 1YMM, then its peptide is transferred into the fixed
experimental HLA and TCR assembly.  It is explicitly not flexible docking,
energy minimization, affinity prediction, or evidence of cross-reactivity.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from assemble_ob1a12_template_models import (
    AA3, OUT as TRANSFER_OUT, ca, global_pairs, heavy_atoms, kabsch,
    nearest_count, residue_records, sequence, transformed_line,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "ob1a12_ebv_scaffold_test"
REFERENCE = TRANSFER_OUT / "1YMM_experimental_reference.pdb"
CASES = {
    "positive_mbp_reconstruction": "HUMAN_MYELIN_114806",
    "ebv_mimic_hypothesis": "EBV_TCELL_63843",
    "negative_ebv_control": "EBV_TCELL_2268933",
}


def find_rank1(folder: Path, candidate: str) -> Path:
    matches = sorted(folder.glob(f"{candidate}_*rank_001*.pdb"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one rank-1 PDB for {candidate}, found {len(matches)}")
    return matches[0]


def chain_by_sequence(target: dict[str, list[dict]], reference: list[dict]) -> str:
    """Identify target DRB1 by maximal exact aligned positions."""
    return max(target, key=lambda name: len(global_pairs(sequence(reference), sequence(target[name]))))


def ca_rmsd(first: list[dict], second: list[dict]) -> float | None:
    points = [(ca(a), ca(b)) for a, b in zip(first, second)]
    points = [(a, b) for a, b in points if a is not None and b is not None]
    if not points:
        return None
    return float(np.sqrt(np.mean([np.sum((a - b) ** 2) for a, b in points])))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Provide the folder containing rank-1 PDB files.")
    folder = Path(sys.argv[1]).expanduser()
    if not REFERENCE.exists():
        raise SystemExit("Run assemble_ob1a12_template_models.py first to fetch 1YMM.")
    OUT.mkdir(parents=True, exist_ok=True)
    reference_lines = REFERENCE.read_text().splitlines()
    reference = residue_records(reference_lines)
    ref_drb, ref_peptide = reference["B"], reference["C"]
    ref_pmhc_lines = [line for line in reference_lines if line.startswith("ATOM") and line[21] in {"A", "B"}]
    ref_tcr_lines = [line for line in reference_lines if line.startswith("ATOM") and line[21] in {"D", "E"}]
    ref_hla_atoms = heavy_atoms(reference["A"] + reference["B"])
    ref_tcr_atoms = heavy_atoms(reference["D"] + reference["E"])
    rows = []
    for label, candidate in CASES.items():
        target_path = find_rank1(folder, candidate)
        target = residue_records(target_path.read_text(errors="ignore").splitlines())
        drb_name = chain_by_sequence(target, ref_drb)
        peptide_name = min(target, key=lambda name: len(target[name]))
        source_points, target_points = [], []
        for ref_i, target_i in global_pairs(sequence(ref_drb), sequence(target[drb_name])):
            ref_point, target_point = ca(ref_drb[ref_i]), ca(target[drb_name][target_i])
            if ref_point is not None and target_point is not None:
                source_points.append(target_point)
                target_points.append(ref_point)
        rotation, translation, drb_rmsd = kabsch(np.array(source_points), np.array(target_points))
        peptide_lines, peptide_atoms = [], []
        for residue in target[peptide_name]:
            for atom, xyz, line in residue["atoms"]:
                peptide_lines.append(transformed_line(line, rotation, translation, "C"))
                if not atom.startswith("H"):
                    peptide_atoms.append(xyz @ rotation + translation)
        # Positive control has the same core as 1YMM, so this validates transfer precision.
        transformed_residues = residue_records(peptide_lines).get("C", [])
        mbp_core_rmsd = ca_rmsd(transformed_residues, ref_peptide[:len(transformed_residues)]) if label.startswith("positive") else None
        out_path = OUT / f"{label}__{candidate}__1ymm_scaffold.pdb"
        out_path.write_text("\n".join(ref_pmhc_lines + ["TER"] + peptide_lines + ["TER"] + ref_tcr_lines + ["TER", "END"]) + "\n")
        rows.append({
            "condition": label,
            "candidate_id": candidate,
            "scaffold_pdb": str(out_path),
            "drb1_alignment_ca_pairs": len(source_points),
            "target_drb1_to_1ymm_drb1_ca_rmsd_angstrom": round(drb_rmsd, 3),
            "positive_control_peptide_ca_rmsd_to_1ymm_angstrom": None if mbp_core_rmsd is None else round(mbp_core_rmsd, 3),
            "peptide_hla_heavy_atom_clashes_le_2_0A": nearest_count(peptide_atoms, ref_hla_atoms, 2.0),
            "tcr_peptide_heavy_atom_contacts_le_4_5A": nearest_count(peptide_atoms, ref_tcr_atoms, 4.5),
            "tcr_peptide_heavy_atom_clashes_le_2_0A": nearest_count(peptide_atoms, ref_tcr_atoms, 2.0),
            "interpretation": "Fixed-template compatibility screen only; no flexible docking or affinity inference.",
        })
    with (OUT / "ob1a12_ebv_scaffold_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
