"""Evaluate every available rank for the three Ob.1A12 ternary conditions."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

import evaluate_ob1a12_ternary_models as single


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "ob1a12_ternary_evaluation"
CASES = {"positive_control": "positive", "ebv_hypothesis": "test", "negative_control": "negative"}


def rank_number(path: Path) -> int:
    match = re.search(r"rank_?0*(\d+)", path.name.lower())
    if not match:
        raise ValueError(f"No rank in {path.name}")
    return int(match.group(1))


def evaluate(path: Path, condition: str, reference: dict[str, list[dict]]) -> dict:
    model = single.parse_pdb(path)
    peptide = next(name for name, residues in model.items() if len(residues) < 30)
    hla = single.identify_hla_chains(model, peptide, reference)
    tcr = [name for name in model if name not in hla + [peptide]]
    tcr = single.best_mapping([reference["D"], reference["E"]], model, tcr)
    rotation, translation, hla_rmsd = single.fit_to_reference(model, reference, hla)
    tcr_rmsds = [single.fitted_chain_rmsd(model[chain], ref, rotation, translation) for chain, ref in zip(tcr, (reference["D"], reference["E"]))]
    peptide_rmsd = None
    if condition == "positive_control":
        peptide_rmsd = single.fitted_chain_rmsd(model[peptide], reference["C"][:len(model[peptide])], rotation, translation)
    tcr_atoms = single.heavy_atoms(model[tcr[0]] + model[tcr[1]])
    peptide_atoms = single.heavy_atoms(model[peptide])
    return {
        "condition": condition,
        "rank": rank_number(path),
        "pdb": str(path),
        "peptide_mean_plddt": single.mean_plddt(model[peptide]),
        "tcr_mean_plddt": round((single.mean_plddt(model[tcr[0]]) + single.mean_plddt(model[tcr[1]])) / 2, 2),
        "hla_fit_to_1ymm_ca_rmsd_angstrom": round(hla_rmsd, 3),
        "tcr_alpha_fit_to_1ymm_ca_rmsd_angstrom": round(tcr_rmsds[0], 3),
        "tcr_beta_fit_to_1ymm_ca_rmsd_angstrom": round(tcr_rmsds[1], 3),
        "positive_control_peptide_fit_to_1ymm_ca_rmsd_angstrom": None if peptide_rmsd is None else round(peptide_rmsd, 3),
        "tcr_peptide_heavy_atom_contacts_le_4_5A": single.neighbor_count(peptide_atoms, tcr_atoms, 4.5),
        "tcr_peptide_heavy_atom_clashes_le_2_0A": single.neighbor_count(peptide_atoms, tcr_atoms, 2.0),
    }


def median(rows: list[dict], key: str) -> float:
    return round(float(np.median([float(row[key]) for row in rows if row[key] not in {None, ""}])), 3)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Provide the TCR PDBs folder.")
    folder = Path(sys.argv[1]).expanduser()
    reference = single.parse_pdb(single.REFERENCE)
    all_rows = []
    for condition, token in CASES.items():
        paths = sorted((path for path in folder.rglob("*.pdb") if token in path.name.lower()), key=rank_number)
        if len(paths) < 3:
            raise FileNotFoundError(f"Expected at least three {condition} PDBs, found {len(paths)}")
        all_rows.extend(evaluate(path, condition, reference) for path in paths)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "ob1a12_ternary_all_ranks.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)
    summary = []
    for condition in CASES:
        rows = [row for row in all_rows if row["condition"] == condition]
        summary.append({
            "condition": condition,
            "n_models": len(rows),
            "median_peptide_plddt": median(rows, "peptide_mean_plddt"),
            "median_tcr_plddt": median(rows, "tcr_mean_plddt"),
            "median_hla_fit_rmsd_angstrom": median(rows, "hla_fit_to_1ymm_ca_rmsd_angstrom"),
            "median_tcr_alpha_fit_rmsd_angstrom": median(rows, "tcr_alpha_fit_to_1ymm_ca_rmsd_angstrom"),
            "median_tcr_beta_fit_rmsd_angstrom": median(rows, "tcr_beta_fit_to_1ymm_ca_rmsd_angstrom"),
            "median_tcr_peptide_contacts": median(rows, "tcr_peptide_heavy_atom_contacts_le_4_5A"),
            "median_tcr_peptide_clashes": median(rows, "tcr_peptide_heavy_atom_clashes_le_2_0A"),
            "interpretation": "Ensemble comparison only. The positive control must recover the experimental TCR orientation before relative contacts can support a structural hypothesis.",
        })
    with (OUT / "ob1a12_ternary_ensemble_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    print(f"Wrote {OUT / 'ob1a12_ternary_all_ranks.csv'}")
    print(f"Wrote {OUT / 'ob1a12_ternary_ensemble_summary.csv'}")


if __name__ == "__main__":
    main()
