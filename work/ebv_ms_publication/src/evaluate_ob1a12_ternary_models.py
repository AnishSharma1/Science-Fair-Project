"""Evaluate the three Ob.1A12 ternary ColabFold outputs against controls.

Usage:
  python3 src/evaluate_ob1a12_ternary_models.py /path/to/three_model_results

It accepts the rank-1 PDB for each of the three FASTA records in
processed/ob1a12_ternary_colabfold_inputs.fasta.  This produces a transparent
comparison table; it never converts a structure-prediction score into proof of
TCR binding or cross-reactivity.
"""

from __future__ import annotations

import csv
import sys
from collections import OrderedDict
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

from assemble_ob1a12_template_models import AA3, global_pairs, kabsch, sequence


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "processed" / "ob1a12_template_transfer" / "1YMM_experimental_reference.pdb"
OUT = ROOT / "processed" / "ob1a12_ternary_evaluation"
CASES = {
    "positive_control": "positive_control__ob1a12__dr15__mbp_85_98",
    "ebv_hypothesis": "test_hypothesis__ob1a12__dr15__ebv_balf5_627_641",
    "negative_control": "negative_control__ob1a12__dr15__ebv_balf4",
}


def parse_pdb(path: Path) -> dict[str, list[dict]]:
    chains: dict[str, OrderedDict[tuple[str, str], dict]] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[16:17] not in {" ", "A"}:
            continue
        name = line[17:20].strip()
        if name not in AA3:
            continue
        chain = line[21].strip() or "_"
        key = (line[22:26].strip(), line[26].strip())
        residues = chains.setdefault(chain, OrderedDict())
        residue = residues.setdefault(key, {"name": name, "atoms": [], "bfactors": []})
        try:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            bfactor = float(line[60:66])
        except ValueError:
            continue
        residue["atoms"].append((line[12:16].strip(), xyz))
        residue["bfactors"].append(bfactor)
    return {chain: list(records.values()) for chain, records in chains.items()}


def ca(residue: dict) -> np.ndarray | None:
    return next((xyz for atom, xyz in residue["atoms"] if atom == "CA"), None)


def mean_plddt(residues: list[dict]) -> float:
    values = [value for residue in residues for value in residue["bfactors"]]
    return round(sum(values) / len(values), 2)


def heavy_atoms(residues: list[dict]) -> np.ndarray:
    return np.asarray([xyz for residue in residues for atom, xyz in residue["atoms"] if not atom.startswith("H")])


def neighbor_count(query: np.ndarray, reference: np.ndarray, cutoff: float) -> int:
    count = 0
    for start in range(0, len(query), 256):
        distances2 = np.sum((query[start:start + 256, None, :] - reference[None, :, :]) ** 2, axis=2)
        count += int(np.sum(np.any(distances2 <= cutoff ** 2, axis=1)))
    return count


def best_mapping(reference: list[list[dict]], target: dict[str, list[dict]], candidates: list[str]) -> list[str]:
    return list(max(
        permutations(candidates),
        key=lambda order: sum(len(global_pairs(sequence(ref), sequence(target[chain]))) for ref, chain in zip(reference, order)),
    ))


def identify_hla_chains(model: dict[str, list[dict]], peptide: str, reference: dict[str, list[dict]]) -> list[str]:
    """Identify HLA by sequence, so partial experimental chains remain usable in QA."""
    candidates = [name for name in model if name != peptide]
    best_pair, best_score = None, -1
    for pair in combinations(candidates, 2):
        mapped = best_mapping([reference["A"], reference["B"]], model, list(pair))
        score = sum(
            len(global_pairs(sequence(ref), sequence(model[chain])))
            for ref, chain in zip((reference["A"], reference["B"]), mapped)
        )
        if score > best_score:
            best_pair, best_score = mapped, score
    assert best_pair is not None
    return best_pair


def fit_to_reference(target: dict[str, list[dict]], reference: dict[str, list[dict]], target_hla: list[str]) -> tuple[np.ndarray, np.ndarray, float]:
    source, destination = [], []
    for ref_chain, target_chain in zip((reference["A"], reference["B"]), target_hla):
        for ref_index, target_index in global_pairs(sequence(ref_chain), sequence(target[target_chain])):
            ref_point, target_point = ca(ref_chain[ref_index]), ca(target[target_chain][target_index])
            if ref_point is not None and target_point is not None:
                source.append(target_point)
                destination.append(ref_point)
    return kabsch(np.asarray(source), np.asarray(destination))


def fitted_chain_rmsd(target: list[dict], reference: list[dict], rotation: np.ndarray, translation: np.ndarray) -> float | None:
    pairs = []
    for ref_index, target_index in global_pairs(sequence(reference), sequence(target)):
        ref_point, target_point = ca(reference[ref_index]), ca(target[target_index])
        if ref_point is not None and target_point is not None:
            pairs.append((target_point @ rotation + translation, ref_point))
    if not pairs:
        return None
    return float(np.sqrt(np.mean([np.sum((a - b) ** 2) for a, b in pairs])))


def find_rank1(folder: Path, identifier: str) -> Path:
    """Accept both ColabFold filenames and the student's organized folders."""
    exact = sorted(folder.rglob(f"{identifier}*rank_001*.pdb"))
    if exact:
        files = exact
    else:
        token = (
            "positive" if identifier.startswith("positive_control")
            else "negative" if identifier.startswith("negative_control")
            else "test"
        )
        files = sorted(
            path for path in folder.rglob("*.pdb")
            if token in path.name.lower()
            and ("rank_001" in path.name.lower() or "rank001" in path.name.lower())
        )
    if len(files) != 1:
        raise FileNotFoundError(f"Expected one rank-1 PDB for {identifier}, found {len(files)}")
    return files[0]


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Provide the folder containing the three ternary rank-1 PDBs.")
    if not REFERENCE.exists():
        raise SystemExit("Experimental reference missing; run the template-transfer script first.")
    OUT.mkdir(parents=True, exist_ok=True)
    reference = parse_pdb(REFERENCE)
    rows = []
    for condition, identifier in CASES.items():
        path = find_rank1(Path(sys.argv[1]).expanduser(), identifier)
        model = parse_pdb(path)
        peptide = next(name for name, residues in model.items() if len(residues) < 30)
        hla_chains = identify_hla_chains(model, peptide, reference)
        tcr_chains = [name for name in model if name not in hla_chains + [peptide]]
        if len(tcr_chains) != 2:
            raise ValueError(f"Could not identify five chains in {path.name}")
        target_hla = hla_chains
        target_tcr = best_mapping([reference["D"], reference["E"]], model, tcr_chains)
        rotation, translation, hla_rmsd = fit_to_reference(model, reference, target_hla)
        tcr_rmsds = [fitted_chain_rmsd(model[chain], ref, rotation, translation) for chain, ref in zip(target_tcr, (reference["D"], reference["E"]))]
        peptide_rmsd = None
        if condition == "positive_control":
            peptide_rmsd = fitted_chain_rmsd(model[peptide], reference["C"][:len(model[peptide])], rotation, translation)
        tcr_atoms, peptide_atoms = heavy_atoms(model[target_tcr[0]] + model[target_tcr[1]]), heavy_atoms(model[peptide])
        rows.append({
            "condition": condition,
            "rank1_pdb": str(path),
            "peptide_mean_plddt": mean_plddt(model[peptide]),
            "tcr_mean_plddt": round((mean_plddt(model[target_tcr[0]]) + mean_plddt(model[target_tcr[1]])) / 2, 2),
            "hla_fit_to_1ymm_ca_rmsd_angstrom": round(hla_rmsd, 3),
            "tcr_alpha_fit_to_1ymm_ca_rmsd_angstrom": None if tcr_rmsds[0] is None else round(tcr_rmsds[0], 3),
            "tcr_beta_fit_to_1ymm_ca_rmsd_angstrom": None if tcr_rmsds[1] is None else round(tcr_rmsds[1], 3),
            "positive_control_peptide_fit_to_1ymm_ca_rmsd_angstrom": None if peptide_rmsd is None else round(peptide_rmsd, 3),
            "tcr_peptide_heavy_atom_contacts_le_4_5A": neighbor_count(peptide_atoms, tcr_atoms, 4.5),
            "tcr_peptide_heavy_atom_clashes_le_2_0A": neighbor_count(peptide_atoms, tcr_atoms, 2.0),
            "interpretation": "Comparative structure-prediction screen only; no result establishes binding or cross-reactivity.",
        })
    with (OUT / "ob1a12_ternary_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT / 'ob1a12_ternary_comparison.csv'}")


if __name__ == "__main__":
    main()
