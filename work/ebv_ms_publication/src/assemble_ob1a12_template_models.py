"""Make *initial*, rigid-body TCR hypotheses from the experimental 1YMM complex.

1YMM is an X-ray structure of the human Ob.1A12 TCR bound to
HLA-DRA*01:01/DRB1*15:01 and the MBP85-99 core ENPVVHFFKNIVTP.
This script aligns the HLA chains of each completed ColabFold pMHC model to
1YMM and transfers only the experimental TCR coordinates.  It does *not* dock
or relax the TCR and its contact/clash values must not be interpreted as an
affinity or cross-reactivity prediction.

The files are useful starting geometries for a later restrained-docking or
MD stage, with an exact experimental positive control.

Usage:
  python3 src/assemble_ob1a12_template_models.py /path/to/PDBs
"""

from __future__ import annotations

import csv
import sys
import urllib.request
from collections import OrderedDict
from itertools import permutations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "processed" / "ob1a12_template_transfer"
REFERENCE_URL = "https://files.rcsb.org/download/1YMM.pdb"
CASES = {
    "positive_mbp": "HUMAN_MYELIN_114806",
    "test_ebv_mimic": "EBV_TCELL_63843",
    "negative_ebv_control": "EBV_TCELL_2268933",
}


def residue_records(lines: list[str]) -> dict[str, list[dict]]:
    """Parse protein residues and retain original ATOM lines for writing."""
    chains: dict[str, OrderedDict[tuple[str, str], dict]] = {}
    for line in lines:
        if not line.startswith("ATOM") or line[16:17] not in {" ", "A"}:
            continue
        chain = line[21].strip() or "_"
        key = (line[22:26].strip(), line[26].strip())
        residues = chains.setdefault(chain, OrderedDict())
        residue = residues.setdefault(key, {"name": line[17:20].strip(), "atoms": []})
        atom = line[12:16].strip()
        try:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        residue["atoms"].append((atom, xyz, line))
    return {chain: list(items.values()) for chain, items in chains.items()}


AA3 = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C", "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P", "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V"}


def sequence(residues: list[dict]) -> str:
    return "".join(AA3.get(record["name"], "X") for record in residues)


def global_pairs(left: str, right: str) -> list[tuple[int, int]]:
    """Needleman-Wunsch index pairs; enough for the nearly identical HLA chains."""
    m, n = len(left), len(right)
    score = np.zeros((m + 1, n + 1), dtype=int)
    score[:, 0] = np.arange(m + 1) * -2
    score[0, :] = np.arange(n + 1) * -2
    trace = np.zeros((m + 1, n + 1), dtype=np.int8)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            options = (score[i - 1, j - 1] + (2 if left[i - 1] == right[j - 1] else -1), score[i - 1, j] - 2, score[i, j - 1] - 2)
            trace[i, j] = int(np.argmax(options))
            score[i, j] = max(options)
    pairs: list[tuple[int, int]] = []
    i, j = m, n
    while i or j:
        direction = trace[i, j] if i and j else (1 if i else 2)
        if direction == 0:
            if left[i - 1] == right[j - 1]:
                pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif direction == 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def matched_hla_order(reference_hla: list[list[dict]], target: dict[str, list[dict]], chains: list[str]) -> list[str]:
    """Map DRA and DRB by sequence identity, never by the PDB chain order."""
    best_order, best_matches = None, -1
    for order in permutations(chains):
        matches = sum(
            len(global_pairs(sequence(ref_chain), sequence(target[target_chain])))
            for ref_chain, target_chain in zip(reference_hla, order)
        )
        if matches > best_matches:
            best_order, best_matches = list(order), matches
    assert best_order is not None
    return best_order


def ca(residue: dict) -> np.ndarray | None:
    return next((xyz for atom, xyz, _ in residue["atoms"] if atom == "CA"), None)


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    source_center, target_center = source.mean(axis=0), target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(covariance)
    # Points are row vectors, so source @ rotation maps onto target.
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1] *= -1
        rotation = u @ vt
    translated = (source - source_center) @ rotation + target_center
    rmsd = float(np.sqrt(np.mean(np.sum((translated - target) ** 2, axis=1))))
    return rotation, target_center - source_center @ rotation, rmsd


def transformed_line(line: str, rotation: np.ndarray, translation: np.ndarray, chain: str) -> str:
    xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    x, y, z = xyz @ rotation + translation
    return f"{line[:21]}{chain}{line[22:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"


def heavy_atoms(residues: list[dict]) -> list[np.ndarray]:
    return [xyz for residue in residues for atom, xyz, _ in residue["atoms"] if not atom.startswith("H")]


def nearest_count(query: list[np.ndarray], reference: list[np.ndarray], cutoff: float) -> int:
    """Count query atoms having a neighbor, in bounded NumPy batches."""
    query_array, reference_array = np.asarray(query), np.asarray(reference)
    count = 0
    for start in range(0, len(query_array), 256):
        chunk = query_array[start:start + 256]
        squared = np.sum((chunk[:, None, :] - reference_array[None, :, :]) ** 2, axis=2)
        count += int(np.sum(np.any(squared <= cutoff ** 2, axis=1)))
    return count


def find_rank1(folder: Path, candidate: str) -> Path:
    matches = sorted(folder.glob(f"{candidate}_*rank_001*.pdb"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one rank-1 PDB for {candidate}, found {len(matches)}")
    return matches[0]


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Provide the folder containing the 86 rank-1 PDBs.")
    folder = Path(sys.argv[1]).expanduser()
    OUT.mkdir(parents=True, exist_ok=True)
    ref_path = OUT / "1YMM_experimental_reference.pdb"
    if not ref_path.exists():
        urllib.request.urlretrieve(REFERENCE_URL, ref_path)
    ref_lines = ref_path.read_text().splitlines()
    ref = residue_records(ref_lines)
    # 1YMM: A=DRA, B=DRB1*15:01, C=MBP, D/E=Ob.1A12 alpha/beta.
    reference_hla = [ref["A"], ref["B"]]
    rows = []
    for label, candidate in CASES.items():
        target_path = find_rank1(folder, candidate)
        target_lines = target_path.read_text(errors="ignore").splitlines()
        target = residue_records(target_lines)
        chain_order = sorted(target, key=lambda name: len(target[name]), reverse=True)
        hla_chains = matched_hla_order(reference_hla, target, chain_order[:2])
        peptide_chain = chain_order[2]
        source_ca, target_ca = [], []
        for reference_chain, target_chain in zip(reference_hla, hla_chains):
            for left, right in global_pairs(sequence(reference_chain), sequence(target[target_chain])):
                source_point, target_point = ca(reference_chain[left]), ca(target[target_chain][right])
                if source_point is not None and target_point is not None:
                    source_ca.append(source_point)
                    target_ca.append(target_point)
        rotation, translation, hla_rmsd = kabsch(np.array(source_ca), np.array(target_ca))
        # Keep the target pMHC intact and append the experimentally observed TCR.
        output_lines = [line for line in target_lines if line.startswith("ATOM")]
        tcr_atoms = []
        for ref_chain, out_chain in (("D", "D"), ("E", "E")):
            for residue in ref[ref_chain]:
                for atom, xyz, line in residue["atoms"]:
                    transformed = transformed_line(line, rotation, translation, out_chain)
                    output_lines.append(transformed)
                    if not atom.startswith("H"):
                        tcr_atoms.append(xyz @ rotation + translation)
            output_lines.append("TER")
        output_lines.append("END")
        output_path = OUT / f"{label}__{candidate}__ob1a12_template_transfer.pdb"
        output_path.write_text("\n".join(output_lines) + "\n")
        peptide_atoms = heavy_atoms(target[peptide_chain])
        hla_atoms = heavy_atoms(target[hla_chains[0]] + target[hla_chains[1]])
        rows.append({
            "condition": label,
            "candidate_id": candidate,
            "input_pmhc_rank1_pdb": str(target_path),
            "output_initial_complex_pdb": str(output_path),
            "template": "1YMM (experimental Ob.1A12--DRB1*15:01--MBP)",
            "hla_alignment_ca_pairs": len(source_ca),
            "template_to_target_hla_ca_rmsd_angstrom": round(hla_rmsd, 3),
            "eligible_as_initial_tcr_geometry": hla_rmsd <= 6.0,
            "tcr_peptide_heavy_atom_contacts_le_4_5A": nearest_count(tcr_atoms, peptide_atoms, 4.5),
            "tcr_pmhc_heavy_atom_clashes_le_2_0A": nearest_count(tcr_atoms, hla_atoms + peptide_atoms, 2.0),
            "interpretation": "Rigid-body template transfer only; not a docking score, affinity estimate, or cross-reactivity result. A model with HLA RMSD >6 A is retained for audit but excluded from docking.",
        })
    with (OUT / "ob1a12_template_transfer_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
