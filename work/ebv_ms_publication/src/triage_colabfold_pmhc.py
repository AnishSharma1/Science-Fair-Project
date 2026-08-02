"""Coordinate-level triage for completed ColabFold HLA-DRB1*15:01 pMHC models.

This is deliberately a *pMHC quality and comparability* screen.  It does not
infer TCR cross-reactivity.  It reports peptide-level confidence and how much
of each peptide is buried against HLA before any TCR-docking claim is made.

Usage:
  python3 src/triage_colabfold_pmhc.py /path/to/ebv_ms_pmhc_batch
"""

from __future__ import annotations

import csv
import math
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "processed" / "pmhc_candidate_manifest.csv"
OUT = ROOT / "processed"
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# These pairs were selected before structure inspection from assay evidence,
# explicit candidate provenance, and transparent exploratory resemblance.
PAIRS = [
    ("EBV_TCELL_63843", "HUMAN_MYELIN_114806"),
    ("EBV_TCELL_2268933", "HUMAN_MYELIN_5516"),
    ("EBV_TCELL_950", "HUMAN_MYELIN_112214"),
    ("EBV_TCELL_2268934", "HUMAN_MYELIN_112226"),
]


def load_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return {row["candidate_id"]: row for row in csv.DictReader(handle)}


def parse_pdb(path: Path) -> dict[str, list[dict]]:
    """Read residues and atoms from a PDB without trusting the filename."""
    chains: dict[str, OrderedDict[tuple[str, str], dict]] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[16:17] not in {" ", "A"}:
            continue
        residue_name = line[17:20].strip()
        if residue_name not in AA3:
            continue
        chain = line[21].strip() or "_"
        key = (line[22:26].strip(), line[26].strip())
        residues = chains.setdefault(chain, OrderedDict())
        residue = residues.setdefault(key, {"aa": AA3[residue_name], "atoms": [], "bfactors": []})
        atom = line[12:16].strip()
        element = (line[76:78].strip() or atom[:1]).upper()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            bfactor = float(line[60:66])
        except ValueError:
            continue
        residue["atoms"].append((atom, element, xyz))
        residue["bfactors"].append(bfactor)
    return {chain: list(residues.values()) for chain, residues in chains.items()}


def sequence(residues: list[dict]) -> str:
    return "".join(r["aa"] for r in residues)


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def peptide_metrics(peptide: list[dict], hla: list[dict]) -> dict[str, float | int | str]:
    hla_atoms = [atom for residue in hla for atom in residue["atoms"] if atom[1] != "H"]
    per_residue_contacts = []
    per_residue_confidence = []
    for residue in peptide:
        atoms = [atom for atom in residue["atoms"] if atom[1] != "H"]
        contacts = 0
        for _, _, xyz in atoms:
            if any(distance(xyz, hla_atom[2]) <= 4.0 for hla_atom in hla_atoms):
                contacts += 1
        per_residue_contacts.append(contacts)
        per_residue_confidence.append(sum(residue["bfactors"]) / len(residue["bfactors"]))
    # Lower contact counts are an exposure proxy only; a real SASA calculation
    # is a later refinement, not something to substitute with a label here.
    return {
        "peptide_sequence_from_pdb": sequence(peptide),
        "peptide_residues": len(peptide),
        "peptide_mean_plddt": round(sum(per_residue_confidence) / len(per_residue_confidence), 2),
        "peptide_min_plddt": round(min(per_residue_confidence), 2),
        "mean_hla_contacting_atoms_per_residue": round(sum(per_residue_contacts) / len(per_residue_contacts), 2),
        "lowest_contact_positions_1based": ";".join(
            str(i + 1) for i, value in enumerate(per_residue_contacts)
            if value == min(per_residue_contacts)
        ),
    }


def residue_detail_rows(candidate: str, peptide: list[dict], hla: list[dict]) -> list[dict]:
    """Return transparent per-residue values; no residue is called TCR-facing."""
    hla_atoms = [atom for residue in hla for atom in residue["atoms"] if atom[1] != "H"]
    rows = []
    for position, residue in enumerate(peptide, start=1):
        atoms = [atom for atom in residue["atoms"] if atom[1] != "H"]
        contacts = sum(
            any(distance(xyz, hla_atom[2]) <= 4.0 for hla_atom in hla_atoms)
            for _, _, xyz in atoms
        )
        rows.append({
            "candidate_id": candidate,
            "peptide_position_1based": position,
            "residue": residue["aa"],
            "residue_mean_plddt": round(sum(residue["bfactors"]) / len(residue["bfactors"]), 2),
            "hla_contacting_heavy_atoms": contacts,
        })
    return rows


def ca_coordinates(residues: list[dict]) -> np.ndarray:
    coords = []
    for residue in residues:
        ca = next((xyz for atom, _, xyz in residue["atoms"] if atom == "CA"), None)
        if ca is None:
            raise ValueError("A residue lacks a CA atom")
        coords.append(ca)
    return np.asarray(coords, dtype=float)


def kabsch(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit source to target; return rotation, translation, and fitted RMSD."""
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (source - source_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(covariance)
    # Points are row vectors, so source @ rotation maps onto target.
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    translation = target_center - source_center @ rotation
    fitted = source @ rotation + translation
    return rotation, translation, float(np.sqrt(np.mean(np.sum((fitted - target) ** 2, axis=1))))


def local_alignment_indices(a: str, b: str) -> list[tuple[int, int]]:
    """Documented Smith-Waterman alignment with residue positions retained."""
    score = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    trace = [[None] * (len(b) + 1) for _ in range(len(a) + 1)]
    best = (0, 0, 0)
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            choices = [
                (0, None),
                (score[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1), "diag"),
                (score[i - 1][j] - 2, "up"),
                (score[i][j - 1] - 2, "left"),
            ]
            score[i][j], trace[i][j] = max(choices, key=lambda item: item[0])
            if score[i][j] > best[0]:
                best = (score[i][j], i, j)
    _, i, j = best
    pairs = []
    while i > 0 and j > 0 and score[i][j] > 0:
        step = trace[i][j]
        if step == "diag":
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif step == "up":
            i -= 1
        elif step == "left":
            j -= 1
        else:
            break
    return list(reversed(pairs))


def pair_geometry(ebv: dict, human: dict) -> dict:
    """Compare homologous HLA frames and only the transparent local alignment."""
    ebv_hla = ebv["hla"]
    human_hla = human["hla"]
    # Only the N-terminal peptide-binding platform of each class-II chain is
    # used as the frame.  Distal immunoglobulin-like domains can move relative
    # to the groove and would otherwise obscure a peptide-pose comparison.
    ebv_ca = np.vstack([ca_coordinates(chain[:85]) for chain in ebv_hla])
    human_ca = np.vstack([ca_coordinates(chain[:85]) for chain in human_hla])
    if ebv_ca.shape != human_ca.shape:
        return {"status": "HLA_CHAIN_LENGTH_MISMATCH"}
    rotation, translation, hla_rmsd = kabsch(ebv_ca, human_ca)
    pairs = local_alignment_indices(sequence(ebv["peptide"]), sequence(human["peptide"]))
    ebv_peptide_ca = ca_coordinates(ebv["peptide"])
    human_peptide_ca = ca_coordinates(human["peptide"])
    fitted_ebv = ebv_peptide_ca @ rotation + translation
    if not pairs:
        return {"status": "NO_LOCAL_ALIGNMENT", "hla_groove_ca_rmsd_after_fit": round(hla_rmsd, 3)}
    distances = [distance(tuple(fitted_ebv[i]), tuple(human_peptide_ca[j])) for i, j in pairs]
    aligned = ";".join(
        f"{i + 1}{sequence(ebv['peptide'])[i]}:{j + 1}{sequence(human['peptide'])[j]}"
        for i, j in pairs
    )
    return {
        "status": "PASS",
        "hla_groove_ca_rmsd_after_fit": round(hla_rmsd, 3),
        "locally_aligned_peptide_residues": len(pairs),
        "aligned_positions_ebv_to_human": aligned,
        "local_peptide_ca_rmsd_after_hla_fit": round(float(np.sqrt(np.mean(np.square(distances)))), 3),
    }


def locate_rank1_models(batch_dir: Path, requested: set[str]) -> dict[str, Path]:
    models: dict[str, Path] = {}
    # The Drive-synced batch is very large.  Do not recursively enumerate or
    # hydrate all 86 models when the present decision concerns four pairs.
    for candidate in requested:
        hits = sorted(batch_dir.glob(f"{candidate}_unrelaxed_rank_001_*.pdb"))
        if len(hits) == 1:
            models[candidate] = hits[0]
    return models


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Pass the downloaded ebv_ms_pmhc_batch folder as the only argument.")
    batch_dir = Path(sys.argv[1]).expanduser().resolve()
    if not batch_dir.is_dir():
        raise SystemExit(f"Not a folder: {batch_dir}")
    manifest = load_manifest()
    # Screen every completed rank-1 pMHC before any biological shortlist is
    # chosen.  PAIRS is only used later to retain provenance for the initial
    # hypothesis pairs; it must never limit the quality-control population.
    requested = set(manifest)
    models = locate_rank1_models(batch_dir, requested)
    if not models:
        raise SystemExit("No rank-1 ColabFold PDB files found in that folder.")

    rows = []
    residue_rows = []
    parsed = {}
    for candidate, path in sorted(models.items()):
        chains = parse_pdb(path)
        expected = manifest[candidate]["peptide"]
        peptide_chain = next((chain for chain, residues in chains.items() if sequence(residues) == expected), None)
        # ColabFold can emit the alpha and beta chains in different record
        # orders.  Their lengths distinguish DRA from DRB here; normalize the
        # order before any inter-model coordinate comparison.
        long_chains = sorted(
            (residues for chain, residues in chains.items() if chain != peptide_chain and len(residues) >= 150),
            key=lambda residues: (len(residues), sequence(residues)),
        )
        row = {"candidate_id": candidate, "pdb_path": str(path), "expected_peptide": expected}
        if peptide_chain is None or len(long_chains) != 2:
            row.update({"status": "FAILED_LAYOUT_OR_SEQUENCE_QA"})
        else:
            row.update({"status": "PASS", "peptide_chain": peptide_chain})
            row.update(peptide_metrics(chains[peptide_chain], long_chains[0] + long_chains[1]))
            residue_rows.extend(residue_detail_rows(candidate, chains[peptide_chain], long_chains[0] + long_chains[1]))
            parsed[candidate] = {"peptide": chains[peptide_chain], "hla": long_chains}
        rows.append(row)

    fieldnames = sorted({key for row in rows for key in row})
    out = OUT / "colabfold_pmhc_peptide_qa.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    residue_out = OUT / "colabfold_pmhc_peptide_residue_qa.csv"
    with residue_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(residue_rows[0]))
        writer.writeheader()
        writer.writerows(residue_rows)

    available = {row["candidate_id"] for row in rows if row["status"] == "PASS"}
    pair_rows = []
    for ebv, human in PAIRS:
        pair_rows.append({
            "ebv_candidate_id": ebv,
            "human_candidate_id": human,
            "both_rank1_pdbs_present_and_sequence_validated": ebv in available and human in available,
            "next_decision": "inspect aligned peptide poses and TCR-facing residues" if ebv in available and human in available else "obtain or repair missing PDB",
            **(pair_geometry(parsed[ebv], parsed[human]) if ebv in available and human in available else {}),
        })
    pair_out = OUT / "colabfold_pmhc_pair_triage.csv"
    with pair_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    # Full, unbiased candidate matrix: every EBV peptide with direct T-cell
    # evidence against every modeled human-myelin peptide.  This is the input
    # for biological prioritization; the four PAIRS above are not privileged.
    qa_by_id = {row["candidate_id"]: row for row in rows}
    tier1_ebv = sorted(candidate for candidate in available if candidate.startswith("EBV_TCELL_"))
    myelin = sorted(candidate for candidate in available if candidate.startswith("HUMAN_MYELIN_"))
    matrix_rows = []
    for ebv in tier1_ebv:
        for human in myelin:
            geometry = pair_geometry(parsed[ebv], parsed[human])
            matrix_rows.append({
                "ebv_candidate_id": ebv,
                "human_candidate_id": human,
                "ebv_peptide_mean_plddt": qa_by_id[ebv].get("peptide_mean_plddt"),
                "human_peptide_mean_plddt": qa_by_id[human].get("peptide_mean_plddt"),
                **geometry,
            })
    matrix_out = OUT / "colabfold_tier1_ebv_myelin_geometry_matrix.csv"
    with matrix_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in matrix_rows for key in row}))
        writer.writeheader()
        writer.writerows(matrix_rows)
    print(f"Wrote {out}")
    print(f"Wrote {pair_out}")
    print(f"Wrote {residue_out}")
    print(f"Wrote {matrix_out}")
    print(f"Validated rank-1 pMHC models: {sum(row['status'] == 'PASS' for row in rows)}/{len(rows)}")


if __name__ == "__main__":
    main()
