"""Coordinate-level QA for legacy MHC-II pMHC/PDB files.

The script infers the peptide from the short standard-amino-acid chain and
never trusts filename labels as sequence evidence.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
PROJECT_DATA = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/"
    "Downloads/Projects and Data"
)
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

import sys
sys.path.insert(0, str(ROOT / "src"))
from explore_similarity import smith_waterman


def bucket(path: Path) -> str:
    s = str(path)
    if "MHC Docking Analysis" in s:
        return "legacy_docking_complex"
    if "Tetramer Docking Analysis" in s:
        return "legacy_tetramer_export"
    return "legacy_predicted_pmhc"


def parse_sequences(path: Path) -> tuple[dict[str, str], dict[str, list[float]]]:
    residues = defaultdict(list)
    b_factors = defaultdict(list)
    for line in path.read_text(errors="ignore").splitlines():
        if line[:6].strip() != "ATOM":
            continue
        chain = line[21].strip() or "_"
        residue = (line[22:26].strip(), line[26].strip(), line[17:20].strip())
        if residue not in residues[chain]:
            residues[chain].append(residue)
        try:
            b_factors[chain].append(float(line[60:66]))
        except ValueError:
            pass
    seqs = {chain: "".join(AA3.get(r[2], "X") for r in vals) for chain, vals in residues.items()}
    return seqs, b_factors


def parse(path: Path, template_alpha: str, template_beta: str) -> dict:
    seqs, b_factors = parse_sequences(path)
    peptide_candidates = [c for c, seq in seqs.items() if 9 <= len(seq) <= 40 and set(seq) <= set(AA3.values())]
    peptide_chain = sorted(peptide_candidates, key=lambda c: (len(seqs[c]), c))[0] if peptide_candidates else ""
    peptide = seqs.get(peptide_chain, "")
    other_lengths = sorted(len(seq) for c, seq in seqs.items() if c != peptide_chain)
    clean_layout = len(seqs) >= 3 and len(other_lengths) >= 2 and all(200 <= n <= 350 for n in other_lengths[:2])
    role_counts = {"DRA_like": 0, "DRB_like": 0}
    for chain, sequence in seqs.items():
        if chain == peptide_chain or len(sequence) < 40:
            continue
        _, _, alpha_aligned, alpha_identity = smith_waterman(sequence, template_alpha)
        _, _, beta_aligned, beta_identity = smith_waterman(sequence, template_beta)
        if alpha_aligned >= 150 and alpha_identity >= 0.8 and alpha_identity > beta_identity:
            role_counts["DRA_like"] += 1
        if beta_aligned >= 150 and beta_identity >= 0.8 and beta_identity > alpha_identity:
            role_counts["DRB_like"] += 1
    proper_class_ii_layout = role_counts["DRA_like"] == 1 and role_counts["DRB_like"] == 1 and bool(peptide_chain)
    manifest = pd.read_csv(PROC / "pmhc_candidate_manifest.csv")
    hits = manifest[manifest.peptide == peptide]
    if clean_layout and proper_class_ii_layout and hits.shape[0] and bucket(path) == "legacy_predicted_pmhc":
        status = "validated_reference_structure"
    elif clean_layout and hits.shape[0] and role_counts["DRB_like"] >= 2 and role_counts["DRA_like"] == 0:
        status = "peptide_match_but_missing_alpha_chain"
    elif len(seqs) < 3:
        status = "incomplete_structure"
    elif bucket(path) != "legacy_predicted_pmhc":
        status = "legacy_ambiguous_complex"
    elif not hits.shape[0]:
        status = "unmapped_peptide"
    else:
        status = "layout_requires_review"
    return {
        "structure_path": str(path),
        "source_bucket": bucket(path),
        "filename": path.name,
        "chain_ids": ":".join(sorted(seqs)),
        "chain_lengths": json.dumps({c: len(s) for c, s in sorted(seqs.items())}, sort_keys=True),
        "peptide_chain": peptide_chain,
        "peptide_sequence": peptide,
        "candidate_ids": ";".join(hits.candidate_id.astype(str).tolist()),
        "hla": ";".join(sorted(set(hits.hla.astype(str)))) if len(hits) else "",
        "molecule_role_counts": json.dumps(role_counts, sort_keys=True),
        "proper_class_ii_layout": proper_class_ii_layout,
        "mean_peptide_bfactor": round(sum(b_factors.get(peptide_chain, [])) / len(b_factors.get(peptide_chain, [])), 3) if b_factors.get(peptide_chain) else None,
        "mean_all_bfactor": round(sum(v for vals in b_factors.values() for v in vals) / sum(len(v) for v in b_factors.values()), 3) if b_factors else None,
        "status": status,
    }


def main() -> None:
    paths = sorted(p for p in PROJECT_DATA.rglob("*.pdb") if "MHCII_" in p.name)
    template_sequences, _ = parse_sequences(ROOT / "raw/rcsb_templates/8TBP.pdb")
    template_alpha = template_sequences["A"]
    template_beta = template_sequences["B"]
    rows = [parse(path, template_alpha, template_beta) for path in paths]
    df = pd.DataFrame(rows)
    df.to_csv(PROC / "pmhc_structure_qa.csv", index=False)
    summary = df.groupby(["source_bucket", "status"]).size().reset_index(name="n")
    summary.to_csv(PROC / "pmhc_structure_qa_summary.csv", index=False)
    lines = ["# pMHC structure QA summary", "", "| source_bucket | status | n |", "|---|---|---:|"]
    lines.extend(f"| {r.source_bucket} | {r.status} | {int(r.n)} |" for r in summary.itertuples())
    (PROC / "pmhc_structure_qa_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
