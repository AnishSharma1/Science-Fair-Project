#!/usr/bin/env python3
"""Template-free CA-atom calibration comparison for downloaded AF3 models."""
import json
from pathlib import Path
import shlex
import numpy as np

ROOT = Path("/Users/anishsharma/Documents/New project")
METRICS = ROOT / "outputs/ebv_ms_model_package/results_analysis/af3_seed_job_metrics.json"
REFS = ROOT / "outputs/ebv_ms_model_package/reference_structures"
OUT = ROOT / "outputs/ebv_ms_model_package/results_analysis"

AA = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}
MAPS = {
    "CAL_1YMM": {"reference":"1YMM", "pairs":{"A":"D","B":"E","C":"A","D":"B","E":"C"}},
    "CAL_1ZGL": {"reference":"1ZGL", "pairs":{"A":"D","B":"E","C":"A","D":"B","E":"C"}},
    "CAL_2WBJ": {"reference":"2WBJ", "pairs":{"A":"C","B":"D","C":"A","D":"B","E":"D"}},
}

def atom_columns(lines):
    for start, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers = []
        pos = start + 1
        while pos < len(lines) and lines[pos].startswith("_atom_site."):
            headers.append(lines[pos].strip())
            pos += 1
        if headers:
            return headers, pos
    raise ValueError("No atom_site loop found")

def ca_chains(cif_path):
    lines = Path(cif_path).read_text().splitlines()
    headers, pos = atom_columns(lines)
    index = {key: idx for idx, key in enumerate(headers)}
    needed = ["_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
    if any(key not in index for key in needed):
        raise ValueError(f"Missing atom fields in {cif_path}")
    chains = {}
    seen = set()
    while pos < len(lines):
        line = lines[pos].strip()
        if not line or line.startswith("#"):
            break
        values = shlex.split(line)
        if len(values) == len(headers):
            atom = values[index["_atom_site.label_atom_id"]]
            comp = values[index["_atom_site.label_comp_id"]]
            chain = values[index["_atom_site.label_asym_id"]]
            seq_id = values[index["_atom_site.label_seq_id"]]
            if atom == "CA" and comp in AA and (chain, seq_id) not in seen:
                seen.add((chain, seq_id))
                chains.setdefault(chain, []).append((AA[comp], np.array([float(values[index["_atom_site.Cartn_x"]]), float(values[index["_atom_site.Cartn_y"]]), float(values[index["_atom_site.Cartn_z"]])])))
        pos += 1
    return chains

def align(seq_a, seq_b):
    m, n = len(seq_a), len(seq_b)
    score = np.zeros((m + 1, n + 1), dtype=int)
    trace = np.zeros((m + 1, n + 1), dtype=np.int8)
    score[1:, 0] = -np.arange(1, m + 1)
    score[0, 1:] = -np.arange(1, n + 1)
    trace[1:, 0] = 1
    trace[0, 1:] = 2
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            opts = (score[i-1, j-1] + (2 if seq_a[i-1] == seq_b[j-1] else -1), score[i-1, j] - 1, score[i, j-1] - 1)
            trace[i, j] = int(np.argmax(opts))
            score[i, j] = opts[trace[i, j]]
    pairs = []
    i, j = m, n
    while i or j:
        direction = trace[i, j]
        if i and j and direction == 0:
            if seq_a[i-1] == seq_b[j-1]: pairs.append((i-1, j-1))
            i, j = i - 1, j - 1
        elif i and (not j or direction == 1):
            i -= 1
        else:
            j -= 1
    return list(reversed(pairs))

def matched_coords(pred, ref, mapping, selected):
    a, b = [], []
    for pred_chain in selected:
        ref_chain = mapping[pred_chain]
        left, right = pred[pred_chain], ref[ref_chain]
        pairs = align("".join(x[0] for x in left), "".join(x[0] for x in right))
        for i, j in pairs:
            a.append(left[i][1]); b.append(right[j][1])
    return np.asarray(a), np.asarray(b)

def transform_fit(mobile, target):
    center_m, center_t = mobile.mean(axis=0), target.mean(axis=0)
    x, y = mobile - center_m, target - center_t
    u, _, vt = np.linalg.svd(x.T @ y)
    d = np.eye(3); d[2, 2] = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ d @ u.T
    return rotation, center_m, center_t

def rmsd(mobile, target, transform=None):
    if transform is None: transform = transform_fit(mobile, target)
    rotation, center_m, center_t = transform
    moved = (mobile - center_m) @ rotation + center_t
    return float(np.sqrt(np.mean(np.sum((moved - target) ** 2, axis=1)))), transform

def main():
    rows = json.loads(METRICS.read_text())
    output = []
    for condition, spec in MAPS.items():
        ref = ca_chains(REFS / f"{spec['reference']}.cif")
        for row in rows:
            if row.get("condition") != condition or row.get("status") != "completed":
                continue
            candidates = list(Path(row["dir"]).glob(f"*_model_{row['bestModel']}.cif"))
            if len(candidates) != 1:
                raise ValueError(f"Expected one best-model CIF for {row['name']}, found {len(candidates)}")
            model = candidates[0]
            pred = ca_chains(model)
            mapping = spec["pairs"]
            pmhc_pred, pmhc_ref = matched_coords(pred, ref, mapping, "CDE")
            tcr_pred, tcr_ref = matched_coords(pred, ref, mapping, "AB")
            all_pred, all_ref = matched_coords(pred, ref, mapping, "ABCDE")
            pmhc_rmsd, pmhc_fit = rmsd(pmhc_pred, pmhc_ref)
            tcr_placement, _ = rmsd(tcr_pred, tcr_ref, pmhc_fit)
            tcr_fold, _ = rmsd(tcr_pred, tcr_ref)
            complex_rmsd, _ = rmsd(all_pred, all_ref)
            output.append({"seed_job":row["name"], "condition":condition, "reference":spec["reference"], "best_model":row["bestModel"], "pMHC_CA_atoms":len(pmhc_pred), "TCR_CA_atoms":len(tcr_pred), "all_CA_atoms":len(all_pred), "pMHC_CA_RMSD_A":pmhc_rmsd, "TCR_placement_CA_RMSD_A":tcr_placement, "TCR_fold_CA_RMSD_A":tcr_fold, "whole_complex_CA_RMSD_A":complex_rmsd, "best_ranking_score":row["bestRanking"], "best_iptm":row["bestIptm"]})
    header = list(output[0])
    def val(x): return f"{x:.3f}" if isinstance(x, float) else str(x)
    (OUT / "af3_calibrator_structural_metrics.tsv").write_text("\t".join(header) + "\n" + "\n".join("\t".join(val(row[key]) for key in header) for row in output) + "\n")
    lines = ["# AF3 calibrator structural comparison", "", "CA-atom comparisons use template-excluded AF3 best-ranked models. pMHC RMSD aligns MHC-alpha, MHC-beta, and peptide; TCR placement RMSD then evaluates TCR alpha/beta in that pMHC-aligned frame.", ""]
    for condition in MAPS:
        use = [r for r in output if r["condition"] == condition]
        if not use:
            lines.append(f"## {condition}\n\nNo completed seeds found.\n")
            continue
        lines.append(f"## {condition}\n\n- Completed seeds: {len(use)}\n- Mean pMHC CA RMSD: {np.mean([r['pMHC_CA_RMSD_A'] for r in use]):.2f} Å\n- Mean TCR placement CA RMSD: {np.mean([r['TCR_placement_CA_RMSD_A'] for r in use]):.2f} Å\n- Mean TCR fold CA RMSD: {np.mean([r['TCR_fold_CA_RMSD_A'] for r in use]):.2f} Å\n- Mean whole-complex CA RMSD: {np.mean([r['whole_complex_CA_RMSD_A'] for r in use]):.2f} Å\n")
    lines.append("\nThese are structural calibration metrics, not evidence of binding or cross-reactivity.")
    (OUT / "AF3_CALIBRATOR_STRUCTURAL_COMPARISON.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(output, indent=2))

if __name__ == "__main__": main()
