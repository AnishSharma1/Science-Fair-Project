#!/usr/bin/env python3
"""Compare AF3 Hy.2E11 TCR-peptide contacts across independent seed jobs.

This is a hypothesis-generating interface-consensus analysis only. It does not
validate AF3 docking or establish binding/cross-reactivity.
"""
import json
import shlex
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT = Path("/Users/anishsharma/Documents/New project")
METRICS = ROOT / "outputs/ebv_ms_model_package/results_analysis/af3_seed_job_metrics.json"
OUT = ROOT / "outputs/ebv_ms_model_package/results_analysis"
CONDITIONS = [
    "HY_MBP_DRB1", "HY_BALF5_FULL15", "HY_BALF5_CORE14",
    "DECOY_01_HY_MBP_DRB5", "DECOY_02_HY_ENGA_DRB1",
]
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def atom_loop(path):
    lines = Path(path).read_text().splitlines()
    for start, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        headers, pos = [], start + 1
        while pos < len(lines) and lines[pos].startswith("_atom_site."):
            headers.append(lines[pos].strip()); pos += 1
        if headers:
            index = {x:i for i,x in enumerate(headers)}
            needed = ["_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_seq_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
            if any(x not in index for x in needed):
                raise ValueError(f"Missing atom columns in {path}")
            records = []
            while pos < len(lines):
                line = lines[pos].strip()
                if not line or line.startswith("#"): break
                values = shlex.split(line)
                if len(values) == len(headers):
                    comp = values[index["_atom_site.label_comp_id"]]
                    atom = values[index["_atom_site.label_atom_id"]]
                    if comp in AA3 and not atom.startswith("H"):
                        records.append((values[index["_atom_site.label_asym_id"]], int(values[index["_atom_site.label_seq_id"]]), AA3[comp], np.array([float(values[index["_atom_site.Cartn_x"]]), float(values[index["_atom_site.Cartn_y"]]), float(values[index["_atom_site.Cartn_z"]])])) )
                pos += 1
            return records
    raise ValueError(f"No atom-site loop in {path}")

def contacts(path, cutoff=4.5):
    residues = defaultdict(lambda: {"aa": None, "atoms": []})
    for chain, pos, aa, xyz in atom_loop(path):
        residues[(chain,pos)]["aa"] = aa
        residues[(chain,pos)]["atoms"].append(xyz)
    peptide = [(key, val) for key, val in residues.items() if key[0] == "E"]
    out = set()
    for tcr_chain in ("A", "B"):
        for key, val in residues.items():
            if key[0] != tcr_chain: continue
            t_atoms = np.asarray(val["atoms"])
            for pep_key, pep_val in peptide:
                p_atoms = np.asarray(pep_val["atoms"])
                if np.any(np.sum((t_atoms[:,None,:] - p_atoms[None,:,:]) ** 2, axis=2) <= cutoff ** 2):
                    out.add((key[0], key[1], val["aa"], pep_key[1], pep_val["aa"]))
    return out

def fmt_residue(item):
    chain, pos, aa = item
    return f"{chain}{pos}{aa}"

def main():
    rows = json.loads(METRICS.read_text())
    all_contacts, per_model = defaultdict(list), defaultdict(list)
    for row in rows:
        if row.get("condition") not in CONDITIONS or row.get("status") != "completed": continue
        files = list(Path(row["dir"]).glob(f"*_model_{row['bestModel']}.cif"))
        if len(files) != 1: raise ValueError(f"Missing best model for {row['name']}")
        value = contacts(files[0])
        all_contacts[row["condition"]].append(value)
        per_model[row["condition"]].append((row["name"], value))

    summaries = {}
    for condition, models in all_contacts.items():
        tcr_counter = Counter()
        peptide_counter = Counter()
        pair_counter = Counter()
        for model_contacts in models:
            tcr_counter.update({(a,b,c) for a,b,c,_,_ in model_contacts})
            peptide_counter.update({(d,e) for _,_,_,d,e in model_contacts})
            pair_counter.update(model_contacts)
        summaries[condition] = {
            "models": len(models),
            "mean_tcr_peptide_contacts": round(float(np.mean([len(x) for x in models])), 1),
            "stable_tcr_residues": [fmt_residue(x) for x,n in sorted(tcr_counter.items()) if n == len(models)],
            "stable_peptide_residues": [f"P{x[0]}{x[1]}" for x,n in sorted(peptide_counter.items()) if n == len(models)],
            "stable_contact_pairs": [f"{fmt_residue(x[:3])}-P{x[3]}{x[4]}" for x,n in sorted(pair_counter.items()) if n == len(models)],
            "tcr_contact_frequency": {fmt_residue(x): n for x,n in sorted(tcr_counter.items())},
        }

    def tcr_set(condition):
        return set(summaries[condition]["stable_tcr_residues"])
    comparisons = {}
    for left, right in [("HY_MBP_DRB1", "HY_BALF5_FULL15"), ("HY_MBP_DRB1", "HY_BALF5_CORE14"), ("HY_BALF5_FULL15", "HY_BALF5_CORE14")]:
        a,b = tcr_set(left),tcr_set(right)
        comparisons[f"{left}__{right}"] = {"shared_stable_TCR_residues": sorted(a & b), "Jaccard": round(len(a & b) / len(a | b), 3) if a | b else None}

    payload = {"contact_cutoff_A": 4.5, "summaries": summaries, "Hy_comparisons": comparisons}
    (OUT / "hy_tcr_peptide_contact_consensus.json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = ["# AF3 Hy.2E11 TCR-peptide contact consensus", "", "Heavy-atom contacts use a 4.5 Å cutoff in the best-ranked model from each of five independent seed jobs. A stable residue/contact occurs in all five seed jobs. This is a hypothesis-generating consistency analysis; the template-excluded calibrators did not recover known ternary geometry, so these contacts cannot establish binding, specificity, or cross-reactivity.", ""]
    for condition in CONDITIONS:
        s = summaries[condition]
        lines += [f"## {condition}", "", f"- Models: {s['models']}", f"- Mean TCR-peptide residue-pair contacts/model: {s['mean_tcr_peptide_contacts']}", f"- Stable TCR residues: {', '.join(s['stable_tcr_residues']) or 'none'}", f"- Stable peptide residues: {', '.join(s['stable_peptide_residues']) or 'none'}", f"- Stable TCR–peptide pairs: {', '.join(s['stable_contact_pairs']) or 'none'}", ""]
    lines += ["## Hy.2E11 shared-contact comparison", ""]
    for label, value in comparisons.items():
        lines += [f"- {label}: stable-TCR-residue Jaccard = {value['Jaccard']}; shared = {', '.join(value['shared_stable_TCR_residues']) or 'none'}"]
    (OUT / "AF3_HY_TCR_PEPTIDE_CONTACT_CONSENSUS.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(payload, indent=2))

if __name__ == "__main__": main()
