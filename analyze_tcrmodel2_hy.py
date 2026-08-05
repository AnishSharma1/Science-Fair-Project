#!/usr/bin/env python3
"""Residue-contact consistency in downloaded TCRmodel2 Hy.2E11 models."""
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT = Path("/Users/anishsharma/Documents/New project")
RESULTS = ROOT / "outputs/ebv_ms_model_package/tcrmodel2_results"
OUT = ROOT / "outputs/ebv_ms_model_package/results_analysis"
CONDITIONS = {"HY_MBP_DRB1": ["HY_MBP_DRB1"], "HY_BALF5": ["HY_BALF5_REPLICATE_1", "HY_BALF5_REPLICATE_2"]}
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def residues(path):
    out = defaultdict(lambda: {"aa":None, "atoms":[]})
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"): continue
        atom, comp, chain, pos = line[12:16].strip(), line[17:20].strip(), line[21], int(line[22:26])
        if comp not in AA3 or atom.startswith("H"): continue
        out[(chain,pos)]["aa"] = AA3[comp]
        out[(chain,pos)]["atoms"].append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    return out

def contacts(path, cutoff=4.5):
    r = residues(path); out = set()
    for (tc,ti), tv in r.items():
        if tc not in {"D","E"}: continue
        x=np.asarray(tv["atoms"])
        for (pc,pi), pv in r.items():
            if pc != "C": continue
            y=np.asarray(pv["atoms"])
            if np.any(np.sum((x[:,None,:]-y[None,:,:])**2, axis=2) <= cutoff**2): out.add((tc,ti,tv["aa"],pi,pv["aa"]))
    return out

def label(x): return f"{x[0]}{x[1]}{x[2]}"

def main():
    models={condition:[contacts(RESULTS / run / "ranked_0.pdb") for run in runs] for condition,runs in CONDITIONS.items()}
    summary={}
    for condition, sets in models.items():
        residue_counts=Counter(); pair_counts=Counter()
        for one in sets:
            residue_counts.update({x[:3] for x in one}); pair_counts.update(one)
        stable={label(x) for x,n in residue_counts.items() if n==len(sets)}
        summary[condition]={"models":len(sets), "mean_pairs":round(float(np.mean([len(x) for x in sets])),1), "stable_TCR_peptide_residues":sorted(stable), "stable_pairs":[f"{label(x[:3])}-P{x[3]}{x[4]}" for x,n in sorted(pair_counts.items()) if n==len(sets)]}
    a,b=set(summary['HY_MBP_DRB1']['stable_TCR_peptide_residues']),set(summary['HY_BALF5']['stable_TCR_peptide_residues'])
    summary['comparison']={"shared_stable_TCR_residues":sorted(a&b), "Jaccard":round(len(a&b)/len(a|b),3)}
    (OUT/'tcrmodel2_hy_contact_consensus.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# TCRmodel2 Hy.2E11 contact consistency','', 'Heavy-atom peptide contacts use a 4.5 A cutoff in ranked_0 models. This summarizes internal model consistency only; it is not evidence of binding or cross-reactivity.','']
    for key in ['HY_MBP_DRB1','HY_BALF5']:
        s=summary[key]; lines += [f'## {key}','',f"- Models: {s['models']}",f"- Mean TCR-peptide residue-pair contacts/model: {s['mean_pairs']}",f"- Stable TCR residues: {', '.join(s['stable_TCR_peptide_residues'])}",f"- Stable pairs: {', '.join(s['stable_pairs'])}",'']
    lines += ['## MBP vs BALF5', '', f"- Shared stable TCR residues: {', '.join(summary['comparison']['shared_stable_TCR_residues'])}", f"- Jaccard overlap: {summary['comparison']['Jaccard']}"]
    (OUT/'TCRMODEL2_HY_CONTACT_CONSENSUS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))

if __name__ == '__main__': main()
