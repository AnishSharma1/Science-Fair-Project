#!/usr/bin/env python3
"""CA structural comparison of downloaded TCRmodel2 calibrator predictions."""
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from compare_af3_calibrators import REFS, OUT, AA, align, matched_coords, rmsd

ROOT = Path("/Users/anishsharma/Documents/New project")
RESULTS = ROOT / "outputs/ebv_ms_model_package/tcrmodel2_results"
SPECS = {
    "CAL_1YMM": {"ref":"1YMM", "pairs":{"A":"A","B":"B","C":"C","D":"D","E":"E"}},
    "CAL_1ZGL_EXCLUDED": {"ref":"1ZGL", "pairs":{"A":"A","B":"B","C":"C","D":"D","E":"E"}},
    "CAL_2WBJ_EXCLUDED": {"ref":"2WBJ", "pairs":{"A":"A","B":"B","C":"D","D":"C","E":"D"}},
}

def pdb_ca_chains(path):
    chains, seen = {}, set()
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM"): continue
        if line[12:16].strip() != "CA": continue
        comp, chain, res = line[17:20].strip(), line[21], line[22:26].strip()
        if comp not in AA or (chain,res) in seen: continue
        seen.add((chain,res))
        chains.setdefault(chain, []).append((AA[comp], np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return chains

def cif_ca_chains(path):
    from compare_af3_calibrators import ca_chains
    return ca_chains(path)

def main():
    output = []
    for condition, spec in SPECS.items():
        pred = pdb_ca_chains(RESULTS / condition / "ranked_0.pdb")
        ref = cif_ca_chains(REFS / f"{spec['ref']}.cif")
        pmhc_pred, pmhc_ref = matched_coords(pred, ref, spec["pairs"], "ABC")
        tcr_pred, tcr_ref = matched_coords(pred, ref, spec["pairs"], "DE")
        all_pred, all_ref = matched_coords(pred, ref, spec["pairs"], "ABCDE")
        pmhc_value, pmhc_fit = rmsd(pmhc_pred, pmhc_ref)
        placement, _ = rmsd(tcr_pred, tcr_ref, pmhc_fit)
        fold, _ = rmsd(tcr_pred, tcr_ref)
        whole, _ = rmsd(all_pred, all_ref)
        output.append((condition, spec["ref"], len(pmhc_pred), pmhc_value, len(tcr_pred), placement, fold, whole))
    header = "condition\treference\tpMHC_CA_atoms\tpMHC_CA_RMSD_A\tTCR_CA_atoms\tTCR_placement_CA_RMSD_A\tTCR_fold_CA_RMSD_A\twhole_complex_CA_RMSD_A"
    rows = [header] + ["\t".join([a,b,str(c),f"{d:.3f}",str(e),f"{f:.3f}",f"{g:.3f}",f"{h:.3f}"]) for a,b,c,d,e,f,g,h in output]
    (OUT / "tcrmodel2_calibrator_structural_metrics.tsv").write_text("\n".join(rows) + "\n")
    print("\n".join(rows))

if __name__ == "__main__": main()
