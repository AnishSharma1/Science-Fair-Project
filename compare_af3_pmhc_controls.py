#!/usr/bin/env python3
"""Structural calibration of template-excluded AF3 pMHC-only controls."""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from compare_af3_calibrators import METRICS, REFS, OUT, ca_chains, matched_coords, rmsd

CONTROLS = {
    "PMHC_MBP_DRB1": {"reference": "1YMM", "pairs": {"A":"A", "B":"B", "C":"C"}},
    "PMHC_MBP_DRB5": {"reference": "1ZGL", "pairs": {"A":"A", "B":"B", "C":"C"}},
    "PMHC_ENGA_2WBJ": {"reference": "2WBJ", "pairs": {"A":"A", "B":"B", "C":"D"}},
}

def main():
    rows = json.loads(METRICS.read_text())
    results = []
    for condition, spec in CONTROLS.items():
        reference = ca_chains(REFS / f"{spec['reference']}.cif")
        for row in rows:
            if row.get("condition") != condition or row.get("status") != "completed": continue
            files = list(Path(row["dir"]).glob(f"*_model_{row['bestModel']}.cif"))
            if len(files) != 1: raise ValueError(f"Expected model for {row['name']}")
            predicted = ca_chains(files[0])
            mobile, target = matched_coords(predicted, reference, spec["pairs"], "ABC")
            results.append({"seed_job": row["name"], "condition": condition, "reference": spec["reference"], "CA_atoms": len(mobile), "pMHC_CA_RMSD_A": rmsd(mobile, target)[0], "best_ranking_score": row["bestRanking"], "best_iptm": row["bestIptm"]})
    header = list(results[0])
    def render(value): return f"{value:.3f}" if isinstance(value, float) else str(value)
    (OUT / "af3_pmhc_control_structural_metrics.tsv").write_text("\t".join(header) + "\n" + "\n".join("\t".join(render(row[key]) for key in header) for row in results) + "\n")
    print(json.dumps(results, indent=2))

if __name__ == "__main__": main()
