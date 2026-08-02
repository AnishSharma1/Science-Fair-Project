"""Write the explicit pMHC modeling input specification."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"


def main() -> None:
    candidates = pd.read_csv(PROC / "pmhc_candidate_manifest.csv")
    candidates = candidates.to_dict(orient="records")
    spec = {
        "template": {
            "pdb_id": "8TBP",
            "source_url": "https://www.rcsb.org/structure/8TBP",
            "alpha_chains": ["A", "C"],
            "beta_chains": ["B", "D"],
            "peptide_chains": ["E", "F"],
            "hla": "HLA-DRB1*15:01",
            "method": "X-ray diffraction",
        },
        "target_hla": "HLA-DRB1*15:01",
        "target_mhc_class": "II",
        "modeling_policy": {
            "peptide_sequence_is_full_iedb_sequence": True,
            "do_not_trim_to_class_i_9mer": True,
            "preserve_candidate_id_and_iedb_id": True,
            "predicted_models_are_hypothesis_generating": True,
            "tcr_docking_is_out_of_scope_for_this_stage": True,
        },
        "structure_generation": {
            "input_format": "ColabFold multimer FASTA",
            "all_candidate_inputs": "processed/pmhc_colabfold_inputs.fasta",
            "per_candidate_inputs": "processed/colabfold_inputs/{candidate_id}.fasta",
            "chain_order": "mature HLA-DRA : mature HLA-DRB1*15:01 : full IEDB peptide",
            "template_chain_source": "8TBP chains A and B",
            "coordinates_generated": False,
        },
        "candidates": candidates,
    }
    (PROC / "pmhc_modeling_inputs.json").write_text(json.dumps(spec, indent=2), encoding="utf-8")
    print("candidates", len(candidates), "template", spec["template"]["pdb_id"])


if __name__ == "__main__":
    main()
