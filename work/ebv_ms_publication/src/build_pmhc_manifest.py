"""Build the candidate manifest for the upcoming pMHC modeling phase."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"


def main() -> None:
    rows = []
    ebv_t = pd.read_csv(PROC / "tcell_ebv_drb1501.csv").drop_duplicates("peptide")
    positive = {"Positive", "Positive-Low", "Positive-Intermediate", "Positive-High"}
    for _, r in ebv_t[ebv_t.outcome.isin(positive)].iterrows():
        rows.append({
            "candidate_id": f"EBV_TCELL_{int(r.iedb_epitope_id)}",
            "arm": "EBV",
            "evidence_tier": "Tier 1: DRB1*15:01 T-cell assay",
            "peptide": r.peptide,
            "peptide_length": int(r.peptide_length),
            "source_antigen": r.source_antigen_name,
            "source_accession": r.source_antigen_accession,
            "iedb_assay_id": r.iedb_assay_id,
            "iedb_epitope_id": r.iedb_epitope_id,
            "pubmed_id": r.pubmed_id,
            "hla": r.mhc_allele,
            "mhc_class": r.mhc_class,
            "modeling_status": "eligible_after input/structure QA",
        })
    ebv_m = pd.read_csv(PROC / "mhc_ebv_drb1501.csv").drop_duplicates("peptide")
    for _, r in ebv_m.iterrows():
        rows.append({
            "candidate_id": f"EBV_MHC_{int(r.iedb_epitope_id)}",
            "arm": "EBV",
            "evidence_tier": "Tier 2: DRB1*15:01 MHC-ligand assay",
            "peptide": r.peptide,
            "peptide_length": int(r.peptide_length),
            "source_antigen": r.source_antigen_name,
            "source_accession": r.source_antigen_accession,
            "iedb_assay_id": r.iedb_assay_id,
            "iedb_epitope_id": r.iedb_epitope_id,
            "pubmed_id": r.pubmed_id,
            "hla": r.mhc_allele,
            "mhc_class": r.mhc_class,
            "modeling_status": "eligible after input/structure QA",
        })
    human = pd.read_csv(PROC / "human_drb1501_mhc_ii_iedb_enriched.csv")
    human = human[(human.candidate_class == "myelin_candidate") & (human.provenance_status == "coordinate_validated")]
    for _, r in human.drop_duplicates("peptide").iterrows():
        rows.append({
            "candidate_id": f"HUMAN_MYELIN_{int(r.iedb_epitope_id)}",
            "arm": "Human myelin",
            "evidence_tier": "Tier 3: HLA-DRB1*15:01 aggregate epitope record",
            "peptide": r.peptide,
            "peptide_length": int(r.peptide_length),
            "source_antigen": r.source_antigen_name,
            "source_accession": "see human_epitope_accession_map.csv",
            "iedb_assay_id": "",
            "iedb_epitope_id": r.iedb_epitope_id,
            "pubmed_id": "",
            "hla": r.mhc_allele,
            "mhc_class": r.mhc_class,
            "modeling_status": "eligible after input/structure QA",
        })
    out = pd.DataFrame(rows).drop_duplicates(subset=["arm", "peptide", "evidence_tier"])
    out.to_csv(PROC / "pmhc_candidate_manifest.csv", index=False)
    with (PROC / "pmhc_candidate_peptides.fasta").open("w", encoding="utf-8") as handle:
        for _, r in out.iterrows():
            handle.write(f">{r.candidate_id}|{r.arm.replace(' ', '_')}|{r.hla}|IEDB_EPITOPE:{r.iedb_epitope_id}\n{r.peptide}\n")
    print(out.groupby(["arm", "evidence_tier"]).size().to_string())


if __name__ == "__main__":
    main()
