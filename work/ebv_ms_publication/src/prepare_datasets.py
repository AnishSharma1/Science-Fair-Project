"""Create auditable, normalized datasets for the EBV-MS MHC-II rebuild.

This script never edits the source workbook or legacy project. It writes only
derived tables into processed/. Every row retains its IEDB identifier and the
source URL used for retrieval.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
PROCESSED = ROOT / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

IEDB_XLSX = Path(
    "/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/"
    "Downloads/Review Later/Unsorted Files/downloads/IEDB data.xlsx"
)
IEDB_API = "https://query-api.iedb.org"
RETRIEVED = "2026-08-01"


def clean_antigen(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def clean_source_antigen(item: object) -> dict:
    if not isinstance(item, dict):
        return {}
    return {
        "accession": item.get("accession", ""),
        "name": item.get("name", ""),
        "source_organism_name": item.get("source_organism_name", ""),
        "source_organism_iri": item.get("source_organism_iri", ""),
        "starting_position": item.get("starting_position", ""),
        "ending_position": item.get("ending_position", ""),
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_api(path: Path, assay_type: str) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for item in data:
        antigen = clean_source_antigen(item.get("curated_source_antigen"))
        rows.append(
            {
                "assay_type": assay_type,
                "iedb_assay_id": item.get("tcell_id", item.get("mhc_id", "")),
                "iedb_epitope_id": item.get("structure_id", ""),
                "peptide": str(item.get("linear_sequence") or "").strip().upper(),
                "peptide_length": len(str(item.get("linear_sequence") or "").strip()),
                "source_antigen_accession": antigen.get("accession", ""),
                "source_antigen_name": antigen.get("name", ""),
                "source_organism": antigen.get("source_organism_name", item.get("source_organism_name", "")),
                "source_organism_iri": antigen.get("source_organism_iri", ""),
                "source_start": antigen.get("starting_position", ""),
                "source_end": antigen.get("ending_position", ""),
                "mhc_class": item.get("mhc_class", ""),
                "mhc_allele": item.get("mhc_allele_name", ""),
                "mhc_resolution": item.get("mhc_allele_resolution", ""),
                "mhc_restriction": item.get("mhc_restriction", ""),
                "outcome": item.get("qualitative_measure", ""),
                "assay_description": re.sub(r"<[^>]+>", " ", str(item.get("assay_description", ""))).strip(),
                "host_organism": item.get("host_organism_name", ""),
                "pubmed_id": item.get("pubmed_id", ""),
                "reference_id": item.get("reference_id", ""),
                "source_url": f"{IEDB_API}/{path.stem.split('_')[0]}_search",
                "retrieved_date": RETRIEVED,
            }
        )
    return rows


def main() -> None:
    # Human self-antigen arm from the user's IEDB workbook. The workbook has
    # two header rows; the DRB1*15:01 block begins at column 22.
    sheet = pd.read_excel(IEDB_XLSX, sheet_name="Sheet1", header=None)
    human = sheet.iloc[2:, 22:28].copy()
    human.columns = [
        "iedb_epitope_id",
        "peptide",
        "source_antigen_name",
        "source_organism",
        "n_references",
        "n_assays",
    ]
    human = human.dropna(subset=["iedb_epitope_id", "peptide"])
    human["iedb_epitope_id"] = human["iedb_epitope_id"].astype(int)
    human["peptide"] = human["peptide"].astype(str).str.strip().str.upper()
    human["source_antigen_name"] = human["source_antigen_name"].map(clean_antigen)
    human["source_organism"] = human["source_organism"].map(clean_antigen)
    human["peptide_length"] = human["peptide"].str.len()
    human["mhc_class"] = "II"
    human["mhc_allele"] = "HLA-DRB1*15:01"
    human["evidence_type"] = "IEDB workbook aggregate"
    human["source_file"] = str(IEDB_XLSX)
    human["source_url"] = "https://www.iedb.org/"
    human["retrieved_date"] = RETRIEVED
    human["candidate_class"] = human["source_antigen_name"].str.contains(
        r"Myelin basic protein|Myelin proteolipid protein|Myelin-oligodendrocyte glycoprotein",
        case=False,
        regex=True,
        na=False,
    ).map({True: "myelin_candidate", False: "human_background"})
    human.to_csv(PROCESSED / "human_drb1501_mhc_ii_iedb.csv", index=False)

    api_tables = [
        ("tcell_ebv_drb1501.json", "T-cell response"),
        ("mhc_ebv_drb1501.json", "MHC ligand"),
        ("tcell_human_drb1501.json", "Human comparator T-cell response"),
    ]
    for filename, assay_type in api_tables:
        path = RAW / filename
        if not path.exists():
            continue
        rows = normalize_api(path, assay_type)
        if not rows:
            continue
        fields = list(rows[0])
        out_name = filename.replace(".json", ".csv")
        write_csv(PROCESSED / out_name, rows, fields)

    # Keep a machine-readable provenance record for reproducibility.
    legacy = Path("/Users/anishsharma/Documents/New project/molecular_mimicry_pipeline_v3_fixed.py")
    manifest = {
        "project": "EBV-MS molecular mimicry rebuild",
        "scope": "MHC-II / HLA-DRB1*15:01 primary analysis",
        "retrieved_date": RETRIEVED,
        "iedb_workbook": str(IEDB_XLSX),
        "iedb_api_base": IEDB_API,
        "legacy_code": str(legacy),
        "legacy_code_sha256": hashlib.sha256(legacy.read_bytes()).hexdigest() if legacy.exists() else None,
        "source_tables": sorted(p.name for p in RAW.glob("*.json")),
        "derived_tables": sorted(p.name for p in PROCESSED.glob("*.csv")),
        "status": "rebuild in progress; no publication claims authorized",
    }
    (PROCESSED / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
