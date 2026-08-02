"""Validate human epitope provenance using IEDB source-antigen metadata."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "processed"
RAW = ROOT / "raw"


def main() -> None:
    human = pd.read_csv(PROC / "human_drb1501_mhc_ii_iedb.csv")
    wanted = set(human.iedb_epitope_id.astype(int))
    metadata = json.loads((RAW / "human_epitope_metadata.json").read_text())
    rows = []
    for record in metadata:
        eid = int(record["structure_id"])
        if eid not in wanted:
            continue
        peptide = (record.get("structure_descriptions") or [""])[0]
        sources = record.get("curated_source_antigens") or []
        human_sources = [s for s in sources if "Homo sapiens" in str(s.get("source_organism_name"))]
        for source in human_sources:
            start = source.get("starting_position")
            end = source.get("ending_position")
            coord_len = int(end - start + 1) if start is not None and end is not None else None
            rows.append(
                {
                    "iedb_epitope_id": eid,
                    "peptide": peptide,
                    "accession": source.get("accession", ""),
                    "source_antigen_name": source.get("name", ""),
                    "source_organism": source.get("source_organism_name", ""),
                    "start": start,
                    "end": end,
                    "coordinate_length": coord_len,
                    "peptide_length": len(peptide),
                    "coordinate_length_matches": coord_len == len(peptide),
                    "mapping_status": "coordinate_validated" if coord_len == len(peptide) else "quarantine_missing_coordinates",
                    "accession_type": "UniProt" if re.fullmatch(r"[A-Z][0-9][A-Z0-9]{3}[0-9](?:\.\d+)?", str(source.get("accession", ""))) else "sequence database",
                    "source_url": f"https://query-api.iedb.org/epitope_search?structure_id=eq.{eid}",
                }
            )
    fields = list(rows[0])
    with (PROC / "human_epitope_accession_map.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    status_by_id = {}
    for eid in wanted:
        matching = [r for r in rows if r["iedb_epitope_id"] == eid]
        status_by_id[eid] = (
            "coordinate_validated"
            if any(r["coordinate_length_matches"] for r in matching)
            else "quarantine_missing_coordinates"
        )
    enriched = human.copy()
    enriched["provenance_status"] = enriched["iedb_epitope_id"].map(status_by_id)
    enriched.to_csv(PROC / "human_drb1501_mhc_ii_iedb_enriched.csv", index=False)

    counts = Counter(r["iedb_epitope_id"] for r in rows)
    report = [
        "# Human epitope provenance validation",
        "",
        f"- Workbook epitopes checked: **{len(wanted)}**",
        f"- Epitopes with at least one Homo sapiens source mapping: **{len(counts)}**",
        f"- Human source mappings retained: **{len(rows)}**",
        f"- Coordinate-length mismatches: **{sum(not r['coordinate_length_matches'] for r in rows)}**",
        f"- Epitopes quarantined for missing coordinates: **{sum(v == 'quarantine_missing_coordinates' for v in status_by_id.values())}**",
        f"- Epitopes with multiple human source mappings: **{sum(v > 1 for v in counts.values())}**",
        "",
        "Every retained mapping has a source organism of Homo sapiens. Most have",
        "a coordinate span equal to the reported peptide length; the mismatches",
        "must be reviewed rather than silently discarded. Multiple mappings",
        "are preserved because the same epitope may be represented by several",
        "sequence accessions or isoform records. The analysis should use the",
        "mapping table rather than silently choosing one accession.",
        "",
        "This validates provenance and coordinates; it does not prove that the",
        "peptide is disease-specific, pathogenic, or naturally presented in MS.",
    ]
    (PROC / "human_epitope_validation.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("epitopes", len(wanted), "human mappings", len(rows), "mismatches", sum(not r["coordinate_length_matches"] for r in rows))


if __name__ == "__main__":
    main()
