#!/usr/bin/env python3
"""Attach biologically readable protein/residue labels to pMHC audit tables."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path


PROCESSED = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent
AUDIT_DIR = PROCESSED / "complete_model_pipeline_audit_2026-08-12"


PROTEIN_NAMES = {
    "AAA45887.1": ("EBV membrane protein", "EBV membrane protein"),
    "CAD53462.1": ("BALF5 DNA polymerase", "BALF5"),
    "CEQ37348.1": ("BALF4 glycoprotein B", "BALF4 gB"),
    "P03182.1": ("BHRF1 apoptosis regulator", "BHRF1"),
    "P03186.1": ("BPLF1 large tegument protein deneddylase", "BPLF1"),
    "P03187.1": ("Triplex capsid protein 1", "Triplex capsid protein 1"),
    "P03200.1": ("BLLF1 envelope glycoprotein gp350", "gp350"),
    "Q3KSQ3.1": ("BXLF2 envelope glycoprotein H", "gH"),
    "Q3KST1.1": ("EBNA4 (EBNA3B)", "EBNA4"),
    "Q777A4": ("Latent membrane protein 1", "LMP1"),
    "Q777E1": ("Epstein-Barr nuclear antigen 1", "EBNA1"),
    "YP_401677.1": ("Epstein-Barr nuclear antigen 1", "EBNA1"),
    "P02686": ("Myelin basic protein", "MBP"),
    "P60201": ("Myelin proteolipid protein", "PLP1"),
    "Q16653": ("Myelin-oligodendrocyte glycoprotein", "MOG"),
    "AAB08089.1": ("Myelin-oligodendrocyte glycoprotein", "MOG"),
    "P01889.3": ("HLA class I histocompatibility antigen B alpha chain", "HLA-B"),
    "O95167.1": ("NADH dehydrogenase subunit NDUFA3", "NDUFA3"),
    "P23528.3": ("Cofilin-1", "CFL1"),
    "P68871.2": ("Hemoglobin subunit beta", "HBB"),
    "EAW90170.1": ("Sentrin-specific protease 3", "SENP3"),
    "P07197.3": ("Neurofilament medium polypeptide", "NEFM"),
    "P69905.2": ("Hemoglobin subunit alpha", "HBA"),
    "P05090.1": ("Apolipoprotein D", "APOD"),
    "Q14CZ8": ("Glial cell adhesion molecule", "GlialCAM"),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = ""
    sequence: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if header:
                records[header] = "".join(sequence)
            header = line[1:]
            sequence = []
        elif line:
            sequence.append(line)
    if header:
        records[header] = "".join(sequence)
    return records


def header_has_accession(header: str, accession: str) -> bool:
    bare = accession.split(".")[0]
    first = header.split()[0]
    return first == accession or first.split("|")[1:2] == [accession] or bare in first.split("|")


def exact_starts(sequence: str, peptide: str) -> list[int]:
    return [match.start() + 1 for match in re.finditer(f"(?={re.escape(peptide)})", sequence)]


def parent_name(accession: str, fallback: str) -> tuple[str, str]:
    if accession in PROTEIN_NAMES:
        return PROTEIN_NAMES[accession]
    bare = accession.split(".")[0]
    if bare in PROTEIN_NAMES:
        return PROTEIN_NAMES[bare]
    return fallback, fallback.replace(" ", "-")


def build_annotations() -> list[dict[str, str]]:
    manifest = read_csv(PROCESSED / "pmhc_candidate_manifest.csv")
    human_map = read_csv(PROCESSED / "human_epitope_accession_map.csv")
    human_by_epitope: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in human_map:
        human_by_epitope[row["iedb_epitope_id"]].append(row)

    ebv_fasta = read_fasta(OUT_DIR / "source_records" / "ebv_parent_proteins.fasta")
    human_fasta = read_fasta(OUT_DIR / "source_records" / "human_parent_proteins_uniprot.fasta")
    annotations: list[dict[str, str]] = []

    for row in manifest:
        peptide = row["peptide"]
        recorded_accession = row["source_accession"]
        accession = recorded_accession
        start = end = None
        match_count = 0
        coordinate_source = ""
        mapping_status = ""

        if row["arm"] == "EBV":
            matching_sequences = [
                sequence
                for header, sequence in ebv_fasta.items()
                if header_has_accession(header, accession)
            ]
            starts = [position for sequence in matching_sequences for position in exact_starts(sequence, peptide)]
            match_count = len(starts)
            if match_count == 1:
                start = starts[0]
                end = start + len(peptide) - 1
                mapping_status = "unique_exact_match_to_recorded_accession"
            else:
                mapping_status = "multiple_exact_matches" if match_count > 1 else "not_found"
            coordinate_source = "NCBI Protein accession sequence"
        else:
            canonical_match = re.search(r"UniProt:([A-Z0-9]+)", row["source_antigen"])
            if not canonical_match:
                raise ValueError(f"No human UniProt accession in {row['candidate_id']}")
            accession = canonical_match.group(1)
            matching_sequences = [
                sequence
                for header, sequence in human_fasta.items()
                if header_has_accession(header, accession)
            ]
            starts = [position for sequence in matching_sequences for position in exact_starts(sequence, peptide)]
            match_count = len(starts)
            if match_count == 1:
                start = starts[0]
                end = start + len(peptide) - 1
                mapping_status = "unique_exact_match_to_canonical_uniprot"
                coordinate_source = "UniProt canonical sequence"
            else:
                # One MOG peptide belongs to a shorter accession/isoform rather than canonical Q16653.
                candidates = [
                    item
                    for item in human_by_epitope[row["iedb_epitope_id"]]
                    if item["peptide"] == peptide and item["mapping_status"] == "coordinate_validated"
                ]
                preferred = next((item for item in candidates if item["accession"] == "AAB08089.1"), None)
                if preferred is None and candidates:
                    preferred = candidates[0]
                if preferred is None:
                    mapping_status = "not_found"
                else:
                    accession = preferred["accession"]
                    start = int(preferred["start"])
                    end = int(preferred["end"])
                    match_count = 1
                    mapping_status = "IEDB_coordinate_validated_noncanonical_accession"
                    coordinate_source = "IEDB human epitope accession map"

        if start is None or end is None:
            raise ValueError(f"Unresolved protein coordinates for {row['candidate_id']}")
        protein, short = parent_name(accession, row["source_antigen"].split(" (UniProt:")[0])
        region = f"{short} residues {start}-{end}"
        readable = f"{region} ({peptide})"
        annotations.append(
            {
                "candidate_id": row["candidate_id"],
                "readable_candidate_name": readable,
                "protein_region_label": region,
                "parent_protein_name": protein,
                "short_protein_name": short,
                "resolved_parent_accession": accession,
                "parent_residue_start_1_based": str(start),
                "parent_residue_end_1_based": str(end),
                "peptide": peptide,
                "arm": row["arm"],
                "source_antigen_recorded": row["source_antigen"],
                "source_accession_recorded": recorded_accession,
                "mapping_status": mapping_status,
                "exact_match_count": str(match_count),
                "coordinate_source": coordinate_source,
                "stable_id": row["candidate_id"],
            }
        )

    if len(annotations) != 86 or len({row["candidate_id"] for row in annotations}) != 86:
        raise ValueError("Expected exactly 86 unique candidate annotations")

    # Add the matched-background candidates that appear in the AF3 audit inventory.
    background_manifest = read_csv(PROCESSED / "expanded_background" / "background_pmhc_candidate_manifest.csv")
    for row in background_manifest:
        peptide = row["peptide"]
        canonical = re.search(r"UniProt:([A-Z0-9]+)", row["source_antigen"])
        candidates = [
            item
            for item in human_by_epitope[row["iedb_epitope_id"]]
            if item["peptide"] == peptide and item["mapping_status"] == "coordinate_validated"
        ]
        preferred = []
        if canonical:
            preferred = [item for item in candidates if item["accession"].split(".")[0] == canonical.group(1)]
        selected_pool = preferred or [item for item in candidates if item["accession_type"] == "UniProt"] or candidates
        if not selected_pool:
            raise ValueError(f"No validated background coordinate for {row['candidate_id']}")
        selected = selected_pool[0]
        accession = selected["accession"]
        start, end = int(selected["start"]), int(selected["end"])
        protein_fallback = row["source_antigen"].split(" (UniProt:")[0]
        protein, short = parent_name(accession, protein_fallback)
        region = f"{short} residues {start}-{end}"
        annotations.append(
            {
                "candidate_id": row["candidate_id"],
                "readable_candidate_name": f"{region} ({peptide})",
                "protein_region_label": region,
                "parent_protein_name": protein,
                "short_protein_name": short,
                "resolved_parent_accession": accession,
                "parent_residue_start_1_based": str(start),
                "parent_residue_end_1_based": str(end),
                "peptide": peptide,
                "arm": row["arm"],
                "source_antigen_recorded": row["source_antigen"],
                "source_accession_recorded": row["source_accession"],
                "mapping_status": "IEDB_coordinate_validated_background_accession",
                "exact_match_count": "1",
                "coordinate_source": "IEDB human epitope accession map",
                "stable_id": row["candidate_id"],
            }
        )

    # The four special controls use one GlialCAM segment with different phosphoserine sites.
    for phosphosite in (376, 377, 383, 384):
        candidate_id = f"GLIALCAM_370_389_PS{phosphosite}"
        peptide = "ATGRTHSSPPRAPSSPGRSR"
        region = f"GlialCAM residues 370-389, pSer{phosphosite}"
        annotations.append(
            {
                "candidate_id": candidate_id,
                "readable_candidate_name": f"{region} ({peptide})",
                "protein_region_label": region,
                "parent_protein_name": "Glial cell adhesion molecule",
                "short_protein_name": "GlialCAM",
                "resolved_parent_accession": "Q14CZ8",
                "parent_residue_start_1_based": "370",
                "parent_residue_end_1_based": "389",
                "peptide": peptide,
                "arm": "Special GlialCAM control",
                "source_antigen_recorded": "GlialCAM phosphopeptide control",
                "source_accession_recorded": "Q14CZ8",
                "mapping_status": "predefined_control_coordinates",
                "exact_match_count": "1",
                "coordinate_source": "Predefined AlphaFold request label",
                "stable_id": candidate_id,
            }
        )

    annotations.append(
        {
            "candidate_id": "DECOY_02_HY_ENGA_DRB1_S101",
            "readable_candidate_name": "Excluded non-pMHC five-chain decoy (no peptide region)",
            "protein_region_label": "No single protein/peptide region",
            "parent_protein_name": "",
            "short_protein_name": "",
            "resolved_parent_accession": "",
            "parent_residue_start_1_based": "",
            "parent_residue_end_1_based": "",
            "peptide": "",
            "arm": "Excluded control",
            "source_antigen_recorded": "",
            "source_accession_recorded": "",
            "mapping_status": "excluded_non_pmhc_chain_layout",
            "exact_match_count": "0",
            "coordinate_source": "AF3 inventory classification",
            "stable_id": "DECOY_02_HY_ENGA_DRB1_S101",
        }
    )

    if len({row["candidate_id"] for row in annotations}) != len(annotations):
        raise ValueError("Duplicate candidate IDs in combined annotation table")
    return annotations


def annotate_candidate_table(source: Path, target: Path, by_id: dict[str, dict[str, str]]) -> None:
    output: list[dict[str, str]] = []
    for row in read_csv(source):
        ann = by_id[row["candidate_id"]]
        output.append(
            {
                "candidate_id": row["candidate_id"],
                "readable_candidate_name": ann["readable_candidate_name"],
                "protein_region_label": ann["protein_region_label"],
                "parent_protein_name": ann["parent_protein_name"],
                "resolved_parent_accession": ann["resolved_parent_accession"],
                **{key: value for key, value in row.items() if key != "candidate_id"},
            }
        )
    write_csv(target, output)


def annotate_pair_table(source: Path, target: Path, by_id: dict[str, dict[str, str]]) -> None:
    output: list[dict[str, str]] = []
    for row in read_csv(source):
        ebv = by_id[row["ebv_candidate_id"]]
        human = by_id[row["human_candidate_id"]]
        pair_name = f"{ebv['protein_region_label']} vs {human['protein_region_label']}"
        leading: dict[str, str] = {}
        if "discovery_priority_rank" in row:
            leading["discovery_priority_rank"] = row["discovery_priority_rank"]
        leading.update(
            {
                "pair_id": row["pair_id"],
                "readable_pair_name": pair_name,
                "ebv_readable_candidate_name": ebv["readable_candidate_name"],
                "human_readable_candidate_name": human["readable_candidate_name"],
                "ebv_parent_protein": ebv["parent_protein_name"],
                "ebv_parent_accession": ebv["resolved_parent_accession"],
                "ebv_parent_residues_1_based": f"{ebv['parent_residue_start_1_based']}-{ebv['parent_residue_end_1_based']}",
                "human_parent_protein": human["parent_protein_name"],
                "human_parent_accession": human["resolved_parent_accession"],
                "human_parent_residues_1_based": f"{human['parent_residue_start_1_based']}-{human['parent_residue_end_1_based']}",
            }
        )
        skip = set(leading)
        output.append({**leading, **{key: value for key, value in row.items() if key not in skip}})
    write_csv(target, output)


def main() -> None:
    annotations = build_annotations()
    by_id = {row["candidate_id"]: row for row in annotations}
    write_csv(OUT_DIR / "candidate_protein_region_annotations.csv", annotations)

    annotate_pair_table(
        AUDIT_DIR / "01_COMPLETE_32_PAIR_SCORECARD.csv",
        AUDIT_DIR / "01_COMPLETE_32_PAIR_SCORECARD_WITH_PROTEIN_REGIONS.csv",
        by_id,
    )
    annotate_pair_table(
        AUDIT_DIR / "02_ALL_1000_STRUCTURE_COMPARISONS.csv",
        AUDIT_DIR / "02_ALL_1000_STRUCTURE_COMPARISONS_WITH_PROTEIN_REGIONS.csv",
        by_id,
    )
    annotate_candidate_table(
        AUDIT_DIR / "03_ALL_150_SAVED_AF3_JOB_FOLDERS.csv",
        AUDIT_DIR / "03_ALL_150_SAVED_AF3_JOB_FOLDERS_WITH_PROTEIN_REGIONS.csv",
        by_id,
    )
    annotate_candidate_table(
        AUDIT_DIR / "04_UNIQUE_AF3_JOB_QUALITY_SUMMARY.csv",
        AUDIT_DIR / "04_UNIQUE_AF3_JOB_QUALITY_SUMMARY_WITH_PROTEIN_REGIONS.csv",
        by_id,
    )
    annotate_candidate_table(
        AUDIT_DIR / "05_INDIVIDUAL_AF3_MODEL_QUALITY_METRICS.csv",
        AUDIT_DIR / "05_INDIVIDUAL_AF3_MODEL_QUALITY_METRICS_WITH_PROTEIN_REGIONS.csv",
        by_id,
    )


if __name__ == "__main__":
    main()
